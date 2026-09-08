"""Audit saved scores and reproduce final assembly. No predictor inference or training.

The Borda/MMR/apportion functions are unchanged extractions of the consolidation.
The average-rank helper replaces scipy.stats.rankdata for this fixed environment.
Clustering is frozen: fitted scaler, centers and labels are supplied as artifacts.
"""
from pathlib import Path
import csv,gzip,hashlib,json
import numpy as np
from rapidfuzz import process
from rapidfuzz.distance import Levenshtein as RFLev
ROOT=Path(__file__).resolve().parent
def rankdata(values,method='average'):
    if method!='average':raise ValueError(method)
    a=np.asarray(values);order=np.argsort(a,kind='stable');sorted_values=a[order]
    starts=np.r_[0,np.flatnonzero(sorted_values[1:]!=sorted_values[:-1])+1];ends=np.r_[starts[1:],len(a)]
    result=np.empty(len(a),dtype=float)
    result[order]=np.repeat((starts+ends+1)/2,ends-starts)
    return result
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))

def borda(values,specs):
    values=np.asarray(values,dtype=float)
    if values.ndim!=2 or values.shape[1]!=len(specs) or not np.isfinite(values).all():raise ValueError('Missing/nonfinite/inconsistent panel')
    n=len(values);ranks=np.full(values.shape,np.nan);bs=np.full(values.shape,np.nan);active=[];omitted=[];seen={}
    for j,e in enumerate(specs):
        key=(e.get('implementation',e['id']),e.get('checkpoint'),e['direction'])
        if key in seen:omitted.append({'id':e['id'],'reason':'identical evaluator implementation/checkpoint copy','copy_of':seen[key]});continue
        seen[key]=e['id']
        if n>1 and np.all(values[:,j]==values[0,j]):omitted.append({'id':e['id'],'reason':'constant evaluator; no ranking information'});continue
        ranks[:,j]=rankdata(values[:,j] if e['direction']=='minimize' else -values[:,j],method='average')
        bs[:,j]=1. if n==1 else (n-ranks[:,j])/(n-1);active.append(j)
    families={}
    for name in sorted({specs[j]['family'] for j in active}):
        ids=[j for j in active if specs[j]['family']==name];families[name]=bs[:,ids].mean(axis=1)
    q=np.mean(list(families.values()),axis=0) if families else np.zeros(n)
    return ranks,bs,families,q,active,omitted

def mmr(indices,seqs,q,count,lam):
    ids=np.asarray(indices,dtype=np.int64);strings=[seqs[i] for i in ids];qsub=q[ids];maximum=np.zeros(len(ids));available=np.ones(len(ids),bool);chosen=[];events=[]
    for step in range(min(count,len(ids))):
        utility=lam*qsub-(1-lam)*maximum
        # IDs are already ordered by Q then sequence; argmax's first match implements the prescribed secondary ties.
        utility[~available]=-np.inf;best=int(np.argmax(qsub if step==0 else utility))
        assert available[best]
        chosen.append(int(ids[best]));events.append((int(ids[best]),step+1,float(maximum[best]),float(utility[best])))
        available[best]=False
        sims=process.cdist([strings[best]],strings,scorer=RFLev.normalized_similarity,dtype=np.float64,workers=1)[0]
        maximum=np.maximum(maximum,sims)
    return chosen,events

def apportion(capacities,slots):
    cap=np.asarray(capacities,dtype=np.int64);slots=min(int(slots),int(cap.sum()))
    if not len(cap) or not cap.sum():return np.zeros_like(cap)
    # Integer arithmetic avoids floating tie ambiguity in largest remainder.
    numerator=cap*slots;alloc=numerator//cap.sum();remainder=numerator%cap.sum();left=slots-int(alloc.sum())
    order=sorted(range(len(cap)),key=lambda i:(-int(remainder[i]),i))
    for i in order:
        if left and alloc[i]<cap[i]:alloc[i]+=1;left-=1
    assert left==0 and int(alloc.sum())==slots and np.all(alloc<=cap)
    return alloc

def main(generated_sequences=None):
    for rel,h in read(ROOT/'data_manifest.json').items():
        if sha(ROOT/rel)!=h:raise ValueError('Changed score/selection evidence: '+rel)
    with gzip.open(ROOT/'data/pool_sequences.txt.gz','rt',encoding='ascii') as f:seqs=f.read().splitlines()
    assert seqs==sorted(set(seqs)) and len(seqs)==499998
    if generated_sequences is not None and seqs!=generated_sequences:
        raise ValueError('Frozen score-table population differs from regenerated source pools')
    raw=np.load(ROOT/'data/raw_scores.npy');cache=np.load(ROOT/'data/rank_cache.npz')
    r,b,f,q,active,omitted=borda(raw,read(ROOT/'panel_config.json')['evaluators'])
    order=np.lexsort((np.array(seqs),-q));ranks=np.empty(len(q),dtype=np.int64);ranks[order]=np.arange(1,len(q)+1)
    assert np.array_equal(q,cache['Q']) and np.array_equal(order,cache['quality_order']) and np.array_equal(ranks,cache['rank'])
    med=np.lexsort((np.array(seqs),np.median(r[:,active],axis=1)))
    assert np.array_equal(med,cache['median_order'])
    print('All 499,998 consensus scores and ranks match the frozen ranking.',flush=True)
    cfg=read(ROOT/'selection_config.json');screen=np.load(ROOT/'data/mmr_screen.npy');labels=np.load(ROOT/'data/cluster_labels.npy')
    assert np.array_equal(screen,np.array(sorted(screen,key=lambda i:(-q[i],seqs[i]))))
    chosen,events=mmr(screen,seqs,q,200,cfg['mmr_lambda']);expected=read(ROOT/'selection_record.json')
    assert chosen==expected['mmr_insertion_pool_indices']
    top=sorted(chosen[:100],key=lambda i:ranks[i]);mask=np.zeros(len(seqs),bool);mask[top]=True;k=int(labels.max())+1
    allocation=apportion(np.bincount(labels[~mask],minlength=k),50000-len(top));taken=np.zeros(k,dtype=np.int64)
    for i in order:
        c=labels[i]
        if not mask[i] and taken[c]<allocation[c]:mask[i]=True;taken[c]+=1
    lib=order[mask[order]].tolist();assert len(lib)==50000 and len(top)==100 and set(top)<=set(lib)
    output=ROOT/'replayed';output.mkdir(exist_ok=True)
    hashes={}
    for name,ids in [('library',lib),('top',top)]:
        payload=''.join(f'>amp_cons_q{ranks[i]:06d}\n{seqs[i]}\n' for i in ids).encode('ascii')
        hashes[name]=hashlib.sha256(payload).hexdigest();assert hashes[name]==expected[name+'_sha256']
        (output/(name+'.fasta')).write_bytes(payload)
    # Persist per-evaluator tied ranks/Borda and family contributions on demand.
    np.savez_compressed(output/'expanded_rank_table.npz',raw=raw,evaluator_ranks=r,evaluator_borda=b,
                        **{'family_'+name:value for name,value in f.items()},Q=q,quality_rank=ranks)
    report={'all_scores_and_quality_ranks_exact':True,'all_mmr_insertion_indices_exact':True,
            'cluster_allocation_exact':True,'output_sha256':hashes,'omitted_evaluators':omitted,
            'predictors_executed':False,'clustering_refitted':False}
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('MMR, proportional selection and both artifact bytes match.',flush=True)
if __name__=='__main__':main()
