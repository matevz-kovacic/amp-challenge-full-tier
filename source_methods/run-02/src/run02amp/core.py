"""Seeded conditional sequence generation and explicitly provisional ranking."""
import argparse,csv,hashlib,json,math,os,time
from collections import Counter
from pathlib import Path
for _key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[_key]='2'
import numpy as np
from rapidfuzz import process
from rapidfuzz.distance import Indel
from .features import ALPHABET,HYDROPHOBIC,esm_embeddings

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read_fasta(path):
    rows=[];seq=[]
    for line in Path(path).read_text('utf-8').splitlines():
        line=line.strip()
        if not line:continue
        if line.startswith('>'):
            if seq:rows.append(''.join(seq).upper());seq=[]
        else:seq.append(line)
    if seq:rows.append(''.join(seq).upper())
    return rows
def properties(s):
    n=len(s);c=Counter(s);h=[int(a in HYDROPHOBIC) for a in s]
    longest=run=samerun=samemax=0;prev=None
    for a,v in zip(s,h):
        run=run+1 if v else 0;longest=max(longest,run)
        samerun=samerun+1 if a==prev else 1;samemax=max(samemax,samerun);prev=a
    return {'length':n,'charge_proxy':c['K']+c['R']-c['D']-c['E'],
      'hydrophobic_fraction':sum(h)/n,'basic_fraction':(c['K']+c['R'])/n,
      'proline_count':c['P'],'glycine_fraction':c['G']/n,'hydrophobic_run':longest,
      'identical_run':samemax,'entropy':-sum((v/n)*math.log(v/n) for v in c.values())/math.log(20),
      'binary_helix_moment':abs(sum(v*complex(math.cos(math.radians(100*i)),math.sin(math.radians(100*i))) for i,v in enumerate(h)))/n}
def prior_pass(s,p):
    return (12<=p['length']<=32 and 'C' not in s and 3<=p['charge_proxy']<=8
        and .30<=p['hydrophobic_fraction']<=.55 and p['basic_fraction']<=.45
        and p['proline_count']<=1 and p['glycine_fraction']<=.20
        and p['hydrophobic_run']<=4 and p['identical_run']<=3)

class Operator:
    def __init__(self,model,seed,regime):
        self.rng=np.random.Generator(np.random.PCG64(seed));self.model=model;self.regime=regime
        self.base=np.array(model['unigram'],dtype=np.float64)
        self.cdf=np.cumsum(np.array(model['conditional'],dtype=np.float64),axis=-1)
        self.lengthp=np.array(model['length_probabilities'],dtype=np.float64)
    def batch(self,size=2048):
        rng=self.rng;lengths=rng.choice(np.arange(12,33),size=size,p=self.lengthp)
        families=np.full(size,self.regime,dtype='U12')
        if self.regime=='mixture':families=np.where(rng.random(size)<.5,'markov','grammar')
        phases=rng.uniform(0,2*np.pi,size);strength=rng.uniform(1.5,3,size)
        tokens=np.zeros((size,32),dtype=np.int16);before=np.full(size,20);last=before.copy()
        hyd=np.array([a in HYDROPHOBIC for a in ALPHABET],float)
        basic=np.array([a in 'KR' for a in ALPHABET],float)
        for pos in range(32):
            cdf=self.cdf[before,last].copy()
            comp=families=='composition'
            if comp.any():cdf[comp]=np.cumsum(self.base)
            gram=families=='grammar'
            if gram.any():
                face=np.cos(phases[gram]+math.radians(100)*pos)
                weights=self.base[None,:]*np.exp(strength[gram,None]*face[:,None]*(hyd[None,:]-.7*basic[None,:]))
                cdf[gram]=np.cumsum(weights/weights.sum(1,keepdims=True),axis=1)
            cdf[:,-1]=1.0
            chosen=(rng.random(size)[:,None]>cdf).sum(axis=1)
            tokens[:,pos]=chosen;before,last=last,chosen
        return [(''.join(ALPHABET[a] for a in tokens[i,:n]),str(families[i])) for i,n in enumerate(lengths)]

def generate_library(model,regime,seed,count,excluded):
    operator=Operator(model,seed,regime);kept=[];seen=set();draws=0;fail=Counter();families=Counter()
    while len(kept)<count:
        if draws>=max(100000,count*200):raise RuntimeError('proposal cap reached; no implicit relaxation')
        for s,f in operator.batch():
            draws+=1;p=properties(s)
            if not prior_pass(s,p):fail['structural_prior']+=1;continue
            if s in excluded:fail['known_exact']+=1;continue
            if s in seen:fail['duplicate']+=1;continue
            seen.add(s);kept.append((s,f));families[f]+=1
            if len(kept)==count:break
    return kept,dict(raw_draws=draws,accepted=len(kept),rejections=dict(fail),family_counts=dict(families),batch_size=2048)

def bin_index(s):return min((len(s)-12)//5,3)
def select_top(library,reference,training,assets,regime):
    import joblib
    families=['markov','grammar'] if regime=='mixture' else [regime]
    capacity=1024//(4*len(families));counts=Counter();short=[];novelty_rejected=0
    for s,f in library:
        key=(f,bin_index(s))
        if counts[key]>=capacity:continue
        match=process.extractOne(s,reference,scorer=Indel.normalized_similarity)
        if match[1]>.8:novelty_rejected+=1;continue
        short.append(dict(sequence=s,family=f,length_bin=key[1],max_reference_ratio=match[1],nearest_reference=match[0]))
        counts[key]+=1
        if len(short)==1024:break
    if len(short)!=1024:raise RuntimeError(f'insufficient stratified novelty shortlist: {len(short)}')
    sequences=[r['sequence'] for r in short]
    embeddings=esm_embeddings(sequences,assets/'esm2')
    model=joblib.load(assets/'hemo-model.joblib');model.n_jobs=1
    risk=model.predict_proba(embeddings)[:,list(model.classes_).index(1)]
    train=np.load(assets/'hemo-train-embeddings.npy',allow_pickle=False)
    norm=lambda x:x/np.maximum(np.linalg.norm(x,axis=1,keepdims=True),1e-12)
    cosine=norm(embeddings)@norm(train).T
    for i,row in enumerate(short):
        row['hemo_probability']=round(float(risk[i]),8)
        row['max_train_cosine']=round(float(cosine[i].max()),8)
        row['max_train_indel']=process.extractOne(row['sequence'],training,scorer=Indel.normalized_similarity)[1]
    targets={(f,b):(25 if len(families)==1 else ([13,12,13,12] if f=='markov' else [12,13,12,13])[b]) for f in families for b in range(4)}
    selected=[];used=Counter();remaining=list(range(len(short)));nearest=np.zeros(len(short))
    for rank in range(1,101):
        valid=[i for i in remaining if used[(short[i]['family'],short[i]['length_bin'])]<targets[(short[i]['family'],short[i]['length_bin'])] and nearest[i]<=.8]
        if not valid:raise RuntimeError('diverse ranking quota infeasible; no fallback')
        index=min(valid,key=lambda i:(round(short[i]['hemo_probability']+.25*nearest[i],8),short[i]['hemo_probability'],short[i]['sequence']))
        row=dict(short[index],rank=rank,selection_utility=round(short[index]['hemo_probability']+.25*nearest[index],8),max_previously_selected_ratio=float(nearest[index]))
        selected.append(row);remaining.remove(index);used[(row['family'],row['length_bin'])]+=1
        for i in remaining:nearest[i]=max(nearest[i],Indel.normalized_similarity(short[i]['sequence'],row['sequence']))
    return selected,short,dict(novelty_shortlist_rejections=novelty_rejected,shortlist_size=len(short),quota_counts={f'{k[0]}:{k[1]}':v for k,v in used.items()})

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2]);args=parser.parse_args()
    base=args.root.resolve();assets=base/'assets';config=json.loads((base/'config.json').read_text())
    identity=json.loads((assets/'identity.json').read_text())
    for rel,digest in identity['files'].items():
        if sha(assets/rel)!=digest:raise ValueError('runtime asset mismatch: '+rel)
    model=json.loads((assets/'generator.json').read_text());reference=read_fasta(assets/'antibacterial.fasta')
    training=(assets/'hemo-train-sequences.txt').read_text().splitlines()
    excluded=set(reference)|set((assets/'known-assay-sequences.txt').read_text().splitlines())
    started=time.perf_counter()
    library,stats=generate_library(model,config['regime'],config['seed'],50000,excluded)
    selected,short,ranking=select_top(library,reference,training,assets,config['regime'])
    output=base/'generate';output.mkdir(exist_ok=True)
    with (output/'library.fasta').open('w',encoding='ascii',newline='\n') as f:
        for i,(s,family) in enumerate(library,1):f.write(f'>run02_{i:05d}_{family}\n{s}\n')
    with (output/'top.fasta').open('w',encoding='ascii',newline='\n') as f:
        for row in selected:f.write(f">rank_{row['rank']:03d}_risk_{row['hemo_probability']:.8f}\n{row['sequence']}\n")
    audit=dict(protocol_id=config['protocol_id'],config=config,proposal=stats,ranking=ranking,
      interpretation='Hemo risk is a research preference, not observed safety or potency. No efficacy claim.',selected=selected,shortlist=short)
    (output/'audit.json').write_text(json.dumps(audit,sort_keys=True,indent=2)+'\n',encoding='utf-8')
    elapsed=time.perf_counter()-started
    print(json.dumps(dict(generated=50000,top=100,wall_seconds=elapsed,raw_draws=stats['raw_draws'],library_sha256=sha(output/'library.fasta'),top_sha256=sha(output/'top.fasta'))))
if __name__=='__main__':main()
