"""Run-01 empirical and physical generative hypotheses. Original code: MIT.

Prior likelihood and physical heuristics are not validated potency scores.
"""
import argparse, collections, functools, hashlib, json, math, random, re, time
from pathlib import Path
from rapidfuzz.distance import Indel
from rapidfuzz import process

AA='ACDEFGHIKLMNPQRSTVWY';VOCAB=AA.replace('C','');HYDRO=set('AILMFWVY')

def read_fasta(path):
    headers=[];seqs=[];header=None;parts=[]
    for line in Path(path).read_text().splitlines():
        line=line.strip()
        if not line:continue
        if line.startswith('>'):
            if header is not None:headers.append(header);seqs.append(''.join(parts))
            header=line[1:];parts=[]
        else:parts.append(line.upper())
    if header is not None:headers.append(header);seqs.append(''.join(parts))
    return headers,seqs

def properties(s):
    n=len(s);charge=s.count('K')+s.count('R')-s.count('D')-s.count('E');hydro=sum(a in HYDRO for a in s)/n
    freq=collections.Counter(s)
    return {'length':n,'charge':charge,'charge_density':charge/n,'hydrophobic_fraction':hydro,
      'entropy':-sum((v/n)*math.log(v/n) for v in freq.values()),
      'binary_hydrophobic_moment':abs(sum((a in HYDRO)*complex(math.cos(i*math.pi/1.8),math.sin(i*math.pi/1.8)) for i,a in enumerate(s)))/n}

def invalid_reason(s):
    if not 12<=len(s)<=40:return 'length'
    if not set(s)<=set(VOCAB):return 'alphabet_or_cysteine'
    p=properties(s)
    if not 2<=p['charge']<=12 or p['charge_density']>.4:return 'charge'
    if not .25<=p['hydrophobic_fraction']<=.65:return 'hydrophobic_fraction'
    if re.search(r'(.)\1\1\1',s):return 'homopolymer'
    return None

def train_prior(records):
    parents=[];weights=[];counts={'':{a:.25 for a in VOCAB}};lengths={}
    for record in records:
        s=record['sequence'];rate=record['activity_fraction']['16'];w=1+4*rate
        if not 12<=len(s)<=40 or not set(s)<=set(AA):continue
        parents.append(s);weights.append(w);lengths[len(s)]=lengths.get(len(s),0)+w
        for i,a in enumerate(s):
            if a not in VOCAB:continue
            counts[''][a]+=w
            for k in [1,2]:
                if i>=k:
                    context=s[i-k:i];bucket=counts.setdefault(context,{})
                    bucket[a]=bucket.get(a,0)+w
    assert parents
    return {'schema':'run01-empirical-prior-v1','vocabulary':VOCAB,'parents':parents,'parent_weights':weights,
      'length_values':sorted(lengths),'length_weights':[lengths[n] for n in sorted(lengths)],'counts':counts,
      'smoothing_order1':15.0,'smoothing_order2':5.0,'fitness_weight':'1+4*observed eleven-strain activity fraction at16uM',
      'potency_validation':'NOT_ESTABLISHED: generic and direct regressors failed preregistered criteria.'}

class Operator:
    def __init__(self,model,seed):
        self.model=model;self.rng=random.Random(seed);self.seed=seed
        total=sum(model['counts'][''].values());self.unigram=tuple(model['counts'][''].get(a,0)/total for a in VOCAB)
    @functools.lru_cache(maxsize=1000)
    def probabilities(self,context):
        if not context:return self.unigram
        base=self.probabilities(context[1:]);bucket=self.model['counts'].get(context,{})
        alpha=self.model['smoothing_order1' if len(context)==1 else 'smoothing_order2'];total=sum(bucket.values())+alpha
        return tuple((bucket.get(a,0)+alpha*p)/total for a,p in zip(VOCAB,base))
    def propose(self,kind):
        if kind=='markov':
            n=self.rng.choices(self.model['length_values'],self.model['length_weights'])[0];s=''
            for _ in range(n):s+=self.rng.choices(VOCAB,self.probabilities(s[-2:]))[0]
            return s
        if kind=='mutation':
            parent=self.rng.choices(self.model['parents'],self.model['parent_weights'])[0];s=list(parent)
            count=max(math.ceil(.25*len(s)),int(self.rng.uniform(.25,.45)*len(s)))
            positions=set(self.rng.sample(range(len(s)),count)) | {i for i,a in enumerate(s) if a=='C'}
            for i in sorted(positions):
                choices=[a for a in VOCAB if a!=s[i]];w=[self.unigram[VOCAB.index(a)] for a in choices]
                s[i]=self.rng.choices(choices,w)[0]
            return ''.join(s)
        if kind=='helix':
            n=self.rng.randint(16,32);phase=self.rng.random()*2*math.pi;s=''
            for i in range(n):
                hydrophobic=math.cos(i*math.pi/1.8+phase)>.1
                if self.rng.random()<.25:hydrophobic=not hydrophobic
                letters='AILVFMFWY' if hydrophobic else 'KKRRQNSTAG'
                s+=self.rng.choice(letters)
            return s
        raise ValueError(kind)
    def prior_score(self,s):
        likelihood=sum(math.log(self.probabilities(s[max(0,i-2):i])[VOCAB.index(a)]) for i,a in enumerate(s))/len(s)
        p=properties(s)
        penalty=.015*abs(len(s)-24)+.025*abs(p['charge']-6)+.4*abs(p['hydrophobic_fraction']-.45)
        return round(likelihood-penalty,12)

def maximum_ratio(sequence,reference):
    return float(process.extractOne(sequence,reference,scorer=Indel.normalized_similarity)[1])

def generate(model,reference,config):
    seed=config.get('seed',0);n=config.get('library_size',50000);top_n=config.get('top_size',100)
    kinds=config.get('kinds',['markov','mutation']);operator=Operator(model,seed);refset=set(reference)
    seen=set();items=[];rejections=collections.Counter();proposals=collections.Counter();accepted=collections.Counter()
    quota={kind:n//len(kinds)+(i<n%len(kinds)) for i,kind in enumerate(kinds)}
    for kind in kinds:
        while accepted[kind]<quota[kind]:
            if proposals[kind]>=config.get('max_proposals_per_operator',1000000):raise RuntimeError('Proposal limit reached: '+kind)
            s=operator.propose(kind);proposals[kind]+=1;reason=invalid_reason(s)
            if reason:rejections[kind+':'+reason]+=1;continue
            if s in seen:rejections[kind+':duplicate']+=1;continue
            if s in refset:rejections[kind+':reference_exact']+=1;continue
            seen.add(s);accepted[kind]+=1;items.append({'sequence':s,'operator':kind,'score':operator.prior_score(s)})
    ranked={kind:sorted([r for r in items if r['operator']==kind],key=lambda r:(-r['score'],r['sequence'])) for kind in kinds}
    cursors={k:0 for k in kinds};chosen=[];checks=0;reject_top=collections.Counter()
    for rank in range(top_n):
        kind=kinds[rank%len(kinds)]
        while cursors[kind]<len(ranked[kind]):
            item=ranked[kind][cursors[kind]];cursors[kind]+=1;s=item['sequence'];checks+=1
            if any(Indel.normalized_similarity(s,r['sequence'])>.7 for r in chosen):reject_top['top_redundancy']+=1;continue
            maximum=maximum_ratio(s,reference)
            if maximum>.8:reject_top['reference_similarity']+=1;continue
            chosen.append(dict(item,maximum_reference_ratio=maximum));break
        else:raise RuntimeError('Insufficient novel balanced top candidates: '+kind)
    stats={'seed':seed,'proposals':dict(proposals),'accepted':dict(accepted),'rejections':dict(rejections),
      'top_checks':checks,'top_rejections':dict(reject_top),'top_max_reference_ratio':max(r['maximum_reference_ratio'] for r in chosen),
      'ranking':'Alternating operator balance; descending prior-density score within each; novelty and pairwise diversity constraints.',
      'scientific_status':'Unvalidated generative hypothesis; no efficacy or selectivity score.'}
    return items,chosen,stats

def main():
    p=argparse.ArgumentParser();p.add_argument('--seed',type=int);a=p.parse_args()
    root=Path.cwd();model=json.loads((root/'assets/model.json').read_text());config=json.loads((root/'config.json').read_text())
    if a.seed is not None:config['seed']=a.seed
    reference_path=root/'assets/antibacterial.fasta'
    if hashlib.sha256(reference_path.read_bytes()).hexdigest()!='cbbeac64ba95746d87961e8ad9dd0849ae8058d15a300b2e7f6990730ca521e9':raise ValueError('Pinned reference missing or corrupt')
    _,reference=read_fasta(reference_path);start=time.monotonic();items,top,stats=generate(model,reference,config)
    out=root/'generate';out.mkdir(exist_ok=True)
    identifiers={item['sequence']:f'run01_{item["operator"]}_{i+1:05d}' for i,item in enumerate(items)}
    for name,rows in [('library',items),('top',top)]:
        (out/(name+'.fasta')).write_bytes(''.join(f'>{identifiers[r["sequence"]]}\n{r["sequence"]}\n' for r in rows).encode('ascii'))
    (out/'generation.json').write_text(json.dumps(stats,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(json.dumps({'library_count':len(items),'top_count':len(top),'wall_seconds':time.monotonic()-start,'proposals':stats['proposals']}))
if __name__=='__main__':main()
