"""CPU generative models with explicit, unvalidated centrality/diversity priors."""
import argparse
import hashlib
import json
import math
import random
import statistics
import time
from collections import Counter,defaultdict
from pathlib import Path
from rapidfuzz import process
from rapidfuzz.distance import Indel

ALPHABET='ADEFGHIKLMNPQRSTVWY' # cysteine-free primary design hypothesis
REFERENCE_SHA256='cbbeac64ba95746d87961e8ad9dd0849ae8058d15a300b2e7f6990730ca521e9'

def digest(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def descriptors(s):
    n=len(s);c=Counter(s)
    return [n,(c['K']+c['R']-c['D']-c['E'])/n,sum(c[a] for a in 'AILMFWVY')/n,(c['G']+c['P'])/n]

def fit(sequences):
    sequences=sorted(set(sequences)); counts=defaultdict(Counter)
    for s in sequences:
        prefix='^^'
        for a in s:
            for n in range(3):counts[prefix[-n:] if n else ''][a]+=1
            prefix+=a
    ds=list(zip(*(descriptors(s) for s in sequences)))
    med=[statistics.median(x) for x in ds]
    scale=[max(.05 if j else 4.,statistics.median(abs(v-med[j]) for v in x)*1.4826) for j,x in enumerate(ds)]
    return {'schema':'run04-generative-weights-v1','sequences':sequences,
            'lengths':dict(Counter(map(len,sequences))),
            'contexts':{k:[v[a] for a in ALPHABET] for k,v in sorted(counts.items())},
            'descriptor_medians':med,'descriptor_scales':scale,'backoff_mass':8.0}

class Model:
    def __init__(self,weights):
        self.w=weights;self.cache={}
        globalcounts=weights['contexts'][''];total=sum(globalcounts)+len(ALPHABET)
        self.prior=[(x+1)/total for x in globalcounts]

    def probabilities(self,context):
        if context in self.cache:return self.cache[context]
        if not context:return self.prior
        previous=self.probabilities(context[1:]);counts=self.w['contexts'].get(context)
        if counts is None:return previous
        tau=self.w['backoff_mass'];den=sum(counts)+tau
        result=[(c+tau*p)/den for c,p in zip(counts,previous)]
        self.cache[context]=result;return result

    def draw(self,rng,probabilities):
        u=rng.random();running=0.
        for a,p in zip(ALPHABET,probabilities):
            running+=p
            if u<running:return a
        return ALPHABET[-1]

    def sample(self,rng,method):
        templates=self.w['sequences']
        if method=='recombine':
            a,b=rng.choice(templates),rng.choice(templates)
            cut1=rng.randint(max(1,len(a)//3),max(1,2*len(a)//3))
            cut2=rng.randint(max(1,len(b)//3),max(1,2*len(b)//3))
            s=a[:cut1]+b[cut2:]
            return ''.join(self.draw(rng,self.prior) if rng.random()<.25 else c for c in s)
        length=len(rng.choice(templates));prefix='^^'
        for _ in range(length):
            context=prefix[-2:] if method=='markov' else ''
            prefix+=self.draw(rng,self.probabilities(context))
        return prefix[2:]

    def log_likelihood(self,s,order=2):
        total=0.;prefix='^^'
        for a in s:
            context=prefix[-order:] if order else ''
            total+=math.log(self.probabilities(context)[ALPHABET.index(a)])
            prefix+=a
        return total/len(s)

    def score(self,s):
        # Versioned custom search score only. It is not MIC, HC50 or organizer grading.
        z=[(v-m)/d for v,m,d in zip(descriptors(s),self.w['descriptor_medians'],self.w['descriptor_scales'])]
        return self.log_likelihood(s)-.15*sum(min(4,abs(v))**2 for v in z)

def plausible(s):
    c=Counter(s);n=len(s);charge=c['K']+c['R']-c['D']-c['E']
    h=sum(c[a] for a in 'AILMFWVY')/n
    entropy=-sum(v/n*math.log2(v/n) for v in c.values())
    return (12<=n<=35 and 2<=charge<=10 and .25<=h<=.65 and entropy>=2.
            and all(a*4 not in s for a in ALPHABET))

def fasta_sequences(path):
    seqs=[];parts=[]
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            if parts:seqs.append(''.join(parts).upper());parts=[]
        elif line.strip():parts.append(line.strip())
    if parts:seqs.append(''.join(parts).upper())
    return seqs

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seed',type=int,default=0)
    p.add_argument('--n-sequences',type=int,default=50000)
    p.add_argument('--top-k',type=int,default=100)
    a=p.parse_args();start=time.monotonic();root=Path.cwd()
    config=json.loads((root/'config.json').read_text())
    if digest(root/'assets/antibacterial.fasta')!=REFERENCE_SHA256:raise ValueError('Wrong reference')
    w=json.loads((root/'assets/weights.json').read_text());model=Model(w)
    reference=fasta_sequences(root/'assets/antibacterial.fasta')
    known=set(reference)|set(w['sequences']);rng=random.Random(a.seed)
    seen=set();library=[];counts=Counter()
    while len(library)<a.n_sequences:
        if counts['proposals']>a.n_sequences*100:raise RuntimeError('Proposal budget exhausted')
        counts['proposals']+=1;s=model.sample(rng,config['method'])
        if not 8<=len(s)<=50:counts['invalid_length']+=1;continue
        if s in known:counts['exact_known']+=1;continue
        if s in seen:counts['duplicate']+=1;continue
        library.append(s);seen.add(s)
    generation_seconds=time.monotonic()-start
    ranking_namespace='custom/run04-centrality-v1'
    eligible=[i for i,s in enumerate(library) if plausible(s)]
    if config.get('ranking_mode')=='esm':
        from amp_run04.esm_rank import pool_scores
        def centrality(s):
            z=[(v-m)/d for v,m,d in zip(descriptors(s),w['descriptor_medians'],w['descriptor_scales'])]
            return -sum(min(4,abs(v))**2 for v in z)
        pool=sorted(eligible,key=lambda i:(-centrality(library[i]),library[i]))[:512]
        if len(pool)!=512:raise RuntimeError('Fewer than 512 eligible reranking candidates')
        neural=pool_scores([library[i] for i in pool],root/'assets/esm2-8m')
        scores={i:value+.15*centrality(library[i]) for i,value in zip(pool,neural)}
        order=sorted(pool,key=lambda i:(-scores[i],library[i]))
        counts['neural_scored']=len(pool)
        ranking_namespace='custom/run04-neural-selection-v1'
    else:
        order=sorted(eligible,key=lambda i:(-model.score(library[i]),library[i]))
    selected=[];selection_indices=[]
    for i in order:
        counts['top_considered']+=1;s=library[i]
        if any(Indel.normalized_similarity(s,t)>.70 for t in selected):
            counts['top_internal_similarity']+=1;continue
        nearest=process.extractOne(s,reference,scorer=Indel.normalized_similarity,score_cutoff=.8)
        if nearest is not None and nearest[1]>.8:
            counts['top_reference_similarity']+=1;continue
        nearest=process.extractOne(s,w['sequences'],scorer=Indel.normalized_similarity,score_cutoff=.8)
        if nearest is not None and nearest[1]>.8:
            counts['top_training_similarity']+=1;continue
        selected.append(s);selection_indices.append(i)
        if len(selected)==a.top_k:break
    if len(selected)!=a.top_k:raise RuntimeError(f'Only {len(selected)} eligible diverse top sequences')
    out=root/'generate';out.mkdir(exist_ok=True)
    for name,indices in [('library.fasta',range(len(library))),('top.fasta',selection_indices)]:
        data=''.join(f'>run04_{config["method"]}_{i+1:06d}\n{library[i]}\n' for i in indices).encode()
        tmp=out/(name+'.tmp');tmp.write_bytes(data);tmp.replace(out/name)
    print(json.dumps({'method':config['method'],'seed':a.seed,'library_count':len(library),'top_count':len(selected),
                      'counts':counts,'generation_seconds':generation_seconds,'wall_seconds':time.monotonic()-start,
                      'proposal_sequences_per_second':counts['proposals']/generation_seconds,
                      'library_sha256':digest(out/'library.fasta'),'top_sha256':digest(out/'top.fasta'),
                      'ranking_namespace':ranking_namespace,'config_sha256':digest(root/'config.json'),
                      'weights_sha256':digest(root/'assets/weights.json')},sort_keys=True))

if __name__=='__main__':main()
