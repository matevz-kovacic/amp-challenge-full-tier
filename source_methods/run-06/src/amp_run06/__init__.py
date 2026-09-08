"""Seeded stochastic peptide grammars. Desirability is not a biological predictor."""
import argparse
import json
import math
import random
import time
from collections import Counter
from pathlib import Path

from rapidfuzz import process
from rapidfuzz.distance import Indel

ALPHABET='ACDEFGHIKLMNPQRSTVWY'
HYDRO=set('AFILMVWY')
POSITIVE=set('KR')
GROUPS={'hydrophobic':'AAIIILLLVVFMWY','positive':'KKKKRR','polar':'SSSTTNNQQGG','acidic':'DE'}

def descriptors(s):
    n=len(s);counts=Counter(s)
    charge=sum(counts[a] for a in 'KR')-sum(counts[a] for a in 'DE')
    hydro=sum(a in HYDRO for a in s)/n
    real=sum((a in HYDRO)*math.cos(i*math.pi*100/180) for i,a in enumerate(s))
    imag=sum((a in HYDRO)*math.sin(i*math.pi*100/180) for i,a in enumerate(s))
    moment=math.hypot(real,imag)/n
    longest=run=repeat=longest_repeat=0;last=None
    for a in s:
        run=run+1 if a in HYDRO else 0;longest=max(longest,run)
        repeat=repeat+1 if a==last else 1;longest_repeat=max(longest_repeat,repeat);last=a
    return {'length':n,'charge':charge,'charge_fraction':charge/n,'hydrophobic_fraction':hydro,
            'binary_hydrophobic_moment':moment,'longest_hydrophobic_run':longest,
            'maximum_positive_5window':max(sum(x in POSITIVE for x in s[i:i+5]) for i in range(n-4)),
            'longest_repeat':longest_repeat,'entropy':-sum((v/n)*math.log2(v/n) for v in counts.values())}

def desirability(d):
    # A bounded, hand-chosen physicochemical preference; no assay fitting or calibrated units.
    return -((d['hydrophobic_fraction']-.45)/.14)**2-((d['charge_fraction']-.20)/.13)**2-((d['length']-21)/10)**2-((d['binary_hydrophobic_moment']-.20)/.18)**2

def allowed(d):
    return (2<=d['charge']<=8 and .30<=d['hydrophobic_fraction']<=.62 and d['longest_hydrophobic_run']<=4
            and d['longest_repeat']<=2 and d['maximum_positive_5window']<=3 and d['entropy']>=2.5)

def propose(rng,mode):
    n=rng.randrange(12,33);phase=rng.random()*2*math.pi
    hinge=rng.randrange(max(5,n//3),min(n-4,2*n//3)+1)
    chars=[]
    for i in range(n):
        if mode=='asymmetric' and i==hinge:
            chars.append(rng.choice('PGS'));continue
        if mode=='iid':
            weights=(.44,.24,.27,.05)
        else:
            offset=phase+(1.3 if mode=='asymmetric' and i>hinge else 0)
            face=math.cos(i*math.pi*100/180+offset)
            h=.44+.27*face
            weights=(h,.24-.08*face,.27-.19*face,.05)
        group=rng.choices(tuple(GROUPS),weights=weights,k=1)[0]
        chars.append(rng.choice(GROUPS[group]))
    return ''.join(chars)

def sample_context(rng,model):
    lengths=[int(x) for x in model['length_counts']]
    n=rng.choices(lengths,weights=list(model['length_counts'].values()),k=1)[0]
    context='^^';s=[]
    for _ in range(n):
        key=context
        while key not in model['distributions']:key=key[1:]
        a=rng.choices(model['alphabet'],weights=model['distributions'][key],k=1)[0]
        s.append(a);context=(context+a)[-2:]
    return ''.join(s)

def read_reference(path):
    seqs=[];parts=[]
    for line in path.read_text(encoding='utf-8').splitlines():
        line=line.strip()
        if line.startswith('>'):
            if parts:seqs.append(''.join(parts));parts=[]
        elif line:parts.append(line.upper())
    if parts:seqs.append(''.join(parts))
    return sorted(set(seqs))

def generate(config,reference):
    import hashlib
    if hashlib.sha256(reference.read_bytes()).hexdigest()!=config['reference_sha256']:
        raise ValueError('Pinned reference hash mismatch')
    refs=read_reference(reference);refset=set(refs);rng=random.Random(config['seed'])
    modes=config['modes'];model=json.loads(Path('assets/context_model.json').read_text(encoding='utf-8'));refset.update(read_reference(Path('assets/teacher_library.fasta')));library=[];seen=set();records={};trials=0;rejections=Counter()
    while len(library)<config['library_size']:
        if trials>=config['max_proposals']:
            partial=Path('generate/partial_library.fasta');partial.parent.mkdir(exist_ok=True)
            partial.write_bytes(''.join(f'>partial_{i+1}\n{x}\n' for i,x in enumerate(library)).encode())
            print(json.dumps({'partial_library_count':len(library),'proposal_trials':trials,'rejections':dict(rejections)}))
            raise RuntimeError('Proposal cap exhausted; incomplete output is explicitly named partial_library.fasta')
        mode=modes[len(library)%len(modes)];s=sample_context(rng,model) if mode=='distilled' else propose(rng,mode);trials+=1
        if s in seen or s in refset or 'C' in s or not 8<=len(s)<=50:rejections['duplicate_or_exact_reference']+=1;continue
        d=descriptors(s)
        if not allowed(d):rejections['design_constraints']+=1;continue
        library.append(s);seen.add(s);records[s]={'mode':mode,'descriptors':d,'score':desirability(d)}
    ranked={m:sorted((s for s in library if records[s]['mode']==m),key=lambda s:(-records[s]['score'],s)) for m in dict.fromkeys(modes)}
    cursors={m:0 for m in ranked};top=[];evaluated=0
    # Round-robin coverage across regimes. Within each regime, descending fixed desirability.
    while len(top)<config['top_size']:
        m=list(ranked)[len(top)%len(ranked)]
        if cursors[m]>=len(ranked[m]):raise RuntimeError('Insufficient novel top candidates in regime '+m)
        s=ranked[m][cursors[m]];cursors[m]+=1;evaluated+=1
        if any(Indel.normalized_similarity(s,t)>config['top_internal_max_similarity'] for t in top):
            rejections['top_internal_similarity']+=1;continue
        nearest=process.extractOne(s,refs,scorer=Indel.normalized_similarity,score_cutoff=config['top_reference_max_similarity'])
        if nearest is not None and nearest[1]>config['top_reference_max_similarity']:
            rejections['top_reference_similarity']+=1;continue
        top.append(s)
    stats={'proposal_trials':trials,'accepted_library':len(library),'top_evaluated':evaluated,
           'rejections':dict(rejections),'modes':Counter(records[s]['mode'] for s in library),
           'top_records':[{'rank':i+1,'sequence':s,**records[s]} for i,s in enumerate(top)]}
    return library,top,stats

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--config',default='config.json');args=parser.parse_args()
    config=json.loads(Path(args.config).read_text(encoding='utf-8'));start=time.monotonic()
    library,top,stats=generate(config,Path('assets/antibacterial.fasta'))
    directory=Path('generate');directory.mkdir(exist_ok=True)
    # Stable output serialization; headers and ranking order are part of the reconstruction contract.
    identifiers={s:i+1 for i,s in enumerate(library)}
    for name,seqs in [('library',library),('top',top)]:
        data=''.join(f'>run06_{identifiers[s]:05d}\n{s}\n' for s in seqs).encode('ascii')
        temporary=directory/(name+'.fasta.tmp');temporary.write_bytes(data);temporary.replace(directory/(name+'.fasta'))
    stats['wall_seconds']=time.monotonic()-start
    print(json.dumps(stats,sort_keys=True))

if __name__=='__main__':main()
