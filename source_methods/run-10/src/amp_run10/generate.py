"""A generative stochastic grammar with fixed, disclosed selection rules.

All peptide strings are sampled at runtime. No output sequence list is stored as
model input. A default-stream reference-screening exclusion certificate is an
offline filter asset, not an activity model. Only the sealed default is certified.
"""
import hashlib
import json
import math
import random
from collections import Counter, defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
HYDRO=set('AVILMFWY')
AA='ACDEFGHIKLMNPQRSTVWY'
FAMILIES=['iid','helix','turn','proline']
def digest(s): return hashlib.sha256(s.encode()).hexdigest()
def identity(obj): return digest(json.dumps(obj,sort_keys=True,separators=(',',':')))
def weighted(rng,letters,weights): return rng.choices(letters,weights=weights,k=1)[0]

def propose(rng,family):
    n=rng.choice(list(range(16,33))); phase=rng.randrange(18)
    split=rng.randrange(n//2-2,n//2+3); second_phase=rng.randrange(18)
    seq=[]
    for i in range(n):
        if family=='iid':
            a=weighted(rng,'ADEFGHIKLMNPQRSTVWY',[10,2,2,4,5,1,7,12,12,1,4,2,4,11,5,4,7,2,3])
        elif family in ['helix','turn']:
            if family=='turn' and split<=i<split+2: a=weighted(rng,'GSP',[5,3,2])
            else:
                position=(5*i+(second_phase if family=='turn' and i>=split+2 else phase))%18
                hydrophobic_probability=.78 if position<8 else .18
                a=weighted(rng,'AVILMFWY',[14,12,14,22,2,6,2,4]) if rng.random()<hydrophobic_probability else weighted(rng,'KRQNSTDEHG',[24,20,8,6,8,6,2,2,1,3])
        elif family=='proline':
            a=weighted(rng,'PRKQNSTGAFWILYV',[22,18,12,5,4,5,4,4,6,4,2,4,4,3,3])
        else: raise ValueError(family)
        seq.append(a)
    return ''.join(seq)

def descriptors(s):
    c=Counter(s); n=len(s); q=c['K']+c['R']-c['D']-c['E']; h=sum(c[a] for a in HYDRO)/n
    entropy=-sum(v/n*math.log2(v/n) for v in c.values()); run=hydro_run=0; max_repeat=max_hydro=0; prev=''
    for a in s:
        run=run+1 if a==prev else 1; prev=a; max_repeat=max(max_repeat,run)
        hydro_run=hydro_run+1 if a in HYDRO else 0; max_hydro=max(max_hydro,hydro_run)
    return {'length':n,'charge_proxy':q,'hydrophobic_fraction':h,'entropy_bits':entropy,'max_repeat':max_repeat,'max_hydrophobic_run':max_hydro,'proline_fraction':c['P']/n,'acidic_fraction':(c['D']+c['E'])/n}

def plausible(s,family,top=False):
    d=descriptors(s)
    if not (2<=d['charge_proxy']<=9 and d['max_repeat']<=3 and d['max_hydrophobic_run']<=4 and d['entropy_bits']>=2.65 and d['acidic_fraction']<=.15): return False
    if family=='proline':
        if not (.15<=d['hydrophobic_fraction']<=.45 and .12<=d['proline_fraction']<=.35): return False
    elif not .30<=d['hydrophobic_fraction']<=.60: return False
    if top and (d['length']>30 or not 3<=d['charge_proxy']<=8): return False
    return True

def ratio(a,b):
    # Bit-parallel longest-common-subsequence implementation of exact normalized
    # indel similarity. Cross-checked against pinned Levenshtein in research.
    masks={}
    for i,c in enumerate(b): masks[c]=masks.get(c,0)|(1<<i)
    row=0
    for c in a:
        x=masks.get(c,0)|row; row=x & ~(x-((row<<1)|1))
    return 2*row.bit_count()/(len(a)+len(b))

def generate(config,filters=None):
    filters=filters or {'excluded_library':[],'excluded_top':[]}
    if filters.get('configuration_sha256',identity(config))!=identity(config): raise ValueError('filter certificate belongs to another configuration')
    banned=set(filters['excluded_library']); banned_top=set(filters['excluded_top']); rng=random.Random(config['seed']); sequences=[]; families=[]; seen=set(); proposals=0; counters=defaultdict(int)
    while len(sequences)<config['library_size']:
        family=config['families'][len(sequences)%len(config['families'])]
        s=propose(rng,family); proposals+=1; counters[family+'_proposed']+=1
        if proposals>config['library_size']*100: raise RuntimeError('proposal budget exhausted')
        if s in seen or digest(s) in banned or not plausible(s,family): continue
        seen.add(s); sequences.append(s); families.append(family); counters[family+'_accepted']+=1
    # Per-family queues preserve seeded generation order. Round-robin traversal
    # makes both the first 50 and all 100 cover the configured families.
    queues={family:[i for i,(s,f) in enumerate(zip(sequences,families)) if f==family and plausible(s,f,top=True) and digest(s) not in banned_top] for family in config['families']}
    positions={f:0 for f in queues}; chosen=[]; ranked_families=[]
    while len(chosen)<config['top_size']:
        family=config['families'][len(chosen)%len(config['families'])]
        queue=queues[family]; selected=None
        while positions[family]<len(queue):
            index=queue[positions[family]]; positions[family]+=1; s=sequences[index]
            if all(ratio(s,sequences[j])<=config['internal_similarity_max'] for j in chosen): selected=index; break
        if selected is None: raise RuntimeError('top selection exhausted family '+family)
        chosen.append(selected); ranked_families.append(family)
    return sequences,[sequences[i] for i in chosen],{'proposals':proposals,'family_counts':dict(counters),'library_families':families,'top_families':ranked_families,'top_library_indices':chosen}

def write_outputs(library,top):
    out=Path('generate'); out.mkdir(exist_ok=True)
    index={s:i+1 for i,s in enumerate(library)}
    for name,seqs in [('library.fasta',library),('top.fasta',top)]:
        data=''.join(f'>run10_{index[s]:05d}\n{s}\n' for s in seqs).encode('ascii')
        temporary=out/(name+'.tmp'); temporary.write_bytes(data); temporary.replace(out/name)

def main():
    config=json.loads((ROOT/'config.json').read_text()); filters=json.loads((ROOT/'filter_certificate.json').read_text())
    expected=digest(Path(__file__).read_text())
    if filters['generator_text_sha256']!=expected: raise ValueError('generator text differs from certified default')
    library,top,stats=generate(config,filters); write_outputs(library,top)
    print(json.dumps({'library':len(library),'top':len(top),'proposals':stats['proposals'],'seed':config['seed'],'scope':'de novo generated; biological activity unmeasured'}))
if __name__=='__main__': main()
