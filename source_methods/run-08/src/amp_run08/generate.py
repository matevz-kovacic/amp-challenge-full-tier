"""Reconstruct a deterministic library by sampling learned probability distributions."""
import bisect
import hashlib
from importlib.resources import files
import json
from pathlib import Path
import random

ALPHABET='ACDEFGHIKLMNPQRSTVWY'

def digest(data):
    return hashlib.sha256(data).hexdigest()

def cumulative(weights):
    total=0.; result=[]
    for weight in weights:
        if weight<0: raise ValueError('Negative probability weight')
        total+=weight; result.append(total)
    if not total>0: raise ValueError('Zero probability mass')
    return result

def draw(rng,cdf):
    return bisect.bisect_right(cdf,rng.random()*cdf[-1])

def sample_library(model,seed,count,excluded=()):
    rng=random.Random(seed); excluded=set(excluded)
    lengths=sorted(map(int,model['length_counts']))
    length_cdf=cumulative([model['length_counts'].get(str(n),model['length_counts'].get(n)) for n in lengths])
    kind=model['kind']
    def emission_cdf(weights):
        return cumulative([weight if aa!='C' else 0. for aa,weight in zip(ALPHABET,weights)])
    if kind=='hmm':
        starts=cumulative(model['pi']); transitions=[cumulative(row) for row in model['transition']]
        emissions=[emission_cdf(row) for row in model['emission']]
    elif kind=='composition':
        emission=emission_cdf(model['probabilities'])
    elif kind=='context':
        contexts={ctx:emission_cdf(p) for ctx,p in model['probabilities'].items()}
    else: raise ValueError('Unknown model kind')
    result=[]; seen=set(); attempts=0; duplicate_rejections=0; known_rejections=0
    while len(result)<count:
        attempts+=1
        if attempts>max(1000000,count*100): raise RuntimeError('Generation trial bound exceeded')
        length=lengths[draw(rng,length_cdf)]
        if not 8<=length<=50: raise ValueError('Invalid empirical length')
        chars=[]
        if kind=='hmm': state=draw(rng,starts)
        for position in range(length):
            if kind=='hmm':
                c=ALPHABET[draw(rng,emissions[state])]
                state=draw(rng,transitions[state])
            elif kind=='composition': c=ALPHABET[draw(rng,emission)]
            else:
                ctx=('^'*model['order']+''.join(chars))[-model['order']:]
                while ctx not in contexts: ctx=ctx[1:]
                c=ALPHABET[draw(rng,contexts[ctx])]
            chars.append(c)
        seq=''.join(chars)
        if seq in seen: duplicate_rejections+=1; continue
        if digest(seq.encode()) in excluded: known_rejections+=1; continue
        seen.add(seq); result.append(seq)
    return result,{'proposal_count':attempts,'accepted_unique':len(result),'duplicate_rejections':duplicate_rejections,'known_rejections':known_rejections}

def write_outputs(library,selection,out):
    out=Path(out); out.mkdir(parents=True,exist_ok=True)
    (out/'library.fasta').write_text(''.join(f'>run08_{i+1:06d}\n{s}\n' for i,s in enumerate(library)),encoding='utf-8',newline='\n')
    lines=[]
    for rank,item in enumerate(selection['ranked'],1):
        i=item['library_index']; seq=library[i]
        if digest(seq.encode())!=item['sequence_sha256']: raise ValueError('Frozen ranking does not match regenerated library')
        lines.append(f'>rank_{rank:03d}|library=run08_{i+1:06d}\n{seq}\n')
    if len(lines)!=100: raise ValueError('Top selection must contain 100 records')
    (out/'top.fasta').write_text(''.join(lines),encoding='utf-8',newline='\n')

def main():
    assets=files('amp_run08').joinpath('assets')
    config=json.loads(assets.joinpath('config.json').read_text())
    for name,expected in config['asset_sha256'].items():
        if digest(assets.joinpath(name).read_bytes())!=expected: raise ValueError('Asset hash mismatch: '+name)
    model=json.loads(assets.joinpath('model.json').read_text())
    excluded=assets.joinpath('excluded_sequence_sha256.txt').read_text().splitlines()
    selection=json.loads(assets.joinpath('selection.json').read_text())
    library,counts=sample_library(model,config['seed'],50000,excluded)
    write_outputs(library,selection,'generate')
    print(json.dumps({'seed':config['seed'],'model':model['kind'],**counts},sort_keys=True))

if __name__=='__main__': main()
