"""Standalone integer-weight stochastic peptide generator (stdlib only).

This source is copied into each candidate package. It has no research-framework
dependency, network call, learned biological claim, or stored output library.
"""
import bisect,collections,hashlib,json,random,time
from pathlib import Path

ALPHABET='ACDEFGHIKLMNPQRSTVWY'
HYDROPHOBIC=set('AILMFWVY')
def digest_sequences(sequences):return hashlib.sha256(('\n'.join(sequences)+'\n').encode('ascii')).hexdigest()
def cumulative(weights):
    total=0;out=[]
    for x in weights:total+=x;out.append(total)
    return out
def draw(rng,weights):return bisect.bisect_right(weights,rng.randrange(weights[-1]))
def violation(s):
    n=len(s);charge=sum(s.count(a) for a in 'KR')-sum(s.count(a) for a in 'DE')
    if not 2<=charge<=int(0.4*n):return 'charge_filter'
    frac=sum(a in HYDROPHOBIC for a in s)/n
    if not 0.25<=frac<=0.65:return 'hydrophobic_fraction_filter'
    run=0;same=0;last=None
    for a in s:
        run=run+1 if a in HYDROPHOBIC else 0
        same=same+1 if a==last else 1;last=a
        if run>6:return 'hydrophobic_run_filter'
        if same>3:return 'homopolymer_filter'
    return None
def propose(config,blocked):
    rng=random.Random(config['seed']);seen=set();seqs=[];families=[];rejects=collections.Counter()
    init=cumulative(config['markov']['initial'])
    transitions=[cumulative(row) for row in config['markov']['transitions']]
    length_weights=cumulative(config['markov']['length_weights'])
    helix=[[cumulative(row) for row in phase] for phase in config['helix']['phase_weights']]
    n=config['library_size'];attempted=0
    while len(seqs)<n and attempted<config['maximum_proposals']:
        attempted+=1
        kind=config['mode']
        if kind=='mixture':kind='markov' if rng.randrange(2)==0 else 'helix'
        if kind=='markov':
            length=config['length_min']+draw(rng,length_weights)
            prev=draw(rng,init);chars=[ALPHABET[prev]]
            for _ in range(length-1):prev=draw(rng,transitions[prev]);chars.append(ALPHABET[prev])
        elif kind=='helix':
            length=rng.randrange(config['length_min'],config['length_max']+1)
            phase=rng.randrange(len(helix));reverse=bool(rng.randrange(2))
            chars=[ALPHABET[draw(rng,helix[phase][i])] for i in range(length)]
            if reverse:chars.reverse()
        else:raise ValueError('Unknown generative family')
        s=''.join(chars);bad=violation(s)
        if bad:rejects[bad]+=1;continue
        if s in seen:rejects['duplicate']+=1;continue
        if hashlib.sha256(s.encode('ascii')).hexdigest() in blocked:rejects['exact_reference_overlap']+=1;continue
        seen.add(s);seqs.append(s);families.append(kind)
    if len(seqs)!=n:raise RuntimeError(f'Proposal budget exhausted: {len(seqs)} of {n}; {dict(rejects)}')
    return seqs,families,{'attempted':attempted,'accepted':len(seqs),'rejections':dict(sorted(rejects.items())),'family_counts':dict(collections.Counter(families))}
def write_outputs(sequences,selection,directory):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    if digest_sequences(sequences)!=selection['library_sequence_sha256']:raise ValueError('Library identity differs from selection; rebuild and reverify the whole candidate')
    indices=selection['ordered_library_indices_zero_based']
    if len(indices)!=100 or len(set(indices))!=100 or not all(0<=i<len(sequences) for i in indices):raise ValueError('Invalid selected indices')
    library=''.join(f'>run09_{i+1:05d}\n{s}\n' for i,s in enumerate(sequences))
    top=''.join(f'>rank_{rank:03d}_library_{i+1:05d}\n{sequences[i]}\n' for rank,i in enumerate(indices,1))
    for name,text in [('library.fasta',library),('top.fasta',top)]:
        p=directory/name;temp=p.with_suffix('.fasta.tmp');temp.write_bytes(text.encode('ascii'));temp.replace(p)
def main():
    start=time.perf_counter();assets=Path(__file__).with_name('assets')
    config=json.loads((assets/'model.json').read_text())
    blocked=set((assets/'reference-sha256.txt').read_text().splitlines())
    expected=config['reference_hash_list_sha256']
    if hashlib.sha256((assets/'reference-sha256.txt').read_bytes()).hexdigest()!=expected:raise ValueError('Reference exclusion asset corrupt')
    sequences,families,stats=propose(config,blocked)
    selection=json.loads((assets/'selection.json').read_text())
    write_outputs(sequences,selection,Path.cwd()/'generate')
    print(json.dumps({'generator':config['mode'],'seed':config['seed'],'library_sequence_sha256':digest_sequences(sequences),**stats,'wall_seconds':time.perf_counter()-start}))
if __name__=='__main__':main()
