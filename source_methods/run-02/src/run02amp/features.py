"""Explicit sequence descriptors and CPU-only frozen ESM2 representations."""
import math
import os
import hashlib
from collections import Counter
from pathlib import Path
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[key]='2'
import numpy as np

ALPHABET='ACDEFGHIKLMNPQRSTVWY'
HYDROPHOBIC=set('AILMFWVY')
def validate(sequence):
    if not 8<=len(sequence)<=50 or set(sequence)-set(ALPHABET): raise ValueError('noncanonical or out-of-range peptide')
def descriptors(sequences,dipeptides=False):
    vectors=[]
    for sequence in sequences:
        validate(sequence);n=len(sequence);counts=Counter(sequence)
        fractions=[counts[a]/n for a in ALPHABET]
        h=[float(a in HYDROPHOBIC) for a in sequence]
        run=longest=0
        for v in h:
            run=run+1 if v else 0;longest=max(longest,run)
        moment=lambda theta:abs(sum(v*complex(math.cos(i*theta),math.sin(i*theta)) for i,v in enumerate(h)))/n
        features=fractions+[n/50,(counts['K']+counts['R']-counts['D']-counts['E'])/n,
            sum(h)/n,-sum(p*math.log(p) for p in fractions if p)/math.log(20),longest/n,
            moment(math.radians(100)),moment(math.pi)]
        for terminus in [sequence[:3],sequence[-3:]]:
            features.extend(terminus.count(a)/3 for a in ALPHABET)
        if dipeptides:
            pairs=Counter(zip(sequence,sequence[1:]))
            features.extend(pairs[(a,b)]/(n-1) for a in ALPHABET for b in ALPHABET)
        vectors.append(features)
    return np.asarray(vectors,dtype=np.float64)

def load_esm(asset_directory):
    import torch
    import esm
    torch.set_num_threads(2)
    if torch.get_num_interop_threads()!=1:torch.set_num_interop_threads(1)
    torch.use_deterministic_algorithms(True);torch.manual_seed(0)
    directory=Path(asset_directory)
    expected={'esm2_t6_8M_UR50D.pt':'46f002a9870c9bdecd0ea887acb1f9a38a6b561e8f8bf8a6990b679b9d31b928',
              'esm2_t6_8M_UR50D-contact-regression.pt':'8f7a4557d57713b97ba0e484303007efb7230d25299c0ac47a0a1b12a87bbb9d'}
    for name,digest in expected.items():
        if hashlib.sha256((directory/name).read_bytes()).hexdigest()!=digest:raise ValueError('ESM checkpoint hash mismatch')
    # These are the directly retrieved official FAIR checkpoints. Their hashes
    # are checked and recorded before loading; no user-supplied pickle is used.
    data=torch.load(directory/'esm2_t6_8M_UR50D.pt',map_location='cpu',weights_only=False)
    regression=torch.load(directory/'esm2_t6_8M_UR50D-contact-regression.pt',map_location='cpu',weights_only=False)
    model,alphabet=esm.pretrained.load_model_and_alphabet_core('esm2_t6_8M_UR50D',data,regression)
    model=model.to('cpu').eval()
    return model,alphabet

def esm_embeddings(sequences,asset_directory,batch_size=64):
    import torch
    model,alphabet=load_esm(asset_directory);converter=alphabet.get_batch_converter();result=[]
    for begin in range(0,len(sequences),batch_size):
        batch=sequences[begin:begin+batch_size]
        for sequence in batch:validate(sequence)
        _,_,tokens=converter([(str(i),s) for i,s in enumerate(batch)])
        with torch.inference_mode():
            output=model(tokens.to('cpu'),repr_layers=[6],return_contacts=False)['representations'][6]
        result.extend(output[i,1:len(s)+1].mean(0).numpy() for i,s in enumerate(batch))
        print(f'CPU embeddings {min(begin+batch_size,len(sequences))}/{len(sequences)}',flush=True)
    return np.asarray(result,dtype=np.float32)
