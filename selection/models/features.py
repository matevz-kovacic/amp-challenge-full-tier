"""Unchanged descriptor function from the fitted HemoPI2 screening implementation."""
import numpy as np
ALPHABET="ACDEFGHIKLMNPQRSTVWY"
INDEX={c:i for i,c in enumerate(ALPHABET)}
KD=np.array([1.8,2.5,-3.5,-3.5,2.8,-.4,-3.2,4.5,-3.9,3.8,1.9,-3.5,-1.6,-3.5,-4.5,-.8,-.7,4.2,-.9,-1.3])

def features(sequences):
    rows=[]
    for s in sequences:
        aac=np.bincount([INDEX[c] for c in s],minlength=20)/len(s)
        di=np.zeros(400)
        for a,b in zip(s,s[1:]): di[INDEX[a]*20+INDEX[b]]+=1/(len(s)-1)
        extras=[len(s)/50,float(aac@KD),(s.count('K')+s.count('R')-s.count('D')-s.count('E'))/len(s),sum(s.count(c) for c in 'FWY')/len(s)]
        rows.append(np.concatenate([aac,extras,di]))
    return np.array(rows)
