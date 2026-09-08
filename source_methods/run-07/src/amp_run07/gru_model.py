"""NumPy reset-after GRU inference, independently checked against TensorFlow 2.2.1.

Equations/reference: TensorFlow v2.2.1 keras/layers/recurrent.py GRUCell,
Apache-2.0; no TensorFlow runtime or third-party weights required.
"""
from __future__ import annotations
import math
from collections import Counter
import numpy as np
ALPHABET='ACDEFGHIKLMNPQRSTVWY'

class GRUModel:
    def __init__(self,weights,profile_model=None):
        self.w={k:np.asarray(v,dtype=np.float32) for k,v in weights.items()}
        self.width=self.w['recurrent_kernel'].shape[0]
        self.profile_model=profile_model

    def step(self,tokens,state):
        x=self.w['embedding'][tokens]@self.w['kernel']+self.w['bias'][0]
        v=state@self.w['recurrent_kernel']+self.w['bias'][1]
        h=self.width
        update=1/(1+np.exp(-(x[:,:h]+v[:,:h])))
        reset=1/(1+np.exp(-(x[:,h:2*h]+v[:,h:2*h])))
        candidate=np.tanh(x[:,2*h:]+reset*v[:,2*h:])
        state=update*state+(1-update)*candidate
        return state@self.w['output_kernel']+self.w['output_bias'],state

    def teacher_logits(self,sequences):
        width=max(map(len,sequences))+1
        tokens=np.zeros((len(sequences),width),dtype=np.int32);tokens[:,0]=21
        for i,s in enumerate(sequences):
            tokens[i,1:len(s)+1]=[ALPHABET.index(a)+1 for a in s]
        state=np.zeros((len(sequences),self.width),dtype=np.float32);out=[]
        for j in range(width):
            logits,state=self.step(tokens[:,j],state);out.append(logits)
        return np.stack(out,axis=1)

    def nll_batch(self,sequences):
        logits=self.teacher_logits(sequences).astype(np.float64)
        logits[:,:8,0]=-np.inf
        m=np.max(logits,axis=2,keepdims=True)
        logp=logits-m-np.log(np.sum(np.exp(logits-m),axis=2,keepdims=True))
        result=[]
        for i,s in enumerate(sequences):
            ids=[ALPHABET.index(a)+1 for a in s]
            ll=sum(logp[i,j,t] for j,t in enumerate(ids))
            if len(s)<50: ll+=logp[i,len(s),0]
            result.append(-float(ll)/(len(s)*math.log(2)))
        return result

    def nll_per_residue(self,seq): return self.nll_batch([seq])[0]

    def component(self,seq):
        if self.profile_model is None: return 0
        return self.profile_model.component(seq)

    def propose_batch(self,rng,batch_size=256):
        state=np.zeros((batch_size,self.width),dtype=np.float32)
        tokens=np.full(batch_size,21,dtype=np.int32)
        alive=np.ones(batch_size,dtype=bool);sequences=['']*batch_size
        for j in range(50):
            logits,state=self.step(tokens,state)
            # Float64 sampling probabilities are part of this fixed implementation.
            logits=logits.astype(np.float64)
            if j<8: logits[:,0]=-np.inf
            p=np.exp(logits-np.max(logits,axis=1,keepdims=True));p/=p.sum(axis=1,keepdims=True)
            u=rng.random(batch_size)
            tokens=np.minimum(np.sum(u[:,None]>np.cumsum(p,axis=1),axis=1),20).astype(np.int32)
            for i,t in enumerate(tokens):
                if alive[i] and t: sequences[i]+=ALPHABET[t-1]
            alive &= tokens!=0
        return sequences

    def library(self,n,seed=0,excluded=None,no_cysteine=True,max_proposals=1000000):
        rng=np.random.default_rng(seed);excluded=set(excluded or []);seen=set();out=[];rejected=Counter();proposed=0
        while proposed<max_proposals:
            for s in self.propose_batch(rng,256):
                proposed+=1
                if no_cysteine and 'C' in s: rejected['cysteine']+=1
                elif s in seen: rejected['duplicate']+=1
                elif s in excluded: rejected['reference_exact']+=1
                else:
                    seen.add(s);out.append(s)
                    if len(out)==n: return out,{'seed':seed,'batch_size':256,'proposals':proposed,'rejected':dict(rejected)}
                if proposed>=max_proposals: break
        raise RuntimeError('Proposal cap reached')
