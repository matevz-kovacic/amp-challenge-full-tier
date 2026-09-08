"""Original data-free sequence sampler and documented design descriptors.

No activity, toxicity or synthesis probability is inferred by these descriptors.
The default outputs are computed from the configuration and seed, not read back
from stored sequences. All optional research parameters have fixed defaults.
"""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import random
import time

ALPHABET="ACDEFGHIKLMNPQRSTVWY"
HYDRO=set("AILMFWVY")
CONFIG=Path(__file__).with_name("config.json")

def descriptors(sequence):
    n=len(sequence);c=Counter(sequence)
    charge=c["K"]+c["R"]-c["D"]-c["E"]
    hydro=sum(c[a] for a in HYDRO)/n
    # Explicit binary amphipathic-pattern statistic, not an empirical energy scale.
    x=sum((int(a in HYDRO)-hydro)*math.cos(math.radians(100*i)) for i,a in enumerate(sequence))
    y=sum((int(a in HYDRO)-hydro)*math.sin(math.radians(100*i)) for i,a in enumerate(sequence))
    moment=math.hypot(x,y)/n
    entropy=-sum((v/n)*math.log2(v/n) for v in c.values())
    run=1;max_run=1
    for a,b in zip(sequence,sequence[1:]):run=run+1 if a==b else 1;max_run=max(max_run,run)
    return {"length":n,"formal_sidechain_charge":charge,"charge_density":charge/n,"hydrophobic_fraction":hydro,"binary_helical_segregation":moment,"residue_entropy_bits":entropy,"proline_fraction":c["P"]/n,"glycine_fraction":c["G"]/n,"aromatic_fraction":sum(c[a] for a in "FWY")/n,"max_identical_run":max_run}

def design_score(sequence):
    d=descriptors(sequence)
    # A declared plausibility heuristic. Higher values are not measured potency.
    score=(-abs(d["charge_density"]-0.25)-abs(d["hydrophobic_fraction"]-0.45)
           +0.5*min(d["binary_helical_segregation"],0.25)
           -0.04*max(0,d["max_identical_run"]-2)-0.05*sequence.count("P")
           -0.02*abs(d["length"]-24))
    return round(score,12)

def propose(rng,config):
    n=rng.randint(config["length_min"],config["length_max"])
    family=config["family"]
    if family=="composition":
        return "".join(rng.choices(ALPHABET,weights=[config["composition_weights"][a] for a in ALPHABET],k=n))
    if family in ["amphipathic","scrambled_amphipathic"]:
        phase=rng.randrange(360);parts=[]
        for i in range(n):
            angle=(100*i+phase)%360
            if angle<80 or angle>280:
                chars="ALIVFWKRTSGNQ";weights=[22,23,14,14,4,2,6,4,3,3,1,2,2]
            else:
                chars="KRQNSTAGLIVF";weights=[28,18,10,6,8,6,6,6,6,2,2,2]
            parts.append(rng.choices(chars,weights=weights,k=1)[0])
        if family=="scrambled_amphipathic":rng.shuffle(parts)
        return "".join(parts)
    if family=="block_amphiphile":
        # Orthogonal, nonperiodic charge/hydrophobe arrangement control.
        left=rng.randint(3,6);right=rng.randint(3,6)
        core=n-left-right
        return ("".join(rng.choices("KRQST",weights=[4,3,1,1,1],k=left))+
                "".join(rng.choices("ALIVFGS",weights=[3,3,2,2,1,1,1],k=core))+
                "".join(rng.choices("KRQST",weights=[4,3,1,1,1],k=right)))
    raise ValueError(f"Unknown family {family}")

def generate(config):
    rng=random.Random(config["seed"]);library=[];seen=set();proposals=0
    while len(library)<config["count"]:
        proposals+=1
        if proposals>config["max_proposals"]:raise RuntimeError("Explicit proposal cap reached")
        s=propose(rng,config)
        if not 8<=len(s)<=50 or set(s)-set(ALPHABET):raise ValueError("Configuration generated an invalid string")
        if s not in seen:seen.add(s);library.append(s)
    if config.get("selected_indices") is not None:
        indices=config["selected_indices"]
        if len(indices)!=config["top_count"] or len(set(indices))!=len(indices):raise ValueError("Invalid fixed selection indices")
        if any(type(i) is not int or i<0 or i>=len(library) for i in indices):raise ValueError("Selection outside library")
    elif config["rank_method"]=="hash_control":
        indices=sorted(range(len(library)),key=lambda i:(hashlib.sha256(library[i].encode()).hexdigest(),library[i]))[:config["top_count"]]
    elif config["rank_method"]=="design_heuristic":
        indices=sorted(range(len(library)),key=lambda i:(-design_score(library[i]),library[i]))[:config["top_count"]]
    else:raise ValueError("Unknown ranking method")
    return library,indices,proposals

def write_outputs(directory,library,indices):
    directory.mkdir(parents=True,exist_ok=True)
    (directory/"library.fasta").write_bytes("".join(f">run05_{i+1:05d}\n{s}\n" for i,s in enumerate(library)).encode("ascii"))
    (directory/"top.fasta").write_bytes("".join(f">rank_{r+1:03d}_run05_{i+1:05d}\n{library[i]}\n" for r,i in enumerate(indices)).encode("ascii"))

def main():
    p=argparse.ArgumentParser();p.add_argument("--config",type=Path,default=CONFIG);p.add_argument("--output",type=Path,default=Path("generate"));p.add_argument("--seed",type=int,default=None)
    args=p.parse_args();config=json.loads(args.config.read_text(encoding="utf-8"))
    if args.seed is not None:config["seed"]=args.seed
    start=time.monotonic();library,indices,proposals=generate(config);write_outputs(args.output,library,indices)
    elapsed=time.monotonic()-start
    print(json.dumps({"family":config["family"],"seed":config["seed"],"proposals":proposals,"unique_sequences":len(library),"top":len(indices),"wall_seconds":elapsed,"proposals_per_second":proposals/elapsed,"library_sha256":hashlib.sha256((args.output/"library.fasta").read_bytes()).hexdigest(),"top_sha256":hashlib.sha256((args.output/"top.fasta").read_bytes()).hexdigest()}))
if __name__=="__main__":main()
