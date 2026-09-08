"""Reconstruct all source libraries from frozen generative code and assets.

Source top-list ranking is superseded by frozen consolidation selection. For
run03/04 monolithic entry points, an AST adapter executes the unchanged prefix
through the library-generation loop and returns before source shortlist scoring.
No fitting function or saved candidate FASTA is used to produce a library.
"""
import os,sys
sys.dont_write_bytecode=True
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import ast,contextlib,hashlib,importlib.util,json,time
import numpy as np
from threadpoolctl import threadpool_limits

def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def sequence_hash(strings):return hashlib.sha256(('\n'.join(strings)+'\n').encode('ascii')).hexdigest()
def fasta(p):
    strings=[];parts=[]
    for line in Path(p).read_text(encoding='utf-8-sig').splitlines():
        if line.startswith('>'):
            if parts:strings.append(''.join(parts).upper());parts=[]
        elif line.strip():parts.append(line.strip())
    if parts:strings.append(''.join(parts).upper())
    return strings
def module(name,p):
    spec=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
@contextlib.contextmanager
def context(p):
    cwd=Path.cwd();argv=sys.argv.copy();path=sys.path.copy()
    try:
        os.chdir(p);sys.argv=['generate'];sys.path[:0]=[str(p),str(p/'src')];yield
    finally:os.chdir(cwd);sys.argv=argv;sys.path[:]=path
def main_prefix(m,p):
    tree=ast.parse(p.read_text(encoding='utf-8-sig'));fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='main')
    stop=next(i for i,n in enumerate(fn.body) if isinstance(n,ast.While) and 'len(library)' in ast.unparse(n.test))
    fn.name='consolidation_library_prefix';fn.body=fn.body[:stop+1]+[ast.Return(value=ast.Name(id='library',ctx=ast.Load()))]
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),str(p),'exec'),m.__dict__)
    return m.consolidation_library_prefix()
def generate_one(run,root):
    p=Path(root)/run
    with context(p),threadpool_limits(limits=1):
        if run=='run-01':
            m=module('source01',p/'generator.py');lib,_,_=m.generate(read(p/'assets/model.json'),fasta(p/'assets/antibacterial.fasta'),read(p/'config.json'));return [r['sequence'] for r in lib]
        if run=='run-02':
            import run02amp.core as m
            assets=p/'assets';cfg=read(p/'config.json');excluded=set(fasta(assets/'antibacterial.fasta'))|set((assets/'known-assay-sequences.txt').read_text().splitlines())
            lib,_=m.generate_library(read(assets/'generator.json'),cfg['regime'],cfg['seed'],50000,excluded);return [s for s,_ in lib]
        if run=='run-03':
            module('sequence_model',p/'sequence_model.py');m=module('source03',p/'generator.py');return main_prefix(m,p/'generator.py')
        if run=='run-04':
            src=p/'src/amp_run04/generate.py';m=module('source04',src);return main_prefix(m,src)
        if run=='run-05':
            m=module('source05',p/'src/amp_run05/generate.py');return m.generate(read(p/'src/amp_run05/config.json'))[0]
        if run=='run-06':
            m=module('source06',p/'src/amp_run06/__init__.py');return m.generate(read(p/'config.json'),p/'assets/antibacterial.fasta')[0]
        if run=='run-07':
            from amp_run07.gru_model import GRUModel
            model=GRUModel(dict(np.load(p/'assets/gru-weights.npz',allow_pickle=False)))
            return model.library(50000,seed=0,excluded=set(fasta(p/'assets/antibacterial.fasta')),no_cysteine=True,max_proposals=1000000)[0]
        if run=='run-08':
            m=module('source08',p/'src/amp_run08/generate.py');assets=p/'src/amp_run08/assets';cfg=read(assets/'config.json')
            return m.sample_library(read(assets/'model.json'),cfg['seed'],50000,(assets/'excluded_sequence_sha256.txt').read_text().splitlines())[0]
        if run=='run-09':
            m=module('source09',p/'src/amp_run09/generate.py');assets=p/'src/amp_run09/assets';return m.propose(read(assets/'model.json'),set((assets/'reference-sha256.txt').read_text().splitlines()))[0]
        if run=='run-10':
            m=module('source10',p/'src/amp_run10/generate.py');return m.generate(read(p/'config.json'),read(p/'filter_certificate.json'))[0]
    raise ValueError(run)
def main():
    base=Path(__file__).resolve().parent
    # Packaged mode: source_methods is alongside this file. Root research mode uses a prepared package.
    if not (base/'source_methods').exists():base=base/'submission'
    manifest=read(base/'replay_manifest.json')
    for rel,h in manifest['runtime_file_hashes'].items():
        if sha(base/rel)!=h:raise ValueError('Frozen runtime asset changed: '+rel)
    pools={};evidence=[];start=time.perf_counter()
    for run,info in manifest['source_pools'].items():
        ts=time.perf_counter();strings=generate_one(run,base/'source_methods');actual=sequence_hash(strings)
        if len(strings)!=info['count'] or actual!=info['ordered_sequence_sha256']:raise ValueError('Regenerated source pool differs: '+run)
        pools[run]=strings;evidence.append({'run':run,'count':len(strings),'ordered_sequence_sha256':actual,'wall_seconds':time.perf_counter()-ts})
        print(run,'reconstructed',len(strings),'hash matches',flush=True)
    selection=read(base/'frozen_selection.json');out=Path.cwd()/'generate';out.mkdir(exist_ok=True)
    for kind in ['library','top']:
        lines=[];seen=set()
        for row in selection[kind]:
            s=pools[row['source_run']][row['source_index_zero_based']]
            if hashlib.sha256(s.encode()).hexdigest()!=row['sequence_sha256']:raise ValueError('Selection identity changed')
            if s in seen:raise ValueError('Duplicate output')
            seen.add(s);lines.append(f">amp_cons_q{row['quality_rank']:06d}\n{s}\n")
        b=''.join(lines).encode('ascii')
        if hashlib.sha256(b).hexdigest()!=manifest['output_sha256'][kind]:raise ValueError('Reconstructed artifact differs from frozen selection')
        (out/(kind+'.fasta')).write_bytes(b)
    report={'source_pools':evidence,'all_pool_hashes_match':True,'read_saved_source_candidate_output_fastas':False,'reference_or_teacher_filter_assets_read':True,'fitted_any_model':False,'output_sha256':manifest['output_sha256'],'wall_seconds':time.perf_counter()-start}
    (out/'generative_replay.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('Reconstructed consolidation from all ten generative source pools.',flush=True)
if __name__=='__main__':main()
