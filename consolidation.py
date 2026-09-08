"""Default reviewer entry point: ten generative source pools, then consolidation."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[key]='1'
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import argparse,gzip,hashlib,importlib.util,json,time
from replay_sources import generate_one,read,sha,sequence_hash
from reference_assets import ensure_reference

ROOT=Path(__file__).resolve().parent

def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    result=importlib.util.module_from_spec(spec);sys.modules[name]=result;spec.loader.exec_module(result)
    return result

def verify_outputs(output):
    official=module('official_submission_checks',ROOT/'scripts/portable_sequence_checks.py')
    library=official._verify_sequences(output/'library.fasta')
    official._verify_top(output/'top.fasta',library,100)
    reference=set(official._read_fasta(ROOT/'source_methods/run-01/assets/antibacterial.fasta')[1])
    official._verify_no_overlap(library,reference)
    official._veritfy_max_simularity(set(official._read_fasta(output/'top.fasta')[1]),reference)
    return {'library_unique_count':len(library),'top_unique_count':100,'top_is_library_subset':True,
            'official_sequence_and_reference_checks':'PASS'}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--export-source-pools',action='store_true',help='Also write regenerated source libraries for inspection; default off.')
    args=parser.parse_args()
    manifest=read(ROOT/'replay_manifest.json');start=time.perf_counter();output=Path.cwd()/'generate'
    reference_report=ensure_reference(ROOT)
    for rel,expected in manifest['runtime_file_hashes'].items():
        if sha(ROOT/rel)!=expected:raise ValueError('Frozen runtime asset changed: '+rel)
    print('[1/3] Regenerating all ten source libraries from code, weights and fixed seeds.',flush=True)
    unique=set();records=0;source_reports=[]
    for run,info in manifest['source_pools'].items():
        ts=time.perf_counter();strings=generate_one(run,ROOT/'source_methods');actual=sequence_hash(strings)
        if len(strings)!=info['count'] or actual!=info['ordered_sequence_sha256']:
            raise ValueError('Regenerated source pool differs: '+run)
        records+=len(strings);unique.update(strings)
        source_reports.append({'run':run,'count':len(strings),'ordered_sequence_sha256':actual,'wall_seconds':time.perf_counter()-ts})
        print(f'{run}: {len(strings):,} generated; ordered sequence SHA-256 {actual}; PASS',flush=True)
        if args.export_source_pools:
            folder=output/'source_libraries';folder.mkdir(parents=True,exist_ok=True)
            (folder/(run+'.fasta')).write_bytes(''.join(f'>{run}_{i:06d}\n{s}\n' for i,s in enumerate(strings,1)).encode('ascii'))
    sequences=sorted(unique)
    if records!=500000 or len(sequences)!=499998:raise ValueError('Source comparison population changed')
    print(f'Pooled {records:,} source records into {len(sequences):,} distinct N-to-C identities.',flush=True)
    print('[2/3] Recomputing common-panel Borda ranks, MMR and proportional cluster selection from frozen scores.',flush=True)
    selection=module('frozen_consolidation_selection',ROOT/'selection/replay_selection.py')
    selection.main(generated_sequences=sequences)
    output.mkdir(exist_ok=True)
    for kind,expected in manifest['output_sha256'].items():
        payload=(ROOT/'selection/replayed'/(kind+'.fasta')).read_bytes()
        if hashlib.sha256(payload).hexdigest()!=expected:raise ValueError('Consolidated artifact differs: '+kind)
        (output/(kind+'.fasta')).write_bytes(payload)
    print('[3/3] Checking the final files with the disclosed BSD-derived sequence checks using equivalent normalized-indel similarity.',flush=True)
    result=verify_outputs(output)
    report={'entry_point':'uv run generate','source_pools':source_reports,'source_records':records,
            'unique_comparison_population':len(sequences),'common_score_table_population_matches_regenerated_pools':True,
            'ranking_and_assembly_recomputed':True,'predictor_inference_repeated':False,'model_training_performed':False,
            'clustering_refitted':False,'read_saved_source_candidate_output_fastas':False,
            'source_pool_exports_requested':args.export_source_pools,'output_sha256':manifest['output_sha256'],
            'reference_acquisition':reference_report,'validation':result,'wall_seconds':time.perf_counter()-start}
    (output/'generation_report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    for kind,digest in manifest['output_sha256'].items():print(f'{kind}.fasta SHA-256 {digest}',flush=True)
    print('PASS: 50,000 unique eligible library members; ranked top 100; fixed submitted bytes reproduced.',flush=True)

if __name__=='__main__':main()
