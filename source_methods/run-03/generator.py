"""Reconstruct the chemistry-restricted statistical candidate under fixed seed 0."""
import argparse
import json
import random
import time
from pathlib import Path

from rapidfuzz.distance import Indel
from sequence_model import Model


def read_reference(path):
    sequences, parts = [], []
    for line in path.read_text(encoding='utf-8').splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if parts:
                sequences.append(''.join(parts))
            parts = []
        else:
            parts.append(line.upper())
    if parts:
        sequences.append(''.join(parts))
    return sequences


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed', type=int, default=0)
    args = parser.parse_args()
    start = time.monotonic()
    root = Path.cwd()
    weights = json.loads((root / 'assets/weights.json').read_text(encoding='utf-8'))
    model = Model(weights)
    references = read_reference(root / 'assets/antibacterial.fasta')
    reference_set = set(references)
    rng, library, seen = random.Random(args.seed), [], set()
    proposals = duplicates = exact_overlaps = 0
    while len(library) < 50000:
        proposals += 1
        if proposals > 500000:
            raise RuntimeError('Proposal bound exceeded without 50000 valid unique sequences')
        seq = model.sample(rng)
        if seq in seen:
            duplicates += 1
            continue
        if seq in reference_set:
            exact_overlaps += 1
            continue
        library.append(seq)
        seen.add(seq)
    means, sds = weights['typicality_mean'], weights['typicality_sd']
    ranked = []
    for i, seq in enumerate(library):
        features = model.features(seq)
        score = sum(weight * abs((x - mean) / sd) for weight, x, mean, sd in zip([1, 0.5, 0.5], features, means, sds))
        ranked.append((score, i, seq))
    top, similarity_rejects, internal_rejects = [], 0, 0
    for score, i, seq in sorted(ranked):
        if any(Indel.normalized_similarity(seq, ref) > 0.8 for ref in references):
            similarity_rejects += 1
            continue
        if any(Indel.normalized_similarity(seq, existing) > 0.8 for existing in top):
            internal_rejects += 1
            continue
        top.append(seq)
        if len(top) == 100:
            break
    if len(top) != 100:
        raise RuntimeError('Insufficient top-list candidates after disclosed constraints')
    out = root / 'generate'
    out.mkdir(exist_ok=True)
    ids = {seq: i for i, seq in enumerate(library, 1)}
    for name, sequences in [('library', library), ('top', top)]:
        payload = ''.join(f'>run03_{ids[seq]:05d}\n{seq}\n' for seq in sequences)
        (out / (name + '.fasta')).write_bytes(payload.encode('ascii'))
    elapsed = time.monotonic() - start
    print(json.dumps({'seed': args.seed, 'model_order': model.order, 'proposals': proposals,
                      'duplicates_rejected': duplicates, 'exact_reference_rejected': exact_overlaps,
                      'top_reference_rejected': similarity_rejects, 'top_internal_rejected': internal_rejects,
                      'library_count': len(library), 'top_count': len(top), 'wall_seconds': elapsed,
                      'proposals_per_second': proposals / elapsed,
                      'score_scope': 'Uncalibrated training-distribution typicality, not MIC or HC50'}))


if __name__ == '__main__':
    main()
