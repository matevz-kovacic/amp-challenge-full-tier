# Consensus AMP library

Full-tier AMP Challenge submission package: ten generative source methods, a fixed family-weighted Borda ranking, MMR top-list assembly and proportional cluster-stratified library selection. Project-owned code, weights and generated artifacts are MIT licensed; third-party data retains its attributed licenses. Antimicrobial efficacy and benchmark performance are unknown.

```sh
git clone https://github.com/matevz-kovacic/amp-challenge-full-tier.git
cd amp-challenge-full-tier
uv sync
uv run generate
```

Run from the repository root. Python 3.12.3 and inference dependencies are pinned. CPU inference uses one numerical thread and fixed per-source seeds; no training, GPU, API key or LLM is required. The output is:

- generate/library.fasta: exactly 50,000 unique eligible peptides.
- generate/top.fasta: exactly 100 library members, ordered by the original common quality rank.

```text
c7fe0d42580856f9531afe8be20dca22e6f275dbe466eba0f3dcffcb239ef0a3  generate/library.fasta
31cb2b9b81c89dd2a261cb5bbf7aeb208924577b1f5e2fd61c9c9b33ebfa3273  generate/top.fasta
```

## Public reference and offline use

The first generation downloads the public organizer reference from the immutable commit URL in reference_manifest.json and verifies its exact byte length and SHA-256. No candidate sequence is transmitted. The reference is not mirrored or relicensed in this repository; upstream terms apply. Subsequent runs use .cache/reference/antibacterial.fasta. For offline use, set AMP_REFERENCE_FASTA to an existing copy with the required hash before generation. The normal default needs no arguments or environment variables. Installation may also need public Python/package downloads.

## What generation actually does

consolidation.py regenerates all ten complete 50,000-record source pools from their code, weights and fixed seeds. It verifies each ordered pool hash and the exact 499,998-identity union. It then recomputes full-population Borda ranks and MMR/stratified selection from the frozen common scores and fitted cluster assignments. Existing output FASTAs can be deleted. The source pools are generated, not copied from cached candidate FASTAs. Frozen predictor scores and cluster assignments are disclosed method assets; predictor inference and model fitting are not repeated.

The program validates count, alphabet, length, uniqueness, top membership, exact reference exclusion and top normalized-indel similarity. scripts/portable_sequence_checks.py is an attributed BSD derivative of the official validator, using the same RapidFuzz Indel operation underlying the original Levenshtein.ratio. Diversity uses the separate unit-cost edit metric. No reverse-equivalence is applied to identity or official novelty.

Optional commands:

```sh
uv run generate --export-source-pools
uv run python selection/replay_selection.py
uv run python selection/test_selection.py
```

The first also exports the ten intermediate libraries. The second replays only ranking and assembly from cached scores. The third runs meaningful rank, tie, missingness, identity and allocation tests. These are one consolidated submission, not ten independent entries. Source-run old shortlist branches are superseded; use the root entry point.

## Evidence and disclosure

Read ABSTRACT.md, METHOD.md, DECISION.md, RUNS_AND_CONSOLIDATION.md, TRAINING_DATA.md, DATA_AND_FILTERS.md, DISCLOSURE.md, LICENSE_SCOPE.md and THIRD_PARTY_NOTICES.md. schemas, raw scores, fixed ranks, source provenance, selected/reserve lists and reverse relationships are under selection/. All required generative weights are under source_methods/. Optional screening checkpoints and their audit environment are under selection/models/.

REVIEWER_DRY_RUN.md records actual full-tier verification and the unchanged output hashes. FULL_REQUIREMENTS.md maps the organizer requirements. KAGGLE_WRITEUP.md is prepared submission text; publishing this repository does not post it to Kaggle. Public read access requires no organizer invitation.
