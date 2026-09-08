# Full-tier reviewer verification: PASS

The runtime and environment at commit `a18834ddcf2e3fa0d849cff10f28bdd6ad1a4c58` passed fresh Windows and Ubuntu checkouts. Later commits add documentation and regression-test evidence; runtime assets remain identical.

```sh
uv sync
uv run generate
uv run generate
```

| Check | Observed result |
|---|---|
| Windows, untouched full organizer CLI | PASS; fresh clone/venv, plain uv sync, two generation calls, original sequence/reference checks and byte equality; 115.4 seconds total |
| Ubuntu WSL2, commands above | PASS; generation 50.1 and 47.5 seconds |
| Public reference acquisition | Both first runs obtained the exact upstream file by HTTPS, with no local override; repeated runs used the verified cache |
| Linux output FASTAs removed before generation | Both recreated with the exact original SHA-256 hashes |
| Default runtime GPL package check | Levenshtein is absent on both platforms |
| Original official sequence checks | Passed using the untouched validator in separate audit interpreters |
| Ratio replacement | All 3,944,800 top100/reference values match Levenshtein 0.27.3 exactly as float64 |
| Rank/assembly regressions | Six selection tests passed on both platforms |
| Public-release regressions | Three additional Windows tests passed: 0.8 boundary, reversal identity, reference cache and hash failure |

Actual transcripts and reports are in validation/. Local machine paths are replaced for privacy. Dependency wheels were reused from existing caches; the locked versions were not upgraded. A dependency-resolution attempt initially lacked offline index metadata; refreshing public index metadata produced the pinned lock with Levenshtein removed. There were no failed generator attempts or model changes in the full-tier tests.

Every generation regenerated 500,000 source records, verified 499,998 distinct source identities, recomputed all frozen quality ranks and final selection, and reproduced the exact original 50,000/100 FASTA hashes. No predictor inference, training or LLM ran. Runtime dependence on the pinned public reference is explicit; it is not a hidden local input.

The untouched organizer validator is scripts/verify_submission.py. Its pinned SHA-256 is 3f2eb1bd61200abfccf07d86e9c226d57f3d12abcf25715af1d90f41531942cf. It can be run in a separate environment containing Levenshtein 0.27.3, for example after obtaining the reference:

```sh
uv run --isolated --with Levenshtein==0.27.3 python scripts/verify_submission.py https://github.com/matevz-kovacic/amp-challenge-full-tier --dir ../organizer-check --antibacterial-fasta .cache/reference/antibacterial.fasta
```

That audit tool clones, installs, runs generation twice, validates the files/reference, and checks repeated bytes. It is optional tooling with its own dependency license; it does not change the default inference environment.
