# Full-tier reviewer verification

Fresh-checkout verification of this release is being recorded before publication. Commands are uv sync followed by uv run generate twice. The submitted library and top100 hashes remain those in README.md. Default inference must not import the GPL Levenshtein package.

The untouched organizer validator is scripts/verify_submission.py. Its pinned SHA-256 is 3f2eb1bd61200abfccf07d86e9c226d57f3d12abcf25715af1d90f41531942cf. It can be run in a separate environment containing Levenshtein 0.27.3, for example after obtaining the reference:

```sh
uv run --isolated --with Levenshtein==0.27.3 python scripts/verify_submission.py https://github.com/matevz-kovacic/amp-challenge-full-tier --dir ../organizer-check --antibacterial-fasta .cache/reference/antibacterial.fasta
```

That audit tool clones, installs, runs generation twice, validates the files/reference, and checks repeated bytes. It is optional tooling with its own dependency license; it does not change the default inference environment.
