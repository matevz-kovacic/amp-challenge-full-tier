# Full-tier requirements mapping

Source: https://github.com/szczurek-lab/amp-challenge-2027, checked 2026-09-08 UTC. The user reports accepting the Kaggle conditions. This release targets the full submission tier; organizer acceptance, scientific qualification and eventual co-authorship decisions are not claimed.

| Requirement | Artifact/evidence |
|---|---|
| Abstract | ABSTRACT.md |
| 50,000 library and ranked top100 | generate/library.fasta and generate/top.fasta; exact immutable hashes in README.md |
| Ranking and set assembly documentation | METHOD.md, DECISION.md and selection/; individual Q remains separate from MMR insertion order |
| Training data, external databases, filters and intervention | TRAINING_DATA.md, DATA_AND_FILTERS.md, DISCLOSURE.md and attributed data/weights |
| Public repository with inference code and weights | Public release URL in README.md; source_methods/ plus consolidation.py and replay_sources.py |
| Permissive OSI-approved project license | MIT in LICENSE and pyproject.toml; scope and third-party exceptions in LICENSE_SCOPE.md |
| uv and defined Python version | pyproject.toml, uv.lock and .python-version |
| Default generate entry point and fixed random seeds | uv run generate; all source seeds are fixed and disclosed; optional source export defaults off |
| Identical repeated output | REVIEWER_DRY_RUN.md and validation/; compares against the originally selected FASTA bytes |
| Full training-data disclosure and non-public-data release | All biological data sources are public; exact available/processed rows and provenance are disclosed. Project-created weights, candidate pools and audit artifacts are included under the stated license. |

The full library must have no exact matches in the organizer reference. The top100 additionally reject Levenshtein.ratio strictly greater than 0.8; equality passes. The known official wording conflict over broad library similarity and first50 versus top100 wet-lab sampling remains documented in METHOD.md and DECISION.md. Both shortlist sizes have reported metrics. Hidden competition weights are unknown.

The organizer's sample directory layout illustrates an entry point, not a required neural architecture. This repository implements the documented uv entry-point contract and retains the unchanged validator for external checking. Its separate GPL-dependent verification environment is not part of default inference. No newly fitted model or candidate experiment was added for this release.
