# Ten source methods and the consolidation program

One repository supplies everything needed for one consolidated AMP Challenge submission. Each `source_methods/run-01` through `run-10` directory contains its necessary inference implementation, frozen weights/configuration and filters. Original research-framework files, agent instructions, prompts, logs, unused predictors and old per-run shortlist experiments are not needed for this runtime.

| Source | Packaged generative method | Fixed seed | Regenerated records |
|---|---|---:|---:|
| run-01 | Empirical Markov/mutation model | 0 | 50,000 |
| run-02 | HemoPI2-derived count mixture | 20260908 | 50,000 |
| run-03 | Chemistry-restricted statistical sequence model | 0 | 50,000 |
| run-04 | DRAMP empirical sequence model | 0 | 50,000 |
| run-05 | Procedural composition grammar | 50317 | 50,000 |
| run-06 | Synthetic-teacher-derived context model | 20260908 | 50,000 |
| run-07 | Frozen GRU with NumPy inference | 0 | 50,000 |
| run-08 | DRAMP hidden Markov model | 20260908 | 50,000 |
| run-09 | Empirical/Markov mixture | 0 | 50,000 |
| run-10 | Procedural sequence grammar | 271828 | 50,000 |

`consolidation.py` coordinates all ten source methods via the inspected adapters in `replay_sources.py`. Their full union must exactly match the 499,998 N-to-C strings defining the saved common score table. The program then invokes `selection/replay_selection.py` to recompute per-evaluator average ranks, family-weighted Borda Q, MMR insertion order and proportional cluster-stratified membership. Individual quality ranking remains separate from diversity-aware set assembly.

The common raw scores and fitted cluster labels are frozen method assets. Reproduction uses them consistently for every pooled identity; it does not repeat predictor inference, model training or clustering fitting. Saved scaler/centers, clustering parameters, predictor weights, feature code, source provenance and expanded-rank replay are retained for review. Frozen selection indices are retained as provenance, but the default command now calculates selection from the common scores and frozen cluster assignments.

The final output contains 50,000 library members and a quality-ordered top 100. Only those two FASTAs are the submitted peptide collections. `--export-source-pools` is an optional inspection facility for the ten intermediate libraries, not ten independent Kaggle entries. DATA_AND_FILTERS.md and METHOD.md describe data, chemistry, eligibility and selection limitations.

Full-tier v1.1 release: source01/source03 use the existing RapidFuzz normalized-indel operation directly; the coordinator uses an attributed portable validator and fetches the exact pinned reference. Original input/output sequence hashes, fixed scores/ranks and cluster assignments remain unchanged. See LICENSE_SCOPE.md and TRAINING_DATA.md.
