# Complete training-data disclosure

All training sources are public. No proprietary or otherwise non-public biological training dataset was used. New project-created sequence selections, fitted parameters and generated candidate artifacts are disclosed here. No training runs are required for inference. DATA_AND_FILTERS.md contains source authors, DOI, versions, original file hashes and preprocessing.

| Source/component | Exact data disclosure |
|---|---|
| run01 and run09 | source_methods/run-01/assets/study-mic.csv: attributed 46-peptide study table; model.json files contain empirical fitted parameters |
| run02 | source_methods/run-02/assets/generator.json includes training_rows and original HemoPI2 training hash; training_data/hemopi2_cross_val_dataset.csv supplies the public source |
| run03 | source_methods/run-03/assets/training_records.json and dramp-attribution.json contain the selected DRAMP rows and split/provenance |
| run04 | source_methods/run-04/assets/training_rows.json contains the selected DRAMP source rows |
| run05 and run10 | Procedural grammars; no fitted training data |
| run06 | source_methods/run-06/assets/teacher_library.fasta supplies the public synthetic teacher; its source method implements the frozen hash split/filter and context_model.json retains the fitted counts |
| run07 | source_methods/run-07/assets/training.json discloses exact sequences, DRAMP IDs, groups and split assignments; gru-weights.npz contains trained weights |
| run08 | training_data/run08_hmm_training.json discloses all 292 cysteine-free training strings and DRAMP IDs |
| Common hemolysis panel | Both original public HemoPI2 CSVs are in training_data/. hemolysis_model_training.json supplies the exact canonical, deduplicated, non-conflicting development rows and labels. The independent CSV was evaluated, not fitted. selection/models contains frozen checkpoints and feature code. |

DRAMP/HemoPI2 source data and attributed subsets retain CC BY 4.0 (licenses/CC-BY-4.0.txt). The historical table uses the final article's CC BY 4.0 basis, with original attribution retained. Synthetic teacher data retains its starter-kit MIT notice. Source-specific transformation details are in DATA_AND_FILTERS.md and the existing source files. Training_data/manifest.json pins the additional data files.

The organizer's aggregated antibacterial reference is used only for rejection/novelty checks. It is fetched from a fixed public upstream commit and hash-verified, not redistributed or claimed as MIT training data. Every source still uses the exact same reference bytes as in the verified original generation. No candidate is sent in that request.
