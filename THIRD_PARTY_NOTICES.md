# Third-party notices

Project-owned components are released under MIT; original source-specific MIT notices remain. The NumPy GRU implementation retains its TensorFlow-equation Apache-2.0 attribution. Model weights are fitted project artifacts, not third-party foundation-model checkpoints.

Default inference dependencies: Python (PSF), NumPy (BSD-3-Clause with its bundled notices), RapidFuzz (MIT), threadpoolctl (BSD-3-Clause). Build tooling: setuptools (MIT); uv (MIT OR Apache-2.0). Optional checkpoint-audit packages scikit-learn, SciPy and joblib retain their BSD-family notices. Dependency binaries are obtained through the pinned package environment.

scripts/verify_submission.py is the untouched organizer BSD-3-Clause validator. scripts/portable_sequence_checks.py is a BSD-3-Clause derivative that changes only its similarity import/calls to MIT-licensed RapidFuzz. The untouched validator can be run in a separate verification environment requiring Levenshtein (GPL-2.0-or-later); this optional tool is not a default generation dependency. Its license text is retained under licenses/.

DRAMP and HemoPI2 datasets/subsets retain CC BY 4.0 with authors, DOI, source versions, hashes and preprocessing in DATA_AND_FILTERS.md and TRAINING_DATA.md. The historical table is attributed to the final CC BY 4.0 article. The synthetic teacher retains the starter-kit MIT notice. The organizer reference is fetched from its pinned upstream repository and is not redistributed in this Git history. Third-party data is not relabeled MIT.

Full notice texts are retained under licenses/ and source_methods/. LICENSE_SCOPE.md states the scope of the root grant.
