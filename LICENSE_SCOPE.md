# License scope

SPDX-License-Identifier: MIT

The root MIT license applies to project-owned source code, inference adapters, documentation, trained model parameters, generated source candidate pools, frozen prediction/rank/selection artifacts and final peptide collections, unless a more specific retained notice applies. Existing MIT, BSD and Apache notices and copyright holders are preserved. This license is granted only for rights held by the submission contributors; it makes no claim over third-party rights or peptide patents.

DRAMP and HemoPI2 source data and their attributed subsets retain CC BY 4.0. The historical experimental table retains its published CC BY 4.0 basis. The public synthetic teacher file retains its starter-kit MIT notice. See DATA_AND_FILTERS.md, TRAINING_DATA.md and THIRD_PARTY_NOTICES.md.

The organizer reference is not mirrored or relicensed in this public repository. Inference obtains the exact public organizer file described by reference_manifest.json and verifies its SHA-256. Cached files retain upstream terms. No undisclosed non-public biological training data is used. All project-created training/filter selections and artifacts needed to understand the method are disclosed in this release.

The default runtime uses permissively licensed dependencies. The untouched organizer validator is BSD-3-Clause code; its separate optional verification environment uses the GPL-2.0-or-later Levenshtein package, with that dependency's notice retained. GPL tooling is not imported by the default generation program.
