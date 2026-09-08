# Training data, external databases and filters

All ten source methods are frozen; no training occurred during consolidation or package preparation. Source code, generative weights and configuration hashes appear in `replay_manifest.json`. No private biological training dataset, private scientific service, or undisclosed pretrained model is required by the final generator.

| Source | Generative data and method | Relevant preprocessing/filtering |
|---|---|---|
| 01 | Small empirical model from 46 Torres et al. historical study peptides | Empirical counts and original design filters; organizer-reference exclusion |
| 02 | HemoPI2-derived sequence count mixture | Canonical length filtering; exclude organizer reference and known assay sequences |
| 03 | DRAMP chemistry-restricted sequence model | 341 rows / 340 strings; grouped split 272/25/43; only training weights; explicit linear/free-N/free-C/L/unmodified annotations |
| 04 | DRAMP empirical sequence model | 293 rows / 292 unique strings, same explicit chemistry restriction plus no cysteine; frozen count/recombination configuration |
| 05 | Original procedural stochastic grammar | No training dataset or pretrained weights; fixed composition/amphipathic configuration |
| 06 | Context model fitted to public HydrAMP synthetic teacher outputs | Fixed hash split and feasibility filters leave 1,561 training strings; exclude teacher/reference exact matches |
| 07 | Newly fitted GRU, exported as NumPy arrays | DRAMP grouped training subset of 2,761 rows; source chemistry support restrictions and no-cysteine generation; epoch-10 weights |
| 08 | DRAMP hidden Markov model | 292 chemistry-curated cysteine-free strings; exact known-sequence digest exclusion |
| 09 | Empirical/Markov mixture based on the 46-peptide historical study | Fixed design grammar, length restrictions and reference digest exclusion |
| 10 | Original stochastic grammar | No learned data; fixed configuration and rejection-digest certificate |

The original optional design filters differ by source and are reproduced only within their originating generator. Consolidation imposes the common hard rules and common panel described in METHOD.md; it does not impose any run's optional design criteria on the entire pooled library. Old source ESM-based shortlist ordering and other historical predictors are superseded and are not required for this submitted inference.

DRAMP source: Ma Tianyue and Heng Zheng (2024), *Peptides data of DRAMP database*, Figshare 27233508 version 2, DOI https://doi.org/10.6084/m9.figshare.27233508.v2, CC BY 4.0. general_amps.xlsx SHA-256: 2196e49100cae4553045873ce640eadb6168555dec5097686842e316ba462f2b. Changes include chemistry/activity/alphabet/length filtering, deduplication and source-specific grouped splitting. Retained records/IDs are packaged for sources 03, 04 and 07. Database annotations are not independent verification of physical chemistry.

HemoPI2 source: Anand Singh Rathore, Nishant Kumar, Shubham Choudhury, Naman Kumar Mehta and Gajendra Raghava, *Prediction of Hemolytic Peptides and their Hemolytic Concentration*, https://doi.org/10.5281/zenodo.14676712, dataset deposit CC BY 4.0. cross_val_dataset.csv SHA-256: 7bdaf3ede499d1eda2712585d2e52d7700f3f138776d1a0a46e2ca88e8152da0; independent_dataset.csv: 500013c2244219762ff3ff4a03401c7419790c83d0ef0c3aeebfdbea426b3eb5. Author GitHub revision 2b67a5c85422b25ae847100ebaa81ad586950928 supplied matching archival files. Source02 generation uses training-derived counts; the source08-derived common screening models were independently fitted with scikit-learn on canonical, deduplicated training records. Original HemoPI2 GPL predictor code and weights are not used. Published strong-hemolysis labels correspond to HC50 <=100 micromolar with heterogeneous or incompletely specified assay/chemistry context. Independent data support diagnostics and exclusion, not fitting.

Historical study: Torres, Chen, Wan, Chatterjee and de la Fuente-Nunez, *Generative latent diffusion language modeling yields anti-infective synthetic peptides*, Cell Biomaterials 1, 100183 (2025), https://doi.org/10.1016/j.celbio.2025.100183. The final publication is recorded as CC BY 4.0; the earlier preprint has different terms and is not the license basis. Table source: szczurek-lab/ampdiffusion-starter-kit commit 1a862af9078e6b55c87d1fa576f3da81851ba94b, experimental/mic.csv, SHA-256 3da5adb7814f00b5f331eedd3dc909c55ee66cc1561e1bf33d5c1753d30eb82d. The selected 46-peptide cohort is not an unbiased AMP population. Historical results are not measurements of these candidates.

Synthetic teacher: szczurek-lab/hydramp-starter-kit commit 7804df862872ccc6d09fe01c41bafbca194cfa31, generate_broad_spectrum/library.fasta, SHA-256 8a717e059d3b6d63cf0ca2a58a83d143a6202f97387b2180d9f7f906f3955c25, original MIT notice retained. These published synthetic strings train source06's count model and serve as an exclusion filter. Original HydrAMP neural weights and experimental labels are not part of the submitted model.

Novelty-only organizer reference: szczurek-lab/amp-challenge-2027 commit 5c8a5d8e2551c8cf572d3d3bfcfe7633b109d91e, data/antibacterial.fasta, SHA-256 cbbeac64ba95746d87961e8ad9dd0849ae8058d15a300b2e7f6990730ca521e9. It contains 39,448 normalized canonical strings with aggregated upstream database attribution. The public release fetches this exact file from its pinned upstream URL instead of mirroring it. Original headers are preserved in the fetched cache; the organizer BSD notice remains. It is not a quantitative potency training set. Organizer repository licensing does not relicense every upstream database.

The final library and top 100 have no duplicates, reverse pairs or palindromes. Neither forward nor reversed strings exactly match the pinned organizer reference or the 7,795 comparable canonical 8–50 strings extracted from the additional DRAMP snapshot. Top-100 maximum forward reference Levenshtein ratio is 0.6923076923; the separate reversed-string diagnostic maximum is 0.6486486486. These checks cover identified frozen snapshots, not every current public AMP database. No all-repository or patent-novelty claim is made.

Full-tier release: TRAINING_DATA.md maps all source fitting data and processed screening labels. LICENSE_SCOPE.md distinguishes project MIT grants from retained upstream licenses.
