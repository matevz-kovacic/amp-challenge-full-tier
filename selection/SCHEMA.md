# Frozen selection data

All arrays use zero-based row indices into data/pool_sequences.txt.gz, which contains 499,998 unique normalized N-to-C strings sorted lexicographically. Chemistry is linear_L_free_N_free_C_unmodified for this assessed pool. No source-run ordering affects these indices.

- raw_scores.npy: float64 N x 3; columns hemo_logistic, hemo_trees, physicochemical_desirability.
- rank_cache.npz: Q (float64), rank (1-based quality rank), quality_order (row-index permutation), median_order (diagnostic permutation).
- cluster_labels.npy: integer cluster per pooled row; cluster_model.npz: saved fitted scaler/centers and representation data. clustering_report.json and selection_config.json pin the fit.
- mmr_screen.npy: ordered 10,000 top-eligible pooled row indices.
- top_eligibility.csv.gz: exact stricter reference checks for all screened contenders, including rejection reason and reference ratio.
- source_occurrences.csv.gz: one row per original input record, including identity digest, source run, record number, original header, raw joined sequence, previous top rank, chemistry and rejection reasons.
- selected_top100.csv: quality export position, original Q/rank, MMR insertion order/utility, source provenance, cluster and novelty evidence.
- reserve_candidates.csv: subsequent MMR reserve priority, original Q/rank and membership flag; reserve candidates are separate from the chosen 100.
- data_manifest.json: SHA-256 pins for the selection inputs. Root frozen_selection.json: sequence hash and source-run/zero-based source index plus quality rank for every output member; it stores no output strings.

Running replay_selection.py expands raw scores, per-evaluator average ranks/Borda, family contributions, Q and quality rank into replayed/expanded_rank_table.npz, in the same pool row order. It reproduces the submitted FASTA hashes without predictor inference. Saved score/rank files remain read-only to the replay.

Some copied historical configuration records retain hashes naming original evidence files that are not duplicated here; original_pin_hashes identifies those historical pins. The root replay_manifest and this directory's data_manifest define the actual packaged inference/audit assets. The scalar source code and checkpoint are pinned independently. No claim is made that old per-run standalone top-list branches are executable from this reduced package.

reverse_pairs.csv records the one exact reverse relationship in the full source pool. Both strings remain separate N-to-C identities. Neither appears in the final library.
