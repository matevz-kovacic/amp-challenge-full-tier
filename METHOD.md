# Method and fixed selection procedure

Ten completed runs supply 500,000 records. Ordinary FASTA formatting is normalized; identity is the exact N-to-C sequence with compatible canonical linear L chemistry and free, unmodified termini. Two cross-run duplicate records leave 499,998 candidates. Reverse sequences remain different peptides; equal composition is not identity. Source occurrences and original identifiers are retained in `selection/data/source_occurrences.csv.gz`.

The fixed common panel scores all 499,998 hard-eligible candidates identically:

| Evaluator | Direction | Meaning and dependence |
|---|---|---|
| HemoPI2-derived logistic classifier | Lower | Uncalibrated probability of the published strong-hemolysis label; composition plus four descriptors |
| HemoPI2-derived ExtraTrees classifier | Lower | Same endpoint and training labels; composition, descriptors and dipeptides |
| Physicochemical desirability | Higher | Explicit heuristic, negative squared deviations from chosen length, charge fraction, hydrophobic fraction and binary helical-moment targets |

The two classifiers form one family; the heuristic forms another. For evaluator j, oriented average-tie ranks r on the full fixed population N become b=(N-r)/(N-1). Identical evaluator copies collapse; globally constant evaluators are omitted; nonfinite scores fail. N=1 is explicitly handled. Within-family Borda means are averaged with equal family weights. Thus Q=0.25 b_logistic + 0.25 b_trees + 0.50 b_heuristic in this panel. Exact Q ties resolve lexicographically by sequence. Q is an internal ordinal preference, not MIC, HC50, a probability, or the organizer aggregate. Multiple source occurrences give no additional vote.

The classifier ensemble's existing homology-filtered historical validation has n=100 and AUROC 0.748 (bootstrap interval 0.655–0.843). It does not establish transfer to newly generated peptides or the challenge assays. No common efficacy predictor with adequate reproducible evidence was available within the bounded existing panel. The heuristic has no validated activity direction. The panel is explicitly MIXED_PREDICTIVE_HEURISTIC. A median-rank diagnostic using the same three columns has 0/100 leading overlap; its different family emphasis and preference for short sequences are disclosed, without changing the Borda default.

Hard library eligibility: canonical alphabet, length 8–50, exact uniqueness, declared allowed chemistry, and no exact match to the pinned organizer reference. The top screen additionally rejects Levenshtein.ratio strictly greater than 0.8 against that reference; equality passes. This is normalized-indel similarity, not unit-cost normalized edit similarity. Ten thousand eligible candidates were found after examining 10,001 quality leaders, with one reference rejection. All compared candidates retain their original full-population Q.

Greedy MMR selects 100 candidates plus 100 reserves from that fixed screen, with lambda=0.90. The first choice maximizes Q; subsequent choices maximize 0.90 Q(x) minus 0.10 times the largest similarity to a previously selected candidate. Assembly similarity is 1 minus unit-cost Levenshtein distance divided by maximum sequence length, in the normal N-to-C direction. Ties use Q then sequence. Export order is original Q, with insertion order separately recorded in `selection/selected_top100.csv`. Reserves are separate and are not necessarily library members.

MMR retains 79 of the quality-only top 100. For all 100, mean Q is 0.960944 versus 0.962126 and mean pairwise normalized edit distance is 0.832243 versus 0.816748. For the first 50 in quality order, mean Q is 0.966800 and diversity 0.813031. Both sizes are reported because public sampling descriptions have differed.

The library reserves the selected 100, then uses full-pool standardized length plus 20 amino-acid fractions as limited coverage coordinates. scikit-learn 1.7.2 MiniBatchKMeans uses seed 0 and 256 clusters (498,071 distinct feature vectors). n_init=3, batch_size=4096, init_size=12288, max_iter=100, max_no_improvement=10, tol=0, reassignment_ratio=0.01, one numerical thread. Remaining slots are proportional to each cluster's remaining population using capacity-limited integer largest remainders; cluster-ID ties are deterministic. Within each cluster, fixed library Q determines selection. This is neither a biological-family classifier nor a run quota.

The final library's mean Q is 0.756913. Five deterministic SHA-256-seeded 1,024-member diversity samples give mean 0.819976 (sample SD 0.001303); the reference comparison is 0.854845. Mean length is 19.3385 versus 21.5586 in the full eligible pool. These describe the actual selected set, without claiming that greater diversity or reference conformity guarantees qualification. Full metrics, distribution changes, cluster allocation and source contributions are under `selection/`.

Default inference regenerates every originating source pool before applying the frozen index map. Historical per-run top-list scoring is superseded; run03/run04 adapters execute the unchanged original generation prefix and stop before old shortlist ranking. All source-pool hashes and final byte hashes are checked. Replaying the fixed selection does not retrain models or refit the comparison population.

Attribution: Borda/mean-rank aggregation — Li, Wang and Xiao, Briefings in Bioinformatics 20(1), 178–189, DOI 10.1093/bib/bbx101. MMR — Carbonell and Goldstein, SIGIR 1998, DOI 10.1145/290941.291025. These established methods are adapted here; family grouping, equal weights, lambda, feature representation and proportional clustering policy are declared engineering choices.
