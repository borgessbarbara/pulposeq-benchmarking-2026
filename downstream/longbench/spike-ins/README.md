# Synthetic RNA controls (Sequins + SIRVs + ERCCs) analysis

LongBench spike-ins, per protocol (ONT cDNA, ONT DRS, PacBio), on the oarfish runs.

## `01_sequin_edger.R` - sequin ground truth (Fig. 1d)
- Contrast Mix B (SCLC-A) vs Mix A; expected log2FC from Anaquin's `RnaQuinIsoformMixture`.
- Counts: oarfish, via `read_oarfish()` (`../../helpers/oarfish_counts.R`):
  transcript level `catchOarfish()`, gene level `catchOarfishGene()`, both divided by their
  RTA overdispersion (`count_mode = "divided"`, primary). Gene level is also run on undivided
  summed counts (`count_mode = "undivided"`) as a sensitivity check.
- Filter >= 10 reads in >= 2 samples on the undivided counts; library sizes and TMM factors
  from the ratio-1 control sequins; edgeR glmQLFit (robust) -> glmTreat (|log2FC| > 0.5), FDR 0.05.
- Scoring: truly DE = |expected log2FC| > 0.5; wrong-sign calls are false discoveries.
- Identification: detected = >= 1 read in the curated oarfish table; also reported for every
  model oarfish quantified and for Bambu's raw table; novel models on `chrIS` (raw vs curated).
- Outputs (`results/sequins/`): `edger/<platform>/<level>/<count_mode>/sequin_de.csv`,
  `edger/sequin_de_results.csv`, `edger/fold_change_summary_edger.csv` (with
  `median_overdispersion`), `sequin_identification.csv`, `novel_models_on_sequins.csv`.

## `02_sirv_ercc.R` - SIRV-Set 4 and ERCC (Fig. 1b, c)
- SIRV E0 isoform accuracy, long-SIRV length bias, ERCC dose response and detection limit,
  novel models on spike-in sequences, SIRV identification.
- Counts: curated oarfish counts, **undivided** (accuracy, nothing is tested);
  `sirv_e0_per_isoform.csv` adds each isoform's RTA overdispersion and the summary its median.
- Outputs (`results/sirv_ercc/`): `sirv_ercc_summary.csv`, `length_bias.csv`,
  `ercc_per_transcript.csv`, `novel_models_on_spikeins.csv`, `<platform>/sirv_e0_per_isoform.csv`.
