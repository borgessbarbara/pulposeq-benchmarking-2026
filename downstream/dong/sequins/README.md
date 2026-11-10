# Synthetic RNA controls (Sequins) analysis with Anaquin

`02_sequin_edger.R`: sequin ground truth for the Dong DE (Fig. 2c). H1975 + Mix A,
HCC827 + Mix B; expected log2FC = log2(Mix A / Mix B).

- Same test, filter and scoring as `diff-exp/01_de_edger.R` and LongBench's sequin script.
- Counts: oarfish, divided by the RTA overdispersion (`count_mode = "divided"`, primary);
  gene level also undivided (`count_mode = "undivided"`) as a sensitivity check. The filter uses
  the undivided reads, so both count modes test the same sequins.
- Outputs (`results/`): `<level>/<count_mode>/sequin_de.csv`, `sequin_de_results.csv`,
  `summary.csv` (one row per level x count mode).
