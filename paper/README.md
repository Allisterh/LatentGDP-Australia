# Working paper

`latentgdp_wp.qmd` is a Quarto working paper on the LatentGDP model family. It
renders to `latentgdp_wp.pdf` (and `latentgdp_wp.html`).

## Files

| File | What it is |
|---|---|
| `latentgdp_wp.qmd` | The paper. Every number in the text is computed by inline R when it renders. |
| `references.bib` | Bibliography. Each entry was checked against Crossref, the publisher or an official page. |
| `build_inputs.R` | Extracts what the paper needs from the local fit objects into `data/`. |
| `fit_rdp_sample.R` | Fits m2, m4, m5, m9, m12, m13 on RDP 2014-12's sample, 1980 Q1 - 2013 Q2 (`--end`, `--out` for another end date; output only under `output/sample_rdp/`). |
| `ppc_m13.R` | Prior and posterior predictive checks of the production model (writes `data/ppc_*.csv` and `output/ppc/ppc_draws.rds`). |
| `data/*.csv` | Small extracted inputs (each well under 1 MB). Safe to commit. |
| `figures/*.png` | Figures, written by the .qmd when it renders. |
| `latentgdp_wp.pdf` | Rendered paper. |

The .qmd reads only `paper/data/*.csv`, `output/tables/*.csv`,
`output/sample_1980/tables/*.csv` and `docs/*.json|csv`. It never reads a fit
object.

Three CSVs in `data/` are transcriptions, not extractions:
- `rdp2014_reported.csv`: values published in RBA RDP 2014-12;
- `record_facts.csv`: facts recorded only in NOTES.md or Stan file headers,
  each with its source;
- `test_counts.csv`: parsed from `output/sample_1980/logs/test_*.log`.

## Rebuild

From the repository root:

```bash
# 1. Only if the fits changed: re-extract the inputs (about 1 minute; reads
#    output/fit_production.rds, output/fit_m*.rds, output/sample_1980/fit_m*.rds,
#    output/fit_m13_seed2.rds, output/fit_m13_omegag1.rds,
#    output/sample_1980/fit_m13_omegag1.rds, output/sample_rdp/**/fit_m*.rds and
#    data/processed/gdp_growth.rds; writes only paper/data/).
Rscript paper/build_inputs.R

# 1b. The analyses added for the paper (about 45 minutes of sampling on 4 cores;
#     sequential and resumable): RDP-sample fits, fits ending in 2019 Q4, and the
#     omega_g ~ half-N(0, 1) refit of the production model on 1980 Q1 - 2026 Q2.
nohup bash output/sample_rdp/run_new_analyses.sh >> output/sample_rdp/logs/driver.log 2>&1 &
nohup bash output/sample_rdp/run_end2019.sh >> output/sample_rdp/logs/driver.log 2>&1 &
Rscript paper/ppc_m13.R          # seconds; needs output/fit_production.rds

# 2. Render (PDF and HTML).
quarto render paper/latentgdp_wp.qmd
# or one format only:
quarto render paper/latentgdp_wp.qmd --to pdf
```

Requirements:
- Quarto (tested with 1.9) and tectonic (the PDF engine; xelatex/pdflatex are
  not needed);
- R with ggplot2, patchwork, jsonlite and knitr; build_inputs.R also needs
  cmdstanr and posterior;
- the PDF uses the STIX Two Text and STIX Two Math fonts. On macOS they are
  system fonts. Elsewhere install them (they are free) or change `mainfont` /
  `mathfont` in the .qmd header.

To regenerate the fits and tables themselves, see the repository README
("Local reproduction") and the paper's appendix on reproduction.
