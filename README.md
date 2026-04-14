# emoDAG_demo

Reanalysis of Lin et al. (2026), "The Dynamic Interplay Between Goal Setting, Performance, and Emotions in Self-Regulated Learning" (*Journal of Educational Psychology*). 

## What this project does

I reanalyze Study 1 of Lin et al. (2026) and show that:

1. The reported negative emotion-to-performance effect is not robust to the inclusion of a log(trial) learning curve.
2. WAIC and PSIS-LOO, as used in the original study, are diagnostically unreliable and theoretically biased for this class of autoregressive models. Exact leave-future-out cross-validation (LFO-CV) provides an unbiased alternative.
3. Under LFO-CV, the original model selection conclusion (Alt1 over Hypothesized) is no longer supported.

## Data

The data are from Lin et al. (2026) and are publicly available at [https://osf.io/qc6hm/](https://osf.io/qc6hm/). Download `data.csv` and place it in `data/` before running the scripts.

## Project structure

```
├── stan/                      # Stan model files
│   ├── hyp_optimized.stan     # Hypothesized model (symmetric GPD)
│   ├── alt1_optimized.stan    # Alternate 1 (split success/failure)
│   └── mine1_logtrial.stan    # LogTrial (Alt1 + log(t) effects)
├── run_script/
│   ├── run.R                  # Fit all models, WAIC, PSIS-LOO, parameter estimates
│   ├── run_lfo.R              # Exact LFO-CV (one model at a time)
│   └── compare_lfo.R          # Aggregate and compare LFO results
├── figure_script/
│   ├── fig2_badk.R            # Pareto k diagnostic plots
│   ├── fig3_LFOCV.R           # LFO per-block ΔELPD bar charts
│   └── fig4_forest.R          # Standardized posterior forest plot
├── models/                    # Saved model fits and LFO results (.rds, .csv)
├── figure/                    # Output figures (.pdf, .png)
├── paper/                     # Manuscript (Quarto/apaquarto)
├── csf/                       # HPC job scripts (University of Manchester CSF)
└── data/                      # Place data.csv here (not tracked; see above)
```

## How to reproduce

**Prerequisites:** R (≥ 4.3), [CmdStan](https://mc-stan.org/cmdstanr/), and the following R packages: `cmdstanr`, `posterior`, `loo`, `tidyverse`, `ggdist`, `ggtext`, `patchwork`.

**Step 1: Fit models and compute WAIC/LOO.**

```r
source("run_script/run.R")
```

This fits all three models and prints WAIC, PSIS-LOO, stacking weights, and parameter estimates. Key results are recorded in the comments of the script.

**Step 2: Run exact LFO-CV.**

```r
# One model at a time (each takes ~8 hours):
Rscript run_script/run_lfo.R hyp
Rscript run_script/run_lfo.R alt1
Rscript run_script/run_lfo.R mine1

# Then aggregate:
source("run_script/compare_lfo.R")
```

**Step 3: Generate figures.**

```r
# After run.R and compare_lfo.R have been run:
source("figure_script/fig2_badk.R")
source("figure_script/fig3_LFOCV.R")
source("figure_script/fig4_forest.R")
```

Most numerical results are available as comments in `run_script/run.R` and `run_script/run_lfo.R` without needing to rerun the models. Pre-computed model fits are saved in `models/`.

## Reference

Lin, W. M., FitzGibbon, L., Theobald, M., Breitwieser, J., Brod, G., Murayama, K., & Sakaki, M. (2026). The dynamic interplay between goal setting, performance, and emotions in self-regulated learning: A computational modeling approach. *Journal of Educational Psychology*. https://doi.org/10.1037/edu0001022