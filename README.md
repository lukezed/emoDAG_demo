# emoDAG_demo

Reanalysis of Lin et al. (2026), "The Dynamic Interplay Between Goal Setting, Performance, and Emotions in Self-Regulated Learning" (*Journal of Educational Psychology*). 

## What this project does

I reanalyze Lin et al. (2026), Study 1 in full and Study 2 as an extension, and show that:

1. Study 1's reported negative emotion-to-performance effect is not robust to a log(trial) learning curve; a practice effect accounts for it.
2. The apparent sign reversal of the emotion-to-performance effect across the two studies reflects practice confounding rather than opposing emotional mechanisms. Once a time trend is modelled, the effect is near-null in Study 1 and a small positive in Study 2.
3. WAIC and PSIS-LOO, as used in the original study, are diagnostically unreliable for this class of autoregressive models; exact leave-future-out cross-validation (LFO-CV) is used instead. Under LFO-CV the goal-revision asymmetry (Alt1) is supported, while the emotion-to-performance conclusion is not.

## Data

Study 1 data are from Lin et al. (2026), publicly available at [https://osf.io/qc6hm/](https://osf.io/qc6hm/); download `data.csv` into `data/`. Study 2 data are from Breitwieser et al. (2022) (the source dataset for Lin et al.'s Study 2); place the raw file as `data2.csv` in `data/`, then run `s2_prep.R`. Neither dataset is redistributed in this repository.

## Project structure

```
├── stan/                      # Stan model files (2x2: symmetric/split GPD x no/with log-trial)
│   ├── hyp_optimized.stan     # Hyp   — symmetric GPD, no log-trial
│   ├── alt1_optimized.stan    # Alt1  — split success/failure, no log-trial
│   ├── mine1_logtrial.stan    # LogT  — Alt1 + log(t) effects
│   └── mine2_logtrial_sym.stan# Hyp + log(t) effects (completes the 2x2)
├── run_script/                # named <study>_<step>: s1_/s2_ by study, lfo_ shared
│   ├── s1_fit.R               # Study 1: fit the 2x2 (hyp/alt1/mine1/mine2), WAIC, PSIS-LOO, parameters
│   ├── s2_prep.R              # Study 2: build analysis data (scaling + emotion index)
│   ├── s2_fit.R               # Study 2: fit the 2x2 for pre/post emotions
│   ├── lfo_fit.R              # Exact LFO-CV, one model at a time (either study)
│   └── lfo_compare.R          # Aggregate and compare LFO results
├── figure_script/
│   ├── fig2_badk.R            # Pareto k diagnostic plots
│   ├── fig3_LFOCV.R           # LFO per-block ΔELPD bar charts
│   └── fig4_forest.R          # Standardized posterior forest plot
├── models/                    # Saved model fits and LFO results (.rds, .csv)
├── figure/                    # Output figures (.pdf, .png)
├── paper/                     # Manuscript (Quarto/apaquarto)
├── csf/                       # HPC job scripts (University of Manchester CSF)
└── data/                      # Place data.csv (Study 1) / data2.csv (Study 2) here (not tracked)
```

## How to reproduce

**Prerequisites:** R (≥ 4.3), [CmdStan](https://mc-stan.org/cmdstanr/), and the following R packages: `cmdstanr`, `posterior`, `loo`, `tidyverse`, `ggdist`, `ggtext`, `patchwork`.

**Step 1: Fit Study 1 models and compute WAIC/LOO.**

```r
source("run_script/s1_fit.R")   # fits the 2x2: hyp, alt1, mine1, mine2
```

This prints WAIC, PSIS-LOO, stacking weights, and parameter estimates. Key results are recorded in the comments of the script.

**Step 2: Run exact LFO-CV (Study 1).**

```r
# One model at a time (each takes ~8 hours):
Rscript run_script/lfo_fit.R hyp
Rscript run_script/lfo_fit.R alt1
Rscript run_script/lfo_fit.R mine1
Rscript run_script/lfo_fit.R mine2

# Then aggregate:
source("run_script/lfo_compare.R")
```

**Step 3: Study 2 reanalysis.**

```r
# Build Study 2 data (place data2.csv in data/ first), then fit the 2x2:
Rscript run_script/s2_prep.R pre      # -> data/data2_pre.csv
Rscript run_script/s2_prep.R post     # -> data/data2_post.csv
Rscript run_script/s2_fit.R pre       # fit the pre-emotion 2x2
Rscript run_script/s2_fit.R post alt1 mine1   # post-emotion control
```

The Study 2 pipeline recovers Lin et al.'s reported Study 2 baseline; results are recorded in the comments of `s2_fit.R`.

**Step 4: Generate figures.**

```r
# After s1_fit.R and lfo_compare.R have been run:
source("figure_script/fig2_badk.R")
source("figure_script/fig3_LFOCV.R")
source("figure_script/fig4_forest.R")
```

Most numerical results are available as comments in `run_script/s1_fit.R`, `lfo_fit.R`, and `s2_fit.R` without rerunning the models. Pre-computed model fits are saved in `models/`.

## Reference

Lin, W. M., FitzGibbon, L., Theobald, M., Breitwieser, J., Brod, G., Murayama, K., & Sakaki, M. (2026). The dynamic interplay between goal setting, performance, and emotions in self-regulated learning: A computational modeling approach. *Journal of Educational Psychology*. https://doi.org/10.1037/edu0001022