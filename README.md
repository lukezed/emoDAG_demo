# emoDAG_demo

Code and results for "Does Positive Affect Impair Subsequent Performance? A Practice Confound in Learning Studies" (Zhang, under review).

The manuscript refits the computational model of goals, performance, and affect of Lin et al. (2026) to both of their datasets, adds practice terms (log of the trial number) to each equation, compares the models by leave-one-group-out (LOGO) cross-validation, and uses simulation to test whether omitting the practice terms biases the affect-to-performance path.

## Main results

- **Study 1** (arithmetic task, 357 adults, 20 blocks). The affect-to-performance path in Alt1 was −0.104 [−0.174, −0.036] without practice terms and −0.043 [−0.114, 0.028] with them.
- **Study 2** (examination preparation, 303 medical students, up to 40 sessions). The path stayed positive: 0.017 [0.001, 0.033] without practice terms and 0.020 [0.004, 0.037] with them.
- **Prediction.** The practice terms improved prediction for held-out participants in both studies (LOGO ΔELPD 34 to 49, each more than three standard errors).
- **Simulation.** With a true path of zero, the model without practice terms gave a negative estimate in 90 of 100 datasets and a credibly negative estimate in 18. The model with practice terms recovered both zero and a true path of −0.10.
- **Time coding.** The results were the same with a linear trend or a separate indicator for each block.

## Data

Neither dataset is included in this repository.

- **Study 1** is openly available from the repository of Lin et al. at <https://osf.io/qc6hm/>. Download `data.csv` into `data/`.
- **Study 2** comes from Breitwieser et al. (2022) and was shared by its authors under a license that does not permit redistribution. With the raw file saved as `data/data2.csv`, `run_script/s2_prep.R` builds the analysis data.

## Structure

```
stan/              Stan models
  hyp_optimized.stan       Hyp: symmetric goal-performance coefficient
  alt1_optimized.stan      Alt1: separate coefficients for goal success and failure
  mine2_logtrial_sym.stan  Hyp + log t
  mine1_logtrial.stan      Alt1 + log t
  mine1_timeX.stan         Alt1 with a linear or block-indicator coding of time
  sim_alt1.stan            Alt1 for the simulation fits
  sim_mine1.stan           Alt1 + log t for the simulation fits
run_script/
  s1_fit.R                 Study 1: fit the four models; WAIC, PSIS-LOO, parameters
  s2_prep.R                Study 2: build the analysis data
  s2_fit.R                 Study 2: fit the four models (pre- or post-session affect)
  logo.R                   LOGO cross-validation (PSIS with exact refits)
  time_coding.R            alternative codings of time
  s1_sim_check.R           simulation with a known affect-to-performance path
  s1_sim_describe.R        simulated vs real Study 1 data (Appendix D)
figure_script/     figure code (R and TikZ)
figure/            figures used in the manuscript
models/            result tables (.csv); model fits (.rds) are not tracked
archive/lfo/       earlier leave-future-out analysis, not used in the manuscript
```

In the scripts, `mine1` is Alt1 + log t and `mine2` is Hyp + log t.

## Reproducing the results

Requirements: R (≥ 4.3), CmdStan via [`cmdstanr`](https://mc-stan.org/cmdstanr/), and the packages `posterior`, `loo`, `tidyverse`, `ggdist`, `ggtext`, and `patchwork`. Run all scripts from the repository root. Each script's header documents its arguments, and key results are recorded in comments.

```sh
# Study 1: fit the four models (saved to models/fit_*.rds)
Rscript run_script/s1_fit.R

# Study 2: build the data, then fit the four models
Rscript run_script/s2_prep.R pre
Rscript run_script/s2_prep.R post
Rscript run_script/s2_fit.R pre
Rscript run_script/s2_fit.R post alt1 mine1

# LOGO cross-validation (Table 2); "refit" runs the exact refits
Rscript run_script/logo.R s1 refit
Rscript run_script/logo.R s2pre refit

# Alternative codings of time (Table 4)
Rscript run_script/time_coding.R s1 linear
Rscript run_script/time_coding.R s1 dummies
Rscript run_script/time_coding.R s2pre linear
Rscript run_script/time_coding.R s2pre dummies

# Simulation (Table 5, Appendix D): reps, fitted model, true path
Rscript run_script/s1_sim_check.R 1:100 alt1 0
Rscript run_script/s1_sim_check.R 1:100 mine1 0
Rscript run_script/s1_sim_check.R 1:100 mine1 -0.10
Rscript run_script/s1_sim_describe.R
```

Dataset `r` of the simulation uses seed `2026 + r`. Results append to the CSV files in `models/`, so an interrupted run can be resumed with the missing range.

Figures: `figure_script/fig4_forest.R`, `fig2_badk.R`, and `fig6_s2_badk.R` read the saved fits; the two diagrams are compiled with `xelatex figure_script/fig5_dag.tex` and `xelatex figure_script/fig5b_dag_estimates.tex`.

## References

Breitwieser, J., Neubauer, A. B., Schmiedek, F., & Brod, G. (2022). Self-regulation prompts promote the achievement of learning goals – but only briefly: Uncovering hidden dynamics in the effects of a psychological intervention. *Learning and Instruction, 80*, 101560. https://doi.org/10.1016/j.learninstruc.2021.101560

Lin, W. M., FitzGibbon, L., Theobald, M., Breitwieser, J., Brod, G., Murayama, K., & Sakaki, M. (2026). The dynamic interplay between goal setting, performance, and emotions in self-regulated learning: A computational modeling approach. *Journal of Educational Psychology*. https://doi.org/10.1037/edu0001022
