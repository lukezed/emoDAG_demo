# prep_study2.R — build Study 2 analysis data in Study-1 format
#
# Usage: Rscript run_script/prep_study2.R [pre|post]   (default: pre)
# Output: data/data2_<group>.csv (subj, block, goals, perfs, emots),
# so every existing Study-1 Stan model and run script works unchanged.
#
# Follows Lin et al. Study 2 cleaning (Supplement, "Additional Considerations"):
#   - performance = n_unique_qs_shifted (questions attempted)
#   - drop days with missing/zero performance ("no studying"); re-sequence
#   - emotion index valence-aggregated (verified vs Table 4: pre 0.66, post 0.43)
#   - keep participants with >= 20 complete days
#
# ponytail: participants with a missing GOAL (pre) or missing EMOTION (post,
# 163 ppl) on some days are dropped row-wise rather than imputed as Stan
# parameters (Lin's approach). Keeps the Stan model identical to Study 1.
# Upgrade path: port Lin's missing-data imputation if a reviewer asks.

library(tidyverse)

group <- (commandArgs(trailingOnly = TRUE)[1] %||% "pre")
stopifnot(group %in% c("pre", "post"))

raw <- read_csv("data/data2.csv", show_col_types = FALSE) |>
  mutate(day = as.integer(day))

# ---- Emotion index (Lin's recipe: reverse-code negatives 6-x, centre at 3) ----
# pre  items: enjoyment_b, anger_b, tension_b
# post items: enjoyment_a, pride_a, anger_a, tension_a
emo_index <- if (group == "pre") {
  function(d) rowMeans(cbind(d$enjoyment_b - 3, (6 - d$anger_b) - 3, (6 - d$tension_b) - 3))
} else {
  function(d) rowMeans(cbind(d$enjoyment_a - 3, d$pride_a - 3,
                             (6 - d$anger_a) - 3, (6 - d$tension_a) - 3))
}

# Goals and performance are question COUNTS (~90); the emotion index is ~1.
# Lin scaled counts down so all three variables share a comparable scale (their
# reported sigma_p = 0.388, gamma_p = 0.269, gamma_g = 0.608 are only recoverable
# with counts / 100). The factor is undocumented; 100 is recovered by matching
# those three quantities. Emotion index kept raw (Table 4 SD 0.77). Without this,
# the N(0,3) priors and NCP fail to converge (emotion coefficients blow up).
SCALE <- 100

# ---- Clean + re-sequence -----------------------------------------------------
clean <- raw |>
  mutate(
    perf  = n_unique_qs_shifted / SCALE,
    goal  = goal_qs / SCALE,
    emot  = emo_index(raw)
  ) |>
  # drop "no studying" days and rows with any missing modelled variable
  filter(!is.na(perf), perf > 0, !is.na(goal), !is.na(emot)) |>
  arrange(record_id, day)

# Pre-learning emotion is measured in the MORNING, before that day's task, so it
# should predict the SAME day's goal/performance (mood-as-information), and be a
# result of the PREVIOUS day's outcome (Supplement: "E_B at t+1 estimated from
# goal/performance at t"). Straight alignment puts the emotion one step early and
# collapses beta_E->P to ~0; leading the emotion by one kept session recovers
# Lin's +0.018. Post-learning emotion is measured after the task and needs no
# shift (it is already a result of the same day's GPD).
if (group == "pre") {
  clean <- clean |>
    group_by(record_id) |>
    mutate(emot = lead(emot)) |>          # emots[s] <- pre-emotion(s+1)
    filter(!is.na(emot)) |>               # drops each person's last session
    ungroup()
}

clean <- clean |>
  arrange(record_id, day) |>
  group_by(record_id) |>
  mutate(block = row_number()) |>          # re-sequence: gaps treated as consecutive
  ungroup()

# ---- Exclude < 20 complete days ----------------------------------------------
keep_ids <- clean |> count(record_id) |> filter(n >= 20) |> pull(record_id)

study2 <- clean |>
  filter(record_id %in% keep_ids) |>
  mutate(subj = as.integer(factor(record_id))) |>
  transmute(subj, block, goals = goal, perfs = perf, emots = emot) |>
  arrange(subj, block)

outfile <- sprintf("data/data2_%s.csv", group)
write_csv(study2, outfile)

# ---- Report / validation -----------------------------------------------------
cat(sprintf("Group: %s  ->  %s\n", group, outfile))
cat(sprintf("Participants kept: %d  (Lin reports 305)\n", n_distinct(study2$subj)))
cat(sprintf("Total modelled rows: %d\n", nrow(study2)))
cat(sprintf("Blocks/person: min %d median %.0f max %d\n",
            min(table(study2$subj)), median(table(study2$subj)), max(table(study2$subj))))

cat(sprintf("\nEmotion index (%s) mean = %.2f  (Lin Table 4: pre 0.66 / post 0.43)\n",
            group, mean(study2$emots)))
cat(sprintf("cor(goal, perf) = %.2f  (Lin Table 4: .46)\n",
            cor(study2$goals, study2$perfs)))

# --- THE KEY QUESTION: is there a practice trend in Study 2? -------------------
# regress each outcome on log(session), within-person centred, pooled slope
trend <- study2 |>
  group_by(subj) |>
  mutate(lt = log(block),
         perf_c = perfs - mean(perfs),
         goal_c = goals - mean(goals),
         emot_c = emots - mean(emots),
         lt_c   = lt - mean(lt)) |>
  ungroup()

slope <- function(y, x) sum(x * y) / sum(x * x)   # within-person OLS slope, no intercept
cat("\n--- Within-person trend on log(session) (pooled) ---\n")
cat(sprintf("performance ~ log(session): slope = %+.3f   <-- practice trend?\n",
            slope(trend$perf_c, trend$lt_c)))
cat(sprintf("goal        ~ log(session): slope = %+.3f   (Lin: goals decrease)\n",
            slope(trend$goal_c, trend$lt_c)))
cat(sprintf("emotion     ~ log(session): slope = %+.3f\n",
            slope(trend$emot_c, trend$lt_c)))

cat(sprintf("\nWrote %s\n", outfile))
