// ============================================================================
// mine1 — Alt1 + fixed log(trial) effect
//
// vs Alt1: +3 population params only, 0 individual params added
// Goal:  beta_gp_g_s / beta_gp_g_f (from alt1) + beta_lt_g * log(t)
// Perf:  beta_gp_p (unchanged)                  + beta_lt_p * log(t)
// Emot:  beta_gp_e (unchanged)                  + beta_lt_e * log(t)
// ============================================================================

data {
  int<lower=1> N;
  int<lower=1> T;
  array[T] int<lower=1,upper=N> id;
  array[T] int<lower=2> trial;
  vector[N] g1;
  vector[N] p1;
  vector[N] e1;
  vector[T] g;
  vector[T] p;
  vector[T] e;
}

transformed data {
  array[T] int is_first;
  vector[T] log_trial;
  for (i in 1:T) {
    is_first[i] = (trial[i] == 2);
    log_trial[i] = log(trial[i]);
  }
}

parameters {
  // Population means (same 12 as alt1)
  real pop_gamma_g;   real pop_alpha_g;
  real pop_beta_s;    real pop_beta_f;    real pop_beta_e_g;
  real pop_gamma_p;   real pop_alpha_p;   real pop_beta_gp_p;  real pop_beta_e_p;
  real pop_gamma_e;   real pop_alpha_e;   real pop_beta_gp_e;

  // Fixed log(trial) effects (+3)
  real beta_lt_g;
  real beta_lt_p;
  real beta_lt_e;

  // Population SDs (same 12 as alt1)
  real<lower=0> pop_gamma_g_sd;   real<lower=0> pop_alpha_g_sd;
  real<lower=0> pop_beta_s_sd;    real<lower=0> pop_beta_f_sd;    real<lower=0> pop_beta_e_g_sd;
  real<lower=0> pop_gamma_p_sd;   real<lower=0> pop_alpha_p_sd;
  real<lower=0> pop_beta_gp_p_sd; real<lower=0> pop_beta_e_p_sd;
  real<lower=0> pop_gamma_e_sd;   real<lower=0> pop_alpha_e_sd;   real<lower=0> pop_beta_gp_e_sd;

  // NCP offsets (N x 12, same as alt1)
  matrix[N, 12] z;

  // Residual SDs
  real<lower=0> sigma_g;
  real<lower=0> sigma_p;
  real<lower=0> sigma_e;
}

transformed parameters {
  vector[N] gamma_g   = pop_gamma_g   + pop_gamma_g_sd   * z[,1];
  vector[N] alpha_g   = pop_alpha_g   + pop_alpha_g_sd   * z[,2];
  vector[N] beta_gp_s = pop_beta_s    + pop_beta_s_sd    * z[,3];
  vector[N] beta_gp_f = pop_beta_f    + pop_beta_f_sd    * z[,4];
  vector[N] beta_e_g  = pop_beta_e_g  + pop_beta_e_g_sd  * z[,5];
  vector[N] gamma_p   = pop_gamma_p   + pop_gamma_p_sd   * z[,6];
  vector[N] alpha_p   = pop_alpha_p   + pop_alpha_p_sd   * z[,7];
  vector[N] beta_gp_p = pop_beta_gp_p + pop_beta_gp_p_sd * z[,8];
  vector[N] beta_e_p  = pop_beta_e_p  + pop_beta_e_p_sd  * z[,9];
  vector[N] gamma_e   = pop_gamma_e   + pop_gamma_e_sd   * z[,10];
  vector[N] alpha_e   = pop_alpha_e   + pop_alpha_e_sd   * z[,11];
  vector[N] beta_gp_e = pop_beta_gp_e + pop_beta_gp_e_sd * z[,12];
}

model {
  pop_gamma_g ~ normal(0, 3);  pop_alpha_g ~ normal(0, 3);
  pop_beta_s  ~ normal(0, 3);  pop_beta_f  ~ normal(0, 3);  pop_beta_e_g ~ normal(0, 3);
  pop_gamma_p ~ normal(0, 3);  pop_alpha_p ~ normal(0, 3);
  pop_beta_gp_p ~ normal(0, 3); pop_beta_e_p ~ normal(0, 3);
  pop_gamma_e ~ normal(0, 3);  pop_alpha_e ~ normal(0, 3);  pop_beta_gp_e ~ normal(0, 3);

  beta_lt_g ~ normal(0, 3);
  beta_lt_p ~ normal(0, 3);
  beta_lt_e ~ normal(0, 3);

  pop_gamma_g_sd ~ normal(0, 3);  pop_alpha_g_sd ~ normal(0, 3);
  pop_beta_s_sd  ~ normal(0, 3);  pop_beta_f_sd  ~ normal(0, 3);  pop_beta_e_g_sd ~ normal(0, 3);
  pop_gamma_p_sd ~ normal(0, 3);  pop_alpha_p_sd ~ normal(0, 3);
  pop_beta_gp_p_sd ~ normal(0, 3); pop_beta_e_p_sd ~ normal(0, 3);
  pop_gamma_e_sd ~ normal(0, 3);  pop_alpha_e_sd ~ normal(0, 3);  pop_beta_gp_e_sd ~ normal(0, 3);

  to_vector(z) ~ std_normal();

  sigma_g ~ normal(0, 3);
  sigma_p ~ normal(0, 3);
  sigma_e ~ normal(0, 3);

  for (i in 1:T) {
    int s = id[i];
    real pg, pp, pe, disc, mg, mp;
    real lt = log_trial[i];

    if (is_first[i]) {
      pg = g1[s]; pp = p1[s]; pe = e1[s];
    } else {
      pg = g[i-1]; pp = p[i-1]; pe = e[i-1];
    }

    disc = pp - pg;
    mg = gamma_g[s] + alpha_g[s] * pg
         + beta_gp_s[s] * fmax(disc, 0) + beta_gp_f[s] * fmin(disc, 0)
         + beta_e_g[s] * pe + beta_lt_g * lt;
    g[i] ~ normal(mg, sigma_g);

    mp = gamma_p[s] + alpha_p[s] * pp
         + beta_gp_p[s] * (g[i] - pp)
         + beta_e_p[s] * pe + beta_lt_p * lt;
    p[i] ~ normal(mp, sigma_p);

    e[i] ~ normal(gamma_e[s] + alpha_e[s] * pe
                  + beta_gp_e[s] * (p[i] - g[i])
                  + beta_lt_e * lt,
                  sigma_e);
  }
}

generated quantities {
  vector[T] log_lik;
  vector[T] sam_g;
  vector[T] sam_p;
  vector[T] sam_e;

  for (i in 1:T) {
    int s = id[i];
    real pg, pp, pe, disc, mg, mp, me;
    real lt = log_trial[i];

    if (is_first[i]) {
      pg = g1[s]; pp = p1[s]; pe = e1[s];
    } else {
      pg = g[i-1]; pp = p[i-1]; pe = e[i-1];
    }

    disc = pp - pg;
    mg = gamma_g[s] + alpha_g[s] * pg
         + beta_gp_s[s] * fmax(disc, 0) + beta_gp_f[s] * fmin(disc, 0)
         + beta_e_g[s] * pe + beta_lt_g * lt;
    mp = gamma_p[s] + alpha_p[s] * pp
         + beta_gp_p[s] * (g[i] - pp)
         + beta_e_p[s] * pe + beta_lt_p * lt;
    me = gamma_e[s] + alpha_e[s] * pe
         + beta_gp_e[s] * (p[i] - g[i])
         + beta_lt_e * lt;

    log_lik[i] = normal_lpdf(g[i] | mg, sigma_g)
               + normal_lpdf(p[i] | mp, sigma_p)
               + normal_lpdf(e[i] | me, sigma_e);
    sam_g[i] = normal_rng(mg, sigma_g);
    sam_p[i] = normal_rng(mp, sigma_p);
    sam_e[i] = normal_rng(me, sigma_e);
  }
}
