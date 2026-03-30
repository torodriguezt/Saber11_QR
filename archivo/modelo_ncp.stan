/*
 * Modelo completo con Non-Centered Parameterization (NCP)
 * Para MCMC (HMC/NUTS) con datos completos
 *
 * Cambios vs modelo_completo.stan original:
 *   1. Non-centered parameterization para alpha y mu
 *   2. Priors mejorados (half-normal en vez de inv_gamma)
 *   3. Mantiene skew_double_exponential directamente
 */

data {
  int<lower=1> N;
  int<lower=1> J;
  int<lower=1> K;
  int<lower=1> H;
  int<lower=1> E;
  int<lower=1> N2;
  int<lower=1> I;
  int<lower=1> L;
  int<lower=1> M;
  int<lower=1> P;

  int<lower=1, upper=L> naturaleza[J];

  int<lower=1, upper=J> cole[N];
  int<lower=1, upper=K> muni[N];
  int<lower=1, upper=H> horas[N];
  int<lower=1, upper=E> estrato[N];
  int<lower=1, upper=I> internet[N];
  int<lower=1, upper=M> educ_madre[N];
  int<lower=1, upper=P> educ_padre[N];
  vector[N] y;
  real<lower=0, upper=1> tau;
  int<lower=1, upper=J> col2[N2];
  int<lower=1, upper=K> muni2[N2];
}

parameters {
  // Non-centered: z-scores
  vector[K]   mu_raw;
  vector[J]   alpha_raw;

  real        mu_global;
  real<lower=0> sigma_global;
  real<lower=0> sigma_cole;
  real<lower=0> sigma;

  vector[H]   horas_effect;
  vector[E]   estrato_effect;
  vector[I]   internet_effect;
  vector[M]   educ_madre_effect;
  vector[P]   educ_padre_effect;
  vector[L]   naturaleza_effect;
}

transformed parameters {
  // Reconstruir desde z-scores
  vector[K] mu = mu_global + sigma_global * mu_raw;

  vector[J] alpha;
  for (j in 1:N2)
    alpha[col2[j]] = mu[muni2[j]] + sigma_cole * alpha_raw[col2[j]];

  // Centrado sum-to-zero
  vector[H] horas_cen       = horas_effect       - mean(horas_effect);
  vector[E] estrato_cen     = estrato_effect     - mean(estrato_effect);
  vector[I] internet_cen    = internet_effect    - mean(internet_effect);
  vector[M] educ_madre_cen  = educ_madre_effect  - mean(educ_madre_effect);
  vector[P] educ_padre_cen  = educ_padre_effect  - mean(educ_padre_effect);
  vector[L] naturaleza_cen  = naturaleza_effect  - mean(naturaleza_effect);
}

model {
  // Priors
  mu_global    ~ normal(247, 80);
  sigma_global ~ normal(0, 50);
  sigma_cole   ~ normal(0, 50);
  sigma        ~ normal(0, 30);

  mu_raw    ~ std_normal();
  alpha_raw ~ std_normal();

  horas_effect        ~ normal(0, 30);
  estrato_effect      ~ normal(0, 30);
  internet_effect     ~ normal(0, 30);
  educ_madre_effect   ~ normal(0, 30);
  educ_padre_effect   ~ normal(0, 30);
  naturaleza_effect   ~ normal(0, 30);

  // Verosimilitud
  for (i in 1:N)
    y[i] ~ skew_double_exponential(
              alpha[cole[i]]
            + estrato_cen[estrato[i]]
            + internet_cen[internet[i]]
            + horas_cen[horas[i]]
            + educ_madre_cen[educ_madre[i]]
            + educ_padre_cen[educ_padre[i]]
            + naturaleza_cen[ naturaleza[ cole[i] ] ],
            2 * sigma, tau);
}
