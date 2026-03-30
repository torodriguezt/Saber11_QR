/*
 * Modelo completo reparametrizado para Variational Bayes (ADVI / Pathfinder)
 *
 * La ALD (skew_double_exponential) tiene gradientes no-suaves que hacen
 * fallar a ADVI. Solucion: representar la ALD como mezcla normal-exponencial
 * (Kozumi & Kobayashi, 2011):
 *
 *   y_i | w_i ~ Normal(loc_i + delta * w_i,  sqrt(gamma_sq * w_i))
 *   w_i       ~ Exponential(sigma_ald)
 *
 * donde delta = sigma_ald * (1-2*tau) / (tau*(1-tau))
 *       gamma_sq = 2 * sigma_ald^2 / (tau*(1-tau))
 *
 * Esto es matematicamente equivalente a y_i ~ ALD(loc_i, sigma_ald, tau)
 * pero es completamente diferenciable -> VB funciona.
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

transformed data {
  // Constantes de la mezcla normal-exponencial para la ALD
  real tau_comp = tau * (1.0 - tau);
  real delta_unit = (1.0 - 2.0 * tau) / tau_comp;    // se multiplica por sigma_ald
  real gamma_sq_unit = 2.0 / tau_comp;                // se multiplica por sigma_ald^2
}

parameters {
  // --- Non-centered parameterization ---
  vector[K]   mu_raw;           // z-scores para municipios
  vector[J]   alpha_raw;        // z-scores para colegios

  real        mu_global;
  real<lower=0> sigma_global;
  real<lower=0> sigma_cole;
  real<lower=0> sigma_ald;      // escala de la ALD (reemplaza sigma)

  // Efectos covariables
  vector[H]   horas_effect;
  vector[E]   estrato_effect;
  vector[I]   internet_effect;
  vector[M]   educ_madre_effect;
  vector[P]   educ_padre_effect;
  vector[L]   naturaleza_effect;

  // Variables latentes de la mezcla (una por observacion)
  vector<lower=0>[N] w;
}

transformed parameters {
  // Reconstruir mu y alpha desde z-scores (non-centered)
  vector[K] mu = mu_global + sigma_global * mu_raw;

  vector[J] alpha;
  for (j in 1:N2)
    alpha[col2[j]] = mu[muni2[j]] + sigma_cole * alpha_raw[col2[j]];

  // Centrado de efectos (sum-to-zero)
  vector[H] horas_cen       = horas_effect       - mean(horas_effect);
  vector[E] estrato_cen     = estrato_effect     - mean(estrato_effect);
  vector[I] internet_cen    = internet_effect    - mean(internet_effect);
  vector[M] educ_madre_cen  = educ_madre_effect  - mean(educ_madre_effect);
  vector[P] educ_padre_cen  = educ_padre_effect  - mean(educ_padre_effect);
  vector[L] naturaleza_cen  = naturaleza_effect  - mean(naturaleza_effect);

  // Parametros de la mezcla (dependen de sigma_ald)
  real delta = sigma_ald * delta_unit;
  real gamma_sq = sigma_ald * sigma_ald * gamma_sq_unit;
}

model {
  // --- Priors ---
  mu_global    ~ normal(247, 80);
  sigma_global ~ normal(0, 50);
  sigma_cole   ~ normal(0, 50);
  sigma_ald    ~ normal(0, 30);

  mu_raw    ~ std_normal();
  alpha_raw ~ std_normal();

  horas_effect        ~ normal(0, 30);
  estrato_effect      ~ normal(0, 30);
  internet_effect     ~ normal(0, 30);
  educ_madre_effect   ~ normal(0, 30);
  educ_padre_effect   ~ normal(0, 30);
  naturaleza_effect   ~ normal(0, 30);

  // --- Variables latentes de la mezcla ---
  w ~ exponential(1.0 / sigma_ald);

  // --- Verosimilitud (mezcla normal-exponencial = ALD suavizada) ---
  for (i in 1:N) {
    real loc_i = alpha[cole[i]]
               + estrato_cen[estrato[i]]
               + internet_cen[internet[i]]
               + horas_cen[horas[i]]
               + educ_madre_cen[educ_madre[i]]
               + educ_padre_cen[educ_padre[i]]
               + naturaleza_cen[ naturaleza[ cole[i] ] ];

    y[i] ~ normal(loc_i + delta * w[i], sqrt(gamma_sq * w[i]));
  }
}
