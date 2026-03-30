data {
  int<lower=1> N;
  int<lower=1> J;
  int<lower=1> K;
  int<lower=1> H;
  int<lower=1> E;
  int<lower=1> N2;
  int<lower=1> I;
  int<lower=1> L;                       // <── NUEVO (2 niveles)
  int<lower=1> M;
  int<lower=1> P;

  // naturaleza ahora es un índice 1..L para cada colegio
  int<lower=1, upper=L> naturaleza[J];  // <── cambiado

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
  vector[J]   alpha;
  real<lower=0> sigma_cole;
  vector[K]   mu;
  real        mu_global;
  real<lower=0> sigma;
  real<lower=0> sigma_global;

  //vector[H]   horas_effect;
  //vector[E]   estrato_effect;
  //vector[I]   internet_effect;

  //vector[M]   educ_madre_effect;
  //vector[P]   educ_padre_effect;

  //vector[L]   naturaleza_effect;        // <── reemplaza beta_naturaleza
}

transformed parameters {
  //vector[H] horas_cen       = horas_effect       - mean(horas_effect);
  //vector[E] estrato_cen     = estrato_effect     - mean(estrato_effect);
  //vector[I] internet_cen    = internet_effect    - mean(internet_effect);
  //vector[M] educ_madre_cen  = educ_madre_effect  - mean(educ_madre_effect);
  //vector[P] educ_padre_cen  = educ_padre_effect  - mean(educ_padre_effect);

  // nuevo centrado para naturaleza
  //vector[L] naturaleza_cen  = naturaleza_effect  - mean(naturaleza_effect);
}

model {
  // Priors generales
  mu_global     ~ normal(247, 82);
  sigma_global  ~ inv_gamma(0.001, 0.001);
  mu            ~ normal(mu_global, sigma_global);
  sigma         ~ inv_gamma(0.001, 0.001);
  sigma_cole    ~ cauchy(0, 10);

  // Priors covariables
  //horas_effect        ~ normal(0, 1000);
  //estrato_effect      ~ normal(0, 1000);
  //internet_effect     ~ normal(0, 1000);
  //educ_madre_effect   ~ normal(0, 1000);
  //educ_padre_effect   ~ normal(0, 1000);
  //naturaleza_effect   ~ normal(0, 1000);     // <── nuevo

  // Jerarquía colegios
  for (j in 1:N2)
    alpha[col2[j]] ~ normal(mu[muni2[j]], sigma_cole);

  // Verosimilitud
  for (i in 1:N)
    y[i] ~ skew_double_exponential(
              alpha[cole[i]],
            2 * sigma, tau);
}
