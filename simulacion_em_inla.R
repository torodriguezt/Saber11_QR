# ==============================================================================
# SCRIPT CORREGIDO: REGRESIÓN CUANTÍLICA (ALD) CON INLA Y ALGORITMO EM
# ==============================================================================
library(INLA)

# 1. SIMULACIÓN DE DATOS 
set.seed(123)
N <- 3000
tau <- 0.10
J_sch <- 50
J_mun <- 20

school_idx <- sample(1:J_sch, N, replace = TRUE)
muni_idx <- sample(1:J_mun, N, replace = TRUE)

x1 <- rnorm(N, 0, 1)
x2 <- rnorm(N, 0, 1)

beta0_true <- 250.0
beta1_true <- 10.0
beta2_true <- -5.0
sigma_true <- 15.0
sd_sch_true <- 15.0
sd_mun_true <- 8.0

u_sch <- rnorm(J_sch, 0, sd_sch_true)
u_mun <- rnorm(J_mun, 0, sd_mun_true)

theta <- (1 - 2 * tau) / (tau * (1 - tau))
kappa2 <- 2 / (tau * (1 - tau))

z <- rexp(N, rate = 1 / sigma_true)
error_ald <- theta * z + sqrt(sigma_true * kappa2 * z) * rnorm(N)
y <- beta0_true + beta1_true * x1 + beta2_true * x2 + u_sch[school_idx] + u_mun[muni_idx] + error_ald

data <- data.frame(y, x1, x2, school_idx, muni_idx)

# ==============================================================================
# FASE 0: INICIALIZACIÓN
# ==============================================================================
cat(sprintf("\n############### tau = %.2f ###############\n", tau))
cat(sprintf("N = %d \n", N))

formula_free <- y ~ x1 + x2 + 
  f(school_idx, model = "iid", constr = TRUE) + 
  f(muni_idx, model = "iid", constr = TRUE)

fit0 <- inla(formula_free, data = data, family = "gaussian", 
             control.compute = list(config = TRUE))

hyp0 <- fit0$summary.hyperpar
prec_gauss0 <- hyp0[1, "mean"]
sigma_est <- sqrt(1 / prec_gauss0)

prec_school <- 1 / (15^2)
prec_muni <- 1 / (15^2)

cat(sprintf("[DEBUG Fase 0] Gauss_prec = %.4f | sigma inicial = %.3f\n", prec_gauss0, sigma_est))

mu_est <- fit0$summary.fitted.values$mean

# ==============================================================================
# ALGORITMO EM
# ==============================================================================
max_iter <- 80
tol <- 1e-4
diff <- 1

for (iter in 1:max_iter) {
  
  gamma2 <- sigma_est * kappa2
  
  # 2. PASO E: Calcular variables latentes esperadas
  res <- data$y - mu_est
  
  chi <- pmax((res^2) / gamma2, 1e-10) 
  psi <- (theta^2) / gamma2 + (2 / sigma_est)
  
  inv_w <- sqrt(psi / chi)
  w <- sqrt(chi / psi) + (1 / psi)
  
  # [CORRECCIÓN 1]: Proteger contra valores atípicos que generan 1/0 en inv_w
  inv_w <- pmax(inv_w, 1e-6)
  
  # [CORRECCIÓN 2]: La pseudo-respuesta matemáticamente exacta de la función Q
  data$y_tilde <- data$y - (theta / inv_w)
  data$scale_vec <- inv_w / gamma2
  
  # 3. PASO M: Maximización con INLA fijando las varianzas previas
  formula_fixed <- y_tilde ~ x1 + x2 + 
    f(school_idx, model = "iid", 
      hyper = list(prec = list(initial = log(prec_school), fixed = TRUE)),
      constr = TRUE) + 
    f(muni_idx, model = "iid", 
      hyper = list(prec = list(initial = log(prec_muni), fixed = TRUE)),
      constr = TRUE)
  
  fit_em <- inla(formula_fixed, data = data, family = "gaussian", 
                 scale = data$scale_vec,
                 control.family = list(hyper = list(prec = list(initial = 0, fixed = TRUE))))
  
  # 4. Actualizar parámetros
  sch_mean <- fit_em$summary.random$school_idx$mean
  sch_var_post <- fit_em$summary.random$school_idx$sd^2
  new_var_sch <- mean(sch_mean^2 + sch_var_post)
  
  mun_mean <- fit_em$summary.random$muni_idx$mean
  mun_var_post <- fit_em$summary.random$muni_idx$sd^2
  new_var_mun <- mean(mun_mean^2 + mun_var_post)
  
  prec_school <- 1 / (0.5 * (1/prec_school) + 0.5 * new_var_sch)
  prec_muni   <- 1 / (0.5 * (1/prec_muni)   + 0.5 * new_var_mun)
  
  # Obtener residuos usando los datos ORIGINALES, no la pseudo-respuesta
  mu_est <- fit_em$summary.fitted.values$mean
  res_new <- data$y - mu_est
  
  # [CORRECCIÓN 3]: Actualización de sigma usando la derivada analítica EXACTA del Paso M
  term1 <- (inv_w / (2 * kappa2)) * res_new^2
  term2 <- (theta / kappa2) * res_new
  term3 <- ( (theta^2) / (2 * kappa2) + 1 ) * w
  
  sigma_new <- mean(term1 - term2 + term3)
  
  # Criterio de parada
  diff <- abs(sigma_new - sigma_est) / sigma_est
  sigma_est <- sigma_new
  
  cat(sprintf(" Iter %2d: sigma=%.3f SD_sch=%.2f SD_mun=%.2f rel=%.5f\n", 
              iter, sigma_est, sqrt(1/prec_school), sqrt(1/prec_muni), diff))
  
  if (diff < tol) break
}

# ==============================================================================
# AJUSTE FINAL 
# ==============================================================================
cat("\n[INFO] Convergencia alcanzada. Ejecutando ajuste final...\n")

formula_final <- y_tilde ~ x1 + x2 + 
  f(school_idx, model = "iid", 
    hyper = list(prec = list(initial = log(prec_school), fixed = TRUE)),
    constr = TRUE) + 
  f(muni_idx, model = "iid", 
    hyper = list(prec = list(initial = log(prec_muni), fixed = TRUE)),
    constr = TRUE)

fit_final <- inla(formula_final, data = data, family = "gaussian", 
                  scale = data$scale_vec,
                  control.family = list(hyper = list(prec = list(initial = 0, fixed = TRUE))),
                  control.compute = list(config = TRUE))

# Resultados
cat(sprintf("\n --- Resultados tau = %.2f ---\n", tau))
cat(sprintf(" beta0:    True= %5.1f | Est= %5.1f [%5.1f, %5.1f]\n", beta0_true, fit_final$summary.fixed["(Intercept)", "mean"], fit_final$summary.fixed["(Intercept)", "0.025quant"], fit_final$summary.fixed["(Intercept)", "0.975quant"]))
cat(sprintf(" beta1:    True= %5.1f | Est= %5.1f [%5.1f, %5.1f]\n", beta1_true, fit_final$summary.fixed["x1", "mean"], fit_final$summary.fixed["x1", "0.025quant"], fit_final$summary.fixed["x1", "0.975quant"]))
cat(sprintf(" beta2:    True= %5.1f | Est= %5.1f [%5.1f, %5.1f]\n", beta2_true, fit_final$summary.fixed["x2", "mean"], fit_final$summary.fixed["x2", "0.025quant"], fit_final$summary.fixed["x2", "0.975quant"]))
cat(sprintf(" sigma:    True= %5.1f | Est= %.3f\n", sigma_true, sigma_est))
cat(sprintf(" SD_sch:   True= %5.1f | Est= %.2f\n", sd_sch_true, sqrt(1/prec_school)))
cat(sprintf(" SD_muni:  True= %5.1f | Est= %.2f\n", sd_mun_true, sqrt(1/prec_muni)))