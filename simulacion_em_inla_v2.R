# ==============================================================================
# SCRIPT CORREGIDO v2: REGRESIÓN CUANTÍLICA (ALD) CON INLA Y ALGORITMO EM
# Un solo Dataset con error Normal, Múltiples Cuantiles
# El intercepto verdadero se desplaza: beta0 + qnorm(tau)*sigma_error
# ==============================================================================
library(INLA)

# ==============================================================================
# 1. SIMULACIÓN DE DATOS (UN SOLO CONJUNTO DE DATOS)
# ==============================================================================
set.seed(123)
N <- 3000
J_sch <- 50
J_mun <- 20

school_idx <- sample(1:J_sch, N, replace = TRUE)
muni_idx <- sample(1:J_mun, N, replace = TRUE)

x1 <- rnorm(N, 0, 1)
x2 <- rnorm(N, 0, 1)

beta0_base <- 250.0
beta1_true <- 10.0
beta2_true <- -5.0
sd_sch_true <- 15.0
sd_mun_true <- 8.0
sigma_error <- 15.0

u_sch <- rnorm(J_sch, 0, sd_sch_true)
u_mun <- rnorm(J_mun, 0, sd_mun_true)

# Error Normal: el cuantil tau-esimo del intercepto es beta0 + qnorm(tau)*sigma
error_sim <- rnorm(N, mean = 0, sd = sigma_error)
y <- beta0_base + beta1_true * x1 + beta2_true * x2 +
  u_sch[school_idx] + u_mun[muni_idx] + error_sim

data <- data.frame(y, x1, x2, school_idx, muni_idx)

taus <- c(0.10, 0.25, 0.50, 0.75, 0.90)
resultados_finales <- list()

cat("========== CONFIGURACION ==========\n")
cat("  N:           ", N, "\n")
cat("  J_sch:       ", J_sch, "\n")
cat("  J_mun:       ", J_mun, "\n")
cat("  beta0_base:  ", beta0_base, "\n")
cat("  beta1_true:  ", beta1_true, "\n")
cat("  beta2_true:  ", beta2_true, "\n")
cat("  sigma_error: ", sigma_error, "\n")
cat("  sd_sch_true: ", sd_sch_true, "\n")
cat("  sd_mun_true: ", sd_mun_true, "\n")
cat("  Cuantiles:   ", taus, "\n\n")

# ==============================================================================
# 2. BUCLE PRINCIPAL SOBRE CADA CUANTIL
# ==============================================================================
for (tau in taus) {

  cat(sprintf("\n===============================================================\n"))
  cat(sprintf("################# EJECUTANDO tau = %.2f #################\n", tau))
  cat(sprintf("===============================================================\n"))

  theta <- (1 - 2 * tau) / (tau * (1 - tau))
  kappa2 <- 2 / (tau * (1 - tau))

  # --- FASE 0: INICIALIZACIÓN ---
  formula_free <- y ~ x1 + x2 +
    f(school_idx, model = "iid", constr = TRUE) +
    f(muni_idx, model = "iid", constr = TRUE)

  fit0 <- inla(formula_free, data = data, family = "gaussian",
               control.compute = list(config = TRUE), verbose = FALSE)

  prec_gauss0 <- fit0$summary.hyperpar[1, "mean"]
  sigma_est <- sqrt(1 / prec_gauss0)

  prec_school <- 1 / (15^2)
  prec_muni <- 1 / (15^2)
  mu_est <- fit0$summary.fitted.values$mean

  cat(sprintf("[DEBUG Fase 0] Gauss_prec = %.4f | sigma inicial = %.3f\n",
              prec_gauss0, sigma_est))

  # --- ALGORITMO EM ---
  max_iter <- 80
  tol <- 1e-4
  diff <- 1

  for (iter in 1:max_iter) {

    gamma2 <- sigma_est * kappa2

    # PASO E: Calcular variables latentes esperadas
    res <- data$y - mu_est
    chi <- pmax((res^2) / gamma2, 1e-10)
    psi <- (theta^2) / gamma2 + (2 / sigma_est)

    inv_w <- pmax(sqrt(psi / chi), 1e-6)
    w <- sqrt(chi / psi) + (1 / psi)

    data$y_tilde <- data$y - (theta / inv_w)
    data$scale_vec <- inv_w / gamma2

    # PASO M: INLA con hiperparametros fijos
    formula_fixed <- y_tilde ~ x1 + x2 +
      f(school_idx, model = "iid",
        hyper = list(prec = list(initial = log(prec_school), fixed = TRUE)),
        constr = TRUE) +
      f(muni_idx, model = "iid",
        hyper = list(prec = list(initial = log(prec_muni), fixed = TRUE)),
        constr = TRUE)

    fit_em <- inla(formula_fixed, data = data, family = "gaussian",
                   scale = data$scale_vec,
                   control.family = list(
                     hyper = list(prec = list(initial = 0, fixed = TRUE))
                   ), verbose = FALSE)

    # Actualizar precisiones de efectos aleatorios
    new_var_sch <- mean(fit_em$summary.random$school_idx$mean^2 +
                        fit_em$summary.random$school_idx$sd^2)
    new_var_mun <- mean(fit_em$summary.random$muni_idx$mean^2 +
                        fit_em$summary.random$muni_idx$sd^2)

    prec_school <- 1 / (0.5 * (1/prec_school) + 0.5 * new_var_sch)
    prec_muni   <- 1 / (0.5 * (1/prec_muni)   + 0.5 * new_var_mun)

    # Residuos con datos ORIGINALES
    mu_est <- fit_em$summary.fitted.values$mean
    res_new <- data$y - mu_est

    # Actualización de sigma (derivada analítica)
    term1 <- (inv_w / (2 * kappa2)) * res_new^2
    term2 <- (theta / kappa2) * res_new
    term3 <- ((theta^2) / (2 * kappa2) + 1) * w

    sigma_new <- mean(term1 - term2 + term3)

    # Criterio de parada
    diff <- abs(sigma_new - sigma_est) / sigma_est
    sigma_est <- sigma_new

    cat(sprintf(" Iter %2d: sigma=%.3f SD_sch=%.2f SD_mun=%.2f rel=%.5f\n",
                iter, sigma_est, sqrt(1/prec_school), sqrt(1/prec_muni), diff))

    if (diff < tol) break
  }

  # --- AJUSTE FINAL ---
  cat("\n[INFO] Convergencia alcanzada. Ejecutando ajuste final...\n")
  fit_final <- inla(formula_fixed, data = data, family = "gaussian",
                    scale = data$scale_vec,
                    control.family = list(
                      hyper = list(prec = list(initial = 0, fixed = TRUE))
                    ),
                    control.compute = list(config = TRUE), verbose = FALSE)

  # Beta0 verdadero para este cuantil
  beta0_true_tau <- beta0_base + qnorm(tau, mean = 0, sd = sigma_error)

  resultados_finales[[as.character(tau)]] <- list(
    tau = tau,
    beta0_est  = fit_final$summary.fixed["(Intercept)", "mean"],
    beta0_true = beta0_true_tau,
    beta1 = fit_final$summary.fixed["x1", "mean"],
    beta2 = fit_final$summary.fixed["x2", "mean"],
    sigma = sigma_est,
    sd_sch = sqrt(1/prec_school),
    sd_mun = sqrt(1/prec_muni),
    iteraciones = iter
  )
}

# ==============================================================================
# 3. TABLA RESUMEN
# ==============================================================================
cat("\n\n")
cat("========================================================================================\n")
cat("                     TABLA RESUMEN DE SIMULACIÓN POR CUANTIL                            \n")
cat("========================================================================================\n")
cat(sprintf("%-5s | %10s %10s | %8s %8s | %8s | %8s %8s | %s\n",
            "tau", "B0_True", "B0_Est", "beta1", "beta2",
            "sigma", "SD_sch", "SD_mun", "Iters"))
cat(paste(rep("-", 88), collapse = ""), "\n")

for (res in resultados_finales) {
  cat(sprintf("%5.2f | %10.1f %10.1f | %8.1f %8.1f | %8.3f | %8.2f %8.2f | %4d\n",
              res$tau,
              res$beta0_true, res$beta0_est,
              res$beta1, res$beta2,
              res$sigma,
              res$sd_sch, res$sd_mun,
              res$iteraciones))
}
cat(paste(rep("-", 88), collapse = ""), "\n")
cat(sprintf("REF   | %10s %10s | %8.1f %8.1f | %8s | %8.1f %8.1f |\n",
            "varies", "", beta1_true, beta2_true,
            "", sd_sch_true, sd_mun_true))
cat("========================================================================================\n")

# ==============================================================================
# 4. GUARDAR
# ==============================================================================
save(resultados_finales, beta0_base, beta1_true, beta2_true,
     sigma_error, sd_sch_true, sd_mun_true,
     file = "archivo/resultados_simulacion_v2.RData")
cat("\n>>> Resultados: archivo/resultados_simulacion_v2.RData\n")
