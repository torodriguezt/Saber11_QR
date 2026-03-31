###############################################################################
# Estudio de simulacion: regresion cuantilica bayesiana via EM-INLA
#
# Objetivo: verificar que el algoritmo EM con INLA recupera los parametros
# verdaderos (efectos fijos, varianzas de efectos aleatorios, sigma_ald)
# bajo datos generados desde el modelo ALD con mezcla normal-exponencial.
#
# Parametrizacion correcta (Kozumi & Kobayashi, 2011):
#   Y = mu + delta*W + sqrt(gamma_sq * W) * Z,  W ~ Exp(1), Z ~ N(0,1)
#   donde delta = sigma*(1-2p)/(p*(1-p)), gamma_sq = sigma^2 * 2/(p*(1-p))
#   => Y ~ ALD(mu, sigma, p)
#
# Diseno:
#   - n_sim replicas independientes
#   - Datos con estructura jerarquica (grupos) y covariables continuas/binarias
#   - Se evalua sesgo, RMSE e intervalos de cobertura al 95%
###############################################################################

library(INLA)
library(lme4)

# =============================================================================
# 1. Configuracion del estudio
# =============================================================================
set.seed(2024)

n_sim       <- 10       # numero de replicas (reducir a 5-10 para prueba rapida)
N           <- 2000     # observaciones por replica
n_grupos    <- 30       # numero de grupos (efecto aleatorio)
TAU_vec     <- c(0.25, 0.5, 0.75)  # cuantiles a evaluar

# Parametros verdaderos
beta_true      <- c(50, 3, -2, 5)  # intercepto, x1, x2, x3
sigma_u_true   <- 8                # SD del efecto aleatorio de grupo
sigma_ald_true <- 10               # escala ALD

# =============================================================================
# 2. Funcion generadora de datos ALD jerarquicos
# =============================================================================
#' Genera datos desde la mezcla ALD:
#'   Y_i = mu_i + delta*W_i + sqrt(gamma_sq * W_i) * Z_i
#'   W_i ~ Exp(1),  Z_i ~ N(0,1)
#'   mu_i = X_i * beta + u_{g[i]}
#'   u_j ~ N(0, sigma_u^2)
#'
#' Esto produce Y_i ~ ALD(mu_i, sigma_ald, tau)
generar_datos_ald <- function(N, n_grupos, beta, sigma_u, sigma_ald, tau) {

  tau_comp   <- tau * (1 - tau)
  delta_unit <- (1 - 2 * tau) / tau_comp
  gamma_unit <- 2 / tau_comp

  # Escalar por sigma_ald
  delta    <- sigma_ald * delta_unit
  gamma_sq <- sigma_ald^2 * gamma_unit

  # Covariables
  x1 <- rnorm(N, 0, 1)
  x2 <- rbinom(N, 1, 0.5)
  x3 <- rnorm(N, 2, 0.5)

  # Efecto aleatorio de grupo
  grupo <- sample(1:n_grupos, N, replace = TRUE)
  u     <- rnorm(n_grupos, 0, sigma_u)

  # Predictor lineal (cuantil condicional tau)
  mu <- beta[1] + beta[2] * x1 + beta[3] * x2 + beta[4] * x3 + u[grupo]

  # Mezcla ALD con W ~ Exp(1)
  w <- rexp(N, rate = 1)
  z <- rnorm(N)
  y <- mu + delta * w + sqrt(gamma_sq * w) * z

  data.frame(y = y, x1 = x1, x2 = x2, x3 = x3,
             grupo = grupo, mu_true = mu)
}

# =============================================================================
# 3. Funcion EM-INLA (version limpia para simulacion)
# =============================================================================
em_inla_qr <- function(datos, tau, max_iter = 50, tol = 0.005, damp = 0.5,
                       sigma_init = 5, verbose = FALSE) {

  N <- nrow(datos)
  tau_comp   <- tau * (1 - tau)
  delta_unit <- (1 - 2 * tau) / tau_comp
  gamma_unit <- 2 / tau_comp

  sigma_ald <- sigma_init
  w     <- rep(1, N)       # E[W_i]: para la pseudo-respuesta y_tilde
  inv_w <- rep(1, N)       # E[1/W_i]: para el scale de INLA (precision)
  clip_w <- function(v) pmax(pmin(v, 1e4), 1e-6)

  # --- Inicializacion de sigma_u via lme4 ---
  lmer_init <- tryCatch(
    lmer(y ~ x1 + x2 + x3 + (1 | grupo), data = datos, REML = TRUE),
    error = function(e) NULL
  )
  if (!is.null(lmer_init)) {
    vc <- as.data.frame(VarCorr(lmer_init))
    sigma_u_init <- vc$sdcor[vc$grp == "grupo"]
    if (verbose) cat("  Init lme4: sigma_u =", round(sigma_u_init, 2), "\n")
  } else {
    sigma_u_init <- 5
  }
  prec_u_fixed <- 1 / max(sigma_u_init, 0.1)^2

  # =====================================================================
  # FASE 1: EM con precision del efecto aleatorio FIJA (de lme4)
  # Esto permite que sigma_ald y w converjan sin que los w absorban
  # la variacion entre grupos (sigma_u no participa en la optimizacion).
  # =====================================================================
  formula_fixed_u <- y_tilde ~ 1 + x1 + x2 + x3 +
    f(grupo, model = "iid",
      hyper = list(prec = list(initial = log(prec_u_fixed), fixed = TRUE)))

  converged <- FALSE
  if (verbose) cat("  --- Fase 1: EM con sigma_u fija ---\n")

  for (iter in 1:max_iter) {
    delta    <- sigma_ald * delta_unit
    gamma_sq <- sigma_ald^2 * gamma_unit

    datos$y_tilde <- datos$y - delta * w        # usa E[W]
    scale_vec     <- clip_w(inv_w)               # usa E[1/W], NO 1/E[W]

    fit <- tryCatch(
      inla(formula_fixed_u,
           family    = "gaussian",
           data      = datos,
           scale     = scale_vec,
           control.compute   = list(config = FALSE),
           control.predictor = list(compute = TRUE, link = 1),
           control.inla      = list(strategy = "gaussian",
                                    int.strategy = "eb"),
           num.threads = "1:1",
           verbose = FALSE,
           safe = FALSE),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      if (verbose) cat("  Iter", iter, ": INLA fallo.\n")
      break
    }

    mu_fitted <- fit$summary.linear.predictor[, "mean"]

    prec_gauss    <- fit$summary.hyperpar[
      "Precision for the Gaussian observations", "mean"]
    sigma_ald_new <- sqrt(1 / (prec_gauss * gamma_unit))

    # M-step: actualizar w desde GIG(1/2, a, b)
    resid <- datos$y - mu_fitted
    a_gig <- delta^2 / gamma_sq + 2
    b_gig <- pmax(resid^2 / gamma_sq, 1e-10)

    sqrt_ab <- sqrt(a_gig * b_gig)

    # E[W] para la pseudo-respuesta (y_tilde = y - delta*E[W])
    e_w     <- sqrt(b_gig / a_gig) * (1 + 1 / pmax(sqrt_ab, 1e-10))
    e_w     <- clip_w(e_w)

    # E[1/W] para el scale de INLA (la Q-function requiere esto, no 1/E[W])
    # Para GIG(1/2, a, b): E[1/W] = sqrt(a/b)
    e_inv_w <- clip_w(sqrt(a_gig / b_gig))

    sigma_ald_damped <- damp * sigma_ald_new + (1 - damp) * sigma_ald
    w_damped         <- damp * e_w     + (1 - damp) * w
    inv_w_damped     <- damp * e_inv_w + (1 - damp) * inv_w

    rel_change <- abs(sigma_ald_damped - sigma_ald) / max(abs(sigma_ald), 1e-10)
    sigma_ald  <- sigma_ald_damped
    w          <- w_damped
    inv_w      <- inv_w_damped

    if (verbose) {
      cat(sprintf("  F1 Iter %2d: sigma_ald=%.3f, cambio=%.5f\n",
                  iter, sigma_ald, rel_change))
    }

    if (rel_change < tol && iter > 3) {
      converged <- TRUE
      break
    }
  }
  n_iter_f1 <- iter

  # =====================================================================
  # FASE 2: Ajuste final con sigma_u LIBRE
  # Con w estabilizados, INLA puede separar la variacion entre grupos
  # (sigma_u) del ruido ALD (capturado por w) sin riesgo de colapso.
  # =====================================================================
  if (verbose) cat("  --- Fase 2: ajuste final con sigma_u libre ---\n")

  delta    <- sigma_ald * delta_unit
  gamma_sq <- sigma_ald^2 * gamma_unit
  datos$y_tilde <- datos$y - delta * w

  # Recalcular E[1/W] con los parametros finales de Fase 1
  resid_final   <- datos$y - mu_fitted
  a_gig_final   <- delta^2 / gamma_sq + 2
  b_gig_final   <- pmax(resid_final^2 / gamma_sq, 1e-10)
  e_inv_w_final <- clip_w(sqrt(a_gig_final / b_gig_final))
  scale_vec     <- e_inv_w_final

  # PC prior centrada en la estimacion lme4
  pc_param <- c(2 * max(sigma_u_init, 1), 0.01)
  formula_free_u <- y_tilde ~ 1 + x1 + x2 + x3 +
    f(grupo, model = "iid",
      hyper = list(prec = list(prior = "pc.prec", param = pc_param)))

  fit_final <- tryCatch(
    inla(formula_free_u,
         family    = "gaussian",
         data      = datos,
         scale     = scale_vec,
         control.compute   = list(config = TRUE),
         control.predictor = list(compute = TRUE, link = 1),
         num.threads = "1:1",
         verbose = FALSE),
    error = function(e) NULL
  )

  if (is.null(fit_final)) return(NULL)

  # Extraer resultados
  beta_hat  <- fit_final$summary.fixed[, "mean"]
  beta_lo   <- fit_final$summary.fixed[, "0.025quant"]
  beta_hi   <- fit_final$summary.fixed[, "0.975quant"]

  prec_grupo <- fit_final$summary.hyperpar[
    grep("grupo", rownames(fit_final$summary.hyperpar)), "mean"]
  sigma_u_hat <- 1 / sqrt(prec_grupo)

  if (verbose) {
    cat(sprintf("  Final: sigma_ald=%.3f, sigma_u=%.3f\n",
                sigma_ald, sigma_u_hat))
  }

  list(beta_hat    = beta_hat,
       beta_lo     = beta_lo,
       beta_hi     = beta_hi,
       sigma_u_hat = sigma_u_hat,
       sigma_ald   = sigma_ald,
       converged   = converged,
       n_iter      = n_iter_f1,
       fit         = fit_final)
}

# =============================================================================
# 4. Verificacion rapida (1 replica, verbose)
# =============================================================================
cat("========== Verificacion rapida (1 replica, tau=0.5) ==========\n")
dat_test <- generar_datos_ald(N, n_grupos, beta_true,
                               sigma_u_true, sigma_ald_true, tau = 0.5)
fit_test <- em_inla_qr(dat_test, tau = 0.5, verbose = TRUE)

if (!is.null(fit_test)) {
  cat("\nResultados vs verdaderos:\n")
  cat(sprintf("  %-12s %8s %8s\n", "Param", "Verdad", "Estimado"))
  noms <- c("intercept", "x1", "x2", "x3")
  for (j in seq_along(beta_true)) {
    cat(sprintf("  %-12s %8.2f %8.2f\n", noms[j], beta_true[j], fit_test$beta_hat[j]))
  }
  cat(sprintf("  %-12s %8.2f %8.2f\n", "sigma_u", sigma_u_true, fit_test$sigma_u_hat))
  cat(sprintf("  %-12s %8.2f %8.2f\n", "sigma_ald", sigma_ald_true, fit_test$sigma_ald))
}

# =============================================================================
# 5. Ejecutar simulacion completa
# =============================================================================
resultados <- list()

for (tau in TAU_vec) {
  cat("\n============================================================\n")
  cat("  Cuantil tau =", tau, "\n")
  cat("============================================================\n")

  res_tau <- data.frame(
    sim       = integer(0),
    param     = character(0),
    verdadero = numeric(0),
    estimado  = numeric(0),
    lo95      = numeric(0),
    hi95      = numeric(0),
    stringsAsFactors = FALSE
  )

  n_exito <- 0
  n_conv  <- 0

  for (s in 1:n_sim) {
    t0 <- Sys.time()
    dat <- generar_datos_ald(N, n_grupos, beta_true,
                             sigma_u_true, sigma_ald_true, tau)

    fit <- em_inla_qr(dat, tau = tau, verbose = FALSE)
    dt  <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)

    if (is.null(fit)) {
      cat(sprintf("  [sim %3d/%d] FALLO (%.1fs)\n", s, n_sim, dt))
      next
    }

    n_exito <- n_exito + 1
    if (fit$converged) n_conv <- n_conv + 1

    # Guardar betas
    nombres_beta <- c("intercept", "x1", "x2", "x3")
    for (j in seq_along(beta_true)) {
      res_tau <- rbind(res_tau, data.frame(
        sim       = s,
        param     = nombres_beta[j],
        verdadero = beta_true[j],
        estimado  = fit$beta_hat[j],
        lo95      = fit$beta_lo[j],
        hi95      = fit$beta_hi[j]
      ))
    }

    # sigma_u
    res_tau <- rbind(res_tau, data.frame(
      sim       = s,
      param     = "sigma_u",
      verdadero = sigma_u_true,
      estimado  = fit$sigma_u_hat,
      lo95      = NA,
      hi95      = NA
    ))

    # sigma_ald
    res_tau <- rbind(res_tau, data.frame(
      sim       = s,
      param     = "sigma_ald",
      verdadero = sigma_ald_true,
      estimado  = fit$sigma_ald,
      lo95      = NA,
      hi95      = NA
    ))

    cat(sprintf("  [sim %3d/%d] conv=%s, iters=%d, sigma_ald=%.2f, sigma_u=%.2f (%.1fs)\n",
                s, n_sim, fit$converged, fit$n_iter, fit$sigma_ald,
                fit$sigma_u_hat, dt))
  }

  cat(sprintf("\n  Resumen: %d/%d exitosas, %d/%d convergieron\n",
              n_exito, n_sim, n_conv, n_exito))

  resultados[[as.character(tau)]] <- res_tau
}

# =============================================================================
# 6. Resumen de resultados
# =============================================================================
cat("\n\n################################################################\n")
cat("                RESUMEN DEL ESTUDIO DE SIMULACION\n")
cat("################################################################\n")

for (tau in TAU_vec) {
  res <- resultados[[as.character(tau)]]
  if (nrow(res) == 0) next

  cat(sprintf("\n--- tau = %.2f (n_sim efectivas = %d) ---\n",
              tau, length(unique(res$sim))))
  cat(sprintf("  %-12s %8s %8s %8s %8s %8s\n",
              "Parametro", "Verdad", "Media", "Sesgo", "RMSE", "Cob95%"))
  cat(paste(rep("-", 62), collapse = ""), "\n")

  for (p in unique(res$param)) {
    sub <- res[res$param == p, ]
    verdad  <- sub$verdadero[1]
    media   <- mean(sub$estimado, na.rm = TRUE)
    sesgo   <- media - verdad
    rmse    <- sqrt(mean((sub$estimado - verdad)^2, na.rm = TRUE))

    if (all(!is.na(sub$lo95))) {
      cobertura <- mean(sub$lo95 <= verdad & verdad <= sub$hi95, na.rm = TRUE)
      cob_str   <- sprintf("%.1f%%", cobertura * 100)
    } else {
      cob_str <- "  ---"
    }

    cat(sprintf("  %-12s %8.2f %8.2f %8.3f %8.3f %8s\n",
                p, verdad, media, sesgo, rmse, cob_str))
  }
}

# =============================================================================
# 7. Guardar resultados
# =============================================================================
save(resultados, beta_true, sigma_u_true, sigma_ald_true,
     file = "archivo/resultados_simulacion_inla_ald.RData")
cat("\n>>> Resultados guardados en archivo/resultados_simulacion_inla_ald.RData\n")

# =============================================================================
# 8. Graficos diagnosticos (requiere ggplot2)
# =============================================================================
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  todos <- do.call(rbind, lapply(names(resultados), function(tau) {
    df <- resultados[[tau]]
    df$tau <- as.numeric(tau)
    df
  }))

  # Sesgo por parametro y cuantil
  resumen <- do.call(rbind, lapply(split(todos, list(todos$param, todos$tau)), function(sub) {
    data.frame(
      param    = sub$param[1],
      tau      = sub$tau[1],
      verdad   = sub$verdadero[1],
      media    = mean(sub$estimado, na.rm = TRUE),
      sd_est   = sd(sub$estimado, na.rm = TRUE)
    )
  }))

  p1 <- ggplot(resumen, aes(x = factor(tau), y = media - verdad)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = (media - verdad) - 1.96 * sd_est,
                      ymax = (media - verdad) + 1.96 * sd_est),
                  width = 0.1) +
    facet_wrap(~ param, scales = "free_y") +
    labs(title = "Sesgo por parametro y cuantil (IC 95% de Monte Carlo)",
         x = "Cuantil (tau)", y = "Sesgo (estimado - verdadero)") +
    theme_minimal()

  ggsave("archivo/simulacion_sesgo.png", p1, width = 10, height = 6)
  cat(">>> Grafico guardado en archivo/simulacion_sesgo.png\n")

  # Boxplots de estimaciones
  p2 <- ggplot(todos, aes(x = factor(tau), y = estimado)) +
    geom_boxplot(fill = "lightblue", alpha = 0.7) +
    geom_hline(aes(yintercept = verdadero), color = "red",
               linetype = "dashed", linewidth = 0.8) +
    facet_wrap(~ param, scales = "free_y") +
    labs(title = "Distribucion de estimaciones por cuantil",
         subtitle = "Linea roja = valor verdadero",
         x = "Cuantil (tau)", y = "Estimacion") +
    theme_minimal()

  ggsave("archivo/simulacion_boxplots.png", p2, width = 10, height = 6)
  cat(">>> Grafico guardado en archivo/simulacion_boxplots.png\n")
}
