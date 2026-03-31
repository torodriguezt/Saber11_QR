###############################################################################
# Regresion cuantilica bayesiana jerarquica via EM-INLA (v2)
# Parametrizacion v_i ~ Exp(sigma) (Kozumi & Kobayashi, 2011, Sec 3.1)
#
# Algoritmo EM con hiperparametros fijos en INLA:
#   a) INLA resuelve regresion ponderada (prec fija, RE fijos)
#   b) Actualizar v_i desde GIG(1/2, chi, psi)
#   c) Actualizar sigma desde derivada analitica del Q-function (factor 2/3)
#   d) Actualizar precisiones de RE desde E[alpha^2 + Var(alpha)]
###############################################################################

library(INLA)

# =============================================================================
# 1. Carga de datos
# =============================================================================
cat("Leyendo dataset completo SB11_20232.TXT ...\n")
datos <- as.data.frame(readr::read_delim("SB11_20232.TXT", delim = "\u00AC",
                                         quote = "", show_col_types = FALSE))
cat("Filas leidas:", nrow(datos), "\n")

# =============================================================================
# 2. Limpieza y filtrado
# =============================================================================
vars_requeridas <- c("FAMI_TIENEINTERNET", "FAMI_ESTRATOVIVIENDA",
                     "ESTU_HORASSEMANATRABAJA", "FAMI_EDUCACIONMADRE",
                     "FAMI_EDUCACIONPADRE", "PUNT_GLOBAL",
                     "COLE_CODIGO_ICFES", "COLE_COD_MCPIO_UBICACION",
                     "COLE_NATURALEZA")

for (v in vars_requeridas) {
  datos <- datos[!is.na(datos[[v]]) & datos[[v]] != "", ]
}

datos <- datos[!datos$FAMI_EDUCACIONMADRE %in% c("No sabe", "No Aplica") &
               !datos$FAMI_EDUCACIONPADRE %in% c("No sabe", "No Aplica"), ]
datos <- datos[datos$FAMI_ESTRATOVIVIENDA != "Sin Estrato", ]
cat("Datos limpios:", nrow(datos), "\n")

# =============================================================================
# 3. Recodificaciones
# =============================================================================
datos$internet <- factor(ifelse(datos$FAMI_TIENEINTERNET == "Si", "Si", "No"))

datos$estrato <- factor(ifelse(datos$FAMI_ESTRATOVIVIENDA == "Estrato 1", 1L,
                 ifelse(datos$FAMI_ESTRATOVIVIENDA == "Estrato 2", 2L,
                 ifelse(datos$FAMI_ESTRATOVIVIENDA == "Estrato 3", 3L,
                 ifelse(datos$FAMI_ESTRATOVIVIENDA == "Estrato 4", 4L,
                 ifelse(datos$FAMI_ESTRATOVIVIENDA == "Estrato 5", 5L, 6L))))))

datos$horas <- factor(ifelse(datos$ESTU_HORASSEMANATRABAJA == "0",
                              "No_trabaja", "Trabaja"))

agrupar_educacion <- function(nivel) {
  nivel <- trimws(nivel)
  if (nivel == "Ninguno") return("Ninguno")
  else if (nivel %in% c("Primaria completa", "Primaria incompleta")) return("Primaria")
  else if (nivel %in% c("Secundaria (Bachillerato) completa",
                         "Secundaria (Bachillerato) incompleta")) return("Secundaria")
  else if (nivel %in% c("T\u00e9cnica o tecnol\u00f3gica completa",
                         "T\u00e9cnica o tecnol\u00f3gica incompleta")) return("Tecnica")
  else if (nivel %in% c("Educaci\u00f3n profesional completa",
                         "Educaci\u00f3n profesional incompleta")) return("Profesional")
  else if (nivel == "Postgrado") return("Postgrado")
  else return(NA_character_)
}

niveles_educ <- c("Ninguno", "Primaria", "Secundaria",
                  "Tecnica", "Profesional", "Postgrado")

datos$educ_madre <- factor(sapply(datos$FAMI_EDUCACIONMADRE, agrupar_educacion),
                           levels = niveles_educ)
datos$educ_padre <- factor(sapply(datos$FAMI_EDUCACIONPADRE, agrupar_educacion),
                           levels = niveles_educ)

datos$naturaleza <- factor(datos$COLE_NATURALEZA)

datos <- datos[complete.cases(datos[, c("internet", "estrato", "horas",
                                        "educ_madre", "educ_padre",
                                        "naturaleza", "PUNT_GLOBAL")]), ]

# Indices para efectos aleatorios
munis_unicos <- sort(unique(datos$COLE_COD_MCPIO_UBICACION))
coles_unicos <- sort(unique(datos$COLE_CODIGO_ICFES))
datos$muni_idx <- as.integer(factor(datos$COLE_COD_MCPIO_UBICACION,
                                     levels = munis_unicos))
datos$cole_idx <- as.integer(factor(datos$COLE_CODIGO_ICFES,
                                     levels = coles_unicos))

datos$y <- as.numeric(datos$PUNT_GLOBAL)
N <- nrow(datos)
n_coles <- length(coles_unicos)
n_munis <- length(munis_unicos)

cat("\n========== Dimensiones ==========\n")
cat("  Estudiantes: ", N, "\n")
cat("  Colegios:    ", n_coles, "\n")
cat("  Municipios:  ", n_munis, "\n")

# =============================================================================
# 4. Tipo de contrastes
# =============================================================================
CONTRASTE <- "sum"

if (CONTRASTE == "sum") {
  cat("\n>>> Usando contrastes SUM-TO-ZERO (ANOVA bayesiana)\n")
  contrasts(datos$estrato)    <- contr.sum(levels(datos$estrato))
  contrasts(datos$internet)   <- contr.sum(levels(datos$internet))
  contrasts(datos$horas)      <- contr.sum(levels(datos$horas))
  contrasts(datos$educ_madre) <- contr.sum(levels(datos$educ_madre))
  contrasts(datos$educ_padre) <- contr.sum(levels(datos$educ_padre))
  contrasts(datos$naturaleza) <- contr.sum(levels(datos$naturaleza))
} else {
  cat("\n>>> Usando contrastes TREATMENT (nivel de referencia)\n")
}

# =============================================================================
# 5. Cuantiles a estimar
# =============================================================================
TAUS <- c(0.1, 0.25, 0.5, 0.75, 0.9)

# =============================================================================
# 6. Algoritmo EM-INLA (hiperparametros fijos)
# =============================================================================
resultados_todos <- list()

for (TAU in TAUS) {

  cat(sprintf("\n\n############### tau = %.2f ###############\n", TAU))

  theta  <- (1 - 2 * TAU) / (TAU * (1 - TAU))
  kappa2 <- 2 / (TAU * (1 - TAU))

  cat(sprintf("  theta = %.4f, kappa2 = %.4f\n", theta, kappa2))

  # --- Fase 0: ajuste libre para inicializar ---
  datos$y_tilde <- datos$y

  formula_free <- y_tilde ~ 1 +
    estrato + internet + horas + educ_madre + educ_padre + naturaleza +
    f(cole_idx, model = "iid", constr = TRUE) +
    f(muni_idx, model = "iid", constr = TRUE)

  fit0 <- inla(formula_free, data = datos, family = "gaussian",
               control.compute = list(config = TRUE),
               num.threads = parallel::detectCores(),
               verbose = FALSE)

  hyp0 <- fit0$summary.hyperpar
  sigma_est  <- sqrt(1 / hyp0[1, "mean"])
  mu_est     <- fit0$summary.fitted.values$mean

  # Inicializar precisiones de RE con valores conservadores (SD=15)
  # NO usar hyp0 porque INLA libre puede colapsar los RE a 0
  # (especialmente con la pseudo-respuesta y pesos del EM)
  prec_cole  <- 1 / (15^2)
  prec_muni  <- 1 / (15^2)

  cat(sprintf("  Fase 0: sigma=%.3f, SD_cole=%.2f (init), SD_muni=%.2f (init)\n",
              sigma_est, 1/sqrt(prec_cole), 1/sqrt(prec_muni)))

  # --- Loop EM ---
  MAX_ITER <- 80
  TOL <- 1e-4

  for (iter in 1:MAX_ITER) {
    t_start <- Sys.time()

    gamma2 <- sigma_est * kappa2

    # PASO E: variables latentes GIG(1/2, chi, psi)
    res <- datos$y - mu_est
    chi <- pmax((res^2) / gamma2, 1e-10)
    psi <- (theta^2) / gamma2 + (2 / sigma_est)

    inv_w <- pmax(sqrt(psi / chi), 1e-6)
    w     <- sqrt(chi / psi) + (1 / psi)

    # Pseudo-respuesta y pesos
    datos$y_tilde  <- datos$y - (theta / inv_w)
    scale_vec      <- inv_w / gamma2

    # PASO M: INLA con todos los hiperparametros fijos
    formula_fixed <- eval(bquote(
      y_tilde ~ 1 +
        estrato + internet + horas + educ_madre + educ_padre + naturaleza +
        f(cole_idx, model = "iid",
          hyper = list(prec = list(initial = .(log(prec_cole)), fixed = TRUE)),
          constr = TRUE) +
        f(muni_idx, model = "iid",
          hyper = list(prec = list(initial = .(log(prec_muni)), fixed = TRUE)),
          constr = TRUE)
    ))

    fit_em <- tryCatch(
      inla(formula_fixed, data = datos, family = "gaussian",
           scale = scale_vec,
           control.family = list(
             hyper = list(prec = list(initial = 0, fixed = TRUE))
           ),
           control.predictor = list(compute = TRUE, link = 1),
           num.threads = parallel::detectCores(),
           verbose = FALSE, safe = FALSE),
      error = function(e) {
        cat("  [!] INLA error:", e$message, "\n")
        return(NULL)
      }
    )

    if (is.null(fit_em)) break

    # Actualizar precisiones de efectos aleatorios
    cole_mean     <- fit_em$summary.random$cole_idx$mean
    cole_var_post <- fit_em$summary.random$cole_idx$sd^2
    new_var_cole  <- mean(cole_mean^2 + cole_var_post)

    muni_mean     <- fit_em$summary.random$muni_idx$mean
    muni_var_post <- fit_em$summary.random$muni_idx$sd^2
    new_var_muni  <- mean(muni_mean^2 + muni_var_post)

    prec_cole <- 1 / (0.5 * (1/prec_cole) + 0.5 * new_var_cole)
    prec_muni <- 1 / (0.5 * (1/prec_muni) + 0.5 * new_var_muni)

    # Residuos con datos originales
    mu_est  <- fit_em$summary.fitted.values$mean
    res_new <- datos$y - mu_est

    # Actualizar sigma (factor 2/3 del Q-function)
    term1 <- (inv_w / (2 * kappa2)) * res_new^2
    term2 <- (theta / kappa2) * res_new
    term3 <- ((theta^2) / (2 * kappa2) + 1) * w

    sigma_new <- (2 / 3) * mean(term1 - term2 + term3)

    # Convergencia
    rel_change <- abs(sigma_new - sigma_est) / max(abs(sigma_est), 1e-10)
    sigma_est  <- sigma_new

    t_elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

    cat(sprintf("  Iter %2d: sigma=%.3f SD_cole=%.2f SD_muni=%.2f rel=%.5f (%.1fs)\n",
                iter, sigma_est, sqrt(1/prec_cole), sqrt(1/prec_muni),
                rel_change, t_elapsed))

    if (rel_change < TOL && iter > 2) {
      cat("  >>> Convergencia alcanzada!\n")
      break
    }
  }

  # --- Ajuste final ---
  cat("  >>> Ajuste final...\n")

  formula_final <- eval(bquote(
    y_tilde ~ 1 +
      estrato + internet + horas + educ_madre + educ_padre + naturaleza +
      f(cole_idx, model = "iid",
        hyper = list(prec = list(initial = .(log(prec_cole)), fixed = TRUE)),
        constr = TRUE) +
      f(muni_idx, model = "iid",
        hyper = list(prec = list(initial = .(log(prec_muni)), fixed = TRUE)),
        constr = TRUE)
  ))

  fit_final <- inla(formula_final, data = datos, family = "gaussian",
                    scale = scale_vec,
                    control.family = list(
                      hyper = list(prec = list(initial = 0, fixed = TRUE))
                    ),
                    control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
                    control.predictor = list(compute = TRUE, link = 1),
                    num.threads = parallel::detectCores(),
                    verbose = FALSE)

  # --- Resultados de este cuantil ---
  cat(sprintf("\n========== Efectos fijos (cuantil %.2f) ==========\n", TAU))
  cat("Contrastes:", CONTRASTE, "\n\n")
  print(round(fit_final$summary.fixed, 3))

  if (CONTRASTE == "sum") {
    cat("\n--- Efectos completos sum-to-zero (incluyendo cat. omitida) ---\n")

    reconstruir_stz <- function(fit, prefijo, niveles) {
      fijos <- fit$summary.fixed
      idx <- grep(paste0("^", prefijo), rownames(fijos))
      medias <- fijos[idx, "mean"]
      ultima <- -sum(medias)
      efectos <- c(medias, ultima)
      names(efectos) <- niveles
      return(round(efectos, 3))
    }

    cat("Estrato:\n")
    print(reconstruir_stz(fit_final, "estrato", levels(datos$estrato)))
    cat("\nInternet:\n")
    print(reconstruir_stz(fit_final, "internet", levels(datos$internet)))
    cat("\nHoras:\n")
    print(reconstruir_stz(fit_final, "horas", levels(datos$horas)))
    cat("\nEducacion madre:\n")
    print(reconstruir_stz(fit_final, "educ_madre", niveles_educ))
    cat("\nEducacion padre:\n")
    print(reconstruir_stz(fit_final, "educ_padre", niveles_educ))
    cat("\nNaturaleza:\n")
    print(reconstruir_stz(fit_final, "naturaleza", levels(datos$naturaleza)))
  }

  cat(sprintf("\n  sigma (escala ALD): %.3f\n", sigma_est))
  cat(sprintf("  SD_cole: %.2f\n", sqrt(1/prec_cole)))
  cat(sprintf("  SD_muni: %.2f\n", sqrt(1/prec_muni)))

  # Guardar
  resultados_todos[[as.character(TAU)]] <- list(
    fit        = fit_final,
    tau        = TAU,
    sigma      = sigma_est,
    prec_cole  = prec_cole,
    prec_muni  = prec_muni,
    sd_cole    = sqrt(1/prec_cole),
    sd_muni    = sqrt(1/prec_muni),
    contraste  = CONTRASTE
  )
}

# =============================================================================
# 7. Resumen comparativo entre cuantiles
# =============================================================================
cat("\n\n")
cat("==========================================================================\n")
cat("             RESUMEN COMPARATIVO ENTRE CUANTILES                          \n")
cat("==========================================================================\n")
cat(sprintf("%-5s | %10s | %8s %8s | %s\n",
            "tau", "sigma", "SD_cole", "SD_muni", "Intercepto"))
cat(paste(rep("-", 60), collapse = ""), "\n")

for (r in resultados_todos) {
  b0 <- r$fit$summary.fixed["(Intercept)", "mean"]
  cat(sprintf("%5.2f | %10.3f | %8.2f %8.2f | %8.2f\n",
              r$tau, r$sigma, r$sd_cole, r$sd_muni, b0))
}
cat("==========================================================================\n")

# =============================================================================
# 8. Guardar
# =============================================================================
mapeo <- list(
  municipios = data.frame(idx = seq_along(munis_unicos), cod_mcpio = munis_unicos),
  colegios   = data.frame(idx = seq_along(coles_unicos), cod_icfes = coles_unicos)
)

save(resultados_todos, file = "inla_ald_v2_todos_cuantiles.RData")
save(mapeo, file = "mapeo_indices_inla.RData")

cat("\n>>> Resultados: inla_ald_v2_todos_cuantiles.RData\n")
cat(">>> Mapeo: mapeo_indices_inla.RData\n")
