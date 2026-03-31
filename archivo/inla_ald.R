###############################################################################
# Regresion cuantilica bayesiana jerarquica via INLA
# Aproximacion Gaussiana de la ALD (Kozumi & Kobayashi, 2011)
#
# La ALD se escribe como mezcla normal-exponencial:
#   y_i | w_i ~ N(mu_i + delta*w_i, sigma_ald * sqrt(tau_comp * w_i))
#   w_i       ~ Exp(1/sigma_ald)
#
# Algoritmo EM con INLA en el E-step:
#   1. Dado w_i: ajustar Gaussiano ponderado en INLA
#   2. Dado el ajuste: actualizar w_i desde la GIG posterior
#   3. Repetir hasta convergencia
#
# Ref: Yue & Rue (2011), Kozumi & Kobayashi (2011)
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

cat("\n========== Dimensiones ==========\n")
cat("  Estudiantes: ", N, "\n")
cat("  Colegios:    ", length(coles_unicos), "\n")
cat("  Municipios:  ", length(munis_unicos), "\n")

# =============================================================================
# 4. Tipo de contrastes
# =============================================================================
# "treatment" = nivel de referencia (default R, para interpretacion directa)
# "sum"       = sum-to-zero (para ANOVA bayesiana, sin nivel de referencia)
CONTRASTE <- "sum"   # <--- cambiar a "treatment" si se quiere nivel de referencia

if (CONTRASTE == "sum") {
  cat("\n>>> Usando contrastes SUM-TO-ZERO (ANOVA bayesiana)\n")
  cat(">>> Cada coeficiente es la desviacion del efecto de esa categoria\n")
  cat(">>> respecto a la media global. La ultima categoria se omite\n")
  cat(">>> pero se calcula como: -(suma de las demas).\n\n")
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
# 5. Parametros de la mezcla ALD
# =============================================================================
TAU <- 0.8

tau_comp <- TAU * (1 - TAU)                     # tau*(1-tau)
delta_unit <- (1 - 2 * TAU) / tau_comp          # (1-2*tau)/(tau*(1-tau))
gamma_unit <- 2 / tau_comp                      # 2/(tau*(1-tau))

cat("\nParametros ALD (tau =", TAU, "):\n")
cat("  delta_unit:", delta_unit, "\n")
cat("  gamma_unit:", gamma_unit, "\n")

# =============================================================================
# 5. Formula INLA (igual para todas las iteraciones)
# =============================================================================
# y_tilde_i ~ N(mu_i, sigma_ald * sqrt(gamma_unit * w_i))
# donde y_tilde_i = y_i - delta * w_i (respuesta ajustada)
# y las precisiones de INLA = 1 / (gamma_unit * w_i)

formula_mod <- y_tilde ~ 1 +
  estrato + internet + horas + educ_madre + educ_padre + naturaleza +
  f(muni_idx, model = "iid",
    hyper = list(prec = list(prior = "pc.prec", param = c(50, 0.01)))) +
  f(cole_idx, model = "iid",
    hyper = list(prec = list(prior = "pc.prec", param = c(50, 0.01))))

# =============================================================================
# 6. Algoritmo EM: INLA + actualizacion de w
# =============================================================================
MAX_ITER <- 50
TOL <- 0.01          # tolerancia en cambio relativo de sigma_ald
DAMP <- 0.3          # factor de amortiguamiento (0 = sin cambio, 1 = cambio total)
sigma_ald <- 10       # estimacion inicial conservadora

# Inicializar w_i = E[Exp(1)] = 1  (exponencial unitaria)
w     <- rep(1, N)       # E[W_i]: para la pseudo-respuesta y_tilde
inv_w <- rep(1, N)       # E[1/W_i]: para el scale de INLA (precision)

# Clipear para estabilidad numerica
clip_w <- function(w_vec) pmax(pmin(w_vec, 1e4), 1e-6)

cat("\n========== Iniciando EM con INLA ==========\n")

for (iter in 1:MAX_ITER) {
  t_start <- Sys.time()

  # --- E-step: ajustar INLA con w fijos ---
  delta <- sigma_ald * delta_unit
  gamma_sq <- sigma_ald^2 * gamma_unit

  # Respuesta ajustada (E[W]) y pesos de precision (E[1/W])
  datos$y_tilde <- datos$y - delta * w
  scale_vec <- clip_w(inv_w)

  fit <- inla(
    formula_mod,
    family = "gaussian",
    data = datos,
    scale = scale_vec,
    control.compute = list(config = FALSE),
    control.predictor = list(compute = TRUE, link = 1),
    control.inla = list(strategy = "simplified.laplace", int.strategy = "ccd"),
    num.threads = parallel::detectCores(),
    verbose = FALSE,
    safe = FALSE
  )

  # Fitted values
  mu_fitted <- fit$summary.linear.predictor[, "mean"]

  # Actualizar sigma_ald desde la precision estimada
  prec_gauss <- fit$summary.hyperpar[
    "Precision for the Gaussian observations", "mean"]
  sigma_ald_new <- sqrt(1 / (prec_gauss * gamma_unit))

  # --- M-step: actualizar w_i desde la posterior GIG(1/2, a, b) ---
  # E[w] = sqrt(b/a) * (1 + 1/sqrt(a*b))  (formula cerrada para p=1/2)
  resid <- datos$y - mu_fitted
  # Prior w ~ Exp(1): contribuye exp(-w), coeficiente en a = 2
  a_gig <- delta^2 / gamma_sq + 2
  b_gig <- resid^2 / gamma_sq

  # Evitar division por cero
  b_gig <- pmax(b_gig, 1e-10)

  # GIG(1/2, a, b): formulas cerradas
  sqrt_ab <- sqrt(a_gig * b_gig)

  # E[W] para la pseudo-respuesta (y_tilde = y - delta*E[W])
  e_w     <- sqrt(b_gig / a_gig) * (1 + 1 / pmax(sqrt_ab, 1e-10))
  e_w     <- clip_w(e_w)

  # E[1/W] para el scale de INLA (Q-function requiere E[1/W], no 1/E[W])
  # Para GIG(1/2, a, b): E[1/W] = sqrt(a/b)
  e_inv_w <- clip_w(sqrt(a_gig / b_gig))

  # Actualizaciones amortiguadas para evitar oscilaciones
  sigma_ald_damped <- DAMP * sigma_ald_new + (1 - DAMP) * sigma_ald
  w_damped         <- DAMP * e_w     + (1 - DAMP) * w
  inv_w_damped     <- DAMP * e_inv_w + (1 - DAMP) * inv_w

  # Convergencia: cambio relativo en sigma_ald
  rel_change <- abs(sigma_ald_damped - sigma_ald) / max(abs(sigma_ald), 1e-10)
  sigma_ald <- sigma_ald_damped

  # Actualizar w y inv_w
  w     <- w_damped
  inv_w <- inv_w_damped

  t_elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

  cat(sprintf("  Iter %2d: sigma_ald = %.3f, cambio_rel = %.5f (%.1fs)\n",
              iter, sigma_ald, rel_change, t_elapsed))

  if (rel_change < TOL && iter > 2) {
    cat(">>> Convergencia alcanzada!\n")
    break
  }
}

# =============================================================================
# 7. Ajuste final con estrategia completa (no simplificada)
# =============================================================================
cat("\n>>> Ajuste final con estrategia completa...\n")
delta_final <- sigma_ald * delta_unit
gamma_sq_final <- sigma_ald^2 * gamma_unit
datos$y_tilde <- datos$y - delta_final * w

# Recalcular E[1/W] con los parametros finales del EM
resid_final   <- datos$y - mu_fitted
a_gig_final   <- delta_final^2 / gamma_sq_final + 2
b_gig_final   <- pmax(resid_final^2 / gamma_sq_final, 1e-10)
scale_vec_final <- clip_w(sqrt(a_gig_final / b_gig_final))

fit_final <- inla(
  formula_mod,
  family = "gaussian",
  data = datos,
  scale = scale_vec_final,
  control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE, link = 1),
  num.threads = parallel::detectCores(),
  verbose = TRUE
)

# =============================================================================
# 8. Resultados
# =============================================================================
cat("\n========== Efectos fijos (cuantil", TAU, ") ==========\n")
cat("Contrastes:", CONTRASTE, "\n\n")
print(round(fit_final$summary.fixed, 3))

# Si usamos sum-to-zero, calcular la categoria omitida (ultima de cada factor)
if (CONTRASTE == "sum") {
  cat("\n========== Efectos completos sum-to-zero (incluyendo cat. omitida) ==========\n")
  cat("La ultima categoria de cada factor = -(suma de las demas)\n\n")

  # Funcion para reconstruir efectos completos
  reconstruir_stz <- function(fit, prefijo, niveles) {
    fijos <- fit$summary.fixed
    idx <- grep(paste0("^", prefijo), rownames(fijos))
    medias <- fijos[idx, "mean"]
    ultima <- -sum(medias)
    efectos <- c(medias, ultima)
    names(efectos) <- niveles
    return(round(efectos, 3))
  }

  cat("Estrato (1-6):\n")
  print(reconstruir_stz(fit_final, "estrato", levels(datos$estrato)))

  cat("\nInternet (No, Si):\n")
  print(reconstruir_stz(fit_final, "internet", levels(datos$internet)))

  cat("\nHoras (No_trabaja, Trabaja):\n")
  print(reconstruir_stz(fit_final, "horas", levels(datos$horas)))

  cat("\nEducacion madre:\n")
  print(reconstruir_stz(fit_final, "educ_madre", niveles_educ))

  cat("\nEducacion padre:\n")
  print(reconstruir_stz(fit_final, "educ_padre", niveles_educ))

  cat("\nNaturaleza:\n")
  print(reconstruir_stz(fit_final, "naturaleza", levels(datos$naturaleza)))
}

cat("\n========== Hiperparametros ==========\n")
print(round(fit_final$summary.hyperpar, 3))

cat("\n========== SD efectos aleatorios ==========\n")
hyper_names <- rownames(fit_final$summary.hyperpar)
for (h in hyper_names) {
  if (grepl("Precision", h)) {
    prec_mean <- fit_final$summary.hyperpar[h, "mean"]
    sd_est <- 1 / sqrt(prec_mean)
    cat(" ", h, "-> SD =", round(sd_est, 2), "\n")
  }
}

cat("\n  sigma_ald (escala ALD):", round(sigma_ald, 3), "\n")

# =============================================================================
# 9. Guardar
# =============================================================================
resultados <- list(
  fit = fit_final,
  sigma_ald = sigma_ald,
  w = w,
  tau = TAU,
  delta = sigma_ald * delta_unit,
  gamma_sq = sigma_ald^2 * gamma_unit,
  contraste = CONTRASTE
)
save(resultados, file = "inla_ald_modelo_completo.RData")

mapeo <- list(
  municipios = data.frame(idx = seq_along(munis_unicos), cod_mcpio = munis_unicos),
  colegios   = data.frame(idx = seq_along(coles_unicos), cod_icfes = coles_unicos)
)
save(mapeo, file = "mapeo_indices_inla.RData")

cat("\n>>> Resultados guardados en inla_ald_modelo_completo.RData\n")
cat(">>> Mapeo en mapeo_indices_inla.RData\n")