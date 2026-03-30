###############################################################################
# Regresion cuantilica bayesiana jerarquica con R-INLA
#
# Usa el enfoque de Padellini & Rue (2018): link = "quantile"
# El predictor lineal modela DIRECTAMENTE el cuantil tau de la distribucion,
# no la media. Esto reemplaza la ALD por una reparametrizacion de la Gaussiana.
#
# Ref: https://arxiv.org/abs/1804.03714
#      https://github.com/tulliapadellini/INLA-quantreg
#
# Instalar R-INLA (no esta en CRAN):
#   install.packages("INLA",
#     repos = c(getOption("repos"),
#               INLA = "https://inla.r-inla-download.org/R/stable"),
#     dep = TRUE)
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

# Indices numericos para efectos aleatorios
munis_unicos <- sort(unique(datos$COLE_COD_MCPIO_UBICACION))
coles_unicos <- sort(unique(datos$COLE_CODIGO_ICFES))
datos$muni_idx <- as.integer(factor(datos$COLE_COD_MCPIO_UBICACION,
                                     levels = munis_unicos))
datos$cole_idx <- as.integer(factor(datos$COLE_CODIGO_ICFES,
                                     levels = coles_unicos))

datos$y <- as.numeric(datos$PUNT_GLOBAL)

cat("\n========== Dimensiones ==========\n")
cat("  Estudiantes: ", nrow(datos), "\n")
cat("  Colegios:    ", length(coles_unicos), "\n")
cat("  Municipios:  ", length(munis_unicos), "\n")

# =============================================================================
# 4. Cuantil a estimar
# =============================================================================
TAU <- 0.8

# =============================================================================
# 5. Formula del modelo
# =============================================================================
# Estructura jerarquica:
#   - f(muni_idx, model="iid") : intercepto aleatorio por municipio
#   - f(cole_idx, model="iid") : intercepto aleatorio por colegio
#   Cada colegio pertenece a un solo municipio (anidamiento implicito).
#
# Efectos fijos: estrato, internet, horas, educ_madre, educ_padre, naturaleza

formula_mod <- y ~ 1 +
  estrato + internet + horas + educ_madre + educ_padre + naturaleza +
  f(muni_idx, model = "iid",
    hyper = list(prec = list(prior = "pc.prec", param = c(50, 0.01)))) +
  f(cole_idx, model = "iid",
    hyper = list(prec = list(prior = "pc.prec", param = c(50, 0.01))))

# Priors PC (Penalized Complexity):
#   pc.prec con param = c(50, 0.01) -> P(sigma > 50) = 0.01

# =============================================================================
# 6. Ajuste: Gaussiana con link quantile (Padellini & Rue, 2018)
# =============================================================================
# El predictor lineal eta_i modela directamente Q_{tau}(y_i | x_i),
# el cuantil tau de y_i. INLA internamente mapea eta_i al parametro
# canonico (media) de la Gaussiana via la funcion cuantil inversa.

# Verificar version de INLA y links disponibles
cat("Version de INLA:", as.character(packageVersion("INLA")), "\n")

# Intentar con quantile link; si no existe, usar modelo gaussiano directo
# y derivar cuantiles desde la predictiva posterior
use_quantile_link <- tryCatch({
  # Test rapido con subset minimo para verificar si el link existe
  test_df <- datos[1:100, ]
  test_fit <- inla(
    y ~ 1,
    family = "gaussian",
    data = test_df,
    control.family = list(
      control.link = list(model = "quantile", quantile = TAU)
    )
  )
  TRUE
}, error = function(e) {
  cat(">>> link='quantile' no disponible para Gaussiana en esta version.\n")
  cat(">>> Usando modelo Gaussiano jerarquico + cuantiles desde la predictiva.\n\n")
  FALSE
})

if (use_quantile_link) {
  cat("\n>>> Ajustando con link='quantile', tau =", TAU, "...\n")
  fit <- inla(
    formula_mod,
    family = "gaussian",
    data = datos,
    control.family = list(
      control.link = list(model = "quantile", quantile = TAU)
    ),
    control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
    control.predictor = list(compute = TRUE),
    verbose = TRUE
  )
} else {
  # Modelo Gaussiano jerarquico estandar (link = identity, modela la media)
  # Los cuantiles se derivan en la seccion 9
  cat("\n>>> Ajustando modelo Gaussiano jerarquico (modela la media)...\n")
  cat(">>> El cuantil 0.8 se derivara de la predictiva posterior.\n\n")
  fit <- inla(
    formula_mod,
    family = "gaussian",
    data = datos,
    control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
    control.predictor = list(compute = TRUE),
    verbose = TRUE
  )
}

# =============================================================================
# 7. Resultados
# =============================================================================
cat("\n========== Efectos fijos (sobre el cuantil", TAU, ") ==========\n")
print(round(fit$summary.fixed, 3))

cat("\n========== Hiperparametros ==========\n")
print(round(fit$summary.hyperpar, 3))

# Varianzas de los efectos aleatorios
cat("\n========== SD efectos aleatorios ==========\n")
hyper_names <- rownames(fit$summary.hyperpar)
for (h in hyper_names) {
  if (grepl("Precision", h)) {
    prec_mean <- fit$summary.hyperpar[h, "mean"]
    sd_est <- 1 / sqrt(prec_mean)
    cat(" ", h, "-> SD =", round(sd_est, 2), "\n")
  }
}

# DIC y WAIC
cat("\n========== Criterios de informacion ==========\n")
cat("  DIC:  ", fit$dic$dic, "\n")
cat("  WAIC: ", fit$waic$waic, "\n")

# =============================================================================
# 8. Derivar cuantil 0.8 (si se uso modelo Gaussiano de media)
# =============================================================================
if (!use_quantile_link) {
  cat("\n========== Cuantil", TAU, "desde la predictiva posterior ==========\n")
  cat("Para un modelo Gaussiano: Q_tau(y_i) = mu_i + sigma * qnorm(tau)\n\n")

  # sigma residual = 1/sqrt(precision)
  prec_lik <- fit$summary.hyperpar["Precision for the Gaussian observations", "mean"]
  sigma_lik <- 1 / sqrt(prec_lik)
  ajuste_cuantil <- sigma_lik * qnorm(TAU)

  cat("  sigma residual estimada:", round(sigma_lik, 2), "\n")
  cat("  Ajuste para cuantil", TAU, ": +", round(ajuste_cuantil, 2), "puntos\n\n")

  # Efectos fijos: en el modelo de media, el efecto sobre el cuantil
  # es el MISMO que sobre la media (homocedasticidad Gaussiana).
  # Solo el intercepto cambia: intercepto_q80 = intercepto_media + sigma*qnorm(0.8)
  cat("  Interpretacion de efectos fijos:\n")
  cat("  - Los coeficientes de covariables son los MISMOS para media y cuantil 0.8\n")
  cat("    (bajo homocedasticidad Gaussiana, un efecto que sube la media\n")
  cat("     sube todos los cuantiles en la misma magnitud)\n")
  cat("  - Solo el intercepto se ajusta: intercepto_q80 =",
      round(fit$summary.fixed["(Intercept)", "mean"] + ajuste_cuantil, 2), "\n")

  # Efectos aleatorios sobre el cuantil: tambien son los mismos
  # alpha_j en media = alpha_j en cuantil (misma varianza entre colegios)
  cat("\n  Efectos aleatorios:\n")
  cat("  - SD entre municipios y entre colegios son los mismos para\n")
  cat("    media y cuantil 0.8 (homocedasticidad)\n")
  cat("  - Esto es una LIMITACION: si la varianza entre colegios\n")
  cat("    es diferente en la cola superior (cuantil 0.8) que en la media,\n")
  cat("    el modelo Gaussiano no lo captura.\n")
  cat("  - El modelo Stan con ALD SI permite heterocedasticidad implicita.\n")
}

# =============================================================================
# 9. Guardar
# =============================================================================
save(fit, file = "inla_modelo_completo.RData")

mapeo <- list(
  municipios = data.frame(idx = seq_along(munis_unicos), cod_mcpio = munis_unicos),
  colegios   = data.frame(idx = seq_along(coles_unicos), cod_icfes = coles_unicos)
)
save(mapeo, file = "mapeo_indices_inla.RData")

cat("\n>>> Modelo guardado en inla_modelo_completo.RData\n")
cat(">>> Mapeo en mapeo_indices_inla.RData\n")
