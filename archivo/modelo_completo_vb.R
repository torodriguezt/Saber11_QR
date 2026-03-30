###############################################################################
# Modelo completo - Variational Bayes (ADVI + Pathfinder)
# Dataset COMPLETO: SB11_20232.TXT (~551K estudiantes)
#
# Usa el modelo reparametrizado (non-centered) en modelo_completo_vb.stan
# No necesita puntos iniciales, clusterizacion, ni particiones.
#
# Dos metodos disponibles:
#   1. rstan::vb()         - ADVI, rapido y simple
#   2. cmdstanr Pathfinder - mas robusto, recomendado
###############################################################################

# =============================================================================
# 0. Configuracion
# =============================================================================
USE_CMDSTANR <- FALSE   # TRUE = cmdstanr (Pathfinder), FALSE = rstan (ADVI)

# Subir limite de memoria (default 24GB puede no alcanzar con ~566K parametros)
Sys.setenv("R_MAX_VSIZE" = "48Gb")

# =============================================================================
# 1. Librerias
# =============================================================================
if (USE_CMDSTANR) {
  library(cmdstanr)
} else {
  library(rstan)
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
}

# =============================================================================
# 2. Carga de datos completos desde TXT
# =============================================================================
cat("Leyendo dataset completo SB11_20232.TXT ...\n")
datos <- as.data.frame(readr::read_delim("SB11_20232.TXT", delim = "\u00AC",
                                         quote = "", show_col_types = FALSE))
cat("Filas leidas:", nrow(datos), "\n")

# =============================================================================
# 3. Limpieza y filtrado
# =============================================================================

# Variables necesarias
vars_requeridas <- c("FAMI_TIENEINTERNET", "FAMI_ESTRATOVIVIENDA",
                     "ESTU_HORASSEMANATRABAJA", "FAMI_EDUCACIONMADRE",
                     "FAMI_EDUCACIONPADRE", "PUNT_GLOBAL",
                     "COLE_CODIGO_ICFES", "COLE_COD_MCPIO_UBICACION",
                     "COLE_NATURALEZA")

# Eliminar filas con NA en variables clave
for (v in vars_requeridas) {
  datos <- datos[!is.na(datos[[v]]) & datos[[v]] != "", ]
}
cat("Despues de eliminar NA/vacios:", nrow(datos), "\n")

# Eliminar valores no validos en educacion
datos <- datos[!datos$FAMI_EDUCACIONMADRE %in% c("No sabe", "No Aplica") &
               !datos$FAMI_EDUCACIONPADRE %in% c("No sabe", "No Aplica"), ]

# Eliminar "Sin Estrato"
datos <- datos[datos$FAMI_ESTRATOVIVIENDA != "Sin Estrato", ]
cat("Despues de filtros de calidad:", nrow(datos), "\n")

# =============================================================================
# 4. Recodificaciones
# =============================================================================

# Internet: Si -> 2, No -> 1
datos$internet <- ifelse(datos$FAMI_TIENEINTERNET == "Si", 2L,
                  ifelse(datos$FAMI_TIENEINTERNET == "No", 1L, NA_integer_))

# Estrato: texto -> numerico 1-6
datos$estrato <- ifelse(datos$FAMI_ESTRATOVIVIENDA == "Estrato 1", 1L,
                 ifelse(datos$FAMI_ESTRATOVIVIENDA == "Estrato 2", 2L,
                 ifelse(datos$FAMI_ESTRATOVIVIENDA == "Estrato 3", 3L,
                 ifelse(datos$FAMI_ESTRATOVIVIENDA == "Estrato 4", 4L,
                 ifelse(datos$FAMI_ESTRATOVIVIENDA == "Estrato 5", 5L,
                 ifelse(datos$FAMI_ESTRATOVIVIENDA == "Estrato 6", 6L, NA_integer_))))))

# Horas de trabajo: "0" -> 1 (no trabaja), otro -> 2 (trabaja)
datos$horas <- ifelse(datos$ESTU_HORASSEMANATRABAJA == "0", 1L, 2L)

# Educacion: agrupar completa/incompleta, luego a factor numerico
agrupar_educacion <- function(nivel) {
  nivel <- trimws(nivel)
  if (nivel == "Ninguno") return("Ninguno")
  else if (nivel %in% c("Primaria completa", "Primaria incompleta")) return("Primaria")
  else if (nivel %in% c("Secundaria (Bachillerato) completa", "Secundaria (Bachillerato) incompleta")) return("Secundaria")
  else if (nivel %in% c("T\u00e9cnica o tecnol\u00f3gica completa", "T\u00e9cnica o tecnol\u00f3gica incompleta")) return("T\u00e9cnica/Tecnol\u00f3gica")
  else if (nivel %in% c("Educaci\u00f3n profesional completa", "Educaci\u00f3n profesional incompleta")) return("Profesional")
  else if (nivel == "Postgrado") return("Postgrado")
  else return(NA_character_)
}

niveles_educ <- c("Ninguno", "Primaria", "Secundaria",
                  "T\u00e9cnica/Tecnol\u00f3gica", "Profesional", "Postgrado")

datos$educ_madre <- as.integer(factor(
  sapply(datos$FAMI_EDUCACIONMADRE, agrupar_educacion),
  levels = niveles_educ))
datos$educ_padre <- as.integer(factor(
  sapply(datos$FAMI_EDUCACIONPADRE, agrupar_educacion),
  levels = niveles_educ))

# Naturaleza del colegio: OFICIAL -> 1, NO OFICIAL -> 2
datos$naturaleza_cod <- ifelse(datos$COLE_NATURALEZA == "OFICIAL", 1L,
                        ifelse(datos$COLE_NATURALEZA == "NO OFICIAL", 2L, NA_integer_))

# Eliminar filas con NA despues de recodificacion
datos <- datos[complete.cases(datos[, c("internet", "estrato", "horas",
                                        "educ_madre", "educ_padre",
                                        "naturaleza_cod", "PUNT_GLOBAL")]), ]
cat("Despues de recodificacion y limpieza final:", nrow(datos), "\n")

# =============================================================================
# 5. Construir indices jerarquicos (municipio -> colegio)
#    Usamos los codigos reales, renumerados a 1..K y 1..J
# =============================================================================

# Municipios: COLE_COD_MCPIO_UBICACION -> indice 1..K
munis_unicos <- sort(unique(datos$COLE_COD_MCPIO_UBICACION))
muni_map <- setNames(seq_along(munis_unicos), munis_unicos)
datos$muni_idx <- as.integer(muni_map[as.character(datos$COLE_COD_MCPIO_UBICACION)])

# Colegios: COLE_CODIGO_ICFES -> indice 1..J
coles_unicos <- sort(unique(datos$COLE_CODIGO_ICFES))
cole_map <- setNames(seq_along(coles_unicos), coles_unicos)
datos$cole_idx <- as.integer(cole_map[as.character(datos$COLE_CODIGO_ICFES)])

# Tabla colegio-nivel: un registro por colegio con su municipio y naturaleza
datos_cole <- datos[!duplicated(datos$cole_idx), c("cole_idx", "muni_idx", "naturaleza_cod")]
datos_cole <- datos_cole[order(datos_cole$cole_idx), ]

# Verificar que todos los colegios tienen asignacion
stopifnot(nrow(datos_cole) == length(coles_unicos))
stopifnot(all(datos_cole$cole_idx == 1:nrow(datos_cole)))

# Dimensiones
N  <- nrow(datos)
J  <- length(coles_unicos)
K  <- length(munis_unicos)
N2 <- J
M  <- length(unique(datos$educ_madre))
P  <- length(unique(datos$educ_padre))
L  <- 2L

cat("\n========== Dimensiones del modelo ==========\n")
cat("  Estudiantes (N):  ", N, "\n")
cat("  Colegios (J):     ", J, "\n")
cat("  Municipios (K):   ", K, "\n")
cat("  Niveles educ madre (M):", M, "\n")
cat("  Niveles educ padre (P):", P, "\n")

# =============================================================================
# 6. Preparar datos para Stan
# =============================================================================
stan_data <- list(
  N  = N,
  J  = J,
  K  = K,
  H  = 2L,
  E  = 6L,
  I  = 2L,
  L  = L,
  N2 = N2,
  M  = M,
  P  = P,
  cole        = datos$cole_idx,
  muni        = datos$muni_idx,
  horas       = datos$horas,
  estrato     = datos$estrato,
  internet    = datos$internet,
  educ_madre  = datos$educ_madre,
  educ_padre  = datos$educ_padre,
  naturaleza  = datos_cole$naturaleza_cod,   # vector de largo J (1 por colegio)
  y           = as.numeric(datos$PUNT_GLOBAL),
  tau         = 0.8,
  col2        = datos_cole$cole_idx,         # vector 1:J
  muni2       = datos_cole$muni_idx          # municipio de cada colegio
)

# Validaciones rapidas
stopifnot(length(stan_data$naturaleza) == J)
stopifnot(length(stan_data$col2) == N2)
stopifnot(length(stan_data$muni2) == N2)
stopifnot(all(stan_data$cole >= 1 & stan_data$cole <= J))
stopifnot(all(stan_data$muni >= 1 & stan_data$muni <= K))

cat("\nDatos Stan preparados. Iniciando ajuste variacional...\n")

# =============================================================================
# 7. Ajuste con Variational Bayes
# =============================================================================

if (USE_CMDSTANR) {
  # ---- Opcion A: cmdstanr con Pathfinder (recomendado) ----
  mod <- cmdstan_model("modelo_completo_vb.stan")

  fit_vb <- mod$pathfinder(
    data = stan_data,
    num_paths = 4,
    draws = 4000,
    seed = 42
  )

  cat("\n========== Resumen Pathfinder ==========\n")
  print(fit_vb$summary())

  # Guardar
  fit_vb$save_object(file = "vb_modelo_completo_full.rds")

} else {
  # ---- Opcion B: rstan::vb() con ADVI ----
  stan_mod <- stan_model(file = "modelo_completo_vb.stan")

  fit_vb <- vb(
    stan_mod,
    data = stan_data,
    algorithm = "meanfield",
    iter = 50000,
    output_samples = 500,      # reducido: 566K params x 4000 = 18GB, no cabe
    tol_rel_obj = 0.01,        # tolerancia mas permisiva (0.001 era muy estricto)
    eta = 0.01,                # step-size mas bajo para evitar oscilaciones
    adapt_engaged = TRUE,
    adapt_iter = 500,          # mas iteraciones de adaptacion (default 50)
    seed = 42
  )

  cat("\n========== Resumen ADVI (rstan) ==========\n")
  print(fit_vb, pars = c("mu_global", "sigma_global", "sigma_cole", "sigma_ald",
                          "horas_cen", "estrato_cen", "internet_cen",
                          "educ_madre_cen", "educ_padre_cen", "naturaleza_cen"))

  # Guardar
  save(fit_vb, file = "vb_modelo_completo_full.RData")
}

cat("\n>>> Ajuste variacional completado exitosamente.\n")

# =============================================================================
# 8. Diagnosticos
# =============================================================================

if (USE_CMDSTANR) {
  draws <- fit_vb$draws(format = "df")
} else {
  draws <- extract(fit_vb)
}

cat("\n========== Parametros globales (media +/- sd) ==========\n")
cat("  mu_global:    ", mean(draws$mu_global), " (", sd(draws$mu_global), ")\n")
cat("  sigma_global: ", mean(draws$sigma_global), " (", sd(draws$sigma_global), ")\n")
cat("  sigma_cole:   ", mean(draws$sigma_cole), " (", sd(draws$sigma_cole), ")\n")
cat("  sigma_ald:    ", mean(draws$sigma_ald), " (", sd(draws$sigma_ald), ")\n")

if (!USE_CMDSTANR) {
  cat("\n========== Efectos covariables (media VB) ==========\n")
  cat("  horas_cen:       ", colMeans(draws$horas_cen), "\n")
  cat("  estrato_cen:     ", colMeans(draws$estrato_cen), "\n")
  cat("  internet_cen:    ", colMeans(draws$internet_cen), "\n")
  cat("  educ_madre_cen:  ", colMeans(draws$educ_madre_cen), "\n")
  cat("  educ_padre_cen:  ", colMeans(draws$educ_padre_cen), "\n")
  cat("  naturaleza_cen:  ", colMeans(draws$naturaleza_cen), "\n")
}

# =============================================================================
# 9. Guardar mapeo de indices para interpretacion posterior
# =============================================================================
mapeo <- list(
  municipios = data.frame(idx = seq_along(munis_unicos), cod_mcpio = munis_unicos),
  colegios   = data.frame(idx = seq_along(coles_unicos), cod_icfes = coles_unicos)
)
save(mapeo, file = "mapeo_indices_vb.RData")
cat("\nMapeo de indices guardado en mapeo_indices_vb.RData\n")
cat("  (para saber que alpha[j] corresponde a cual colegio ICFES)\n")
