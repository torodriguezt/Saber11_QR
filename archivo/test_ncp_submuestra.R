###############################################################################
# Test rapido: NCP + MCMC con submuestra de 10K observaciones
# Objetivo: verificar si la non-centered parameterization resuelve
# la convergencia antes de invertir horas en el dataset completo.
###############################################################################

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

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

# =============================================================================
# 3. Recodificaciones
# =============================================================================
datos$internet <- ifelse(datos$FAMI_TIENEINTERNET == "Si", 2L,
                  ifelse(datos$FAMI_TIENEINTERNET == "No", 1L, NA_integer_))

datos$estrato <- ifelse(datos$FAMI_ESTRATOVIVIENDA == "Estrato 1", 1L,
                 ifelse(datos$FAMI_ESTRATOVIVIENDA == "Estrato 2", 2L,
                 ifelse(datos$FAMI_ESTRATOVIVIENDA == "Estrato 3", 3L,
                 ifelse(datos$FAMI_ESTRATOVIVIENDA == "Estrato 4", 4L,
                 ifelse(datos$FAMI_ESTRATOVIVIENDA == "Estrato 5", 5L,
                 ifelse(datos$FAMI_ESTRATOVIVIENDA == "Estrato 6", 6L, NA_integer_))))))

datos$horas <- ifelse(datos$ESTU_HORASSEMANATRABAJA == "0", 1L, 2L)

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
  sapply(datos$FAMI_EDUCACIONMADRE, agrupar_educacion), levels = niveles_educ))
datos$educ_padre <- as.integer(factor(
  sapply(datos$FAMI_EDUCACIONPADRE, agrupar_educacion), levels = niveles_educ))

datos$naturaleza_cod <- ifelse(datos$COLE_NATURALEZA == "OFICIAL", 1L,
                        ifelse(datos$COLE_NATURALEZA == "NO OFICIAL", 2L, NA_integer_))

datos <- datos[complete.cases(datos[, c("internet", "estrato", "horas",
                                        "educ_madre", "educ_padre",
                                        "naturaleza_cod", "PUNT_GLOBAL")]), ]
cat("Datos limpios:", nrow(datos), "\n")

# =============================================================================
# 4. Submuestra de 10K
# =============================================================================
set.seed(42)
idx_sub <- sample(nrow(datos), min(10000, nrow(datos)))
datos_sub <- datos[idx_sub, ]

# Reindexar colegios y municipios para la submuestra (1..J_sub, 1..K_sub)
coles_sub <- sort(unique(datos_sub$COLE_CODIGO_ICFES))
munis_sub <- sort(unique(datos_sub$COLE_COD_MCPIO_UBICACION))
cole_remap <- setNames(seq_along(coles_sub), coles_sub)
muni_remap <- setNames(seq_along(munis_sub), munis_sub)

datos_sub$cole_idx <- as.integer(cole_remap[as.character(datos_sub$COLE_CODIGO_ICFES)])
datos_sub$muni_idx <- as.integer(muni_remap[as.character(datos_sub$COLE_COD_MCPIO_UBICACION)])

# Tabla colegio: un registro por colegio
datos_cole <- datos_sub[!duplicated(datos_sub$cole_idx),
                         c("cole_idx", "muni_idx", "naturaleza_cod")]
datos_cole <- datos_cole[order(datos_cole$cole_idx), ]

J_sub <- length(coles_sub)
K_sub <- length(munis_sub)
M <- length(unique(datos_sub$educ_madre))
P <- length(unique(datos_sub$educ_padre))

cat("\n========== Submuestra ==========\n")
cat("  Estudiantes (N): ", nrow(datos_sub), "\n")
cat("  Colegios (J):    ", J_sub, "\n")
cat("  Municipios (K):  ", K_sub, "\n")

# =============================================================================
# 5. Datos Stan
# =============================================================================
stan_data <- list(
  N = nrow(datos_sub), J = J_sub, K = K_sub,
  H = 2L, E = 6L, I = 2L, L = 2L, N2 = J_sub, M = M, P = P,
  cole = datos_sub$cole_idx, muni = datos_sub$muni_idx,
  horas = datos_sub$horas, estrato = datos_sub$estrato,
  internet = datos_sub$internet,
  educ_madre = datos_sub$educ_madre, educ_padre = datos_sub$educ_padre,
  naturaleza = datos_cole$naturaleza_cod,
  y = as.numeric(datos_sub$PUNT_GLOBAL), tau = 0.8,
  col2 = datos_cole$cole_idx, muni2 = datos_cole$muni_idx
)

# =============================================================================
# 6. MCMC test
# =============================================================================
cat("\nCompilando modelo NCP...\n")
stan_mod <- stan_model(file = "modelo_ncp.stan")

cat("\n>>> Corriendo 1000 iter, 2 cadenas en submuestra de 10K...\n")
fit_test <- sampling(
  stan_mod,
  data = stan_data,
  chains = 2,
  iter = 1000,
  warmup = 500,
  cores = 2,
  refresh = 100
)

# =============================================================================
# 7. Diagnosticos
# =============================================================================
cat("\n========== Resultados ==========\n")
params_clave <- c("mu_global", "sigma_global", "sigma_cole", "sigma",
                   "horas_cen", "estrato_cen", "internet_cen",
                   "educ_madre_cen", "educ_padre_cen", "naturaleza_cen")
print(fit_test, pars = params_clave)

# Rhat
summ <- summary(fit_test, pars = params_clave)$summary
rhats <- summ[, "Rhat"]
cat("\nRango de Rhat:", range(rhats, na.rm = TRUE), "\n")

# Divergencias
n_div <- sum(get_sampler_params(fit_test, inc_warmup = FALSE)[[1]][, "divergent__"])
n_div2 <- sum(get_sampler_params(fit_test, inc_warmup = FALSE)[[2]][, "divergent__"])
cat("Divergencias: cadena 1 =", n_div, ", cadena 2 =", n_div2, "\n")

if (all(rhats < 1.1, na.rm = TRUE) && (n_div + n_div2) < 50) {
  cat("\n>>> VEREDICTO: NCP funciona! Rhat OK y pocas divergencias.\n")
  cat(">>> Se puede proceder con el dataset completo.\n")
} else if (all(rhats < 1.1, na.rm = TRUE)) {
  cat("\n>>> VEREDICTO: Rhat OK pero hay divergencias.\n")
  cat(">>> Puede funcionar con adapt_delta mas alto (0.9 o 0.95).\n")
} else {
  cat("\n>>> VEREDICTO: NCP no resuelve el problema por si solo.\n")
  cat(">>> Hay que explorar alternativas.\n")
}

save(fit_test, file = "test_ncp_submuestra.RData")
cat("\nGuardado en test_ncp_submuestra.RData\n")
