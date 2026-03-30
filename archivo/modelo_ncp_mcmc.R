###############################################################################
# Modelo NCP (Non-Centered Parameterization) + MCMC
# Dataset COMPLETO: SB11_20232.TXT
#
# La NCP deberia resolver los problemas de convergencia del modelo original.
# Empezamos con pocas iteraciones para verificar antes de correr largo.
###############################################################################

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# =============================================================================
# 1. Carga de datos completos desde TXT
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
cat("Despues de eliminar NA/vacios:", nrow(datos), "\n")

datos <- datos[!datos$FAMI_EDUCACIONMADRE %in% c("No sabe", "No Aplica") &
               !datos$FAMI_EDUCACIONPADRE %in% c("No sabe", "No Aplica"), ]
datos <- datos[datos$FAMI_ESTRATOVIVIENDA != "Sin Estrato", ]
cat("Despues de filtros de calidad:", nrow(datos), "\n")

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
cat("Despues de limpieza final:", nrow(datos), "\n")

# =============================================================================
# 4. Indices jerarquicos
# =============================================================================
munis_unicos <- sort(unique(datos$COLE_COD_MCPIO_UBICACION))
muni_map <- setNames(seq_along(munis_unicos), munis_unicos)
datos$muni_idx <- as.integer(muni_map[as.character(datos$COLE_COD_MCPIO_UBICACION)])

coles_unicos <- sort(unique(datos$COLE_CODIGO_ICFES))
cole_map <- setNames(seq_along(coles_unicos), coles_unicos)
datos$cole_idx <- as.integer(cole_map[as.character(datos$COLE_CODIGO_ICFES)])

datos_cole <- datos[!duplicated(datos$cole_idx), c("cole_idx", "muni_idx", "naturaleza_cod")]
datos_cole <- datos_cole[order(datos_cole$cole_idx), ]

N  <- nrow(datos)
J  <- length(coles_unicos)
K  <- length(munis_unicos)
N2 <- J
M  <- length(unique(datos$educ_madre))
P  <- length(unique(datos$educ_padre))
L  <- 2L

cat("\n  N:", N, "  J:", J, "  K:", K, "  M:", M, "  P:", P, "\n")

# =============================================================================
# 5. Datos Stan
# =============================================================================
stan_data <- list(
  N = N, J = J, K = K, H = 2L, E = 6L, I = 2L, L = L, N2 = N2, M = M, P = P,
  cole = datos$cole_idx, muni = datos$muni_idx,
  horas = datos$horas, estrato = datos$estrato, internet = datos$internet,
  educ_madre = datos$educ_madre, educ_padre = datos$educ_padre,
  naturaleza = datos_cole$naturaleza_cod,
  y = as.numeric(datos$PUNT_GLOBAL), tau = 0.8,
  col2 = datos_cole$cole_idx, muni2 = datos_cole$muni_idx
)

# =============================================================================
# 6. MCMC - Test con submuestra primero, luego datos completos
# =============================================================================
cat("\nCompilando modelo NCP...\n")
stan_mod <- stan_model(file = "modelo_ncp.stan")

# --- Paso 0: Test rapido con submuestra (10K obs) ---
# Esto tarda ~minutos y nos dice si NCP resuelve la convergencia
# antes de invertir horas en el dataset completo.
set.seed(42)
idx_sub <- sample(N, min(10000, N))
datos_sub_test <- datos[idx_sub, ]

# Reconstruir indices para la submuestra
coles_sub <- sort(unique(datos_sub_test$cole_idx))
munis_sub <- sort(unique(datos_sub_test$muni_idx))
cole_remap <- setNames(seq_along(coles_sub), coles_sub)
muni_remap <- setNames(seq_along(munis_sub), munis_sub)

datos_sub_test$cole_new <- as.integer(cole_remap[as.character(datos_sub_test$cole_idx)])
datos_sub_test$muni_new <- as.integer(muni_remap[as.character(datos_sub_test$muni_idx)])

datos_cole_sub <- datos_sub_test[!duplicated(datos_sub_test$cole_new),
                                  c("cole_new", "muni_new", "naturaleza_cod")]
datos_cole_sub <- datos_cole_sub[order(datos_cole_sub$cole_new), ]

stan_data_sub <- list(
  N = nrow(datos_sub_test),
  J = length(coles_sub), K = length(munis_sub),
  H = 2L, E = 6L, I = 2L, L = L,
  N2 = length(coles_sub), M = M, P = P,
  cole = datos_sub_test$cole_new, muni = datos_sub_test$muni_new,
  horas = datos_sub_test$horas, estrato = datos_sub_test$estrato,
  internet = datos_sub_test$internet,
  educ_madre = datos_sub_test$educ_madre, educ_padre = datos_sub_test$educ_padre,
  naturaleza = datos_cole_sub$naturaleza_cod,
  y = as.numeric(datos_sub_test$PUNT_GLOBAL), tau = 0.8,
  col2 = datos_cole_sub$cole_new, muni2 = datos_cole_sub$muni_new
)

cat("\n>>> Test con submuestra: N =", stan_data_sub$N,
    ", J =", stan_data_sub$J, ", K =", stan_data_sub$K, "\n")
cat(">>> 1000 iter, 2 cadenas...\n")

fit_test <- sampling(
  stan_mod,
  data = stan_data_sub,
  chains = 2,
  iter = 1000,
  warmup = 500,
  cores = 2,
  refresh = 100,
  thin = 1
)

# Verificar convergencia con Rhat
cat("\n========== Diagnostico rapido ==========\n")
params_clave <- c("mu_global", "sigma_global", "sigma_cole", "sigma",
                   "horas_cen", "estrato_cen", "internet_cen",
                   "educ_madre_cen", "educ_padre_cen", "naturaleza_cen")
print(fit_test, pars = params_clave)

# Extraer Rhat
summ <- summary(fit_test, pars = params_clave)$summary
rhats <- summ[, "Rhat"]
cat("\nRango de Rhat:", range(rhats, na.rm = TRUE), "\n")
if (all(rhats < 1.1, na.rm = TRUE)) {
  cat(">>> Rhat OK (<1.1). El modelo NCP converge!\n")
  cat(">>> Puedes proceder con la corrida larga.\n")
} else {
  cat(">>> ADVERTENCIA: Algunos Rhat > 1.1. Revisar antes de correr largo.\n")
}

save(fit_test, file = "ncp_test_corto.RData")

# --- Paso 2: corrida larga (descomentar cuando el test pase) ---
# cat("\n>>> Corrida larga (5000 iter, 4 cadenas)...\n")
# fit_full <- sampling(
#   stan_mod,
#   data = stan_data,
#   chains = 4,
#   iter = 5000,
#   warmup = 2500,
#   cores = 4,
#   refresh = 100,
#   thin = 2
# )
#
# print(fit_full, pars = params_clave)
# save(fit_full, file = "ncp_modelo_completo_full.RData")

# =============================================================================
# 7. Mapeo de indices
# =============================================================================
mapeo <- list(
  municipios = data.frame(idx = seq_along(munis_unicos), cod_mcpio = munis_unicos),
  colegios   = data.frame(idx = seq_along(coles_unicos), cod_icfes = coles_unicos)
)
save(mapeo, file = "mapeo_indices_ncp.RData")
cat("\nMapeo guardado en mapeo_indices_ncp.RData\n")
