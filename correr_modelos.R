###############################################################################
# Regresion cuantilica jerarquica bayesiana via INLA-EM
#
# Modelo completo: estrato + educ_madre + educ_padre + internet + horas + naturaleza
# Con efectos aleatorios: municipio + colegio
#
# Dos tipos de contraste:
#   - sum-to-zero (ANOVA bayesiana)
#   - treatment   (nivel de referencia)
###############################################################################

library(INLA)
library(ggplot2)

theme_set(theme_minimal(base_size = 13))

# =============================================================================
# 1. Carga y limpieza de datos
# =============================================================================
cat("========================================\n")
cat("Leyendo dataset completo SB11_20232.TXT\n")
cat("========================================\n")
datos_raw <- as.data.frame(readr::read_delim("SB11_20232.TXT", delim = "\u00AC",
                                              quote = "", show_col_types = FALSE))
cat("Filas leidas:", nrow(datos_raw), "\n")

vars_requeridas <- c("FAMI_TIENEINTERNET", "FAMI_ESTRATOVIVIENDA",
                     "ESTU_HORASSEMANATRABAJA", "FAMI_EDUCACIONMADRE",
                     "FAMI_EDUCACIONPADRE", "PUNT_GLOBAL",
                     "COLE_CODIGO_ICFES", "COLE_COD_MCPIO_UBICACION",
                     "COLE_NATURALEZA")

datos <- datos_raw
for (v in vars_requeridas) {
  datos <- datos[!is.na(datos[[v]]) & datos[[v]] != "", ]
}
datos <- datos[!datos$FAMI_EDUCACIONMADRE %in% c("No sabe", "No Aplica") &
               !datos$FAMI_EDUCACIONPADRE %in% c("No sabe", "No Aplica"), ]
datos <- datos[datos$FAMI_ESTRATOVIVIENDA != "Sin Estrato", ]

# Recodificaciones
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

munis_unicos <- sort(unique(datos$COLE_COD_MCPIO_UBICACION))
coles_unicos <- sort(unique(datos$COLE_CODIGO_ICFES))
datos$muni_idx <- as.integer(factor(datos$COLE_COD_MCPIO_UBICACION, levels = munis_unicos))
datos$cole_idx <- as.integer(factor(datos$COLE_CODIGO_ICFES, levels = coles_unicos))
datos$y <- as.numeric(datos$PUNT_GLOBAL)

cat("Datos limpios:", nrow(datos), "\n")
cat("Colegios:", length(coles_unicos), "  Municipios:", length(munis_unicos), "\n\n")

# =============================================================================
# 2. Graficos descriptivos
# =============================================================================
cat("========================================\n")
cat("Generando graficos descriptivos...\n")
cat("========================================\n")

dir.create("figuras", showWarnings = FALSE)

# --- 2a. Distribucion de estudiantes por estrato (barras + conteo) ---
tab_estrato <- as.data.frame(table(datos$estrato))
names(tab_estrato) <- c("Estrato", "n")
tab_estrato$pct <- round(100 * tab_estrato$n / sum(tab_estrato$n), 1)
tab_estrato$label <- paste0(format(tab_estrato$n, big.mark = ","), "\n(", tab_estrato$pct, "%)")

cat("\nDistribucion por estrato:\n")
print(tab_estrato[, c("Estrato", "n", "pct")])

p1 <- ggplot(tab_estrato, aes(x = Estrato, y = n)) +
  geom_col(fill = "#3B82F6", width = 0.7) +
  geom_text(aes(label = label), vjust = -0.3, size = 3.5) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Numero de estudiantes por estrato",
       subtitle = paste0("N = ", format(nrow(datos), big.mark = ",")),
       x = "Estrato", y = "Estudiantes") +
  theme(plot.title = element_text(face = "bold"))
ggsave("figuras/01_n_por_estrato.png", p1, width = 7, height = 5, dpi = 150)

# --- 2b. Boxplot puntaje por estrato ---
p2 <- ggplot(datos, aes(x = estrato, y = y)) +
  geom_boxplot(fill = "#3B82F6", alpha = 0.6, outlier.size = 0.3, outlier.alpha = 0.1) +
  geom_hline(yintercept = quantile(datos$y, 0.8), linetype = "dashed", color = "red") +
  annotate("text", x = 6.4, y = quantile(datos$y, 0.8) + 5,
           label = "P80 global", color = "red", size = 3.5) +
  labs(title = "Distribucion de PUNT_GLOBAL por estrato",
       subtitle = "Linea roja = percentil 80 global",
       x = "Estrato", y = "Puntaje global")
ggsave("figuras/02_boxplot_estrato.png", p2, width = 7, height = 5, dpi = 150)

# --- 2c. P80 por estrato (lo que el modelo estima) ---
p80_estrato <- aggregate(y ~ estrato, data = datos, FUN = function(x) quantile(x, 0.8))
names(p80_estrato)[2] <- "P80"
p80_estrato$n <- tab_estrato$n

p3 <- ggplot(p80_estrato, aes(x = estrato, y = P80)) +
  geom_col(fill = "#10B981", width = 0.7) +
  geom_text(aes(label = round(P80, 0)), vjust = -0.5, size = 4) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(title = "Percentil 80 observado de PUNT_GLOBAL por estrato",
       x = "Estrato", y = "P80 puntaje global")
ggsave("figuras/03_p80_por_estrato.png", p3, width = 7, height = 5, dpi = 150)

# --- 2d. Estrato vs educacion madre (tabla de calor) ---
tab_cross <- as.data.frame(prop.table(table(datos$estrato, datos$educ_madre), margin = 1))
names(tab_cross) <- c("Estrato", "Educ_madre", "Prop")
tab_cross$Educ_madre <- factor(tab_cross$Educ_madre, levels = niveles_educ)

p4 <- ggplot(tab_cross, aes(x = Educ_madre, y = Estrato, fill = Prop)) +
  geom_tile(color = "white") +
  geom_text(aes(label = paste0(round(100 * Prop), "%")), size = 3.5) +
  scale_fill_gradient(low = "white", high = "#3B82F6", labels = scales::percent) +
  labs(title = "Educacion de la madre por estrato",
       subtitle = "Porcentaje dentro de cada estrato (fila)",
       x = "Educacion madre", y = "Estrato", fill = "Proporcion") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("figuras/04_heatmap_estrato_educ_madre.png", p4, width = 8, height = 5, dpi = 150)

# --- 2e. Violin plot puntaje por educacion madre ---
p5 <- ggplot(datos, aes(x = educ_madre, y = y)) +
  geom_violin(fill = "#8B5CF6", alpha = 0.5, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  labs(title = "Distribucion de PUNT_GLOBAL por educacion de la madre",
       x = "Educacion madre", y = "Puntaje global") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("figuras/05_violin_educ_madre.png", p5, width = 8, height = 5, dpi = 150)

# --- 2f. P80 por estrato y naturaleza del colegio ---
p80_est_nat <- aggregate(y ~ estrato + naturaleza, data = datos,
                          FUN = function(x) quantile(x, 0.8))
names(p80_est_nat)[3] <- "P80"

p6 <- ggplot(p80_est_nat, aes(x = estrato, y = P80, fill = naturaleza)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = round(P80, 0), group = naturaleza),
            position = position_dodge(width = 0.7), vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("NO OFICIAL" = "#3B82F6", "OFICIAL" = "#EF4444")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(title = "P80 de PUNT_GLOBAL por estrato y naturaleza del colegio",
       x = "Estrato", y = "P80 puntaje global", fill = "Naturaleza")
ggsave("figuras/06_p80_estrato_naturaleza.png", p6, width = 8, height = 5, dpi = 150)

# --- 2g. Internet y horas por estrato ---
tab_inet <- as.data.frame(prop.table(table(datos$estrato, datos$internet), margin = 1))
names(tab_inet) <- c("Estrato", "Internet", "Prop")

p7 <- ggplot(tab_inet, aes(x = Estrato, y = Prop, fill = Internet)) +
  geom_col(width = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("No" = "#EF4444", "Si" = "#10B981")) +
  labs(title = "Acceso a internet por estrato",
       x = "Estrato", y = "Proporcion", fill = "Internet")
ggsave("figuras/07_internet_por_estrato.png", p7, width = 7, height = 5, dpi = 150)

cat(">>> Graficos guardados en figuras/\n\n")

# =============================================================================
# 3. Funcion EM-INLA reutilizable
# =============================================================================
ajustar_inla_ald <- function(datos, formula_mod, tau = 0.8,
                              max_iter = 50, tol = 0.01, damp = 0.3,
                              sigma_init = 10, verbose_em = TRUE) {

  N <- nrow(datos)
  tau_comp <- tau * (1 - tau)
  delta_unit <- (1 - 2 * tau) / tau_comp
  gamma_unit <- 2 / tau_comp

  sigma_ald <- sigma_init
  w <- rep(sigma_ald, N)
  clip_w <- function(w_vec) pmax(pmin(w_vec, 1e4), 1e-6)
  w <- clip_w(w)

  for (iter in 1:max_iter) {
    t0 <- Sys.time()

    delta <- sigma_ald * delta_unit
    gamma_sq <- sigma_ald^2 * gamma_unit

    datos$y_tilde <- datos$y - delta * w
    scale_vec <- clip_w(1 / w)

    fit <- inla(
      formula_mod,
      family = "gaussian",
      data = datos,
      scale = scale_vec,
      control.compute = list(config = FALSE),
      control.predictor = list(compute = TRUE, link = 1),
      control.inla = list(strategy = "gaussian", int.strategy = "eb"),
      num.threads = parallel::detectCores(),
      verbose = FALSE,
      safe = FALSE
    )

    mu_fitted <- fit$summary.linear.predictor[, "mean"]
    prec_gauss <- fit$summary.hyperpar["Precision for the Gaussian observations", "mean"]
    sigma_ald_new <- sqrt(1 / (prec_gauss * gamma_unit))

    resid <- datos$y - mu_fitted
    a_gig <- delta^2 / gamma_sq + 2 / sigma_ald
    b_gig <- pmax(resid^2 / gamma_sq, 1e-10)
    sqrt_ab <- sqrt(a_gig * b_gig)
    w_new <- clip_w(sqrt(b_gig / a_gig) * (1 + 1 / pmax(sqrt_ab, 1e-10)))

    sigma_ald_d <- damp * sigma_ald_new + (1 - damp) * sigma_ald
    w_d <- damp * w_new + (1 - damp) * w

    rel_change <- abs(sigma_ald_d - sigma_ald) / max(abs(sigma_ald), 1e-10)
    sigma_ald <- sigma_ald_d
    w <- w_d

    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (verbose_em) {
      cat(sprintf("  Iter %2d: sigma_ald=%.3f, cambio=%.5f (%.1fs)\n",
                  iter, sigma_ald, rel_change, elapsed))
    }

    if (rel_change < tol && iter > 2) {
      if (verbose_em) cat("  >>> Convergencia!\n")
      break
    }
  }

  # Ajuste final completo
  datos$y_tilde <- datos$y - sigma_ald * delta_unit * w
  scale_final <- clip_w(1 / w)

  fit_final <- inla(
    formula_mod,
    family = "gaussian",
    data = datos,
    scale = scale_final,
    control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
    control.predictor = list(compute = TRUE, link = 1),
    num.threads = parallel::detectCores(),
    verbose = FALSE
  )

  return(list(
    fit = fit_final,
    sigma_ald = sigma_ald,
    w = w,
    tau = tau,
    n_iter = min(iter, max_iter)
  ))
}

# =============================================================================
# 4. Funcion para imprimir resultados
# =============================================================================
imprimir_resultados <- function(res, nombre, contraste, factores_info) {
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("MODELO:", nombre, " | Contraste:", contraste, "\n")
  cat(strrep("=", 70), "\n")

  cat("\nEfectos fijos (cuantil", res$tau, "):\n")
  print(round(res$fit$summary.fixed[, c("mean", "sd", "0.025quant", "0.975quant")], 3))

  if (contraste == "sum") {
    cat("\nEfectos sum-to-zero completos (incl. categoria omitida):\n")
    for (info in factores_info) {
      fijos <- res$fit$summary.fixed
      idx <- grep(paste0("^", info$prefijo), rownames(fijos))
      if (length(idx) > 0) {
        medias <- fijos[idx, "mean"]
        ultima <- -sum(medias)
        efectos <- c(medias, ultima)
        names(efectos) <- info$niveles
        cat("\n ", info$nombre, ":\n")
        print(round(efectos, 3))
      }
    }
  }

  cat("\nSD efectos aleatorios:\n")
  for (h in rownames(res$fit$summary.hyperpar)) {
    if (grepl("Precision", h)) {
      sd_est <- 1 / sqrt(res$fit$summary.hyperpar[h, "mean"])
      cat("  ", h, "-> SD =", round(sd_est, 2), "\n")
    }
  }
  cat("  sigma_ald:", round(res$sigma_ald, 3), "\n")
  cat("  DIC:", round(res$fit$dic$dic, 1), "\n")
  cat("  EM convergio en", res$n_iter, "iteraciones\n")
}

# =============================================================================
# 5. Definir formula del modelo completo
# =============================================================================
re_terms <- "+ f(muni_idx, model='iid', hyper=list(prec=list(prior='pc.prec', param=c(50,0.01)))) + f(cole_idx, model='iid', hyper=list(prec=list(prior='pc.prec', param=c(50,0.01))))"

formula_completo <- as.formula(paste(
  "y_tilde ~ 1 + estrato + internet + horas + educ_madre + educ_padre + naturaleza",
  re_terms))

factores_info <- list(
  list(prefijo = "estrato",    nombre = "Estrato",     niveles = levels(datos$estrato)),
  list(prefijo = "internet",   nombre = "Internet",    niveles = levels(datos$internet)),
  list(prefijo = "horas",      nombre = "Horas",       niveles = levels(datos$horas)),
  list(prefijo = "educ_madre", nombre = "Educ madre",  niveles = niveles_educ),
  list(prefijo = "educ_padre", nombre = "Educ padre",  niveles = niveles_educ),
  list(prefijo = "naturaleza", nombre = "Naturaleza",  niveles = levels(datos$naturaleza))
)

# =============================================================================
# 6. Correr modelo completo con ambos contrastes
# =============================================================================
resultados <- list()

for (tipo_contraste in c("sum", "treatment")) {

  cat("\n\n")
  cat(strrep("#", 70), "\n")
  cat("# MODELO COMPLETO - CONTRASTES:", toupper(tipo_contraste), "\n")
  cat(strrep("#", 70), "\n")

  if (tipo_contraste == "sum") {
    contrasts(datos$estrato)    <- contr.sum(levels(datos$estrato))
    contrasts(datos$internet)   <- contr.sum(levels(datos$internet))
    contrasts(datos$horas)      <- contr.sum(levels(datos$horas))
    contrasts(datos$educ_madre) <- contr.sum(levels(datos$educ_madre))
    contrasts(datos$educ_padre) <- contr.sum(levels(datos$educ_padre))
    contrasts(datos$naturaleza) <- contr.sum(levels(datos$naturaleza))
  } else {
    contrasts(datos$estrato)    <- contr.treatment(levels(datos$estrato))
    contrasts(datos$internet)   <- contr.treatment(levels(datos$internet))
    contrasts(datos$horas)      <- contr.treatment(levels(datos$horas))
    contrasts(datos$educ_madre) <- contr.treatment(levels(datos$educ_madre))
    contrasts(datos$educ_padre) <- contr.treatment(levels(datos$educ_padre))
    contrasts(datos$naturaleza) <- contr.treatment(levels(datos$naturaleza))
  }

  key <- paste0("completo_", tipo_contraste)
  cat("\n>>> Ajustando modelo completo (", tipo_contraste, ")\n")
  resultados[[key]] <- ajustar_inla_ald(datos, formula_completo)
  imprimir_resultados(resultados[[key]], "Completo", tipo_contraste, factores_info)
}

# =============================================================================
# 7. Tabla comparativa
# =============================================================================
cat("\n\n")
cat(strrep("=", 70), "\n")
cat("TABLA COMPARATIVA\n")
cat(strrep("=", 70), "\n\n")

cat(sprintf("%-30s %10s %10s %10s\n", "Modelo", "DIC", "sigma_ald", "Iters"))
cat(strrep("-", 62), "\n")
for (nm in names(resultados)) {
  r <- resultados[[nm]]
  cat(sprintf("%-30s %10.1f %10.3f %10d\n",
              nm, r$fit$dic$dic, r$sigma_ald, r$n_iter))
}

# =============================================================================
# 8. Graficos de resultados del modelo
# =============================================================================
cat("\n>>> Generando graficos de resultados...\n")

# --- 8a. Efectos sum-to-zero de todas las variables ---
res_sum <- resultados[["completo_sum"]]
fijos <- res_sum$fit$summary.fixed

reconstruir_stz <- function(fit, prefijo, niveles) {
  fijos <- fit$summary.fixed
  idx <- grep(paste0("^", prefijo), rownames(fijos))
  medias_mean <- fijos[idx, "mean"]
  medias_lo   <- fijos[idx, "0.025quant"]
  medias_hi   <- fijos[idx, "0.975quant"]
  # Ultima categoria
  ult_mean <- -sum(medias_mean)
  ult_lo   <- NA  # no hay IC exacto para la suma
  ult_hi   <- NA
  data.frame(
    nivel = niveles,
    mean  = c(medias_mean, ult_mean),
    lo    = c(medias_lo, ult_lo),
    hi    = c(medias_hi, ult_hi)
  )
}

vars_plot <- list(
  list(prefijo = "estrato",    nombre = "Estrato",       niveles = levels(datos$estrato)),
  list(prefijo = "educ_madre", nombre = "Educ. madre",   niveles = niveles_educ),
  list(prefijo = "educ_padre", nombre = "Educ. padre",   niveles = niveles_educ),
  list(prefijo = "internet",   nombre = "Internet",      niveles = levels(datos$internet)),
  list(prefijo = "horas",      nombre = "Horas trabajo",  niveles = levels(datos$horas)),
  list(prefijo = "naturaleza", nombre = "Naturaleza",    niveles = levels(datos$naturaleza))
)

all_effects <- do.call(rbind, lapply(vars_plot, function(v) {
  df <- reconstruir_stz(res_sum$fit, v$prefijo, v$niveles)
  df$variable <- v$nombre
  df$nivel <- factor(df$nivel, levels = v$niveles)
  df
}))
all_effects$variable <- factor(all_effects$variable,
                                levels = sapply(vars_plot, function(v) v$nombre))

p8 <- ggplot(all_effects, aes(x = nivel, y = mean)) +
  geom_col(fill = "#3B82F6", width = 0.7) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.3, na.rm = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  facet_wrap(~ variable, scales = "free_x", ncol = 3) +
  labs(title = "Efectos sum-to-zero sobre el P80 de PUNT_GLOBAL",
       subtitle = "Modelo completo con RE de municipio y colegio | Barras de error = IC 95%",
       x = NULL, y = "Efecto (puntos respecto a la media global)") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        strip.text = element_text(face = "bold"))
ggsave("figuras/08_efectos_stz_completo.png", p8, width = 12, height = 8, dpi = 150)

# --- 8b. Efectos treatment (diferencias respecto a referencia) ---
res_trt <- resultados[["completo_treatment"]]
fijos_trt <- res_trt$fit$summary.fixed
# Excluir intercepto
fijos_trt <- fijos_trt[rownames(fijos_trt) != "(Intercept)", ]

df_trt <- data.frame(
  coef = rownames(fijos_trt),
  mean = fijos_trt[, "mean"],
  lo   = fijos_trt[, "0.025quant"],
  hi   = fijos_trt[, "0.975quant"]
)
df_trt$coef <- factor(df_trt$coef, levels = rev(df_trt$coef))
df_trt$signo <- ifelse(df_trt$mean > 0, "positivo", "negativo")

p9 <- ggplot(df_trt, aes(x = coef, y = mean, color = signo)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c("positivo" = "#10B981", "negativo" = "#EF4444"), guide = "none") +
  coord_flip() +
  labs(title = "Efectos fijos (treatment) sobre el P80 de PUNT_GLOBAL",
       subtitle = "Referencia: Estrato 1, Ninguno (educ), No internet, No trabaja, No oficial",
       x = NULL, y = "Efecto (puntos) con IC 95%")
ggsave("figuras/09_forest_treatment.png", p9, width = 10, height = 7, dpi = 150)

cat(">>> Graficos de resultados guardados en figuras/\n")

# =============================================================================
# 9. Guardar todo
# =============================================================================
save(resultados, file = "resultados_modelo_completo.RData")

mapeo <- list(
  municipios = data.frame(idx = seq_along(munis_unicos), cod_mcpio = munis_unicos),
  colegios   = data.frame(idx = seq_along(coles_unicos), cod_icfes = coles_unicos)
)
save(mapeo, file = "mapeo_indices.RData")

cat("\n>>> Resultados guardados en resultados_modelo_completo.RData\n")
cat(">>> Mapeo en mapeo_indices.RData\n")
cat(">>> Graficos en figuras/\n")
