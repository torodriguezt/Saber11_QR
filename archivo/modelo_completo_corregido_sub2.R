library(rstan)

# Leer el archivo CSV
datos_sub <- read.csv("johnatan3/SUB_BASE2.csv", header = TRUE, sep = ",")

# Filtrar valores NA
datos_sub <- datos_sub[!is.na(datos_sub$FAMI_TIENEINTERNET) & 
                         !is.na(datos_sub$FAMI_ESTRATOVIVIENDA) & 
                         !is.na(datos_sub$ESTU_HORASSEMANATRABAJA) &
                         !is.na(datos_sub$EDAD_2023) &
                         !is.na(datos_sub$FAMI_EDUCACIONMADRE) &
                         !is.na(datos_sub$FAMI_EDUCACIONPADRE), ]

# Eliminar valores no válidos en educación
datos_sub <- datos_sub[!datos_sub$FAMI_EDUCACIONMADRE %in% c("No sabe", "No Aplica") & 
                         !datos_sub$FAMI_EDUCACIONPADRE %in% c("No sabe", "No Aplica"), ]

# Eliminar registros con "Sin Estrato"
datos_sub <- datos_sub[datos_sub$FAMI_ESTRATOVIVIENDA != "Sin Estrato", ]

# Recodificar FAMI_TIENEINTERNET: "Si" -> 2, "No" -> 1
datos_sub$FAMI_TIENEINTERNET <- ifelse(datos_sub$FAMI_TIENEINTERNET == "Si", 2,
                                       ifelse(datos_sub$FAMI_TIENEINTERNET == "No", 1, NA))

# Recodificar FAMI_ESTRATOVIVIENDA a numérico
datos_sub$FAMI_ESTRATOVIVIENDA <- ifelse(datos_sub$FAMI_ESTRATOVIVIENDA == "Estrato 1", 1,
                                         ifelse(datos_sub$FAMI_ESTRATOVIVIENDA == "Estrato 2", 2,
                                                ifelse(datos_sub$FAMI_ESTRATOVIVIENDA == "Estrato 3", 3,
                                                       ifelse(datos_sub$FAMI_ESTRATOVIVIENDA == "Estrato 4", 4,
                                                              ifelse(datos_sub$FAMI_ESTRATOVIVIENDA == "Estrato 5", 5,
                                                                     ifelse(datos_sub$FAMI_ESTRATOVIVIENDA == "Estrato 6", 6, NA))))))

# Recodificar ESTU_HORASSEMANATRABAJA a binario: trabaja = 1, no trabaja = 0
datos_sub$ESTU_HORASSEMANATRABAJA <- ifelse(datos_sub$ESTU_HORASSEMANATRABAJA == "0", 1, 2)

# Eliminar espacios y normalizar texto en niveles educativos
datos_sub$FAMI_EDUCACIONMADRE <- trimws(as.character(datos_sub$FAMI_EDUCACIONMADRE))
datos_sub$FAMI_EDUCACIONPADRE <- trimws(as.character(datos_sub$FAMI_EDUCACIONPADRE))

# Función para agrupar los niveles educativos (completos e incompletos juntos)
agrupar_educacion <- function(nivel) {
  if (nivel == "Ninguno") {
    return("Ninguno")
  } else if (nivel %in% c("Primaria completa", "Primaria incompleta")) {
    return("Primaria")
  } else if (nivel %in% c("Secundaria (Bachillerato) completa", "Secundaria (Bachillerato) incompleta")) {
    return("Secundaria")
  } else if (nivel %in% c("Técnica o tecnológica completa", "Técnica o tecnológica incompleta")) {
    return("Técnica/Tecnológica")
  } else if (nivel %in% c("Educación profesional completa", "Educación profesional incompleta")) {
    return("Profesional")
  } else if (nivel == "Postgrado") {
    return("Postgrado")
  } else {
    return(NA)  # Por si hay valores fuera de los previstos
  }
}

# Aplicar función de agrupación
datos_sub$FAMI_EDUCACIONMADRE_AGRUPADA <- sapply(datos_sub$FAMI_EDUCACIONMADRE, agrupar_educacion)
datos_sub$FAMI_EDUCACIONPADRE_AGRUPADA <- sapply(datos_sub$FAMI_EDUCACIONPADRE, agrupar_educacion)

# Convertir agrupaciones a factores con niveles ordenados
niveles_educ <- c("Ninguno", "Primaria", "Secundaria", "Técnica/Tecnológica", "Profesional", "Postgrado")
datos_sub$FAMI_EDUCACIONMADRE <- as.integer(factor(datos_sub$FAMI_EDUCACIONMADRE_AGRUPADA, levels = niveles_educ))
datos_sub$FAMI_EDUCACIONPADRE <- as.integer(factor(datos_sub$FAMI_EDUCACIONPADRE_AGRUPADA, levels = niveles_educ))

# Cambiar nombre de columna colegio
colnames(datos_sub)[colnames(datos_sub) == "colegio"] <- "Colegio"

# Naturaleza del colegio: recodificar de "OFICIAL" -> 1 y "NO OFICIAL" -> 0  
datos_sub_distinct <- datos_sub[!duplicated(datos_sub[c("ClusterMuni", "Colegio")]), ]
datos_sub_distinct$COLE_NATURALEZA <- ifelse(datos_sub_distinct$COLE_NATURALEZA == "OFICIAL", 1,
                                             ifelse(datos_sub_distinct$COLE_NATURALEZA == "NO OFICIAL", 2, NA))

# Eliminar filas con NA en las variables clave
datos_sub <- datos_sub[!is.na(datos_sub$FAMI_TIENEINTERNET) & 
                         !is.na(datos_sub$FAMI_ESTRATOVIVIENDA) & 
                         !is.na(datos_sub$ESTU_HORASSEMANATRABAJA) &
                         !is.na(datos_sub$EDAD_2023) &
                         !is.na(datos_sub$FAMI_EDUCACIONMADRE) &
                         !is.na(datos_sub$FAMI_EDUCACIONPADRE), ]

# Convertir variables a enteros
datos_sub$FAMI_TIENEINTERNET <- as.integer(datos_sub$FAMI_TIENEINTERNET)
datos_sub$FAMI_ESTRATOVIVIENDA <- as.integer(datos_sub$FAMI_ESTRATOVIVIENDA)
datos_sub$ESTU_HORASSEMANATRABAJA <- as.integer(datos_sub$ESTU_HORASSEMANATRABAJA)
datos_sub$FAMI_EDUCACIONMADRE <- as.integer(datos_sub$FAMI_EDUCACIONMADRE)
datos_sub$FAMI_EDUCACIONPADRE <- as.integer(datos_sub$FAMI_EDUCACIONPADRE)

# Asegurarse de que la variable Colegio sea numérica
datos_sub_distinct$Colegio <- as.integer(datos_sub_distinct$Colegio)
datos_sub$Colegio <- as.integer(datos_sub$Colegio)

# Definir el número total de estudiantes, colegios y municipios ----------------
N <- nrow(datos_sub)
K <- length(unique(datos_sub$ClusterMuni))

# Para los colegios se usa el dataframe sin duplicados
datos_unicos <- datos_sub[!duplicated(datos_sub[c("ClusterMuni", "Colegio")]), ]
conteo_colegios_por_cluster <- aggregate(Colegio ~ ClusterMuni, data = datos_unicos, FUN = length)
total_colegios <- sum(conteo_colegios_por_cluster$Colegio)
J <- total_colegios
N2 <- total_colegios

# Definir número de niveles para las variables nuevas
M <- length(unique(datos_sub$FAMI_EDUCACIONMADRE))
P <- length(unique(datos_sub$FAMI_EDUCACIONPADRE))
L <- 2  


stan_data <- list(
  N  = N,
  J  = J,
  K  = K,
  H  = 2,
  N2 = N2,
  E  = 6,
  I  = 2,
  L  = L,                              # <── NUEVO
  M  = M,
  P  = P,
  cole        = datos_sub$Colegio,
  muni        = datos_sub$ClusterMuni,
  horas       = datos_sub$ESTU_HORASSEMANATRABAJA,
  estrato     = datos_sub$FAMI_ESTRATOVIVIENDA,
  internet    = datos_sub$FAMI_TIENEINTERNET,
  educ_madre  = datos_sub$FAMI_EDUCACIONMADRE,
  educ_padre  = datos_sub$FAMI_EDUCACIONPADRE,
  naturaleza  = datos_sub_distinct$COLE_NATURALEZA,
  y           = datos_sub$PUNT_GLOBAL,
  tau         = 0.8,
  col2        = datos_sub_distinct$Colegio,
  muni2       = datos_sub_distinct$ClusterMuni
)

# Ajuste del modelo en Stan --------------------------------------------------
# Función para generar los puntos iniciales de todos los parámetros basados en las medias posteriores
# Función para generar los puntos iniciales de todos los parámetros basados en las medias posteriores
init_fn <- function() {
  list(
    # Inicialización de los parámetros alpha con sus medias posteriores
    alpha = c(262.20, 312.04, 345.54, 284.02, 271.08, 243.56, 292.65, 
              248.29, 229.60, 282.13, 269.55, 300.36, 366.31, 329.06),
    
    # Inicialización de los parámetros mu con sus medias posteriores
    mu = c(291.69, 280.09, 274.82, 297.43),
    
    # Inicialización de mu_global con su media posterior
    mu_global = 285.07,
    
    # Inicialización de los parámetros sigma con sus medias posteriores
    sigma = 12.37,
    sigma_cole = 36.94,
    sigma_global = 18.46,
    
    # Inicialización de las covariables
    horas_cen = c(4.32, -4.32),  # Asegúrate de que estas sean las medias de tus covariables
    estrato_cen = c(8.58, 5.23, 1.88, -6.80, -0.88, -8.02),  # Ejemplo de valores
    internet_cen = c(-5.72, 5.72),  # Valores para la covariable de internet
    educ_madre_cen = c(-18.22, -11.80, -3.94, 7.83, 9.75, 16.38),  # Media de las covariables de educación madre
    educ_padre_cen = c(-11.39, -6.64, -4.03, 3.82, 4.31, 13.93),  # Media de las covariables de educación padre
    naturaleza_cen = c(-17.70, 17.70),  # Para las covariables de naturaleza, si son dos, agrega las medias
    lp__ = 0  # Solo un valor arbitrario para lp__ (log posterior)
  )
}

# Ajuste del modelo en Stan con los puntos iniciales
fit <- stan(data = stan_data, 
            file = "johnatan3/modelo_completo.stan",
            chains = 4,
            iter = 300000,
            cores = 4,
            refresh = 3000,
            thin = 4,
            init = init_fn)  # Pasar la función de inicialización
fit

options(max.print = 1050)
print(fit)
save(fit, file = "johnatan3/particion_modelo_completo_sub2.RData")
