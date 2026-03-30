library(rstan)


# Leer el archivo CSV
datos_sub <- read.csv("johnatan3/SUB_BASE4.csv", header = TRUE, sep = ",")

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
  tau         = 0.3,
  col2        = datos_sub_distinct$Colegio,
  muni2       = datos_sub_distinct$ClusterMuni
)


init_fn <- function() {
  list(
    # Inicialización de los parámetros alpha con sus medias posteriores
    alpha = c(200.67, 246.39, 286.98, 223.01, 211.19, 188.29, 227.40, 
              190.63, 178.06, 222.21, 205.34, 238.42, 326.52, 276.54),
    
    # Inicialización de los parámetros mu con sus medias posteriores
    mu = c(233.06, 223.20, 219.39, 240.88),
    
    # Inicialización de mu_global con su media posterior
    mu_global = 229.54,
    
    # Inicialización de los parámetros sigma con sus medias posteriores
    sigma = 14.20,
    sigma_cole = 40.35,
    sigma_global = 18.19,
    
    # Inicialización de las covariables
    horas_cen = c(4.66, -4.66),  # Asegúrate de que estas sean las medias de tus covariables
    estrato_cen = c(8.35, 7.14, 4.75, -0.38, -5.59, -14.27),  # Ejemplo de valores
    internet_cen = c(-2.72, 2.72),  # Valores para la covariable de internet
    educ_madre_cen = c(-12.84, -6.74, -3.34, 4.90, 6.63, 11.39),  # Media de las covariables de educación madre
    educ_padre_cen = c(-15.14, -9.71, -4.73, 2.46, 3.54, 23.57),  # Media de las covariables de educación padre
    naturaleza_cen = c(-10.70, 10.70)  # Para las covariables de naturaleza, si son dos, agrega las medias
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
save(fit, file = "johnatan3/particion03_modelo_completo_sub4.RData")