# Regresion Cuantilica Bayesiana Jerarquica - SABER 11

Estimacion del efecto de covariables socioeconomicas sobre el **percentil 80** del puntaje global del examen SABER 11 (Colombia, 2023-2), usando regresion cuantilica bayesiana jerarquica con estructura estudiante -> colegio -> municipio.

## Metodo

Se utiliza un algoritmo **EM con INLA** basado en la representacion de la Asymmetric Laplace Distribution (ALD) como mezcla normal-exponencial (Kozumi & Kobayashi, 2011; Yue & Rue, 2011):

1. **E-step**: dado los pesos latentes w_i, INLA ajusta un modelo Gaussiano ponderado jerarquico (~20s con 415K obs)
2. **M-step**: dados los parametros, se actualizan los w_i con formula cerrada desde la distribucion GIG
3. Converge en ~20 iteraciones (~10 min total)

Ver [metodologia.md](metodologia.md) para la descripcion completa del algoritmo, modelo y referencias.

## Estructura del proyecto

```
.
├── correr_modelos.R       # Script principal: graficos descriptivos + modelo INLA-EM
├── inla_ald.R             # Script standalone del algoritmo EM-INLA (una especificacion)
├── metodologia.md         # Documentacion del metodo, modelo y referencias
├── figuras/               # Graficos generados por correr_modelos.R
└── archivo/               # Intentos previos (Stan MCMC, Variational Bayes) - no se usan
```

## Datos

El dataset `SB11_20232.TXT` contiene los resultados del examen SABER 11 del segundo semestre de 2023, publicados por el ICFES. No se incluye en el repositorio por su tamano (~450MB). Debe colocarse en la raiz del proyecto.

**Fuente**: [ICFES - Datos Abiertos](https://www.icfes.gov.co/)

## Covariables

| Variable | Descripcion |
|----------|-------------|
| `estrato` | Estrato socioeconomico (1-6) |
| `educ_madre` | Educacion de la madre (Ninguno, Primaria, Secundaria, Tecnica, Profesional, Postgrado) |
| `educ_padre` | Educacion del padre (mismos niveles) |
| `internet` | Acceso a internet en el hogar (Si/No) |
| `horas` | Trabaja durante la semana (Trabaja/No_trabaja) |
| `naturaleza` | Naturaleza del colegio (Oficial/No oficial) |

## Efectos aleatorios

- **Municipio** (`muni_idx`): intercepto aleatorio por municipio de ubicacion del colegio
- **Colegio** (`cole_idx`): intercepto aleatorio por colegio

Priors: PC priors con P(sigma > 50) = 0.01.

## Como correr

```r
# Requiere R-INLA instalado:
# install.packages("INLA",
#   repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"),
#   dep = TRUE)

# Correr modelo completo + graficos
source("correr_modelos.R")
```

Genera:
- Graficos descriptivos y de resultados en `figuras/`
- Resultados del modelo en `resultados_modelo_completo.RData`

## Contrastes

Se ajustan dos versiones del modelo:
- **Sum-to-zero**: cada coeficiente es la desviacion respecto a la media global (ANOVA bayesiana)
- **Treatment**: cada coeficiente es la diferencia respecto a la categoria de referencia

## Referencias principales

- Kozumi, H. & Kobayashi, G. (2011). Gibbs sampling methods for Bayesian quantile regression. *JSCS*.
- Yue, Y. R. & Rue, H. (2011). Bayesian inference for additive mixed quantile regression models. *CSDA*.
- Rue, H., Martino, S. & Chopin, N. (2009). Approximate Bayesian inference for latent Gaussian models using INLA. *JRSS-B*.
- Yu, K. & Moyeed, R. A. (2001). Bayesian quantile regression. *Statistics & Probability Letters*.
