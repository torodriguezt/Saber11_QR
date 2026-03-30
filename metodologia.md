# Regresion Cuantilica Bayesiana Jerarquica via INLA-EM

## 1. Objetivo

Estimar el efecto de covariables socioeconomicas (estrato, educacion de los
padres, acceso a internet, horas de trabajo, naturaleza del colegio) sobre
el **cuantil 0.8** (percentil 80) de los puntajes PUNT_GLOBAL del examen
SABER 11, respetando la estructura jerarquica:

```
Estudiante (i) → Colegio (j) → Municipio (k)
```

## 2. El modelo

### 2.1 Regresion cuantilica y la distribucion Laplace asimetrica (ALD)

En regresion cuantilica clasica (Koenker & Bassett, 1978), se minimiza la
funcion de perdida "check":

```
rho_tau(u) = u * (tau - I(u < 0))
```

La interpretacion bayesiana (Yu & Moyeed, 2001) utiliza la **Asymmetric
Laplace Distribution (ALD)** como likelihood:

```
y_i | mu_i, sigma, tau ~ ALD(mu_i, sigma, tau)
```

donde mu_i es el predictor lineal (intercepto + efectos fijos + efectos
aleatorios) y tau = 0.8 es el cuantil de interes.

### 2.2 Estructura jerarquica

```
mu_i = beta_0 + X_i * beta + alpha_{j[i]} + gamma_{k[i]}

alpha_j ~ Normal(0, sigma_cole^2)   (efecto aleatorio colegio)
gamma_k ~ Normal(0, sigma_muni^2)   (efecto aleatorio municipio)
```

X_i incluye las covariables del estudiante i (estrato, internet, horas de
trabajo, educacion madre, educacion padre, naturaleza del colegio).

### 2.3 Priors

- sigma_cole, sigma_muni: Penalized Complexity priors (PC priors,
  Simpson et al., 2017): P(sigma > 50) = 0.01
- beta: priors implicitos de INLA (Gaussianos difusos)

## 3. El problema computacional

La ALD no es diferenciable (tiene un "quiebre" en su moda), lo que causa
problemas para:

- **MCMC (HMC/NUTS)**: gradientes discontinuos causan divergencias o
  convergencia extremadamente lenta. Con parametrizacion centrada y 300,000
  iteraciones no se logro convergencia.
- **Variational Bayes (ADVI)**: falla en la adaptacion del step-size
  ("All proposed step-sizes failed").
- **INLA directamente**: R-INLA removio la familia ALD. El link "quantile"
  solo esta disponible para familias discretas (Poisson, Binomial).

## 4. Por que necesitamos un algoritmo EM

El modelo tiene dos bloques de cantidades desconocidas:

1. **Parametros del modelo**: beta (efectos fijos), alpha_j (efecto colegio),
   gamma_k (efecto municipio), sigma_ald, sigma_cole, sigma_muni.
2. **Variables latentes**: w_i (una por cada una de las ~551,000 observaciones).

El problema es que ninguno de los metodos estandar puede estimar ambos
bloques simultaneamente de forma eficiente:

- **INLA** es excelente para el bloque 1 (modelos Gaussianos latentes
  jerarquicos), pero necesita que la likelihood sea de una familia
  conocida. La ALD no esta disponible, y con las w_i desconocidas la
  likelihood no es Gaussiana.

- **Si conocieramos las w_i**, el modelo se reduce a un Gaussiano
  ponderado — que INLA resuelve en segundos.

- **Si conocieramos los parametros del modelo**, las w_i se pueden
  calcular con formula cerrada desde su distribucion condicional
  (una GIG).

El algoritmo EM (Expectation-Maximization, Dempster et al., 1977)
resuelve exactamente este tipo de problema: cuando hay variables
latentes (las w_i) que impiden la estimacion directa, EM alterna
entre:

- **E-step**: dado el ajuste actual de los parametros, calcular la
  esperanza de las variables latentes (actualizar w_i).
- **M-step**: dadas las w_i, maximizar la verosimilitud respecto a
  los parametros (ajustar INLA).

En nuestro caso, cada iteracion del EM:
1. Fija los w_i y le pasa a INLA un modelo Gaussiano ponderado
   (que INLA resuelve en ~20 segundos con 551K observaciones).
2. Usa los resultados de INLA para recalcular los w_i con formula
   cerrada (instantaneo).
3. Repite hasta que los parametros se estabilicen (~25 iteraciones).

Esto convierte un problema intratable (ALD jerarquica con 551K obs)
en una secuencia de problemas faciles (Gaussianos ponderados en INLA).

## 5. La solucion: Mezcla normal-exponencial + EM con INLA

### 5.1 Representacion como mezcla (Kozumi & Kobayashi, 2011)

La ALD se puede escribir como una mezcla de normales:

```
y_i | w_i ~ Normal(mu_i + delta * w_i,  sqrt(gamma * w_i))
w_i       ~ Exponencial(1 / sigma_ald)
```

donde:
- delta = sigma_ald * (1 - 2*tau) / (tau * (1 - tau))
- gamma = 2 * sigma_ald^2 / (tau * (1 - tau))

Esta representacion es **matematicamente equivalente** a la ALD, pero
condicionando en w_i el modelo es Gaussiano — exactamente lo que INLA
sabe resolver de forma eficiente.

### 5.2 Algoritmo EM

**Idea**: alternar entre estimar los parametros del modelo (dados los w_i)
y actualizar los w_i (dados los parametros).

**Inicializacion**:
- sigma_ald = 10 (estimacion conservadora)
- w_i = sigma_ald para todo i (media de la Exponencial)

**Iteracion t = 1, 2, ...**:

**E-step** (ajustar INLA):
1. Calcular delta y gamma desde sigma_ald
2. Crear respuesta ajustada: y_tilde_i = y_i - delta * w_i
3. Crear pesos de precision: scale_i = 1 / w_i
4. Ajustar modelo Gaussiano jerarquico ponderado en INLA:
   ```
   y_tilde ~ efectos_fijos + f(muni_idx, iid) + f(cole_idx, iid)
   family = "gaussian", scale = 1/w
   ```
5. Extraer valores ajustados mu_i y precision del likelihood

**M-step** (actualizar w_i):
1. Calcular residuos: r_i = y_i - mu_i
2. La posterior de w_i | r_i, sigma es una distribucion **GIG(1/2, a, b)**
   (Generalized Inverse Gaussian) con:
   ```
   a = delta^2 / gamma + 2 / sigma_ald
   b = r_i^2 / gamma
   ```
3. Actualizar w_i con la esperanza de la GIG (formula cerrada):
   ```
   E[w_i] = sqrt(b/a) * (1 + 1/sqrt(a*b))
   ```
4. Actualizar sigma_ald desde la precision estimada por INLA

**Amortiguamiento** (damping):
Para evitar oscilaciones, las actualizaciones se amortiguan:
```
sigma_ald = 0.3 * sigma_nuevo + 0.7 * sigma_anterior
w_i       = 0.3 * w_nuevo     + 0.7 * w_anterior
```

**Convergencia**: cuando el cambio relativo en sigma_ald < 0.01.

**Ajuste final**: una vez convergido, se corre INLA una ultima vez sin
simplificaciones para obtener las distribuciones posteriores completas.

### 5.3 Por que funciona

- Cada iteracion de INLA resuelve un modelo Gaussiano jerarquico con
  pesos de precision heterogeneos — algo para lo que INLA esta optimizado
  (segundos con 551K observaciones).
- Los w_i se actualizan con formula cerrada (no requieren sampling).
- El algoritmo EM converge en ~25 iteraciones (~10 minutos total).
- La estructura jerarquica (municipio -> colegio -> estudiante) se
  preserva intacta en cada iteracion.

## 6. Tipos de contrastes

### 6.1 Treatment (nivel de referencia)

```
estrato2 = efecto de estrato 2 RELATIVO a estrato 1
```

Util para comparaciones directas contra una categoria base.

### 6.2 Sum-to-zero (ANOVA bayesiana)

```
estrato_k = desviacion del estrato k respecto a la MEDIA GLOBAL
sum(estrato_k) = 0
```

Util para ANOVA bayesiana: cada efecto se interpreta como cuanto se devia
esa categoria de la media. No hay nivel de referencia privilegiado.

Equivalente al centrado que se usaba en el modelo Stan original:
`estrato_cen = estrato_effect - mean(estrato_effect)`.

## 7. Comparacion de especificaciones

| Modelo | Proposito |
|--------|-----------|
| Completo | Modelo principal con todas las covariables |
| Sin educ padres | Ver si estrato cambia de signo (confounding) |
| Sin estrato | Ver efecto "limpio" de educacion sin confounding |

Si el efecto de estrato cambia de signo entre "Completo" y "Sin educ
padres", confirma confounding: estrato y educacion de padres estan
fuertemente correlacionados, y al incluir ambos, estrato pierde su
efecto marginal porque la educacion lo "absorbe".

## 8. Limitaciones

1. **Aproximacion EM**: no propaga la incertidumbre de w_i a los
   parametros del modelo. Los intervalos de credibilidad pueden ser
   ligeramente optimistas.

2. **Mezcla ALD como likelihood**: Sriram et al. (2013) muestran que
   la ALD como working likelihood subestima la varianza posterior.
   Los intervalos de credibilidad son orientativos, no exactos.

3. **Amortiguamiento**: el factor de damping (0.3) es heuristico.
   Valores menores dan convergencia mas suave pero mas lenta.

4. **Homocedasticidad de los efectos aleatorios**: el modelo asume
   que la varianza entre colegios y municipios es la misma en todos
   los cuantiles. Si la varianza cambia con el cuantil, esto no se
   captura.

## 9. Referencias

### Regresion cuantilica bayesiana

- **Yu, K., & Moyeed, R. A. (2001)**. Bayesian quantile regression.
  *Statistics & Probability Letters*, 54(4), 437-447.
  La formulacion bayesiana original usando la ALD como likelihood.

- **Koenker, R., & Bassett, G. (1978)**. Regression quantiles.
  *Econometrica*, 46(1), 33-50.
  El paper fundacional de regresion cuantilica.

### Mezcla normal-exponencial de la ALD

- **Kozumi, H., & Kobayashi, G. (2011)**. Gibbs sampling methods for
  Bayesian quantile regression. *Journal of Statistical Computation
  and Simulation*, 81(11), 1565-1578.
  La representacion de la ALD como mezcla que usamos en el algoritmo EM.

- **Reed, C., & Yu, K. (2011)**. An efficient Gibbs sampler for
  Bayesian quantile regression. Technical report.
  Extension del Gibbs sampler basado en la mezcla.

### Critica a la ALD como likelihood

- **Sriram, K., Ramamoorthi, R. V., & Ghosh, P. (2013)**. Posterior
  consistency of Bayesian quantile regression based on the misspecified
  asymmetric Laplace density. *Bayesian Analysis*, 8(4), 879-898.
  Muestra que la ALD subestima la varianza posterior.

- **Yang, Y., Wang, H. J., & He, X. (2016)**. Posterior inference in
  Bayesian quantile regression with asymmetric Laplace likelihood.
  *International Statistical Review*, 84(3), 327-344.

### INLA y regresion cuantilica

- **Rue, H., Martino, S., & Chopin, N. (2009)**. Approximate Bayesian
  inference for latent Gaussian models by using integrated nested Laplace
  approximations. *Journal of the Royal Statistical Society: Series B*,
  71(2), 319-392.
  El paper fundamental de INLA.

- **Yue, Y. R., & Rue, H. (2011)**. Bayesian inference for additive
  mixed quantile regression models. *Computational Statistics & Data
  Analysis*, 55(1), 84-96.
  EM con INLA para regresion cuantilica — la base directa de nuestro
  algoritmo.

- **Padellini, T., & Rue, H. (2019)**. Model-aware quantile regression
  for discrete data. arXiv:1804.03714.
  Propone el quantile link como alternativa a la ALD.

### Priors

- **Simpson, D., Rue, H., Riebler, A., Martins, T. G., & Sorbye, S. H.
  (2017)**. Penalising model component complexity: A principled, practical
  approach to constructing priors. *Statistical Science*, 32(1), 1-28.
  Los PC priors usados para los efectos aleatorios.

### Modelos jerarquicos y parametrizacion

- **Betancourt, M., & Girolami, M. (2015)**. Hamiltonian Monte Carlo
  for hierarchical models. In *Current Trends in Bayesian Methodology
  with Applications*. CRC Press.
  Explica el problema del "funnel" y la non-centered parameterization.

- **Papaspiliopoulos, O., Roberts, G. O., & Skold, M. (2007)**. A
  general framework for the parametrization of hierarchical models.
  *Statistical Science*, 22(1), 59-73.
  Marco teorico para elegir entre parametrizacion centrada y no centrada.
