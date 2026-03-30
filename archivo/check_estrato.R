###############################################################################
# Verificacion rapida: efecto de estrato con y sin controlar por educacion
# Si estrato se invierte al controlar por educacion -> confounding
###############################################################################

library(INLA)

# Cargar datos ya limpios del modelo INLA (si esta en memoria, si no recargar)
# Asumimos que 'datos' ya esta en el environment con las recodificaciones

# Modelo 1: solo estrato (sin educacion de padres)
cat(">>> Modelo 1: y ~ estrato + internet + horas + naturaleza + RE\n")
fit1 <- inla(
  y ~ 1 + estrato + internet + horas + naturaleza +
    f(muni_idx, model = "iid") + f(cole_idx, model = "iid"),
  family = "gaussian", data = datos, verbose = FALSE
)
cat("\nEfectos de estrato SIN controlar por educacion:\n")
print(round(fit1$summary.fixed[grep("estrato", rownames(fit1$summary.fixed)),
                                 c("mean", "0.025quant", "0.975quant")], 2))

# Modelo 2: estrato + educacion (modelo completo)
cat("\n>>> Modelo 2: y ~ estrato + internet + horas + educ_madre + educ_padre + naturaleza + RE\n")
fit2 <- inla(
  y ~ 1 + estrato + internet + horas + educ_madre + educ_padre + naturaleza +
    f(muni_idx, model = "iid") + f(cole_idx, model = "iid"),
  family = "gaussian", data = datos, verbose = FALSE
)
cat("\nEfectos de estrato CONTROLANDO por educacion:\n")
print(round(fit2$summary.fixed[grep("estrato", rownames(fit2$summary.fixed)),
                                 c("mean", "0.025quant", "0.975quant")], 2))

cat("\nSi los signos se invierten entre Modelo 1 y 2 -> confounding confirmado\n")
