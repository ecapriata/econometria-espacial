rm(list=ls())

# Ejemplo con datos de fecundidad departamental en Argentina, 2001

# Declaramos el directorio de trabajo
setwd("C:/Users/estef/Desktop/San Andrés/2022/3er trimestre/Econometría Espacial/HW2-20221124T153945Z-001/HW2")

# cargamos las librerias espaciales
library(R.matlab)
library(spdep)
library(spatialreg)
library(rgdal)
options("rgdal_show_exportToProj4_warnings"="none")

# leemos el shapefile
data <- readMat("C:/Users/estef/Desktop/San Andrés/2022/3er trimestre/Econometría Espacial/HW2-20221124T153945Z-001/HW2/mystery_process.mat")
names(data)

w <- data$W
class(w)

library(spdep)
w_lista <- mat2listw(w, style = "W")
summary(w_lista)

# regresión no-espacial
# Comenzamos en el modelo 8.
reg.ols <- lm(y ~ x1 + x2, data=data)
summary(reg.ols)

# Test I de Moran
lm.morantest(reg.ols,w_lista, alternative = "two.sided")
# Encontramos estructura espacial.

# Tests LM's
lms <- lm.LMtests(reg.ols, w_lista, test = "all")
tests <- t(sapply(lms, function(x) c(x$statistic, x$parameter, x$p.value)))
colnames(tests) <- c("Test","df","p-valor")
printCoefmat(tests)
# Agrupamos los tests.
# Trabajando al 5%. Todos los tests son significativos. Si los dos me dan significativo paso directamente a un SARAR (M2), pero es difícil confiar en sus estimaciones.

# los resultados encuentran autocorrelación espacial: estimar modelo espacial

# 1. SLX
# Hacemos un spacial lag.
w_x1 <- lag.listw(w_lista, data$x1)
w_x2 <- lag.listw(w_lista, data$x2)
reg.slx <- lm(y ~ x1 + x2 + w_x1 + w_x2, data=data)
summary(reg.slx)
# Los lags son significativos.

lm.morantest(reg.slx,w_lista, alternative = "two.sided")
# El i de Moran me da significativo (tengo problema de especificación).

lms <- lm.LMtests(reg.slx, w_lista, test = "all")
tests <- t(sapply(lms, function(x) c(x$statistic, x$parameter, x$p.value)))
colnames(tests) <- c("Test","df","p-valor")
printCoefmat(tests)
# Leo los LM robustos, porque los simples son casi siempre significativos. El RLMlag no da significativo. Y el RLMerr da significativo (tipo SEM).

# 2. SEM
reg.sem <- spatialreg::errorsarlm(y ~ x1 + x2, data=data, w_lista)
summary(reg.sem)
(LR_lambda=2*(reg.sem$LL - logLik(reg.ols)))

# 3. SLM
reg.slm <- spatialreg::lagsarlm(y ~ x1 + x2, data=data, w_lista)
summary(reg.slm)
(LR_rho=2*(reg.slm$LL - logLik(reg.ols)))
# El LR es equivalente a la diferencia que estamos calculando.

# 4. SDM
# Puedo estimar usando lagsar. Es una estimación de un SLM pero en las x tengo las Wx.
# si en lasarlm ponemos Durbin, podemos elegir qué variables rezagar. El "mixed" es el tipo de modelo.
reg.sdm <- spatialreg::lagsarlm(y ~ x1 + x2, Durbin=~x1 + x2,data=data, w_lista, type = "mixed")
summary(reg.sdm)

# Partimos de arriba. Las dos estrategias usualmente tienen que llevar al mismo modelo final. La estrategia a la "Henry" es más robusta.
(LR_rho=2*(reg.sdm$LL - logLik(reg.slx)))
LR.Sarlm(reg.sdm,reg.slx)

# test de factores comunes
LR.Sarlm(reg.sdm,reg.sem)

# Test entre SDM vs SLM
LR.Sarlm(reg.sdm,reg.slm)

# 5. SDEM
reg.sdem <- spatialreg::errorsarlm(y ~ x1 + x2, Durbin=~x1 + x2,data=data, w_lista, etype = "emixed")
summary(reg.sdem)
(LR_lambda1=2*(reg.sdem$LL - logLik(reg.slx)))
# Usa Durbin al final para ver qué variables rezagar.
# El akaike (AIC) sale por defecto en todos los modelos.VER EL VALOR MÁS PRÓXIMO A CERO.

LR.Sarlm(reg.sdem,reg.slx)
LR.Sarlm(reg.sdem,reg.sem)
# se reduce a un SEM

# 6. SARAR model. ESTIMARLO CUANDO LOS DOS LM ROBUSTOS SEAN SIGNIFICATIVOS.
reg.sarar <- spatialreg::sacsarlm(y ~ x1 + x2, data=data, w_lista)
summary(reg.sarar)
(LR_lambda2=2*(reg.sarar$LL - reg.slm$LL))
(LR_rho2=2*(reg.sarar$LL - reg.sdem$LL))
LR.Sarlm(reg.sarar,reg.sem)
LR.Sarlm(reg.sarar,reg.slm)
