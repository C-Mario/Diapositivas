
library(rucm)

#Para el modelo de nivel local

#| output-location: fragment
Nilo <- datasets::Nile
# ajustar modelo 
modelNilo <- rucm::ucm(formula = Nile~0, data = Nile, level = TRUE)
#
plot(Nilo)
lines(modelNilo$s.level, col = "cyan",lwd=1.5)

#Para el modelo de tendencia local lineal

modelNilo2 <- rucm::ucm(formula=Nile~0, data = Nile,level = TRUE)







