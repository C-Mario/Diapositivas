rm(list = ls())
graphics.off()

library(tsm)
library(dlm)
library(dplyr)
library(zoo)
library(forecast)

#Leer y transformar los datos

## EXPORTACIONES----

#Cargar datos
load(file = "_environment.RData")

#Ajuste del formato de los datos
options(scipen = 12)
Exportaciones$Total <- round(Exportaciones$Total)
Exportaciones$Mes <- as.yearmon(Exportaciones$Mes)

# Cambiar la ventana de datos
Exportaciones <- Exportaciones[97:378,]

# Crear dataframe
exportaciones <- data.frame(as.Date(Exportaciones$Mes),Exportaciones$Total)
colnames(exportaciones) <- c("Fecha", "Dato")

#Crear objeto de tipo ts indicandole la fecha de inicio y la frecuencia
exportaciones_ts <- ts(exportaciones$Dato,start = c(2000,01),frequency = 12)

#Transformación logarítmica con lambda = 0.45
lexportaciones_ts <- (1/0.45)*((exportaciones_ts^(0.45))-1)

plot(lexportaciones_ts);inf<-lexportaciones_ts
#inf<-diff(lexportaciones_ts, lag=1);plot(inf)

## QUEJAS EQUIPAJE ####

quejas<- read.csv("baggagecomplaints.csv")
quejas<- quejas %>% filter(Airline=="American Eagle")
quejas<- data.frame(quejas$Date,as.numeric(quejas$Baggage))
colnames(quejas)<-c("Fecha","Total")

quejas_ts<-ts(quejas$Total,start=c(2004,01),frequency=12)

plot(quejas_ts);inf<-quejas_ts


### Dividir en entrenamiento y prueba ####

lserie <- length(inf)
ntrain <- trunc(lserie*0.8)
train <- window(inf, end = time(inf)[ntrain])
test <- window(inf, start = time(inf)[ntrain]+1/12) 
ntest <- length(test)
paste("Número de datos en el conjunto de entrenamiento:", ntrain)
paste("Número de datos en el conjunto de entrenamiento:", ntest)

## QUEJAS EQUIPAJE ####

#### LOCAL LEVEL MODEL ####

# Incluye una pendiente estocástica, en el paquete glm
#  se construyen estos modelos de acuerdo con el orden del
#  polinomio, donde un polinomio de primer orden es un modelo
#  de nivel local.
#  
#  Dos parámetros para error, uno en la ecuación de medición y
#  otro en la ecuación de estado.

fn <- function(parm) {
  dlmModPoly(order = 1, dV = exp(parm[1]), dW = exp(parm[2]))
}

# Se asumen los valores iniciales de los parámetros como cero
# y se procede a ajustar el modelo con métodos de máxima verosim

fit <- dlmMLE(train, rep(0, 2), build = fn, hessian = TRUE)
(conv <- fit$convergence) # Hay convergencia

# Función de verosimilitud y criterios de información

loglik <- dlmLL(train, dlmModPoly(1))
n.coef <- 2
r.aic <- (2 * (loglik)) + 2 * (sum(n.coef))  #dlmLL caculates the neg. LL
r.bic <- (2 * (loglik)) + (log(length(train))) * (n.coef)

# Variance - Covariance Matrix

mod <- fn(fit$par)
obs.error.var <- V(mod)
state.error.var <- W(mod)

# Filtro de Kalman y suavizado

filtered <- dlmFilter(train, mod = mod)
smoothed <- dlmSmooth(filtered)

# Resultados del filtro de Kalman

resids <- residuals(filtered, sd = FALSE)
mu <- dropFirst(smoothed$s)
mu.1 <- mu[1]
mu.end <- mu[length(mu)]

# Tendencia estocástica Ajustada y resuduales

par(mfrow = c(2, 1), mar = c(2.2, 2.2, 1, 1), cex = 0.8)
plot.ts(train, col = "darkgrey", xlab = "", ylab = "", lwd = 1.5)
lines(mu, col = "black")
legend("topright", legend = c("Observed Deflator", "Stochastic level"), 
       lwd = c(2, 1), col = c("darkgrey", "black"), bty = "n")

plot.ts(resids, ylab = "", xlab = "", col = "darkgrey", 
        lwd = 1.5)
abline(h = 0)
legend("topright", legend = "Residuals", lwd = 1.5, col = "darkgrey", 
       bty = "n")

# Criterios de información y varianza del error

cat("AIC", r.aic)

cat("BIC", r.bic)

cat("V.variance", obs.error.var)

cat("W.variance", state.error.var)

# Diagnóstico de los residuales

par(mfrow = c(1, 1), mar = c(2.2, 2.2, 1, 1), cex = 0.8)
hist(resids, prob = TRUE, col = "grey", main = "", breaks = seq(-12,6, length.out = 50)
)

ac(resids)  # acf

Box.test(resids, lag = 12, type = "Box-Pierce", fitdf = 2)  # joint autocorrelation

#Test de normalidad
tseries::jarque.bera.test(resids)

### Pronostico

alpha <- mu

comb.state <- alpha

forecast <- dlmForecast(filtered, nAhead = 48)
var.2 <- unlist(forecast$Q)
wid.2 <- qnorm(0.05, lower = FALSE) * sqrt(var.2)
comb.fore <- forecast$f

result <- ts(c(comb.state, comb.fore), start = c(2017,1), frequency = 12)


par(mfrow = c(1, 1), mar = c(2.2, 2.2, 1, 1), cex = 0.8)
plot.ts(inf[1:273], col = "darkgrey", plot.type = "single", 
        xlab = "", ylab = "", lty = 1.5)
lines(result[1:273], col = "red", lwd = 0.8)
abline(v = 225, col = "blue", lwd = 1, lty = 3)
legend("topleft", legend = c("Observed Deflator", "Stochastic level"), 
       lwd = c(1.5, 1), col = c("red", "darkgrey"), bty = "n")

### Intervalos de confianza---------

conf.tmp <- unlist(dlmSvd2var(smoothed$U.S, smoothed$D.S))
conf <- ts(as.numeric(conf.tmp)[-1], start = c(2000, 1), 
           frequency = 12)
wid <- qnorm(0.05, lower = FALSE) * sqrt(conf)

conf.pos <- mu + wid
conf.neg <- mu - wid

mu.f <- dropFirst(filtered$a)
cov.tmp <- unlist(dlmSvd2var(filtered$U.R, filtered$D.R))

if (sum(dim(mod$FF)) == 2) {
  variance <- cov.tmp + as.numeric(V(mod))
} else {
  variance <- (sapply(cov.tmp, function(x) mod$FF %*% 
                        x %*% t(mod$FF))) + V(mod)
}

par(mfrow = c(1, 1), mar = c(2.2, 2.2, 1, 1), cex = 0.8)
plot.ts(train, col = "darkgrey", xlab = "", ylab = "", lwd = 1.5)
lines(mu, col = "black")
lines(conf.pos, col = "red")
lines(conf.neg, col = "red")
legend("topright", legend = c("Observed Deflator", "Stochastic level", 
                              "Confidence Interval"), lwd = c(1.5, 1, 1), col = c("darkgrey", 
                                                                                  "black", "red"), bty = "n")

### Pronostico-----------

comb.state <- cbind(mu, conf.pos, conf.neg)

forecast <- dlmForecast(filtered, nAhead = 12)
var.2 <- unlist(forecast$Q)
wid.2 <- qnorm(0.05, lower = FALSE) * sqrt(var.2)
comb.fore <- cbind(forecast$f, forecast$f + wid.2, forecast$f - 
                     wid.2)

result <- ts(rbind(comb.state, comb.fore), start = c(2000,1), frequency = 12)

par(mfrow = c(1, 1), mar = c(2.2, 2.2, 1, 1), cex = 0.8)
plot.ts(result, col = c("black", "red", "red"), plot.type = "single", 
        xlab = "", ylab = "", lty = c(1, 2, 2))
lines(train, col = "darkgrey", lwd = 1.5)
legend("topleft", legend = c("Observed Deflator", "Stochastic level"), 
       lwd = c(1.5, 1), col = c("darkgrey", "black"), bty = "n")


### Modelo de tendencia lineal local ####

fn <- function(parm) {
  dlmModPoly(order = 2, dV = exp(parm[1]), dW = exp(parm[2:3]))
}

fit <- dlmMLE(train, rep(0, 3), build = fn, hessian = TRUE)
conv <- fit$convergence  # zero for converged

loglik <- dlmLL(train, dlmModPoly(2))
n.coef <- 3
r.aic <- (2 * (loglik)) + 2 * (sum(n.coef))  #dlmLL caculates the neg. LL
r.bic <- (2 * (loglik)) + (log(length(train))) * (n.coef)

mod <- fn(fit$par)
obs.error.var <- V(mod)
state.error.var <- diag(W(mod))

filtered <- dlmFilter(train, mod = mod)
smoothed <- dlmSmooth(filtered)

resids <- residuals(filtered, sd = FALSE)

mu <- dropFirst(smoothed$s[, 1])
upsilon <- dropFirst(smoothed$s[, 2])
mu.1 <- mu[1]
mu.end <- mu[length(mu)]
ups.1 <- upsilon[1]
ups.end <- upsilon[length(mu)]

par(mfrow = c(3, 1), mar = c(2.2, 2.2, 1, 1), cex = 0.8)
plot.ts(train, col = "darkgrey", xlab = "", ylab = "", lwd = 2)
lines(mu, col = "black")
legend("topright", legend = c("Observed Deflator", "Stochastic level"), 
       lwd = c(2, 1), col = c("darkgrey", "black"), bty = "n")

plot.ts(upsilon, col = "darkgrey", xlab = "", ylab = "", 
        lwd = 2)
legend("topright", legend = "Slope", lwd = 2, col = "darkgrey", 
       bty = "n")

plot.ts(resids, ylab = "", xlab = "", col = "darkgrey", 
        lwd = 2)
abline(h = 0)
legend("topright", legend = "Residuals", lwd = 2, col = "darkgrey", 
       bty = "n")

cat("AIC", r.aic)

cat("BIC", r.bic)

cat("V.variance", obs.error.var)

cat("W.variance", state.error.var)


par(mfrow = c(1, 1), mar = c(2.2, 2.2, 1, 1), cex = 0.8)
hist(resids, prob = TRUE, col = "grey", main = "", breaks = seq(-12,6, length.out = 50)
)

ac(resids)  # acf

Box.test(resids, lag = 12, type = "Box-Pierce", fitdf = 2)  # joint autocorrelation

#Test de normalidad
tseries::jarque.bera.test(resids)

### Pronostico-------------

alpha <- mu+upsilon

comb.state <- alpha

forecast <- dlmForecast(filtered, nAhead = 48)
var.2 <- unlist(forecast$Q)
wid.2 <- qnorm(0.05, lower = FALSE) * sqrt(var.2)
comb.fore <- forecast$f

result <- ts(c(comb.state, comb.fore), start = c(2017,1), frequency = 12)


par(mfrow = c(1, 1), mar = c(2.2, 2.2, 1, 1), cex = 0.8)
plot.ts(inf[1:273], col = "darkgrey", plot.type = "single", 
        xlab = "", ylab = "", lty = 1.5)
lines(result[1:273], col = "red", lwd = 0.8)
abline(v = 225, col = "blue", lwd = 1, lty = 3)
legend("topleft", legend = c("Observed Deflator", "Stochastic level"), 
       lwd = c(1.5, 1), col = c("red", "darkgrey"), bty = "n")

### Modelo de nivel local con componente estacional ####

fn <- function(parm) {
  mod <- dlmModPoly(order = 1) + dlmModSeas(frequency = 12)
  V(mod) <- exp(parm[1])
  diag(W(mod))[1:2] <- exp(parm[2:3])
  return(mod)
}

fit <- dlmMLE(train, rep(0, 3), build = fn, hessian = TRUE)
conv <- fit$convergence  # zero for converged

loglik <- dlmLL(train, dlmModPoly(1) + dlmModSeas(12))
n.coef <- 3
r.aic <- (2 * (loglik)) + 2 * (sum(n.coef))  #dlmLL caculates the neg. LL
r.bic <- (2 * (loglik)) + (log(length(train))) * (n.coef)

mod <- fn(fit$par)
obs.error.var <- V(mod)
state.error.var <- diag(W(mod))

filtered <- dlmFilter(train, mod = mod)
smoothed <- dlmSmooth(filtered)
resids <- residuals(filtered, sd = FALSE)
mu <- dropFirst(smoothed$s[, 1])
gammas <- dropFirst(smoothed$s[, 2])
mu.1 <- mu[1]
mu.end <- mu[length(mu)]
gammas.1 <- gammas[1]
gammas.end <- gammas[length(mu)]


par(mfrow = c(3, 1), mar = c(2.2, 2.2, 1, 1), cex = 0.8)
plot.ts(train, col = "darkgrey", xlab = "", ylab = "", lwd = 2)
lines(mu, col = "black")
legend("topright", legend = c("Observed Deflator", "Stochastic level"), 
       lwd = c(2, 1), col = c("darkgrey", "black"), bty = "n")

plot.ts(gammas, col = "darkgrey", xlab = "", ylab = "", 
        lwd = 2)
legend("topright", legend = "Seasonal", lwd = 2, col = "darkgrey", 
       bty = "n")

plot.ts(resids, ylab = "", xlab = "", col = "darkgrey", 
        lwd = 2)
abline(h = 0)
legend("topright", legend = "Residuals", lwd = 2, col = "darkgrey", 
       bty = "n")


### Ajuste de las componentes combinadas

par(mfrow = c(1, 1), mar = c(2.2, 2.2, 1, 1), cex = 0.8)
alpha <- mu + gammas
par(mfrow = c(1, 1), mar = c(2.2, 2.2, 1, 1), cex = 0.8)
plot.ts(train, col = "darkgrey", xlab = "", ylab = "", lwd = 2)
lines(alpha, col = "black")
legend("topright", legend = c("Observed Deflator", "State Components"), 
       lwd = c(2, 1), col = c("darkgrey", "black"), bty = "n")

cat("AIC", r.aic)

cat("BIC", r.bic)

cat("V.variance", obs.error.var)

cat("W.variance", state.error.var)

par(mfrow = c(1, 1), mar = c(2.2, 2.2, 1, 1), cex = 0.8)
hist(resids, prob = TRUE, col = "grey", main = "", breaks = seq(-12,6, length.out = 50)
)

ac(resids)  # acf

Box.test(resids, lag = 12, type = "Box-Pierce", fitdf = 2)  # joint autocorrelation

#Test de normalidad
tseries::jarque.bera.test(resids)


### Pronostico

alpha <- mu + gammas

comb.state <- alpha

forecast <- dlmForecast(filtered, nAhead = 48)
var.2 <- unlist(forecast$Q)
wid.2 <- qnorm(0.05, lower = FALSE) * sqrt(var.2)
comb.fore <- forecast$f

result <- ts(c(comb.state, comb.fore), start = c(2017,1), frequency = 12)


par(mfrow = c(1, 1), mar = c(2.2, 2.2, 1, 1), cex = 0.8)
plot.ts(inf[1:273], col = "darkgrey", plot.type = "single", 
        xlab = "", ylab = "", lty = 1.5)
lines(result[1:273], col = "red", lwd = 0.8)
abline(v = 225, col = "blue", lwd = 1, lty = 3)
legend("topleft", legend = c("Observed Deflator", "Stochastic level"), 
       lwd = c(1.5, 1), col = c("red", "darkgrey"), bty = "n")


### Modelo de nivel local con Tendencia y  componente estacional ####

fn <- function(parm) {
  mod <- dlmModPoly(order = 2) + dlmModSeas(frequency = 12)
  V(mod) <- exp(parm[1])
  diag(W(mod))[1:2] <- exp(parm[2:3])
  return(mod)
}

fit <- dlmMLE(train, rep(0, 3), build = fn, hessian = TRUE)
conv <- fit$convergence  # zero for converged

loglik <- dlmLL(train, dlmModPoly(1) + dlmModSeas(12))
n.coef <- 3
r.aic <- (2 * (loglik)) + 2 * (sum(n.coef))  #dlmLL caculates the neg. LL
r.bic <- (2 * (loglik)) + (log(length(train))) * (n.coef)

mod <- fn(fit$par)

obs.error.var <- V(mod)
state.error.var <- diag(W(mod))

filtered <- dlmFilter(train, mod = mod)
smoothed <- dlmSmooth(filtered)

resids <- residuals(filtered, sd = FALSE)
mu <- dropFirst(smoothed$s[, 1])
upsilon <- dropFirst(smoothed$s[, 2])
gammas <- dropFirst(smoothed$s[, 3])
mu.1 <- mu[1]
mu.end <- mu[length(mu)]
upsilon.1 <- upsilon[1]
upsilon.end <- upsilon[length(mu)]
gammas.1 <- gammas[1]
gammas.end <- gammas[length(mu)]

par(mfrow = c(4, 1), mar = c(2.2, 2.2, 1, 1), cex = 0.8)
plot.ts(train, col = "darkgrey", xlab = "", ylab = "", lwd = 2)
lines(mu, col = "black")
legend("topright", legend = c("Observed Deflator", "Stochastic level"), 
       lwd = c(2, 1), col = c("darkgrey", "black"), bty = "n")

plot.ts(upsilon, col = "darkgrey", xlab = "", ylab = "", 
        lwd = 2)
legend("topright", legend = "Slope", lwd = 2, col = "darkgrey", 
       bty = "n")

plot.ts(gammas, col = "darkgrey", xlab = "", ylab = "", 
        lwd = 2)
legend("topright", legend = "Seasonal", lwd = 2, col = "darkgrey", 
       bty = "n")

plot.ts(resids, ylab = "", xlab = "", col = "darkgrey", 
        lwd = 2)
abline(h = 0)
legend("topright", legend = "Residuals", lwd = 2, col = "darkgrey", 
       bty = "n")


### Ajuste de las componentes combinadas

par(mfrow = c(1, 1), mar = c(2.2, 2.2, 1, 1), cex = 0.8)
alpha <- mu + upsilon + gammas
par(mfrow = c(1, 1), mar = c(2.2, 2.2, 1, 1), cex = 0.8)
plot.ts(train, col = "darkgrey", xlab = "", ylab = "", lwd = 2)
lines(alpha, col = "black")
legend("topright", legend = c("Observed Deflator", "State Components"), 
       lwd = c(2, 1), col = c("darkgrey", "black"), bty = "n")

cat("AIC", r.aic)

cat("BIC", r.bic)

cat("V.variance", obs.error.var)

cat("W.variance", state.error.var)

par(mfrow = c(1, 1), mar = c(2.2, 2.2, 1, 1), cex = 0.8)
hist(resids, prob = TRUE, col = "grey", main = "", breaks = seq(-12,6, length.out = 50)
)

ac(resids)  # acf

Box.test(resids, lag = 12, type = "Ljung", fitdf = 2)  # joint autocorrelation

shapiro.test(resids)  # normality

### Pronostico

alpha <- mu + upsilon + gammas

comb.state <- alpha

forecast <- dlmForecast(filtered, nAhead = 48)

var.2 <- unlist(forecast$Q)
wid.2 <- qnorm(0.05, lower = FALSE) * sqrt(var.2)
comb.fore <- forecast$f

result <- ts(c(comb.state, comb.fore), start = c(2017,1), frequency = 12)


par(mfrow = c(1, 1), mar = c(2.2, 2.2, 1, 1), cex = 0.8)
plot.ts(inf[1:273], col = "darkgrey", plot.type = "single", 
        xlab = "", ylab = "", lty = 1.5)
lines(result[1:273], col = "red", lwd = 0.8)
abline(v = 225, col = "blue", lwd = 1, lty = 3)
legend("topleft", legend = c("Observed Deflator", "Stochastic level"), 
       lwd = c(1.5, 1), col = c("red", "darkgrey"), bty = "n")


### ROLLING ####

h <- 1
n <- ntest - h + 1

fc <- ts(numeric(n), start=c(2018,10), freq=12)

for(i in 1:n){  
  x <- window(lexportaciones_ts, end=c(2018, 9+(i-1)))
  filtered <- dlmFilter(x, mod = mod)
  fc[i] <-dlmForecast(filtered, nAhead = h)$f
}

test1<-InvBoxCox(test,lambda=0.45)
fc1<-InvBoxCox(fc,lambda=0.45)

dife=(test1-fc1)^2
ecm=(1/(ntest))*sum(dife)
ecm
recm1pasoModelo1=sqrt(ecm)
recm1pasoModelo1

result <- ts(c(comb.state, fc), start = c(2000,1), frequency = 12)

result<-InvBoxCox(result, lambda=0.45)

par(mfrow = c(1, 1), mar = c(2.2, 2.2, 1, 1), cex = 0.8)
plot.ts(exportaciones_ts, col = "darkgrey", plot.type = "single", 
        xlab = "", ylab = "", lty = 1.5)
lines(result, col = "red", lwd = 0.8)
abline(v = c(2018,10), col = "blue", lwd = 1, lty = 3)
legend("topleft", legend = c("Observed Deflator", "Stochastic level"), 
       lwd = c(1.5, 1), col = c("red", "darkgrey"), bty = "n")



