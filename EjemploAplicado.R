library(astsa)
library(tidyverse)
library(GGally)
library(gridExtra)
library(gifski)
library(forecast)
library(tseries)
library(xts)

##### plot de acf #####
ggacf <- function(x, type = c("correlation", "covariance", "partial"), ci = 0.95, ...) {
  type <- match.arg(type)
  x.acf <- acf(x, type = type, plot = FALSE, ...)
  x.df <- with(x.acf, data.frame(lag, acf))
  x.ci <- qnorm((1 + ci)/2) / sqrt(length(x))
  print(x.ci)
  ggplot(data = x.df, aes(x = lag, y = acf)) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(aes(xend = lag, yend = 0)) +
    geom_hline(yintercept = c(x.ci, -x.ci), color = "blue", linetype = "dashed") +
    scale_y_continuous(breaks = seq(0.0, 1.0, 0.2))
}

###### Leer Bitcoin ######

BTC_Daily <- read.csv("BTC-Daily.csv")
Data <- data.frame(BTC_Daily$date,BTC_Daily$close)
colnames(Data) <- c("FechaTiempo", "Valor")

Data$FechaTiempo <- strftime(Data$FechaTiempo, format="%Y-%m-%d")
Data$FechaTiempo <- as.Date(Data$FechaTiempo)

Bitcoin <- Data %>%
  filter(FechaTiempo >= as.Date("2017-01-01"),
         FechaTiempo <= as.Date("2021-12-31")) |> 
  arrange(FechaTiempo)

head(Bitcoin); tail(Bitcoin)


## Diferencias sin tranformación logarítmica

returnPlot1 <- Bitcoin %>%
  select(ends_with("Valor")) %>%
  apply(.,2,function(x) diff(x)) %>% 
  as_tibble %>% 
  mutate(date = Bitcoin$FechaTiempo[-1]) %>% 
  gather(coin, usd,-date) %>% 
  ggplot(aes(date,usd)) +
  geom_point() +
  ylab("Price Difference (USD)") +
  facet_wrap(~coin,nrow=3,scales="free_y")

returnPlot1

## Diferencias con tranformación logarítmica

returnPlot2 <- Bitcoin %>%
  select(ends_with("Valor")) %>%
  apply(.,2,function(x) log10(x) %>% diff) %>% 
  as_tibble %>% 
  mutate(date = Bitcoin$FechaTiempo[-1]) %>% 
  gather(coin, log_return,-date) %>% 
  ggplot(aes(date,log_return)) +
  geom_point() +
  ylab("Daily Log Return (log10)") +
  facet_wrap(~coin,nrow=3,scales="free_y")

returnPlot2

## ACF y PACF

Bitcoin %>%
  pull("Valor") %>%
  log10 %>%
  ggacf +
  ggtitle("Valor")


Bitcoin %>%
  pull("Valor") %>%
  log10 %>%
  ggacf(type="partial") +
  ylab("Partial ACF") +
  ggtitle("Valor")


## Normality

qqplot <- Bitcoin %>%
  select(ends_with("Valor")) %>%
  apply(.,2,function(x) log10(x) %>% diff) %>% 
  as_tibble %>% 
  mutate(date = Bitcoin$FechaTiempo[-1]) %>% 
  gather(coin, log_return,-date) %>% 
  ggplot(aes(sample = log_return)) +
  geom_qq() + 
  geom_qq_line() + 
  facet_wrap(~coin,nrow=3,scales="free_y")

qqplot
