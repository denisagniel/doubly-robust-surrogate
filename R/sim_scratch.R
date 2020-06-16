library(tidyverse)
library(here)
library(glue)
library(bama)
library(mvnfast)

n <- 1000
p <- 500
rho <- 0.5
x <- rnorm(n)
a <- rbinom(n, prob = plogis(x), size = 1)
s_1 <- 0.5 + 1.9*x^2 + rnorm(n)
s_0 <- 1.7 + rnorm(n)
s <- s_1*a + (1-a)*s_0

lm(s ~ a + x)
lm(s ~ a*x)
mean(s_1) - mean(s_0)

Delta_s <- 1
y_1 <- Delta_s + x + s_1 + rnorm(n)
y_0 <- -x + s_0 + rnorm(n)
y <- y_1*a + (1-a)*y_0

ds <- tibble(
  x, a, s, y
)
smod <- lm(s ~ a*x, data = ds)
ymod <- lm(y ~ a*(x + s), data = ds)
psmod <- glm(a ~ s*x, data = ds)

ymod0 <- lm(y ~ a*(x), data = ds)
psmod0 <- glm(a ~ x, data = ds)
Delta <- mean(y_1) - mean(y_0)
R <- 1 - Delta_s/Delta
ds <- ds %>%
  mutate(pi = predict(psmod, type= 'response'),
         yhat1 = predict(ymod, newdata = ds %>%
                          mutate(a = 1)),
         yhat0 = predict(ymod, newdata = ds %>%
                          mutate(a = 0)),
         u_s = yhat1 - yhat0 + 
           a*(y - yhat1)/pi -
           (1-a)*(y-yhat0)/(1-pi),
         
         pi_0 = predict(psmod0, type= 'response'),
         yhat1_0 = predict(ymod0, newdata = ds %>%
                           mutate(a = 1)),
         yhat0_0 = predict(ymod0, newdata = ds %>%
                           mutate(a = 0)),
         u_0 = yhat1_0 - yhat0_0 + 
           a*(y - yhat1_0)/pi_0 -
           (1-a)*(y-yhat0_0)/(1-pi_0))
dr_delta_s <- mean(ds$u_s)
dr_delta <- mean(ds$u_0)
dr_r <- 1 - dr_delta_s/dr_delta
bama_tst <- bama(Y = as.vector(y),
                 A = a,
                 M = s %>% matrix(n, 1),
                 C1 = cbind(1, x),
                 C2 = cbind(1, x),
                 beta.m = rep(0, p),
                 alpha.a = rep(0, p),
                 burnin = 1000,
                 ndraws = 500)
bama_nie <- colMeans(bama_tst$alpha.a* bama_tst$beta.m) %>% sum
bama_nde <- mean(bama_tst$beta.a)
bama_te <- bama_nie + bama_nde

bama_delta_s <- bama_nde
bama_delta <- bama_te
bama_r <- bama_nie/bama_te

