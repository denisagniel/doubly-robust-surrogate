library(tidyverse)
library(here)
library(glue)
library(bama)
library(mvnfast)
library(SuperLearner)
library(crossurr)
library(HIMA)
library(clustermq)

simfn <- function(n, p, q, R = 0.8, linear = TRUE, run = 0) {
  library(tidyverse)
  library(here)
  library(glue)
  library(bama)
  library(mvnfast)
  library(SuperLearner)
  library(crossurr)
  library(HIMA)
  library(clustermq)
  
  set.seed(0)
  ax_beta <- rnorm(q)
  beta_s <- matrix(seq(-1, 1, length = q), q, p, byrow = TRUE)
  beta_s0 <- matrix(seq(-2, 0, length = q), q, p, byrow = TRUE)
  alpha_s1 <- matrix(c(0.75, 0.25, runif(p-2)), n, p, byrow = TRUE)
  alpha_s0 <- matrix(c(0, 0, runif(p- 2) - 0.5), n, p, byrow = TRUE)
  delta <- alpha_s1[1,1] + alpha_s1[1,2] - alpha_s0[1,1] - alpha_s0[1,2]
  Delta_s <- delta*(1-R)/R
  Delta <- Delta_s + delta
  
  set.seed(run)
  x <- rnorm(n*q) %>% matrix(n, q)
  a <- rbinom(n, prob = plogis(x %*% ax_beta), size = 1)
  
  s_1 <- alpha_s1 + x %*% beta_s + rnorm(n*p) %>% matrix(n, p)
  s_0 <- alpha_s0 + x %*% beta_s0 + rnorm(n*p) %>% matrix(n, p)
  s <- s_1*a + (1-a)*s_0
  
  if (linear) {
    y_1 <- Delta_s + rowMeans(x) + s_1[,1] + s_1[,2] + rnorm(n)
    y_0 <- rowMeans(x) + s_0[,1] + s_0[,2] + rnorm(n)
  } else {
    y_1 <- Delta_s + rowMeans(exp(x)) + s_1[,1] + s_1[,2] + rnorm(n)
    y_0 <- rowMeans(exp(x)) + s_0[,1] + s_0[,2] + rnorm(n)
  }
  y <- y_1*a + (1-a)*y_0
  
  
  dsi <- tibble(
    id = 1:n,
    a, y
  ) 
  sds <- tibble(
    id = rep(1:n, p),
    s = c(s),
    sn = glue('s.{rep(1:p, each = n)}')
  )
  xds <- tibble(
    id = rep(1:n, q),
    x = c(x),
    xn = glue('x.{rep(1:q, each = n)}')
  ) %>%
    spread(xn, x)
  
  ds <- dsi %>%
    inner_join(sds)
  wds <- ds %>%
    spread(sn, s) %>%
    inner_join(xds)
  xf_delta_s <- xfit_dr(ds = wds,
                        xvars = c(paste('s.', 1:p, sep =''),
                                  paste('x.', 1:q, sep ='')),
                        yname = y,
                        aname = a,
                        K = 4,
                        outcome_learners = c("SL.mean", "SL.glmnet", 
                                             "SL.ridge", "SL.lm", "SL.svm", 'SL.ranger'),
                        ps_learners = c("SL.mean", "SL.glmnet", 
                                        "SL.glm", "SL.lda", "SL.qda",
                                        "SL.svm", 'SL.ranger'),
                        trim_at = 0.01,
                        mthd = 'superlearner')
  
  xf_delta <- xfit_dr(ds = wds,
                      xvars = c(paste('x.', 1:q, sep ='')),
                      yname = y,
                      aname = a,
                      K = 4,
                      outcome_learners = c("SL.mean", "SL.glmnet", 
                                           "SL.ridge", "SL.lm", "SL.svm", 'SL.ranger'),
                      ps_learners = c("SL.mean", "SL.glmnet", 
                                      "SL.glm", "SL.lda", "SL.qda",
                                      "SL.svm", 'SL.ranger'),
                      trim_at = 0.01,
                      method = 'superlearner')
  
  dr_delta_s <- xf_delta_s$estimate
  dr_delta <- xf_delta$estimate
  dr_r <- 1 - dr_delta_s/dr_delta
  
  bama_tst <- bama(Y = as.vector(y),
                   A = a,
                   M = s,
                   C1 = cbind(1, x),
                   C2 = cbind(1, x),
                   beta.m = rep(0, p),
                   alpha.a = rep(0, p),
                   burnin = 2000,
                   ndraws = 1000)
  bama_nie <- colMeans(bama_tst$alpha.a* bama_tst$beta.m) %>% sum
  bama_nde <- mean(bama_tst$beta.a)
  bama_te <- bama_nie + bama_nde
  
  bama_delta_s <- bama_nde
  bama_delta <- bama_te
  bama_r <- bama_nie/bama_te
  
  
  hima_tst <- hima(X = a,
                   Y = y,
                   M = s,
                   COV.XM = x,
                   COV.MY = x,
                   verbose = TRUE)
  hima_nie <- sum(hima_tst$`alpha*beta`)
  hima_r <- sum(hima_tst$`% total effect`/100)
  hima_nde <- (1-hima_r)*hima_nie/hima_r
  hima_delta_s <- hima_nde
  hima_delta <- hima_nde + hima_nie
  
  
  out_ds <- tibble(
    type = c('true', 'xfdr', 'bama', 'him'),
    Delta = c(Delta, dr_delta, bama_delta, hima_delta),
    Delta_s = c(Delta_s, dr_delta_s, bama_delta_s, hima_delta_s),
    R = c(R, dr_r, bama_r, hima_r)
  )
  out_ds
}

sim_params <- expand.grid(n = c(50, 500, 1000),
                          p = c(5, 100, 250),
                          q = c(5, 100, 250),
                          linear = c(TRUE, FALSE),
                          R = 0.5,
                          run = 1:1000)
# tst <- sim_params %>% filter(p < 5000, q < 5000, n < 1000) %>% sample_n(2)
# tst
# Q_rows(tst, simfn)

options(
  clustermq.defaults = list(ptn="short",
                            log_file="Rout/log%a.log",
                            time_amt = "12:00:00"
  )
)
sim_res <- Q_rows(sim_params, simfn, 
                  fail_on_error = FALSE,
                  n_jobs = 750)
saveRDS(sim_res, here('results/05_lnl-sim-results.rds'))
