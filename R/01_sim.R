library(tidyverse)
library(here)
library(glue)
library(bama)
library(mvnfast)
library(SuperLearner)
library(crossurr)
library(HIMA)
library(clustermq)

simfn <- function(n, p, R = 0.5, rho = 0.4, run = 0) {
  Delta <- 2.25
  Delta_s <- Delta*(1-R)
  
  Sigma <- matrix(rho, p, p) + (1-rho)*diag(p)
  x1 <- runif(n, -2, 5)
  x2 <- rbinom(n, prob = 0.5, size = 1)
  x <- cbind(x1, x2)
  a <- rbinom(n, prob = plogis(-x1 + 2*x1*x2), size = 1)
  
  s_spm <- matrix(rep(1:0, c(10, p - 10)), n, p, byrow = TRUE)
  s_1 <-  1.5 + (x1 + x2)*s_spm + rmvn(n, mu = rep(0, p), sigma = Sigma)
  s_0 <- 2 + x2*s_spm - x1*x2 + rmvn(n, mu = rep(0, p), sigma = Sigma)
  s <- s_1*a + (1-a)*s_0
  
  y_1 <- Delta_s + x[,1] + x[,2] + rowMeans(s_1[,1:15]) + rnorm(n)
  y_0 <- x[,1] + x[,2] + rowMeans(s_0[,1:15]) + rnorm(n)
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
    id = rep(1:n, 2),
    x = c(x),
    xn = glue('x.{rep(1:2, each = n)}')
  ) %>%
    spread(xn, x)
  
  ds <- dsi %>%
    inner_join(sds)
  wds <- ds %>%
    spread(sn, s) %>%
    inner_join(xds)
  
  xf_delta_s <- xfit_dr(ds = wds,
                        xvars = c(paste('s.', 1:p, sep =''),
                                  'x.1', 'x.2'),
                        yname = y,
                        aname = a,
                        K = 4,
                        outcome_learners = c("SL.mean", "SL.glmnet", 
                                             "SL.ridge", "SL.lm", "SL.svm", 'SL.ranger'),
                        ps_learners = c("SL.mean", "SL.glmnet", 
                                        "SL.glm", "SL.lda", "SL.qda",
                                        "SL.svm", 'SL.ranger'),
                        cvControl = list(V = 4))
  xf_delta <- xfit_dr(ds = wds,
                      xvars = c('x.1', 'x.2'),
                      yname = y,
                      aname = a,
                      K = 4,
                      outcome_learners = c("SL.mean", "SL.glmnet", 
                                           "SL.ridge", "SL.lm", "SL.svm", 'SL.ranger', 'SL.glm.interaction', 'SL.step.interaction', 'SL.rpart'),
                      ps_learners = c("SL.mean", "SL.glmnet", 
                                      "SL.glm", "SL.lda", "SL.qda",
                                      "SL.svm", "SL.glm.interaction",'SL.bayesglm', 'SL.rpart', 'SL.knn'),
                      cvControl = list(V = 4))
  
  
  
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
                   burnin = 3000,
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

sim_params <- expand.grid(n = c(100, 500, 1000),
                          p = c(50, 100, 500),
                          R = c(0.2, 0.9),
                          run = 1:3)

options(
  clustermq.defaults = list(ptn="short",
                            log_file="Rout/log%a.log",
                            time_amt = "12:00:00"
  )
)
sim_res <- Q_rows(sim_params, simfn, 
                  fail_on_error = FALSE,
                  n_jobs = 20)
saveRDS(sim_res, here('results/tst_sim-results.rds'))
