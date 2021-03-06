library(tidyverse)
library(here)
library(glue)
library(bama)
library(mvnfast)
library(SuperLearner)
library(crossurr)
library(HIMA)
library(clustermq)
library(freebird)

simfn <- function(n, p, q, sig = 1, R = 0.8, linear = TRUE, run = 0, write = TRUE) {
  library(tidyverse)
  library(here)
  library(glue)
  library(bama)
  library(mvnfast)
  library(SuperLearner)
  library(crossurr)
  library(HIMA)
  library(clustermq)
  library(freebird)
  
  # browser()
  set.seed(0)
  ax_beta <- rnorm(q)
  beta_s <- matrix(c(seq(-1, 1, length = 5),
                     rep(0, q - 5)),
                     q, p)
  beta_s0 <- matrix(c(seq(-2, 0, length = 5), 
                      rep(0, q - 5)),
                    q, p)
  alpha_s1 <- matrix(c(0.75, 0.25, runif(p-2)), n, p, byrow = TRUE)
  alpha_s0 <- matrix(c(0, 0, runif(p- 2) - 0.5), n, p, byrow = TRUE)
  if (!linear) {
    delta <- exp(alpha_s1[1,1]) + alpha_s1[1,2] - exp(alpha_s0[1,1]) - alpha_s0[1,2]
  } else {
    delta <- alpha_s1[1,1] + alpha_s1[1,2] - alpha_s0[1,1] - alpha_s0[1,2]
  }
  
  Delta_s <- delta*(1-R)/R
  Delta <- Delta_s + delta
  
  set.seed(run)
  x <- rnorm(n*q, sd = 1) %>% matrix(n, q)
  a <- rbinom(n, prob = plogis(x %*% ax_beta), size = 1)
  
  s_1 <- alpha_s1 + x %*% beta_s + rnorm(n*p, sd = sig) %>% matrix(n, p)
  s_0 <- alpha_s0 + x %*% beta_s0 + rnorm(n*p, sd = sig) %>% matrix(n, p)
  s <- s_1*a + (1-a)*s_0
  
  if (linear) {
    y_1 <- Delta_s + rowMeans(x[,1:pmin(q, 25)]) + s_1[,1] + s_1[,2] + rnorm(n, sd = sig)
    y_0 <- rowMeans(x[,1:pmin(q, 25)]) + s_0[,1] + s_0[,2] + rnorm(n, sd = sig)
  } else {
    y_1 <- Delta_s + rowMeans(x) + exp(s_1[,1]) + s_1[,2] + rnorm(n, sd = sig)
    y_0 <- rowMeans(x) + exp(s_0[,1]) + s_0[,2] + rnorm(n, sd = sig)
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
  
  xf_surr <- xf_surrogate(ds = wds,
                          x = paste('x.', 1:q, sep =''),
                          s = paste('s.', 1:p, sep =''),
                          a = 'a',
                          y = 'y',
                          K = 4,
                          outcome_learners = c("SL.mean", "SL.glmnet", 
                                               "SL.ridge", "SL.lm", "SL.svm", 'SL.ranger'),
                          ps_learners = c("SL.mean", "SL.glmnet", 
                                          "SL.glm", "SL.lda", "SL.qda",
                                          "SL.svm", 'SL.ranger'),
                          trim_at = 0.01,
                          mthd = 'superlearner')
  xfl_surr <- xf_surrogate(ds = wds,
                          x = paste('x.', 1:q, sep =''),
                          s = paste('s.', 1:p, sep =''),
                          a = 'a',
                          y = 'y',
                          K = 4,
                          trim_at = 0.01,
                          mthd = 'lasso')
  
  dr_delta <- xf_surr$deltahat
  drl_delta <- xfl_surr$deltahat
  dr_delta_s <- xf_surr$deltahat_s
  drl_delta_s <- xfl_surr$deltahat_s
  dr_r <- 1 - dr_delta_s/dr_delta
  drl_r <- 1 - drl_delta_s/drl_delta
  dr_r_cil <- dr_r - 1.96*xf_surr$R_se
  dr_r_cih <- dr_r + 1.96*xf_surr$R_se
  drl_r_cil <- drl_r - 1.96*xfl_surr$R_se
  drl_r_cih <- drl_r + 1.96*xfl_surr$R_se
  
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
  bama_r <- 1 - mean(bama_tst$beta.a/(bama_tst$beta.a + rowSums(bama_tst$alpha.a* bama_tst$beta.m)))
  
  bama_ds_cil <- quantile(bama_tst$beta.a, 0.025)
  bama_ds_cih <- quantile(bama_tst$beta.a, 0.975)
  bama_d_cil <- quantile(bama_tst$beta.a + rowSums(bama_tst$alpha.a* bama_tst$beta.m), 0.025)
  bama_d_cih <- quantile(bama_tst$beta.a + rowSums(bama_tst$alpha.a* bama_tst$beta.m), 0.975)
  bama_r_cil <- quantile(bama_tst$beta.a/(bama_tst$beta.a + rowSums(bama_tst$alpha.a* bama_tst$beta.m)), 0.025)
  bama_r_cih <- quantile(bama_tst$beta.a/(bama_tst$beta.a + rowSums(bama_tst$alpha.a* bama_tst$beta.m)), 0.975)
  # browser()
  
  hima_tst <- try(hima(X = a,
                   Y = y,
                   M = s,
                   COV.XM = x,
                   COV.MY = x,
                  verbose = TRUE), silent = TRUE)
  if (inherits(hima_tst, 'try-error')) {
    hima_nie <- NA
    hima_r <- NA
    hima_nde <- NA
    hima_delta_s <- NA
    hima_delta <- NA
  } else {
    hima_nie <- sum(hima_tst$`alpha*beta`)
    hima_r <- sum(hima_tst$`% total effect`/100)
    hima_nde <- (1-hima_r)*hima_nie/hima_r
    hima_delta_s <- hima_nde
    hima_delta <- hima_nde + hima_nie
  }
  
  fb_tst <- freebird::hilma(Y = y,
                            G = s,
                            S = a %>% matrix(n, 1))
  fb_nie <- fb_tst$beta_hat
  fb_nde <- fb_tst$alpha1_hat
  fb_delta_s <- fb_nde
  fb_delta <- fb_nde + fb_nie
  fb_r <- fb_nie/fb_delta
  
  Delta_0 <- mean(y_1) - mean(y_0)
  R_0 <- 1 - Delta_s/Delta_0
  
  # browser()
  
  out_ds <- data.frame(
    type = c('true', 'xfdr', 'xfl', 'bama', 'hima', 'fb'),
    # Delta = c(Delta, dr_delta, drl_delta, bama_delta, hima_delta, fb_delta),
    Delta = c(Delta_0, dr_delta, drl_delta, bama_delta, hima_delta, fb_delta),
    Delta_s = c(Delta_s, dr_delta_s, drl_delta_s, bama_delta_s, hima_delta_s, fb_delta_s),
    # R = c(R, dr_r, drl_r, bama_r, hima_r)
    R = c(R_0, dr_r, drl_r, bama_r, hima_r, fb_r),
    R_cil = c(R, dr_r_cil, drl_r_cil, bama_r_cil, NA, NA),
    R_cih = c(R, dr_r_cih, drl_r_cih, bama_r_cih, NA, NA)
  ) %>% as_tibble
  if (write) {
    write_csv(out_ds, glue('/n/scratch3/users/d/dma12/doubly-robust-surrogate/res_n{n}-p{p}-q{q}-l{linear}-{run}.csv'))
  }
  out_ds
}
sim_params <- expand.grid(n = c(100, 500),
                          p = 100,
                          q = 100,
                          linear = TRUE,
                          sig = c(0.1, 0.5),
                          R = 0.5,
                          run = 1:1000)
# tst <- sim_params %>% filter(p < 5000, q < 5000, n < 1000) %>% sample_n(1)
# tst
# with(tst, simfn(n = n,
#                 p = p,
#                 q = q,
#                 sig = 0.0001,
#                 R = 0.5,
#                 run = run-1,
#                 write = FALSE))

options(
  clustermq.defaults = list(ptn="short",
                            log_file="Rout/log%a.log",
                            time_amt = "12:00:00"
  )
)
sim_res <- Q_rows(sim_params, simfn, 
                  fail_on_error = FALSE,
                  n_jobs = 200)
saveRDS(sim_res, here('results/05_lnl-sim-results.rds'))
