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

simfn <- function(n, p, dim_s, beta_sparse, R = 0.5, rho = 0.4, run = 0, write = TRUE) {
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
  
  set.seed(run)
  Delta <- 2.25
  Delta_s <- Delta*(1-R)
  q <- 2
  sig <- 1/sqrt(p)
  sig_s<- 3
  
  Sigma <- matrix(rho, p, p) + (1-rho)*diag(p)
  x1 <- runif(n, -2, 5)
  x2 <- rbinom(n, prob = 0.5, size = 1)
  x <- cbind(x1, x2)
  a <- rbinom(n, prob = plogis(-x1 + 2*x1*x2), size = 1)
  
  s_spm <- matrix(rep(1:0, c(10, p - 10)), n, p, byrow = TRUE)
  m_1 <-  1.5 + (x1 + x2)*s_spm + rmvn(n, mu = rep(0, p), sigma = Sigma)*sig_s
  m_0 <- 2 + x2*s_spm - x1*x2 + rmvn(n, mu = rep(0, p), sigma = Sigma)*sig_s
  # eta_1 <- 1.5 + 10*(x1 + x2)*s_spm + rmvn(n, mu = rep(0, p), sigma = Sigma)
  # eta_0 <- 2 + 10*x2*s_spm - 10*x1*x2 + rmvn(n, mu = rep(0, p), sigma = Sigma)
  # s_1 <- rbinom(n*p, size = 1, prob = plogis(eta_1)) %>% matrix(n, p)
  # s_0 <- rbinom(n*p, size = 1, prob = plogis(eta_0)) %>% matrix(n, p)
  m <- m_1*a + (1-a)*m_0
  
  beta_1 <- abs(rnorm(p*dim_s))/p %>% matrix(p, dim_s)
  beta_1[beta_sparse:p,] <- 0
  
  s_1 <- m_1 %*% beta_1 + rnorm(n*dim_s, sd = 0.1) %>% matrix(n, dim_s)
  s_0 <- m_0 %*% beta_1 + rnorm(n*dim_s, sd = 0.1) %>% matrix(n, dim_s)
  s <- s_1*a + (1-a)*s_0
  
  # browser()
  
  # browser()
  y_1 <- Delta_s + x[,1] + x[,2] + rowMeans(m_1[,1:15]) + rnorm(n, sd = sig)
  y_0 <- x[,1] + x[,2] + rowMeans(m_0[,1:15]) + rnorm(n, sd = sig)
  y <- y_1*a + (1-a)*y_0
  
  dsi <- tibble(
    id = 1:n,
    a, y
  ) 
  sds <- tibble(
    id = rep(1:n, dim_s),
    s = c(s),
    sn = glue('s.{rep(1:dim_s, each = n)}')
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
  
  xf_surr <- xf_surrogate(ds = wds,
                          x = paste('x.', 1:q, sep =''),
                          s = paste('s.', 1:dim_s, sep =''),
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
  # xfr_surr <- xfr_surrogate(ds = wds,
  #                           x = paste('x.', 1:q, sep =''),
  #                           s = paste('s.', 1:p, sep =''),
  #                           a = 'a',
  #                           y = 'y',
  #                           K = 4,
  #                           outcome_learners = c("SL.mean", "SL.glmnet", 
  #                                                "SL.ridge", "SL.lm", "SL.svm", 'SL.ranger'),
  #                           ps_learners = c("SL.mean", "SL.glmnet", 
  #                                           "SL.glm", "SL.lda", "SL.qda",
  #                                           "SL.svm", 'SL.ranger'),
  #                           trim_at = 0.01,
  #                           mthd = 'superlearner',
  #                           splits = 5)
  xfl_surr <- xf_surrogate(ds = wds,
                           x = paste('x.', 1:q, sep =''),
                           s = paste('s.', 1:dim_s, sep =''),
                           a = 'a',
                           y = 'y',
                           K = 4,
                           trim_at = 0.01,
                           mthd = 'lasso')
  dr_delta <- xf_surr$deltahat
  drl_delta <- xfl_surr$deltahat
  # drr_delta <- xfr_surr$Dm
  dr_delta_s <- xf_surr$deltahat_s
  drl_delta_s <- xfl_surr$deltahat_s
  # drr_delta_s <- xfr_surr$Dsm
  dr_r <- 1 - dr_delta_s/dr_delta
  drl_r <- 1 - drl_delta_s/drl_delta
  # drr_r <- xfr_surr$Rm
  dr_r_cil <- dr_r - 1.96*xf_surr$R_se
  dr_r_cih <- dr_r + 1.96*xf_surr$R_se
  drl_r_cil <- drl_r - 1.96*xfl_surr$R_se
  drl_r_cih <- drl_r + 1.96*xfl_surr$R_se
  # drr_r_cil <- xfr_surr$R_cil0
  # drr_r_cih <- xfr_surr$R_cih0
  
  
  
  Delta_0 <- mean(y_1) - mean(y_0)
  R_0 <- 1 - Delta_s/Delta_0
  
  # browser()
  
  out_ds <- data.frame(
    type = c('true', 'xfdr',  'xfl'),
    # Delta = c(Delta, dr_delta, drl_delta, bama_delta, hima_delta, fb_delta),
    Delta = c(Delta_0, dr_delta,  drl_delta),
    Delta_s = c(Delta_s, dr_delta_s,   drl_delta_s),
    # R = c(R, dr_r, drl_r, bama_r, hima_r)
    R = c(R_0, dr_r, drl_r),
    R_cil = c(R, dr_r_cil, drl_r_cil),
    R_cih = c(R, dr_r_cih, drl_r_cih)
  ) %>% as_tibble
  if (write) {
    write_csv(out_ds, glue('/n/scratch3/users/d/dma12/doubly-robust-surrogate/res-surr_dim{dim_s}-{run}.csv'))
  }
  out_ds
}

sim_params <- expand.grid(n = 500,
                          p = 50,
                          R = 0.5,
                          beta_sparse = c(5, 15, 25, 50),
                          dim_s = c(1, 10, 50, 75),
                          run = 1:1000)
tst <- sim_params %>% filter(n < 1000) %>% sample_n(1)
tst
with(tst, simfn(n = n,
                p = p,
                R = R,
                beta_sparse = beta_sparse,
                dim_s = dim_s,
                run = run,
                write = TRUE))
# # 
# simfn(n = 500, p = 50, R = 0.9, write = FALSE, run = 850)


options(
  clustermq.defaults = list(ptn="short",
                            log_file="Rout/log%a.log",
                            time_amt = "12:00:00"
  )
)
sim_res <- Q_rows(sim_params, simfn, 
                  fail_on_error = FALSE,
                  n_jobs = 200)
saveRDS(sim_res, here('results/31_surrogate-only-results.rds'))
