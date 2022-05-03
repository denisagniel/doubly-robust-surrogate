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

simfn <- function(n, dim_s, zeta_3, run = 0, write = TRUE) {
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
  set.seed(run)

  x <- rnorm(n)
  xx <- matrix(c(x, rnorm(n)), n, 2)
  a <- rbinom(n, prob = plogis(x), size = 1)
  
  beta_1 <- 2
  beta_2 <- 0.1
  eps_u <- rnorm(n)
  u <- beta_1*a + beta_2*x*a + eps_u
  u1 <- beta_1 + beta_2*x + eps_u
  u0 <- eps_u
  
  gamma_1 <- 0.4
  gamma_2 <- 2
  eps_s <- rnorm(n)
  s <- gamma_1*a + gamma_2*u + eps_s
  s1 <- gamma_1 + gamma_2*u1 + eps_s
  s0 <- gamma_2*u0 + eps_s
  
  ss <- matrix(c(s, rnorm(n*(dim_s - 1))), n, dim_s)
  
  zeta_1 <- 0.1
  zeta_2 <- 1
  eps_z <- rnorm(n)
  y <- pnorm(x^2) + x^2*a + zeta_1*a + zeta_2*u + zeta_3*u*a + eps_z
  y1 <- pnorm(x^2) + x^2 + zeta_1 + zeta_2*u1 + zeta_3*u1 + eps_z
  y0 <- pnorm(x^2) +  zeta_2*u0 + eps_z
  delta_0 <- mean(y1 - y0)
  
  e_u1 <- mean(beta_1 + beta_2*x + gamma_2/(gamma_2^2 + 1)*(s - gamma_1 - gamma_2*beta_1 - gamma_2*beta_2*x))
  e_u0 <- mean(gamma_2/(gamma_2^2 + 1)*s)
  
  delta_s0 <- zeta_1 + (zeta_2 + zeta_3)*e_u1 - zeta_2*e_u0 + mean(x^2)
  r0 <- 1 - delta_s0/delta_0
  
  dsi <- tibble(
    id = 1:n,
    a, y
  ) 
  sds <- tibble(
    id = rep(1:n, dim_s),
    s = c(ss),
    sn = glue('s.{rep(1:dim_s, each = n)}')
  )
  xds <- tibble(
    id = rep(1:n, 2),
    x = c(xx),
    xn = glue('x.{rep(1:2, each = n)}')
  ) %>%
    spread(xn, x)
  
  ds <- dsi %>%
    inner_join(sds)
  wds <- ds %>%
    spread(sn, s) %>%
    inner_join(xds)
  
  xf_surr <- xf_surrogate(ds = wds,
                          x = paste('x.', 1:2, sep =''),
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
  # browser()
  xfl_surr <- xf_surrogate(ds = wds,
                           x = paste('x.', 1:2, sep =''),
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
  
  # browser()
  
  out_ds <- data.frame(
    type = c('true', 'xfdr',  'xfl'),
    # Delta = c(Delta, dr_delta, drl_delta, bama_delta, hima_delta, fb_delta),
    Delta = c(delta_0, dr_delta,  drl_delta),
    Delta_s = c(delta_s0, dr_delta_s,   drl_delta_s),
    # R = c(R, dr_r, drl_r, bama_r, hima_r)
    R = c(r0, dr_r, drl_r),
    R_cil = c(r0, dr_r_cil, drl_r_cil),
    R_cih = c(r0, dr_r_cih, drl_r_cih)
  ) %>% as_tibble
  if (write) {
    write_csv(out_ds, glue('/n/scratch3/users/d/dma12/doubly-robust-surrogate/res-surr_dim{dim_s}-{run}.csv'))
  }
  out_ds
}

sim_params <- expand.grid(n = c(50, 1000),
                          dim_s = c(1, 10, 50),
                          zeta_3 = c(0, 10),
                          run = 1:3)
# tst <- sim_params %>% sample_n(1)
# tst
# with(tst, simfn(n = n,
#                 dim_s = dim_s,
#                 zeta_3 = zeta_3,
#                 run = run,
#                 write = FALSE))


options(
  clustermq.defaults = list(ptn="short",
                            log_file="Rout/log%a.log",
                            time_amt = "1:00:00"
  )
)
sim_res <- Q_rows(sim_params, simfn, 
                  fail_on_error = FALSE,
                  n_jobs = 4)
saveRDS(sim_res, here('results/31_surrogate-only-results.rds'))
