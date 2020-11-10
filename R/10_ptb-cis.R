library(tidyverse)
library(here)
library(glue)
library(mvnfast)
library(SuperLearner)
library(crossurr)
library(clustermq)


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
  # browser()
  xf_surr <- map(1:20, function(i) {
    xf_surrogate(ds = wds,
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
                          mthd = 'superlearner',
                          n_ptb = 1000,
                 seed = i)
  }) %>%
    bind_rows %>%
    summarise_all(median)
  xfl_surr <- map(1:20, function(i) {
    xf_surrogate(ds = wds,
                           x = paste('x.', 1:q, sep =''),
                           s = paste('s.', 1:p, sep =''),
                           a = 'a',
                           y = 'y',
                           K = 4,
                           trim_at = 0.01,
                           mthd = 'lasso',
                           n_ptb = 1000,
                 seed = i)
  }) %>%
    bind_rows %>%
    summarise_all(median)
  
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
  
  # browser()
  
  out_ds <- data.frame(
    type = c('xfdr', 'xfl'),
    R = c(dr_r, drl_r),
    R_cil = c(dr_r_cil,drl_r_cil),
    R_cih = c(dr_r_cih,drl_r_cih),
    R_bci_l = c(xf_surr$R_qci_l, xfl_surr$R_qci_l),
    R_bci_h = c(xf_surr$R_qci_h, xfl_surr$R_qci_h),
    R_se = c(xf_surr$R_se, xfl_surr$R_se),
    R_ptb_se = c(xf_surr$R_ptb_se, xfl_surr$R_ptb_se)
  ) %>% as_tibble
  if (write) {
    write_csv(out_ds, glue('/n/scratch3/users/d/dma12/doubly-robust-surrogate/ptb_n{n}-p{p}-R{R}-{run}.csv'))
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
tst <- sim_params %>% filter(n < 1000) %>% sample_n(1)
tst
with(tst, simfn(n = n,
                p = p,
                q = q,
                R = R,
                run = runif(1),
                write = FALSE))

# simfn(n = 500, p = 50, R = 0.9, write = FALSE, run = -12)


options(
  clustermq.defaults = list(ptn="medium",
                            log_file="Rout/log%a.log",
                            time_amt = "72:00:00"
  )
)
sim_res <- Q_rows(sim_params, simfn, 
                  fail_on_error = FALSE,
                  n_jobs = 500)
saveRDS(sim_res, here('results/10_ptb-cis.rds'))
