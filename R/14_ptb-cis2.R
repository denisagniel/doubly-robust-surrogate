library(tidyverse)
library(here)
library(glue)
library(mvnfast)
library(SuperLearner)
library(crossurr)
library(clustermq)


simfn <- function(n, p, R = 0.5, rho = 0.4, run = 0, write = TRUE) {
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
  
  Delta <- 2.25
  Delta_s <- Delta*(1-R)
  q <- 2
  
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

sim_params <- expand.grid(n = 500,
                          p = c(50, 100, 500),
                          R = c(0.2, 0.9),
                          run = 1:1000)
tst <- sim_params %>% filter(n < 1000) %>% sample_n(1)
tst
with(tst, simfn(n = n,
                p = p,
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
saveRDS(sim_res, here('results/14_ptb-cis2.rds'))
