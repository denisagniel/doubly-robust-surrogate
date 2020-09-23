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
  library(mvnfast)
  library(SuperLearner)
  library(crossurr)
  library(clustermq)
  
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
  
  xf_delta_s <- xfit_dr(ds = wds,
                        x = c(paste('s.', 1:p, sep =''),
                              paste('x.', 1:q, sep ='')),
                        y = 'y',
                        a = 'a',
                        K = 4,
                        outcome_learners = c("SL.mean", "SL.glmnet", 
                                             "SL.ridge", "SL.lm", "SL.svm", 'SL.ranger'),
                        ps_learners = c("SL.mean", "SL.glmnet", 
                                        "SL.glm", "SL.lda", "SL.qda",
                                        "SL.svm", 'SL.ranger'),
                        trim_at = 0.01,
                        mthd = 'superlearner')
  
  xf_delta <- xfit_dr(ds = wds,
                      x = c(paste('x.', 1:q, sep ='')),
                      y = 'y',
                      a = 'a',
                      K = 4,
                      outcome_learners = c("SL.mean", "SL.glmnet", 
                                           "SL.ridge", "SL.lm", "SL.svm", 'SL.ranger'),
                      ps_learners = c("SL.mean", "SL.glmnet", 
                                      "SL.glm", "SL.lda", "SL.qda",
                                      "SL.svm", 'SL.ranger'),
                      trim_at = 0.01,
                      mthd = 'superlearner')
  
  dr_delta_s <- xf_delta_s$estimate
  dr_delta <- xf_delta$estimate
  dr_r <- 1 - dr_delta_s/dr_delta
  
  dr_ds_cil <- dr_delta_s - 1.96*xf_delta_s$se
  dr_ds_cih <- dr_delta_s + 1.96*xf_delta_s$se
  dr_d_cil <- dr_delta - 1.96*xf_delta$se
  dr_d_cih <- dr_delta + 1.96*xf_delta$se
  
  phi_ds <- xf_delta_s$observation_data[[1]] %>% select(u_i)
  phi_d <- xf_delta$observation_data[[1]] %>% select(u_i)
  phi <- cbind(phi_ds, phi_d) %>% as.matrix
  gdot <- c(dr_delta^(-1), dr_delta_s/dr_delta^2)
  Sigma <- t(phi) %*% phi / n
  sigma <- t(gdot) %*% Sigma %*% gdot
  r_se <- sqrt(sigma)/sqrt(n)
  
  dr_r_cil <- dr_r - 1.96*r_se
  dr_r_cih <- dr_r + 1.96*r_se
  
  deltastar_s <- map(1:500, function(s) {
    phi_ds %>%
      mutate(g = rexp(n),
             ustar_i = u_i*g) %>%
      summarise(mean(ustar_i)) %>%
      unlist
  }) %>% unlist
  deltastar <- map(1:500, function(s) {
    phi_d %>%
      mutate(g = rexp(n),
             ustar_i = u_i*g) %>%
      summarise(mean(ustar_i)) %>%
      unlist
  }) %>% unlist
  R_bci_l <- quantile(1 - deltastar_s/deltastar, 0.025)
  R_bci_h <- quantile(1 - deltastar_s/deltastar, 0.975)
  
  Delta_0 <- mean(y_1) - mean(y_0)
  R_0 <- 1 - Delta_s/Delta_0
  
  # browser()
  
  out_ds <- data.frame(
    R = dr_r,
    R_cil = dr_r_cil,
    R_cih = dr_r_cih,
    R_bci_l = R_bci_l,
    R_bci_h = R_bci_h
  ) %>% as_tibble
  if (write) {
    write_csv(out_ds, glue('/n/scratch3/users/d/dma12/doubly-robust-surrogate/ptb_n{n}-p{p}-R{R}-{run}.csv'))
  }
  out_ds
}

sim_params <- expand.grid(n = c(500, 1000),
                          p = c(50, 100, 500),
                          R = c(0.2, 0.9),
                          run = 1:1000)
# tst <- sim_params %>% filter(n < 1000) %>% sample_n(1)
# tst
# with(tst, simfn(n = 500,
#                 p = 50,
#                 R = 0.9,
#                 run = runif(1),
#                 write = FALSE))

# simfn(n = 500, p = 50, R = 0.9, write = FALSE, run = -12)


options(
  clustermq.defaults = list(ptn="medium",
                            log_file="Rout/log%a.log",
                            time_amt = "72:00:00"
  )
)
sim_res <- Q_rows(sim_params, simfn, 
                  fail_on_error = FALSE,
                  n_jobs = 100)
saveRDS(sim_res, here('results/10_ptb-cis.rds'))
