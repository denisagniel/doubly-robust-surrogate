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
  
  xfl_delta_s <- xfit_dr(ds = wds,
                         xvars = c(paste('s.', 1:p, sep =''),
                                   paste('x.', 1:q, sep ='')),
                         yname = y,
                         aname = a,
                         K = 4,
                         trim_at = 0.01,
                         mthd = 'lasso')
  
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
                      mthd = 'superlearner')
  xfl_delta <- xfit_dr(ds = wds,
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
                       mthd = 'lasso')
  
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
  
  drl_delta_s <- xfl_delta_s$estimate
  drl_delta <- xfl_delta$estimate
  drl_r <- 1 - drl_delta_s/drl_delta
  
  drl_ds_cil <- drl_delta_s - 1.96*xfl_delta_s$se
  drl_ds_cih <- drl_delta_s + 1.96*xfl_delta_s$se
  drl_d_cil <- drl_delta - 1.96*xfl_delta$se
  drl_d_cih <- drl_delta + 1.96*xfl_delta$se
  
  # browser()
  
  phi_dsl <- xfl_delta_s$observation_data[[1]] %>% select(u_i)
  phi_dl <- xfl_delta$observation_data[[1]] %>% select(u_i)
  phil <- cbind(phi_dsl, phi_dl) %>% as.matrix
  gdotl <- c(drl_delta^(-1), drl_delta_s/drl_delta^2)
  Sigmal <- t(phil) %*% phil / n
  sigmal <- t(gdotl) %*% Sigmal %*% gdotl
  rl_se <- sqrt(sigmal)/sqrt(n)
  
  drl_r_cil <- drl_r - 1.96*rl_se
  drl_r_cih <- drl_r + 1.96*rl_se
  
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
  # browser()
  bama_ds_cil <- quantile(bama_tst$beta.a, 0.025)
  bama_ds_cih <- quantile(bama_tst$beta.a, 0.975)
  bama_d_cil <- quantile(bama_tst$beta.a + rowSums(bama_tst$alpha.a* bama_tst$beta.m), 0.025)
  bama_d_cih <- quantile(bama_tst$beta.a + rowSums(bama_tst$alpha.a* bama_tst$beta.m), 0.975)
  bama_r_cih <- 1 - quantile(bama_tst$beta.a/(bama_tst$beta.a + rowSums(bama_tst$alpha.a* bama_tst$beta.m)), 0.025)
  bama_r_cil <- 1 - quantile(bama_tst$beta.a/(bama_tst$beta.a + rowSums(bama_tst$alpha.a* bama_tst$beta.m)), 0.975)
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
    Delta = c(Delta, dr_delta, drl_delta, bama_delta, hima_delta, fb_delta),
    # Delta = c(Delta_0, dr_delta, drl_delta, bama_delta, hima_delta, fb_delta),
    Delta_cil = c(Delta_0, dr_d_cil, drl_d_cil, bama_d_cil, NA, NA),
    Delta_cih = c(Delta_0, dr_d_cih, drl_d_cih, bama_d_cih, NA, NA),
    Delta_s = c(Delta_s, dr_delta_s, drl_delta_s, bama_delta_s, hima_delta_s, fb_delta_s),
    Delta_s_cil = c(Delta_s, dr_ds_cil, drl_ds_cil, bama_ds_cil, NA, NA),
    Delta_s_cih = c(Delta_s, dr_ds_cih, drl_ds_cih, bama_ds_cih, NA, NA),
    R = c(R, dr_r, drl_r, bama_r, hima_r, fb_r),
    # R = c(R_0, dr_r, drl_r, bama_r, hima_r, fb_r),
    R_cil = c(R, dr_r_cil, drl_r_cil, bama_r_cil, NA, NA),
    R_cih = c(R, dr_r_cih, drl_r_cih, bama_r_cih, NA, NA)
  ) %>% as_tibble
  if (write) {
    write_csv(out_ds, glue('/n/scratch3/users/d/dma12/doubly-robust-surrogate/res1_n{n}-p{p}-R{R}-{run}.csv'))
  }
  out_ds
}

sim_params <- expand.grid(n = 500,
                          p = c(50, 100, 500),
                          R = c(0.2, 0.9),
                          run = 1:3)
tst <- sim_params %>% filter(n < 1000) %>% sample_n(1)
tst
with(tst, simfn(n = n,
                p = p,
                R = R,
                run = run,
                write = FALSE))

simfn(n = 500, p = 50, R = 0.9, write = FALSE, run = -12)


options(
  clustermq.defaults = list(ptn="short",
                            log_file="Rout/log%a.log",
                            time_amt = "12:00:00"
  )
)
sim_res <- Q_rows(sim_params, simfn, 
                  fail_on_error = FALSE,
                  n_jobs = 200)
saveRDS(sim_res, here('results/01_sim-results.rds'))
