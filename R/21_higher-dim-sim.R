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
library(RCAL)
library(SIS)

simfn <- function(n, K, r, est, run = 0, write = TRUE) {
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
  library(RCAL)
  library(SIS)
  
  
  set.seed(run)
  p <- 1000
  q <- 10
  sig <- sqrt(1/p)
  Delta_s <- 5
  x <- rnorm(n*q, sd = 1) %>% matrix(n, q)
  a <- rbinom(n, prob = 0.5, size = 1)
  
  s_00 <- rbinom(n*4, size = 1, prob = 0.2) %>% matrix(n, 4)
  s_10 <- rbinom(n*4, size = 1, prob = 0.8) %>% matrix(n, 4)
  s_11 <- s_01 <- runif(n*(p-4), 0, 6) %>% matrix(n, p-4)
  
  s_0 <- cbind(s_00, s_01)
  s_1 <- cbind(s_10, s_11)
  s <- s_1*a + (1-a)*s_0
  
  y_1 <- Delta_s + rowMeans(x[,1:pmin(q, 25)]) + 5*(s_1[,1]*s_1[,2]*s_1[,10]- s_1[,3])  - s_1[,2]  + rnorm(n, sd = sig)
  y_0 <- rowMeans(x[,1:pmin(q, 25)]) + s_0[,1] + s_0[,2] + s_0[,3] + rnorm(n, sd = sig)
  y <- y_1*a + (1-a)*y_0
  # mean(y_1) - mean(y_0)
  
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
  
  if (est == 'xfdr') {
    print('fit xfdr')
    xf_surr <- map(1:r, function(xx) {
      xf_surrogate(ds = wds,
                   x = paste('x.', 1:q, sep =''),
                   s = paste('s.', 1:p, sep =''),
                   a = 'a',
                   y = 'y',
                   K = K,
                   outcome_learners = c("SL.mean", "SL.glmnet",
                                        "SL.ridge", "SL.lm", "SL.svm", 'SL.ranger'),
                   ps_learners = c("SL.mean", "SL.glmnet",
                                   "SL.glm", "SL.lda", "SL.qda",
                                   "SL.svm", 'SL.ranger'),
                   trim_at = 0.01,
                   mthd = 'superlearner', ncores = 1)
    })
    
    out_ds <- tibble(
      type = 'xfdr',
      Delta = map(xf_surr, pluck, 'deltahat') %>% unlist %>% median,
      Delta_s = map(xf_surr, pluck, 'deltahat_s') %>% unlist %>% median,
      R = 1 - Delta_s/Delta,
      R_cil = R - 1.96*map(xf_surr, pluck, 'R_se') %>% unlist %>% median,
      R_cih = R + 1.96*map(xf_surr, pluck, 'R_se') %>% unlist %>% median
    )
  } else if (est == 'xfl') {
    print('fit xfl')
    xfl_surr <- map(1:r, function(xx) {
      xf_surrogate(ds = wds,
                   x = paste('x.', 1:q, sep =''),
                   s = paste('s.', 1:p, sep =''),
                   a = 'a',
                   y = 'y',
                   K = K,
                   trim_at = 0.01,
                   mthd = 'lasso',
                   ncores = 1)
    })
    
    out_ds <- tibble(
      type = 'xfl',
      Delta = map(xfl_surr, pluck, 'deltahat') %>% unlist %>% median,
      Delta_s = map(xfl_surr, pluck, 'deltahat_s') %>% unlist %>% median,
      R = 1 - Delta_s/Delta,
      R_cil = R - 1.96*map(xfl_surr, pluck, 'R_se') %>% unlist %>% median,
      R_cih = R + 1.96*map(xfl_surr, pluck, 'R_se') %>% unlist %>% median
    )
  } else if (est == 'xfs') {
    print('fit xfs')
    xfs_surr <- map(1:r, function(xx) {
      xf_surrogate(ds = wds,
                   x = paste('x.', 1:q, sep =''),
                   s = paste('s.', 1:p, sep =''),
                   a = 'a',
                   y = 'y',
                   K = K,
                   trim_at = 0.01,
                   mthd = 'sis',
                   ncores = 1)
    })
    
    out_ds <- tibble(
      type = 'xfs',
      Delta = map(xfs_surr, pluck, 'deltahat') %>% unlist %>% median,
      Delta_s = map(xfs_surr, pluck, 'deltahat_s') %>% unlist %>% median,
      R = 1 - Delta_s/Delta,
      R_cil = R - 1.96*map(xfs_surr, pluck, 'R_se') %>% unlist %>% median,
      R_cih = R + 1.96*map(xfs_surr, pluck, 'R_se') %>% unlist %>% median)
    
  } else if (est == 'xfc') {
    print('fit xfc')
    xfc_surr <- map(1:r, function(xx) {
      xf_surrogate(ds = wds,
                   x = paste('x.', 1:q, sep =''),
                   s = paste('s.', 1:p, sep =''),
                   a = 'a',
                   y = 'y',
                   K = K,
                   trim_at = 0.01,
                   mthd = 'cal',
                   ncores = 1)
    })
    
    out_ds <- tibble(
      type = 'xfc',
      Delta = map(xfc_surr, pluck, 'deltahat') %>% unlist %>% median,
      Delta_s = map(xfc_surr, pluck, 'deltahat_s') %>% unlist %>% median,
      R = 1 - Delta_s/Delta,
      R_cil = R - 1.96*map(xfc_surr, pluck, 'R_se') %>% unlist %>% median,
      R_cih = R + 1.96*map(xfc_surr, pluck, 'R_se') %>% unlist %>% median)
    
  } else if (est == 'hima') {
    print('fit hima')
    hima_tst <- try(hima(X = a,
                         Y = y,
                         M = s,
                         COV.XM = x,
                         COV.MY = x,
                         verbose = TRUE), silent = TRUE)
    if (inherits(hima_tst, 'try-error')) {
      hima_ds <- NULL
    } else {
      hima_nie <- sum(hima_tst$`alpha*beta`)
      hima_r <- sum(hima_tst$`% total effect`/100)
      hima_nde <- (1-hima_r)*hima_nie/hima_r
      hima_delta_s <- hima_nde
      hima_delta <- hima_nde + hima_nie
      out_ds <- tibble(
        type = 'hima',
        Delta = hima_delta,
        Delta_s = hima_delta_s,
        R = hima_r,
        R_cil = NA,
        R_cih = NA
      )
    }
  } else if (est == 'fb') {
    print('fit fb')
    fb_tst <- freebird::hilma(Y = y,
                              G = s,
                              S = a %>% matrix(n, 1))
    fb_nie <- fb_tst$beta_hat
    fb_nde <- fb_tst$alpha1_hat
    fb_delta_s <- fb_nde
    fb_delta <- fb_nde + fb_nie
    fb_r <- fb_nie/fb_delta
    
    out_ds <- tibble(
      type = 'fb',
      Delta = fb_delta,
      Delta_s = fb_delta_s,
      R = fb_r,
      R_cil = NA,
      R_cih = NA
    )
  }
  Delta_0 <- mean(y_1) - mean(y_0)
  y_11 <- Delta_s + rowMeans(x[,1:pmin(q, 25)]) + 5*(s_1[,1]*s_1[,2]*s_1[,10]- s_1[,3])  - s_1[,2]  + rnorm(n, sd = sig)
  y_00 <- rowMeans(x[,1:pmin(q, 25)]) + s_0[,1] + s_0[,2] + s_0[,3] + rnorm(n, sd = sig)
  y_10 <- Delta_s + rowMeans(x[,1:pmin(q, 25)]) + 5*(s_0[,1]*s_0[,2]*s_0[,10]- s_0[,3])  - s_0[,2]  + rnorm(n, sd = sig)
  y_01 <- rowMeans(x[,1:pmin(q, 25)]) + s_1[,1] + s_1[,2] + s_1[,3] + rnorm(n, sd = sig)
  Delta_s0 <- mean(y_11 - y_01)/2 + mean(y_10 - y_00)/2
  R_0 <- 1 - Delta_s0/Delta_0
  

  
  ds_0 <- tibble(
    type = 'true',
    Delta = Delta_0,
    Delta_s = Delta_s0,
    R = R_0,
    R_cil = NA,
    R_cih = NA
  )
  
  
 
  
  
  
  
  
  
  # out_ds <- full_join(
  #   ds_0,
  #   xf_ds
  # ) %>%
  #   full_join(xfl_ds) %>%
  #   full_join(xfs_ds) %>%
  #   full_join(xfc_ds) %>%
  #   full_join(fb_ds)
  # if (!is.null(hima_ds)) {
  #   out_ds <- out_ds %>%
  #     full_join(hima_ds)
  # }
  
  out_ds <- out_ds %>%
    full_join(ds_0) %>%
    mutate(n = n,
           K = K,
           r = r,
           run = run)
  
  # browser()
  if (write) {
    write_csv(out_ds, glue('/n/scratch3/users/d/dma12/doubly-robust-surrogate/res1_n{n}-k{K}-r{r}_{est}_{run}.csv'))
  }
  out_ds
}

sim_params <- expand.grid(n = 200,
                          K = c(2, 10),
                          r = c(1, 5, 11),
                          est = c('xfdr',
                                  'xfc',
                                  'xfl',
                                  'xfs',
                                  'hima',
                                  'fb'),
                          run = 1:10) %>%
  filter(!(K == 20 & r == 11))
# tst <- sim_params %>% sample_n(1)
# tst
# with(tst, simfn(n = 200,
#                 K = 2,
#                 r = 1,
#                 est = est,
#                 run = run,
#                 write = TRUE))

options(
  clustermq.defaults = list(ptn="short",
                            log_file="Rout/log%a.log",
                            time_amt = "12:00:00"
  )
)
sim_res <- Q_rows(sim_params, simfn, 
                  fail_on_error = FALSE,
                  n_jobs = 10)
saveRDS(sim_res, here('results/21_higher-dim-sim_results.rds'))
