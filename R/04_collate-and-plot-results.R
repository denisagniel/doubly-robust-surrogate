library(tidyverse)
library(here)
library(glue)

res <- read_rds(here('results/01_sim-results.rds'))
sim_params <- expand.grid(n = c(100, 500, 1000),
                          p = c(50, 100, 500),
                          R_0n = c(0.2, 0.9),
                          run = 1:1000)
sim_params <- sim_params %>%
  mutate(sim = as.character(1:nrow(sim_params)))
resds <- res %>%
  bind_rows(.id = 'sim') %>%
  inner_join(sim_params) %>%
  mutate(Delta = case_when(
    type == 'true' & R_0n == 0.2 ~ 3.05,
    type == 'true' & R_0n == 0.9 ~ 1.475,
    TRUE ~ Delta),
    R = ifelse(type == 'true', 1 - Delta_s/Delta, R))

true_vals <- resds %>%
  filter(type == 'true') %>%
  select(R_0n, Delta_0 = Delta, 
         Delta_s0 = Delta_s, 
         R_0 = R) %>%
  unique

resds %>% sample_n(1) %>% select(sim) %>% inner_join(resds)


ggplot(resds, aes(x = R, y = type)) +
  ggridges::geom_density_ridges() +
  scale_x_continuous(limits = c(0, 1)) +
  facet_wrap(~ n + p + R_0n)

res_sum <- resds %>%
  group_by(n, p, R_0n, type) %>%
  summarise_at(vars(Delta, Delta_s, R), list(m = ~median(.),
                                             lo = ~quantile(., 0.025),
                                             hi = ~quantile(., 0.975))) %>%
  group_by(n, p, R_0n) %>%
  ungroup %>%
  inner_join(true_vals) %>%
  filter(type != 'true')

res_sum %>% sample_n(1) %>% select(n, p, R_0n) %>% inner_join(res_sum)

delta_pl <- 
  ggplot(res_sum,
         aes(x = type,
             y = Delta_m)) +
  geom_point() +
  geom_linerange(aes(ymin = Delta_lo,
                     ymax = Delta_hi)) +
  facet_wrap(~ n + p + R_0n, scales = 'free') +
  geom_hline(aes(yintercept = Delta_0), linetype = 2) +
  theme_bw()
delta_pl

delta_spl <- 
  ggplot(res_sum,
         aes(x = type,
             y = Delta_s_m)) +
  geom_point() +
  geom_linerange(aes(ymin = Delta_s_lo,
                     ymax = Delta_s_hi)) +
  facet_wrap(~ n + p + R_0n, scales = 'free') +
  geom_hline(aes(yintercept = Delta_s0), linetype = 2) +
  theme_bw()
delta_spl

r_pl <- 
  ggplot(res_sum %>%
           filter(n == 500) %>%
           mutate(r0 = glue('R = {round(R_0, 2)}'),
                  pp = glue('p = {p}'),
                  ppp = fct_relevel(pp, 'p = 50', 'p = 100', 'p = 500')),
         aes(x = type,
             y = R_m)) +
  geom_point() +
  geom_linerange(aes(ymin = R_lo,
                     ymax = R_hi)) +
  facet_wrap(~ ppp + r0, scales = 'free') +
  geom_hline(aes(yintercept = R_0), linetype = 2) +
  theme_bw() +
  ggtitle(expression(hat(R)),
          subtitle = ('n = 500'))
r_pl
ggsave(here('results/initial-r-plot.png'),
       height = 6, width = 6)
  ggplot(res_sum %>%
           filter(type == 'xfdr'),
         aes(x = type,
             y = R_m)) +
  geom_point() +
  geom_linerange(aes(ymin = R_lo,
                     ymax = R_hi)) +
  facet_wrap(~ n + p + R_0n, scales = 'free') +
  geom_hline(aes(yintercept = R_0), linetype = 2) +
  theme_bw()
