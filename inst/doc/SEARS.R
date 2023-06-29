## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, warning = FALSE---------------------------------------------------
library(SEARS)

## -----------------------------------------------------------------------------
SEARS(p.p = 0.2, p.d = c(0.2, 0.2, 0.2, 0.2, 0.2), p.tox = c(0.03, 0.06, 0.17, 0.3, 0.5),
      k1 = 0.95, k2 = 0.8, pi_t = 0.17, pi_e = 0.2, pT = 0.17, csize = 3, csize2 = 3, 
      d.cs = 36, p.cs = 36, phase1_size = 30, n_earlystop = 18, Nsim = 1000, 
      n_catchup = 3, control_arm = "fixed", power_c = 0.5, lower_bound = 0.05, 
      weight1 = 0.5, weight2 = 0.5, seed = 100)

## -----------------------------------------------------------------------------
SEARS(p.p = 0.2, p.d = c(0.1, 0.5, 0.6, 0.7, 0.8), p.tox = c(0.03, 0.06, 0.17, 0.3, 0.5),
      k1 = 0.95, k2 = 0.8, pi_t = 0.17, pi_e = 0.2, pT = 0.17, csize = 3, csize2 = 3, 
      d.cs = 36, p.cs = 36, phase1_size = 30, n_earlystop = 18, Nsim = 1000, 
      n_catchup = 3, control_arm = "fixed", power_c = 0.5, lower_bound = 0.05, 
      weight1 = 0, weight2 = 0, seed = 100)

## -----------------------------------------------------------------------------
SEARS(p.p = 0.2, p.d = c(0.1, 0.5, 0.6, 0.7, 0.8), p.tox = c(0.03, 0.06, 0.17, 0.3, 0.5),
      k1 = 0.95, k2 = 0.8, pi_t = 0.17, pi_e = 0.2, pT = 0.17, csize = 3, csize2 = 3, 
      d.cs = 36, p.cs = 36, phase1_size = 30, n_earlystop = 18, Nsim = 1000, 
      n_catchup = 3, control_arm = "fixed", power_c = 0.5, lower_bound = 0.05, 
      weight1 = 0.5, weight2 = 0.5, seed = 100)

## -----------------------------------------------------------------------------
SEARS(p.p = 0.2, p.d = c(0.01, 0.05, 0.6, 0.7, 0.8), p.tox = c(0.03, 0.06, 0.12, 0.3, 0.5),
      k1 = 0.95, k2 = 0.8, pi_t = 0.17, pi_e = 0.2, pT = 0.17, csize = 3, csize2 = 3, 
      d.cs = 36, p.cs = 36, phase1_size = 30, n_earlystop = 18, Nsim = 1000, n_catchup = 3, 
      control_arm = "fixed", power_c = 0.5, lower_bound = 0.05, weight1 = 0.5, 
      weight2 = 0.5, seed = 100)

## -----------------------------------------------------------------------------
SEARS(p.p = 0.2, p.d = c(0.01, 0.6, 0.65, 0.7, 0.8), p.tox = c(0.03, 0.06, 0.12, 0.3, 0.5),
      k1 = 0.95, k2 = 0.8, pi_t = 0.17, pi_e = 0.2, pT = 0.17, csize = 3, csize2 = 3, 
      d.cs = 36, p.cs = 36, phase1_size = 30, n_earlystop = 18, Nsim = 1000, 
      n_catchup = 10, control_arm = "", power_c = 0.5, lower_bound = 0.05, weight1 = 0.5, 
      weight2 = 0.5, seed = 100)

## -----------------------------------------------------------------------------
SEARS(p.p = 0.2, p.d = c(0.01, 0.6, 0.65, 0.7, 0.8), p.tox = c(0.03, 0.06, 0.12, 0.3, 0.5),
      k1 = 0.95, k2 = 0.8, pi_t = 0.17, pi_e = 0.2, pT = 0.17, csize = 3, csize2 = 3, 
      d.cs = 36, p.cs = 36, phase1_size = 30, n_earlystop = 18, Nsim = 1000, 
      n_catchup = 10, control_arm = "fixed", power_c = 0.5, lower_bound = 0.05, 
      weight1 = 0.5, weight2 = 0.5, seed = 100)

## -----------------------------------------------------------------------------
SEARS(p.p = 0.2, p.d = c(0.01, 0.6, 0.65, 0.7, 0.8), p.tox = c(0.03, 0.06, 0.12, 0.3, 0.5),
      k1 = 0.95, k2 = 0.8, pi_t = 0.17, pi_e = 0.2, pT = 0.17, csize = 3, csize2 = 3, 
      d.cs = 36, p.cs = 36, phase1_size = 30, n_earlystop = 18, Nsim = 1000, 
      n_catchup = 10, control_arm = "fixed", power_c = 0, lower_bound = 0, 
      weight1 = 0.5, weight2 = 0.5, seed = 100)

