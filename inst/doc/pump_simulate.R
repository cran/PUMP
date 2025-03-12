## ----initialize, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  warning = FALSE,
  message = FALSE
)
library( knitr )
library( PUMP )

set.seed( 524235326 )

## ----gen.data-----------------------------------------------------------------
pp <- pump_power( "d3.1_m3rr2rr", MDES = 0.2, 
                  M = 5, rho = 0.8,
                  MTP = "BH",
                  nbar = 30, J = 7, K = 5, Tbar = 0.5 )
sim.data <- gen_sim_data( pp )

## ----first.dataset------------------------------------------------------------
head( sim.data[[1]] )

## ----single.outcome, warning=FALSE--------------------------------------------
pp.one <- update( pp, M = 1 )
sim3 <- gen_sim_data( pp.one )
head( sim3 )

## ----sep.data-----------------------------------------------------------------
sim.data.v2 <- gen_sim_data( pp, return.as.dataframe = FALSE )
names( sim.data.v2 )

## ----model.params-------------------------------------------------------------
model.params.list <- list(
  M = 3                             # number of outcomes
  , J = 7                           # number of schools
  , K = 5                           # number of districts
                                    # (for two-level model, set K = 1)
  , nbar = 30                       # number of individuals per school
  , rho.default = 0.5               # default rho value (optional)
  ################################################## impact
  , MDES = 0.125                    # minimum detectable effect size      
  ################################################## level 3: districts
  , R2.3 = 0.1                      # percent of district variation
                                      # explained by district covariates
  , ICC.3 = 0.2                     # district intraclass correlation
  , omega.3 = 0.1                   # ratio of district effect size variability
                                      # to random effects variability
  ################################################## level 2: schools
  , R2.2 = 0.1                      # percent of school variation
                                    # explained by school covariates
  , ICC.2 = 0.2                     # school intraclass correlation	
  , omega.2 = 0.1                   # ratio of school effect size variability
                                      # to random effects variability
  ################################################## level 1: individuals
  , R2.1 = 0.1                      # percent of indiv variation explained
                                      # by indiv covariates
)

## ----model.params.full, eval = FALSE------------------------------------------
# M <- 3
# rho.default <- 0.5
# default.rho.matrix <- gen_corr_matrix(M = M, rho.scalar = rho.default)
# default.kappa.matrix <- matrix(0, M, M)
# 
# model.params.list <- list(
#   M = 3                             # number of outcomes
#   , J = 7                           # number of schools
#   , K = 5                           # number of districts
#                                     # (for two-level model, set K = 1)
#   , nbar = 30                       # number of individuals per school
#   , S.id = NULL                     # N-length vector of school assignments
#   , D.id = NULL                     # N-length vector of district assignments
#   ################################################## grand mean outcome and impact
#   , Xi0 = 0                         # scalar grand mean outcome under no treatment
#   , MDES = rep(0.125, M)            # minimum detectable effect size
#   ################################################## level 3: districts
#   , R2.3 = rep(0.1, M)              # percent of district variation
#                                       # explained by district covariates
#   , rho.V = default.rho.matrix      # MxM correlation matrix of district covariates
#   , ICC.3 = rep(0.2, M)             # district intraclass correlation
#   , omega.3 = rep(0.1, M)           # ratio of district effect size variability
#                                       # to random effects variability
#   , rho.w0 = default.rho.matrix     # MxM matrix of correlations for district random effects
#   , rho.w1 = default.rho.matrix     # MxM matrix of correlations for district impacts
#   , kappa.w =  default.kappa.matrix # MxM matrix of correlations between district
#                                       # random effects and impacts
#   ################################################## level 2: schools
#   , R2.2 = rep(0.1, M)              # percent of school variation
#                                       # explained by school covariates
#   , rho.X = default.rho.matrix      # MxM correlation matrix of school covariates
#   , ICC.2 = rep(0.2, M)             # school intraclass correlation	
#   , omega.2 = rep(0.1, M)           # ratio of school effect size variability
#                                       # to random effects variability
#   , rho.u0 = default.rho.matrix     # MxM matrix of correlations for school random effects
#   , rho.u1 = default.rho.matrix     # MxM matrix of correlations for school impacts
#   , kappa.u = default.kappa.matrix  # MxM matrix of correlations between school
#                                       # random effects and impacts
#   ################################################## level 1: individuals
#   , R2.1 = rep(0.1, M)              # percent of indiv variation explained
#                                       # by indiv covariates
#   , rho.C = default.rho.matrix      # MxM correlation matrix of individual covariates
#   , rho.r = default.rho.matrix      # MxM matrix of correlations for individual residuals
# )

## ----gen.sim.data-------------------------------------------------------------
sim.data <- gen_sim_data(d_m = 'd3.3_m3rc2rc', model.params.list, Tbar = 0.5)

## ----convert.params-----------------------------------------------------------
dgp.params.list <- convert_params(model.params.list)

## ----gen.full.data------------------------------------------------------------
sim.data <- gen_base_sim_data(dgp.params.list, 
                              dgp.params = TRUE,
                              return.as.dataframe = FALSE )

## ----tx-----------------------------------------------------------------------
d_m <- 'd3.3_m3rc2rc'
sim.data$T.x <- gen_T.x(
    d_m = d_m,
    S.id = sim.data$ID$S.id,
    D.id = sim.data$ID$D.id,
    Tbar = 0.5
)
sim.data$Yobs <- gen_Yobs(sim.data, T.x = sim.data$T.x)

## ----convert.dataframe--------------------------------------------------------
sim.data <- PUMP:::makelist_samp( sim.data )

