## ----initialize, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  warning = FALSE,
  message = FALSE
)
library( knitr )
library( PUMP )

set.seed( 524235326 )

## -----------------------------------------------------------------------------
model.params.list <- list(
  M = 3                             # number of outcomes
  , J = 30                          # number of schools
  , K = 10                          # number of districts
                                      # (for two-level model, set K = 1)
  , nbar = 50                       # number of individuals per school
  , rho.default = 0.5               # default rho value (optional)
  ################################################## impact
  , MDES = 0.125                    # minimum detectable effect size      
  ################################################## level 3: districts
  , numCovar.3 = 1                  # number of district covariates
  , R2.3 = 0.1                      # percent of district variation
                                      # explained by district covariates
  , ICC.3 = 0.2                     # district intraclass correlation
  , omega.3 = 0.1                   # ratio of district effect size variability
                                      # to random effects variability
  ################################################## level 2: schools
  , numCovar.2 = 1                  # number of school covariates
  , R2.2 = 0.1                      # percent of school variation
                                    # explained by school covariates
  , ICC.2 = 0.2                     # school intraclass correlation	
  , omega.2 = 0.1                   # ratio of school effect size variability
                                      # to random effects variability
  ################################################## level 1: individuals
  , numCovar.1 = 1                  # number of individual covariates
  , R2.1 = 0.1                      # percent of indiv variation explained
                                      # by indiv covariates
)

## ---- eval = FALSE------------------------------------------------------------
#  M <- 3
#  rho.default <- 0.5
#  default.rho.matrix <- gen_corr_matrix(M = M, rho.scalar = rho.default)
#  default.kappa.matrix <- matrix(0, M, M)
#  
#  model.params.list <- list(
#    M = 3                             # number of outcomes
#    , J = 30                          # number of schools
#    , K = 10                          # number of districts
#                                        # (for two-level model, set K = 1)
#    , nbar = 50                       # number of individuals per school
#    , S.id = NULL                     # N-length vector of school assignments
#    , D.id = NULL                     # N-length vector of district assignments
#    ################################################## grand mean outcome and impact
#    , Xi0 = 0                         # scalar grand mean outcome under no treatment
#    , MDES = rep(0.125, M)            # minimum detectable effect size
#    ################################################## level 3: districts
#    , numCovar.3 = 1                  # number of district covariates
#    , R2.3 = rep(0.1, M)              # percent of district variation
#                                        # explained by district covariates
#    , rho.V = default.rho.matrix      # MxM correlation matrix of district covariates
#    , ICC.3 = rep(0.2, M)             # district intraclass correlation
#    , omega.3 = rep(0.1, M)           # ratio of district effect size variability
#                                        # to random effects variability
#    , rho.w0 = default.rho.matrix     # MxM matrix of correlations for district random effects
#    , rho.w1 = default.rho.matrix     # MxM matrix of correlations for district impacts
#    , kappa.w =  default.kappa.matrix # MxM matrix of correlations between district
#                                        # random effects and impacts
#    ################################################## level 2: schools
#    , numCovar.2 = 1                  # number of school covariates
#    , R2.2 = rep(0.1, M)              # percent of school variation
#                                        # explained by school covariates
#    , rho.X = default.rho.matrix      # MxM correlation matrix of school covariates
#    , ICC.2 = rep(0.2, M)             # school intraclass correlation	
#    , omega.2 = rep(0.1, M)           # ratio of school effect size variability
#                                        # to random effects variability
#    , rho.u0 = default.rho.matrix     # MxM matrix of correlations for school random effects
#    , rho.u1 = default.rho.matrix     # MxM matrix of correlations for school impacts
#    , kappa.u = default.kappa.matrix  # MxM matrix of correlations between school
#                                        # random effects and impacts
#    ################################################## level 1: individuals
#    , numCovar.1 = 1                  # number of individual covariates
#    , R2.1 = rep(0.1, M)              # percent of indiv variation explained
#                                        # by indiv covariates
#    , rho.C = default.rho.matrix      # MxM correlation matrix of individual covariates
#    , rho.r = default.rho.matrix      # MxM matrix of correlations for individual residuals
#  )

## -----------------------------------------------------------------------------
d_m <- 'd3.3_m3rc2rc'
dgp.params.list <- convert_params(model.params.list)
sim.data.full <- gen_full_data(dgp.params.list)

## -----------------------------------------------------------------------------
T.x <- gen_T.x(
    d_m = d_m, S.id = sim.data.full$ID$S.id, D.id = sim.data.full$ID$D.id,
    nbar = model.params.list$nbar, Tbar = 0.5
)
sim.data.obs <- sim.data.full
sim.data.obs$Yobs <- gen_Yobs(sim.data.full, T.x)

## ---- echo = FALSE------------------------------------------------------------
str(sim.data.obs)

