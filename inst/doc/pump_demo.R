## ----initialize, include = FALSE----------------------------------------------
library( PUMP )
recompile <- FALSE

knitr::opts_chunk$set(
  cache = FALSE,
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE,
  fig.width = 7,
  fig.height = 4,
  fig.align = "center"
)
options(knitr.kable.NA = '')
library( dplyr )
library( ggplot2 )
library( here )
library( kableExtra )
library( knitr )
set.seed( 0440044 )

## ---- echo = FALSE------------------------------------------------------------
info <- pump_info(comment = FALSE)
kable( info$Context, format = 'latex', booktabs = TRUE ) %>%
  kable_styling(position = "center")

## ---- echo = FALSE------------------------------------------------------------
kable(info$Parameters, format = 'latex', booktabs = TRUE ) %>%
  kable_styling(position = "center")

## ----MDEScalc, eval = FALSE---------------------------------------------------
#  m <- pump_mdes(
#    d_m = "d3.2_m3fc2rc",         # choice of design and analysis strategy
#    MTP = "HO",                   # multiple testing procedure
#    target.power = 0.80,          # desired power
#    power.definition = "D1indiv", # power type
#    M = 5,                        # number of outcomes
#    J = 3,                        # number of schools per block
#    K = 21,                       # number districts
#    nbar = 258,                   # average number of students per school
#    Tbar = 0.50,                  # prop treated
#    alpha = 0.05,                 # significance level
#    numCovar.1 = 5,               # number of covariates at level 1
#    numCovar.2 = 3,               # number of covariates at level 2
#    R2.1 = 0.1, R2.2 = 0.7,       # explanatory power of covariates for each level
#    ICC.2 = 0.05, ICC.3 = 0.4,    # intraclass correlation coefficients
#    rho = 0.4 )                   # how correlated outcomes are

## ----MDEScalc-compile, echo = FALSE-------------------------------------------
if(recompile)
{
    m <- pump_mdes(
      d_m = "d3.2_m3fc2rc",         # choice of design and analysis strategy
      MTP = "HO",                   # multiple testing procedure
      target.power = 0.80,          # desired power
      power.definition = "D1indiv", # power type
      M = 5,                        # number of outcomes
      J = 3,                        # number of schools per block
      K = 21,                       # number districts
      nbar = 258,                   # average number of students per school
      Tbar = 0.50,                  # prop treated
      alpha = 0.05,                 # significance level
      numCovar.1 = 5,               # number of covariates at level 1
      numCovar.2 = 3,               # number of covariates at level 2
      R2.1 = 0.1, R2.2 = 0.7,       # explanatory power of covariates for each level
      ICC.2 = 0.05, ICC.3 = 0.4,    # intraclass correlation coefficients
      rho = 0.4 ) 
    
    saveRDS(m, here::here("vignettes/output", "MDEScalc.RDS"))
} else
{
    m <- readRDS(here::here("vignettes/output", "MDEScalc.RDS"))
}

## ----echo = TRUE--------------------------------------------------------------
knitr::kable( m, digits = 3 ) %>%
  kableExtra::kable_styling( position = "center" )

## ----MDEScalcmin1, echo = TRUE, eval = FALSE----------------------------------
#  m2 <- update( m, power.definition = "min1" )

## ----MDEScalcmin1-compile, echo = FALSE---------------------------------------
if(recompile)
{
    m2 <- update( m, power.definition = "min1" )
    saveRDS(m2, here::here("vignettes/output", "MDEScalcmin1.RDS"))
} else
{
    m2 <- readRDS(here::here("vignettes/output", "MDEScalcmin1.RDS"))
}
print( m2 )

## ----MDESwithNumZero, echo = TRUE, eval = FALSE-------------------------------
#  m3 <- update( m2, numZero = 2 )

## ----MDESwithNumZero-compile, echo = FALSE------------------------------------
if(recompile)
{
    m3 <- update( m2, numZero = 2 )
    saveRDS(m3, here::here("vignettes/output", "MDESwithNumZero.RDS"))
} else
{
    m3 <- readRDS(here::here("vignettes/output", "MDESwithNumZero.RDS"))
}
print( m3 )

## ----samplesizecalc, eval = FALSE---------------------------------------------
#  smp <- pump_sample(
#    d_m = "d3.2_m3fc2rc",
#    MTP = "HO",
#    typesample = "K",
#    target.power = 0.80, power.definition = "min1", tol = 0.01,
#    MDES = 0.10, M = 5, nbar = 258, J = 3,
#    Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 3,
#    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.40, rho = 0.4 )
#  
#  print( smp )

## ----samplesizecalc-compile, echo = FALSE-------------------------------------
if(recompile)
{
    smp <- pump_sample(
      d_m = "d3.2_m3fc2rc",
      MTP = "HO",
      typesample = "K",
      target.power = 0.80, power.definition = "min1", tol = 0.01,
      MDES = 0.10, M = 5, nbar = 258, J = 3,
      Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 3,
      R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.40, rho = 0.4 )
    saveRDS(smp, here::here("vignettes/output", "samplesizecalc.RDS"))
} else
{
    smp <- readRDS(here::here("vignettes/output", "samplesizecalc.RDS"))
}
print( smp )

## ----samplesizeverify, eval = FALSE-------------------------------------------
#  p_check <- update( smp, type = "power", tnum = 20000,
#                     long.table = TRUE )

## ----samplesizeverify-compile, echo = FALSE-----------------------------------
if(recompile)
{
    p_check <- update( smp, type = "power", tnum = 20000,
                   long.table = TRUE )
    saveRDS(p_check, here::here("vignettes/output", "samplesizeverify.RDS"))
} else
{
    p_check<- readRDS(here::here("vignettes/output", "samplesizeverify.RDS"))
}
knitr::kable( p_check, digits = 2 ) %>%
  kableExtra::kable_styling(position = "center")

## ----plotsamplepowercurve, fig.height = 3.5, fig.width=5, echo=TRUE-----------
plot( smp )

## ----powbase, eval = FALSE----------------------------------------------------
#  pow <- update( p_check, tnum = 10000 )

## ----powbase-compile, echo = FALSE--------------------------------------------
if(recompile)
{
    pow <- update( p_check, tnum = 10000 )
    saveRDS(pow, here::here("vignettes/output", "powbase.RDS"))
} else
{
    pow <- readRDS(here::here("vignettes/output", "powbase.RDS"))
}

## ----othercorrections, eval = FALSE-------------------------------------------
#  p2 <- update( pow,
#                MTP = c( "BF", "HO", "WY-SD" ) )
#  plot( p2 )

## ----othercorrections-compile, echo = FALSE-----------------------------------
if(recompile)
{
    p2 <- update( pow,
      MTP = c( "BF", "HO", "WY-SD" ) )
    saveRDS(p2, here::here("vignettes/output", "othercorrections.RDS"))
} else
{
    p2 <- readRDS(here::here("vignettes/output", "othercorrections.RDS"))
}
plot( p2 )

## ----powICC, eval = FALSE-----------------------------------------------------
#  p_b <- update( pow, ICC.2 = 0.20, ICC.3 = 0.25 )
#  print( p_b )

## ----powICC-compile, echo = FALSE---------------------------------------------
if(recompile)
{
    p_b <- update( pow, ICC.2 = 0.20, ICC.3 = 0.25 )
    saveRDS(p_b, here::here("vignettes/output", "powICC.RDS"))
} else
{
    p_b <- readRDS(here::here("vignettes/output", "powICC.RDS"))
}
print( p_b )

## ----powR2, eval = FALSE------------------------------------------------------
#  p_d <- update( pow,
#            	   R2.1 = c( 0.1, 0.3, 0.1, 0.2, 0.2 ),
#            	   R2.2 = c( 0.4, 0.8, 0.3, 0.2, 0.2 ) )
#  print( p_d )

## ----powR2-compile, echo = FALSE----------------------------------------------
if(recompile)
{
    p_d <- update( pow,
          	   R2.1 = c( 0.1, 0.3, 0.1, 0.2, 0.2 ),
          	   R2.2 = c( 0.4, 0.8, 0.3, 0.2, 0.2 ) )
    saveRDS(p_d, here::here("vignettes/output", "powR2.RDS"))
} else
{
    p_d <- readRDS(here::here("vignettes/output", "powR2.RDS"))
}
print( p_d )

## -----------------------------------------------------------------------------
summary( p_d )

## ----ICCgrid, eval = FALSE----------------------------------------------------
#  grid <- update_grid( pow,
#          	ICC.2 = seq( 0, 0.3, 0.05 ),
#          	ICC.3 = seq( 0, 0.60, 0.20 ) )
#  
#  plot( grid, power.definition = "min1" )

## ----ICCgrid-compile, echo = FALSE--------------------------------------------
if(recompile)
{
    grid <- update_grid( pow,
        	ICC.2 = seq( 0, 0.3, 0.05 ),
        	ICC.3 = seq( 0, 0.60, 0.20 ) )
    saveRDS(grid, here::here("vignettes/output", "ICCgrid.RDS"))
} else
{
    grid <- readRDS(here::here("vignettes/output", "ICCgrid.RDS"))
}
plot( grid, power.definition = "min1" )

## ----rhogrid, eval = FALSE----------------------------------------------------
#  gridRho <- update_grid( pow,
#          	  MTP = c( "BF", "WY-SD" ),
#          	  rho = seq( 0, 0.9, by = 0.15 ),
#          	  tnum = 1000,
#          	  B = 3000 )

## ----rhogrid-compile, echo = FALSE--------------------------------------------
if(recompile)
{
    gridRho <- update_grid( pow,
        	  MTP = c( "BF", "WY-SD" ),
        	  rho = seq( 0, 0.9, by = 0.15 ),
        	  tnum = 1000,
        	  B = 3000 )
    saveRDS(gridRho, here::here("vignettes/output", "rhogrid.RDS"))
} else
{
    gridRho <- readRDS(here::here("vignettes/output", "rhogrid.RDS"))
}

## -----------------------------------------------------------------------------
plot( gridRho )

## ----numzerogrid, fig.height = 2.5, eval = FALSE------------------------------
#  gridZero <- update_grid( pow,
#          	             numZero = 0:4,
#                           M = 5 )
#  plot( gridZero, nrow = 1 )

## ----numzerogrid-compile, echo = FALSE, fig.height = 2.5----------------------
if(recompile)
{
    gridZero <- update_grid( pow,
                             numZero = 0:4,
                             M = 5 )
    saveRDS(gridZero, here::here("vignettes/output", "numzerogrid.RDS"))
} else
{
   gridZero <- readRDS(here::here("vignettes/output", "numzerogrid.RDS"))
}
plot( gridZero, nrow = 1 )

