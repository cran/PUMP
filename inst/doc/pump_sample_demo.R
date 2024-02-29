## ----initialize, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  cache = FALSE,
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE
)
library( knitr )
library( PUMP )

set.seed( 524235326 )

## ----echo = FALSE-------------------------------------------------------------
set.seed( 101010 )

## -----------------------------------------------------------------------------
p <- pump_power(
    d_m = "d3.1_m3rr2rr",
    MTP = "HO",
    nbar = 50,
    K = 15,
    J = 20,
    M = 3,
    MDES = rep(0.125, 3),
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.1,
    ICC.2 = 0.2, ICC.3 = 0.2,
    omega.2 = 0.2, omega.3 = 0.2,
    rho = 0.5, tnum = 100000
)

## ----echo = FALSE-------------------------------------------------------------
knitr::kable(p)
target_power <- round( p$indiv.mean[2], digits = 3 )

## ----echo = FALSE-------------------------------------------------------------
set.seed( 40444040 )

## -----------------------------------------------------------------------------
K <- pump_sample(
  d_m = "d3.1_m3rr2rr",
  typesample = "K",
  MTP = "HO",
  target.power = target_power,
  power.definition = "D1indiv",
  J = 20,
  nbar = 50,
  M = 3,
  MDES = 0.125,
  Tbar = 0.5, alpha = 0.05,
  numCovar.1 = 1, numCovar.2 = 1,
  R2.1 = 0.1, R2.2 = 0.1,
  ICC.2 = 0.2, ICC.3 = 0.2,
  omega.2 = 0.2, omega.3 = 0.2, rho = 0.5
)

## ----echo = FALSE-------------------------------------------------------------
knitr::kable(K)

## -----------------------------------------------------------------------------
p <- update(K, type = "power", tnum = 100000)

## ----echo = FALSE-------------------------------------------------------------
knitr::kable(p)

## ----echo = FALSE-------------------------------------------------------------
set.seed( 333434445 )

## -----------------------------------------------------------------------------
J1 <- pump_sample(
    d_m = "d3.1_m3rr2rr",
    typesample = "J",
    MTP = "HO",
    target.power = target_power,
    power.definition = "D1indiv",
    K = 15,
    nbar = 50,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.1,
    ICC.2 = 0.2, ICC.3 = 0.2,
    omega.2 = 0.2, omega.3 = 0.2,
    rho = 0.5
)

## ----include=FALSE------------------------------------------------------------
sp <- search_path(J1)
sp$dx[nrow(sp)] * 4

## ----echo = FALSE-------------------------------------------------------------
knitr::kable(J1)

## ----echo = TRUE--------------------------------------------------------------
search_path(J1)

## ----fig.width=7, fig.align="center"------------------------------------------
plot(J1, type = "search")

## ----echo = FALSE-------------------------------------------------------------
set.seed( 333434447 )

## -----------------------------------------------------------------------------
power_curve(J1)

## ----fig.width=5, fig.align="center"------------------------------------------
plot(J1)

## ----echo = FALSE-------------------------------------------------------------
set.seed(2344)

## -----------------------------------------------------------------------------
pp1 <- pump_power(
    d_m = "d3.3_m3rc2rc",
    MTP = "HO",
    nbar = 50,
    K = 20,
    J = 40,
    M = 3,
    MDES = rep(0.25, 3),
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1, numCovar.3 = 1,
    R2.1 = 0.1, R2.2 = 0.1, R2.3 = 0.1,
    ICC.2 = 0.1, ICC.3 = 0.1,
    omega.2 = 0, omega.3 = 0, rho = 0.5
)

## ----echo = FALSE-------------------------------------------------------------
knitr::kable(pp1)

## ----warning = TRUE-----------------------------------------------------------
nbar1 <- pump_sample(
        d_m = "d3.3_m3rc2rc",
        power.definition = "D1indiv",
        target.power = 0.2594,
        typesample = "nbar",
        MTP = "HO",
        K = 20,
        J = 40,
        M = 3,
        MDES = rep(0.25, 3),
        Tbar = 0.5, alpha = 0.05,
        numCovar.1 = 1, numCovar.2 = 1, numCovar.3 = 1,
        R2.1 = 0.1, R2.2 = 0.1, R2.3 = 0.1,
        ICC.2 = 0.1, ICC.3 = 0.1,
        omega.2 = 0, omega.3 = 0, rho = 0.5
)

## ----echo = FALSE-------------------------------------------------------------
knitr::kable(nbar1)

## ----fig.width=5, fig.align="center"------------------------------------------
plot( nbar1 )

## ----echo = FALSE-------------------------------------------------------------
set.seed(2344)

## -----------------------------------------------------------------------------
pp2 <- pump_power(
    d_m = "d3.3_m3rc2rc",
    MTP = "HO",
    nbar = 10,
    K = 20,
    J = 40,
    M = 3,
    MDES = rep(0.25, 3),
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1, numCovar.3 = 1,
    R2.1 = 0.1, R2.2 = 0.1, R2.3 = 0.1,
    ICC.2 = 0.1, ICC.3 = 0.1,
    omega.2 = 0, omega.3 = 0, rho = 0.5
)

## ----echo = FALSE-------------------------------------------------------------
knitr::kable(pp2)

## ----echo = FALSE-------------------------------------------------------------
set.seed(2344)

## -----------------------------------------------------------------------------
nbar2 <- pump_sample(
    d_m = "d3.3_m3rc2rc",
    typesample = "nbar",
    MTP = "HO",
    target.power = pp1$D1indiv[2],
    power.definition = "D1indiv",
    K = 20,
    J = 40,
    M = 3,
    MDES = rep(0.25, 3),
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1, numCovar.3 = 1,
    R2.1 = 0.1, R2.2 = 0.1, R2.3 = 0.1,
    ICC.2 = 0.1, ICC.3 = 0.1,
    omega.2 = 0, omega.3 = 0, rho = 0.5,
    max_sample_size_nbar = 100
)

## ----echo = FALSE-------------------------------------------------------------
knitr::kable(nbar2)

## -----------------------------------------------------------------------------
nbar3 <- pump_sample(
        d_m = "d3.3_m3rc2rc",
        power.definition = "D1indiv",
        target.power = 0.4,
        typesample = "nbar",
        MTP = "HO",
        K = 20,
        J = 40,
        M = 3,
        MDES = rep(0.25, 3),
        Tbar = 0.5, alpha = 0.05,
        numCovar.1 = 1, numCovar.2 = 1, numCovar.3 = 1,
        R2.1 = 0.1, R2.2 = 0.1, R2.3 = 0.1,
        ICC.2 = 0.1, ICC.3 = 0.1,
        omega.2 = 0, omega.3 = 0, rho = 0.5
)

## ----echo = FALSE-------------------------------------------------------------
knitr::kable(nbar3)

## ----fig.width=7, fig.align="center"------------------------------------------
plot(nbar3, type = "search")

