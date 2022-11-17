##
## R package reThReg by Yingfa Xie and Jun Yan
## Copyright (C) 2022
##
## This file is part of the R package reThReg.
##
## The R package reThReg is free software: You can redistribute it and/or
## modify it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or any later
## version (at your option). See the GNU General Public License at
## <https://www.gnu.org/licenses/> for details.
##
## The R package reThReg is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##

#' @useDynLib reThReg, .registration = TRUE
#' @importFrom nimble nimbleFunction
#' @importFrom nimble nimbleCode
#' @importFrom nimble nimbleMCMC
#' @importFrom nimble rexp_nimble
#' @importFrom nimble inprod
#' @importFrom nimble dinvgamma
#' @importFrom stats dnorm
#' @importFrom reda Recur
NULL

## cdf of FHT distribution of reflected Brownian motion
pFHT_nimble <- nimble::nimbleFunction(
  run = function(x = double(0),
                 x0 = double(0), nu = double(0), kappa = double(0), 
                 sigma = double(0), log = integer(0, default = 0)){
    returnType(double(0))
    if (x <= 0){
      if(log) return(log(0))
      else return(0)
    }

    n <- 1; smt <- 0; temp <- 1
    ## cdf
    while ( abs(temp) > 0) {
      lambda_n <- (2 * n - 1)^2 * sigma^2 * pi^2 / (8 * (kappa - nu)^2)
      c_n <- (-1)^(n + 1) * 4 * cos((2 * n - 1) * pi * (kappa - x0) / 
                                      (2 * (kappa - nu))) /((2 * n - 1) * pi)
      temp <- c_n * exp(-lambda_n * x)
      smt <- smt + temp
      n <- n + 1
    }
    cprob <- 1 - smt

    if (cprob < 0) cprob = 0
    if (cprob > 1) cprob = 1
    logProb = log(cprob)
    if(log) return(logProb)
    else return(exp(logProb))
  })


dFHT_lklh_nimble <- nimble::nimbleFunction(
  run = function(x = double(0), event = double(0),
                 x0 = double(0), nu = double(0), kappa = double(0), 
                 sigma = double(0), log = integer(0, default = 0)){
    returnType(double(0))

    unit <- 1
    x_intv_end <- x + unit/2
    x_intv_start <- x - unit/2

    if (event == 1){
      ## interval censored
      lklh <- pFHT_nimble(x_intv_end, x0, nu, kappa, sigma) -
        pFHT_nimble(x_intv_start, x0, nu, kappa, sigma)
    } else {
      ## survival probability
      lklh <- 1 - pFHT_nimble(x_intv_end, x0, nu, kappa, sigma)
    }


    logProb = log(lklh)
    if(log) return(logProb)
    else return(exp(logProb))
  })


rFHT_lklh_nimble <- nimble::nimbleFunction(
  run = function(n = double(0), event = double(0),
                 x0 = double(0), nu = double(0),
                 kappa = double(0), sigma = double(0)){
    returnType(double(0))
    # if(n != 1) print("rmyexp only allows n = 1; using n = 1.")
    cand <- rexp_nimble(1, 0.04) + 0.3
    # cand <- rexp_nimble(1, 0.04) + 3
    return(cand)
  })


runMCMC_reThReg <- function(modeltype, fit_data, fit_const, fit_init, 
                            coefName, thin, nburnin, niter, nchains){
  if (modeltype == "1"){
    if (is.null(fit_init)){
      fit_init <- list(alpha = rep(0, fit_const$n_para),
                       beta = rep(0, fit_const$n_para),
                       gamma = 0,
                       theta = c(0.2, 0.2),
                       z1 = rnorm(fit_const$n_subj, 0, 0.1),
                       z2 = rnorm(fit_const$n_subj, 0, 0.1))
    }

    fit_code <- nimbleCode({
      for(i in 1:fit_const$n_evts) {
        x[i] ~ dFHT_lklh_nimble(event = evt[i], x0 = fit_const$x0, 
                                nu = nu, kappa = kappaa[id[i]], 
                                sigma = sigmaa[id[i]])
      }
      for (i in 1:fit_const$n_subj){
        z1[i] ~ dnorm(0, sd = sqrt(theta[1]))
        z2[i] ~ dnorm(0, sd = sqrt(theta[2]))
        kappaa[i] <- fit_const$x0 + exp(alphaa[i] + gamma * z1[i] + z2[i])
        sigmaa[i] <- exp(betaa[i] + z1[i])
        alphaa[i] <- inprod(X[i, ], alpha[1:fit_const$n_para])
        betaa[i] <- inprod(X[i, ], beta[1:fit_const$n_para])
      }
      for (i in 1:fit_const$n_para){
        beta[i] ~ dnorm(0, sd = 10)
      }
      for (i in 1:fit_const$n_para){
        alpha[i] ~ dnorm(0, sd = 10)
      }
      for(i in 1:2){
        theta[i] ~ dinvgamma(1, 1)
      }
      gamma ~ dnorm(0, sd = 10)
    })

  } else if (modeltype == "2") {
    if (is.null(fit_init)){
      fit_init <- list(alpha = rep(0, fit_const$n_para),
                       beta = rep(0, fit_const$n_para),
                       theta = c(0.2, 0.2),
                       z1 = rnorm(fit_const$n_subj, 0, 0.1),
                       z2 = rnorm(fit_const$n_subj, 0, 0.1))
    }

    fit_code <- nimbleCode({
      for(i in 1:fit_const$n_evts) {
        x[i] ~ dFHT_lklh_nimble(event = evt[i], x0 = fit_const$x0, 
                                nu = nu, kappa = kappaa[id[i]], 
                                sigma = sigmaa[id[i]])
      }
      for (i in 1:fit_const$n_subj){
        z1[i] ~ dnorm(0, sd = sqrt(theta[1]))
        z2[i] ~ dnorm(0, sd = sqrt(theta[2]))
        kappaa[i] <- fit_const$x0 + exp(alphaa[i]  + z2[i])
        sigmaa[i] <- exp(betaa[i] + z1[i])
        alphaa[i] <- inprod(X[i, ], alpha[1:fit_const$n_para])
        betaa[i] <- inprod(X[i, ], beta[1:fit_const$n_para])
      }
      for (i in 1:fit_const$n_para){
        beta[i] ~ dnorm(0, sd = 10)
      }
      for (i in 1:fit_const$n_para){
        alpha[i] ~ dnorm(0, sd = 10)
      }
      for(i in 1:2){
        theta[i] ~ dinvgamma(1, 1)
      }
    })
  } else {
    if (is.null(fit_init)){
      fit_init <- list(alpha = rep(0, fit_const$n_para),
                       beta = rep(0, fit_const$n_para),
                       gamma = 0,
                       theta = 0.2,
                       z1 = rnorm(fit_const$n_subj, 0, 0.1))
    }

    fit_code <- nimbleCode({
      for(i in 1:fit_const$n_evts) {
        x[i] ~ dFHT_lklh_nimble(event = evt[i], x0 = fit_const$x0, 
                                nu = nu, kappa = kappaa[id[i]], 
                                sigma = sigmaa[id[i]])
      }
      for (i in 1:fit_const$n_subj){
        z1[i] ~ dnorm(0, sd = sqrt(theta))
        # z2[i] ~ dnorm(0, sd = sqrt(theta[2]))
        kappaa[i] <- fit_const$x0 + exp(alphaa[i]  + gamma * z1[i])
        sigmaa[i] <- exp(betaa[i] + z1[i])
        alphaa[i] <- inprod(X[i, ], alpha[1:fit_const$n_para])
        betaa[i] <- inprod(X[i, ], beta[1:fit_const$n_para])
      }
      for (i in 1:fit_const$n_para){
        beta[i] ~ dnorm(0, sd = 10)
      }
      for (i in 1:fit_const$n_para){
        alpha[i] ~ dnorm(0, sd = 10)
      }
      theta ~ dinvgamma(1, 1)
      gamma ~ dnorm(0, sd = 10)
    })
  }
  print(names(fit_init))
  mcmc <- nimbleMCMC(data = fit_data,
                     constants = fit_const,
                     inits = fit_init,
                     code = fit_code,
                     monitors = names(fit_init),
                     thin = thin,
                     niter = niter,
                     nburnin = nburnin,
                     nchains = nchains,
                     setSeed = TRUE)
  colnames(mcmc)[grepl("alpha|beta", colnames(mcmc))] <- 
    c(paste0(coefName, "_Kappa"), paste0(coefName, "_Sigma"))
  colnames(mcmc)[grepl("theta", colnames(mcmc))] <- 
    c("fVarSigma", "fVarKappa")[1:length(fit_init$theta)]
  # print(names(mcmc))
  return(mcmc)
}




#' @title Recurrent Events Threshold Regression
#'
#' @description The function reThReg return mcmc samples.
#'
#' @param formula   A formula object, with the response on the left
#'     of a "~" operator, and the predictors on the right. The response must be
#'     a recurrent event object as returned by function Recur.
#' @param data      A data frame includes all the variables in "formula".
#'     Data frame should not include intercept
#' @param model     Model type indicator. Three types: "1", "2", and "3"
#' @param x0        The starting point of Brownian motion, which is a known
#'     parameter in model.
#' @param nu        The lower boundary of Brownian motion, which is a known
#'     parameter in model.
#' @param initial   The initial value for MCMC, which is optional.
#'     Default is Null.
#' @param thin      Thinning number of MCMC chain.
#' @param nburnin   Burn-in number of MCMC chain.
#' @param niter     Total number of iteration of MCMC chain before thin
#'     and burn-in.
#' @param nchains   Number of chains
#' @export
reThReg <- function(formula, data, model = "1", x0 = 10, nu = 3.9,
                    initial = NULL,
                    thin = 10, nburnin = 1000, niter = 3000, nchains = 1){
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$data <- data
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mm <- stats::model.matrix(formula, data = mf)
  coefName <- colnames(mm)
  obj <- stats::model.extract(mf, "response")
  obj <- as.data.frame(obj)
  obj$gapt <- obj$time2 - obj$time1
  obj <- obj[, c("gapt", "id", "event")]
  DF <- cbind(obj, as.data.frame(mm))
  # DF <- DF[,colnames(DF) != "(Intercept)"]

  # set up "data" for nimble
  fit_data <- list(x = DF$gapt, evt = DF$event)

  # set up "constant" for nimble
  X <- DF[!duplicated(DF$id), c(4:ncol(DF))] # get covariate matrix
  fit_const = list(n_evts = nrow(DF),
                   n_subj = length(unique(DF$id)),
                   n_para = ncol(X),
                   id = DF$id,
                   X = X,
                   x0 = x0,
                   nu = nu)

  # set up nimble model code
  mcmc <- runMCMC_reThReg(modeltype = model,
                          fit_data = fit_data,
                          fit_const = fit_const,
                          fit_init = initial,
                          coefName = coefName,
                          thin = thin,
                          nburnin = nburnin,
                          niter = niter,
                          nchains = nchains)
  return(mcmc)
}





