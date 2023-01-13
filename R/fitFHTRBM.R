##
## R package reflbrown by Author and Author
## Copyright (C) 2022
##
## This file is part of the R package reflbrown.
##
## The R package reflbrown is free software: You can redistribute it and/or
## modify it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or any later
## version (at your option). See the GNU General Public License at
## <https://www.gnu.org/licenses/> for details.
##
## The R package reflbrown is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##


#' @importFrom nimble nimbleFunction nimbleCode nimbleMCMC rexp_nimble
#' @importFrom nimble inprod dinvgamma getNimbleOption
#' @importFrom stats rnorm dnorm model.matrix
#' @importFrom reda Recur
#' @importFrom Formula Formula
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


runMCMC_FHTRBM <- function(frailty,
                           obj_df,
                           XX_sigma_df,
                           XX_kappa_df,
                           fit_init, x0, nu,
                           thin, nburnin, niter, nchains){
  # 1: global environment
  pos <- 1
  envir = as.environment(pos)
  assign('pFHT_nimble', pFHT_nimble, envir = envir)
  assign('rFHT_lklh_nimble', rFHT_lklh_nimble, envir = envir)
  assign('dFHT_lklh_nimble', dFHT_lklh_nimble, envir = envir)


  # alpha <- beta <- gamma <- theta <- z1 <- z2 <- NULL
  # x <- evt <- x0 <- nu <- NULL
  z1 <- z2 <- beta <- alpha <- NULL
  n_evts <- n_subj <- n_para_sigma <- n_para_kappa <- NULL

  coef_sigma_name <- colnames(XX_sigma_df)
  coef_kappa_name <- colnames(XX_kappa_df)
  obj_sigma_df <- cbind(obj_df, XX_sigma_df)
  obj_kappa_df <- cbind(obj_df, XX_kappa_df)

  # set up "data" for nimble
  fit_data <- list(x = obj_df$gapt, evt = obj_df$event)

  # set up "constant" for nimble
  # get covariate matrix
  X_sigma <- obj_sigma_df[!duplicated(obj_sigma_df$id), c(4:ncol(obj_sigma_df))]
  X_kappa <- obj_kappa_df[!duplicated(obj_kappa_df$id), c(4:ncol(obj_kappa_df))]
  # print(head(X_sigma))
  # print(head(X_kappa))
  fit_const <- list(n_evts = nrow(obj_df),
                    n_subj = length(unique(obj_df$id)),
                    id = obj_df$id,
                    n_para_sigma = ncol(X_sigma),
                    n_para_kappa = ncol(X_kappa),
                    X_sigma = X_sigma,
                    X_kappa = X_kappa,
                    x0 = x0,
                    nu = nu)


  if (frailty == "correlated"){
    if (is.null(fit_init)){
      fit_init <- list(alpha = rep(0, fit_const$n_para_kappa),
                       beta = rep(0, fit_const$n_para_sigma),
                       gamma = 0,
                       theta = c(0.2, 0.2),
                       z1 = rnorm(fit_const$n_subj, 0, 0.1),
                       z2 = rnorm(fit_const$n_subj, 0, 0.1))
    }

    fit_code <- nimble::nimbleCode({
      for(i in 1:n_evts) {
        x[i] ~ dFHT_lklh_nimble(event = evt[i], x0 = x0,
                                nu = nu, kappa = kappaa[id[i]],
                                sigma = sigmaa[id[i]])
      }
      for (i in 1:n_subj){
        z1[i] ~ dnorm(0, sd = sqrt(theta[1]))
        z2[i] ~ dnorm(0, sd = sqrt(theta[2]))
        kappaa[i] <- x0 + exp(alphaa[i] + gamma * z1[i] + z2[i])
        sigmaa[i] <- exp(betaa[i] + z1[i])
        alphaa[i] <- inprod(X_kappa[i, ], alpha[1:n_para_kappa])
        betaa[i] <- inprod(X_sigma[i, ], beta[1:n_para_sigma])
      }
      for (i in 1:n_para_sigma){
        beta[i] ~ dnorm(0, sd = 10)
      }
      for (i in 1:n_para_kappa){
        alpha[i] ~ dnorm(0, sd = 10)
      }
      for(i in 1:2){
        theta[i] ~ dinvgamma(1, 1)
      }
      gamma ~ dnorm(0, sd = 10)
    })

  } else if (frailty == "independent") {
    if (is.null(fit_init)){
      fit_init <- list(alpha = rep(0, fit_const$n_para_kappa),
                       beta = rep(0, fit_const$n_para_sigma),
                       theta = c(0.2, 0.2),
                       z1 = rnorm(fit_const$n_subj, 0, 0.1),
                       z2 = rnorm(fit_const$n_subj, 0, 0.1))
    }

    fit_code <- nimble::nimbleCode({
      for(i in 1:n_evts) {
        x[i] ~ dFHT_lklh_nimble(event = evt[i], x0 = x0,
                                nu = nu, kappa = kappaa[id[i]],
                                sigma = sigmaa[id[i]])
      }
      for (i in 1:n_subj){
        z1[i] ~ dnorm(0, sd = sqrt(theta[1]))
        z2[i] ~ dnorm(0, sd = sqrt(theta[2]))
        kappaa[i] <- x0 + exp(alphaa[i]  + z2[i])
        sigmaa[i] <- exp(betaa[i] + z1[i])
        alphaa[i] <- inprod(X_kappa[i, ], alpha[1:n_para_kappa])
        betaa[i] <- inprod(X_sigma[i, ], beta[1:n_para_sigma])
      }
      for (i in 1:n_para_sigma){
        beta[i] ~ dnorm(0, sd = 10)
      }
      for (i in 1:n_para_kappa){
        alpha[i] ~ dnorm(0, sd = 10)
      }
      for(i in 1:2){
        theta[i] ~ dinvgamma(1, 1)
      }
    })
  } else if (frailty == "shared"){
    if (is.null(fit_init)){
      fit_init <- list(alpha = rep(0, fit_const$n_para_kappa),
                       beta = rep(0, fit_const$n_para_kappa),
                       gamma = 0,
                       theta = 0.2,
                       z1 = rnorm(fit_const$n_subj, 0, 0.1))
    }

    fit_code <- nimble::nimbleCode({
      for(i in 1:n_evts) {
        x[i] ~ dFHT_lklh_nimble(event = evt[i], x0 = x0,
                                nu = nu, kappa = kappaa[id[i]],
                                sigma = sigmaa[id[i]])
      }
      for (i in 1:n_subj){
        z1[i] ~ dnorm(0, sd = sqrt(theta))
        # z2[i] ~ dnorm(0, sd = sqrt(theta[2]))
        kappaa[i] <- x0 + exp(alphaa[i]  + gamma * z1[i])
        sigmaa[i] <- exp(betaa[i] + z1[i])
        alphaa[i] <- inprod(X_kappa[i, ], alpha[1:n_para_kappa])
        betaa[i] <- inprod(X_sigma[i, ], beta[1:n_para_sigma])
      }
      for (i in 1:n_para_sigma){
        beta[i] ~ dnorm(0, sd = 10)
      }
      for (i in 1:n_para_kappa){
        alpha[i] ~ dnorm(0, sd = 10)
      }
      theta ~ dinvgamma(1, 1)
      gamma ~ dnorm(0, sd = 10)
    })
  } else {
    # none frailty
    if (is.null(fit_init)){
      fit_init <- list(alpha = rep(0, fit_const$n_para_kappa),
                       beta = rep(0, fit_const$n_para_kappa))
    }

    fit_code <- nimble::nimbleCode({
      for(i in 1:n_evts) {
        x[i] ~ dFHT_lklh_nimble(event = evt[i], x0 = x0, nu = nu,
                                kappa = kappaa[id[i]], sigma = sigmaa[id[i]])
      }
      for (i in 1:n_subj){
        # z1[i] ~ dnorm(0, sd = sqrt(theta))
        # z2[i] ~ dnorm(0, sd = sqrt(theta[2]))
        kappaa[i] <- x0 + exp(alphaa[i])
        sigmaa[i] <- exp(betaa[i])
        alphaa[i] <- inprod(X_kappa[i, ], alpha[1:n_para_kappa])
        betaa[i] <- inprod(X_sigma[i, ], beta[1:n_para_sigma])
      }
      for (i in 1:n_para_sigma){
        beta[i] ~ dnorm(0, sd = 10)
      }
      for (i in 1:n_para_kappa){
        alpha[i] ~ dnorm(0, sd = 10)
      }
    })
  }
  # print(names(fit_init))
  mcmc <- nimble::nimbleMCMC(data = fit_data,
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
    c(paste0(coef_kappa_name, "_Kappa"), paste0(coef_sigma_name, "_Sigma"))
  colnames(mcmc)[grepl("theta", colnames(mcmc))] <-
    c("fVarSigma", "fVarKappa")[1:length(fit_init$theta)]
  # print(names(mcmc))
  return(mcmc)
}


#' @title Recurrent Events Threshold Regression
#'
#' @description The function fitFhtrbm return MCMC samples.
#'
#' @param formula   A formula object, with the response on the left
#'     of a "~" operator, and the predictors on the right. The response must be
#'     a recurrent event object as returned by function Recur.
#' @param data      A data frame includes all the variables in "formula".
#'     Data frame should not include intercept.
#' @param frailty   Frailty model type, \code{"correlated"}, \code{"independent"},
#'     \code{"shared"}, and \code{"none"}.
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
#' @examples
#' data(simuEvtDat)
#' # need to run longer MCMC chain
#' # note that simuEvtDat dataset does not include intercept column
#' \dontrun{
#' mcmc <- fitFhtrbm(Recur(time = time, id = id, event = event) ~ x1 + x2|x1 + x2,
#'                  data = simuEvtDat,
#'                  frailty = "correlated",
#'                  thin = 10, nburnin = 3000, niter = 10000, nchains = 1)
#' }
#' @export
fitFhtrbm <- function(formula, data, x0 = 10, nu = 3.9,
                      frailty = c("correlated", "independent", "shared", "none"),
                      initial = NULL,
                      thin = 10, nburnin = 3000, niter = 6000, nchains = 1){
  frailty <- match.arg(frailty)
  cl <- match.call(expand.dots = FALSE)
  indx <- match(c("formula", "data"),
                names(cl), nomatch=0)
  mf<- cl[c(1, indx)]
  f <- Formula(formula)
  f1<-formula(f,lhs=1)
  f1<-Formula(f1)
  mf[[1]] <- as.name("model.frame")

  mf$formula <- f1
  mf <- eval(mf, parent.frame())
  obj <- stats::model.extract(mf, "response")
  f2<-formula(f1, lhs = 0)

  formula_sigma <- formula(f1, lhs = 0, rhs = 1)
  formula_kappa <- formula(f1, lhs = 0, rhs = 2)
  XX_sigma <- as.data.frame(model.matrix(formula_sigma, data))
  XX_kappa <- as.data.frame(model.matrix(formula_kappa, data))
  coef_sigma_name <- colnames(XX_sigma)
  coef_kappa_name <- colnames(XX_kappa)

  obj <- as.data.frame(obj)
  obj$gapt <- obj$time2 - obj$time1
  obj <- obj[, c("gapt", "id", "event")]

  # set up nimble model code
  mcmc <- runMCMC_FHTRBM(frailty = frailty,
                         obj_df = obj,
                         XX_sigma_df = XX_sigma,
                         XX_kappa_df = XX_kappa,
                         fit_init = initial,
                         x0 = x0,
                         nu = nu,
                        thin = thin,
                        nburnin = nburnin,
                        niter = niter,
                        nchains = nchains)
  return(mcmc)
}





