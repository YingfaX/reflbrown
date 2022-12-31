##
## R package refbrown by Author and Author
## Copyright (C) 2022
##
## This file is part of the R package refbrown.
##
## The R package refbrown is free software: You can redistribute it and/or
## modify it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or any later
## version (at your option). See the GNU General Public License at
## <https://www.gnu.org/licenses/> for details.
##
## The R package refbrown is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##


#' @importFrom stats rnorm dnorm model.matrix
#' @importFrom reda Recur
#' @importFrom Formula Formula
NULL


#' @title Model comparison criteria
#'
#' @description The function getModelFit calculates the
#' Deviance Information Criterion (DIC) and
#' the Logarithm of the Pseudo-Marginal Likelihood (LPML).
#' The DIC and LPML are calculated based on observed likelihood.
#' Monte Carlo integration is used for approximating the observed likelihood.
#'
#' @param formula   A formula object, with the response on the left
#'     of a "~" operator, and the predictors on the right. The response must be
#'     a recurrent event object as returned by function Recur.
#' @param data      A data frame includes all the variables in "formula".
#'     Data frame should not include intercept.
#' @param frailty   Frailty model type, \code{"correlated"}, \code{"independent"},
#'     \code{"shared"} and \code{"none"}.
#' @param mcmc_sigma MCMC samples for coefficients in sigma model,
#'     which is a matrix.
#' @param mcmc_kappa MCMC samples for coefficients in kappa model,
#'     which is a matrix.
#' @param mcmc_gamma MCMC samples for coefficients of frailty,
#'     which is a matrix. If using the model with this parameter,
#'     user should provide MCMC samples. Default is NULL.
#' @param mcmc_theta MCMC samples for variance parameter of frailty,
#'     which is a matrix. If using the model with this parameter,
#'     user should provide MCMC samples. Default is NULL.
#' @param x0        The starting point of Brownian motion, which is a known
#'     parameter in model.
#' @param nu        The lower boundary of Brownian motion, which is a known
#'     parameter in model.
#' @param mcsize    Monte Carlo approximation sample size. Default is 500.
#' @param setSeed   Logical or numeric argument. Default is TRUE. The random
#'     seed for Monte Carlo samples.
#' @examples
#' \dontrun{
#' # using the MCMC from fitFhtrbm function
#' ret <- getModelFit(Recur(time = time, id = id, event = event) ~ x1 + x2|x1 + x2,
#'                    data = simuEvtDat,
#'                    frailty = "correlated",
#'                    mcmc_sigma = mcmc[, c(4:6)],
#'                    mcmc_kappa = mcmc[, c(1:3)],
#'                    mcmc_gamma = as.matrix(mcmc[, 7]), # need matrix
#'                    mcmc_theta = mcmc[, c(8:9)])
#' }
#'
#' @export
getModelFit <- function(formula,
                        data,
                        frailty = c("correlated", "independent", "shared", "none"),
                        mcmc_sigma,
                        mcmc_kappa,
                        mcmc_gamma = NULL,
                        mcmc_theta = NULL,
                        x0 = 10,
                        nu = 3.9,
                        mcsize = 500,
                        setSeed = TRUE){
  if (setSeed){
    if (is.numeric(setSeed)){
      set.seed(setSeed)
    } else {
      set.seed(1)
    }
  }

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
  n_para_sigma <- ncol(XX_sigma) # include intercept
  n_para_kappa <- ncol(XX_kappa) # include intercept

  obj <- as.data.frame(obj)
  obj$gapt <- obj$time2 - obj$time1
  obj <- obj[, c("gapt", "id", "event")]

  obj_comb <- cbind(obj, XX_sigma, XX_kappa)

  obj_comb$id <- factor(obj_comb$id)
  obj_comb_ls <- split(obj_comb, f = obj_comb$id)

  # get standard MC samples for MC approximation
  mc_z <- matrix(c(rnorm(mcsize), rnorm(mcsize)),
                 nrow = mcsize, ncol = 2, byrow = F)

  mcmc_niter <- nrow(mcmc_sigma)
  n_subj <- length(unique(obj$id))


  if (is.null(mcmc_gamma)) {
    mcmc_gamma <- matrix(0, 1, 1)
  }
  if (is.null(mcmc_theta)) {
    mcmc_theta <- matrix(0, 1, 1)
  }

  mcmc_post_mean_sigma <- matrix(as.numeric(apply(mcmc_sigma, 2, mean)),
                                    ncol = ncol(mcmc_sigma), byrow = TRUE)
  mcmc_post_mean_kappa <- matrix(as.numeric(apply(mcmc_kappa, 2, mean)),
                                 ncol = ncol(mcmc_kappa), byrow = TRUE)
  mcmc_post_mean_gamma <- matrix(as.numeric(apply(mcmc_gamma, 2, mean)),
                                 ncol = ncol(mcmc_gamma), byrow = TRUE)
  mcmc_post_mean_theta <- matrix(as.numeric(apply(mcmc_theta, 2, mean)),
                                 ncol = ncol(mcmc_theta), byrow = TRUE)

  # last row for saving observed likelihood for posterior mean of para
  loglklh_matrix_subj <- matrix(NA, mcmc_niter+1, n_subj)

  for (i in 1:n_subj){

    Sys.sleep(0.01)
    cat("\rFinished", round(i/n_subj*100), "%")

    gapt_subj <- as.numeric(obj_comb_ls[[i]]$gapt)
    event_subj <- as.numeric(obj_comb_ls[[i]]$event)
    covar_sigma <- as.numeric(obj_comb_ls[[i]][1, c(4:(n_para_sigma+3))])
    covar_kappa <- as.numeric(obj_comb_ls[[i]][1, c((n_para_sigma+4):(n_para_sigma+n_para_kappa+3))])

    loglklh_matrix_subj[c(1:mcmc_niter), i] <- log(obsLklh_c(gapt_subj,
                                                             event_subj,
                                                             covar_sigma,
                                                             covar_kappa,
                                                             mcmc_sigma,
                                                             mcmc_kappa,
                                                             mcmc_gamma,
                                                             mcmc_theta,
                                                             mc_z[, 1],
                                                             mc_z[, 2],
                                                             frailty, x0, nu))
    loglklh_matrix_subj[(1+mcmc_niter), i] <- log(obsLklh_c(gapt_subj,
                                                            event_subj,
                                                            covar_sigma,
                                                            covar_kappa,
                                                            mcmc_post_mean_sigma,
                                                            mcmc_post_mean_kappa,
                                                            mcmc_post_mean_gamma,
                                                            mcmc_post_mean_theta,
                                                            mc_z[, 1],
                                                            mc_z[, 2],
                                                            frailty, x0, nu))
  }
  # DIC
  ## mean of deviance
  mean_of_dev <- (-2) * mean(apply(loglklh_matrix_subj[c(1:mcmc_niter), ], 1, sum))
  ## deviance of (posterior) mean
  dev_of_mean <- (-2) * sum(loglklh_matrix_subj[(mcmc_niter+1), ])
  pd_p <- mean_of_dev - dev_of_mean

  # LPML
  LPML <- sum(log(1/apply(1/exp(loglklh_matrix_subj[c(1:mcmc_niter), ]), 2, mean)))

  ## result
  ret <- list(DIC = (dev_of_mean+2*pd_p),
              DIC_mean = dev_of_mean,
              pd = pd_p,
              LPML = LPML)
  return(ret)
}


