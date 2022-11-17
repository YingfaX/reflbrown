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
#' @importFrom Rcpp sourceCpp
#' @importFrom copula safeUroot
NULL

#' @title The First Hitting Time Distribution
#'
#' @name FHTRBM
#' @description    Density, distribution, quantile, and random number
#'     generation function
#' @param x        random values; vector
#' @param x0       starting point; scalar
#' @param nu       lower barrier; scalar
#' @param kappa    upper reflection barrier; scalar
#' @param sigma    volatility; scalar
#' @param log      default is FALSE
#' @examples
#' dFHTRBM(1, 10, 3.9, 20, 3)
#' dFHTRBM(c(1:5), 10, 3.9, 20, 3)
#'
#' @export
dFHTRBM <- function(x, x0, nu, kappa, sigma, log = FALSE){
  denst <- dFHT_c(t = x, x0, nu, kappa, sigma)
  if (log) denst <- log(denst)
  return(denst)
}

#' @rdname FHTRBM
#'
#' @param q        random values; vector
#' @export
#' @examples
#' pFHTRBM(1, 10, 3.9, 20, 3)
#' pFHTRBM(c(1:5), 10, 3.9, 20, 3)
pFHTRBM <- function(q, x0, nu, kappa, sigma){
  return(pFHT_c(t = q, x0, nu, kappa, sigma))
}


qFHT_scale <- function(p, x0, nu, kappa, sigma){
  obj_quantile <- function(t){
    return(pFHTRBM(t, x0, nu, kappa, sigma) - p)
  }
  t_q <- safeUroot(obj_quantile, c(0, .Machine$integer.max))$root
  return(t_q)
}

#' @rdname FHTRBM
#'
#' @param p        vector of probabilities
#' @export
qFHTRBM <- Vectorize(FUN = qFHT_scale)


#' @rdname FHTRBM
#'
#' @param n       number of observations
#' @param k       number of pieces to be included in proposal, default k = 3
#' @param qTail   the user-defined qth quantile for the tail, default q = 0.95
#' @export
rFHTRBM <- function(n, x0, nu, kappa, sigma, k = 3, qTail = 0.95){
  lower_bound <- 0

  getM <- function(x0, nu, kappa, sigma, k, t_point, ft_point, slopes, 
                   intercepts, lambda_1, upper_bound){
    M <- rep(NA, k)
    for(i in 1:k){
      if (i < k) {
        obj_tri <- function(t){
          return(- dFHT_c(t, x0, nu, kappa, sigma) / 
                   (slopes[i] * t + intercepts[i]))
        }
        M[i] <- -optimize(obj_tri, c(t_point[i], t_point[i+1]))$objective
      } else {
        obj_exp <- function(t){
          return(- dFHT_c(t, x0, nu, kappa, sigma) / 
                   (dexp(t, rate = lambda_1) ))
        }
        M[i] <- -optimize(obj_exp, c(t_point[k], upper_bound))$objective
      }
    }
    return(M)
  }

  # get lambda_1 for tail proposal
  lambda_1 <- sigma^2 * pi^2 / (8 * (kappa - nu)^2)

  # get qth quantile point t_q
  t_q <- qFHTRBM(qTail, x0, nu, kappa, sigma)
  ft_qq <- dFHT_c(t_q, x0, nu, kappa, sigma)
  upper_bound <- qFHTRBM(0.9999, x0, nu, kappa, sigma)

  # get highest density point t_mode
  obj_FHT_mode <- function(t){
    return(-dFHT_c(t, x0, nu, kappa, sigma))
  }
  t_mode <- optimize(obj_FHT_mode, c(lower_bound, t_q))$minimum

  # start point for each piece of proposal
  # length k, with 0
  t_point <- c(lower_bound, seq(t_mode, t_q, length.out = (k-1)))
  # the density function of each t point, which will be used to
  # determine the proposal
  ft_point <- dFHT_c(t_point, x0, nu, kappa, sigma)
  # the probability function of each t point, which will be used to
  # decide sampling from which part
  pt_point <- c(pFHT_c(t_point, x0, nu, kappa, sigma)[-1], 1)

  # determine the line (slope and intercept) of each piece of proposal
  slopes <- (ft_point[2:k] - ft_point[1:(k-1)]) /
    (t_point[2:k] - t_point[1:(k-1)]) # length (k-1)
  intercepts <- ft_point[2:k] - slopes * t_point[2:k]  # length (k-1)

  # calculate the integral for each piece, 
  # the integral will be used for sampling
  I <- slopes * ( t_point[2:k]^2 - t_point[1:(k-1)]^2 ) / 2 +
    intercepts * (t_point[2:k] - t_point[1:(k-1)])
  I <- c(I, 1 - pexp(t_q, rate = lambda_1))

  # determine M of each piece
  # length k
  M <- getM(x0, nu, kappa, sigma, k, t_point, ft_point, slopes, intercepts,
            lambda_1, upper_bound)

  # sampling part implemented in C
  ret <- rFHT_c(n, x0, nu, kappa, sigma, k, lambda_1,
                t_point, slopes, intercepts,
                M, I, pt_point)$samples
  return(ret)
}


