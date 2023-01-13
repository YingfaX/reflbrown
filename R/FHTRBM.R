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

#' @useDynLib reflbrown, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom copula safeUroot
#' @importFrom stats dexp pexp optimize
NULL

#' @title The First Hitting Time Distribution
#'
#' @name FHTRBM
#' @description    Density, distribution, quantile, and random number
#'     generation function of first hitting time (hitting down) distribution
#'     of reflected Brownian motion.
#' @references Hu, Q., Wang, Y., and Yang, X. (2012). The hitting time density
#' for a reflected Brownian motion.
#' \emph{Computational Economics}, \bold{40}, 1–18.

#' @param x,q        vector of quantiles.
#' @param x0         scalar; starting point of the reflected Brownian motion.
#' @param nu         scalar; lower boundary of the reflected Brownian motion.
#' @param kappa      scalar; upper reflection barrier of the reflected Brownian
#'                     motion.
#' @param sigma      scalar; volatility of the reflected Brownian motion.
#' @param log,log.p  logical; default is FALSE; if TRUE, probabilities p
#'                     are given as log(p).
#' @examples
#' dfhtrbm(1, 10, 3.9, 20, 3)
#' dfhtrbm(c(1:5), 10, 3.9, 20, 3)
#' @export
dfhtrbm <- function(x, x0, nu, kappa, sigma, log = FALSE){
  if (kappa < x0)
    stop("The upper reflection barrier should not be smaller than
         starting points.")
  if (x0 < nu)
    stop("The starting point should not be smaller than lower boundary.")

  denst <- dfhtrbm_c(t = x, x0, nu, kappa, sigma)
  if (log) denst <- log(denst)
  return(denst)
}

#' @rdname FHTRBM
#'
#' @param lower.tail logical; if TRUE (default), probabilities
#'                     are \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}.
#' @export
#' @examples
#' pfhtrbm(1, 10, 3.9, 20, 3)
#' pfhtrbm(c(1:5), 10, 3.9, 20, 3)
pfhtrbm <- function(q, x0, nu, kappa, sigma, lower.tail = TRUE, log.p = FALSE){
  if (kappa < x0)
    stop("The upper reflection barrier should not be smaller than
         starting points.")
  if (x0 < nu)
    stop("The starting point should not be smaller than lower boundary.")

  prob <- pfhtrbm_c(t = q, x0, nu, kappa, sigma)
  if (!lower.tail) prob <- 1 - prob
  if (log.p) prob <- log(prob)
  return(prob)
}


qfhtrbm_scalar <- function(p, x0, nu, kappa, sigma,
                           lower.tail = TRUE, log.p = FALSE){
  if (kappa < x0)
    stop("The upper reflection barrier should not be smaller than
         starting points.")
  if (x0 < nu)
    stop("The starting point should not be smaller than lower boundary.")

  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p

  obj_quantile <- function(t){
    return(pfhtrbm(t, x0, nu, kappa, sigma) - p)
  }
  quant <- copula::safeUroot(obj_quantile, c(0, .Machine$integer.max))$root

  return(quant)
}

#' @rdname FHTRBM
#'
#' @param p        vector of probabilities
#' @export
qfhtrbm <- Vectorize(FUN = qfhtrbm_scalar)


#' @rdname FHTRBM
#'
#' @param n       number of observations
#' @export
rfhtrbm <- function(n, x0, nu, kappa, sigma){
  if (kappa < x0)
    stop("The upper reflection barrier should not be smaller than
         starting points.")
  if (x0 < nu)
    stop("The starting point should not be smaller than lower boundary.")

  # three piece; # tq be 95 percentile
  lower_bound <- 0; k = 3; qTail = 0.95


  getM <- function(x0, nu, kappa, sigma, k, t_point, ft_point, slopes,
                   intercepts, lambda_1, upper_bound){
    M <- rep(NA, k)
    for(i in 1:k){
      if (i < k) {
        obj_tri <- function(t){
          return(- dfhtrbm_c(t, x0, nu, kappa, sigma) /
                   (slopes[i] * t + intercepts[i]))
        }
        M[i] <- -optimize(obj_tri, c(t_point[i], t_point[i+1]))$objective
      } else {
        obj_exp <- function(t){
          return(- dfhtrbm_c(t, x0, nu, kappa, sigma) /
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
  t_q <- qfhtrbm(qTail, x0, nu, kappa, sigma)
  ft_qq <- dfhtrbm_c(t_q, x0, nu, kappa, sigma)
  upper_bound <- qfhtrbm(0.9999, x0, nu, kappa, sigma)

  # get highest density point t_mode
  obj_FHT_mode <- function(t){
    return(-dfhtrbm_c(t, x0, nu, kappa, sigma))
  }
  t_mode <- optimize(obj_FHT_mode, c(lower_bound, t_q))$minimum

  # start point for each piece of proposal
  # length k, with 0
  t_point <- c(lower_bound, seq(t_mode, t_q, length.out = (k-1)))
  # the density function of each t point, which will be used to
  # determine the proposal
  ft_point <- dfhtrbm_c(t_point, x0, nu, kappa, sigma)
  # the probability function of each t point, which will be used to
  # decide sampling from which part
  pt_point <- c(pfhtrbm_c(t_point, x0, nu, kappa, sigma)[-1], 1)

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
  ret <- rfhtrbm_c(n, x0, nu, kappa, sigma, k, lambda_1,
                t_point, slopes, intercepts,
                M, I, pt_point)$samples
  return(ret)
}




