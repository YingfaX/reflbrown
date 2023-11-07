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


#' @importFrom stats dexp optimize pexp rnorm
NULL

#' @title Simulated Recurrent Events
#'
#' @description The function simData generates simulated recurrent events from
#' reflected Brownian motion.
#'
#' @param size       Patient sample size. This argument should be a numeric value
#' @param endTime    The end of the follow up times of patients.
#'     This argument should be a numeric vector.
#' @param X          The covariates matrix for patients.
#'     This argument should be a matrix which records the covariate vector of
#'     subject $i$ in row i.
#' @param sigmaCoef  The coefficents of volatility.
#'     This argument should be a numeric vector.
#' @param kappaCoef  The coefficents of upper reflection barrier.
#'     This argument should be a numeric vector.
#' @param gamma      The coefficient of frailty.
#'     This argument should be a numeric value.
#' @param x0         The starting point of the reflected Brownian motion.
#'     This argument should be a numeric value.
#' @param nu         The lower boundary of the reflected Brownian motion.
#'     This argument should be a numeric value.
#' @param theta      The variance of frailties.
#'     This argument should be a numeric vector with length $2$.
#' @param sigmaZ     The frailty in sigma, default is NULL.
#'     This argument should be a numeric vector.
#' @param kappaZ     The frailty in kappa, default is NULL.
#'     This argument should be a numeric vector.
#' @param roundTime  A logical value with default value TRUE indicating whether
#'     to round the generated gap times or not round.
#' @param gapTime    A logical value with default value FALSE indicating whether
#'     to return the generated gap times or calender time.
#' @param event_num  The number of simulated random gap times for one patients.
#'     Only the cumulative random gap times smaller the follow up time are
#'     saved as output. This argument should be a numeric value.
#' @param FrailtyZ  Report frailty in the output dataframe or not. Default is False.
#'
#' @examples
#' data(simuCovDat)
#' size <- nrow(simuCovDat)
#' set.seed(2)
#' Intercept <- rep(1, size)
#' simuDat <- simData(size = size,
#'                    endTime = simuCovDat$censorTime,
#'                    X = cbind(Intercept, as.matrix(simuCovDat[, c("x1", "x2")])),
#'                    sigmaCoef = c(0.9, -0.2, -0.1),
#'                    kappaCoef = c(2.9, 0.2, -0.1),
#'                    gamma = -0.4,
#'                    x0 = 10,
#'                    nu = 3.9,
#'                    theta = c(0.2, 0.3))
#'
#' head(simuDat)
#'
#' @export
simData <- function(size,
                    endTime,
                    X,
                    sigmaCoef,
                    kappaCoef,
                    gamma, x0, nu,
                    theta = NULL,
                    sigmaZ = NULL,
                    kappaZ = NULL,
                    roundTime = TRUE,
                    gapTime = FALSE,
                    event_num = 500,
                    FrailtyZ = FALSE){
  if (is.null(sigmaZ)){
    # print("Generating frailty")
    sigmaZ <- rnorm(size, mean = 0, sd = sqrt(theta[1]))
  }
  # print(sigmaZ[1:5])
  sigmas <- exp(X %*% sigmaCoef + sigmaZ)

  if (is.null(kappaZ)){
    kappaZ <- rnorm(size, mean = 0, sd = sqrt(theta[2]))
  }
  # print(kappaZ[1:5])
  kappas <- x0 + exp(X %*% kappaCoef + gamma * sigmaZ + kappaZ)

  # generate hypo-event time and set status
  dat_evt <- data.frame('id' = NA, 'event' = NA,'time' = NA)
  for (i in 1:size){
    simGapTimes <- rfhtrbm(n = event_num, x0 = x0, nu = nu,
                           kappa = kappas[i], sigma = sigmas[i])

    # find index
    time <- cumsum(simGapTimes)
    time_end_indx <- which(cumsum(simGapTimes) >= endTime[i])[1]
    time <- time[1:time_end_indx]

    if (time_end_indx == 1) {
      time[time_end_indx] <- round(endTime[i])
    } else {
      time <- simGapTimes[1:time_end_indx]
      time[time_end_indx] <- endTime[i] - sum(time[1:(time_end_indx - 1)])

      if (!roundTime) {
        if (!gapTime) {
          time <- cumsum(time)
        }
      } else {
        time <- round(time)
        if (!gapTime) {
          time <- cumsum(time)
        }
      }
    }
    event <- c(rep(1, length(time) - 1), 0)
    dat_evt_new <- data.frame('id' = rep(i, length(time)),
                              'event' = event,
                              'time' = time)
    dat_evt <- rbind(dat_evt, dat_evt_new)
  }

  dat_evt <- dat_evt[!is.na(dat_evt$time), ]
  X <- data.frame(X)
  X$id <- c(1:size)
  dat_sim <- merge(dat_evt, X, by = 'id')

  if (FrailtyZ) {
    dat_frailty_out <- data.frame('id' = c(1:size),
                                  'Fsigma' = sigmaZ,
                                  'Fkappa' = kappaZ)
    dat_sim <- merge(dat_sim, dat_frailty_out, by = 'id', all.x = TRUE)
  }

  return(dat_sim)
}





