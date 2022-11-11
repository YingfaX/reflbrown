

#' @title Simulated Recurrent Events 
#' 
#' @description The function simData generates simulated recurrent events from 
#' reflected Brownian motion.
#' 
#' @param size         Patient sample size. This argument should be a numeric value
#' @param endTime      The end of the follow up times of patients. 
#'    This argument should be a numeric vector.
#' @param X            The covariates matrix for patients. 
#'    This argument should be a matrix which records the covariate vector of subject $i$ in row i.  
#' @param kappaCoef    The coefficents of upper reflection barrier. 
#'    This argument should be a numeric vector.
#' @param sigmaCoef    The coefficents of volatility. 
#'    This argument should be a numeric vector.
#' @param theta        The variance of frailties. 
#'    This argument should be a numeric vector with length $2$.
#' @param gamma        The coefficient of frailty. 
#'    This argument should be a numeric value. 
#' @param x0           The starting point of the reflected Brownian motion. 
#'    This argument should be a numeric value.
#' @param nu           The lower boundary of the reflected Brownian motion. 
#'    This argument should be a numeric value.
#' @param roundTime    A logical value with default value TRUE indicating whether 
#'    to round the generated gap times or not round.
#' @param gapTime      A logical value with default value FALSE indicating whether 
#'    to return the generated gap times or calender time.
#' @param event_num    The number of simulated random gap times for one patients. 
#'    Only the cumulative random gap times smaller the follow up time are saved as output. 
#'    This argument should be a numeric value.
#' @export
simData <- function(size, 
                    endTime,
                    X, 
                    kappaCoef, 
                    sigmaCoef, 
                    theta,  
                    gamma,
                    x0,
                    nu, 
                    roundTime = TRUE,
                    gapTime = FALSE,
                    event_num = 500){
  z1 <- rnorm(size, mean = 0, sd = sqrt(theta[1]))
  sigmas <- exp(X %*% sigmaCoef + z1)
  z2 <- rnorm(size, mean = 0, sd = sqrt(theta[2]))
  kappas <- x0 + exp(X %*% kappaCoef + gamma * z1 + z2)
  
  # generate hypo-event time and set status
  dat_evt <- data.frame('id' = NA, 'event' = NA,'time' = NA)
  for (i in 1:size){
    simGapTimes <- rFHT(n = event_num, x0 = x0, nu = nu, kappa = kappas[i], sigma = sigmas[i])
    if (roundTime) {
      simGapTimes <- round(simGapTimes)
    }
    time <- cumsum(simGapTimes)
    time_end_indx <- which(cumsum(simGapTimes) >= endTime[i])[1]
    time <- time[1:time_end_indx]
    time[time_end_indx] <- endTime[i]
    if (gapTime) {
      time <- simGapTimes[1:time_end_indx]
      time[length(time)] <- endTime[i] - sum(time[1:(time_end_indx - 1)] )
    }
    
    # set event status
    event <- c(rep(1, length(time) - 1), 0)
    dat_evt_new <- data.frame('id' = rep(i, length(time)),
                              'event' = event,
                              'time' = time)
    
    # combine event times for all patients
    dat_evt <- rbind(dat_evt, dat_evt_new)
  }
  
  # remove first row
  dat_evt <- dat_evt[!is.na(dat_evt$time), ]
  # merge hypo-event time and covariates
  X <- data.frame(X)
  colnames(X) <- paste0("X", c(1:ncol(X)))
  X$id <- c(1:size)
  dat_sim <- merge(dat_evt, X, by = 'id')
  return(dat_sim)
}



