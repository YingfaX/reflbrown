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

##
## R package rebrown by Author and Author
## Copyright (C) 2022
##
## This file is part of the R package rebrown.
##
## The R package rebrown is free software: You can redistribute it and/or
## modify it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or any later
## version (at your option). See the GNU General Public License at
## <https://www.gnu.org/licenses/> for details.
##
## The R package rebrown is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##

##
## R package rebrown by Yingfa Xie and Jun Yan
## Copyright (C) 2022
##
## This file is part of the R package rebrown.
##
## The R package rebrown is free software: You can redistribute it and/or
## modify it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or any later
## version (at your option). See the GNU General Public License at
## <https://www.gnu.org/licenses/> for details.
##
## The R package rebrown is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##


#' @useDynLib rebrown, .registration = TRUE
#' @importFrom ggplot2 ggplot geom_line geom_hline theme_minimal labs aes
NULL

#' @title Reflected Brownian Motion Path
#' @description Simulating sample path of reflected Brownian motion.
#'     The reflected Brownian motion begins at starting point. The process
#'     stops when it reach the lower boundary and then reset the process to
#'     be the value of the starting point.
#'
#'
#' @param x0    The starting point of reflected Brownian motion, which is a scalar.
#' @param nu    The lower boundary of reflected Brownian motion, which is a scalar.
#' @param kappa The upper reflection barrier of reflected Brownian motion,
#'               which is a scalar.
#' @param sigma The volatility of reflected Brownian motion, which is a scalar.
#' @param time_idx The discrete time index, which is a vector
#' @param plot_path Plot the path if TRUE. Default is FALSE.
#' @return A matrix contains three column, time_idx, Xt, proc_idx. The column
#'    proc_idx indicates which process the Brownian motion belongs to.
#'
#' @examples
#' set.seed(1)
#' time_idx <- seq(0, 200, 1)
#' sim_list <- simPath(x0 = 10, nu = 3.9, kappa = 20, sigma = 1.5,
#'                     time_idx, plot_path = TRUE)

#' @export
simPath <- function(x0, nu, kappa, sigma, time_idx, plot_path = FALSE){
  dt <- diff(time_idx)
  simPathProc <- simPath_c(x0 = x0, nu = nu, kappa = kappa,
                           sigma = sigma, dt = dt)
  simPathProc_mat <- cbind(time_idx, as.matrix(as.data.frame(simPathProc)))
  if (plot_path){
    Xt <- proc_idx <- NULL
    simPathProc_df <- as.data.frame(simPathProc_mat)
    p <- ggplot(simPathProc_df) +
      geom_line(aes(x = time_idx, y = Xt, group = proc_idx)) +
      geom_hline(yintercept=kappa, linetype="longdash", color = "lightblue") +
      geom_hline(yintercept=nu, linetype="dashed", color = "blue") +
      theme_minimal() +
      labs(x = expression(t), y = expression(X(t)))
    print(p)
  }
  return(simPathProc_mat)
}


