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

##' Simulated Covariate Data Set
##'
##' A simulated data set for variables named
##' \code{id}, \code{censorTime}, \code{x1}, and \code{x2}, where
##' \itemize{
##'     \item \code{id}: Subjects identification;
##'     \item \code{censorTime}: Subject follow-up time;
##'     \item \code{x1}: Continuous covariate x1;
##'     \item \code{x2}: Continuous covariate x2;
##' }
##'
##' @details
##' The sample dataset contains the subject-level covariates generated from
##' the fitted Copula model.
##'
##' @docType data
##' @name simuCovDat
##' @format A data frame with 400 rows and 4 variables.
NULL

##' Simulated Event Time Data Set
##'
##' The simulated event time data set contains the calender time of the
##' recurrent events, the event indicator and the subject covariates.
##' \code{id}, \code{event}, \code{time}, \code{x1} and \code{x2}, where
##' \itemize{
##'     \item \code{id}: The patient id;
##'     \item \code{event}: Event indicator, '1' for event, and '0' for censored;
##'     \item \code{time}: Event or censored time (day);
##'     \item \code{x1}; Continuous covariate x1;
##'     \item \code{x2}; Continuous covariate x2;
##' }
##' @docType data
##' @name simuEvtDat
##' @format A data frame with 3779 rows and 5 variables.
NULL
