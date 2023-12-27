#' Options for estimating stderr method
#'
#' Estimating the stderr may have several options which also depend on
#' outcome type and analysis method.
#'
#' @section \code{stderr_method} options:
#'     \code{jk} is mostly the default for continuous and binary outcomes.
#'     \code{naive} is mostly the default for time-to-event outcomes via
#'     either Greenwood or robust formula whichever is applicable.
#'
#'     \code{jk} using Jackknife method within each stratum,
#'     \code{sjk} using simple Jackknife method for combined estimates
#'     such as point estimates in single arm or treatment effects in RCT, or
#'     \code{cjk} for complex Jackknife method including refitting PS model,
#'     matching, trimming, calculating borrowing parameters, and
#'     combining overall estimates.
#'
#'     Note that \code{sjk} may take a while longer to finish and
#'     \code{cjk} will take even much longer to finish.
#'
#'     The \code{sbs} and \code{cbs} are for simple and complex Bootstrap
#'     methods (\code{n_bootstrap = 200} as default).
#'
#'     \code{none} will skip the stderr for speed up the calculation if
#'     only the point estimate is interested.
#'
#' @docType package
#' @name stderr_method
NULL
