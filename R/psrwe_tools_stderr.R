#' Options for estimating stderr method
#'
#' Estimating the stderr may have several options which also depending on
#' outcome type and analysis methods.
#'
#' @details
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
#'     \code{naive} is mostly for time-to-event outcomes. The default is
#'     Greenwood formula or robust formula whichever is available.
#'
#'     \code{jk} is mostly the default for continuous and binary outcomes.
#'
#'     \code{none} will skip the stderr for speed up the calculation if
#'     only the point estimate is interested.
#'
#'
NULL
