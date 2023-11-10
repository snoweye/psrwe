#' PS-integrated Cox proportional hazard method for
#' comparing time-to-event outcomes (stratified approach)
#'
#' Cox proportional hazard (coxph) method evaluates two-arm RCT via
#' PS-integrated method (stratified approach).
#' Variance can be estimated by Jackknife methods.
#' Apply to the case when there is only one external data source and
#' two-arm RCT.
#'
#' @inheritParams psrwe_survkm
#'
#' @param v_time Column name corresponding to event time
#' @param v_event Column name corresponding to event status
#' @param stderr_method Method for computing StdErr (see Details)
#' @param ... Additional Parameters
#'
#' @details \code{stderr_method} includes \code{naive} as default which
#'     mostly follows the calculation provided by
#'     \code{survival::coxph(..., robust = TURE)},
#'     \code{sjk} using simple Jackknife method for combined estimates
#'     such as point estimates in single arm or treatment effects in RCT, or
#'     \code{cjk} for complex Jackknife method including refitting PS model,
#'     matching, trimming, calculating borrowing parameters, and
#'     combining overall estimates.
#'     Note that \code{sjk} may take a while longer to finish and
#'     \code{cjk} will take even much longer to finish.
#'     The \code{sbs} and \code{cbs} is for simple and complex Bootstrap
#'     methods (\code{n_bootstrap = 200} as default).
#'
#'     The PS-integrated coxph method optimizes the composite partial
#'     likelihood to obtain the point estimate of hazard ratio between
#'     two arms.
#'     The implementation of the optimization is equivalent to estimate the
#'     hazard ratio of the weighted and stratified coxph model
#'     via \code{survival::coxph()}.
#'
#' @return A data frame with class name \code{PSRWE_RST_TESTANA}.
#'     It contains the test statistics of each stratum as well as the
#'     Jackknife estimation. The results can be further
#'     summarized by its S3 method \code{summary}.
#'     The results can be also analyzed by \code{psrwe_outana} for outcome
#'     analysis and inference.
#'
#'
#' @examples
#' data(ex_dta_rct)
#' dta_ps_rct <- psrwe_est(ex_dta_rct,
#'                         v_covs = paste("V", 1:7, sep = ""),
#'                         v_grp = "Group", cur_grp_level = "current",
#'                         v_arm = "Arm", ctl_arm_level = "control",
#'                         ps_method = "logistic", nstrata = 5,
#'                         stra_ctl_only = FALSE)
#' ps_bor_rct <- psrwe_borrow(dta_ps_rct, total_borrow = 30)
#' rst_coxphst <- psrwe_survcoxphst(ps_bor_rct,
#'                                  v_time = "Y_Surv",
#'                                  v_event = "Status")
#' rst_coxphst
#'
#' @export
#'
psrwe_survcoxphst <- function(dta_psbor,
                              v_time        = "time",
                              v_event       = "event",
                              stderr_method = c("naive", "sjk", "cjk",
                                                "sbs", "cbs"), 
                              ...) {

    ## check
    stopifnot(dta_psbor$is_rct)

    stopifnot(inherits(dta_psbor,
                       what = get_rwe_class("PSDIST")))

    stopifnot(all(c(v_event, v_time) %in%
                  colnames(dta_psbor$data)))

    stderr_method <- match.arg(stderr_method)

    ## all time points
    data    <- dta_psbor$data

    ## observed (no need so skip)
    # rst_obs <- get_coxph_observed(data, v_time, v_event)
    rst_obs <- NA

    ## call estimation
    f_get_ps_coxph <- switch(stderr_method[1],
                             sjk = get_ps_coxph_sjk,
                             cjk = get_ps_coxph_cjk,
                             sbs = get_ps_coxph_sbs,
                             cbs = get_ps_coxph_cbs,
                             naive = get_ps_coxph,
                             stop("stderr_method is not implemented."))

    rst <- f_get_ps_coxph(dta_psbor,
                          v_event = v_event, v_time = v_time,
                          ...)
    if (stderr_method[1] %in% c("naive")) {
        rst <- get_ps_coxph(dta_psbor,
                            v_event = v_event, v_time = v_time,
                            stderr_method = stderr_method[1],
                            ...)
    } else if(stderr_method[1] == "sjk") {
        rst <- get_ps_coxph_sjk(dta_psbor,
                                v_event = v_event, v_time = v_time,
                                ...)
    } else if(stderr_method[1] == "cjk") {
        rst <- get_ps_coxph_cjk(dta_psbor,
                                v_event = v_event, v_time = v_time,
                                ...)
    } else if(stderr_method[1] == "sbs") {
        rst <- get_ps_coxph_sbs(dta_psbor,
                                v_event = v_event, v_time = v_time,
                                ...)
    } else if(stderr_method[1] == "cbs") {
        rst <- get_ps_coxph_cbs(dta_psbor,
                                v_event = v_event, v_time = v_time,
                                ...)
    } else {
        stop("stderr_errmethod is not implemented.")
    }

    ## return
    rst$Observed <- rst_obs
    rst$stderr_method <- stderr_method
    rst$Method   <- "ps_coxphst"
    rst$Outcome_type <- "tte"
    class(rst)   <- get_rwe_class("ANARST")
    return(rst)
}


#' Get estimates for Cox proportional hazard between two arms
#' for RCT augmenting control (stratified approach)
#'
#' @noRd
#'
get_ps_coxph <- function(dta_psbor,
                         v_event   = NULL,
                         v_time    = NULL,
                         ...) {

    ## prepare data
    data    <- dta_psbor$data
    data    <- data[!is.na(data[["_strata_"]]), ]

    strata  <- levels(data[["_strata_"]])
    nstrata <- length(strata)
    borrow  <- dta_psbor$Borrow$N_Borrow

    ## rearrange data
    v_covs <- c(v_time, v_event, "_grp_", "_arm_")
    cur_d1  <- NULL
    cur_d0  <- NULL
    cur_d1t <- NULL
    for (i in seq_len(nstrata)) {
        cur_01  <- get_cur_d(data, strata[i], v_covs)

        cur_d1  <- rbind(cur_d1, cur_01$cur_d1)    ## This is "cur_d1c"
        cur_d0  <- rbind(cur_d0, cur_01$cur_d0)
        cur_d1t <- rbind(cur_d1t, cur_01$cur_d1t)
    }
### I am here.

    ## effect with borrowing
    cur_effect   <- get_surv_coxphst(cur_d1, cur_d0, cur_d1t,
                                     n_borrow = borrow, ...)

    ## summary
    rst_effect <- f_overall_est(eff_theta, dta_psbor$Borrow$N_Current)

    ## return
    rst <-  list(Control   = NULL,
                 Treatment = NULL,
                 Effect    = rst_effect,
                 Borrow    = dta_psbor$Borrow,
                 Total_borrow = dta_psbor$Total_borrow,
                 is_rct       = dta_psbor$is_rct)
    return(rst)
}

