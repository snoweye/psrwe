#' PS-integrated Cox proportional hazard method for
#' comparing time-to-event outcomes (same proportion approach)
#'
#' Cox proportional hazard (coxph) method evaluates two-arm RCT via
#' PS-integrated method (same proportion approach).
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
#'     \code{cjk} for complex Jackknife method including refitting PS model,
#'     matching, trimming, calculating borrowing parameters, and
#'     combining overall estimates.
#'     Note that \code{sjk} may take a while longer to finish and
#'     \code{cjk} will take even much longer to finish.
#'     The \code{sbs} and \code{cbs} are for simple and complex Bootstrap
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
#'     It contains the test statistics as well as the
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
#' rst_coxphsp <- psrwe_survcoxphsp(ps_bor_rct,
#'                                  v_time = "Y_Surv",
#'                                  v_event = "Status")
#' rst_coxphsp
#'
#' @export
#'
psrwe_survcoxphsp <- function(dta_psbor,
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
    rst_obs <- NULL

    ## call estimation
    if (stderr_method[1] %in% c("naive")) {
        rst <- get_ps_coxphsp(dta_psbor,
                              v_event = v_event, v_time = v_time,
                              stderr_method = stderr_method[1],
                              ...)
    } else if(stderr_method[1] == "sjk") {
        rst <- get_ps_coxphsp_sjk(dta_psbor,
                                  v_event = v_event, v_time = v_time,
                                  ...)
    } else if(stderr_method[1] == "cjk") {
        rst <- get_ps_coxphsp_cjk(dta_psbor,
                                  v_event = v_event, v_time = v_time,
                                  ...)
    } else if(stderr_method[1] == "sbs") {
        rst <- get_ps_coxphsp_sbs(dta_psbor,
                                  v_event = v_event, v_time = v_time,
                                  ...)
    } else if(stderr_method[1] == "cbs") {
        rst <- get_ps_coxphsp_cbs(dta_psbor,
                                  v_event = v_event, v_time = v_time,
                                  ...)
    } else {
        stop("stderr_errmethod is not implemented.")
    }

    ## return
    rst$Observed <- rst_obs
    rst$stderr_method <- stderr_method
    rst$Method   <- "ps_coxphsp"
    rst$Outcome_type <- "tte"
    class(rst)   <- get_rwe_class("ANARST")
    return(rst)
}

#' Get coxph estimation for each stratum (same proportion approach)
#'
#'
#' @noRd
#'
get_surv_coxphsp <- function(d1, d0 = NULL, d1t, n_borrow = 0,
                             ...) {

    ## treatment or control only
    dta_cur <- d1
    dta_ext <- d0
    dta_cur_trt <- d1t
    ns1     <- nrow(dta_cur)
    ns1_trt <- nrow(dta_cur_trt)

    if (is.null(d0)) {
        ns0 <- 0
    } else {
        ns0 <- nrow(dta_ext)
    }

    ## overall estimate
    overall  <- rwe_coxphwa(dta_cur, dta_ext, dta_cur_trt, n_borrow)

    return(overall)
}

#' The coxph estimation (same proportion approach)
#'
#' Estimate overall coxph estimate with all strata together
#'
#'
#' @param dta_cur Matrix of time and event from a PS stratum in current study
#'                (control arm only)
#' @param dta_ext Matrix of time and event from a PS stratum in external data
#'                source (control arm only)
#' @param dta_cur_trt Matrix of time and event from a PS stratum in current
#'                    study (treatment arm only)
#' @param n_borrow Number of subjects to be borrowed
#' @param stderr_method Method for computing StdErr (available for naive only)
#'
#' @return Estimation of overall coxph estimate
#'
#'
#' @export
#'
rwe_coxphsp <- function(dta_cur, dta_ext, dta_cur_trt, n_borrow = 0,
                        stderr_method = "naive") {

    ## current control and external control if available
    cur_data    <- dta_cur
    ns1         <- nrow(dta_cur)
    cur_weights <- rep(1, ns1)

    if (n_borrow > 0) {
        ns0         <- nrow(dta_ext)
        cur_data    <- rbind(cur_data, dta_ext)
        cur_weights <- c(cur_weights,
                         rep(n_borrow / ns0, ns0))
    }

    ## trt arm
    cur_data_trt    <- dta_cur_trt
    ns1_trt         <- nrow(dta_cur_trt)
    cur_weights_trt <- rep(1, ns1_trt)

    ## Combine data of two arms together
    cur_data <- cbind(cur_data, 0)
    cur_data_trt <- cbind(cur_data_trt, 1)
    cur_data_comb <- rbind(cur_data_trt, cur_data)
    cur_data_comb <- data.frame(cur_data_comb)
    colnames(cur_data_comb) <- c("time", "event", "stratum", "arm")

    ## w_i in the same order of cur_data_comb
    cur_weights_comb <- c(cur_weights_trt, cur_weights)

    ## Cox proportional hazard
    cur_coxph <- coxph(Surv(time, event) ~ arm + strata(stragum),
                       data   = cur_data_comb,
                       weight = cur_weights_comb,
                       ties   = "breslow",
                       robust = TRUE)
    mean_d <- cur_coxph$coefficients

    ## summary.coxph() does not do prediction, see predict.coxph().
    pred_tp <- max(cur_data_comb$time)

    ## for coxph naive stderr
    if (stderr_method == "naive") {
        ## robust se when "robust = TRUE" in coxph()
        stderr_d <- sqrt(cur_coxph$coefficients)
    } else {
        ## for none, jk, sjk, cjk, sbs, or cbs
        stderr_d <- NA
    }

    ## combine coxph estimates
    rst_coxph <- c(mean_d, stderr_d, pred_tp)

    colnames(rst_coxph) <- c("Mean", "StdErr", "T")
    return(rst_coxph)
}

