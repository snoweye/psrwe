#' PS-Integrated Cox Proportional Hazard Method For Comparing Time-to-event Outcomes
#'
#' Cox proportional hazard (coxph) method evaluates two-arm RCT for up to
#' a given time point.
#' Variance can be estimated by Jackknife methods.
#' Apply to the case when there is only one external data source and
#' two-arm RCT.
#'
#' @inheritParams psrwe_survkm
#'
#' @param v_time Column name corresponding to event time
#' @param v_event Column name corresponding to event status
#' @param stderr_method Method for computing StdErr (see Details)
#' @param n_bootstrap Number of bootstrap samples (for bootstrap stderr)
#' @param ... Additional Parameters
#'
#' @details \code{stderr_method} includes \code{naive} as default which
#'     mostly follows Greenwood formula,
#'     \code{jk} using Jackknife method within each stratum,
#'     \code{sjk} using simple Jackknife method for combined estimates
#'     such as point estimates in single arm or treatment effects in RCT, or
#'     \code{cjk} for complex Jackknife method including refitting PS model,
#'     matching, trimming, calculating borrowing parameters, and
#'     combining overall estimates.
#'     Note that \code{sjk} may take a while longer to finish and
#'     \code{cjk} will take even much longer to finish.
#'     The \code{sbs} and \code{cbs} is for simple and complex Bootstrap
#'     methods.
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
#' rst_coxph <- psrwe_survcoxph(ps_bor_rct,
#'                              v_time = "Y_Surv",
#'                              v_event = "Status")
#' rst_coxph
#'
#' @export
#'
psrwe_survcoxph <- function(dta_psbor,
                            v_time        = "time",
                            v_event       = "event",
                            stderr_method = c("naive", "jk", "sjk", "cjk",
                                              "sbs", "cbs", "none"), 
                            n_bootstrap = 200,
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

    ## observed
    rst_obs <- get_coxph_observed(data, v_time, v_event)

    ## call estimation
    f_get_ps_coxph <- switch(stderr_method[1],
                             jk = get_ps_coxph,
                             sjk = get_ps_coxph_sjk,
                             cjk = get_ps_coxph_cjk,
                             sbs = get_ps_coxph_sbs,
                             cbs = get_ps_coxph_cbs,
                             none = get_ps_coxph_none,
                             naive = get_ps_coxph_naive,
                             stop("stderr_method is not implemented."))

    rst <- f_get_ps_coxph(dta_psbor,
                          v_event = v_event, v_time = v_time,
                          f_stratum = get_surv_stratum_coxph,
                          n_bootstrap = n_bootstrap,
                          ...)

    ## return
    rst$Observed <- rst_obs
    rst$stderr_method <- stderr_method
    rst$Method   <- "ps_coxph"
    rst$Outcome_type <- "tte"
    class(rst)   <- get_rwe_class("ANARST")
    return(rst)
}

#' Get coxph estimation for each stratum
#'
#'
#' @noRd
#'
get_surv_stratum_coxph <- function(d1, d0 = NULL, d1t, n_borrow = 0,
                                   stderr_method = "jk", ...) {

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
    overall  <- rwe_coxph(dta_cur, dta_ext, dta_cur_trt, n_borrow,
                          stderr_method)

    ## jackknife stderr
    if (stderr_method == "jk") {
        overall_theta <- overall[, 1, drop = TRUE]

        jk_theta      <- rep(0, length(overall_theta))
        for (j in seq_len(ns1)) {
            cur_jk   <- rwe_coxph(dta_cur[-j, ], dta_ext, dta_cur_trt, n_borrow,
                                  stderr_method)
            jk_theta <- jk_theta + (cur_jk[, 1] - overall_theta)^2
        }

        for (j in seq_len(ns1_trt)) {
            cur_jk   <- rwe_coxph(dta_cur, dta_ext, dta_cur_trt[-j, ], n_borrow,
                                  stderr_method)
            jk_theta <- jk_theta + (cur_jk[, 1] - overall_theta)^2
        }

        if (ns0 > 0) {
            for (j in seq_len(ns0)) {
                ext_jk   <- rwe_coxph(dta_cur, dta_ext[-j, ], dta_cur_trt,
                                      n_borrow, stderr_method)
                jk_theta <- jk_theta + (ext_jk[, 1] - overall_theta)^2
            }
        }

        ## summary
        nc_jk <- (ns1 + ns0 + ns1_trt - 1) / (ns1 + ns0 + ns1_trt)
        stderr_theta <- sqrt(jk_theta * nc_jk)

        overall[, 2] <- stderr_theta
    }

    return(overall)
}

#' coxph Estimation
#'
#' Estimate coxph estimates for a single PS
#' stratum
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
#' @return Estimation of coxph estimates
#'
#'
#' @export
#'
rwe_coxph <- function(dta_cur, dta_ext, dta_cur_trt, n_borrow = 0,
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

    ## KM with stratum weights
    colnames(cur_data) <- c("time", "event")
    cur_data <- data.frame(cur_data)
    cur_surv <- survfit(Surv(time, event) ~ 1,
                        data      = cur_data,
                        weights   = cur_weights,
                        conf.type = "none")

    ## trt arm
    cur_data_trt    <- dta_cur_trt
    ns1_trt         <- nrow(dta_cur_trt)
    cur_weights_trt <- rep(1, ns1_trt)

    ## KM with stratum weights for trt arm
    colnames(cur_data_trt) <- c("time", "event")
    cur_data_trt <- data.frame(cur_data_trt)
    cur_surv_trt <- survfit(Surv(time, event) ~ 1,
                            data      = cur_data_trt,
                            weights   = cur_weights_trt,
                            conf.type = "none")

    ## summary.coxph() does not do prediction, see predict.coxph().
    pred_tps <- max(c(cur_surv$time, cur_surv_trt$time))

    rst <- summary(cur_surv)
    rst_trt <- summary(cur_surv_trt)

    ## info needed for coxph
    n_risk_trt <- rst_trt$n.risk
    n_risk_ctl <- rst$n.risk
    n_event_trt <- rst_trt$n.event
    n_event_ctl <- rst$n.event

    n_risk <- n_risk_trt + n_risk_ctl
    n_event <- n_event_trt + n_event_ctl
    p_event <- ifelse(n_risk == 0, 0, n_event / n_risk)
    E_1_j <- n_risk_trt * p_event

    ## coxph main statistic
    mean_d <- n_event_trt - E_1_j
    mean_d <- cumsum(mean_d)

    ## for coxph naive stderr
    if (stderr_method == "naive") {
        stderr_d <- ifelse(n_risk <= 1, 0,
                           E_1_j * (1 - p_event) *
                           n_risk_ctl / (n_risk - 1))
        stderr_d <- sqrt(cumsum(stderr_d))
    } else {
        ## for none, jk, sjk, cjk, sbs, or cbs
        stderr_d <- rep(NA, length(mean_d))
    }

    ## combine coxph estimates
    rst_coxph <- cbind(mean_d, stderr_d, pred_tps)

    colnames(rst_coxph) <- c("Mean", "StdErr", "T")
    return(rst_coxph)
}

