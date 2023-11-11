#' PS-integrated Cox proportional hazard method for
#' comparing time-to-event outcomes (weighted average approach)
#'
#' Cox proportional hazard (coxph) method evaluates two-arm RCT via
#' PS-integrated method (weighted average approach).
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
#'     \code{survival::coxph(..., robust = TRUE)}, and
#'     \code{jk} using Jackknife method within each stratum.
#'
#'     Naive approach calculates log hazard ratio by each stratum, then
#'     takes weighted average of all stratum-specific estimates as other
#'     PS-integrated methods in \pkg{psrwe}.
#'
#' @return A data frame with class name \code{PSRWE_RST_TESTANA}.
#'     It contains the test statistics of each stratum.
#'     The results can be further
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
#' rst_coxphwa <- psrwe_survcoxphwa(ps_bor_rct,
#'                                  v_time = "Y_Surv",
#'                                  v_event = "Status")
#' rst_coxphwa
#'
#' @export
#'
psrwe_survcoxphwa <- function(dta_psbor,
                              v_time        = "time",
                              v_event       = "event",
                              stderr_method = c("naive", "jk"), 
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
    if (stderr_method[1] %in% c("naive", "jk")) {
        rst <- get_ps_coxphwa(dta_psbor,
                              v_event = v_event, v_time = v_time,
                              stderr_method = stderr_method[1],
                              ...)
    } else {
        stop("stderr_errmethod is not implemented.")
    }

    ## return
    rst$Observed <- rst_obs
    rst$stderr_method <- stderr_method
    rst$Method   <- "ps_coxphwa"
    rst$Outcome_type <- "tte"
    class(rst)   <- get_rwe_class("ANARST")
    return(rst)
}

#' Get estimates for Cox proportional hazard between two arms
#' for RCT augmenting control (weighted average approach)
#'
#' @noRd
#'
get_ps_coxphwa <- function(dta_psbor,
                           v_event       = NULL,
                           v_time        = NULL,
                           f_stratum     = get_surv_stratum_coxphwa,
                           f_overall_est = get_overall_est,
                           ...) {

    ## prepare data
    data    <- dta_psbor$data
    data    <- data[!is.na(data[["_strata_"]]), ]

    strata  <- levels(data[["_strata_"]])
    nstrata <- length(strata)
    borrow  <- dta_psbor$Borrow$N_Borrow

    ## estimate
    eff_theta <- NULL
    for (i in seq_len(nstrata)) {
        cur_01  <- get_cur_d(data,
                             strata[i],
                             c(v_time, v_event))

        cur_d1  <- cur_01$cur_d1    ## This is "cur_d1c"
        cur_d0  <- cur_01$cur_d0
        cur_d1t <- cur_01$cur_d1t

        ## effect with borrowing
        cur_effect   <- f_stratum(cur_d1, cur_d0, cur_d1t,
                                  n_borrow = borrow[i], ...)
        eff_theta <- rbind(eff_theta, cur_effect)
    }

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

#' Get coxph estimation for each stratum (weighted average approach)
#'
#'
#' @noRd
#'
get_surv_stratum_coxphwa <- function(d1, d0 = NULL, d1t, n_borrow = 0,
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
    overall  <- rwe_coxphwa(dta_cur, dta_ext, dta_cur_trt, n_borrow,
                            stderr_method)

    ## jackknife stderr
    if (stderr_method == "jk") {
        overall_theta <- overall[, 1, drop = TRUE]

        jk_theta      <- rep(0, length(overall_theta))
        for (j in seq_len(ns1)) {
            cur_jk   <- rwe_coxphwa(dta_cur[-j, ], dta_ext, dta_cur_trt, n_borrow,
                                    stderr_method)
            jk_theta <- jk_theta + (cur_jk[, 1] - overall_theta)^2
        }

        for (j in seq_len(ns1_trt)) {
            cur_jk   <- rwe_coxphwa(dta_cur, dta_ext, dta_cur_trt[-j, ], n_borrow,
                                    stderr_method)
            jk_theta <- jk_theta + (cur_jk[, 1] - overall_theta)^2
        }

        if (ns0 > 0) {
            for (j in seq_len(ns0)) {
                ext_jk   <- rwe_coxphwa(dta_cur, dta_ext[-j, ], dta_cur_trt,
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

#' The coxph estimation (weighted average approach)
#'
#' Estimate coxph estimates for a single PS stratum
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
rwe_coxphwa <- function(dta_cur, dta_ext, dta_cur_trt, n_borrow = 0,
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
    colnames(cur_data_comb) <- c("time", "event", "arm")

    ## w_i in the same order of cur_data_comb
    cur_weights_comb <- c(cur_weights_trt, cur_weights)

    ## Cox proportional hazard
    cur_coxph <- coxph(Surv(time, event) ~ arm,
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

