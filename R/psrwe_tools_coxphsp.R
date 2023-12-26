#' Get estimates for Cox proportional hazard between two arms
#' for RCT augmenting control (same proportion approach)
#'
#' @noRd
#'
get_ps_coxphsp <- function(dta_psbor,
                           v_outcome     = NULL,
                           v_event       = NULL,
                           v_time        = NULL,
                           f_stratum     = NULL,
                           f_overall_est = get_overall_est_coxphsp,
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

    ## effect with borrowing
    cur_effect   <- get_surv_coxphsp(cur_d1, cur_d0, cur_d1t,
                                     n_borrow = borrow, ...)

    ## summary
    rst_effect <- f_overall_est(eff_theta)

    ## return
    rst <-  list(Control   = NULL,
                 Treatment = NULL,
                 Effect    = rst_effect,
                 Borrow    = dta_psbor$Borrow,
                 Total_borrow = dta_psbor$Total_borrow,
                 is_rct       = dta_psbor$is_rct)
    return(rst)
}

#' Summarize overall theta for coxph (same proportion approach)
#'
#'
#' @noRd
#'
get_overall_est_coxphsp <- function(ts1) {
    o_est <- data.frame(Mean   = ts1[, 1],
                        StdErr = ts1[, 2],
                        T      = ts1[, 3])

    list(Stratum_Estimate = NULL,
         Overall_Estimate = o_est)
}

