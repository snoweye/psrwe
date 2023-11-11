#' Get estimates for Cox proportional hazard between two arms
#' for RCT augmenting control (stratified approach)
#'
#' @noRd
#'
get_ps_coxphst <- function(dta_psbor,
                           v_outcome     = NULL,
                           v_event       = NULL,
                           v_time        = NULL,
                           f_stratum     = NULL,
                           f_overall_est = get_overall_est_coxphst,
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
                             c(v_outcome, v_time, v_event))

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


