## Estimating overall variance for PSCL and PSKM estimations

#' Get overall variance for CL and KM estimations
#' including single arm, RCT, ps_cl and ps_km.
#'
#' @noRd
#'
get_ps_clkmfcn <- function(dta_psbor,
                           v_outcome     = NULL,
                           v_event       = NULL,
                           v_time        = NULL,
                           f_stratum     = get_cl_stratum,
                           f_overall_est = get_overall_est_wostderr,
                           f_ps_clkmfcn  = get_ps_cl_km,
                           stderr_method = c("naive", "jk", "sjk", "cjk",
                                             "sbs", "cbs", "none"), 
                           ...) {
    if (stderr_method[1] %in% c("naive", "jk", "none")) {
        rst <- f_ps_clkmfcn(dta_psbor,
                            v_outcome = v_outcome,
                            v_event = v_event, v_time = v_time,
                            f_stratum = f_stratum,
                            stderr_method = stderr_method[1],
                            ...)
    } else if (stderr_method[1] == "sjk") {
        rst <- get_ps_cl_km_sjk(dta_psbor,
                                v_outcome = v_outcome,
                                v_event = v_event, v_time = v_time,
                                f_stratum = f_stratum,
                                f_overall_est = f_overall_est,
                                f_ps_clkmfcn = f_ps_clkmfcn,
                                ...)
    } else if (stderr_method[1] == "cjk") {
        rst <- get_ps_cl_km_cjk(dta_psbor,
                                v_outcome = v_outcome,
                                v_event = v_event, v_time = v_time,
                                f_stratum = f_stratum,
                                f_overall_est = f_overall_est,
                                f_ps_clkmfcn = f_ps_clkmfcn,
                                ...)
    } else if (stderr_method[1] == "sbs") {
        rst <- get_ps_cl_km_sbs(dta_psbor,
                                v_outcome = v_outcome,
                                v_event = v_event, v_time = v_time,
                                f_stratum = f_stratum,
                                f_overall_est = f_overall_est,
                                f_ps_clkmfcn = f_ps_clkmfcn,
                                ...)
    } else if (stderr_method[1] == "cbs") {
        rst <- get_ps_cl_km_cbs(dta_psbor,
                                v_outcome = v_outcome,
                                v_event = v_event, v_time = v_time,
                                f_stratum = f_stratum,
                                f_overall_est = f_overall_est,
                                f_ps_clkmfcn = f_ps_clkmfcn,
                                ...)
    } else {
        stop("stderr_method is not implemented.")
    }

    return(rst)
}

