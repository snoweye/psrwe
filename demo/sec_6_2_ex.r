### Example of Section 6.2.
suppressMessages(library(psrwe, quietly = TRUE))
options(digits = 3)
data(ex_dta_rct)

### First parts of Data.
head(ex_dta_rct)

### Obtain PSs.
dta_ps_rct <- psrwe_est(ex_dta_rct,
                        v_covs = paste("V", 1:7, sep = ""),
                        v_grp = "Group", cur_grp_level = "current",
                        v_arm = "Arm", ctl_arm_level = "control",
                        ps_method = "logistic", nstrata = 5,
                        stra_ctl_only = FALSE)

### Balance assessment of PS stratification.
### See "sec_4_2_ex" for details.

### Obtain discounting parameters.
### See "sec_4_2_ex" for details.
ps_bor_rct <- psrwe_borrow(dta_ps_rct, total_borrow = 30)

### PSCOXPHWA, two-arm RCT, time-to-event outcome, coxph weighted average.
rst_coxphsp_rct <- psrwe_survcoxphsp(ps_bor_rct,
                                     v_time = "Y_Surv",
                                     v_event = "Status")
rst_coxphsp_rct

### Outcome analysis.
oa_coxphsp_rct <- psrwe_outana(rst_coxphsp_rct, alternative = "less")
oa_coxphsp_rct
summary(oa_coxphsp_rct)

### Use simple Jackknife stderr. This may take a while.
rst_coxphsp_rct_sjk <- psrwe_survcoxphsp(ps_bor_rct,
                                         v_time = "Y_Surv",
                                         v_event = "Status",
                                         stderr_method = "sjk")
oa_coxphsp_rct_sjk <- psrwe_outana(rst_coxphsp_rct_sjk, alternative = "less")
summary(oa_coxphsp_rct_sjk)

