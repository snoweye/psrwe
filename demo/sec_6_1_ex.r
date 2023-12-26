### Example of Section 6.1.
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
rst_coxphwa_rct <- psrwe_survcoxphwa(ps_bor_rct,
                                     v_time = "Y_Surv",
                                     v_event = "Status")
rst_coxphwa_rct

### Outcome analysis.
oa_coxphwa_rct <- psrwe_outana(rst_coxphwa_rct, alternative = "less")
oa_coxphwa_rct
print(oa_coxphwa_rct, show_rct = TRUE)
summary(oa_coxphwa_rct)

### Use simple Jackknife stderr. This may take a while.
rst_coxphwa_rct_sjk <- psrwe_survcoxphwa(ps_bor_rct,
                                         v_time = "Y_Surv",
                                         v_event = "Status",
                                         stderr_method = "sjk")
oa_coxphwa_rct_sjk <- psrwe_outana(rst_coxphwa_rct_sjk, alternative = "less")
summary(oa_coxphwa_rct_sjk)

