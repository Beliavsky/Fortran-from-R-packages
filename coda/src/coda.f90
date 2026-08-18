! coda-fortran: computational translation of the R package coda.
! Original coda license: GPL (>= 2). This translation is GPL-2.0-or-later.
module coda
   use coda_kinds, only : dp
   use coda_types, only : mcmc_chain, mcmc_list, make_mcmc, make_mcmc_list, window_mcmc, pool_chains
   use coda_acf, only : autocorr, autocorr_diag, autocorr_diag_list, crosscorr, rejection_rate, rejection_rate_list
   use coda_spectrum, only : spectrum_ar_result, spectrum0_ar, spectrum0, effective_size, effective_size_list
   use coda_summary, only : mcmc_summary, hpd_interval, batch_se, batch_se_list, summarize_mcmc, summarize_mcmc_list
   use coda_diagnostics, only : geweke_result, gelman_result, heidel_result, raftery_result, &
                                geweke_diag, gelman_diag, heidel_diag, raftery_diag
   use coda_math, only : cramer_cdf
   implicit none
   public
end module coda
