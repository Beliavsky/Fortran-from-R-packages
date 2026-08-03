! SPDX-License-Identifier: Artistic-2.0
module ecd_compat
  use ecd_core, only : dec => ecd_pdf, pec => ecd_cdf, qec => ecd_quantile, &
    rec => ecd_random, discr => ecd_discriminant, jinv => ecd_j_invariant, &
    solve_sym => ecd_solve_sym, solve_trig => ecd_solve_trig, y_slope => ecd_y_slope
  use ecd_processes, only : dlaplace0 => laplace_pdf, rlaplace0 => laplace_random, &
    dstdlap => stdlap_pdf, pstdlap => stdlap_cdf, qstdlap => stdlap_quantile, &
    rstdlap => stdlap_random, cfstdlap => stdlap_cf, kstdlap => stdlap_cumulants, &
    dstablecnt => stable_count_pdf, pstablecnt => stable_count_cdf, &
    qstablecnt => stable_count_quantile, rstablecnt => stable_count_random, &
    cfstablecnt => stable_count_cf, kstablecnt => stable_count_cumulants, &
    dsl => sld_pdf, psl => sld_cdf, qsl => sld_quantile, rsl => sld_random, &
    cfsl => sld_cf, ksl => sld_cumulants, dqsl => sld_pdf, pqsl => sld_cdf, &
    qqsl => sld_quantile, rqsl => sld_random, cfqsl => sld_cf, kqsl => sld_cumulants, &
    k2mnt => k2moments, mnt2k => moments2k
  use ecd_fitting, only : ecd_standardfit, fit_ecld_mle, qsld_fit
  use ecd_options, only : ecop_bs_option_price => bs_option_price, &
    ecop_bs_call_price => bs_call_price, ecop_bs_put_price => bs_put_price, &
    ecop_bs_implied_volatility => bs_implied_volatility, &
    ecop_polyfit_option => polyfit_option
  implicit none
  public
end module ecd_compat
