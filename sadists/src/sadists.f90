! SPDX-License-Identifier: LGPL-3.0-or-later
module sadists
    use sadists_kinds, only : dp
    use sadists_special, only : moments_to_cumulants, cumulants_to_moments, normal_moments, &
        chisq_moments, chisq_log_moment, digamma_dp, polygamma_dp, normal_cdf, normal_quantile
    use sadists_approximations, only : edgeworth_pdf, edgeworth_cdf, cornish_fisher_quantile, as269
    use sadists_distributions
    implicit none
    public
end module sadists
