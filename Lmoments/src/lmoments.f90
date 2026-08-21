! SPDX-License-Identifier: GPL-2.0-only
module lmoments
    use lmoments_utils, only: dp, pi, standard_normal_quantile
    use lmoments_core, only: lmoments_sample, lmoments_matrix, lcoefs_sample, &
        lmom_cov, t1_lmoments, shifted_legendre, hosking_lmoments
    use lmoments_quantile_mixtures, only: lmom2normpoly4, lmom2normpoly6, &
        data2normpoly4, data2normpoly6, qnormpoly, pnormpoly, dnormpoly, &
        rnormpoly, normpoly_inv, normpoly_cdf, normpoly_pdf, normpoly_rnd, &
        covnormpoly4, t1lmom2cauchypoly4, data2cauchypoly4, qcauchypoly, &
        pcauchypoly, dcauchypoly, rcauchypoly, cauchypoly_inv, &
        cauchypoly_cdf, cauchypoly_pdf, cauchypoly_rnd
    implicit none
    public
end module lmoments
