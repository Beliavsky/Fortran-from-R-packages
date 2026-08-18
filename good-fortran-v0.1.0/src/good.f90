! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of computational code from R package good 1.0.2.

module good
    use good_kinds, only : dp
    use good_distribution, only : dgood, pgood, qgood, rgood, goodmean, good_moments, good_logpmf
    use good_glm, only : good_glm_fit, good_prediction, good_glm_summary, glm_good, predict_good, summary_good, good_loglik
    implicit none
    public
end module good
