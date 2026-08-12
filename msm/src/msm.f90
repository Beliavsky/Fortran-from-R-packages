! SPDX-License-Identifier: GPL-2.0-or-later
! Public convenience module for msm-fortran.
module msm
    use msm_kinds, only : dp
    use msm_linalg, only : expm, expm_frechet, inverse_matrix
    use msm_stats
    use msm_ctmc
    use msm_emissions
    use msm_hmm
    use msm_simulation
    use msm_distributions
    use msm_inference
    implicit none
    public
end module msm
