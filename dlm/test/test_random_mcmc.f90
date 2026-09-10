! SPDX-License-Identifier: GPL-2.0-or-later
program test_random_mcmc
    use dlm, only : dp, dlm_success, erg_mean, mcmc_mean, mcmc_sd, rwishart
    use dlm, only : dlm_gibbs_dig, dlm_gibbs_dig_result, dlm_mod_poly, dlm_model
    implicit none

    real(dp) :: x(8, 2), means(2), sd(2), sd_one
    real(dp), allocatable :: erg(:, :), w1(:, :), w2(:, :)
    real(dp) :: y(8)
    type(dlm_model) :: model
    type(dlm_gibbs_dig_result) :: gibbs1, gibbs2, gibbs_mixed_prior
    integer :: i, info

    x(:, 1) = [1.0_dp, 2.0_dp, 1.0_dp, 3.0_dp, 2.0_dp, 4.0_dp, 3.0_dp, 5.0_dp]
    x(:, 2) = 2.0_dp * x(:, 1)
    call mcmc_mean(x, means, sd)
    if (maxval(abs(means - [2.625_dp, 5.25_dp])) > 1.0e-14_dp) error stop "mcmcMean means failed"
    if (any(sd < 0.0_dp)) error stop "mcmcMean SD failed"
    call mcmc_sd(x(:, 1), sd_one)
    if (sd_one < 0.0_dp) error stop "mcmcSD failed"

    call erg_mean(x, erg, info, m=3)
    if (info /= dlm_success) error stop "ergMean failed"
    if (any(shape(erg) /= [6, 2])) error stop "ergMean shape failed"
    if (abs(erg(1, 1) - 4.0_dp / 3.0_dp) > 1.0e-14_dp) error stop "ergMean first value failed"
    if (abs(erg(6, 1) - 2.625_dp) > 1.0e-14_dp) error stop "ergMean last value failed"

    call rwishart(6, 3, w1, info, seed=24680)
    if (info /= dlm_success) error stop "rwishart first draw failed"
    call rwishart(6, 3, w2, info, seed=24680)
    if (info /= dlm_success) error stop "rwishart second draw failed"
    if (maxval(abs(w1 - w2)) > 1.0e-14_dp) error stop "rwishart seed reproducibility failed"
    if (maxval(abs(w1 - transpose(w1))) > 1.0e-12_dp) error stop "rwishart symmetry failed"
    if (any([(w1(i, i) <= 0.0_dp, i = 1, 3)])) error stop "rwishart diagonal failed"

    y = [0.2_dp, 0.1_dp, 0.4_dp, 0.8_dp, 0.7_dp, 1.0_dp, 1.1_dp, 1.3_dp]
    call dlm_mod_poly(1, model, info, d_v=0.4_dp, d_w=[0.2_dp], m0=[0.0_dp], &
                      c0=reshape([1.0_dp], [1, 1]))
    if (info /= dlm_success) error stop "Gibbs model construction failed"
    call dlm_gibbs_dig(y, model, 2, gibbs1, info, thin=1, shape_y=2.0_dp, rate_y=1.0_dp, &
                       shape_theta=[2.0_dp], rate_theta=[1.0_dp], seed=13579)
    if (info /= dlm_success) error stop "dlmGibbsDIG first run failed"
    call dlm_gibbs_dig(y, model, 2, gibbs2, info, thin=1, shape_y=2.0_dp, rate_y=1.0_dp, &
                       shape_theta=[2.0_dp], rate_theta=[1.0_dp], seed=13579)
    if (info /= dlm_success) error stop "dlmGibbsDIG second run failed"
    if (any(gibbs1%d_v <= 0.0_dp) .or. any(gibbs1%d_w <= 0.0_dp)) error stop "dlmGibbsDIG positivity failed"
    if (.not. allocated(gibbs1%theta)) error stop "dlmGibbsDIG state storage failed"
    if (any(shape(gibbs1%theta) /= [9, 1, 2])) error stop "dlmGibbsDIG state shape failed"
    if (maxval(abs(gibbs1%d_v - gibbs2%d_v)) > 1.0e-14_dp) error stop "dlmGibbsDIG seed V failed"
    if (maxval(abs(gibbs1%d_w - gibbs2%d_w)) > 1.0e-14_dp) error stop "dlmGibbsDIG seed W failed"
    if (maxval(abs(gibbs1%theta - gibbs2%theta)) > 1.0e-14_dp) error stop "dlmGibbsDIG seed states failed"

    call dlm_mod_poly(2, model, info, d_v=0.4_dp, d_w=[0.2_dp, 0.1_dp], &
                      m0=[0.0_dp, 0.0_dp], c0=reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]))
    if (info /= dlm_success) error stop "Gibbs mixed-prior model construction failed"
    call dlm_gibbs_dig(y, model, 1, gibbs_mixed_prior, info, a_y=2.0_dp, b_y=1.0_dp, &
                       a_theta=[2.0_dp], b_theta=[1.0_dp, 2.0_dp], save_states=.false., seed=97531)
    if (info /= dlm_success) error stop "dlmGibbsDIG mixed prior recycling failed"
    if (any(shape(gibbs_mixed_prior%d_w) /= [1, 2])) error stop "dlmGibbsDIG mixed prior shape failed"
    if (any(gibbs_mixed_prior%d_w <= 0.0_dp)) error stop "dlmGibbsDIG mixed prior positivity failed"

    print *, "test_random_mcmc: PASS"
end program test_random_mcmc
