program test_otrimle
    use otrimle_mod
    use otrimle_linalg, only : scale_tau2_location_scale
    implicit none

    real(dp) :: x(12, 2)
    real(dp) :: kd_measure
    real(dp) :: tau_location
    real(dp) :: tau_scale
    real(dp), allocatable :: ddpm(:)
    real(dp), allocatable :: generated(:, :)
    integer, allocatable :: generated_labels(:)
    integer, allocatable :: init_auto(:)
    integer :: initial(12)
    type(otrimle_fit) :: fit
    type(otrimle_fit) :: optfit
    type(otrimle_grid_point), allocatable :: opt(:)
    type(kernel_density_result) :: kd
    type(otrimle_grid_result) :: gridres
    type(otrimle_sim_result) :: simres
    type(otrimle_sim_summary) :: simsum

    call make_data(x, initial)

    call assert_close(kmeanfun(10.0_dp), exp(-1.83084_dp) * 10.0_dp ** (-0.55323_dp), 1.0e-14_dp, &
        'kmeanfun small-n branch')
    call assert_close(kmeanfun(50.0_dp), exp(-2.4562342_dp) * 50.0_dp ** (-0.3524332_dp), 1.0e-14_dp, &
        'kmeanfun large-n branch')
    call assert_close(ksdfun(10.0_dp), exp(-2.1482546_dp) * 10.0_dp ** (-0.4639833_dp), 1.0e-14_dp, &
        'ksdfun small-n branch')
    call assert_close(ksdfun(50.0_dp), exp(-2.6575388_dp) * 50.0_dp ** (-0.4547216_dp), 1.0e-14_dp, &
        'ksdfun large-n branch')

    call scale_tau2_location_scale(&
        [-4.0_dp, -1.5_dp, -0.3_dp, 0.0_dp, 0.2_dp, 0.7_dp, 1.1_dp, 2.0_dp, 8.0_dp], &
        tau_location, tau_scale)
    call assert_close(tau_location, 0.3005828985635372_dp, 2.0e-15_dp, 'scaleTau2 robust location')
    call assert_close(tau_scale, 1.624757691871481_dp, 3.0e-15_dp, 'scaleTau2 robust scale')

    call init_clust(x, 2, init_auto, k=2, knnd_trim=0.2_dp)
    call assert_true(size(init_auto) == 12, 'InitClust result length')
    call assert_true(maxval(init_auto) == 2, 'InitClust finds requested cluster count')

    call rimle(x, 2, fit, initial=initial, logicd=-8.0_dp, npr_max=0.30_dp, erc=20.0_dp)
    call assert_true(fit%code == 2, 'rimle converges')
    call assert_true(fit%iter < 20, 'rimle iteration count')
    call assert_close(fit%exproportion(1), 1.0_dp / 6.0_dp, 2.0e-3_dp, 'rimle noise proportion')
    call assert_close(fit%mean(1, 1), -2.0_dp, 2.0e-2_dp, 'rimle cluster 1 x mean')
    call assert_close(fit%mean(1, 2), 2.02_dp, 2.0e-2_dp, 'rimle cluster 2 x mean')
    call assert_true(maxval(fit%cluster(11:12)) == 0, 'rimle classifies extreme points as noise')
    call assert_true(maxval(fit%cluster(1:5)) == 1, 'rimle first group assignment')
    call assert_true(minval(fit%cluster(6:10)) == 2, 'rimle second group assignment')

    call otrimle_fit_grid(x, 2, optfit, opt, initial=initial, &
        logicd_grid=[-12.0_dp, -10.0_dp, -8.0_dp, -6.0_dp], npr_max=0.30_dp, erc=20.0_dp, beta=0.2_dp)
    call assert_true(optfit%code > 0, 'otrimle grid fit succeeds')
    call assert_close(optfit%logicd, -12.0_dp, 1.0e-14_dp, 'otrimle chooses minimum criterion grid point')
    call assert_true(size(opt) == 4, 'otrimle optimization table length')

    call kerndensmeasure([-1.5_dp, -1.0_dp, -0.4_dp, 0.0_dp, 0.2_dp, 0.7_dp, 1.1_dp, 1.8_dp], &
        kd, kernn=20)
    call assert_close(kd%measure, 1.9920831119672024e-3_dp, 3.0e-15_dp, &
        'kerndensmeasure R density compatibility')
    call assert_close(kd%cpx(1), 1.7591648034453230e-3_dp, 3.0e-15_dp, &
        'kerndensmeasure first R density ordinate')
    call assert_true(size(kd%cp) == 10 .and. size(kd%cpx) == 20, 'kerndensmeasure output sizes')

    call kerndenscluster(x, fit, kd_measure, ddpm, kernn=20)
    call assert_true(kd_measure >= 0.0_dp, 'kerndenscluster nonnegative')
    call assert_true(size(ddpm) == 2, 'kerndenscluster component diagnostics')

    call set_otrimle_seed(20260914)
    call generator_otrimle(x, fit, generated, generated_labels)
    call assert_true(all(shape(generated) == [12, 2]), 'generator.otrimle data shape')
    call assert_true(size(generated_labels) == 12, 'generator.otrimle labels length')
    call assert_true(minval(generated_labels) >= 0 .and. maxval(generated_labels) <= 2, &
        'generator.otrimle label range')

    call otrimleg(x, [1, 2], gridres, fixlogicd=[-8.0_dp, -8.0_dp], kernn=20)
    call assert_true(size(gridres%solution) == 2, 'otrimleg candidate count')
    call assert_true(gridres%solution(2)%code > 0, 'otrimleg two-cluster fit')
    call assert_true(gridres%ibic(2) < huge(1.0_dp), 'otrimleg finite information criterion')

    call set_otrimle_seed(701)
    call otrimlesimg(x, [2], 2, simres, sim_est_logicd=.false.)
    call assert_true(simres%simruns == 2, 'otrimlesimg simulation count')
    call summarize_otrimlesimgdens(simres, simsum)
    call assert_true(simsum%best_g == 2, 'summary.otrimlesimgdens one-candidate selection')
    call assert_true(size(simsum%cluster) == 12, 'summary.otrimlesimgdens cluster output')

    print '(a)', 'All otrimle tests passed.'

contains

    subroutine make_data(x, initial)
        real(dp), intent(out) :: x(:, :) !! Synthetic two-cluster observations with two extreme noise points.
        integer, intent(out) :: initial(:) !! Known initial labels for the synthetic sample, with zero denoting noise.

        x(:, 1) = [-2.2_dp, -1.8_dp, -2.0_dp, -2.3_dp, -1.7_dp, &
            2.0_dp, 2.2_dp, 1.7_dp, 2.3_dp, 1.9_dp, 8.0_dp, -8.0_dp]
        x(:, 2) = [-2.0_dp, -2.1_dp, -1.7_dp, -2.2_dp, -1.9_dp, &
            2.1_dp, 1.8_dp, 2.0_dp, 2.2_dp, 1.7_dp, -7.0_dp, 7.0_dp]
        initial = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 0, 0]
    end subroutine make_data

    subroutine assert_true(condition, message)
        logical, intent(in) :: condition !! Boolean condition that must be true for the test to pass.
        character(len=*), intent(in) :: message !! Human-readable assertion description emitted on failure.

        if (.not. condition) then
            write (*, '(a)') 'FAILED: ' // trim(message)
            error stop 1
        end if
    end subroutine assert_true

    subroutine assert_close(actual, expected, tolerance, message)
        real(dp), intent(in) :: actual !! Computed scalar value under test.
        real(dp), intent(in) :: expected !! Reference scalar value expected by the test.
        real(dp), intent(in) :: tolerance !! Maximum allowed absolute difference between actual and expected.
        character(len=*), intent(in) :: message !! Human-readable assertion description emitted on failure.

        if (abs(actual - expected) > tolerance) then
            write (*, '(a,2(1x,es24.16))') 'FAILED: ' // trim(message), actual, expected
            error stop 1
        end if
    end subroutine assert_close

end program test_otrimle
