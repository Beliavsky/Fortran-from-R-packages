program test_parity_methods
    use iso_fortran_env, only : int64
    use mice, only : dp, mice_ok, mice_not_converged, mice_rng_state, rng_seed, proportional_odds_fit, &
                     proportional_odds_probabilities, impute_polr, impute_midastouch, impute_mpmm, &
                     impute_mnar_norm, impute_mnar_logreg, impute_lda
    implicit none

    call test_polr_kernel()
    call test_midastouch_kernel()
    call test_mpmm_kernel()
    call test_mnar_kernels()
    call test_lda_kernel()
    print *, "test_parity_methods: PASS"

contains

    subroutine test_polr_kernel()
        real(dp) :: x(60, 1), weights(60), one(4, 1), probs(4, 3)
        real(dp), allocatable :: beta(:), thresholds(:)
        integer :: category(60), info, i
        integer, allocatable :: imputed(:)
        logical :: observed(60), where(60)
        type(mice_rng_state) :: rng

        do i = 1, 60
            x(i, 1) = real(i - 30, dp) / 10.0_dp
            if (x(i, 1) < -0.5_dp) then
                category(i) = 1
            else if (x(i, 1) < 0.8_dp) then
                category(i) = 2
            else
                category(i) = 3
            end if
        end do
        weights = 1.0_dp
        call proportional_odds_fit(x, category, weights, 3, beta, thresholds, info)
        if (info /= mice_ok .and. info /= mice_not_converged) error stop "polr fit status"
        if (size(beta) /= 1 .or. beta(1) <= 0.0_dp) error stop "polr slope direction"
        if (thresholds(2) <= thresholds(1)) error stop "polr threshold order"
        one(:, 1) = [-2.0_dp, -0.2_dp, 0.7_dp, 2.0_dp]
        call proportional_odds_probabilities(one, beta, thresholds, probs)
        do i = 1, 4
            if (abs(sum(probs(i, :)) - 1.0_dp) > 1.0e-12_dp) error stop "polr probability sum"
            if (any(probs(i, :) < 0.0_dp)) error stop "polr negative probability"
        end do
        observed = .true.
        observed(55:60) = .false.
        where = .not. observed
        call rng_seed(rng, 101_int64)
        call impute_polr(category, observed, x, where, 3, rng, imputed, info)
        if (info /= mice_ok .and. info /= mice_not_converged) error stop "polr impute status"
        if (any(imputed < 1) .or. any(imputed > 3)) error stop "polr impute support"
    end subroutine test_polr_kernel

    subroutine test_midastouch_kernel()
        real(dp) :: y(30), x(30, 1), neff
        real(dp), allocatable :: imputed(:)
        logical :: observed(30), where(30)
        type(mice_rng_state) :: rng
        integer :: i, info, j
        logical :: found

        do i = 1, 30
            x(i, 1) = real(i, dp) / 10.0_dp
            y(i) = 2.0_dp + 1.5_dp * x(i, 1) + 0.1_dp * sin(real(i, dp))
        end do
        observed = .true.
        observed(26:30) = .false.
        where = .not. observed
        call rng_seed(rng, 202_int64)
        call impute_midastouch(y, observed, x, where, rng, imputed, info, neff=neff)
        if (info /= mice_ok) error stop "midastouch status"
        if (neff <= 0.0_dp) error stop "midastouch neff"
        do i = 1, size(imputed)
            found = .false.
            do j = 1, 25
                if (abs(imputed(i) - y(j)) <= epsilon(1.0_dp) * max(1.0_dp, abs(y(j)))) found = .true.
            end do
            if (.not. found) error stop "midastouch donor preservation"
        end do
    end subroutine test_midastouch_kernel

    subroutine test_mpmm_kernel()
        real(dp) :: y(30, 2), x(30, 2)
        real(dp), allocatable :: imputed(:, :)
        logical :: observed(30), where(30)
        type(mice_rng_state) :: rng
        integer :: i, info, j
        logical :: found

        do i = 1, 30
            x(i, 1) = real(i, dp) / 10.0_dp
            x(i, 2) = cos(real(i, dp) / 5.0_dp)
            y(i, 1) = 1.0_dp + x(i, 1) + 0.2_dp * x(i, 2)
            y(i, 2) = -0.5_dp + 2.0_dp * x(i, 1) - 0.1_dp * x(i, 2)
        end do
        observed = .true.
        observed(25:30) = .false.
        where = .not. observed
        call rng_seed(rng, 303_int64)
        call impute_mpmm(y, observed, x, where, rng, imputed, info, donors=4)
        if (info /= mice_ok) error stop "mpmm status"
        do i = 1, size(imputed, 1)
            found = .false.
            do j = 1, 24
                if (maxval(abs(imputed(i, :) - y(j, :))) <= 0.0_dp) found = .true.
            end do
            if (.not. found) error stop "mpmm row donor preservation"
        end do
    end subroutine test_mpmm_kernel

    subroutine test_mnar_kernels()
        real(dp) :: y(20), x(20, 1), u(20, 1), delta0(1), delta2(1)
        real(dp), allocatable :: base(:), shifted(:), binary_base(:), binary_shifted(:)
        logical :: observed(20), where(20)
        type(mice_rng_state) :: rng1, rng2
        integer :: i, info

        do i = 1, 20
            x(i, 1) = real(i, dp) / 10.0_dp
            y(i) = 1.0_dp + 0.5_dp * x(i, 1)
            u(i, 1) = 1.0_dp
        end do
        observed = .true.
        observed(16:20) = .false.
        where = .not. observed
        delta0 = 0.0_dp
        delta2 = 2.0_dp
        call rng_seed(rng1, 404_int64)
        call rng_seed(rng2, 404_int64)
        call impute_mnar_norm(y, observed, x, where, u, delta0, rng1, base, info)
        if (info /= mice_ok) error stop "mnar norm base status"
        call impute_mnar_norm(y, observed, x, where, u, delta2, rng2, shifted, info)
        if (info /= mice_ok) error stop "mnar norm shifted status"
        if (maxval(abs((shifted - base) - 2.0_dp)) > 1.0e-11_dp) error stop "mnar normal offset"

        do i = 1, 20
            if (i <= 10) then
                y(i) = 0.0_dp
            else
                y(i) = 1.0_dp
            end if
        end do
        call rng_seed(rng1, 505_int64)
        call rng_seed(rng2, 505_int64)
        call impute_mnar_logreg(y, observed, x, where, u, delta0, rng1, binary_base, info)
        if (info /= mice_ok .and. info /= mice_not_converged) error stop "mnar logreg base status"
        delta2 = 20.0_dp
        call impute_mnar_logreg(y, observed, x, where, u, delta2, rng2, binary_shifted, info)
        if (info /= mice_ok .and. info /= mice_not_converged) error stop "mnar logreg shifted status"
        if (any(binary_shifted < binary_base)) error stop "mnar logreg positive sensitivity"
    end subroutine test_mnar_kernels

    subroutine test_lda_kernel()
        real(dp) :: x(36, 2)
        integer :: category(36), info, i
        integer, allocatable :: imputed(:)
        logical :: observed(36), where(36)
        type(mice_rng_state) :: rng

        do i = 1, 12
            category(i) = 1
            x(i, :) = [-2.0_dp + 0.05_dp * real(i, dp), -1.0_dp]
            category(12 + i) = 2
            x(12 + i, :) = [0.05_dp * real(i, dp), 1.0_dp]
            category(24 + i) = 3
            x(24 + i, :) = [2.0_dp + 0.05_dp * real(i, dp), -0.5_dp]
        end do
        observed = .true.
        observed([10, 22, 34]) = .false.
        where = .not. observed
        call rng_seed(rng, 707_int64)
        call impute_lda(category, observed, x, where, 3, rng, imputed, info)
        if (info /= mice_ok) error stop "lda status"
        if (any(imputed < 1) .or. any(imputed > 3)) error stop "lda support"
    end subroutine test_lda_kernel

end program test_parity_methods
