program mice_parity_example
    use iso_fortran_env, only : int64
    use mice, only : dp, mice_ok, mice_not_converged, mice_rng_state, rng_seed, impute_polr, impute_midastouch
    implicit none

    real(dp) :: y(18), x(18, 1)
    real(dp), allocatable :: yimp(:)
    integer :: category(18), info, i
    integer, allocatable :: cimp(:)
    logical :: observed(18), where(18)
    type(mice_rng_state) :: rng

    do i = 1, 18
        x(i, 1) = real(i - 9, dp) / 4.0_dp
        y(i) = 3.0_dp + 0.8_dp * x(i, 1) + 0.05_dp * cos(real(i, dp))
        if (i <= 6) then
            category(i) = 1
        else if (i <= 12) then
            category(i) = 2
        else
            category(i) = 3
        end if
    end do
    observed = .true.
    observed([6, 12, 18]) = .false.
    where = .not. observed

    call rng_seed(rng, 20260831_int64)
    call impute_midastouch(y, observed, x, where, rng, yimp, info)
    if (info /= mice_ok) error stop "midastouch example failed"
    print *, "midastouch imputations:", yimp

    call rng_seed(rng, 20260831_int64)
    call impute_polr(category, observed, x, where, 3, rng, cimp, info)
    if (info /= mice_ok .and. info /= mice_not_converged) error stop "polr example failed"
    print *, "polr imputations:", cimp
end program mice_parity_example
