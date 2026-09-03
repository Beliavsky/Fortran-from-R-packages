program test_categorical
    use iso_fortran_env, only : int64
    use mice, only : dp, mice_ok, mice_rng_state, rng_seed, impute_logreg, impute_polyreg
    implicit none
    real(dp) :: y(12), x(12, 1)
    real(dp), allocatable :: imp(:)
    integer :: cat(12)
    integer, allocatable :: cat_imp(:)
    logical :: observed(12), where(12)
    type(mice_rng_state) :: rng
    integer :: i, info

    do i = 1, 12
        x(i, 1) = real(i - 6, dp) / 3.0_dp
        if (i <= 6) then
            y(i) = 0.0_dp
        else
            y(i) = 1.0_dp
        end if
        cat(i) = 1 + mod(i - 1, 3)
    end do
    observed = .true.
    observed(10:12) = .false.
    where = .not. observed
    call rng_seed(rng, 321_int64)
    call impute_logreg(y, observed, x, where, rng, imp, info)
    if (info /= mice_ok) error stop "logreg status"
    if (any(imp < 0.0_dp) .or. any(imp > 1.0_dp)) error stop "logreg support"
    call rng_seed(rng, 777_int64)
    call impute_polyreg(cat, observed, x, where, 3, rng, cat_imp, info)
    if (info /= mice_ok) error stop "polyreg status"
    if (any(cat_imp < 1) .or. any(cat_imp > 3)) error stop "polyreg support"
    print *, "test_categorical: PASS"
end program test_categorical
