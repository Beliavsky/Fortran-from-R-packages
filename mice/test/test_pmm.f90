program test_pmm
    use iso_fortran_env, only : int64
    use mice, only : dp, mice_ok, mice_rng_state, rng_seed, impute_pmm
    implicit none
    real(dp) :: y(8), x(8, 1)
    real(dp), allocatable :: imp(:)
    logical :: observed(8), where(8)
    type(mice_rng_state) :: rng
    integer :: i, info

    do i = 1, 8
        x(i, 1) = real(i, dp)
        y(i) = 10.0_dp * real(i, dp)
    end do
    observed = [.true., .true., .true., .true., .true., .false., .false., .false.]
    where = .not. observed
    call rng_seed(rng, 55_int64)
    call impute_pmm(y, observed, x, where, rng, imp, info, donors=3, matchtype=1)
    if (info /= mice_ok) error stop "pmm status"
    do i = 1, size(imp)
        if (.not. any(abs(imp(i) - y(1:5)) < 1.0e-13_dp)) error stop "PMM did not return donor value"
    end do
    print *, "test_pmm: PASS"
end program test_pmm
