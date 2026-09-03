program test_native
    use iso_fortran_env, only : int64
    use mice, only : dp, mice_ok, mice_rng_state, rng_seed, legendre_basis, matchindex
    implicit none
    real(dp), allocatable :: basis(:, :)
    integer, allocatable :: idx(:)
    type(mice_rng_state) :: rng
    real(dp) :: x(3), donors(3), target(3)
    integer :: info

    x = [0.0_dp, 0.5_dp, 1.0_dp]
    call legendre_basis(x, 3, basis, info)
    if (info /= mice_ok) error stop "legendre status"
    if (abs(basis(2, 1)) > 1.0e-13_dp) error stop "legendre P1"
    if (abs(basis(2, 2) + sqrt(5.0_dp) / 2.0_dp) > 1.0e-13_dp) error stop "legendre P2"
    if (abs(basis(2, 3)) > 1.0e-13_dp) error stop "legendre P3"
    donors = [-5.0_dp, 0.0_dp, 10.0_dp]
    target = [-4.9_dp, 0.1_dp, 9.8_dp]
    call rng_seed(rng, 12345_int64)
    call matchindex(donors, target, 1, rng, idx, info)
    if (info /= mice_ok) error stop "matchindex status"
    if (any(idx /= [1, 2, 3])) error stop "nearest donor mismatch"
    print *, "test_native: PASS"
end program test_native
