program test_twolevel
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use iso_fortran_env, only : int64
    use mice, only : dp, mice_ok, mice_invalid_argument, mice_rng_state, rng_seed, impute_2lonly_mean, impute_2lonly_norm
    implicit none
    real(dp) :: y(6), x(6, 1)
    real(dp), allocatable :: imp(:)
    integer :: cluster(6), info
    logical :: observed(6), where(6)
    type(mice_rng_state) :: rng

    cluster = [1, 1, 2, 2, 3, 3]
    y = [5.0_dp, 5.0_dp, 7.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    x(:, 1) = [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp]
    observed = [.true., .true., .true., .false., .false., .false.]
    where = .not. observed
    call impute_2lonly_mean(y, observed, cluster, where, imp, info)
    if (info /= mice_ok) error stop "2l mean status"
    if (abs(imp(1) - 7.0_dp) > 1.0e-13_dp) error stop "2l mean partial repair"
    if (.not. ieee_is_nan(imp(2))) error stop "2l mean empty cluster"
    call rng_seed(rng, 9_int64)
    call impute_2lonly_norm(y, observed, cluster, x, where, rng, imp, info)
    if (info /= mice_invalid_argument) error stop "2l norm must reject partially missing cluster"
    observed = [.true., .true., .true., .true., .false., .false.]
    where = .not. observed
    call rng_seed(rng, 9_int64)
    call impute_2lonly_norm(y, observed, cluster, x, where, rng, imp, info)
    if (info /= mice_ok .or. size(imp) /= 2) error stop "2l norm systematic missingness"
    print *, "test_twolevel: PASS"
end program test_twolevel
