program test_2l_norm_full
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use iso_fortran_env, only : int64
    use mice, only : dp, mice_ok, mice_rng_state, rng_seed, mice_2l_norm_state, impute_2l_norm
    implicit none

    integer, parameter :: nclass = 5, per_class = 5, n = nclass * per_class
    real(dp) :: y(n)
    real(dp), allocatable :: x(:, :), imputed(:)
    integer :: cluster(n), i, g, info, row
    logical :: observed(n), where(n)
    type(mice_rng_state) :: rng
    type(mice_2l_norm_state) :: state

    allocate(x(n, 0))
    row = 0
    do g = 1, nclass
        do i = 1, per_class
            row = row + 1
            cluster(row) = g
            y(row) = 2.0_dp + 0.4_dp * real(g, dp) + 0.05_dp * real(i - 3, dp)
        end do
    end do
    observed = .true.
    do g = 1, nclass
        observed((g - 1) * per_class + per_class) = .false.
    end do
    where = .not. observed
    call rng_seed(rng, 606_int64)
    call impute_2l_norm(y, observed, cluster, x, where, rng, imputed, info, n_iter=40, state=state)
    if (info /= mice_ok) error stop "2l.norm status"
    if (size(imputed) /= nclass) error stop "2l.norm imputation count"
    if (.not. all(ieee_is_finite(imputed))) error stop "2l.norm finite imputations"
    if (any(state%residual_precision <= 0.0_dp)) error stop "2l.norm residual precision"
    if (state%sigma2_0 <= 0.0_dp .or. state%theta <= 0.0_dp) error stop "2l.norm variance hierarchy"
    if (size(state%random_effects, 1) /= nclass) error stop "2l.norm random effects shape"
    print *, "test_2l_norm_full: PASS"
end program test_2l_norm_full
