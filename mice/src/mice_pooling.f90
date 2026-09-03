! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from R package mice 3.19.0 by Stef van Buuren, Karin Groothuis-Oudshoorn,
! and mice contributors; see NOTICE.md and PROVENANCE.md for attribution.
! Computational translation derived from mice 3.19.0 pooling and D3 routines.
module mice_pooling
    use r_kinds, only : dp
    use r_linalg, only : inverse_matrix
    use mice_status, only : mice_ok, mice_invalid_argument, mice_invalid_shape, mice_singular
    implicit none
    private

    type, public :: pool_scalar_result
        integer :: m = 0
        real(dp) :: qbar = 0.0_dp
        real(dp) :: ubar = 0.0_dp
        real(dp) :: b = 0.0_dp
        real(dp) :: total = 0.0_dp
        real(dp) :: df = 0.0_dp
        real(dp) :: r = 0.0_dp
        real(dp) :: fmi = 0.0_dp
    end type pool_scalar_result

    type, public :: pool_vector_result
        integer :: m = 0
        real(dp), allocatable :: qbar(:)
        real(dp), allocatable :: ubar(:, :)
        real(dp), allocatable :: between(:, :)
        real(dp), allocatable :: total(:, :)
    end type pool_vector_result

    type, public :: d3_result
        real(dp) :: statistic = 0.0_dp
        integer :: numerator_df = 0
        real(dp) :: denominator_df = 0.0_dp
        real(dp) :: relative_increase = 0.0_dp
    end type d3_result

    public :: barnard_rubin
    public :: pool_scalar
    public :: pool_vector
    public :: pooled_wald
    public :: d3_from_deviances

contains

    pure real(dp) function barnard_rubin(m, b, total, df_complete) result(df)
        integer, intent(in), value :: m !! Number of imputations; must exceed one for finite between-imputation degrees of freedom.
        real(dp), intent(in), value :: b !! Between-imputation variance.
        real(dp), intent(in), value :: total !! Total Rubin variance.
        real(dp), intent(in), value :: df_complete !! Complete-data degrees of freedom; a huge value represents infinity.

        real(dp) :: df_old, lambda, tmp

        if (m <= 1 .or. total <= 0.0_dp) then
            df = huge(1.0_dp)
            return
        end if
        lambda = (1.0_dp + 1.0_dp / real(m, dp)) * b / total
        if (lambda <= epsilon(1.0_dp)) then
            df = huge(1.0_dp)
            return
        end if
        df_old = real(m - 1, dp) / (lambda * lambda)
        if (df_complete >= 0.5_dp * huge(1.0_dp)) then
            df = df_old
            return
        end if
        tmp = (1.0_dp - lambda) * (1.0_dp + df_complete) * df_complete
        df = real(m - 1, dp) * tmp / (real((m - 1), dp) * (df_complete + 3.0_dp) + lambda * lambda * tmp)
    end function barnard_rubin

    pure subroutine pool_scalar(q, u, result, info, n, k, reiter)
        real(dp), intent(in) :: q(:) !! Scalar complete-data estimates, one per imputation.
        real(dp), intent(in) :: u(:) !! Complete-data variances corresponding one-to-one with `q`.
        type(pool_scalar_result), intent(out) :: result !! Rubin or Reiter pooled scalar summary.
        integer, intent(out) :: info !! `mice_ok` on success or an argument/shape status code.
        real(dp), intent(in), optional :: n !! Complete-data sample size; omitted means effectively infinite.
        integer, intent(in), optional :: k !! Number of fitted complete-data parameters; default one.
        logical, intent(in), optional :: reiter !! Use Reiter 2003 partially-synthetic pooling instead of Rubin when true.

        real(dp) :: df_complete, denom
        integer :: npar
        logical :: synthetic

        if (size(q) /= size(u) .or. size(q) < 2) then
            info = mice_invalid_shape
            return
        end if
        if (any(u < 0.0_dp)) then
            info = mice_invalid_argument
            return
        end if
        result%m = size(q)
        result%qbar = sum(q) / real(result%m, dp)
        result%ubar = sum(u) / real(result%m, dp)
        result%b = sum((q - result%qbar)**2) / real(result%m - 1, dp)
        if (result%ubar > 0.0_dp) then
            result%r = (1.0_dp + 1.0_dp / real(result%m, dp)) * result%b / result%ubar
        else if (result%b > 0.0_dp) then
            result%r = huge(1.0_dp)
        else
            result%r = 0.0_dp
        end if
        synthetic = .false.
        if (present(reiter)) synthetic = reiter
        if (synthetic) then
            result%total = result%ubar + result%b / real(result%m, dp)
            denom = result%b / real(result%m, dp)
            if (denom > 0.0_dp) then
                result%df = real(result%m - 1, dp) * (1.0_dp + result%ubar / denom)**2
            else
                result%df = huge(1.0_dp)
            end if
            result%fmi = -1.0_dp
        else
            result%total = result%ubar + real(result%m + 1, dp) * result%b / real(result%m, dp)
            df_complete = huge(1.0_dp)
            npar = 1
            if (present(k)) npar = k
            if (present(n)) df_complete = max(n - real(npar, dp), 1.0_dp)
            result%df = barnard_rubin(result%m, result%b, result%total, df_complete)
            if (result%r >= 0.5_dp * huge(1.0_dp)) then
                result%fmi = 1.0_dp
            else
                result%fmi = (result%r + 2.0_dp / (result%df + 3.0_dp)) / (result%r + 1.0_dp)
            end if
        end if
        info = mice_ok
    end subroutine pool_scalar

    pure subroutine pool_vector(q, covariance, result, info)
        real(dp), intent(in) :: q(:, :) !! Complete-data estimates with parameters in rows and imputations in columns.
        real(dp), intent(in) :: covariance(:, :, :) !! Complete-data covariance matrices, one slice per imputation.
        type(pool_vector_result), intent(out) :: result !! Multivariate Rubin pooled means and covariance components.
        integer, intent(out) :: info !! `mice_ok` on success or a shape status code.

        real(dp), allocatable :: delta(:)
        integer :: i, m, p

        p = size(q, 1)
        m = size(q, 2)
        if (m < 2 .or. size(covariance, 1) /= p .or. size(covariance, 2) /= p .or. size(covariance, 3) /= m) then
            info = mice_invalid_shape
            return
        end if
        result%m = m
        allocate(result%qbar(p), result%ubar(p, p), result%between(p, p), result%total(p, p), delta(p))
        result%qbar = sum(q, dim=2) / real(m, dp)
        result%ubar = sum(covariance, dim=3) / real(m, dp)
        result%between = 0.0_dp
        do i = 1, m
            delta = q(:, i) - result%qbar
            result%between = result%between + spread(delta, 2, p) * spread(delta, 1, p)
        end do
        result%between = result%between / real(m - 1, dp)
        result%total = result%ubar + (1.0_dp + 1.0_dp / real(m, dp)) * result%between
        info = mice_ok
    end subroutine pool_vector

    subroutine pooled_wald(estimate, null_value, covariance, statistic, info)
        real(dp), intent(in) :: estimate(:) !! Pooled parameter estimate being tested.
        real(dp), intent(in) :: null_value(:) !! Null-hypothesis parameter vector of the same length as `estimate`.
        real(dp), intent(in) :: covariance(:, :) !! Positive-definite pooled covariance matrix.
        real(dp), intent(out) :: statistic !! Wald chi-square quadratic form `(estimate-null)' T^-1 (estimate-null)`.
        integer, intent(out) :: info !! `mice_ok` on success or a shape/singularity status code.

        real(dp), allocatable :: inverse(:, :), delta(:)
        integer :: la_info, p

        p = size(estimate)
        if (size(null_value) /= p .or. size(covariance, 1) /= p .or. size(covariance, 2) /= p) then
            info = mice_invalid_shape
            statistic = 0.0_dp
            return
        end if
        call inverse_matrix(covariance, inverse, la_info)
        if (la_info /= 0) then
            info = mice_singular
            statistic = 0.0_dp
            return
        end if
        allocate(delta(p))
        delta = estimate - null_value
        statistic = dot_product(delta, matmul(inverse, delta))
        info = mice_ok
    end subroutine pooled_wald

    pure subroutine d3_from_deviances(dev_full_fit, dev_null_fit, dev_full_restricted, dev_null_restricted, &
                                     n_parameters, result, info)
        real(dp), intent(in) :: dev_full_fit(:) !! Deviances of the freely fitted full model for each imputation.
        real(dp), intent(in) :: dev_null_fit(:) !! Deviances of the freely fitted null model for each imputation.
        real(dp), intent(in) :: dev_full_restricted(:) !! Full-model deviances with coefficients restricted to pooled estimates.
        real(dp), intent(in) :: dev_null_restricted(:) !! Null-model deviances with coefficients restricted to pooled estimates.
        integer, intent(in), value :: n_parameters !! Difference in dimensionality between full and null models.
        type(d3_result), intent(out) :: result !! Meng-Rubin D3 statistic, degrees of freedom, and relative increase.
        integer, intent(out) :: info !! `mice_ok` on success or an argument/shape status code.

        real(dp) :: dev_fit, dev_restricted, rm, v
        integer :: m

        m = size(dev_full_fit)
        if (m < 2 .or. size(dev_null_fit) /= m .or. size(dev_full_restricted) /= m .or. &
            size(dev_null_restricted) /= m) then
            info = mice_invalid_shape
            return
        end if
        if (n_parameters < 1) then
            info = mice_invalid_argument
            return
        end if
        dev_fit = sum(dev_null_fit - dev_full_fit) / real(m, dp)
        dev_restricted = sum(dev_null_restricted - dev_full_restricted) / real(m, dp)
        rm = real(m + 1, dp) * (dev_fit - dev_restricted) / &
            (real(n_parameters, dp) * real(m - 1, dp))
        rm = max(rm, epsilon(1.0_dp))
        result%statistic = dev_restricted / (real(n_parameters, dp) * (1.0_dp + rm))
        result%numerator_df = n_parameters
        result%relative_increase = rm
        v = real(n_parameters * (m - 1), dp)
        if (v > 4.0_dp) then
            result%denominator_df = 4.0_dp + (v - 4.0_dp) * &
                (1.0_dp + (1.0_dp - 2.0_dp / v) / rm)**2
        else
            result%denominator_df = 0.5_dp * v * (1.0_dp + 1.0_dp / real(n_parameters, dp)) * &
                (1.0_dp + 1.0_dp / rm)**2
        end if
        info = mice_ok
    end subroutine d3_from_deviances

end module mice_pooling
