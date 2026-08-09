! SPDX-License-Identifier: GPL-2.0-only
module nls2_stats
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
    use nls2_kinds, only : dp
    use nls2_types, only : nls_result
    implicit none
    private
    public :: nls_loglik, nls_df_residual, pearson_residuals

contains

    real(dp) function nls_loglik(result, weights) result(ll)
        type(nls_result), intent(in) :: result
        real(dp), intent(in), optional :: weights(:)
        integer :: n, i
        real(dp) :: sumlogw, pi
        logical, allocatable :: use(:)

        pi = acos(-1.0_dp)
        if (present(weights)) then
            if (size(weights) /= size(result%residuals)) then
                ll = -huge(1.0_dp)
                return
            end if
            allocate(use(size(weights)))
            use = weights > 0.0_dp
            n = count(use)
            sumlogw = 0.0_dp
            do i = 1, size(weights)
                if (use(i)) sumlogw = sumlogw + log(weights(i))
            end do
        else
            n = size(result%residuals)
            sumlogw = 0.0_dp
        end if
        if (n <= 0 .or. result%rss < 0.0_dp) then
            ll = -huge(1.0_dp)
            return
        end if
        if (result%rss <= tiny(1.0_dp)) then
            ll = ieee_value(1.0_dp, ieee_positive_inf)
            return
        end if
        ll = -0.5_dp * real(n,dp) * (log(2.0_dp*pi) + 1.0_dp - log(real(n,dp)) &
            - sumlogw / real(n,dp) + log(result%rss))
    end function nls_loglik

    integer function nls_df_residual(result, weights) result(df)
        type(nls_result), intent(in) :: result
        real(dp), intent(in), optional :: weights(:)
        integer :: n, p
        p = size(result%par)
        if (allocated(result%linear_par)) p = p + size(result%linear_par)
        if (present(weights)) then
            n = count(weights > 0.0_dp)
        else
            n = size(result%residuals)
        end if
        df = n - p
    end function nls_df_residual

    subroutine pearson_residuals(result, r)
        type(nls_result), intent(in) :: result
        real(dp), intent(out) :: r(:)
        if (size(r) /= size(result%residuals) .or. result%sigma <= 0.0_dp) then
            r = 0.0_dp
        else
            r = result%residuals / result%sigma
        end if
    end subroutine pearson_residuals

end module nls2_stats
