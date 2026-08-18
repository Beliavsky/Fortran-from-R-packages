! SPDX-License-Identifier: GPL-3.0-or-later
module nbconv_math
    use nbconv_kinds, only : dp
    implicit none
    private

    real(dp), parameter, public :: nbconv_pi = acos(-1.0_dp)
    real(dp), parameter :: neg_huge = -huge(1.0_dp)

    public :: log_add_exp
    public :: log_sum_exp
    public :: negbin_logpmf_mu
    public :: negbin_pmf_mu

contains

    pure function log_add_exp(a, b) result(c)
        real(dp), intent(in) :: a, b
        real(dp) :: c, m

        if (a <= neg_huge / 2.0_dp) then
            c = b
        else if (b <= neg_huge / 2.0_dp) then
            c = a
        else
            m = max(a, b)
            c = m + log(exp(a - m) + exp(b - m))
        end if
    end function log_add_exp

    pure function log_sum_exp(x) result(ans)
        real(dp), intent(in) :: x(:)
        real(dp) :: ans, m

        if (size(x) == 0) then
            ans = neg_huge
            return
        end if
        m = maxval(x)
        if (m <= neg_huge / 2.0_dp) then
            ans = neg_huge
        else
            ans = m + log(sum(exp(x - m)))
        end if
    end function log_sum_exp

    pure function negbin_logpmf_mu(k, mu, phi) result(lp)
        integer, intent(in) :: k
        real(dp), intent(in) :: mu, phi
        real(dp) :: lp, p

        if (k < 0 .or. mu < 0.0_dp .or. phi <= 0.0_dp) then
            lp = neg_huge
            return
        end if
        if (mu <= 0.0_dp) then
            if (k == 0) then
                lp = 0.0_dp
            else
                lp = neg_huge
            end if
            return
        end if
        p = phi / (phi + mu)
        lp = log_gamma(real(k, dp) + phi) - log_gamma(phi) - log_gamma(real(k + 1, dp)) &
             + phi * log(p) + real(k, dp) * log1m(p)
    contains
        pure function log1m(x) result(y)
            real(dp), intent(in) :: x
            real(dp) :: y
            y = log(1.0_dp - x)
        end function log1m
    end function negbin_logpmf_mu

    pure function negbin_pmf_mu(k, mu, phi) result(pmf)
        integer, intent(in) :: k
        real(dp), intent(in) :: mu, phi
        real(dp) :: pmf, lp

        lp = negbin_logpmf_mu(k, mu, phi)
        if (lp <= neg_huge / 2.0_dp) then
            pmf = 0.0_dp
        else
            pmf = exp(lp)
        end if
    end function negbin_pmf_mu

end module nbconv_math
