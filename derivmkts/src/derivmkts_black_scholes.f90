! SPDX-License-Identifier: MIT
! Copyright (c) 2014-2022 Robert L. McDonald
module derivmkts_black_scholes
    use derivmkts_kinds, only: dp
    use derivmkts_math, only: norm_cdf, norm_pdf
    use derivmkts_types, only: option_pair, greek_result
    implicit none
    private
    public :: bs_d1, bs_d2, bscall, bsput, assetcall, assetput, cashcall, cashput
    public :: bs_prices, bs_call_greeks, bs_put_greeks, bsopt

contains

    pure elemental real(dp) function bs_d1(s, k, v, r, tt, d) result(x)
        real(dp), intent(in) :: s, k, v, r, tt, d
        if (tt <= 0.0_dp .or. v <= 0.0_dp) then
            if (s*exp(-d*max(tt,0.0_dp)) > k*exp(-r*max(tt,0.0_dp))) then
                x = huge(1.0_dp)
            else
                x = -huge(1.0_dp)
            end if
        else
            x = (log(s/k) + (r-d+0.5_dp*v*v)*tt)/(v*sqrt(tt))
        end if
    end function bs_d1

    pure elemental real(dp) function bs_d2(s, k, v, r, tt, d) result(x)
        real(dp), intent(in) :: s, k, v, r, tt, d
        x = bs_d1(s,k,v,r,tt,d) - v*sqrt(max(tt,0.0_dp))
    end function bs_d2

    pure elemental real(dp) function bscall(s, k, v, r, tt, d) result(price)
        real(dp), intent(in) :: s, k, v, r, tt, d
        real(dp) :: sd, kd
        if (tt <= 0.0_dp) then
            price = max(s-k,0.0_dp)
        else if (v <= 0.0_dp) then
            sd = s*exp(-d*tt)
            kd = k*exp(-r*tt)
            price = max(sd-kd,0.0_dp)
        else
            price = s*exp(-d*tt)*norm_cdf(bs_d1(s,k,v,r,tt,d)) - &
                k*exp(-r*tt)*norm_cdf(bs_d2(s,k,v,r,tt,d))
        end if
    end function bscall

    pure elemental real(dp) function bsput(s, k, v, r, tt, d) result(price)
        real(dp), intent(in) :: s, k, v, r, tt, d
        price = bscall(s,k,v,r,tt,d) + k*exp(-r*max(tt,0.0_dp)) - &
            s*exp(-d*max(tt,0.0_dp))
        if (tt <= 0.0_dp) price = max(k-s,0.0_dp)
    end function bsput

    pure elemental real(dp) function assetcall(s, k, v, r, tt, d) result(price)
        real(dp), intent(in) :: s, k, v, r, tt, d
        if (tt <= 0.0_dp) then
            price = merge(s,0.0_dp,s>k)
        else
            price = s*exp(-d*tt)*norm_cdf(bs_d1(s,k,v,r,tt,d))
        end if
    end function assetcall

    pure elemental real(dp) function cashcall(s, k, v, r, tt, d) result(price)
        real(dp), intent(in) :: s, k, v, r, tt, d
        if (tt <= 0.0_dp) then
            price = merge(1.0_dp,0.0_dp,s>k)
        else
            price = exp(-r*tt)*norm_cdf(bs_d2(s,k,v,r,tt,d))
        end if
    end function cashcall

    pure elemental real(dp) function assetput(s, k, v, r, tt, d) result(price)
        real(dp), intent(in) :: s, k, v, r, tt, d
        if (tt <= 0.0_dp) then
            price = merge(s,0.0_dp,s<k)
        else
            price = s*exp(-d*tt)*norm_cdf(-bs_d1(s,k,v,r,tt,d))
        end if
    end function assetput

    pure elemental real(dp) function cashput(s, k, v, r, tt, d) result(price)
        real(dp), intent(in) :: s, k, v, r, tt, d
        if (tt <= 0.0_dp) then
            price = merge(1.0_dp,0.0_dp,s<k)
        else
            price = exp(-r*tt)*norm_cdf(-bs_d2(s,k,v,r,tt,d))
        end if
    end function cashput

    pure function bs_prices(s,k,v,r,tt,d) result(out)
        real(dp), intent(in) :: s,k,v,r,tt,d
        type(option_pair) :: out
        out%call = bscall(s,k,v,r,tt,d)
        out%put = bsput(s,k,v,r,tt,d)
    end function bs_prices

    pure function bs_call_greeks(s,k,v,r,tt,d) result(g)
        real(dp), intent(in) :: s,k,v,r,tt,d
        type(greek_result) :: g
        real(dp) :: x1, x2, discd, discr
        g%premium = bscall(s,k,v,r,tt,d)
        if (tt <= 0.0_dp .or. v <= 0.0_dp) return
        x1 = bs_d1(s,k,v,r,tt,d)
        x2 = x1-v*sqrt(tt)
        discd = exp(-d*tt)
        discr = exp(-r*tt)
        g%delta = discd*norm_cdf(x1)
        g%gamma = discd*norm_pdf(x1)/(s*v*sqrt(tt))
        g%vega = s*discd*norm_pdf(x1)*sqrt(tt)/100.0_dp
        g%rho = k*tt*discr*norm_cdf(x2)/100.0_dp
        g%theta = (-s*discd*norm_pdf(x1)*v/(2.0_dp*sqrt(tt)) - &
            r*k*discr*norm_cdf(x2) + d*s*discd*norm_cdf(x1))/365.0_dp
        g%psi = -s*tt*discd*norm_cdf(x1)/100.0_dp
        if (abs(g%premium) > 1.0e-14_dp) g%elasticity = s*g%delta/g%premium
    end function bs_call_greeks

    pure function bs_put_greeks(s,k,v,r,tt,d) result(g)
        real(dp), intent(in) :: s,k,v,r,tt,d
        type(greek_result) :: g
        real(dp) :: x1, x2, discd, discr
        g%premium = bsput(s,k,v,r,tt,d)
        if (tt <= 0.0_dp .or. v <= 0.0_dp) return
        x1 = bs_d1(s,k,v,r,tt,d)
        x2 = x1-v*sqrt(tt)
        discd = exp(-d*tt)
        discr = exp(-r*tt)
        g%delta = -discd*norm_cdf(-x1)
        g%gamma = discd*norm_pdf(x1)/(s*v*sqrt(tt))
        g%vega = s*discd*norm_pdf(x1)*sqrt(tt)/100.0_dp
        g%rho = -k*tt*discr*norm_cdf(-x2)/100.0_dp
        g%theta = (-s*discd*norm_pdf(x1)*v/(2.0_dp*sqrt(tt)) + &
            r*k*discr*norm_cdf(-x2) - d*s*discd*norm_cdf(-x1))/365.0_dp
        g%psi = s*tt*discd*norm_cdf(-x1)/100.0_dp
        if (abs(g%premium) > 1.0e-14_dp) g%elasticity = s*g%delta/g%premium
    end function bs_put_greeks

    pure subroutine bsopt(s,k,v,r,tt,d,call_result,put_result)
        real(dp), intent(in) :: s,k,v,r,tt,d
        type(greek_result), intent(out) :: call_result, put_result
        call_result = bs_call_greeks(s,k,v,r,tt,d)
        put_result = bs_put_greeks(s,k,v,r,tt,d)
    end subroutine bsopt

end module derivmkts_black_scholes
