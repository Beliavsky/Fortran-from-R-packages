! SPDX-License-Identifier: MIT
! Copyright (c) 2014-2022 Robert L. McDonald
module derivmkts_barriers
    use derivmkts_kinds, only: dp
    use derivmkts_math, only: norm_cdf
    use derivmkts_black_scholes, only: bscall, bsput, assetcall, assetput, cashcall, cashput
    implicit none
    private
    public :: cashdicall, assetdicall, cashdocall, assetdocall
    public :: cashdiput, assetdiput, cashdoput, assetdoput
    public :: cashuicall, assetuicall, cashuocall, assetuocall
    public :: cashuiput, assetuiput, cashuoput, assetuoput
    public :: calldownin, calldownout, putdownin, putdownout
    public :: callupin, callupout, putupin, putupout
    public :: dicall, docall, diput, doput, uicall, uocall, uiput, uoput
    public :: dr, ur, drdeferred, urdeferred
contains

    pure real(dp) function d5(s,v,r,tt,d,h) result(x)
        real(dp), intent(in) :: s,v,r,tt,d,h
        x = (log(s/h)+(r-d+0.5_dp*v*v)*tt)/(v*sqrt(tt))
    end function d5

    pure real(dp) function d6(s,v,r,tt,d,h) result(x)
        real(dp), intent(in) :: s,v,r,tt,d,h
        x = d5(s,v,r,tt,d,h)-v*sqrt(tt)
    end function d6

    pure real(dp) function d7(s,v,r,tt,d,h) result(x)
        real(dp), intent(in) :: s,v,r,tt,d,h
        x = (log(h/s)+(r-d+0.5_dp*v*v)*tt)/(v*sqrt(tt))
    end function d7

    pure real(dp) function d8(s,v,r,tt,d,h) result(x)
        real(dp), intent(in) :: s,v,r,tt,d,h
        x = d7(s,v,r,tt,d,h)-v*sqrt(tt)
    end function d8

    pure real(dp) function cashdicall(s,k,v,r,tt,d,h) result(price)
        real(dp), intent(in) :: s,k,v,r,tt,d,h
        real(dp) :: kh
        if (s <= h) then
            price = cashcall(s,k,v,r,tt,d)
            return
        end if
        kh = max(k,h)
        price = exp(-r*tt)*(norm_cdf((log(s/k)+(r-d-0.5_dp*v*v)*tt)/(v*sqrt(tt))) - &
            norm_cdf(d6(s,v,r,tt,d,kh)) + (h/s)**(2.0_dp*(r-d)/(v*v)-1.0_dp)* &
            norm_cdf(d8(s,v,r,tt,d,h*min(h/k,1.0_dp))))
    end function cashdicall

    pure real(dp) function assetdicall(s,k,v,r,tt,d,h) result(price)
        real(dp), intent(in) :: s,k,v,r,tt,d,h
        price = exp((r-d)*tt)*s*cashdicall(s,k,v,r,tt,d-v*v,h)
    end function assetdicall

    pure real(dp) function cashdocall(s,k,v,r,tt,d,h) result(price)
        real(dp), intent(in) :: s,k,v,r,tt,d,h
        price = cashcall(s,k,v,r,tt,d)-cashdicall(s,k,v,r,tt,d,h)
    end function cashdocall

    pure real(dp) function assetdocall(s,k,v,r,tt,d,h) result(price)
        real(dp), intent(in) :: s,k,v,r,tt,d,h
        price = s*exp((r-d)*tt)*cashdocall(s,k,v,r,tt,d-v*v,h)
    end function assetdocall

    pure real(dp) function cashdoput(s,k,v,r,tt,d,h) result(price)
        real(dp), intent(in) :: s,k,v,r,tt,d,h
        if (s <= h .or. k <= h) then
            price = 0.0_dp
        else
            price = cashdocall(s,h,v,r,tt,d,h)-cashdocall(s,k,v,r,tt,d,h)
        end if
    end function cashdoput

    pure real(dp) function cashdiput(s,k,v,r,tt,d,h) result(price)
        real(dp), intent(in) :: s,k,v,r,tt,d,h
        price = cashput(s,k,v,r,tt,d)-cashdoput(s,k,v,r,tt,d,h)
    end function cashdiput

    pure real(dp) function assetdoput(s,k,v,r,tt,d,h) result(price)
        real(dp), intent(in) :: s,k,v,r,tt,d,h
        price = s*exp((r-d)*tt)*cashdoput(s,k,v,r,tt,d-v*v,h)
    end function assetdoput

    pure real(dp) function assetdiput(s,k,v,r,tt,d,h) result(price)
        real(dp), intent(in) :: s,k,v,r,tt,d,h
        price = s*exp((r-d)*tt)*cashdiput(s,k,v,r,tt,d-v*v,h)
    end function assetdiput

    pure real(dp) function calldownin(s,k,v,r,tt,d,h) result(price)
        real(dp), intent(in) :: s,k,v,r,tt,d,h
        price = assetdicall(s,k,v,r,tt,d,h)-k*cashdicall(s,k,v,r,tt,d,h)
    end function calldownin

    pure real(dp) function calldownout(s,k,v,r,tt,d,h) result(price)
        real(dp), intent(in) :: s,k,v,r,tt,d,h
        price = bscall(s,k,v,r,tt,d)-calldownin(s,k,v,r,tt,d,h)
    end function calldownout

    pure real(dp) function putdownin(s,k,v,r,tt,d,h) result(price)
        real(dp), intent(in) :: s,k,v,r,tt,d,h
        price = k*cashdiput(s,k,v,r,tt,d,h)-assetdiput(s,k,v,r,tt,d,h)
    end function putdownin

    pure real(dp) function putdownout(s,k,v,r,tt,d,h) result(price)
        real(dp), intent(in) :: s,k,v,r,tt,d,h
        price = bsput(s,k,v,r,tt,d)-putdownin(s,k,v,r,tt,d,h)
    end function putdownout

    pure real(dp) function cashuiput(s,k,v,r,tt,d,h) result(price)
        real(dp), intent(in) :: s,k,v,r,tt,d,h
        real(dp) :: kh
        if (s >= h) then
            price = cashput(s,k,v,r,tt,d)
            return
        end if
        kh = min(k,h)
        price = exp(-r*tt)*(1.0_dp - norm_cdf((log(s/k)+(r-d-0.5_dp*v*v)*tt)/(v*sqrt(tt))) - &
            (1.0_dp-norm_cdf(d6(s,v,r,tt,d,kh))) + &
            (h/s)**(2.0_dp*(r-d)/(v*v)-1.0_dp)*(1.0_dp - &
            norm_cdf(d8(s,v,r,tt,d,h*max(h/k,1.0_dp)))))
    end function cashuiput

    pure real(dp) function cashuoput(s,k,v,r,tt,d,h) result(price)
        real(dp), intent(in) :: s,k,v,r,tt,d,h
        price = cashput(s,k,v,r,tt,d)-cashuiput(s,k,v,r,tt,d,h)
    end function cashuoput

    pure real(dp) function cashuicall(s,k,v,r,tt,d,h) result(price)
        real(dp), intent(in) :: s,k,v,r,tt,d,h
        price = cashuiput(s,1.0e15_dp,v,r,tt,d,h)-cashuiput(s,k,v,r,tt,d,h)
    end function cashuicall

    pure real(dp) function cashuocall(s,k,v,r,tt,d,h) result(price)
        real(dp), intent(in) :: s,k,v,r,tt,d,h
        price = cashcall(s,k,v,r,tt,d)-cashuicall(s,k,v,r,tt,d,h)
    end function cashuocall

    pure real(dp) function assetuiput(s,k,v,r,tt,d,h) result(price)
        real(dp), intent(in) :: s,k,v,r,tt,d,h
        price = s*exp((r-d)*tt)*cashuiput(s,k,v,r,tt,d-v*v,h)
    end function assetuiput

    pure real(dp) function assetuoput(s,k,v,r,tt,d,h) result(price)
        real(dp), intent(in) :: s,k,v,r,tt,d,h
        price = s*exp((r-d)*tt)*cashuoput(s,k,v,r,tt,d-v*v,h)
    end function assetuoput

    pure real(dp) function assetuicall(s,k,v,r,tt,d,h) result(price)
        real(dp), intent(in) :: s,k,v,r,tt,d,h
        price = s*exp((r-d)*tt)*cashuicall(s,k,v,r,tt,d-v*v,h)
    end function assetuicall

    pure real(dp) function assetuocall(s,k,v,r,tt,d,h) result(price)
        real(dp), intent(in) :: s,k,v,r,tt,d,h
        price = s*exp((r-d)*tt)*cashuocall(s,k,v,r,tt,d-v*v,h)
    end function assetuocall

    pure real(dp) function callupin(s,k,v,r,tt,d,h) result(price)
        real(dp), intent(in) :: s,k,v,r,tt,d,h
        price = assetuicall(s,k,v,r,tt,d,h)-k*cashuicall(s,k,v,r,tt,d,h)
    end function callupin

    pure real(dp) function callupout(s,k,v,r,tt,d,h) result(price)
        real(dp), intent(in) :: s,k,v,r,tt,d,h
        price = assetuocall(s,k,v,r,tt,d,h)-k*cashuocall(s,k,v,r,tt,d,h)
    end function callupout

    pure real(dp) function putupin(s,k,v,r,tt,d,h) result(price)
        real(dp), intent(in) :: s,k,v,r,tt,d,h
        price = k*cashuiput(s,k,v,r,tt,d,h)-assetuiput(s,k,v,r,tt,d,h)
    end function putupin

    pure real(dp) function putupout(s,k,v,r,tt,d,h) result(price)
        real(dp), intent(in) :: s,k,v,r,tt,d,h
        price = k*cashuoput(s,k,v,r,tt,d,h)-assetuoput(s,k,v,r,tt,d,h)
    end function putupout

    pure real(dp) function dicall(s,k,v,r,tt,d,h) result(price)
        real(dp), intent(in) :: s,k,v,r,tt,d,h
        price = calldownin(s,k,v,r,tt,d,h)
    end function dicall
    pure real(dp) function docall(s,k,v,r,tt,d,h) result(price)
        real(dp), intent(in) :: s,k,v,r,tt,d,h
        price = calldownout(s,k,v,r,tt,d,h)
    end function docall
    pure real(dp) function diput(s,k,v,r,tt,d,h) result(price)
        real(dp), intent(in) :: s,k,v,r,tt,d,h
        price = putdownin(s,k,v,r,tt,d,h)
    end function diput
    pure real(dp) function doput(s,k,v,r,tt,d,h) result(price)
        real(dp), intent(in) :: s,k,v,r,tt,d,h
        price = putdownout(s,k,v,r,tt,d,h)
    end function doput
    pure real(dp) function uicall(s,k,v,r,tt,d,h) result(price)
        real(dp), intent(in) :: s,k,v,r,tt,d,h
        price = callupin(s,k,v,r,tt,d,h)
    end function uicall
    pure real(dp) function uocall(s,k,v,r,tt,d,h) result(price)
        real(dp), intent(in) :: s,k,v,r,tt,d,h
        price = callupout(s,k,v,r,tt,d,h)
    end function uocall
    pure real(dp) function uiput(s,k,v,r,tt,d,h) result(price)
        real(dp), intent(in) :: s,k,v,r,tt,d,h
        price = putupin(s,k,v,r,tt,d,h)
    end function uiput
    pure real(dp) function uoput(s,k,v,r,tt,d,h) result(price)
        real(dp), intent(in) :: s,k,v,r,tt,d,h
        price = putupout(s,k,v,r,tt,d,h)
    end function uoput

    pure real(dp) function ur(s,v,r,tt,d,h,perpetual) result(val)
        real(dp), intent(in) :: s,v,r,tt,d,h
        logical, intent(in), optional :: perpetual
        logical :: perp
        real(dp) :: g,h1,h2,z1,z2
        perp = .false.
        if (present(perpetual)) perp = perpetual
        if (s >= h) then
            val = 1.0_dp
            return
        end if
        g = sqrt(((r-d)/(v*v)-0.5_dp)**2+2.0_dp*r/(v*v))
        h1 = 0.5_dp-(r-d)/(v*v)+g
        h2 = 0.5_dp-(r-d)/(v*v)-g
        if (perp) then
            val = (s/h)**h1
        else
            z1 = (log(h/s)-g*v*v*tt)/(v*sqrt(tt))
            z2 = (log(h/s)+g*v*v*tt)/(v*sqrt(tt))
            val = (s/h)**h1*norm_cdf(-z1)+(s/h)**h2*norm_cdf(-z2)
        end if
    end function ur

    pure real(dp) function dr(s,v,r,tt,d,h,perpetual) result(val)
        real(dp), intent(in) :: s,v,r,tt,d,h
        logical, intent(in), optional :: perpetual
        logical :: perp
        real(dp) :: g,h1,h2,z1,z2
        perp = .false.
        if (present(perpetual)) perp = perpetual
        if (s <= h) then
            val = 1.0_dp
            return
        end if
        g = sqrt(((r-d)/(v*v)-0.5_dp)**2+2.0_dp*r/(v*v))
        h1 = 0.5_dp-(r-d)/(v*v)+g
        h2 = 0.5_dp-(r-d)/(v*v)-g
        if (perp) then
            val = (s/h)**h2
        else
            z1 = (log(h/s)-g*v*v*tt)/(v*sqrt(tt))
            z2 = (log(h/s)+g*v*v*tt)/(v*sqrt(tt))
            val = (s/h)**h1*norm_cdf(z1)+(s/h)**h2*norm_cdf(z2)
        end if
    end function dr

    pure real(dp) function drdeferred(s,v,r,tt,d,h) result(val)
        real(dp), intent(in) :: s,v,r,tt,d,h
        val = cashdicall(s,1.0e-8_dp,v,r,tt,d,h)
    end function drdeferred

    pure real(dp) function urdeferred(s,v,r,tt,d,h) result(val)
        real(dp), intent(in) :: s,v,r,tt,d,h
        val = cashuicall(s,1.0e-8_dp,v,r,tt,d,h)
    end function urdeferred

end module derivmkts_barriers
