! SPDX-License-Identifier: MIT
! Copyright (c) 2014-2022 Robert L. McDonald
module derivmkts_implied
    use derivmkts_kinds, only: dp
    use derivmkts_black_scholes, only: bscall, bsput
    implicit none
    private
    public :: bscallimpvol, bsputimpvol, bscallimps, bsputimps

contains

    real(dp) function bscallimpvol(s,k,r,tt,d,price,lowvol,highvol,ok) result(root)
        real(dp), intent(in) :: s,k,r,tt,d,price
        real(dp), intent(in), optional :: lowvol, highvol
        logical, intent(out), optional :: ok
        real(dp) :: lo, hi, mid, flo, fhi, fm
        integer :: iter
        lo = 1.0e-8_dp
        if (present(lowvol)) lo = max(lowvol,1.0e-12_dp)
        hi = 2.0_dp
        if (present(highvol)) hi = max(highvol,lo*2.0_dp)
        if (present(ok)) ok = .false.
        root = 0.0_dp
        if (price <= max(s*exp(-d*tt)-k*exp(-r*tt),0.0_dp) .or. price > s*exp(-d*tt)) return
        flo = bscall(s,k,lo,r,tt,d)-price
        fhi = bscall(s,k,hi,r,tt,d)-price
        do while (fhi < 0.0_dp .and. hi < 1.0e6_dp)
            hi = hi*2.0_dp
            fhi = bscall(s,k,hi,r,tt,d)-price
        end do
        if (flo*fhi > 0.0_dp) return
        do iter = 1, 200
            mid = 0.5_dp*(lo+hi)
            fm = bscall(s,k,mid,r,tt,d)-price
            if (abs(fm) < 1.0e-12_dp .or. hi-lo < 1.0e-12_dp*(1.0_dp+mid)) exit
            if (flo*fm <= 0.0_dp) then
                hi = mid
            else
                lo = mid
                flo = fm
            end if
        end do
        root = mid
        if (present(ok)) ok = .true.
    end function bscallimpvol

    real(dp) function bsputimpvol(s,k,r,tt,d,price,lowvol,highvol,ok) result(root)
        real(dp), intent(in) :: s,k,r,tt,d,price
        real(dp), intent(in), optional :: lowvol, highvol
        logical, intent(out), optional :: ok
        root = bscallimpvol(s,k,r,tt,d,price+s*exp(-d*tt)-k*exp(-r*tt),lowvol,highvol,ok)
    end function bsputimpvol

    real(dp) function bscallimps(s,k,v,r,tt,d,price,lower,upper,ok) result(root)
        real(dp), intent(in) :: s,k,v,r,tt,d,price
        real(dp), intent(in), optional :: lower,upper
        logical, intent(out), optional :: ok
        real(dp) :: lo,hi,mid,flo,fm
        integer :: iter
        lo = 1.0e-8_dp + 0.0_dp*s
        hi = max(1.0e4_dp, 10.0_dp*k)
        if (present(lower)) lo=lower
        if (present(upper)) hi=upper
        if (present(ok)) ok=.false.
        flo=bscall(lo,k,v,r,tt,d)-price
        do iter=1,200
            mid = 0.5_dp*(lo + hi)
            fm = bscall(mid,k,v,r,tt,d) - price
            if(abs(fm)<1.0e-12_dp .or. hi-lo<1.0e-12_dp*(1.0_dp+mid)) exit
            if (flo*fm <= 0.0_dp) then
                hi = mid
            else
                lo = mid
                flo = fm
            end if
        end do
        root=mid
        if(present(ok)) ok=.true.
    end function bscallimps

    real(dp) function bsputimps(s,k,v,r,tt,d,price,lower,upper,ok) result(root)
        real(dp), intent(in) :: s,k,v,r,tt,d,price
        real(dp), intent(in), optional :: lower,upper
        logical, intent(out), optional :: ok
        real(dp) :: lo,hi,mid,flo,fm
        integer :: iter
        lo = 1.0e-8_dp + 0.0_dp*s
        hi = max(1.0e4_dp, 10.0_dp*k)
        if (present(lower)) lo=lower
        if (present(upper)) hi=upper
        if (present(ok)) ok=.false.
        flo=bsput(lo,k,v,r,tt,d)-price
        do iter=1,200
            mid = 0.5_dp*(lo + hi)
            fm = bsput(mid,k,v,r,tt,d) - price
            if(abs(fm)<1.0e-12_dp .or. hi-lo<1.0e-12_dp*(1.0_dp+mid)) exit
            if (flo*fm <= 0.0_dp) then
                hi = mid
            else
                lo = mid
                flo = fm
            end if
        end do
        root=mid
        if(present(ok)) ok=.true.
    end function bsputimps

end module derivmkts_implied
