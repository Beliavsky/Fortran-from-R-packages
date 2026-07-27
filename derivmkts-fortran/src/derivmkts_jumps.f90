! SPDX-License-Identifier: MIT
! Copyright (c) 2014-2022 Robert L. McDonald
module derivmkts_jumps
    use derivmkts_kinds, only: dp
    use derivmkts_black_scholes, only: bscall, bsput, assetcall, assetput, cashcall, cashput
    use derivmkts_types, only: option_pair
    implicit none
    private
    public :: mertonjump, assetjump, cashjump
contains
    function mertonjump(s,k,v,r,tt,d,lambda,alphaj,vj) result(out)
        real(dp), intent(in) :: s,k,v,r,tt,d,lambda,alphaj,vj
        type(option_pair) :: out
        out%call = jump_price(s,k,v,r,tt,d,lambda,alphaj,vj,1)
        out%put = jump_price(s,k,v,r,tt,d,lambda,alphaj,vj,2)
    end function mertonjump
    function assetjump(s,k,v,r,tt,d,lambda,alphaj,vj) result(out)
        real(dp), intent(in) :: s,k,v,r,tt,d,lambda,alphaj,vj
        type(option_pair) :: out
        out%call = jump_price(s,k,v,r,tt,d,lambda,alphaj,vj,3)
        out%put = jump_price(s,k,v,r,tt,d,lambda,alphaj,vj,4)
    end function assetjump
    function cashjump(s,k,v,r,tt,d,lambda,alphaj,vj) result(out)
        real(dp), intent(in) :: s,k,v,r,tt,d,lambda,alphaj,vj
        type(option_pair) :: out
        out%call = jump_price(s,k,v,r,tt,d,lambda,alphaj,vj,5)
        out%put = jump_price(s,k,v,r,tt,d,lambda,alphaj,vj,6)
    end function cashjump

    real(dp) function jump_price(s,k,v,r,tt,d,lambda,alphaj,vj,which) result(price)
        real(dp), intent(in) :: s,k,v,r,tt,d,lambda,alphaj,vj
        integer, intent(in) :: which
        real(dp) :: lp,mu,p,cum,vp,rp,kappa,term
        integer :: i
        lp = lambda*exp(alphaj)
        mu = lp*tt
        if (mu < 1.0e-14_dp) then
            price = base_price(s,k,v,r,tt,d,which)
            return
        end if
        p = exp(-mu)
        cum = 0.0_dp
        price = 0.0_dp
        kappa = exp(alphaj)-1.0_dp
        i = 0
        do
            vp = sqrt(v*v+real(i,dp)*vj*vj/tt)
            rp = r-lambda*kappa+real(i,dp)*alphaj/tt
            term = p*base_price(s,k,vp,rp,tt,d,which)
            price = price+term
            cum = cum+p
            if (cum >= 1.0_dp-1.0e-12_dp) exit
            i = i+1
            p = p*mu/real(i,dp)
            if (i > 100000) exit
        end do
    end function jump_price

    pure real(dp) function base_price(s,k,v,r,tt,d,which) result(price)
        real(dp), intent(in) :: s,k,v,r,tt,d
        integer, intent(in) :: which
        select case(which)
        case(1)
            price = bscall(s,k,v,r,tt,d)
        case(2)
            price = bsput(s,k,v,r,tt,d)
        case(3)
            price = assetcall(s,k,v,r,tt,d)
        case(4)
            price = assetput(s,k,v,r,tt,d)
        case(5)
            price = cashcall(s,k,v,r,tt,d)
        case default
            price = cashput(s,k,v,r,tt,d)
        end select
    end function base_price
end module derivmkts_jumps
