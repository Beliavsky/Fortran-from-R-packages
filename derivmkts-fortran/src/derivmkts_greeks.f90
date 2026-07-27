! SPDX-License-Identifier: MIT
! Copyright (c) 2014-2022 Robert L. McDonald
module derivmkts_greeks
    use derivmkts_kinds, only: dp
    use derivmkts_types, only: greek_result
    implicit none
    private
    public :: pricing_function, numerical_greeks, greeks, greeks2

    abstract interface
        pure real(kind(1.0d0)) function pricing_function(s,k,v,r,tt,d) result(price)
            real(kind(1.0d0)), intent(in) :: s,k,v,r,tt,d
        end function pricing_function
    end interface

    interface greeks
        module procedure numerical_greeks
    end interface greeks

    interface greeks2
        module procedure numerical_greeks
    end interface greeks2

contains

    pure function numerical_greeks(fn,s,k,v,r,tt,d,include_theta) result(g)
        procedure(pricing_function) :: fn
        real(dp), intent(in) :: s,k,v,r,tt,d
        logical, intent(in), optional :: include_theta
        type(greek_result) :: g
        real(dp), parameter :: e1=1.0e-4_dp,e2=5.0e-4_dp
        logical :: theta_ok
        real(dp) :: dup,ddn
        theta_ok=.true.
        if(present(include_theta))theta_ok=include_theta
        g%premium=fn(s,k,v,r,tt,d)
        g%delta=(fn(s+e1,k,v,r,tt,d)-fn(s-e1,k,v,r,tt,d))/(2.0_dp*e1)
        dup=(fn(s+e2+e1,k,v,r,tt,d)-fn(s+e2-e1,k,v,r,tt,d))/(2.0_dp*e1)
        ddn=(fn(s-e2+e1,k,v,r,tt,d)-fn(s-e2-e1,k,v,r,tt,d))/(2.0_dp*e1)
        g%gamma=(dup-ddn)/(2.0_dp*e2)
        g%vega=(fn(s,k,v+e1,r,tt,d)-fn(s,k,v-e1,r,tt,d))/(2.0_dp*e1)/100.0_dp
        g%rho=(fn(s,k,v,r+e1,tt,d)-fn(s,k,v,r-e1,tt,d))/(2.0_dp*e1)/100.0_dp
        if(theta_ok) then
            g%theta=-(fn(s,k,v,r,tt+e1,d)-fn(s,k,v,r,tt-e1,d))/(2.0_dp*e1)/365.0_dp
        end if
        g%psi=(fn(s,k,v,r,tt,d+e1)-fn(s,k,v,r,tt,d-e1))/(2.0_dp*e1)/100.0_dp
        if(abs(g%premium)>1.0e-6_dp)g%elasticity=s*g%delta/g%premium
    end function numerical_greeks

end module derivmkts_greeks
