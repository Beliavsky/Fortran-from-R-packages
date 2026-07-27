! SPDX-License-Identifier: MIT
! Copyright (c) 2014-2022 Robert L. McDonald
module derivmkts_compound
    use derivmkts_kinds, only: dp
    use derivmkts_math, only: norm_cdf, bivar_norm_cdf
    use derivmkts_black_scholes, only: bscall, bsput, bs_d1
    use derivmkts_implied, only: bscallimps, bsputimps
    use derivmkts_types, only: compound_result
    implicit none
    private
    public :: binormsdist, calloncall, putoncall, callonput, putonput
    public :: optionsoncall, optionsonput
contains
    real(dp) function binormsdist(x1,x2,rho) result(p)
        real(dp), intent(in) :: x1,x2,rho
        p = bivar_norm_cdf(x1,x2,rho)
    end function binormsdist

    function calloncall(s,kuo,kco,v,r,t1,t2,d) result(out)
        real(dp), intent(in) :: s,kuo,kco,v,r,t1,t2,d
        type(compound_result) :: out
        real(dp) :: a1,a2,x1,x2,rho
        out%critical_spot = bscallimps(s,kuo,v,r,t2-t1,d,kco)
        a1 = bs_d1(s,out%critical_spot,v,r,t1,d)
        a2 = a1-v*sqrt(t1)
        x1 = bs_d1(s,kuo,v,r,t2,d)
        x2 = x1-v*sqrt(t2)
        rho = sqrt(t1/t2)
        out%price = s*exp(-d*t2)*binormsdist(a1,x1,rho) - &
            kuo*exp(-r*t2)*binormsdist(a2,x2,rho) - &
            kco*exp(-r*t1)*norm_cdf(a2)
    end function calloncall

    function putoncall(s,kuo,kco,v,r,t1,t2,d) result(out)
        real(dp), intent(in) :: s,kuo,kco,v,r,t1,t2,d
        type(compound_result) :: out
        real(dp) :: a1,a2,x1,x2,rho
        out%critical_spot = bscallimps(s,kuo,v,r,t2-t1,d,kco)
        a1 = bs_d1(s,out%critical_spot,v,r,t1,d)
        a2 = a1-v*sqrt(t1)
        x1 = bs_d1(s,kuo,v,r,t2,d)
        x2 = x1-v*sqrt(t2)
        rho = sqrt(t1/t2)
        out%price = -s*exp(-d*t2)*binormsdist(-a1,x1,-rho) + &
            kuo*exp(-r*t2)*binormsdist(-a2,x2,-rho) + &
            kco*exp(-r*t1)*norm_cdf(-a2)
    end function putoncall

    function callonput(s,kuo,kco,v,r,t1,t2,d) result(out)
        real(dp), intent(in) :: s,kuo,kco,v,r,t1,t2,d
        type(compound_result) :: out
        real(dp) :: a1,a2,x1,x2,rho
        out%critical_spot = bsputimps(s,kuo,v,r,t2-t1,d,kco)
        a1 = bs_d1(s,out%critical_spot,v,r,t1,d)
        a2 = a1-v*sqrt(t1)
        x1 = bs_d1(s,kuo,v,r,t2,d)
        x2 = x1-v*sqrt(t2)
        rho = sqrt(t1/t2)
        out%price = -s*exp(-d*t2)*binormsdist(-a1,-x1,rho) + &
            kuo*exp(-r*t2)*binormsdist(-a2,-x2,rho) - &
            kco*exp(-r*t1)*norm_cdf(-a2)
    end function callonput

    function putonput(s,kuo,kco,v,r,t1,t2,d) result(out)
        real(dp), intent(in) :: s,kuo,kco,v,r,t1,t2,d
        type(compound_result) :: out
        real(dp) :: a1,a2,x1,x2,rho
        out%critical_spot = bsputimps(s,kuo,v,r,t2-t1,d,kco)
        a1 = bs_d1(s,out%critical_spot,v,r,t1,d)
        a2 = a1-v*sqrt(t1)
        x1 = bs_d1(s,kuo,v,r,t2,d)
        x2 = x1-v*sqrt(t2)
        rho = sqrt(t1/t2)
        out%price = s*exp(-d*t2)*binormsdist(a1,-x1,-rho) - &
            kuo*exp(-r*t2)*binormsdist(a2,-x2,-rho) + &
            kco*exp(-r*t1)*norm_cdf(a2)
    end function putonput

    subroutine optionsoncall(s,kuo,kco,v,r,t1,t2,d,call_price,put_price,critical_spot)
        real(dp), intent(in) :: s,kuo,kco,v,r,t1,t2,d
        real(dp), intent(out) :: call_price,put_price,critical_spot
        type(compound_result) :: x
        x = calloncall(s,kuo,kco,v,r,t1,t2,d)
        call_price = x%price
        critical_spot = x%critical_spot
        put_price = call_price-bscall(s,kuo,v,r,t2,d)+kco*exp(-r*t1)
    end subroutine optionsoncall

    subroutine optionsonput(s,kuo,kco,v,r,t1,t2,d,call_price,put_price,critical_spot)
        real(dp), intent(in) :: s,kuo,kco,v,r,t1,t2,d
        real(dp), intent(out) :: call_price,put_price,critical_spot
        type(compound_result) :: x
        x = callonput(s,kuo,kco,v,r,t1,t2,d)
        call_price = x%price
        critical_spot = x%critical_spot
        put_price = call_price-bsput(s,kuo,v,r,t2,d)+kco*exp(-r*t1)
    end subroutine optionsonput
end module derivmkts_compound
