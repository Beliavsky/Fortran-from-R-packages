! Exact Skvoretz biased-net triadic pseudolikelihood kernel translated from
! R/sna src/likelihood.c.  Upstream copyright (C) Carter T. Butts.
! Licensed under GPL-2.0-or-later; see COPYING.
module sna_bn_triad
    use sna_kinds, only : dp
    implicit none
    private
    public :: bn_lpt_triad, bn_triad_stats, bn_nll_triad
contains

    pure real(dp) function bn_lpt_m(k,pi,sigma,rho,d) result(v)
        integer,intent(in)::k
        real(dp),intent(in)::pi,sigma,rho,d
        v=log(max(1.0_dp-(1.0_dp-pi)*(1.0_dp-rho)**k*(1.0_dp-sigma)**k*(1.0_dp-d),tiny(1.0_dp))) + &
          log(max(1.0_dp-(1.0_dp-sigma)**k*(1.0_dp-d),tiny(1.0_dp)))
    end function bn_lpt_m
    pure real(dp) function bn_lpt_a(k,pi,sigma,rho,d) result(v)
        integer,intent(in)::k
        real(dp),intent(in)::pi,sigma,rho,d
        v=log(max(1.0_dp-(1.0_dp-sigma)**k*(1.0_dp-d),tiny(1.0_dp))) + &
          log(max((1.0_dp-pi)*(1.0_dp-rho)**k*(1.0_dp-sigma)**k*(1.0_dp-d),tiny(1.0_dp)))
    end function bn_lpt_a
    pure real(dp) function bn_lpt_n(k,pi,sigma,rho,d) result(v)
        integer,intent(in)::k
        real(dp),intent(in)::pi,sigma,rho,d
        real(dp)::calc
        calc=1.0_dp-exp(bn_lpt_m(k,pi,sigma,rho,d))-2.0_dp*exp(bn_lpt_a(k,pi,sigma,rho,d))
        if(calc<0.0_dp)calc=0.0_dp
        v=log(max(calc,tiny(1.0_dp)))
    end function bn_lpt_n
    pure real(dp) function bn_lpt_mp1(k,pi,sigma,rho,d) result(v)
        integer,intent(in)::k
        real(dp),intent(in)::pi,sigma,rho,d
        v=bn_lpt_m(k+1,pi,sigma,rho,d)
    end function bn_lpt_mp1
    pure real(dp) function bn_lpt_ap1(k,pi,sigma,rho,d) result(v)
        integer,intent(in)::k
        real(dp),intent(in)::pi,sigma,rho,d
        v=bn_lpt_a(k+1,pi,sigma,rho,d)
    end function bn_lpt_ap1
    pure real(dp) function bn_lpt_np1(k,pi,sigma,rho,d) result(v)
        integer,intent(in)::k
        real(dp),intent(in)::pi,sigma,rho,d
        v=bn_lpt_n(k+1,pi,sigma,rho,d)
    end function bn_lpt_np1
    pure real(dp) function bn_lpt_m1(pi,sigma,rho,d) result(v)
        real(dp),intent(in)::pi,sigma,rho,d
        v=log(max(sigma*(1.0_dp-(1.0_dp-sigma)*(1.0_dp-rho)),tiny(1.0_dp)))
    end function bn_lpt_m1
    pure real(dp) function bn_lpt_a1(pi,sigma,rho,d) result(v)
        real(dp),intent(in)::pi,sigma,rho,d
        v=log(max(sigma*(1.0_dp-sigma)*(1.0_dp-rho),tiny(1.0_dp)))
    end function bn_lpt_a1
    pure real(dp) function bn_lpt_n1(pi,sigma,rho,d) result(v)
        real(dp),intent(in)::pi,sigma,rho,d
        v=log(max(1.0_dp-sigma*(1.0_dp+(1.0_dp-sigma)*(1.0_dp-rho)),tiny(1.0_dp)))
    end function bn_lpt_n1
    pure real(dp) function bn_lpt_sr(pi,sigma,rho,d) result(v)
        real(dp),intent(in)::pi,sigma,rho,d
        v=log(max(1.0_dp-(1.0_dp-sigma)*(1.0_dp-rho),tiny(1.0_dp)))
    end function bn_lpt_sr
    pure real(dp) function bn_lpt_1msr(pi,sigma,rho,d) result(v)
        real(dp),intent(in)::pi,sigma,rho,d
        v=log(max((1.0_dp-sigma)*(1.0_dp-rho),tiny(1.0_dp)))
    end function bn_lpt_1msr

    pure real(dp) function bn_lpt_triad(xy,yx,yz,zy,xz,zx,kxy,kyz,kxz,pi,sigma,rho,d) result(lp)
        integer,intent(in)::xy,yx,yz,zy,xz,zx,kxy,kyz,kxz
        real(dp),intent(in)::pi,sigma,rho,d
        integer::code
        code=32*merge(1,0,xy>0)+16*merge(1,0,yx>0)+8*merge(1,0,yz>0)+4*merge(1,0,zy>0)+2*merge(1,0,xz>0)+merge(1,0,zx>0)
        select case(code)
        case(0)
            lp = ( ( bn_lpt_n ( kxy , pi , sigma , rho , d ) + bn_lpt_n ( kyz , pi , sigma , rho , d ) ) + bn_lpt_n &
                & ( kxz , pi , sigma , rho , d ) )
        case(1)
            lp = ( ( bn_lpt_n ( kxy , pi , sigma , rho , d ) + bn_lpt_n ( kyz , pi , sigma , rho , d ) ) + bn_lpt_a &
                & ( kxz , pi , sigma , rho , d ) )
        case(2)
            lp = ( ( bn_lpt_n ( kxy , pi , sigma , rho , d ) + bn_lpt_n ( kyz , pi , sigma , rho , d ) ) + bn_lpt_a &
                & ( kxz , pi , sigma , rho , d ) )
        case(3)
            lp = ( ( bn_lpt_n ( kxy , pi , sigma , rho , d ) + bn_lpt_n ( kyz , pi , sigma , rho , d ) ) + bn_lpt_m &
                & ( kxz , pi , sigma , rho , d ) )
        case(4)
            lp = ( ( bn_lpt_n ( kxy , pi , sigma , rho , d ) + bn_lpt_a ( kyz , pi , sigma , rho , d ) ) + bn_lpt_n &
                & ( kxz , pi , sigma , rho , d ) )
        case(5)
            lp = ( ( ( log ( ( exp ( bn_lpt_np1 ( kxy , pi , sigma , rho , d ) ) + ( 2.0_dp * exp ( ( bn_lpt_n ( &
                & kxy , pi , sigma , rho , d ) + bn_lpt_n1 ( pi , sigma , rho , d ) ) ) ) ) ) + bn_lpt_a ( kxz , pi &
                & , sigma , rho , d ) ) + bn_lpt_a ( kyz , pi , sigma , rho , d ) ) - log ( 3.0_dp ) )
        case(6)
            lp = ( ( bn_lpt_n ( kxy , pi , sigma , rho , d ) + bn_lpt_a ( kyz , pi , sigma , rho , d ) ) + bn_lpt_a &
                & ( kxz , pi , sigma , rho , d ) )
        case(7)
            lp = ( ( ( log ( ( exp ( bn_lpt_np1 ( kxy , pi , sigma , rho , d ) ) + ( 2.0_dp * exp ( ( bn_lpt_n ( &
                & kxy , pi , sigma , rho , d ) + bn_lpt_n1 ( pi , sigma , rho , d ) ) ) ) ) ) + bn_lpt_a ( kyz , pi &
                & , sigma , rho , d ) ) + bn_lpt_m ( kxz , pi , sigma , rho , d ) ) - log ( 3.0_dp ) )
        case(8)
            lp = ( ( bn_lpt_n ( kxy , pi , sigma , rho , d ) + bn_lpt_a ( kyz , pi , sigma , rho , d ) ) + bn_lpt_n &
                & ( kxz , pi , sigma , rho , d ) )
        case(9)
            lp = ( ( bn_lpt_n ( kxy , pi , sigma , rho , d ) + bn_lpt_a ( kyz , pi , sigma , rho , d ) ) + bn_lpt_a &
                & ( kxz , pi , sigma , rho , d ) )
        case(10)
            lp = ( ( bn_lpt_n ( kxy , pi , sigma , rho , d ) + bn_lpt_a ( kyz , pi , sigma , rho , d ) ) + bn_lpt_a &
                & ( kxz , pi , sigma , rho , d ) )
        case(11)
            lp = ( ( bn_lpt_n ( kxy , pi , sigma , rho , d ) + bn_lpt_a ( kyz , pi , sigma , rho , d ) ) + bn_lpt_m &
                & ( kxz , pi , sigma , rho , d ) )
        case(12)
            lp = ( ( bn_lpt_n ( kxy , pi , sigma , rho , d ) + bn_lpt_m ( kyz , pi , sigma , rho , d ) ) + bn_lpt_n &
                & ( kxz , pi , sigma , rho , d ) )
        case(13)
            lp = ( ( ( log ( ( exp ( bn_lpt_np1 ( kxy , pi , sigma , rho , d ) ) + ( 2.0_dp * exp ( ( bn_lpt_n ( &
                & kxy , pi , sigma , rho , d ) + bn_lpt_n1 ( pi , sigma , rho , d ) ) ) ) ) ) + bn_lpt_m ( kyz , pi &
                & , sigma , rho , d ) ) + bn_lpt_a ( kxz , pi , sigma , rho , d ) ) - log ( 3.0_dp ) )
        case(14)
            lp = ( ( bn_lpt_n ( kxy , pi , sigma , rho , d ) + bn_lpt_m ( kyz , pi , sigma , rho , d ) ) + bn_lpt_a &
                & ( kxz , pi , sigma , rho , d ) )
        case(15)
            lp = ( ( ( log ( ( exp ( bn_lpt_np1 ( kxy , pi , sigma , rho , d ) ) + ( 2.0_dp * exp ( ( bn_lpt_n ( &
                & kxy , pi , sigma , rho , d ) + bn_lpt_n1 ( pi , sigma , rho , d ) ) ) ) ) ) + bn_lpt_m ( kyz , pi &
                & , sigma , rho , d ) ) + bn_lpt_m ( kxz , pi , sigma , rho , d ) ) - log ( 3.0_dp ) )
        case(16)
            lp = ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + bn_lpt_n ( kyz , pi , sigma , rho , d ) ) + bn_lpt_n &
                & ( kxz , pi , sigma , rho , d ) )
        case(17)
            lp = ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + bn_lpt_n ( kyz , pi , sigma , rho , d ) ) + bn_lpt_a &
                & ( kxz , pi , sigma , rho , d ) )
        case(18)
            lp = ( ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + log ( ( exp ( bn_lpt_np1 ( kyz , pi , sigma , rho &
                & , d ) ) + ( 2.0_dp * exp ( ( bn_lpt_n ( kyz , pi , sigma , rho , d ) + bn_lpt_n1 ( pi , sigma , &
                & rho , d ) ) ) ) ) ) ) + bn_lpt_a ( kxz , pi , sigma , rho , d ) ) - log ( 3.0_dp ) )
        case(19)
            lp = ( ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + log ( ( exp ( bn_lpt_np1 ( kyz , pi , sigma , rho &
                & , d ) ) + ( 2.0_dp * exp ( ( bn_lpt_n ( kyz , pi , sigma , rho , d ) + bn_lpt_n1 ( pi , sigma , &
                & rho , d ) ) ) ) ) ) ) + bn_lpt_m ( kxz , pi , sigma , rho , d ) ) - log ( 3.0_dp ) )
        case(20)
            lp = ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + bn_lpt_a ( kyz , pi , sigma , rho , d ) ) + bn_lpt_n &
                & ( kxz , pi , sigma , rho , d ) )
        case(21)
            lp = ( ( ( log ( ( ( exp ( bn_lpt_ap1 ( kxy , pi , sigma , rho , d ) ) + ( 2.0_dp * exp ( ( bn_lpt_a ( &
                & kxy , pi , sigma , rho , d ) + bn_lpt_1msr ( pi , sigma , rho , d ) ) ) ) ) + ( 2.0_dp * exp ( ( &
                & bn_lpt_n ( kxy , pi , sigma , rho , d ) + bn_lpt_a1 ( pi , sigma , rho , d ) ) ) ) ) ) + bn_lpt_a &
                & ( kyz , pi , sigma , rho , d ) ) + bn_lpt_a ( kxz , pi , sigma , rho , d ) ) - log ( 3.0_dp ) )
        case(22)
            lp = ( ( ( log ( ( ( exp ( bn_lpt_ap1 ( kxy , pi , sigma , rho , d ) ) + ( 2.0_dp * exp ( ( bn_lpt_a ( &
                & kxy , pi , sigma , rho , d ) + bn_lpt_1msr ( pi , sigma , rho , d ) ) ) ) ) + ( 2.0_dp * exp ( ( &
                & bn_lpt_n ( kxy , pi , sigma , rho , d ) + bn_lpt_a1 ( pi , sigma , rho , d ) ) ) ) ) ) + bn_lpt_a &
                & ( kyz , pi , sigma , rho , d ) ) + bn_lpt_a ( kxz , pi , sigma , rho , d ) ) - log ( 3.0_dp ) )
        case(23)
            lp = ( ( ( log ( ( ( ( ( exp ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + bn_lpt_ap1 ( kyz , pi , &
                & sigma , rho , d ) ) ) + exp ( ( bn_lpt_ap1 ( kxy , pi , sigma , rho , d ) + bn_lpt_a ( kyz , pi , &
                & sigma , rho , d ) ) ) ) + exp ( ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + bn_lpt_a ( kyz , pi &
                & , sigma , rho , d ) ) + bn_lpt_1msr ( pi , sigma , rho , d ) ) ) ) + ( 2.0_dp * exp ( ( ( &
                & bn_lpt_a ( kxy , pi , sigma , rho , d ) + bn_lpt_n ( kyz , pi , sigma , rho , d ) ) + bn_lpt_a1 ( &
                & pi , sigma , rho , d ) ) ) ) ) + ( 2.0_dp * exp ( ( ( bn_lpt_n ( kxy , pi , sigma , rho , d ) + &
                & bn_lpt_a1 ( pi , sigma , rho , d ) ) + bn_lpt_a ( kyz , pi , sigma , rho , d ) ) ) ) ) ) + &
                & bn_lpt_m ( kxz , pi , sigma , rho , d ) ) + bn_lpt_1msr ( pi , sigma , rho , d ) ) - log ( 3.0_dp &
                & ) )
        case(24)
            lp = ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + bn_lpt_a ( kyz , pi , sigma , rho , d ) ) + bn_lpt_n &
                & ( kxz , pi , sigma , rho , d ) )
        case(25)
            lp = ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + bn_lpt_a ( kyz , pi , sigma , rho , d ) ) + bn_lpt_a &
                & ( kxz , pi , sigma , rho , d ) )
        case(26)
            lp = ( ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + log ( ( exp ( bn_lpt_ap1 ( kyz , pi , sigma , rho &
                & , d ) ) + ( 2.0_dp * exp ( ( bn_lpt_a ( kyz , pi , sigma , rho , d ) + bn_lpt_1msr ( pi , sigma , &
                & rho , d ) ) ) ) ) ) ) + bn_lpt_a ( kxz , pi , sigma , rho , d ) ) - log ( 3.0_dp ) )
        case(27)
            lp = ( ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + log ( ( exp ( bn_lpt_ap1 ( kyz , pi , sigma , rho &
                & , d ) ) + ( 2.0_dp * exp ( ( bn_lpt_a ( kyz , pi , sigma , rho , d ) + bn_lpt_1msr ( pi , sigma , &
                & rho , d ) ) ) ) ) ) ) + bn_lpt_m ( kxz , pi , sigma , rho , d ) ) - log ( 3.0_dp ) )
        case(28)
            lp = ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + bn_lpt_m ( kyz , pi , sigma , rho , d ) ) + bn_lpt_n &
                & ( kxz , pi , sigma , rho , d ) )
        case(29)
            lp = ( ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + bn_lpt_m ( kyz , pi , sigma , rho , d ) ) + log ( &
                & ( exp ( bn_lpt_ap1 ( kxz , pi , sigma , rho , d ) ) + ( 2.0_dp * exp ( ( bn_lpt_a ( kxz , pi , &
                & sigma , rho , d ) + bn_lpt_1msr ( pi , sigma , rho , d ) ) ) ) ) ) ) - log ( 3.0_dp ) )
        case(30)
            lp = ( ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + log ( ( ( exp ( bn_lpt_mp1 ( kyz , pi , sigma , &
                & rho , d ) ) + ( 2.0_dp * exp ( bn_lpt_m ( kyz , pi , sigma , rho , d ) ) ) ) + ( 4.0_dp * exp ( ( &
                & bn_lpt_a ( kyz , pi , sigma , rho , d ) + bn_lpt_sr ( pi , sigma , rho , d ) ) ) ) ) ) ) + &
                & bn_lpt_a ( kxz , pi , sigma , rho , d ) ) - log ( 3.0_dp ) )
        case(31)
            lp = log ( ( ( ( ( ( ( exp ( bn_lpt_m ( kyz , pi , sigma , rho , d ) ) * ( ( exp ( ( bn_lpt_ap1 ( kxy , &
                & pi , sigma , rho , d ) + bn_lpt_m ( kxz , pi , sigma , rho , d ) ) ) + exp ( ( ( bn_lpt_a ( kxy , &
                & pi , sigma , rho , d ) + bn_lpt_mp1 ( kxz , pi , sigma , rho , d ) ) + bn_lpt_1msr ( pi , sigma , &
                & rho , d ) ) ) ) + exp ( ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + bn_lpt_m ( kxz , pi , sigma &
                & , rho , d ) ) + bn_lpt_1msr ( pi , sigma , rho , d ) ) ) ) ) / 3.0_dp ) + ( ( 2.0_dp * exp ( ( ( &
                & ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + bn_lpt_m ( kyz , pi , sigma , rho , d ) ) + &
                & bn_lpt_a ( kxz , pi , sigma , rho , d ) ) + bn_lpt_sr ( pi , sigma , rho , d ) ) + bn_lpt_1msr ( &
                & pi , sigma , rho , d ) ) ) ) / 3.0_dp ) ) + ( ( exp ( ( ( bn_lpt_m ( kyz , pi , sigma , rho , d ) &
                & + bn_lpt_a ( kxz , pi , sigma , rho , d ) ) + bn_lpt_sr ( pi , sigma , rho , d ) ) ) * ( ( exp ( &
                & bn_lpt_a ( kxy , pi , sigma , rho , d ) ) + 1.0_dp ) + exp ( ( bn_lpt_a ( kxy , pi , sigma , rho &
                & , d ) + bn_lpt_1msr ( pi , sigma , rho , d ) ) ) ) ) / 3.0_dp ) ) + ( ( 2.0_dp * exp ( ( ( ( &
                & bn_lpt_n ( kxy , pi , sigma , rho , d ) + bn_lpt_m ( kyz , pi , sigma , rho , d ) ) + bn_lpt_m ( &
                & kxz , pi , sigma , rho , d ) ) + bn_lpt_a1 ( pi , sigma , rho , d ) ) ) ) / 3.0_dp ) ) + ( ( ( &
                & 2.0_dp * exp ( bn_lpt_m ( kyz , pi , sigma , rho , d ) ) ) * ( exp ( ( ( ( bn_lpt_n ( kxy , pi , &
                & sigma , rho , d ) + bn_lpt_a ( kxz , pi , sigma , rho , d ) ) + bn_lpt_a1 ( pi , sigma , rho , d &
                & ) ) + bn_lpt_sr ( pi , sigma , rho , d ) ) ) + exp ( ( ( ( bn_lpt_a ( kxy , pi , sigma , rho , d &
                & ) + bn_lpt_n ( kxz , pi , sigma , rho , d ) ) + bn_lpt_m1 ( pi , sigma , rho , d ) ) + &
                & bn_lpt_1msr ( pi , sigma , rho , d ) ) ) ) ) / 3.0_dp ) ) )
        case(32)
            lp = ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + bn_lpt_n ( kyz , pi , sigma , rho , d ) ) + bn_lpt_n &
                & ( kxz , pi , sigma , rho , d ) )
        case(33)
            lp = ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + bn_lpt_n ( kyz , pi , sigma , rho , d ) ) + bn_lpt_a &
                & ( kxz , pi , sigma , rho , d ) )
        case(34)
            lp = ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + bn_lpt_n ( kyz , pi , sigma , rho , d ) ) + bn_lpt_a &
                & ( kxz , pi , sigma , rho , d ) )
        case(35)
            lp = ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + bn_lpt_n ( kyz , pi , sigma , rho , d ) ) + bn_lpt_m &
                & ( kxz , pi , sigma , rho , d ) )
        case(36)
            lp = ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + bn_lpt_a ( kyz , pi , sigma , rho , d ) ) + bn_lpt_n &
                & ( kxz , pi , sigma , rho , d ) )
        case(37)
            lp = ( ( ( log ( ( exp ( bn_lpt_ap1 ( kxy , pi , sigma , rho , d ) ) + ( 2.0_dp * exp ( ( bn_lpt_a ( &
                & kxy , pi , sigma , rho , d ) + bn_lpt_1msr ( pi , sigma , rho , d ) ) ) ) ) ) + bn_lpt_a ( kyz , &
                & pi , sigma , rho , d ) ) + bn_lpt_a ( kxz , pi , sigma , rho , d ) ) - log ( 3.0_dp ) )
        case(38)
            lp = ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + bn_lpt_a ( kyz , pi , sigma , rho , d ) ) + bn_lpt_a &
                & ( kxz , pi , sigma , rho , d ) )
        case(39)
            lp = ( ( ( log ( ( exp ( bn_lpt_ap1 ( kxy , pi , sigma , rho , d ) ) + ( 2.0_dp * exp ( ( bn_lpt_a ( &
                & kxy , pi , sigma , rho , d ) + bn_lpt_1msr ( pi , sigma , rho , d ) ) ) ) ) ) + bn_lpt_a ( kyz , &
                & pi , sigma , rho , d ) ) + bn_lpt_m ( kxz , pi , sigma , rho , d ) ) - log ( 3.0_dp ) )
        case(40)
            lp = ( ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + bn_lpt_a ( kyz , pi , sigma , rho , d ) ) + log ( &
                & ( exp ( bn_lpt_np1 ( kxz , pi , sigma , rho , d ) ) + ( 2.0_dp * exp ( ( bn_lpt_n ( kxz , pi , &
                & sigma , rho , d ) + bn_lpt_n1 ( pi , sigma , rho , d ) ) ) ) ) ) ) - log ( 3.0_dp ) )
        case(41)
            lp = ( ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + bn_lpt_a ( kyz , pi , sigma , rho , d ) ) + log ( &
                & ( exp ( bn_lpt_ap1 ( kxz , pi , sigma , rho , d ) ) + ( 2.0_dp * exp ( ( bn_lpt_a ( kxz , pi , &
                & sigma , rho , d ) + bn_lpt_1msr ( pi , sigma , rho , d ) ) ) ) ) ) ) - log ( 3.0_dp ) )
        case(42)
            lp = ( ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + bn_lpt_a ( kyz , pi , sigma , rho , d ) ) + log ( &
                & ( exp ( bn_lpt_ap1 ( kxz , pi , sigma , rho , d ) ) + ( 2.0_dp * exp ( ( bn_lpt_a ( kxz , pi , &
                & sigma , rho , d ) + bn_lpt_1msr ( pi , sigma , rho , d ) ) ) ) ) ) ) - log ( 3.0_dp ) )
        case(43)
            lp = ( ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + bn_lpt_a ( kxy , pi , sigma , rho , d ) ) + log ( &
                & ( ( exp ( bn_lpt_mp1 ( kxz , pi , sigma , rho , d ) ) + ( 2.0_dp * exp ( bn_lpt_m ( kxz , pi , &
                & sigma , rho , d ) ) ) ) + ( 4.0_dp * exp ( ( bn_lpt_a ( kxz , pi , sigma , rho , d ) + bn_lpt_sr &
                & ( pi , sigma , rho , d ) ) ) ) ) ) ) - log ( 3.0_dp ) )
        case(44)
            lp = ( ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + bn_lpt_m ( kxy , pi , sigma , rho , d ) ) + log ( &
                & ( exp ( bn_lpt_np1 ( kxz , pi , sigma , rho , d ) ) + ( 2.0_dp * exp ( bn_lpt_n ( kxz , pi , &
                & sigma , rho , d ) ) ) ) ) ) - log ( 3.0_dp ) )
        case(45)
            lp = ( ( ( log ( ( ( exp ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + bn_lpt_ap1 ( kxz , pi , sigma , &
                & rho , d ) ) ) + exp ( ( bn_lpt_ap1 ( kxy , pi , sigma , rho , d ) + bn_lpt_a ( kxz , pi , sigma , &
                & rho , d ) ) ) ) + exp ( ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + bn_lpt_a ( kxz , pi , sigma &
                & , rho , d ) ) + bn_lpt_1msr ( pi , sigma , rho , d ) ) ) ) ) + bn_lpt_m ( kyz , pi , sigma , rho &
                & , d ) ) + bn_lpt_1msr ( pi , sigma , rho , d ) ) - log ( 3.0_dp ) )
        case(46)
            lp = ( ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + bn_lpt_m ( kyz , pi , sigma , rho , d ) ) + log ( &
                & ( exp ( bn_lpt_ap1 ( kxz , pi , sigma , rho , d ) ) + ( 2.0_dp * exp ( ( bn_lpt_a ( kxz , pi , &
                & sigma , rho , d ) + bn_lpt_1msr ( pi , sigma , rho , d ) ) ) ) ) ) ) - log ( 3.0_dp ) )
        case(47)
            lp = log ( ( ( ( ( exp ( bn_lpt_m ( kyz , pi , sigma , rho , d ) ) * ( ( exp ( ( bn_lpt_ap1 ( kxy , pi &
                & , sigma , rho , d ) + bn_lpt_m ( kxz , pi , sigma , rho , d ) ) ) + exp ( ( ( bn_lpt_a ( kxy , pi &
                & , sigma , rho , d ) + bn_lpt_mp1 ( kxz , pi , sigma , rho , d ) ) + bn_lpt_1msr ( pi , sigma , &
                & rho , d ) ) ) ) + exp ( ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + bn_lpt_m ( kxz , pi , sigma &
                & , rho , d ) ) + bn_lpt_1msr ( pi , sigma , rho , d ) ) ) ) ) / 3.0_dp ) + ( ( 2.0_dp * exp ( ( ( &
                & ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + bn_lpt_m ( kyz , pi , sigma , rho , d ) ) + &
                & bn_lpt_a ( kxz , pi , sigma , rho , d ) ) + bn_lpt_sr ( pi , sigma , rho , d ) ) + bn_lpt_1msr ( &
                & pi , sigma , rho , d ) ) ) ) / 3.0_dp ) ) + ( ( exp ( ( ( bn_lpt_m ( kyz , pi , sigma , rho , d ) &
                & + bn_lpt_a ( kxz , pi , sigma , rho , d ) ) + bn_lpt_sr ( pi , sigma , rho , d ) ) ) * ( exp ( &
                & bn_lpt_ap1 ( kxy , pi , sigma , rho , d ) ) + exp ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + &
                & bn_lpt_1msr ( pi , sigma , rho , d ) ) ) ) ) / 3.0_dp ) ) )
        case(48)
            lp = ( ( bn_lpt_m ( kxy , pi , sigma , rho , d ) + bn_lpt_n ( kyz , pi , sigma , rho , d ) ) + bn_lpt_n &
                & ( kxz , pi , sigma , rho , d ) )
        case(49)
            lp = ( ( bn_lpt_m ( kxy , pi , sigma , rho , d ) + bn_lpt_n ( kyz , pi , sigma , rho , d ) ) + bn_lpt_a &
                & ( kxz , pi , sigma , rho , d ) )
        case(50)
            lp = ( ( ( bn_lpt_m ( kxy , pi , sigma , rho , d ) + log ( ( exp ( bn_lpt_np1 ( kyz , pi , sigma , rho &
                & , d ) ) + ( 2.0_dp * exp ( bn_lpt_n ( kyz , pi , sigma , rho , d ) ) ) ) ) ) + bn_lpt_a ( kxz , &
                & pi , sigma , rho , d ) ) - log ( 3.0_dp ) )
        case(51)
            lp = ( ( ( bn_lpt_m ( kxy , pi , sigma , rho , d ) + log ( ( exp ( bn_lpt_np1 ( kyz , pi , sigma , rho &
                & , d ) ) + ( 2.0_dp * exp ( bn_lpt_n ( kyz , pi , sigma , rho , d ) ) ) ) ) ) + bn_lpt_a ( kxz , &
                & pi , sigma , rho , d ) ) - log ( 3.0_dp ) )
        case(52)
            lp = ( ( bn_lpt_m ( kxy , pi , sigma , rho , d ) + bn_lpt_a ( kyz , pi , sigma , rho , d ) ) + bn_lpt_n &
                & ( kxz , pi , sigma , rho , d ) )
        case(53)
            lp = ( ( ( log ( ( ( exp ( bn_lpt_mp1 ( kxy , pi , sigma , rho , d ) ) + ( 2.0_dp * exp ( bn_lpt_m ( &
                & kxy , pi , sigma , rho , d ) ) ) ) + ( 4.0_dp * exp ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + &
                & bn_lpt_sr ( pi , sigma , rho , d ) ) ) ) ) ) + bn_lpt_a ( kyz , pi , sigma , rho , d ) ) + &
                & bn_lpt_a ( kxz , pi , sigma , rho , d ) ) - log ( 3.0_dp ) )
        case(54)
            lp = ( ( ( bn_lpt_m ( kxy , pi , sigma , rho , d ) + log ( ( exp ( bn_lpt_ap1 ( kyz , pi , sigma , rho &
                & , d ) ) + ( 2.0_dp * exp ( ( bn_lpt_a ( kyz , pi , sigma , rho , d ) + bn_lpt_1msr ( pi , sigma , &
                & rho , d ) ) ) ) ) ) ) + bn_lpt_a ( kxz , pi , sigma , rho , d ) ) - log ( 3.0_dp ) )
        case(55)
            lp = log ( ( ( ( ( exp ( bn_lpt_m ( kxz , pi , sigma , rho , d ) ) * ( ( exp ( ( bn_lpt_ap1 ( kyz , pi &
                & , sigma , rho , d ) + bn_lpt_m ( kxy , pi , sigma , rho , d ) ) ) + exp ( ( ( bn_lpt_a ( kyz , pi &
                & , sigma , rho , d ) + bn_lpt_mp1 ( kxy , pi , sigma , rho , d ) ) + bn_lpt_1msr ( pi , sigma , &
                & rho , d ) ) ) ) + exp ( ( ( bn_lpt_a ( kyz , pi , sigma , rho , d ) + bn_lpt_m ( kxy , pi , sigma &
                & , rho , d ) ) + bn_lpt_1msr ( pi , sigma , rho , d ) ) ) ) ) / 3.0_dp ) + ( ( 2.0_dp * exp ( ( ( &
                & ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + bn_lpt_m ( kxz , pi , sigma , rho , d ) ) + &
                & bn_lpt_a ( kyz , pi , sigma , rho , d ) ) + bn_lpt_sr ( pi , sigma , rho , d ) ) + bn_lpt_1msr ( &
                & pi , sigma , rho , d ) ) ) ) / 3.0_dp ) ) + ( ( exp ( ( ( bn_lpt_m ( kxz , pi , sigma , rho , d ) &
                & + bn_lpt_a ( kxy , pi , sigma , rho , d ) ) + bn_lpt_sr ( pi , sigma , rho , d ) ) ) * ( exp ( &
                & bn_lpt_ap1 ( kyz , pi , sigma , rho , d ) ) + exp ( ( bn_lpt_a ( kyz , pi , sigma , rho , d ) + &
                & bn_lpt_1msr ( pi , sigma , rho , d ) ) ) ) ) / 3.0_dp ) ) )
        case(56)
            lp = ( ( ( bn_lpt_m ( kxy , pi , sigma , rho , d ) + bn_lpt_a ( kyz , pi , sigma , rho , d ) ) + log ( &
                & ( exp ( bn_lpt_np1 ( kxz , pi , sigma , rho , d ) ) + ( 2.0_dp * exp ( bn_lpt_n ( kxz , pi , &
                & sigma , rho , d ) ) ) ) ) ) - log ( 3.0_dp ) )
        case(57)
            lp = ( ( ( bn_lpt_m ( kxy , pi , sigma , rho , d ) + bn_lpt_a ( kyz , pi , sigma , rho , d ) ) + log ( &
                & ( exp ( bn_lpt_ap1 ( kxz , pi , sigma , rho , d ) ) + ( 2.0_dp * exp ( ( bn_lpt_a ( kxz , pi , &
                & sigma , rho , d ) + bn_lpt_1msr ( pi , sigma , rho , d ) ) ) ) ) ) ) - log ( 3.0_dp ) )
        case(58)
            lp = ( ( ( log ( ( ( exp ( ( bn_lpt_a ( kyz , pi , sigma , rho , d ) + bn_lpt_ap1 ( kxz , pi , sigma , &
                & rho , d ) ) ) + exp ( ( bn_lpt_ap1 ( kxz , pi , sigma , rho , d ) + bn_lpt_a ( kyz , pi , sigma , &
                & rho , d ) ) ) ) + exp ( ( ( bn_lpt_a ( kyz , pi , sigma , rho , d ) + bn_lpt_a ( kxz , pi , sigma &
                & , rho , d ) ) + bn_lpt_1msr ( pi , sigma , rho , d ) ) ) ) ) + bn_lpt_m ( kxy , pi , sigma , rho &
                & , d ) ) + bn_lpt_1msr ( pi , sigma , rho , d ) ) - log ( 3.0_dp ) )
        case(59)
            lp = log ( ( ( ( ( exp ( bn_lpt_m ( kxy , pi , sigma , rho , d ) ) * ( ( exp ( ( bn_lpt_ap1 ( kyz , pi &
                & , sigma , rho , d ) + bn_lpt_m ( kxz , pi , sigma , rho , d ) ) ) + exp ( ( ( bn_lpt_a ( kyz , pi &
                & , sigma , rho , d ) + bn_lpt_mp1 ( kxz , pi , sigma , rho , d ) ) + bn_lpt_1msr ( pi , sigma , &
                & rho , d ) ) ) ) + exp ( ( ( bn_lpt_a ( kyz , pi , sigma , rho , d ) + bn_lpt_m ( kxz , pi , sigma &
                & , rho , d ) ) + bn_lpt_1msr ( pi , sigma , rho , d ) ) ) ) ) / 3.0_dp ) + ( ( 2.0_dp * exp ( ( ( &
                & ( ( bn_lpt_m ( kxy , pi , sigma , rho , d ) + bn_lpt_a ( kyz , pi , sigma , rho , d ) ) + &
                & bn_lpt_a ( kxz , pi , sigma , rho , d ) ) + bn_lpt_sr ( pi , sigma , rho , d ) ) + bn_lpt_1msr ( &
                & pi , sigma , rho , d ) ) ) ) / 3.0_dp ) ) + ( ( exp ( ( ( bn_lpt_m ( kxy , pi , sigma , rho , d ) &
                & + bn_lpt_a ( kxz , pi , sigma , rho , d ) ) + bn_lpt_sr ( pi , sigma , rho , d ) ) ) * ( exp ( &
                & bn_lpt_ap1 ( kyz , pi , sigma , rho , d ) ) + exp ( ( bn_lpt_a ( kyz , pi , sigma , rho , d ) + &
                & bn_lpt_1msr ( pi , sigma , rho , d ) ) ) ) ) / 3.0_dp ) ) )
        case(60)
            lp = ( ( ( bn_lpt_m ( kxy , pi , sigma , rho , d ) + bn_lpt_m ( kyz , pi , sigma , rho , d ) ) + log ( &
                & ( exp ( bn_lpt_np1 ( kxz , pi , sigma , rho , d ) ) + ( 2.0_dp * exp ( bn_lpt_n ( kxz , pi , &
                & sigma , rho , d ) ) ) ) ) ) - log ( 3.0_dp ) )
        case(61)
            lp = log ( ( ( ( ( exp ( bn_lpt_m ( kyz , pi , sigma , rho , d ) ) * ( ( exp ( ( bn_lpt_ap1 ( kxz , pi &
                & , sigma , rho , d ) + bn_lpt_m ( kxy , pi , sigma , rho , d ) ) ) + exp ( ( ( bn_lpt_a ( kxz , pi &
                & , sigma , rho , d ) + bn_lpt_mp1 ( kxy , pi , sigma , rho , d ) ) + bn_lpt_1msr ( pi , sigma , &
                & rho , d ) ) ) ) + exp ( ( ( bn_lpt_a ( kxz , pi , sigma , rho , d ) + bn_lpt_m ( kxy , pi , sigma &
                & , rho , d ) ) + bn_lpt_1msr ( pi , sigma , rho , d ) ) ) ) ) / 3.0_dp ) + ( ( 2.0_dp * exp ( ( ( &
                & ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + bn_lpt_m ( kyz , pi , sigma , rho , d ) ) + &
                & bn_lpt_a ( kxz , pi , sigma , rho , d ) ) + bn_lpt_sr ( pi , sigma , rho , d ) ) + bn_lpt_1msr ( &
                & pi , sigma , rho , d ) ) ) ) / 3.0_dp ) ) + ( ( exp ( ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) &
                & + bn_lpt_m ( kyz , pi , sigma , rho , d ) ) + bn_lpt_sr ( pi , sigma , rho , d ) ) ) * ( exp ( &
                & bn_lpt_ap1 ( kxz , pi , sigma , rho , d ) ) + exp ( ( bn_lpt_a ( kxz , pi , sigma , rho , d ) + &
                & bn_lpt_1msr ( pi , sigma , rho , d ) ) ) ) ) / 3.0_dp ) ) )
        case(62)
            lp = log ( ( ( ( ( exp ( bn_lpt_m ( kxy , pi , sigma , rho , d ) ) * ( ( exp ( ( bn_lpt_ap1 ( kxz , pi &
                & , sigma , rho , d ) + bn_lpt_m ( kyz , pi , sigma , rho , d ) ) ) + exp ( ( ( bn_lpt_a ( kxz , pi &
                & , sigma , rho , d ) + bn_lpt_mp1 ( kyz , pi , sigma , rho , d ) ) + bn_lpt_1msr ( pi , sigma , &
                & rho , d ) ) ) ) + exp ( ( ( bn_lpt_a ( kxz , pi , sigma , rho , d ) + bn_lpt_m ( kyz , pi , sigma &
                & , rho , d ) ) + bn_lpt_1msr ( pi , sigma , rho , d ) ) ) ) ) / 3.0_dp ) + ( ( 2.0_dp * exp ( ( ( &
                & ( ( bn_lpt_m ( kxy , pi , sigma , rho , d ) + bn_lpt_a ( kyz , pi , sigma , rho , d ) ) + &
                & bn_lpt_a ( kxz , pi , sigma , rho , d ) ) + bn_lpt_sr ( pi , sigma , rho , d ) ) + bn_lpt_1msr ( &
                & pi , sigma , rho , d ) ) ) ) / 3.0_dp ) ) + ( ( exp ( ( ( bn_lpt_m ( kxy , pi , sigma , rho , d ) &
                & + bn_lpt_a ( kyz , pi , sigma , rho , d ) ) + bn_lpt_sr ( pi , sigma , rho , d ) ) ) * ( exp ( &
                & bn_lpt_ap1 ( kxz , pi , sigma , rho , d ) ) + exp ( ( bn_lpt_a ( kxz , pi , sigma , rho , d ) + &
                & bn_lpt_1msr ( pi , sigma , rho , d ) ) ) ) ) / 3.0_dp ) ) )
        case(63)
            lp = log ( ( ( ( ( ( ( exp ( ( ( bn_lpt_m ( kxy , pi , sigma , rho , d ) + bn_lpt_m ( kyz , pi , sigma &
                & , rho , d ) ) + bn_lpt_mp1 ( kxz , pi , sigma , rho , d ) ) ) + exp ( ( ( bn_lpt_m ( kxy , pi , &
                & sigma , rho , d ) + bn_lpt_mp1 ( kyz , pi , sigma , rho , d ) ) + bn_lpt_m ( kxz , pi , sigma , &
                & rho , d ) ) ) ) + exp ( ( ( bn_lpt_mp1 ( kxy , pi , sigma , rho , d ) + bn_lpt_m ( kyz , pi , &
                & sigma , rho , d ) ) + bn_lpt_m ( kxz , pi , sigma , rho , d ) ) ) ) / 3.0_dp ) + ( ( ( ( exp ( ( &
                & bn_lpt_a ( kxz , pi , sigma , rho , d ) + bn_lpt_sr ( pi , sigma , rho , d ) ) ) * ( ( exp ( ( &
                & bn_lpt_m ( kxy , pi , sigma , rho , d ) + bn_lpt_mp1 ( kyz , pi , sigma , rho , d ) ) ) + ( &
                & 2.0_dp * exp ( ( bn_lpt_m ( kxy , pi , sigma , rho , d ) + bn_lpt_m ( kyz , pi , sigma , rho , d &
                & ) ) ) ) ) + exp ( ( bn_lpt_mp1 ( kxy , pi , sigma , rho , d ) + bn_lpt_m ( kyz , pi , sigma , rho &
                & , d ) ) ) ) ) + ( exp ( ( bn_lpt_a ( kyz , pi , sigma , rho , d ) + bn_lpt_sr ( pi , sigma , rho &
                & , d ) ) ) * ( ( exp ( ( bn_lpt_m ( kxy , pi , sigma , rho , d ) + bn_lpt_mp1 ( kxz , pi , sigma , &
                & rho , d ) ) ) + ( 2.0_dp * exp ( ( bn_lpt_m ( kxy , pi , sigma , rho , d ) + bn_lpt_m ( kxz , pi &
                & , sigma , rho , d ) ) ) ) ) + exp ( ( bn_lpt_mp1 ( kxy , pi , sigma , rho , d ) + bn_lpt_m ( kxz &
                & , pi , sigma , rho , d ) ) ) ) ) ) + ( exp ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + &
                & bn_lpt_sr ( pi , sigma , rho , d ) ) ) * ( ( exp ( ( bn_lpt_mp1 ( kyz , pi , sigma , rho , d ) + &
                & bn_lpt_m ( kxz , pi , sigma , rho , d ) ) ) + ( 2.0_dp * exp ( ( ( bn_lpt_m ( kyz , pi , sigma , &
                & rho , d ) + bn_lpt_m ( kxy , pi , sigma , rho , d ) ) + bn_lpt_m ( kxz , pi , sigma , rho , d ) ) &
                & ) ) ) + exp ( ( bn_lpt_m ( kyz , pi , sigma , rho , d ) + bn_lpt_mp1 ( kxz , pi , sigma , rho , d &
                & ) ) ) ) ) ) / 3.0_dp ) ) + ( ( 4.0_dp * ( ( exp ( ( ( ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) &
                & + bn_lpt_m ( kyz , pi , sigma , rho , d ) ) + bn_lpt_a ( kxz , pi , sigma , rho , d ) ) + &
                & bn_lpt_sr ( pi , sigma , rho , d ) ) + bn_lpt_sr ( pi , sigma , rho , d ) ) ) + exp ( ( ( ( ( &
                & bn_lpt_m ( kxy , pi , sigma , rho , d ) + bn_lpt_a ( kyz , pi , sigma , rho , d ) ) + bn_lpt_a ( &
                & kxz , pi , sigma , rho , d ) ) + bn_lpt_sr ( pi , sigma , rho , d ) ) + bn_lpt_sr ( pi , sigma , &
                & rho , d ) ) ) ) + exp ( ( ( ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + bn_lpt_a ( kyz , pi , &
                & sigma , rho , d ) ) + bn_lpt_m ( kxz , pi , sigma , rho , d ) ) + bn_lpt_sr ( pi , sigma , rho , &
                & d ) ) + bn_lpt_sr ( pi , sigma , rho , d ) ) ) ) ) / 3.0_dp ) ) + ( ( ( exp ( ( ( ( ( bn_lpt_a ( &
                & kxy , pi , sigma , rho , d ) + bn_lpt_a ( kyz , pi , sigma , rho , d ) ) + bn_lpt_m ( kxz , pi , &
                & sigma , rho , d ) ) + bn_lpt_sr ( pi , sigma , rho , d ) ) + bn_lpt_sr ( pi , sigma , rho , d ) ) &
                & ) + exp ( ( ( ( ( bn_lpt_a ( kxy , pi , sigma , rho , d ) + bn_lpt_m ( kyz , pi , sigma , rho , d &
                & ) ) + bn_lpt_a ( kxz , pi , sigma , rho , d ) ) + bn_lpt_sr ( pi , sigma , rho , d ) ) + &
                & bn_lpt_sr ( pi , sigma , rho , d ) ) ) ) + exp ( ( ( ( ( bn_lpt_m ( kxy , pi , sigma , rho , d ) &
                & + bn_lpt_a ( kyz , pi , sigma , rho , d ) ) + bn_lpt_a ( kxz , pi , sigma , rho , d ) ) + &
                & bn_lpt_sr ( pi , sigma , rho , d ) ) + bn_lpt_sr ( pi , sigma , rho , d ) ) ) ) / 3.0_dp ) ) )
        end select
    end function bn_lpt_triad

    function bn_triad_stats(g) result(stats)
        real(dp),intent(in)::g(:,:)
        integer,allocatable::stats(:,:)
        integer::n,i,j,k
        n=size(g,1)
        allocate(stats(n,n))
        stats=0
        do i=1,n-1
        do j=i+1,n
            do k=1,n
            if(g(k,i)>0.0_dp.and.g(k,j)>0.0_dp)stats(i,j)=stats(i,j)+1
            end do
            stats(j,i)=stats(i,j)
        end do
        end do
    end function bn_triad_stats

    real(dp) function bn_nll_triad(params,g,stats) result(nll)
        real(dp),intent(in)::params(4),g(:,:)
        integer,intent(in)::stats(:,:)
        integer::n,i,j,k,xy,yx,yz,zy,xz,zx
        n=size(g,1)
        nll=0.0_dp
        do i=1,n-2
        do j=i+1,n-1
        do k=j+1,n
            xy=merge(1,0,g(i,j)>0.0_dp)
            yx=merge(1,0,g(j,i)>0.0_dp)
            yz=merge(1,0,g(j,k)>0.0_dp)
            zy=merge(1,0,g(k,j)>0.0_dp)
            xz=merge(1,0,g(i,k)>0.0_dp)
            zx=merge(1,0,g(k,i)>0.0_dp)
            nll=nll-bn_lpt_triad(xy,yx,yz,zy,xz,zx,stats(i,j),stats(j,k),stats(i,k),params(1),params(2),params(3),params(4))
        end do
        end do
        end do
    end function bn_nll_triad
end module sna_bn_triad
