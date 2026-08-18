! SPDX-License-Identifier: GPL-2.0-only
module elliptic_theta
    use elliptic_kinds, only : dp, pi, close_complex
    implicit none
    private
    public :: k_complete, nome, nome_k
    public :: theta1_q, theta2_q, theta3_q, theta4_q
    public :: theta1_m, theta2_m, theta3_m, theta4_m
    public :: theta1dash_q, theta1dashdash_q, theta1dashdashdash_q
    public :: theta1dash_m, theta1dashdash_m, theta1dashdashdash_m
    public :: theta1_dash_zero_q, theta1_dash_zero_m
    public :: theta_s, theta_c, theta_d, theta_n
    public :: sn, cn, dn, ns, nc, nd, sc, sd, cs, cd, ds, dc
    public :: cc, dd, nn, ss, h_fun, h1_fun, theta_big, theta1_big
    public :: e16_28_1, e16_28_2, e16_28_3, e16_28_4, e16_28_5
    public :: e16_37_1, e16_37_2, e16_37_3, e16_37_4
    public :: e16_38_1, e16_38_2, e16_38_3, e16_38_4

contains

    function k_complete(m, tol, maxiter) result(k)
        complex(dp), intent(in) :: m
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        complex(dp) :: k, a, b, an, bn
        real(dp) :: t
        integer :: it, nmax
        t = 32.0_dp * epsilon(1.0_dp)
        if (present(tol)) t = tol
        nmax = 100
        if (present(maxiter)) nmax = maxiter
        a = (1.0_dp, 0.0_dp)
        an = a
        b = sqrt((1.0_dp,0.0_dp) - m)
        do it = 1, nmax
            an = 0.5_dp * (a + b)
            bn = sqrt(a * b)
            if (abs(an - a) <= t * max(1.0_dp, abs(an), abs(a)) .and. it > 3) then
                k = pi / (2.0_dp * an)
                return
            end if
            a = an
            b = bn
        end do
        k = pi / (2.0_dp * an)
    end function k_complete

    function nome(m) result(q)
        complex(dp), intent(in) :: m
        complex(dp) :: q, k, kd
        k = k_complete(m)
        kd = k_complete((1.0_dp,0.0_dp) - m)
        q = exp(-pi * kd / k)
    end function nome

    function nome_k(kmod) result(q)
        complex(dp), intent(in) :: kmod
        complex(dp) :: q, kk, kd
        ! This intentionally follows elliptic::nome.k(), including its
        ! historical parameter convention.
        kk = k_complete(sqrt(kmod))
        kd = k_complete(sqrt((1.0_dp,0.0_dp) - kmod*kmod))
        q = exp(-pi * kd / kk)
    end function nome_k

    pure logical function converged(term, sumv, n, miniter, tol) result(ok)
        complex(dp), intent(in) :: term, sumv
        integer, intent(in) :: n, miniter
        real(dp), intent(in) :: tol
        ok = n >= miniter .and. abs(term) <= tol * max(1.0_dp, abs(sumv))
    end function converged

    function theta1_q(z, q, maxiter, tol) result(ans)
        complex(dp), intent(in) :: z, q
        integer, intent(in), optional :: maxiter
        real(dp), intent(in), optional :: tol
        complex(dp) :: ans, s, term
        real(dp) :: t
        integer :: n, nmax
        nmax = 10000; if (present(maxiter)) nmax = maxiter
        t = 16.0_dp*epsilon(1.0_dp); if (present(tol)) t = tol
        s = (0.0_dp,0.0_dp)
        do n = 0, nmax
            term = (-1.0_dp)**n * q**(n*(n+1)) * sin(real(2*n+1,dp)*z)
            s = s + term
            if (converged(term,s,n,3,t)) exit
        end do
        ans = 2.0_dp * q**0.25_dp * s
    end function theta1_q

    function theta2_q(z, q, maxiter, tol) result(ans)
        complex(dp), intent(in) :: z, q
        integer, intent(in), optional :: maxiter
        real(dp), intent(in), optional :: tol
        complex(dp) :: ans, s, term
        real(dp) :: t
        integer :: n, nmax
        nmax = 10000; if (present(maxiter)) nmax = maxiter
        t = 16.0_dp*epsilon(1.0_dp); if (present(tol)) t = tol
        s = (0.0_dp,0.0_dp)
        do n = 0, nmax
            term = q**(n*(n+1)) * cos(real(2*n+1,dp)*z)
            s = s + term
            if (converged(term,s,n,3,t)) exit
        end do
        ans = 2.0_dp * q**0.25_dp * s
    end function theta2_q

    function theta3_q(z, q, maxiter, tol) result(ans)
        complex(dp), intent(in) :: z, q
        integer, intent(in), optional :: maxiter
        real(dp), intent(in), optional :: tol
        complex(dp) :: ans, s, term
        real(dp) :: t
        integer :: n, nmax
        nmax = 10000; if (present(maxiter)) nmax = maxiter
        t = 16.0_dp*epsilon(1.0_dp); if (present(tol)) t = tol
        s = (0.0_dp,0.0_dp)
        do n = 1, nmax
            term = q**(n*n) * cos(2.0_dp*real(n,dp)*z)
            s = s + term
            if (converged(term,s,n,3,t)) exit
        end do
        ans = (1.0_dp,0.0_dp) + 2.0_dp*s
    end function theta3_q

    function theta4_q(z, q, maxiter, tol) result(ans)
        complex(dp), intent(in) :: z, q
        integer, intent(in), optional :: maxiter
        real(dp), intent(in), optional :: tol
        complex(dp) :: ans, s, term
        real(dp) :: t
        integer :: n, nmax
        nmax = 10000; if (present(maxiter)) nmax = maxiter
        t = 16.0_dp*epsilon(1.0_dp); if (present(tol)) t = tol
        s = (0.0_dp,0.0_dp)
        do n = 1, nmax
            term = (-1.0_dp)**n * q**(n*n) * cos(2.0_dp*real(n,dp)*z)
            s = s + term
            if (converged(term,s,n,3,t)) exit
        end do
        ans = (1.0_dp,0.0_dp) + 2.0_dp*s
    end function theta4_q

    function theta1dash_q(z, q, maxiter, tol) result(ans)
        complex(dp), intent(in) :: z, q
        integer, intent(in), optional :: maxiter
        real(dp), intent(in), optional :: tol
        complex(dp) :: ans, s, term
        real(dp) :: t, a
        integer :: n, nmax
        nmax=10000; if(present(maxiter)) nmax=maxiter
        t=16.0_dp*epsilon(1.0_dp); if(present(tol)) t=tol
        s=(0.0_dp,0.0_dp)
        do n=0,nmax
            a=real(2*n+1,dp)
            term=(-1.0_dp)**n*q**(n*(n+1))*a*cos(a*z)
            s=s+term
            if(converged(term,s,n,3,t)) exit
        end do
        ans=2.0_dp*q**0.25_dp*s
    end function theta1dash_q

    function theta1dashdash_q(z, q, maxiter, tol) result(ans)
        complex(dp), intent(in) :: z, q
        integer, intent(in), optional :: maxiter
        real(dp), intent(in), optional :: tol
        complex(dp) :: ans, s, term
        real(dp) :: t, a
        integer :: n, nmax
        nmax=10000; if(present(maxiter)) nmax=maxiter
        t=16.0_dp*epsilon(1.0_dp); if(present(tol)) t=tol
        s=(0.0_dp,0.0_dp)
        do n=0,nmax
            a=real(2*n+1,dp)
            term=-(-1.0_dp)**n*q**(n*(n+1))*a*a*sin(a*z)
            s=s+term
            if(converged(term,s,n,3,t)) exit
        end do
        ans=2.0_dp*q**0.25_dp*s
    end function theta1dashdash_q

    function theta1dashdashdash_q(z, q, maxiter, tol) result(ans)
        complex(dp), intent(in) :: z, q
        integer, intent(in), optional :: maxiter
        real(dp), intent(in), optional :: tol
        complex(dp) :: ans, s, term
        real(dp) :: t, a
        integer :: n, nmax
        nmax=10000; if(present(maxiter)) nmax=maxiter
        t=16.0_dp*epsilon(1.0_dp); if(present(tol)) t=tol
        s=(0.0_dp,0.0_dp)
        do n=0,nmax
            a=real(2*n+1,dp)
            term=-(-1.0_dp)**n*q**(n*(n+1))*a**3*cos(a*z)
            s=s+term
            if(converged(term,s,n,3,t)) exit
        end do
        ans=2.0_dp*q**0.25_dp*s
    end function theta1dashdashdash_q

    function theta1_m(z,m,maxiter,tol) result(v)
        complex(dp),intent(in)::z,m
        integer,intent(in),optional::maxiter
        real(dp),intent(in),optional::tol
        complex(dp)::v,q
        q=nome(m); v=theta1_q(z,q,maxiter,tol)
    end function
    function theta2_m(z,m,maxiter,tol) result(v)
        complex(dp),intent(in)::z,m
        integer,intent(in),optional::maxiter
        real(dp),intent(in),optional::tol
        complex(dp)::v,q
        q=nome(m); v=theta2_q(z,q,maxiter,tol)
    end function
    function theta3_m(z,m,maxiter,tol) result(v)
        complex(dp),intent(in)::z,m
        integer,intent(in),optional::maxiter
        real(dp),intent(in),optional::tol
        complex(dp)::v,q
        q=nome(m); v=theta3_q(z,q,maxiter,tol)
    end function
    function theta4_m(z,m,maxiter,tol) result(v)
        complex(dp),intent(in)::z,m
        integer,intent(in),optional::maxiter
        real(dp),intent(in),optional::tol
        complex(dp)::v,q
        q=nome(m); v=theta4_q(z,q,maxiter,tol)
    end function
    function theta1dash_m(z,m,maxiter,tol) result(v)
        complex(dp),intent(in)::z,m
        integer,intent(in),optional::maxiter
        real(dp),intent(in),optional::tol
        complex(dp)::v,q
        q=nome(m); v=theta1dash_q(z,q,maxiter,tol)
    end function
    function theta1dashdash_m(z,m,maxiter,tol) result(v)
        complex(dp),intent(in)::z,m
        integer,intent(in),optional::maxiter
        real(dp),intent(in),optional::tol
        complex(dp)::v,q
        q=nome(m); v=theta1dashdash_q(z,q,maxiter,tol)
    end function
    function theta1dashdashdash_m(z,m,maxiter,tol) result(v)
        complex(dp),intent(in)::z,m
        integer,intent(in),optional::maxiter
        real(dp),intent(in),optional::tol
        complex(dp)::v,q
        q=nome(m); v=theta1dashdashdash_q(z,q,maxiter,tol)
    end function

    function theta1_dash_zero_q(q) result(v)
        complex(dp),intent(in)::q
        complex(dp)::v
        v=theta2_q((0.0_dp,0.0_dp),q)*theta3_q((0.0_dp,0.0_dp),q)* &
          theta4_q((0.0_dp,0.0_dp),q)
    end function
    function theta1_dash_zero_m(m) result(v)
        complex(dp),intent(in)::m
        complex(dp)::v
        v=theta1_dash_zero_q(nome(m))
    end function

    function theta_s(u,m) result(v)
        complex(dp),intent(in)::u,m
        complex(dp)::v,k,x,q
        k=k_complete(m); q=nome(m); x=pi*u/(2.0_dp*k)
        v=2.0_dp*k*theta1_q(x,q)/(pi*theta2_q((0.0_dp,0.0_dp),q)* &
          theta3_q((0.0_dp,0.0_dp),q)*theta4_q((0.0_dp,0.0_dp),q))
    end function
    function theta_c(u,m) result(v)
        complex(dp),intent(in)::u,m
        complex(dp)::v,k,x,q
        k=k_complete(m); q=nome(m); x=pi*u/(2.0_dp*k)
        v=theta2_q(x,q)/theta2_q((0.0_dp,0.0_dp),q)
    end function
    function theta_d(u,m) result(v)
        complex(dp),intent(in)::u,m
        complex(dp)::v,k,x,q
        k=k_complete(m); q=nome(m); x=pi*u/(2.0_dp*k)
        v=theta3_q(x,q)/theta3_q((0.0_dp,0.0_dp),q)
    end function
    function theta_n(u,m) result(v)
        complex(dp),intent(in)::u,m
        complex(dp)::v,k,x,q
        k=k_complete(m); q=nome(m); x=pi*u/(2.0_dp*k)
        v=theta4_q(x,q)/theta4_q((0.0_dp,0.0_dp),q)
    end function

    function sn(u,m) result(v); complex(dp),intent(in)::u,m; complex(dp)::v; v=theta_s(u,m)/theta_n(u,m); end function
    function cn(u,m) result(v); complex(dp),intent(in)::u,m; complex(dp)::v; v=theta_c(u,m)/theta_n(u,m); end function
    function dn(u,m) result(v); complex(dp),intent(in)::u,m; complex(dp)::v; v=theta_d(u,m)/theta_n(u,m); end function
    function ns(u,m) result(v); complex(dp),intent(in)::u,m; complex(dp)::v; v=theta_n(u,m)/theta_s(u,m); end function
    function nc(u,m) result(v); complex(dp),intent(in)::u,m; complex(dp)::v; v=theta_n(u,m)/theta_c(u,m); end function
    function nd(u,m) result(v); complex(dp),intent(in)::u,m; complex(dp)::v; v=theta_n(u,m)/theta_d(u,m); end function
    function sc(u,m) result(v); complex(dp),intent(in)::u,m; complex(dp)::v; v=theta_s(u,m)/theta_c(u,m); end function
    function sd(u,m) result(v); complex(dp),intent(in)::u,m; complex(dp)::v; v=theta_s(u,m)/theta_d(u,m); end function
    function cs(u,m) result(v); complex(dp),intent(in)::u,m; complex(dp)::v; v=theta_c(u,m)/theta_s(u,m); end function
    function cd(u,m) result(v); complex(dp),intent(in)::u,m; complex(dp)::v; v=theta_c(u,m)/theta_d(u,m); end function
    function ds(u,m) result(v); complex(dp),intent(in)::u,m; complex(dp)::v; v=theta_d(u,m)/theta_s(u,m); end function
    function dc(u,m) result(v); complex(dp),intent(in)::u,m; complex(dp)::v; v=theta_d(u,m)/theta_c(u,m); end function
    function cc(u,m) result(v)
        complex(dp),intent(in)::u,m
        complex(dp)::v
        v=(1.0_dp,0.0_dp)+0.0_dp*(u+m)
    end function cc
    function dd(u,m) result(v)
        complex(dp),intent(in)::u,m
        complex(dp)::v
        v=(1.0_dp,0.0_dp)+0.0_dp*(u+m)
    end function dd
    function nn(u,m) result(v)
        complex(dp),intent(in)::u,m
        complex(dp)::v
        v=(1.0_dp,0.0_dp)+0.0_dp*(u+m)
    end function nn
    function ss(u,m) result(v)
        complex(dp),intent(in)::u,m
        complex(dp)::v
        v=(1.0_dp,0.0_dp)+0.0_dp*(u+m)
    end function ss

    function h_fun(u,m) result(v)
        complex(dp),intent(in)::u,m; complex(dp)::v,k
        k=k_complete(m); v=theta1_m(pi*u/(2.0_dp*k),m)
    end function
    function h1_fun(u,m) result(v)
        complex(dp),intent(in)::u,m; complex(dp)::v,k
        k=k_complete(m); v=theta2_m(pi*u/(2.0_dp*k),m)
    end function
    function theta_big(u,m) result(v)
        complex(dp),intent(in)::u,m; complex(dp)::v,k
        k=k_complete(m); v=theta4_m(pi*u/(2.0_dp*k),m)
    end function
    function theta1_big(u,m) result(v)
        complex(dp),intent(in)::u,m; complex(dp)::v,k
        k=k_complete(m); v=theta3_m(pi*u/(2.0_dp*k),m)
    end function

    function e16_28_1(z,m) result(v)
        complex(dp),intent(in)::z,m; complex(dp)::v
        v=theta1_m(z,m)**2*theta4_m((0.0_dp,0.0_dp),m)**2 - &
          theta3_m(z,m)**2*theta2_m((0.0_dp,0.0_dp),m)**2 + &
          theta2_m(z,m)**2*theta3_m((0.0_dp,0.0_dp),m)**2
    end function
    function e16_28_2(z,m) result(v)
        complex(dp),intent(in)::z,m; complex(dp)::v
        v=theta2_m(z,m)**2*theta4_m((0.0_dp,0.0_dp),m)**2 - &
          theta4_m(z,m)**2*theta2_m((0.0_dp,0.0_dp),m)**2 + &
          theta1_m(z,m)**2*theta3_m((0.0_dp,0.0_dp),m)**2
    end function
    function e16_28_3(z,m) result(v)
        complex(dp),intent(in)::z,m; complex(dp)::v
        v=theta3_m(z,m)**2*theta4_m((0.0_dp,0.0_dp),m)**2 - &
          theta4_m(z,m)**2*theta3_m((0.0_dp,0.0_dp),m)**2 + &
          theta1_m(z,m)**2*theta2_m((0.0_dp,0.0_dp),m)**2
    end function
    function e16_28_4(z,m) result(v)
        complex(dp),intent(in)::z,m; complex(dp)::v
        v=theta4_m(z,m)**2*theta4_m((0.0_dp,0.0_dp),m)**2 - &
          theta3_m(z,m)**2*theta3_m((0.0_dp,0.0_dp),m)**2 + &
          theta2_m(z,m)**2*theta2_m((0.0_dp,0.0_dp),m)**2
    end function
    function e16_28_5(m) result(v)
        complex(dp),intent(in)::m; complex(dp)::v
        v=theta2_m((0.0_dp,0.0_dp),m)**4+theta4_m((0.0_dp,0.0_dp),m)**4- &
          theta3_m((0.0_dp,0.0_dp),m)**4
    end function

    function e16_37_1(u,m,maxiter) result(v)
        complex(dp),intent(in)::u,m; integer,intent(in),optional::maxiter
        complex(dp)::v,q,k,x,p,pn; integer::n,nmax
        nmax=1000; if(present(maxiter))nmax=maxiter; q=nome(m); k=k_complete(m); x=pi*u/(2*k); p=1
        do n=1,nmax
            pn=p*(1-2*q**(2*n)*cos(2*x)+q**(4*n)); if(close_complex(p,pn,1e-14_dp).and.n>3)exit; p=pn
        end do
        v=(16*q/(m*(1-m)))**(1.0_dp/6.0_dp)*sin(x)*pn
    end function
    function e16_37_2(u,m,maxiter) result(v)
        complex(dp),intent(in)::u,m; integer,intent(in),optional::maxiter
        complex(dp)::v,q,k,x,p,pn; integer::n,nmax
        nmax=1000; if(present(maxiter))nmax=maxiter; q=nome(m); k=k_complete(m); x=pi*u/(2*k); p=1
        do n=1,nmax
            pn=p*(1+2*q**(2*n)*cos(2*x)+q**(4*n)); if(close_complex(p,pn,1e-14_dp).and.n>3)exit; p=pn
        end do
        v=(16*q*sqrt(1-m)/m)**(1.0_dp/6.0_dp)*cos(x)*pn
    end function
    function e16_37_3(u,m,maxiter) result(v)
        complex(dp),intent(in)::u,m; integer,intent(in),optional::maxiter
        complex(dp)::v,q,k,x,p,pn; integer::n,nmax
        nmax=1000; if(present(maxiter))nmax=maxiter; q=nome(m); k=k_complete(m); x=pi*u/(2*k); p=1
        do n=1,nmax
            pn=p*(1+2*q**(2*n-1)*cos(2*x)+q**(4*n-2)); if(close_complex(p,pn,1e-14_dp).and.n>3)exit; p=pn
        end do
        v=(m*(1-m)/(16*q))**(1.0_dp/12.0_dp)*pn
    end function
    function e16_37_4(u,m,maxiter) result(v)
        complex(dp),intent(in)::u,m; integer,intent(in),optional::maxiter
        complex(dp)::v,q,k,x,p,pn; integer::n,nmax
        nmax=1000; if(present(maxiter))nmax=maxiter; q=nome(m); k=k_complete(m); x=pi*u/(2*k); p=1
        do n=1,nmax
            pn=p*(1-2*q**(2*n-1)*cos(2*x)+q**(4*n-2)); if(close_complex(p,pn,1e-14_dp).and.n>3)exit; p=pn
        end do
        v=(m/(16*q*(1-m)**2))**(1.0_dp/12.0_dp)*pn
    end function
    function e16_38_1(u,m) result(v)
        complex(dp),intent(in)::u,m; complex(dp)::v,k,q,x
        k=k_complete(m);q=nome(m);x=pi*u/(2*k)
        v=theta1_q(x,q)/(2*q**0.25_dp)*sqrt(2*pi*sqrt(q)/(sqrt(m)*sqrt(1-m)*k))
    end function
    function e16_38_2(u,m) result(v)
        complex(dp),intent(in)::u,m; complex(dp)::v,k,q,x
        k=k_complete(m);q=nome(m);x=pi*u/(2*k)
        v=theta2_q(x,q)/(2*q**0.25_dp)*sqrt(2*pi*sqrt(q)/(sqrt(m)*k))
    end function
    function e16_38_3(u,m) result(v)
        complex(dp),intent(in)::u,m; complex(dp)::v,k,q,x
        k=k_complete(m);q=nome(m);x=pi*u/(2*k)
        v=theta3_q(x,q)*sqrt(pi/(2*k))
    end function
    function e16_38_4(u,m) result(v)
        complex(dp),intent(in)::u,m; complex(dp)::v,k,q,x
        k=k_complete(m);q=nome(m);x=pi*u/(2*k)
        v=theta4_q(x,q)*sqrt(pi/(2*sqrt(1-m)*k))
    end function
end module elliptic_theta
