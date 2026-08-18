! SPDX-License-Identifier: GPL-2.0-only
module elliptic_weierstrass
    use elliptic_kinds, only : dp, pi, ci
    use elliptic_theta, only : k_complete, theta1_q, theta2_q, theta3_q, theta4_q, &
        theta1dash_q, theta1dashdashdash_q, theta1_dash_zero_q
    implicit none
    private
    public :: elliptic_parameters, parameters_from_g, parameters_from_omega
    public :: equianharmonic_parameters, lemniscatic_parameters, pseudolemniscatic_parameters
    public :: cubic_roots, eee_cardano, half_periods, primitive_periods, is_primitive
    public :: mn_coordinates, fundamental_parallelogram
    public :: wp, wp_prime, weierstrass_sigma, weierstrass_zeta
    public :: coqueraux, ck_coefficients, amn_coefficients
    public :: wp_laurent, wp_prime_laurent, sigma_laurent, sigma_prime_laurent, zeta_laurent
    public :: g2_from_periods, g3_from_periods, g_from_periods, e18_10_9

    type :: elliptic_parameters
        complex(dp) :: omega(2) = (0.0_dp,0.0_dp)
        complex(dp) :: q = (0.0_dp,0.0_dp)
        complex(dp) :: e(3) = (0.0_dp,0.0_dp)
        complex(dp) :: g(2) = (0.0_dp,0.0_dp)
        complex(dp) :: delta = (0.0_dp,0.0_dp)
        complex(dp) :: eta(3) = (0.0_dp,0.0_dp)
    end type elliptic_parameters
contains
    function cubic_roots(g) result(r)
        complex(dp), intent(in) :: g(2)
        complex(dp) :: r(3),old(3),num,den
        complex(dp) :: eps3
        real(dp) :: rad
        integer :: i,j,it
        rad=max(1.0_dp,abs(g(1)),abs(g(2)))
        eps3=exp(2.0_dp*pi*ci/3.0_dp)
        r=[cmplx(rad,0.0_dp,dp),rad*eps3,rad*eps3**2]
        r(1)=r(1)+(0.137_dp,0.071_dp)
        do it=1,200
            old=r
            do i=1,3
                num=4.0_dp*old(i)**3-g(1)*old(i)-g(2)
                den=(4.0_dp,0.0_dp)
                do j=1,3
                    if(j/=i) den=den*(old(i)-old(j))
                end do
                r(i)=old(i)-num/den
            end do
            if(maxval(abs(r-old))<=1.0e-14_dp*max(1.0_dp,maxval(abs(r))))exit
        end do
        if(abs(aimag(g(1)))<1.0e-14_dp .and. abs(aimag(g(2)))<1.0e-14_dp) call order_real_roots(r,g)
    end function cubic_roots

    subroutine order_real_roots(r,g)
        complex(dp),intent(inout)::r(3)
        complex(dp),intent(in)::g(2)
        complex(dp)::t
        real(dp)::disc
        integer::i,j,ip,in,ir
        disc=real(g(1),dp)**3-27.0_dp*real(g(2),dp)**2
        if(disc>=-1.0e-12_dp)then
            do i=1,2
                do j=i+1,3
                    if(real(r(j),dp)>real(r(i),dp))then;t=r(i);r(i)=r(j);r(j)=t;end if
                end do
            end do
            r=cmplx(real(r,dp),0.0_dp,dp)
        else
            ip=maxloc(aimag(r),1);in=minloc(aimag(r),1);ir=6-ip-in
            r=[r(ip),r(ir),r(in)]
            if(abs(aimag(r(2)))<1.0e-12_dp)r(2)=cmplx(real(r(2),dp),0.0_dp,dp)
        end if
    end subroutine order_real_roots

    function eee_cardano(g) result(e)
        complex(dp),intent(in)::g(2)
        complex(dp)::e(3),g2,g3,d,eps3,alpha,beta,gamma_c
        g2=g(1);g3=g(2);eps3=exp(2*pi*ci/3.0_dp)
        if(abs(g2)<=tiny(1.0_dp))then
            e(1)=(g3/4.0_dp)**(1.0_dp/3.0_dp);e(2)=eps3*e(1);e(3)=e(1)/eps3;return
        end if
        if(abs(g3)<=tiny(1.0_dp))then
            e(1)=sqrt(g2)/2.0_dp;e(2)=-e(1);e(3)=0;return
        end if
        d=g3*g3-g2**3/27.0_dp;alpha=sqrt(g2/12.0_dp)
        beta=(g3+sqrt(d))**(1.0_dp/3.0_dp)/(2.0_dp*alpha)
        gamma_c=1.0_dp/beta
        e(1)=alpha*(beta+gamma_c)
        e(2)=alpha*(eps3*beta+gamma_c/eps3)
        e(3)=alpha*(beta/eps3+gamma_c*eps3)
    end function eee_cardano

    function primitive_periods(p,nsearch) result(out)
        complex(dp),intent(in)::p(2)
        integer,intent(in),optional::nsearch
        complex(dp)::out(2),vals(128),tmpz
        integer::co(128,2),n,a,b,i,j,k,tm(2),first(2),second(2),det
        n=3;if(present(nsearch))n=nsearch;k=0
        do a=0,n
            do b=n,-n,-1
                if(a==0.and.b==0)cycle
                k=k+1;co(k,:)=[a,b];vals(k)=a*p(1)+b*p(2)
            end do
        end do
        do i=1,k-1
            do j=i+1,k
                if(abs(vals(j))<abs(vals(i))-1e-15_dp .or. &
                   (abs(abs(vals(j))-abs(vals(i)))<1e-15_dp .and. &
                    (abs(co(j,1)-1)<abs(co(i,1)-1) .or. &
                     (abs(co(j,1)-1)==abs(co(i,1)-1).and.abs(co(j,2))<abs(co(i,2))))))then
                    tmpz=vals(i);vals(i)=vals(j);vals(j)=tmpz;tm=co(i,:);co(i,:)=co(j,:);co(j,:)=tm
                end if
            end do
        end do
        first=co(1,:);second=co(2,:)
        do i=2,k
            det=first(1)*co(i,2)-first(2)*co(i,1)
            if(det/=0)then;second=co(i,:);exit;end if
        end do
        out(1)=first(1)*p(1)+first(2)*p(2)
        out(2)=second(1)*p(1)+second(2)*p(2)
        if(real(out(1),dp)<0.0_dp)out(1)=-out(1)
        if(aimag(out(2)/out(1))<0.0_dp)out(2)=-out(2)
    end function primitive_periods

    logical function is_primitive(p,tol) result(ok)
        complex(dp),intent(in)::p(2)
        real(dp),intent(in),optional::tol
        real(dp)::t
        complex(dp)::q(2)
        t=1.0e-5_dp;if(present(tol))t=tol;q=primitive_periods(p)
        ok=maxval(abs(q-p))<t
    end function is_primitive

    function half_periods(e,g,make_primitive) result(o)
        complex(dp),intent(in),optional::e(3),g(2)
        logical,intent(in),optional::make_primitive
        complex(dp)::o(2),ee(3),m,s
        logical::prim
        if(present(e))then;ee=e;else if(present(g))then;ee=cubic_roots(g);else;error stop 'half_periods: supply e or g';end if
        m=(ee(2)-ee(3))/(ee(1)-ee(3));s=sqrt(ee(1)-ee(3))
        o(1)=k_complete(m)/s;o(2)=ci*k_complete(1.0_dp-m)/s
        prim=.true.;if(present(make_primitive))prim=make_primitive
        if(prim)o=primitive_periods(o)
    end function half_periods

    function mn_coordinates(z,p) result(mn)
        complex(dp),intent(in)::z,p(2)
        real(dp)::mn(2),den1,den2
        den1=real(p(1),dp)*aimag(p(2))-aimag(p(1))*real(p(2),dp)
        den2=real(p(2),dp)*aimag(p(1))-aimag(p(2))*real(p(1),dp)
        mn(1)=(real(z,dp)*aimag(p(2))-aimag(z)*real(p(2),dp))/den1
        mn(2)=(real(z,dp)*aimag(p(1))-aimag(z)*real(p(1),dp))/den2
    end function mn_coordinates

    function fundamental_parallelogram(z,p,mn_integer) result(out)
        complex(dp),intent(in)::z,p(2)
        integer,intent(out),optional::mn_integer(2)
        complex(dp)::out
        real(dp)::x(2)
        integer::k(2)
        x=mn_coordinates(z,p);k=nint(x);out=z-k(1)*p(1)-k(2)*p(2)
        if(present(mn_integer))mn_integer=k
    end function fundamental_parallelogram

    function coqueraux(z,g,nstep,use_fpp) result(v)
        complex(dp),intent(in)::z,g(2)
        integer,intent(in),optional::nstep
        logical,intent(in),optional::use_fpp
        complex(dp)::v,z0,x
        integer::i,n
        logical::uf
        n=5;if(present(nstep))n=nstep;uf=.false.;if(present(use_fpp))uf=use_fpp
        z0=z
        if(uf)z0=fundamental_parallelogram(z0,2.0_dp*half_periods(g=g))
        z0=z0/(2.0_dp**n);x=1.0_dp/z0**2+z0**2*(g(1)/20.0_dp+z0**2*g(2)/28.0_dp)
        do i=1,n
            x=-2.0_dp*x+(6.0_dp*x*x-g(1)/2.0_dp)**2/(4.0_dp*(4.0_dp*x**3-g(1)*x-g(2)))
        end do
        v=x
    end function coqueraux

    function reorder_roots_for_periods(g,o) result(e)
        complex(dp),intent(in)::g(2),o(2)
        complex(dp)::e(3),raw(3),approx(3)
        logical::used(3)
        integer::i,j,k
        raw=cubic_roots(g);approx(1)=coqueraux(o(1),g);approx(3)=coqueraux(o(2),g);approx(2)=-approx(1)-approx(3)
        used=.false.
        do i=1,3
            k=0
            do j=1,3
                if(used(j))cycle
                if(k==0)then
                    k=j
                else if(abs(raw(j)-approx(i))<abs(raw(k)-approx(i)))then
                    k=j
                end if
            end do
            e(i)=raw(k);used(k)=.true.
        end do
    end function reorder_roots_for_periods

    function parameters_from_g(g) result(par)
        complex(dp),intent(in)::g(2)
        type(elliptic_parameters)::par
        complex(dp)::o(2),e0(3),q,eta1,eta2
        e0=cubic_roots(g);o=half_periods(e=e0,make_primitive=.true.);par%e=reorder_roots_for_periods(g,o)
        par%omega=o;par%g=g;par%delta=g(1)**3-27.0_dp*g(2)**2
        q=exp(pi*ci*o(2)/o(1));par%q=q
        eta1=-pi*pi*theta1dashdashdash_q((0.0_dp,0.0_dp),q)/ &
            (12.0_dp*o(1)*theta1_dash_zero_q(q))
        eta2=o(2)/o(1)*eta1-pi*ci/(2.0_dp*o(1))
        par%eta=[eta1,eta2,-eta1-eta2]
    end function parameters_from_g

    function parameters_from_omega(omega) result(par)
        complex(dp),intent(in)::omega(2)
        type(elliptic_parameters)::par
        complex(dp)::o(2),g(2),q,eta1,eta2
        o=primitive_periods(omega)
        g=g_from_periods(o)
        par%omega=o
        par%g=g
        par%e=reorder_roots_for_periods(g,o)
        par%delta=g(1)**3-27.0_dp*g(2)**2
        q=exp(pi*ci*o(2)/o(1))
        par%q=q
        eta1=-pi*pi*theta1dashdashdash_q((0.0_dp,0.0_dp),q)/ &
            (12.0_dp*o(1)*theta1_dash_zero_q(q))
        eta2=o(2)/o(1)*eta1-pi*ci/(2.0_dp*o(1))
        par%eta=[eta1,eta2,-eta1-eta2]
    end function parameters_from_omega

    function g2_from_periods(b) result(g2)
        complex(dp),intent(in)::b(2)
        complex(dp)::g2,p1,tau,q,t2,t3
        p1=b(1);tau=b(2)/p1;q=exp(pi*ci*tau)
        t2=theta2_q((0.0_dp,0.0_dp),q);t3=theta3_q((0.0_dp,0.0_dp),q)
        g2=(pi/p1)**4*(t2**8-t3**4*t2**4+t3**8)/12.0_dp
    end function g2_from_periods

    function g3_from_periods(b) result(g3)
        complex(dp),intent(in)::b(2)
        complex(dp)::g3,p1,tau,q,t2,t3
        p1=b(1);tau=b(2)/p1;q=exp(pi*ci*tau)
        t2=theta2_q((0.0_dp,0.0_dp),q)**4;t3=theta3_q((0.0_dp,0.0_dp),q)**4
        g3=(pi/(2.0_dp*p1))**6*((8.0_dp/27.0_dp)*(t2**3+t3**3) - &
            (4.0_dp/9.0_dp)*(t2+t3)*t2*t3)
    end function g3_from_periods

    function g_from_periods(b) result(g)
        complex(dp),intent(in)::b(2);complex(dp)::g(2)
        g=[g2_from_periods(b),g3_from_periods(b)]
    end function g_from_periods

    function wp(z,g,params,use_fpp) result(v)
        complex(dp),intent(in)::z
        complex(dp),intent(in),optional::g(2)
        type(elliptic_parameters),intent(in),optional::params
        logical,intent(in),optional::use_fpp
        complex(dp)::v,zz,x,o,q,t10
        type(elliptic_parameters)::p
        logical::uf
        if (present(params)) then
            p=params
        else if (present(g)) then
            p=parameters_from_g(g)
        else
            error stop 'wp: supply g or params'
        end if
        zz=z;uf=.true.;if(present(use_fpp))uf=use_fpp;if(uf)zz=fundamental_parallelogram(zz,2.0_dp*p%omega)
        o=p%omega(1);q=p%q;x=pi*zz/(2.0_dp*o);t10=theta1_dash_zero_q(q)
        v=p%e(1)+pi*pi/(4.0_dp*o*o)*(t10*theta2_q(x,q)/theta2_q((0.0_dp,0.0_dp),q)/theta1_q(x,q))**2
    end function wp

    function wp_prime(z,g,params,use_fpp) result(v)
        complex(dp),intent(in)::z
        complex(dp),intent(in),optional::g(2)
        type(elliptic_parameters),intent(in),optional::params
        logical,intent(in),optional::use_fpp
        complex(dp)::v,zz,x,o,q
        type(elliptic_parameters)::p
        logical::uf
        if (present(params)) then
            p=params
        else if (present(g)) then
            p=parameters_from_g(g)
        else
            error stop 'wp_prime: supply g or params'
        end if
        zz=z;uf=.true.;if(present(use_fpp))uf=use_fpp;if(uf)zz=fundamental_parallelogram(zz,2.0_dp*p%omega)
        o=p%omega(1);q=p%q;x=pi*zz/(2.0_dp*o)
        v=-pi**3/(4.0_dp*o**3)*theta2_q(x,q)*theta3_q(x,q)*theta4_q(x,q)* &
          theta1dash_q((0.0_dp,0.0_dp),q)**3/(theta2_q((0.0_dp,0.0_dp),q)* &
          theta3_q((0.0_dp,0.0_dp),q)*theta4_q((0.0_dp,0.0_dp),q)*theta1_q(x,q)**3)
    end function wp_prime

    recursive function weierstrass_sigma(z,g,params,use_theta) result(v)
        complex(dp),intent(in)::z
        complex(dp),intent(in),optional::g(2)
        type(elliptic_parameters),intent(in),optional::params
        logical,intent(in),optional::use_theta
        complex(dp)::v,o,q,zz
        type(elliptic_parameters)::p
        integer::mn(2),sgn
        logical::ut
        if (present(params)) then
            p=params
        else if (present(g)) then
            p=parameters_from_g(g)
        else
            error stop 'sigma: supply g or params'
        end if
        ut=.true.;if(present(use_theta))ut=use_theta
        if(ut)then
            o=p%omega(1);q=p%q
            v=2.0_dp*o/pi*exp(p%eta(1)*z*z/(2.0_dp*o))*theta1_q(pi*z/(2.0_dp*o),q)/theta1dash_q((0.0_dp,0.0_dp),q)
        else
            zz=fundamental_parallelogram(z,2.0_dp*p%omega,mn)
            sgn=mn(1)+mn(2)+mn(1)*mn(2)
            v=(-1.0_dp)**sgn*weierstrass_sigma(zz,params=p,use_theta=.true.)* &
              exp((zz+mn(1)*p%omega(1)+mn(2)*p%omega(2))* &
                  (2.0_dp*mn(1)*p%eta(1)+2.0_dp*mn(2)*p%eta(2)))
        end if
    end function weierstrass_sigma

    recursive function weierstrass_zeta(z,g,params,use_fpp) result(v)
        complex(dp),intent(in)::z
        complex(dp),intent(in),optional::g(2)
        type(elliptic_parameters),intent(in),optional::params
        logical,intent(in),optional::use_fpp
        complex(dp)::v,zz,o,q
        type(elliptic_parameters)::p
        integer::mn(2)
        logical::uf
        if (present(params)) then
            p=params
        else if (present(g)) then
            p=parameters_from_g(g)
        else
            error stop 'zeta: supply g or params'
        end if
        uf=.true.;if(present(use_fpp))uf=use_fpp
        if(uf)then
            zz=fundamental_parallelogram(z,2.0_dp*p%omega,mn)
            v=weierstrass_zeta(zz,params=p,use_fpp=.false.)+2.0_dp*mn(1)*p%eta(1)+2.0_dp*mn(2)*p%eta(2)
        else
            o=p%omega(1);q=exp(pi*ci*p%omega(2)/p%omega(1));zz=pi*z/(2.0_dp*o)
            v=z*p%eta(1)/o+pi*theta1dash_q(zz,q)/(2.0_dp*o*theta1_q(zz,q))
        end if
    end function weierstrass_zeta

    function ck_coefficients(g,n) result(c)
        complex(dp),intent(in)::g(2)
        integer,intent(in)::n
        complex(dp),allocatable::c(:)
        integer::k,j
        allocate(c(n));c=(0.0_dp,0.0_dp);if(n<2)return
        c(2)=g(1)/20.0_dp;if(n>=3)c(3)=g(2)/28.0_dp
        do k=4,n
            do j=2,k-2
                c(k)=c(k)+c(j)*c(k-j)
            end do
            c(k)=3.0_dp*c(k)/real((2*k+1)*(k-3),dp)
        end do
    end function ck_coefficients

    function wp_laurent(z,g,nmax) result(v)
        complex(dp),intent(in)::z,g(2);integer,intent(in),optional::nmax
        complex(dp)::v,zz,term,old
        complex(dp),allocatable::c(:)
        integer::k,n
        n=80;if(present(nmax))n=nmax;c=ck_coefficients(g,n);v=1.0_dp/z**2;zz=z*z
        do k=2,n
            term=c(k)*zz;old=v;v=v+term;if(abs(v-old)<1e-15_dp*max(1.0_dp,abs(v)).and.abs(c(k))>0)exit;zz=zz*z*z
        end do
    end function wp_laurent

    function wp_prime_laurent(z,g,nmax) result(v)
        complex(dp),intent(in)::z,g(2);integer,intent(in),optional::nmax
        complex(dp)::v,zz,term,old
        complex(dp),allocatable::c(:)
        integer::k,n
        n=80;if(present(nmax))n=nmax;c=ck_coefficients(g,n);v=-2.0_dp/z**3;zz=z
        do k=2,n
            term=real(2*k-2,dp)*c(k)*zz;old=v;v=v+term
            if(abs(v-old)<1e-15_dp*max(1.0_dp,abs(v)).and.abs(c(k))>0)exit
            zz=zz*z*z
        end do
    end function wp_prime_laurent

    function zeta_laurent(z,g,nmax) result(v)
        complex(dp),intent(in)::z,g(2);integer,intent(in),optional::nmax
        complex(dp)::v,zz,term,old
        complex(dp),allocatable::c(:)
        integer::k,n
        n=80;if(present(nmax))n=nmax;c=ck_coefficients(g,n);v=1.0_dp/z;zz=z**3
        do k=2,n
            term=-c(k)*zz/real(2*k-1,dp);old=v;v=v+term
            if(abs(v-old)<1e-15_dp*max(1.0_dp,abs(v)).and.abs(c(k))>0)exit
            zz=zz*z*z
        end do
    end function zeta_laurent

    function amn_coefficients(u) result(a)
        integer,intent(in)::u
        real(dp),allocatable::a(:,:)
        integer::s,m,n
        allocate(a(u,u))
        a=0.0_dp
        a(1,1)=1.0_dp
        do s=1,u-1
            do n=0,s
                m=s-n
                a(n+1,m+1)=3.0_dp*real(m+1,dp)*asub(a,m+1,n-1)+ &
                    (16.0_dp/3.0_dp)*real(n+1,dp)*asub(a,m-2,n+1)- &
                    (1.0_dp/3.0_dp)*real((2*m+3*n-1)*(4*m+6*n-1),dp)* &
                    asub(a,m-1,n)
            end do
        end do
    end function amn_coefficients

    pure real(dp) function asub(a,m,n) result(v)
        real(dp),intent(in)::a(:,:);integer,intent(in)::m,n
        if(m<0.or.n<0.or.m>=size(a,2).or.n>=size(a,1))then;v=0.0_dp;else;v=a(n+1,m+1);end if
    end function asub

    function sigma_laurent(z,g,nmax) result(v)
        complex(dp),intent(in)::z,g(2);integer,intent(in),optional::nmax
        complex(dp)::v
        real(dp),allocatable::a(:,:)
        integer::u,m,n
        u=8;if(present(nmax))u=nmax;a=amn_coefficients(u);v=(0.0_dp,0.0_dp)
        do n=0,u-1;do m=0,u-1
            if(abs(a(n+1,m+1))<=tiny(1.0_dp))cycle
            ! Coefficients a are real but invariants may be complex; use complex expression below.
            v=v+a(n+1,m+1)*(g(1)/2.0_dp)**m*(2.0_dp*g(2))**n* &
                z**(4*m+6*n+1)/gamma(real(4*m+6*n+2,dp))
        end do;end do
    end function sigma_laurent

    function sigma_prime_laurent(z,g,nmax) result(v)
        complex(dp),intent(in)::z,g(2);integer,intent(in),optional::nmax
        complex(dp)::v
        real(dp),allocatable::a(:,:)
        integer::u,m,n
        u=8;if(present(nmax))u=nmax;a=amn_coefficients(u);v=(0.0_dp,0.0_dp)
        do n=0,u-1;do m=0,u-1
            if(abs(a(n+1,m+1))<=tiny(1.0_dp))cycle
            v=v+a(n+1,m+1)*(g(1)/2.0_dp)**m*(2.0_dp*g(2))**n* &
                z**(4*m+6*n)/gamma(real(4*m+6*n+1,dp))
        end do;end do
    end function sigma_prime_laurent

    function e18_10_9(p) result(v)
        type(elliptic_parameters),intent(in)::p
        complex(dp)::v(3),q
        q=exp(pi*ci*p%omega(2)/p%omega(1))
        v(1)=12.0_dp*p%omega(1)**2*p%e(1)-pi*pi*(theta3_q((0.0_dp,0.0_dp),q)**4+theta4_q((0.0_dp,0.0_dp),q)**4)
        v(2)=12.0_dp*p%omega(1)**2*p%e(2)-pi*pi*(theta2_q((0.0_dp,0.0_dp),q)**4-theta4_q((0.0_dp,0.0_dp),q)**4)
        v(3)=12.0_dp*p%omega(1)**2*p%e(3)+pi*pi*(theta2_q((0.0_dp,0.0_dp),q)**4+theta3_q((0.0_dp,0.0_dp),q)**4)
    end function e18_10_9

    function equianharmonic_parameters() result(p)
        type(elliptic_parameters)::p
        complex(dp)::eps,eta,etad
        real(dp)::jj
        jj=gamma(1.0_dp/3.0_dp)**3/(4.0_dp*pi)
        p%omega=[cmplx(jj/2.0_dp,-jj*sqrt(3.0_dp)/2.0_dp,dp),cmplx(jj/2.0_dp,jj*sqrt(3.0_dp)/2.0_dp,dp)]
        eps=exp(pi*ci/3.0_dp);p%e=[4.0_dp**(-1.0_dp/3.0_dp)*eps**2,cmplx(4.0_dp**(-1.0_dp/3.0_dp),0.0_dp,dp), &
            4.0_dp**(-1.0_dp/3.0_dp)*eps**(-2)]
        eta=eps*pi/(2.0_dp*p%omega(2)*sqrt(3.0_dp));etad=-eps**(-1)*pi/(2.0_dp*p%omega(2)*sqrt(3.0_dp))
        p%eta=[etad,eta-etad,-eta];p%g=[(0.0_dp,0.0_dp),(1.0_dp,0.0_dp)];p%delta=(-27.0_dp,0.0_dp)
        p%q=exp(pi*ci*p%omega(2)/p%omega(1))
    end function equianharmonic_parameters

    function lemniscatic_parameters() result(p)
        type(elliptic_parameters)::p
        real(dp)::o,jj
        o=gamma(0.25_dp)**2/(4.0_dp*sqrt(pi));p%omega=[cmplx(o,0.0_dp,dp),cmplx(0.0_dp,o,dp)]
        p%e=[(0.5_dp,0.0_dp),(0.0_dp,0.0_dp),(-0.5_dp,0.0_dp)];jj=pi/(4.0_dp*o)
        p%eta=[cmplx(jj,0.0_dp,dp),cmplx(0.0_dp,-jj,dp),cmplx(-jj,jj,dp)]
        p%g=[(1.0_dp,0.0_dp),(0.0_dp,0.0_dp)];p%delta=(1.0_dp,0.0_dp);p%q=exp(pi*ci*p%omega(2)/p%omega(1))
    end function lemniscatic_parameters

    function pseudolemniscatic_parameters() result(p)
        type(elliptic_parameters)::p
        real(dp)::jj
        complex(dp)::ee
        jj=gamma(0.25_dp)**2/(4.0_dp*sqrt(2.0_dp*pi))
        p%omega=[jj*(1.0_dp-ci),jj*(1.0_dp+ci)];ee=0.5_dp*ci;p%e=[ee,(0.0_dp,0.0_dp),-ee]
        ee=pi/(4.0_dp*p%omega(1));p%eta=[ee,-ee*ci,ee*(ci-1.0_dp)]
        p%g=[(-1.0_dp,0.0_dp),(0.0_dp,0.0_dp)];p%delta=(1.0_dp,0.0_dp);p%q=exp(pi*ci*p%omega(2)/p%omega(1))
    end function pseudolemniscatic_parameters
end module elliptic_weierstrass
