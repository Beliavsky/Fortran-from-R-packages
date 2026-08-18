! SPDX-License-Identifier: GPL-2.0-only
module elliptic_modular
    use elliptic_kinds, only : dp, pi, ci
    use elliptic_theta, only : theta2_q, theta3_q, theta4_q
    use elliptic_weierstrass, only : g2_from_periods, g3_from_periods
    use elliptic_arithmetic, only : divisor_sigma
    implicit none
    private
    public :: dedekind_eta, eta_series, j_invariant, modular_lambda
    public :: mobius_transform, g2_direct, g3_direct, g2_divisor, g3_divisor
    public :: g2_fixed, g3_fixed
    public :: g2_lambert, g3_lambert, farey_series, unimodular_matrices
    public :: lattice_points
contains
    function dedekind_eta(z) result(v)
        complex(dp),intent(in)::z
        complex(dp)::v,q
        q=exp(3.0_dp*pi*ci*z)
        v=theta3_q(pi*(0.5_dp+0.5_dp*z),q)*exp(pi*ci*z/12.0_dp)
    end function dedekind_eta

    function eta_series(z,maxiter) result(v)
        complex(dp),intent(in)::z
        integer,intent(in),optional::maxiter
        complex(dp)::v,term
        integer::n,nmax
        nmax=300;if(present(maxiter))nmax=maxiter
        v=(1.0_dp,0.0_dp)
        do n=1,nmax
            term=exp(2.0_dp*pi*ci*real(n,dp)*z)
            v=v*(1.0_dp-term)
            if(abs(term)<epsilon(1.0_dp))exit
        end do
        v=exp(pi*ci*z/12.0_dp)*v
    end function eta_series

    function j_invariant(tau) result(v)
        complex(dp),intent(in)::tau
        complex(dp)::v,q,t2,t3,t4
        q=exp(pi*ci*tau)
        t2=theta2_q((0.0_dp,0.0_dp),q)
        t3=theta3_q((0.0_dp,0.0_dp),q)
        t4=theta4_q((0.0_dp,0.0_dp),q)
        v=(t2**8+t3**8+t4**8)**3/(t2*t3*t4)**8/54.0_dp
    end function j_invariant

    function modular_lambda(tau) result(v)
        complex(dp),intent(in)::tau
        complex(dp)::v,q
        q=exp(pi*ci*tau)
        v=(theta2_q((0.0_dp,0.0_dp),q)/theta3_q((0.0_dp,0.0_dp),q))**4
    end function modular_lambda

    pure function mobius_transform(m,x) result(v)
        integer,intent(in)::m(2,2)
        complex(dp),intent(in)::x
        complex(dp)::v
        ! R matrices are column-major; elliptic::mob uses c(a,c,b,d).
        v=(real(m(1,1),dp)*x+real(m(1,2),dp))/ &
          (real(m(2,1),dp)*x+real(m(2,2),dp))
    end function mobius_transform

    function g2_direct(b,nmax) result(v)
        complex(dp),intent(in)::b(2)
        integer,intent(in),optional::nmax
        complex(dp)::v,z
        integer::n,i,j
        n=50;if(present(nmax))n=nmax;v=(0.0_dp,0.0_dp)
        do i=-n,n
            do j=-n,n
                if(i==0.and.j==0)cycle
                z=2.0_dp*(real(i,dp)*b(1)+real(j,dp)*b(2))
                v=v+60.0_dp/z**4
            end do
        end do
    end function g2_direct

    function g3_direct(b,nmax) result(v)
        complex(dp),intent(in)::b(2)
        integer,intent(in),optional::nmax
        complex(dp)::v,z
        integer::n,i,j
        n=50;if(present(nmax))n=nmax;v=(0.0_dp,0.0_dp)
        do i=-n,n
            do j=-n,n
                if(i==0.and.j==0)cycle
                z=2.0_dp*(real(i,dp)*b(1)+real(j,dp)*b(2))
                v=v+140.0_dp/z**6
            end do
        end do
    end function g3_direct

    function g2_divisor(b,nmax) result(v)
        complex(dp),intent(in)::b(2)
        integer,intent(in),optional::nmax
        complex(dp)::v,q,s,term
        integer::n,nm
        nm=100;if(present(nmax))nm=nmax
        q=exp(pi*ci*b(2)/b(1));s=(0.0_dp,0.0_dp)
        do n=1,nm
            term=real(divisor_sigma(int(n,8),3),dp)*q**(2*n)
            s=s+term
            if(abs(term)<1e-15_dp*max(1.0_dp,abs(s)).and.n>4)exit
        end do
        v=(pi/b(1))**4*(1.0_dp/12.0_dp+20.0_dp*s)
    end function g2_divisor

    function g3_divisor(b,nmax) result(v)
        complex(dp),intent(in)::b(2)
        integer,intent(in),optional::nmax
        complex(dp)::v,q,s,term
        integer::n,nm
        nm=100;if(present(nmax))nm=nmax
        q=exp(pi*ci*b(2)/b(1));s=(0.0_dp,0.0_dp)
        do n=1,nm
            term=real(divisor_sigma(int(n,8),5),dp)*q**(2*n)
            s=s+term
            if(abs(term)<1e-15_dp*max(1.0_dp,abs(s)).and.n>4)exit
        end do
        v=(pi/b(1))**6*(1.0_dp/216.0_dp-(7.0_dp/3.0_dp)*s)
    end function g3_divisor


    function g2_fixed(b,nmax) result(v)
        complex(dp),intent(in)::b(2)
        integer,intent(in),optional::nmax
        complex(dp)::v,q,s,qn
        integer::n,nm
        nm=50;if(present(nmax))nm=nmax
        q=exp(2.0_dp*pi*ci*b(2)/b(1));s=(0.0_dp,0.0_dp);qn=q
        do n=1,nm
            s=s+real(n,dp)**3*qn/(1.0_dp-qn)
            qn=qn*q
        end do
        v=(pi/b(1))**4*(1.0_dp/12.0_dp+20.0_dp*s)
    end function g2_fixed

    function g3_fixed(b,nmax) result(v)
        complex(dp),intent(in)::b(2)
        integer,intent(in),optional::nmax
        complex(dp)::v,q,s,qn
        integer::n,nm
        nm=50;if(present(nmax))nm=nmax
        q=exp(2.0_dp*pi*ci*b(2)/b(1));s=(0.0_dp,0.0_dp);qn=q
        do n=1,nm
            s=s+real(n,dp)**5*qn/(1.0_dp-qn)
            qn=qn*q
        end do
        v=(pi/b(1))**6*(1.0_dp/216.0_dp-(7.0_dp/3.0_dp)*s)
    end function g3_fixed

    function g2_lambert(b,nmax) result(v)
        complex(dp),intent(in)::b(2)
        integer,intent(in),optional::nmax
        complex(dp)::v,q,s,qn,term
        integer::n,nm
        nm=100;if(present(nmax))nm=nmax
        q=exp(2.0_dp*pi*ci*b(2)/b(1));qn=q;s=(0.0_dp,0.0_dp)
        do n=1,nm
            term=real(n,dp)**3*qn/(1.0_dp-qn);s=s+term
            if(abs(term)<1e-15_dp*max(1.0_dp,abs(s)).and.n>4)exit
            qn=qn*q
        end do
        v=(pi/b(1))**4*(1.0_dp/12.0_dp+20.0_dp*s)
    end function g2_lambert

    function g3_lambert(b,nmax) result(v)
        complex(dp),intent(in)::b(2)
        integer,intent(in),optional::nmax
        complex(dp)::v,q,s,qn,term
        integer::n,nm
        nm=100;if(present(nmax))nm=nmax
        q=exp(2.0_dp*pi*ci*b(2)/b(1));qn=q;s=(0.0_dp,0.0_dp)
        do n=1,nm
            term=real(n,dp)**5*qn/(1.0_dp-qn);s=s+term
            if(abs(term)<1e-15_dp*max(1.0_dp,abs(s)).and.n>4)exit
            qn=qn*q
        end do
        v=(pi/b(1))**6*(1.0_dp/216.0_dp-(7.0_dp/3.0_dp)*s)
    end function g3_lambert

    subroutine farey_series(n,num,den)
        integer,intent(in)::n
        integer,allocatable,intent(out)::num(:),den(:)
        integer,allocatable::tn(:),td(:)
        integer::a,b,k,i,j,g
        if(n<=0)then;allocate(num(2),den(2));num=[0,1];den=[1,0];return;end if
        allocate(tn((n+1)*(n+1)),td((n+1)*(n+1)));k=0
        do b=1,n
            do a=0,b
                g=igcd(a,b)
                if(g==1)then;k=k+1;tn(k)=a;td(k)=b;end if
            end do
        end do
        do i=1,k-1
            do j=i+1,k
                if(real(tn(j),dp)/td(j)<real(tn(i),dp)/td(i))then
                    a=tn(i);tn(i)=tn(j);tn(j)=a;b=td(i);td(i)=td(j);td(j)=b
                end if
            end do
        end do
        allocate(num(k),den(k));num=tn(:k);den=td(:k)
    end subroutine farey_series

    pure integer function igcd(a,b) result(g)
        integer,intent(in)::a,b
        integer::x,y,t
        x=abs(a);y=abs(b)
        if(x==0)then;g=max(y,1);return;end if
        do while(y/=0);t=mod(x,y);x=y;y=t;end do;g=x
    end function igcd

    subroutine unimodular_matrices(n,mats)
        integer,intent(in)::n
        integer,allocatable,intent(out)::mats(:,:,:)
        integer,allocatable::num(:),den(:),tmp(:,:,:)
        integer::ord,i,k
        if(n<=1)then
            allocate(mats(2,2,1))
            mats(:,:,1)=reshape([1,0,0,1],[2,2])
            return
        end if
        allocate(tmp(2,2,max(1,4*n*n+4)))
        k=1
        tmp(:,:,1)=reshape([1,0,0,1],[2,2])
        do ord=1,n
            call farey_series(ord,num,den)
            do i=1,size(num)-1
                if(num(i)*den(i+1)-num(i+1)*den(i)==-1)then
                    if(.not.matrix_exists(tmp,k,den(i),den(i+1),num(i),num(i+1)))then
                        k=k+1
                        tmp(:,:,k)=reshape([den(i),num(i),den(i+1),num(i+1)],[2,2])
                    end if
                end if
            end do
        end do
        allocate(mats(2,2,k));if(k>0)mats=tmp(:,:,:k)
    end subroutine unimodular_matrices

    logical function matrix_exists(a,k,a11,a12,a21,a22) result(exists)
        integer,intent(in)::a(:,:,:),k,a11,a12,a21,a22
        integer::j
        exists=.false.
        do j=1,k
            if(all(a(:,:,j)==reshape([a11,a21,a12,a22],[2,2])))then;exists=.true.;return;end if
        end do
    end function matrix_exists

    function lattice_points(p,n) result(z)
        complex(dp),intent(in)::p(2)
        integer,intent(in)::n
        complex(dp),allocatable::z(:,:)
        integer::i,j
        allocate(z(2*n+1,2*n+1))
        do i=-n,n;do j=-n,n;z(i+n+1,j+n+1)=i*p(1)+j*p(2);end do;end do
    end function lattice_points
end module elliptic_modular
