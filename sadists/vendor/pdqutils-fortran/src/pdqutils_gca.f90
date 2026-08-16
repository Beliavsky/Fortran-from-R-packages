! SPDX-License-Identifier: LGPL-3.0-or-later
! Generalized Gram-Charlier approximations following PDQutils R/gram_charlier.r.
module pdqutils_gca
    use pdqutils_kinds, only : dp
    use pdqutils_special, only : factorial_dp, normal_pdf, normal_cdf, gamma_pdf, gamma_cdf, &
        beta_pdf, beta_cdf, log_beta_dp, nan_dp, pos_inf_dp, neg_inf_dp
    implicit none
    private

    integer, parameter, public :: gca_normal=1, gca_gamma=2, gca_beta=3, gca_arcsine=4, gca_wigner=5
    public :: gca_basis_from_name
    public :: dapx_gca, papx_gca, dapx_gca_vec, papx_gca_vec

    type :: gca_setup_t
        integer :: basis=0
        integer :: order_max=0
        real(dp) :: x=0.0_dp
        real(dp) :: support_lo=0.0_dp
        real(dp) :: support_hi=0.0_dp
        real(dp) :: scalby=1.0_dp
        real(dp) :: shape=0.0_dp
        real(dp) :: shape1=0.0_dp
        real(dp) :: shape2=0.0_dp
        real(dp), allocatable :: moments(:) ! zeroth through order_max
    end type gca_setup_t

contains

    pure integer function gca_basis_from_name(name) result(basis)
        character(len=*), intent(in) :: name
        character(len=len(name)) :: low
        integer :: i,c
        low=name
        do i=1,len(name)
            c=iachar(low(i:i))
            if(c>=iachar('A') .and. c<=iachar('Z')) low(i:i)=achar(c+32)
        end do
        select case(trim(adjustl(low)))
        case('normal'); basis=gca_normal
        case('gamma'); basis=gca_gamma
        case('beta'); basis=gca_beta
        case('arcsine'); basis=gca_arcsine
        case('wigner'); basis=gca_wigner
        case default; basis=0
        end select
    end function gca_basis_from_name

    pure subroutine shift_moments(raw,del,out)
        real(dp), intent(in) :: raw(:),del
        real(dp), intent(out) :: out(size(raw))
        integer :: k,j,n
        real(dp) :: tot,comb
        n=size(raw)-1
        out=raw
        if(n>=1) out(2)=raw(2)+del
        do k=2,n
            tot=0.0_dp
            comb=1.0_dp
            do j=0,k
                if(j>0) comb=comb*real(k-j+1,dp)/real(j,dp)
                tot=tot+comb*del**(k-j)*raw(j+1)
            end do
            out(k+1)=tot
        end do
    end subroutine shift_moments

    pure subroutine scale_moments(raw,scale,out)
        real(dp), intent(in) :: raw(:),scale
        real(dp), intent(out) :: out(size(raw))
        integer :: k
        do k=0,size(raw)-1
            out(k+1)=raw(k+1)*scale**k
        end do
    end subroutine scale_moments

    pure subroutine setup_gca(x,raw_moments,basis_in,support_lo_in,support_hi_in,shape_in,scale_in, &
                              shape1_in,shape2_in,s,ok)
        real(dp), intent(in) :: x,raw_moments(:)
        integer, intent(in), optional :: basis_in
        real(dp), intent(in), optional :: support_lo_in,support_hi_in,shape_in,scale_in,shape1_in,shape2_in
        type(gca_setup_t), intent(out) :: s
        logical, intent(out) :: ok
        integer :: basis,n
        real(dp) :: lo,hi,mu,sigma,varv,theta,kshape,tmu,ts2,tmu2,a,b,mid,half
        real(dp), allocatable :: full(:),tmp(:)

        ok=.false.
        n=size(raw_moments)
        if(n<2) return
        basis=gca_normal
        if(present(basis_in)) basis=basis_in
        if(basis<gca_normal .or. basis>gca_wigner) return

        select case(basis)
        case(gca_normal)
            lo=neg_inf_dp(); hi=pos_inf_dp()
        case(gca_gamma)
            lo=0.0_dp; hi=pos_inf_dp()
        case(gca_beta)
            lo=0.0_dp; hi=1.0_dp
        case(gca_arcsine,gca_wigner)
            lo=-1.0_dp; hi=1.0_dp
        end select
        if(present(support_lo_in)) lo=support_lo_in
        if(present(support_hi_in)) hi=support_hi_in
        if(lo>hi) then
            mu=lo; lo=hi; hi=mu
        end if

        allocate(full(n+1),tmp(n+1),s%moments(n+1))
        full(1)=1.0_dp
        full(2:)=raw_moments
        s%basis=basis; s%order_max=n; s%scalby=1.0_dp

        select case(basis)
        case(gca_normal)
            mu=full(2)
            varv=full(3)-mu*mu
            if(varv<=0.0_dp) return
            sigma=sqrt(varv)
            s%x=(x-mu)/sigma
            call shift_moments(full,-mu,tmp)
            call scale_moments(tmp,1.0_dp/sigma,s%moments)
            s%scalby=1.0_dp/sigma
            s%support_lo=(lo-mu)/sigma
            s%support_hi=(hi-mu)/sigma

        case(gca_gamma)
            if(.not.(lo<pos_inf_dp())) return
            s%x=x-lo
            call shift_moments(full,-lo,s%moments)
            s%support_lo=0.0_dp
            s%support_hi=hi-lo
            if(present(shape_in) .or. present(scale_in)) then
                if(.not.(present(shape_in).and.present(scale_in))) return
                kshape=shape_in; theta=scale_in
            else
                if(s%moments(2)<=0.0_dp) return
                theta=s%moments(3)/s%moments(2)-s%moments(2)
                if(theta<=0.0_dp) return
                kshape=s%moments(2)/theta
            end if
            if(kshape<=0.0_dp .or. theta<=0.0_dp) return
            s%x=s%x/theta
            call scale_moments(s%moments,1.0_dp/theta,tmp)
            s%moments=tmp
            s%scalby=1.0_dp/theta
            s%support_lo=s%support_lo/theta
            s%support_hi=s%support_hi/theta
            s%shape=kshape

        case(gca_beta,gca_arcsine,gca_wigner)
            if(.not.(lo>neg_inf_dp() .and. hi<pos_inf_dp()) .or. hi<=lo) return
            mid=0.5_dp*(hi+lo); half=0.5_dp*(hi-lo)
            s%x=(x-mid)/half
            call shift_moments(full,-mid,tmp)
            call scale_moments(tmp,1.0_dp/half,s%moments)
            s%scalby=1.0_dp/half
            s%support_lo=-1.0_dp; s%support_hi=1.0_dp
            if(basis==gca_arcsine) then
                a=0.5_dp; b=0.5_dp
            else if(basis==gca_wigner) then
                a=1.5_dp; b=1.5_dp
            else if(present(shape1_in) .or. present(shape2_in)) then
                if(.not.(present(shape1_in).and.present(shape2_in))) return
                a=shape1_in; b=shape2_in
            else
                mu=s%moments(2)
                varv=s%moments(3)-mu*mu
                tmu=(mu+1.0_dp)/2.0_dp
                ts2=varv/4.0_dp
                if(tmu<=0.0_dp .or. tmu>=1.0_dp .or. ts2<=0.0_dp) return
                tmu2=ts2+tmu*tmu
                b=(tmu-tmu2)*(1.0_dp-tmu)/ts2
                a=b*tmu/(1.0_dp-tmu)
            end if
            if(a<=0.0_dp .or. b<=0.0_dp) return
            s%shape1=a; s%shape2=b
        end select
        ok=.true.
    end subroutine setup_gca

    pure subroutine hermite_coeffs(n,coef)
        integer, intent(in) :: n
        real(dp), intent(out) :: coef(n+1,n+1)
        integer :: k,j
        coef=0.0_dp; coef(1,1)=1.0_dp
        if(n==0) return
        coef(2,2)=1.0_dp
        do k=1,n-1
            do j=0,k
                coef(k+2,j+2)=coef(k+2,j+2)+coef(k+1,j+1)
                coef(k+2,j+1)=coef(k+2,j+1)-real(k,dp)*coef(k,j+1)
            end do
        end do
    end subroutine hermite_coeffs

    pure subroutine laguerre_coeffs(n,alpha,coef)
        integer, intent(in) :: n
        real(dp), intent(in) :: alpha
        real(dp), intent(out) :: coef(n+1,n+1)
        integer :: k,j
        real(dp) :: den,a0
        coef=0.0_dp; coef(1,1)=1.0_dp
        if(n==0) return
        coef(2,1)=alpha+1.0_dp; coef(2,2)=-1.0_dp
        do k=1,n-1
            den=real(k+1,dp)
            a0=2.0_dp*real(k,dp)+alpha+1.0_dp
            do j=0,k
                coef(k+2,j+1)=coef(k+2,j+1)+a0*coef(k+1,j+1)/den
                coef(k+2,j+2)=coef(k+2,j+2)-coef(k+1,j+1)/den
                coef(k+2,j+1)=coef(k+2,j+1)-(real(k,dp)+alpha)*coef(k,j+1)/den
            end do
        end do
    end subroutine laguerre_coeffs

    pure subroutine jacobi_coeffs(n,alpha,beta,coef)
        integer, intent(in) :: n
        real(dp), intent(in) :: alpha,beta
        real(dp), intent(out) :: coef(n+1,n+1)
        integer :: k,j
        real(dp) :: kk,ab,a0,a1,b0,den
        coef=0.0_dp; coef(1,1)=1.0_dp
        if(n==0) return
        coef(2,1)=0.5_dp*(alpha-beta)
        coef(2,2)=0.5_dp*(alpha+beta+2.0_dp)
        ab=alpha+beta
        do k=1,n-1
            kk=real(k,dp)
            den=2.0_dp*(kk+1.0_dp)*(kk+ab+1.0_dp)*(2.0_dp*kk+ab)
            a0=(2.0_dp*kk+ab+1.0_dp)*(alpha*alpha-beta*beta)
            a1=(2.0_dp*kk+ab+1.0_dp)*(2.0_dp*kk+ab)*(2.0_dp*kk+ab+2.0_dp)
            b0=2.0_dp*(kk+alpha)*(kk+beta)*(2.0_dp*kk+ab+2.0_dp)
            do j=0,k
                coef(k+2,j+1)=coef(k+2,j+1)+a0*coef(k+1,j+1)/den
                coef(k+2,j+2)=coef(k+2,j+2)+a1*coef(k+1,j+1)/den
                coef(k+2,j+1)=coef(k+2,j+1)-b0*coef(k,j+1)/den
            end do
        end do
    end subroutine jacobi_coeffs

    pure real(dp) function poly_eval(row,degree,x) result(v)
        real(dp), intent(in) :: row(:),x
        integer, intent(in) :: degree
        integer :: j
        v=row(degree+1)
        do j=degree,1,-1
            v=v*x+row(j)
        end do
    end function poly_eval

    pure real(dp) function poly_expect(row,degree,moments) result(v)
        real(dp), intent(in) :: row(:),moments(:)
        integer, intent(in) :: degree
        integer :: j
        v=0.0_dp
        do j=0,degree
            v=v+row(j+1)*moments(j+1)
        end do
    end function poly_expect

    pure real(dp) function gca_hn(n,s) result(h)
        integer, intent(in) :: n
        type(gca_setup_t), intent(in) :: s
        real(dp) :: alpha,beta,rn
        if(n==0) then
            h=1.0_dp
            return
        end if
        rn=real(n,dp)
        select case(s%basis)
        case(gca_normal)
            h=factorial_dp(n)
        case(gca_gamma)
            alpha=s%shape-1.0_dp
            h=exp(log_gamma(alpha+1.0_dp+rn)-log_gamma(alpha+1.0_dp)-log_gamma(rn+1.0_dp))
        case(gca_beta,gca_arcsine,gca_wigner)
            alpha=s%shape2-1.0_dp; beta=s%shape1-1.0_dp
            h=exp(log_gamma(rn+alpha+1.0_dp)+log_gamma(rn+beta+1.0_dp)-log_gamma(rn+1.0_dp) &
                -log_gamma(rn+alpha+beta+1.0_dp)-log_beta_dp(alpha+1.0_dp,beta+1.0_dp) &
                -log(2.0_dp*rn+alpha+beta+1.0_dp))
        case default
            h=nan_dp()
        end select
    end function gca_hn

    pure real(dp) function gca_weight(y,s) result(w)
        real(dp), intent(in) :: y
        type(gca_setup_t), intent(in) :: s
        select case(s%basis)
        case(gca_normal)
            w=normal_pdf(y)
        case(gca_gamma)
            w=gamma_pdf(y,s%shape)
        case(gca_beta,gca_arcsine,gca_wigner)
            w=0.5_dp*beta_pdf(0.5_dp*(y+1.0_dp),s%shape1,s%shape2)
        case default
            w=nan_dp()
        end select
    end function gca_weight

    pure subroutine make_polynomials(s,coef)
        type(gca_setup_t), intent(in) :: s
        real(dp), intent(out) :: coef(s%order_max+1,s%order_max+1)
        real(dp) :: alpha,beta
        select case(s%basis)
        case(gca_normal)
            call hermite_coeffs(s%order_max,coef)
        case(gca_gamma)
            call laguerre_coeffs(s%order_max,s%shape-1.0_dp,coef)
        case(gca_beta,gca_arcsine,gca_wigner)
            alpha=s%shape2-1.0_dp; beta=s%shape1-1.0_dp
            call jacobi_coeffs(s%order_max,alpha,beta,coef)
        end select
    end subroutine make_polynomials

    pure real(dp) function intpoly(n,y,s) result(v)
        integer, intent(in) :: n
        real(dp), intent(in) :: y
        type(gca_setup_t), intent(in) :: s
        real(dp), allocatable :: c(:,:)
        real(dp) :: alpha,beta,t,ratio
        v=nan_dp()
        if(n==0) then
            select case(s%basis)
            case(gca_normal)
                v=normal_cdf(y)
            case(gca_gamma)
                v=gamma_cdf(y,s%shape)
            case(gca_beta,gca_arcsine,gca_wigner)
                v=beta_cdf(0.5_dp*(y+1.0_dp),s%shape1,s%shape2)
            end select
            return
        end if
        select case(s%basis)
        case(gca_normal)
            allocate(c(n,n))
            call hermite_coeffs(n-1,c)
            v=-normal_pdf(y)*poly_eval(c(n,:),n-1,y)
        case(gca_gamma)
            allocate(c(n,n))
            call laguerre_coeffs(n-1,s%shape,c)
            v=(s%shape/real(n,dp))*gamma_pdf(y,s%shape+1.0_dp)*poly_eval(c(n,:),n-1,y)
        case(gca_beta,gca_arcsine,gca_wigner)
            alpha=s%shape2-1.0_dp; beta=s%shape1-1.0_dp
            allocate(c(n,n))
            call jacobi_coeffs(n-1,alpha+1.0_dp,beta+1.0_dp,c)
            t=0.5_dp*(y+1.0_dp)
            ratio=exp(log_beta_dp(alpha+2.0_dp,beta+2.0_dp)-log_beta_dp(alpha+1.0_dp,beta+1.0_dp))
            ! Mathematically consistent Jacobi/Beta mapping.  The upstream R
            ! closure swaps these two beta-density shape arguments; that typo
            ! is invisible when the expansion exactly equals its beta parent.
            v=-(ratio/real(n,dp))*beta_pdf(t,beta+2.0_dp,alpha+2.0_dp)*poly_eval(c(n,:),n-1,y)
        case default
            v=nan_dp()
        end select
    end function intpoly

    pure real(dp) function dapx_gca(x,raw_moments,basis,support_lo,support_hi,shape,scale,shape1,shape2, &
                                    log_density) result(v)
        real(dp), intent(in) :: x,raw_moments(:)
        integer, intent(in), optional :: basis
        real(dp), intent(in), optional :: support_lo,support_hi,shape,scale,shape1,shape2
        logical, intent(in), optional :: log_density
        type(gca_setup_t) :: s
        real(dp), allocatable :: coef(:,:)
        real(dp) :: wx,ci
        logical :: ok,ld
        integer :: n
        ld=.false.; if(present(log_density)) ld=log_density
        call setup_gca(x,raw_moments,basis,support_lo,support_hi,shape,scale,shape1,shape2,s,ok)
        if(.not.ok) then
            v=nan_dp(); return
        end if
        if(s%x<=s%support_lo .or. s%x>=s%support_hi) then
            v=merge(neg_inf_dp(),0.0_dp,ld); return
        end if
        allocate(coef(s%order_max+1,s%order_max+1)); call make_polynomials(s,coef)
        wx=gca_weight(s%x,s)
        v=0.0_dp
        do n=0,s%order_max
            ci=poly_expect(coef(n+1,:),n,s%moments)/gca_hn(n,s)
            v=v+ci*wx*poly_eval(coef(n+1,:),n,s%x)
        end do
        v=max(0.0_dp,v*s%scalby)
        if(ld) then
            if(v==0.0_dp) then; v=neg_inf_dp(); else; v=log(v); end if
        end if
    end function dapx_gca

    pure real(dp) function papx_gca(q,raw_moments,basis,support_lo,support_hi,shape,scale,shape1,shape2, &
                                    lower_tail,log_p) result(v)
        real(dp), intent(in) :: q,raw_moments(:)
        integer, intent(in), optional :: basis
        real(dp), intent(in), optional :: support_lo,support_hi,shape,scale,shape1,shape2
        logical, intent(in), optional :: lower_tail,log_p
        type(gca_setup_t) :: s
        real(dp), allocatable :: coef(:,:)
        real(dp) :: ci
        logical :: ok,lower,lp
        integer :: n
        lower=.true.; lp=.false.
        if(present(lower_tail)) lower=lower_tail
        if(present(log_p)) lp=log_p
        call setup_gca(q,raw_moments,basis,support_lo,support_hi,shape,scale,shape1,shape2,s,ok)
        if(.not.ok) then
            v=nan_dp(); return
        end if
        if(s%x<=s%support_lo) then
            v=0.0_dp
        else if(s%x>=s%support_hi) then
            v=1.0_dp
        else
            allocate(coef(s%order_max+1,s%order_max+1)); call make_polynomials(s,coef)
            v=0.0_dp
            do n=0,s%order_max
                ci=poly_expect(coef(n+1,:),n,s%moments)/gca_hn(n,s)
                v=v+ci*intpoly(n,s%x,s)
            end do
            v=min(1.0_dp,max(0.0_dp,v))
        end if
        if(.not.lower) v=1.0_dp-v
        if(lp) then
            if(v==0.0_dp) then; v=neg_inf_dp(); else; v=log(v); end if
        end if
    end function papx_gca

    pure function dapx_gca_vec(x,raw_moments,basis,support_lo,support_hi,shape,scale,shape1,shape2, &
                               log_density) result(v)
        real(dp), intent(in) :: x(:),raw_moments(:)
        integer, intent(in), optional :: basis
        real(dp), intent(in), optional :: support_lo,support_hi,shape,scale,shape1,shape2
        logical, intent(in), optional :: log_density
        real(dp), allocatable :: v(:)
        integer :: i
        allocate(v(size(x)))
        do i=1,size(x)
            v(i)=dapx_gca(x(i),raw_moments,basis,support_lo,support_hi,shape,scale,shape1,shape2,log_density)
        end do
    end function dapx_gca_vec

    pure function papx_gca_vec(q,raw_moments,basis,support_lo,support_hi,shape,scale,shape1,shape2, &
                               lower_tail,log_p) result(v)
        real(dp), intent(in) :: q(:),raw_moments(:)
        integer, intent(in), optional :: basis
        real(dp), intent(in), optional :: support_lo,support_hi,shape,scale,shape1,shape2
        logical, intent(in), optional :: lower_tail,log_p
        real(dp), allocatable :: v(:)
        integer :: i
        allocate(v(size(q)))
        do i=1,size(q)
            v(i)=papx_gca(q(i),raw_moments,basis,support_lo,support_hi,shape,scale,shape1,shape2,lower_tail,log_p)
        end do
    end function papx_gca_vec

end module pdqutils_gca
