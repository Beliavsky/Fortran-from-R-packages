! SPDX-License-Identifier: LGPL-3.0-or-later
! Edgeworth expansion following PDQutils R/edgeworth.r.
module pdqutils_edgeworth
    use pdqutils_kinds, only : dp
    use pdqutils_special, only : factorial_dp, normal_pdf, normal_cdf, nan_dp, pos_inf_dp, neg_inf_dp
    implicit none
    private
    public :: dapx_edgeworth, papx_edgeworth
    public :: dapx_edgeworth_vec, papx_edgeworth_vec

contains

    pure real(dp) function hermite_he(n,x) result(v)
        integer, intent(in) :: n
        real(dp), intent(in) :: x
        integer :: j
        real(dp) :: hm2,hm1,h
        if (n < 0) then
            v=nan_dp(); return
        else if (n == 0) then
            v=1.0_dp; return
        else if (n == 1) then
            v=x; return
        end if
        hm2=1.0_dp; hm1=x
        do j=2,n
            h=x*hm1-real(j-1,dp)*hm2
            hm2=hm1; hm1=h
        end do
        v=hm1
    end function hermite_he

    recursive pure subroutine partition_accum(s,idx,remaining,kvals,sn,eta,density,accum)
        integer, intent(in) :: s,idx,remaining
        integer, intent(inout) :: kvals(:)
        real(dp), intent(in) :: sn(:),eta
        logical, intent(in) :: density
        real(dp), intent(inout) :: accum
        integer :: km,m,r,degree
        real(dp) :: coef
        if (idx > s) then
            if (remaining /= 0) return
            r=sum(kvals(1:s))
            coef=1.0_dp
            do m=1,s
                if (kvals(m)>0) then
                    coef=coef*(sn(m)/factorial_dp(m+2))**kvals(m)/factorial_dp(kvals(m))
                end if
            end do
            if (density) then
                degree=s+2*r
                accum=accum+coef*hermite_he(degree,eta)
            else
                degree=s+2*r-1
                accum=accum-coef*hermite_he(degree,eta)
            end if
            return
        end if
        do km=0,remaining/idx
            kvals(idx)=km
            call partition_accum(s,idx+1,remaining-km*idx,kvals,sn,eta,density,accum)
        end do
        kvals(idx)=0
    end subroutine partition_accum

    pure real(dp) function dapx_edgeworth(x,raw_cumulants,support_lo,support_hi,log_density) result(v)
        real(dp), intent(in) :: x,raw_cumulants(:)
        real(dp), intent(in), optional :: support_lo,support_hi
        logical, intent(in), optional :: log_density
        real(dp) :: mu,sigma,eta,phi,lo,hi,nexterm
        real(dp), allocatable :: sn(:)
        integer, allocatable :: kvals(:)
        integer :: order_max,s,m
        logical :: ld
        lo=neg_inf_dp(); hi=pos_inf_dp(); ld=.false.
        if(present(support_lo)) lo=support_lo
        if(present(support_hi)) hi=support_hi
        if(present(log_density)) ld=log_density
        order_max=size(raw_cumulants)
        if(order_max<2 .or. raw_cumulants(2)<=0.0_dp) then
            v=nan_dp(); return
        end if
        if(x<=lo .or. x>=hi) then
            v=merge(neg_inf_dp(),0.0_dp,ld)
            return
        end if
        mu=raw_cumulants(1); sigma=sqrt(raw_cumulants(2))
        eta=(x-mu)/sigma; phi=normal_pdf(eta); v=phi
        if(order_max>2) then
            allocate(sn(order_max-2),kvals(order_max-2))
            do m=1,order_max-2
                sn(m)=raw_cumulants(m+2)/(raw_cumulants(2)**real(m+1,dp))
            end do
            do s=1,order_max-2
                kvals=0; nexterm=0.0_dp
                call partition_accum(s,1,s,kvals,sn,eta,.true.,nexterm)
                v=v+phi*(sigma**s)*nexterm
            end do
        end if
        v=max(0.0_dp,v/sigma)
        if(ld) then
            if(v==0.0_dp) then
                v=neg_inf_dp()
            else
                v=log(v)
            end if
        end if
    end function dapx_edgeworth

    pure real(dp) function papx_edgeworth(q,raw_cumulants,support_lo,support_hi,lower_tail,log_p) result(v)
        real(dp), intent(in) :: q,raw_cumulants(:)
        real(dp), intent(in), optional :: support_lo,support_hi
        logical, intent(in), optional :: lower_tail,log_p
        real(dp) :: qq,mu,sigma,eta,phi,lo,hi,nexterm,tlo,thi
        real(dp), allocatable :: kap(:),sn(:)
        integer, allocatable :: kvals(:)
        integer :: order_max,s,m
        logical :: lower,lp
        lower=.true.; lp=.false.; lo=neg_inf_dp(); hi=pos_inf_dp()
        if(present(lower_tail)) lower=lower_tail
        if(present(log_p)) lp=log_p
        if(present(support_lo)) lo=support_lo
        if(present(support_hi)) hi=support_hi
        order_max=size(raw_cumulants)
        if(order_max<2 .or. raw_cumulants(2)<=0.0_dp) then
            v=nan_dp(); return
        end if
        allocate(kap(order_max)); kap=raw_cumulants; qq=q
        if(.not.lower) then
            qq=-q
            do m=1,order_max
                if(mod(m,2)==1) kap(m)=-kap(m)
            end do
            tlo=-hi; thi=-lo; lo=tlo; hi=thi
        end if
        if(qq<=lo) then
            v=0.0_dp
        else if(qq>=hi) then
            v=1.0_dp
        else
            mu=kap(1); sigma=sqrt(kap(2)); eta=(qq-mu)/sigma
            phi=normal_pdf(eta); v=normal_cdf(eta)
            if(order_max>2) then
                allocate(sn(order_max-2),kvals(order_max-2))
                do m=1,order_max-2
                    sn(m)=kap(m+2)/(kap(2)**real(m+1,dp))
                end do
                do s=1,order_max-2
                    kvals=0; nexterm=0.0_dp
                    call partition_accum(s,1,s,kvals,sn,eta,.false.,nexterm)
                    v=v+phi*(sigma**s)*nexterm
                end do
            end if
            v=min(1.0_dp,max(0.0_dp,v))
        end if
        if(lp) then
            if(v==0.0_dp) then
                v=neg_inf_dp()
            else
                v=log(v)
            end if
        end if
    end function papx_edgeworth

    pure function dapx_edgeworth_vec(x,raw_cumulants,support_lo,support_hi,log_density) result(v)
        real(dp), intent(in) :: x(:),raw_cumulants(:)
        real(dp), intent(in), optional :: support_lo,support_hi
        logical, intent(in), optional :: log_density
        real(dp), allocatable :: v(:)
        integer :: i
        allocate(v(size(x)))
        do i=1,size(x)
            v(i)=dapx_edgeworth(x(i),raw_cumulants,support_lo,support_hi,log_density)
        end do
    end function dapx_edgeworth_vec

    pure function papx_edgeworth_vec(q,raw_cumulants,support_lo,support_hi,lower_tail,log_p) result(v)
        real(dp), intent(in) :: q(:),raw_cumulants(:)
        real(dp), intent(in), optional :: support_lo,support_hi
        logical, intent(in), optional :: lower_tail,log_p
        real(dp), allocatable :: v(:)
        integer :: i
        allocate(v(size(q)))
        do i=1,size(q)
            v(i)=papx_edgeworth(q(i),raw_cumulants,support_lo,support_hi,lower_tail,log_p)
        end do
    end function papx_edgeworth_vec

end module pdqutils_edgeworth
