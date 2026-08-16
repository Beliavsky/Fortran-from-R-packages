! SPDX-License-Identifier: LGPL-3.0-or-later
! Cornish-Fisher approximation following PDQutils R/cornish_fisher.r.
module pdqutils_cf
    use pdqutils_kinds, only : dp
    use pdqutils_special, only : normal_quantile, nan_dp, pos_inf_dp, neg_inf_dp
    implicit none
    private
    public :: as269, as269_orders, as269_vector
    public :: qapx_cf, qapx_cf_vec, rapx_cf

contains

    pure function as269_orders(z,cumul,order_max) result(retval)
        real(dp), intent(in) :: z,cumul(:)
        integer, intent(in), optional :: order_max
        real(dp), allocatable :: retval(:)
        real(dp), allocatable :: a(:),h(:),pp(:),d(:),delta(:)
        real(dp) :: fac,bc,dd,aa,running
        integer :: om,nord,j,jj,kk,ll,ja,jb,jbl,jal
        om=size(cumul)+2
        if(present(order_max)) om=order_max
        if(om<=2) then
            allocate(retval(1)); retval=z; return
        end if
        nord=om-2
        if(nord>size(cumul)) then
            allocate(retval(1)); retval=nan_dp(); return
        end if
        allocate(retval(nord),a(nord),h(3*nord),pp(3*nord*(nord+1)/2),d(3*nord),delta(nord))
        pp=0.0_dp; d=0.0_dp; delta=0.0_dp
        do j=1,nord
            a(j)=(-1.0_dp)**j*cumul(j)/(real(j+1,dp)*real(j+2,dp))
        end do
        h(1)=-z
        if(size(h)>=2) h(2)=z*z-1.0_dp
        do j=3,size(h)
            h(j)=-(z*h(j-1)+real(j-1,dp)*h(j-2))
        end do
        d(1)=-a(1)*h(2)
        delta(1)=d(1)
        pp(1)=d(1)
        if(size(pp)>=3) pp(3)=a(1)
        ja=0; fac=1.0_dp
        if(nord>1) then
            do jj=2,nord
                fac=fac*real(jj,dp)
                ja=ja+3*(jj-1)
                jb=ja; bc=1.0_dp
                do kk=1,jj-1
                    dd=bc*d(kk); aa=bc*a(kk)
                    jb=jb-3*(jj-kk)
                    do ll=1,3*(jj-kk)
                        jbl=jb+ll; jal=ja+ll
                        pp(jal+1)=pp(jal+1)+dd*pp(jbl)
                        pp(jal+kk+2)=pp(jal+kk+2)+aa*pp(jbl)
                    end do
                    bc=bc*real(jj-kk,dp)/real(kk,dp)
                end do
                pp(ja+jj+2)=pp(ja+jj+2)+a(jj)
                d(jj)=0.0_dp
                do ll=2,3*jj
                    d(jj)=d(jj)-pp(ja+ll)*h(ll-1)
                end do
                pp(ja+1)=d(jj)
                delta(jj)=d(jj)/fac
            end do
        end if
        running=0.0_dp
        do j=1,nord
            running=running+delta(j)
            retval(j)=z+running
        end do
    end function as269_orders

    pure real(dp) function as269(z,cumul,order_max) result(retval)
        real(dp), intent(in) :: z,cumul(:)
        integer, intent(in), optional :: order_max
        real(dp), allocatable :: vals(:)
        vals=as269_orders(z,cumul,order_max)
        retval=vals(size(vals))
    end function as269

    pure function as269_vector(z,cumul,order_max) result(retval)
        real(dp), intent(in) :: z(:),cumul(:)
        integer, intent(in), optional :: order_max
        real(dp), allocatable :: retval(:)
        integer :: i
        allocate(retval(size(z)))
        do i=1,size(z)
            retval(i)=as269(z(i),cumul,order_max)
        end do
    end function as269_vector

    pure real(dp) function qapx_cf(p,raw_cumulants,support_lo,support_hi,lower_tail,log_p) result(v)
        real(dp), intent(in) :: p,raw_cumulants(:)
        real(dp), intent(in), optional :: support_lo,support_hi
        logical, intent(in), optional :: lower_tail,log_p
        real(dp) :: z,lo,hi
        real(dp), allocatable :: gammas(:)
        logical :: lower,lp
        integer :: n,j
        lower=.true.; lp=.false.; lo=neg_inf_dp(); hi=pos_inf_dp()
        if(present(lower_tail)) lower=lower_tail
        if(present(log_p)) lp=log_p
        if(present(support_lo)) lo=support_lo
        if(present(support_hi)) hi=support_hi
        n=size(raw_cumulants)
        if(n<2 .or. raw_cumulants(2)<=0.0_dp) then
            v=nan_dp(); return
        end if
        ! Robust endpoint handling.  The upstream R code lets qnorm return an
        ! infinity, after which high-order polynomial corrections can become NaN.
        if(lp) then
            if(p>0.0_dp) then
                v=nan_dp(); return
            else if(p==0.0_dp) then
                v=merge(hi,lo,lower); return
            else if(p==neg_inf_dp()) then
                v=merge(lo,hi,lower); return
            end if
        else
            if(p<0.0_dp .or. p>1.0_dp) then
                v=nan_dp(); return
            else if(p==0.0_dp) then
                v=merge(lo,hi,lower); return
            else if(p==1.0_dp) then
                v=merge(hi,lo,lower); return
            end if
        end if
        z=normal_quantile(p,lower,lp)
        if(n>2) then
            allocate(gammas(n-2))
            do j=3,n
                gammas(j-2)=raw_cumulants(j)/(raw_cumulants(2)**(0.5_dp*real(j,dp)))
            end do
            z=as269(z,gammas)
        end if
        v=raw_cumulants(1)+sqrt(raw_cumulants(2))*z
        v=max(lo,min(hi,v))
    end function qapx_cf

    pure function qapx_cf_vec(p,raw_cumulants,support_lo,support_hi,lower_tail,log_p) result(v)
        real(dp), intent(in) :: p(:),raw_cumulants(:)
        real(dp), intent(in), optional :: support_lo,support_hi
        logical, intent(in), optional :: lower_tail,log_p
        real(dp), allocatable :: v(:)
        integer :: i
        allocate(v(size(p)))
        do i=1,size(p)
            v(i)=qapx_cf(p(i),raw_cumulants,support_lo,support_hi,lower_tail,log_p)
        end do
    end function qapx_cf_vec

    subroutine rapx_cf(raw_cumulants,x,support_lo,support_hi)
        real(dp), intent(in) :: raw_cumulants(:)
        real(dp), intent(out) :: x(:)
        real(dp), intent(in), optional :: support_lo,support_hi
        real(dp) :: u
        integer :: i
        do i=1,size(x)
            call random_number(u)
            x(i)=qapx_cf(u,raw_cumulants,support_lo,support_hi)
        end do
    end subroutine rapx_cf

end module pdqutils_cf
