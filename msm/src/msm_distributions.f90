! SPDX-License-Identifier: GPL-2.0-or-later
! Distribution helpers translated from R/pexp.R, R/phase.R, R/tnorm.R and R/medists.R.
module msm_distributions
    use msm_kinds, only : dp
    use msm_stats, only : normal_pdf, normal_cdf, rand_uniform, rand_normal, rand_exponential, &
        truncated_normal_pdf, me_truncated_normal_pdf, me_uniform_pdf
    implicit none
    private
    public :: dpexp, ppexp, qpexp, rpexp
    public :: d2phase, p2phase, q2phase, h2phase, r2phase
    public :: dtnorm, ptnorm, qtnorm, rtnorm
contains
    function interval_index(x,t) result(k)
        real(dp), intent(in) :: x,t(:)
        integer :: k,i
        k=1
        do i=2,size(t)
            if(x>=t(i)) k=i
        end do
    end function interval_index

    function dpexp(x,rate,t,log_density) result(v)
        real(dp), intent(in) :: x,rate(:),t(:)
        logical, intent(in), optional :: log_density
        real(dp) :: v,logv
        integer :: k,i
        if(size(rate)/=size(t) .or. size(t)<1) error stop "dpexp: size mismatch"
        if(abs(t(1))>1.0e-12_dp) error stop "dpexp: first change time must be zero"
        if(x<0.0_dp) then; v=0.0_dp; return; end if
        k=interval_index(x,t); logv=log(rate(k))-rate(k)*(x-t(k))
        do i=1,k-1; logv=logv-rate(i)*(t(i+1)-t(i)); end do
        if(present(log_density)) then
            if(log_density) then; v=logv; else; v=exp(logv); end if
        else; v=exp(logv); end if
    end function dpexp

    function ppexp(x,rate,t,lower_tail,log_p) result(v)
        real(dp), intent(in) :: x,rate(:),t(:)
        logical, intent(in), optional :: lower_tail,log_p
        real(dp) :: v,h
        integer :: k,i
        logical :: lt,lp
        lt=.true.; lp=.false.; if(present(lower_tail)) lt=lower_tail; if(present(log_p)) lp=log_p
        if(x<=0.0_dp) then; v=merge(0.0_dp,1.0_dp,lt); if(lp) v=log(max(v,tiny(1.0_dp))); return; end if
        k=interval_index(x,t); h=0.0_dp
        do i=1,k-1; h=h+rate(i)*(t(i+1)-t(i)); end do
        h=h+rate(k)*(x-t(k)); v=1.0_dp-exp(-h)
        if(.not.lt) v=1.0_dp-v
        if(lp) v=log(max(v,tiny(1.0_dp)))
    end function ppexp

    function qpexp(p,rate,t,lower_tail,log_p) result(x)
        real(dp), intent(in) :: p,rate(:),t(:)
        logical, intent(in), optional :: lower_tail,log_p
        real(dp) :: x,pp,h,target,next_h
        logical :: lt,lp
        integer :: k
        lt=.true.; lp=.false.; if(present(lower_tail)) lt=lower_tail; if(present(log_p)) lp=log_p
        pp=p; if(lp) pp=exp(pp); if(.not.lt) pp=1.0_dp-pp
        if(pp<=0.0_dp) then; x=0.0_dp; return; end if
        if(pp>=1.0_dp) then; x=huge(1.0_dp); return; end if
        target=-log(1.0_dp-pp); h=0.0_dp
        do k=1,size(t)-1
            next_h=h+rate(k)*(t(k+1)-t(k))
            if(target<=next_h) then; x=t(k)+(target-h)/rate(k); return; end if
            h=next_h
        end do
        k=size(t); x=t(k)+(target-h)/rate(k)
    end function qpexp

    function rpexp(rate,t,start) result(x)
        real(dp), intent(in) :: rate(:),t(:)
        real(dp), intent(in), optional :: start
        real(dp) :: x,st,u,p0
        st=t(1); if(present(start)) st=start
        if(st<t(1)) error stop "rpexp: start before first change time"
        p0=ppexp(st,rate,t); u=p0+(1.0_dp-p0)*rand_uniform(); x=qpexp(u,rate,t)
    end function rpexp

    pure function d2phase(x,l1,mu1,mu2) result(v)
        real(dp), intent(in) :: x,l1,mu1,mu2
        real(dp) :: v,a
        if(x<0.0_dp .or. min(l1,mu1,mu2)<0.0_dp) then; v=0.0_dp; return; end if
        a=l1+mu1
        if(abs(a-mu2)<=epsilon(1.0_dp)*max(1.0_dp,abs(a))) then
            v=exp(-a*x)*(mu1+a*l1*x)
        else
            v=(-a*exp(-a*x)*(-mu1+mu2)+mu2*l1*exp(-mu2*x))/(a-mu2)
        end if
    end function d2phase

    pure function p2phase(x,l1,mu1,mu2) result(v)
        real(dp), intent(in) :: x,l1,mu1,mu2
        real(dp) :: v,a
        if(x<0.0_dp) then; v=0.0_dp; return; end if
        a=l1+mu1
        if(abs(a-mu2)<=epsilon(1.0_dp)*max(1.0_dp,abs(a))) then
            v=1.0_dp-exp(-a*x)*(1.0_dp+l1*x)
        else
            v=1.0_dp+exp(-a*x)*(-mu1+mu2)/(a-mu2)-l1*exp(-mu2*x)/(a-mu2)
        end if
    end function p2phase

    function q2phase(p,l1,mu1,mu2) result(x)
        real(dp), intent(in) :: p,l1,mu1,mu2
        real(dp) :: x,lo,hi,mid
        integer :: it
        if(p<=0.0_dp) then; x=0.0_dp; return; end if
        if(p>=1.0_dp) then; x=huge(1.0_dp); return; end if
        lo=0.0_dp; hi=max(1.0_dp,10.0_dp/max(min(l1+mu1,mu2),tiny(1.0_dp)))
        do while(p2phase(hi,l1,mu1,mu2)<p); hi=2.0_dp*hi; if(hi>1.0e12_dp) exit; end do
        do it=1,100
            mid=0.5_dp*(lo+hi)
            if(p2phase(mid,l1,mu1,mu2)<p) then; lo=mid; else; hi=mid; end if
        end do
        x=0.5_dp*(lo+hi)
    end function q2phase

    pure function h2phase(x,l1,mu1,mu2) result(v)
        real(dp), intent(in) :: x,l1,mu1,mu2
        real(dp) :: v,s
        s=1.0_dp-p2phase(x,l1,mu1,mu2)
        if(s<=0.0_dp) then; v=huge(1.0_dp); else; v=d2phase(x,l1,mu1,mu2)/s; end if
    end function h2phase

    function r2phase(l1,mu1,mu2) result(x)
        real(dp), intent(in) :: l1,mu1,mu2
        real(dp) :: x,ptrans
        ptrans=l1/(l1+mu1); x=rand_exponential(l1+mu1)
        if(rand_uniform()<ptrans) x=x+rand_exponential(mu2)
    end function r2phase

    pure function dtnorm(x,mean,sd,lower,upper) result(v)
        real(dp), intent(in) :: x,mean,sd,lower,upper
        real(dp) :: v
        v=truncated_normal_pdf(x,mean,sd,lower,upper)
    end function dtnorm

    pure function ptnorm(x,mean,sd,lower,upper) result(v)
        real(dp), intent(in) :: x,mean,sd,lower,upper
        real(dp) :: v,den
        if(x<=lower) then; v=0.0_dp; return; end if
        if(x>=upper) then; v=1.0_dp; return; end if
        den=normal_cdf(upper,mean,sd)-normal_cdf(lower,mean,sd)
        v=(normal_cdf(x,mean,sd)-normal_cdf(lower,mean,sd))/den
    end function ptnorm

    function qtnorm(p,mean,sd,lower,upper) result(x)
        real(dp), intent(in) :: p,mean,sd,lower,upper
        real(dp) :: x,lo,hi,mid
        integer :: it
        if(p<=0.0_dp) then; x=lower; return; end if
        if(p>=1.0_dp) then; x=upper; return; end if
        lo=lower; hi=upper
        if(.not.(lo>-huge(1.0_dp))) lo=mean-12.0_dp*sd
        if(.not.(hi<huge(1.0_dp))) hi=mean+12.0_dp*sd
        do it=1,100
            mid=0.5_dp*(lo+hi)
            if(ptnorm(mid,mean,sd,lower,upper)<p) then; lo=mid; else; hi=mid; end if
        end do
        x=0.5_dp*(lo+hi)
    end function qtnorm

    function rtnorm(mean,sd,lower,upper) result(x)
        real(dp), intent(in) :: mean,sd,lower,upper
        real(dp) :: x
        do
            x=mean+sd*rand_normal(); if(x>=lower .and. x<=upper) exit
        end do
    end function rtnorm
end module msm_distributions
