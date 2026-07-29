! SPDX-License-Identifier: GPL-2.0-or-later
module actuar_supplements
  use actuar_kinds, only : dp
  use actuar_special, only : nan_dp, normal_cdf, regularized_gamma_p, regularized_gamma_q, regularized_beta
  implicit none
  private
  public :: mexp, levexp, mgfexp, mnorm, mgfnorm
  public :: mbeta, levbeta, mgamma, levgamma, mgfgamma
  public :: mweibull, levweibull, mlnorm, levlnorm
  public :: munif, levunif, mgfunif, mchisq, levchisq, mgfchisq
contains

  pure function mexp(order,scale) result(m)
    real(dp), intent(in) :: order,scale
    real(dp) :: m
    if(scale<=0.0_dp) then; m=nan_dp()
    else if(order<=-1.0_dp) then; m=huge(1.0_dp)
    else; m=scale**order*gamma(1.0_dp+order)
    end if
  end function mexp

  pure function levexp(limit,scale,order) result(m)
    real(dp), intent(in) :: limit,scale,order
    real(dp) :: m,u,a
    if(scale<=0.0_dp) then; m=nan_dp()
    else if(limit<=0.0_dp) then; m=0.0_dp
    else if(order<=-1.0_dp) then; m=huge(1.0_dp)
    else
      a=1.0_dp+order; u=limit/scale
      m=scale**order*gamma(a)*regularized_gamma_p(a,u)+limit**order*exp(-u)
    end if
  end function levexp

  pure function mgfexp(t,scale) result(m)
    real(dp), intent(in) :: t,scale
    real(dp) :: m
    if(scale<=0.0_dp .or. scale*t>=1.0_dp) then; m=nan_dp()
    else; m=1.0_dp/(1.0_dp-scale*t)
    end if
  end function mgfexp

  pure function mnorm(order,mean,sd) result(m)
    integer, intent(in) :: order
    real(dp), intent(in) :: mean,sd
    real(dp) :: m
    integer :: i
    if(sd<=0.0_dp .or. order<0) then; m=nan_dp(); return; end if
    m=0.0_dp
    do i=0,order/2
      m=m+gamma(real(order+1,dp))*sd**(2*i)*mean**(order-2*i)/ &
        (2.0_dp**i*gamma(real(i+1,dp))*gamma(real(order-2*i+1,dp)))
    end do
  end function mnorm

  pure function mgfnorm(t,mean,sd) result(m)
    real(dp), intent(in) :: t,mean,sd
    real(dp) :: m
    if(sd<=0.0_dp) then; m=nan_dp()
    else; m=exp(t*mean+0.5_dp*t*t*sd*sd)
    end if
  end function mgfnorm

  pure function mbeta(order,shape1,shape2) result(m)
    real(dp), intent(in) :: order,shape1,shape2
    real(dp) :: m
    if(min(shape1,shape2)<=0.0_dp) then; m=nan_dp()
    else if(order<=-shape1) then; m=huge(1.0_dp)
    else
      m=exp(log_gamma(shape1+order)+log_gamma(shape1+shape2)- &
        log_gamma(shape1)-log_gamma(shape1+shape2+order))
    end if
  end function mbeta

  pure function levbeta(limit,shape1,shape2,order) result(m)
    real(dp), intent(in) :: limit,shape1,shape2,order
    real(dp) :: m,full
    if(min(shape1,shape2)<=0.0_dp) then; m=nan_dp()
    else if(limit<=0.0_dp) then; m=0.0_dp
    else if(order<=-shape1) then; m=huge(1.0_dp)
    else if(limit>=1.0_dp) then; m=mbeta(order,shape1,shape2)
    else
      full=mbeta(order,shape1,shape2)
      m=full*regularized_beta(limit,shape1+order,shape2)+ &
        limit**order*(1.0_dp-regularized_beta(limit,shape1,shape2))
    end if
  end function levbeta

  pure function mgamma(order,shape,scale) result(m)
    real(dp), intent(in) :: order,shape,scale
    real(dp) :: m
    if(shape<=0.0_dp .or. scale<=0.0_dp) then; m=nan_dp()
    else if(order<=-shape) then; m=huge(1.0_dp)
    else; m=scale**order*gamma(shape+order)/gamma(shape)
    end if
  end function mgamma

  pure function levgamma(limit,shape,scale,order) result(m)
    real(dp), intent(in) :: limit,shape,scale,order
    real(dp) :: m,a,u
    if(shape<=0.0_dp .or. scale<=0.0_dp) then; m=nan_dp()
    else if(limit<=0.0_dp) then; m=0.0_dp
    else if(order<=-shape) then; m=huge(1.0_dp)
    else
      a=shape+order; u=limit/scale
      m=scale**order*gamma(a)/gamma(shape)*regularized_gamma_p(a,u)+ &
        limit**order*regularized_gamma_q(shape,u)
    end if
  end function levgamma

  pure function mgfgamma(t,shape,scale) result(m)
    real(dp), intent(in) :: t,shape,scale
    real(dp) :: m
    if(shape<=0.0_dp .or. scale<=0.0_dp .or. scale*t>=1.0_dp) then; m=nan_dp()
    else; m=(1.0_dp-scale*t)**(-shape)
    end if
  end function mgfgamma

  pure function mweibull(order,shape,scale) result(m)
    real(dp), intent(in) :: order,shape,scale
    real(dp) :: m
    if(shape<=0.0_dp .or. scale<=0.0_dp) then; m=nan_dp()
    else if(order<=-shape) then; m=huge(1.0_dp)
    else; m=scale**order*gamma(1.0_dp+order/shape)
    end if
  end function mweibull

  pure function levweibull(limit,shape,scale,order) result(m)
    real(dp), intent(in) :: limit,shape,scale,order
    real(dp) :: m,u,a
    if(shape<=0.0_dp .or. scale<=0.0_dp) then; m=nan_dp()
    else if(limit<=0.0_dp) then; m=0.0_dp
    else if(order<=-shape) then; m=huge(1.0_dp)
    else
      a=1.0_dp+order/shape; u=(limit/scale)**shape
      m=scale**order*gamma(a)*regularized_gamma_p(a,u)+limit**order*exp(-u)
    end if
  end function levweibull

  pure function mlnorm(order,logmean,logsd) result(m)
    real(dp), intent(in) :: order,logmean,logsd
    real(dp) :: m
    if(logsd<=0.0_dp) then; m=nan_dp()
    else; m=exp(order*(logmean+0.5_dp*order*logsd**2))
    end if
  end function mlnorm

  pure function levlnorm(limit,logmean,logsd,order) result(m)
    real(dp), intent(in) :: limit,logmean,logsd,order
    real(dp) :: m,u
    if(logsd<=0.0_dp) then; m=nan_dp()
    else if(limit<=0.0_dp) then; m=0.0_dp
    else
      u=(log(limit)-logmean)/logsd
      m=mlnorm(order,logmean,logsd)*normal_cdf(u-order*logsd)+ &
        limit**order*(1.0_dp-normal_cdf(u))
    end if
  end function levlnorm

  pure function munif(order,xmin,xmax) result(m)
    real(dp), intent(in) :: order,xmin,xmax
    real(dp) :: m,a
    if(xmin>=xmax) then; m=nan_dp()
    else if(abs(order+1.0_dp)<1.0e-14_dp) then
      m=(log(abs(xmax))-log(abs(xmin)))/(xmax-xmin)
    else
      a=order+1.0_dp; m=(xmax**a-xmin**a)/((xmax-xmin)*a)
    end if
  end function munif

  pure function levunif(limit,xmin,xmax,order) result(m)
    real(dp), intent(in) :: limit,xmin,xmax,order
    real(dp) :: m,a
    if(xmin>=xmax) then; m=nan_dp()
    else if(limit<=xmin) then; m=limit**order
    else if(limit>=xmax) then; m=munif(order,xmin,xmax)
    else if(abs(order+1.0_dp)<1.0e-14_dp) then
      m=(log(abs(limit))-log(abs(xmin)))/(xmax-xmin)+ &
        (xmax-limit)/(limit*(xmax-xmin))
    else
      a=order+1.0_dp
      m=(limit**a-xmin**a)/((xmax-xmin)*a)+ &
        limit**order*(xmax-limit)/(xmax-xmin)
    end if
  end function levunif

  pure function mgfunif(t,xmin,xmax) result(m)
    real(dp), intent(in) :: t,xmin,xmax
    real(dp) :: m
    if(xmin>=xmax) then; m=nan_dp()
    else if(abs(t)<1.0e-14_dp) then; m=1.0_dp
    else; m=(exp(t*xmax)-exp(t*xmin))/(t*(xmax-xmin))
    end if
  end function mgfunif

  pure function mchisq(order,df,ncp) result(m)
    integer, intent(in) :: order
    real(dp), intent(in) :: df,ncp
    real(dp) :: m
    real(dp), allocatable :: moments(:)
    integer :: i,j
    if(df<=0.0_dp .or. ncp<0.0_dp .or. order<0) then; m=nan_dp(); return; end if
    if(order==0) then; m=1.0_dp; return; end if
    if(ncp==0.0_dp) then
      m=2.0_dp**order*gamma(real(order,dp)+df/2.0_dp)/gamma(df/2.0_dp)
      return
    end if
    allocate(moments(0:order)); moments=0.0_dp; moments(0)=1.0_dp
    moments(1)=df+ncp
    do i=2,order
      moments(i)=2.0_dp**(i-1)*(df+real(i,dp)*ncp)
      do j=1,i-1
        moments(i)=moments(i)+2.0_dp**(j-1)*(df+real(j,dp)*ncp)* &
          moments(i-j)/gamma(real(i-j+1,dp))
      end do
      moments(i)=moments(i)*gamma(real(i,dp))
    end do
    m=moments(order)
  end function mchisq

  pure function levchisq(limit,df,ncp,order) result(m)
    real(dp), intent(in) :: limit,df,ncp,order
    real(dp) :: m,a,u
    if(df<=0.0_dp .or. ncp<0.0_dp) then; m=nan_dp()
    else if(limit<=0.0_dp) then; m=0.0_dp
    else if(ncp>0.0_dp) then; m=nan_dp()
    else if(order<=-df/2.0_dp) then; m=huge(1.0_dp)
    else
      a=order+df/2.0_dp; u=limit/2.0_dp
      m=2.0_dp**order*gamma(a)/gamma(df/2.0_dp)*regularized_gamma_p(a,u)+ &
        limit**order*regularized_gamma_q(df/2.0_dp,u)
    end if
  end function levchisq

  pure function mgfchisq(t,df,ncp) result(m)
    real(dp), intent(in) :: t,df,ncp
    real(dp) :: m
    if(df<=0.0_dp .or. ncp<0.0_dp .or. 2.0_dp*t>=1.0_dp) then; m=nan_dp()
    else; m=exp(ncp*t/(1.0_dp-2.0_dp*t))*(1.0_dp-2.0_dp*t)**(-df/2.0_dp)
    end if
  end function mgfchisq

end module actuar_supplements
