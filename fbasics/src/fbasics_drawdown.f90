! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from the GPL-2-or-later fBasics package.
module fbasics_drawdown
  use fbasics_kinds, only: dp, pi
  use fbasics_rng, only: rnorm_lcg
  implicit none
  private
  public :: maxdd_expectation, dmaxdd, pmaxdd, rmaxdd, realized_max_drawdown
  real(dp), parameter :: qpx(50) = [ &
      0.0005_dp, 0.001_dp, 0.0015_dp, 0.002_dp, 0.0025_dp, &
      0.005_dp, 0.0075_dp, 0.01_dp, 0.0125_dp, 0.015_dp, &
      0.0175_dp, 0.02_dp, 0.0225_dp, 0.025_dp, 0.0275_dp, &
      0.03_dp, 0.0325_dp, 0.035_dp, 0.0375_dp, 0.04_dp, &
      0.0425_dp, 0.045_dp, 0.05_dp, 0.06_dp, 0.07_dp, &
      0.08_dp, 0.09_dp, 0.1_dp, 0.2_dp, 0.3_dp, &
      0.4_dp, 0.5_dp, 1.5_dp, 2.5_dp, 3.5_dp, &
      4.5_dp, 10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, &
      50.0_dp, 150.0_dp, 250.0_dp, 350.0_dp, 450.0_dp, &
      1000.0_dp, 2000.0_dp, 3000.0_dp, 4000.0_dp, 5000.0_dp ]
  real(dp), parameter :: qpy(50) = [ &
      0.01969_dp, 0.027694_dp, 0.033789_dp, 0.038896_dp, 0.043372_dp, &
      0.060721_dp, 0.073808_dp, 0.084693_dp, 0.094171_dp, 0.102651_dp, &
      0.110375_dp, 0.117503_dp, 0.124142_dp, 0.130374_dp, 0.136259_dp, &
      0.141842_dp, 0.147162_dp, 0.152249_dp, 0.157127_dp, 0.161817_dp, &
      0.166337_dp, 0.170702_dp, 0.179015_dp, 0.194248_dp, 0.207999_dp, &
      0.220581_dp, 0.232212_dp, 0.24305_dp, 0.325071_dp, 0.382016_dp, &
      0.426452_dp, 0.463159_dp, 0.668992_dp, 0.775976_dp, 0.849298_dp, &
      0.905305_dp, 1.088998_dp, 1.253794_dp, 1.351794_dp, 1.42186_dp, &
      1.476457_dp, 1.747485_dp, 1.874323_dp, 1.958037_dp, 2.02063_dp, &
      2.219765_dp, 2.392826_dp, 2.494109_dp, 2.565985_dp, 2.621743_dp ]
  real(dp), parameter :: qnx(50) = [ &
      0.0005_dp, 0.001_dp, 0.0015_dp, 0.002_dp, 0.0025_dp, &
      0.005_dp, 0.0075_dp, 0.01_dp, 0.0125_dp, 0.015_dp, &
      0.0175_dp, 0.02_dp, 0.0225_dp, 0.025_dp, 0.0275_dp, &
      0.03_dp, 0.0325_dp, 0.035_dp, 0.0375_dp, 0.04_dp, &
      0.0425_dp, 0.045_dp, 0.0475_dp, 0.05_dp, 0.055_dp, &
      0.06_dp, 0.065_dp, 0.07_dp, 0.075_dp, 0.08_dp, &
      0.085_dp, 0.09_dp, 0.095_dp, 0.1_dp, 0.15_dp, &
      0.2_dp, 0.25_dp, 0.3_dp, 0.35_dp, 0.4_dp, &
      0.5_dp, 1.0_dp, 1.5_dp, 2.0_dp, 2.5_dp, &
      3.0_dp, 3.5_dp, 4.0_dp, 4.5_dp, 5.0_dp ]
  real(dp), parameter :: qny(50) = [ &
      0.019965_dp, 0.028394_dp, 0.034874_dp, 0.040369_dp, 0.045256_dp, &
      0.064633_dp, 0.079746_dp, 0.092708_dp, 0.104259_dp, 0.114814_dp, &
      0.124608_dp, 0.133772_dp, 0.142429_dp, 0.150739_dp, 0.158565_dp, &
      0.166229_dp, 0.173756_dp, 0.180793_dp, 0.187739_dp, 0.194489_dp, &
      0.201094_dp, 0.207572_dp, 0.213877_dp, 0.220056_dp, 0.231797_dp, &
      0.243374_dp, 0.254585_dp, 0.265472_dp, 0.27607_dp, 0.286406_dp, &
      0.296507_dp, 0.306393_dp, 0.316066_dp, 0.325586_dp, 0.413136_dp, &
      0.491599_dp, 0.564333_dp, 0.633007_dp, 0.698849_dp, 0.762455_dp, &
      0.884593_dp, 1.44552_dp, 1.97074_dp, 2.48396_dp, 2.99094_dp, &
      3.49252_dp, 3.99519_dp, 4.49238_dp, 4.99043_dp, 5.49882_dp ]
contains
  real(dp) function maxdd_expectation(mean,sd,horizon) result(v)
    real(dp),intent(in)::mean,sd,horizon
    real(dp)::x,q,gamma0
    gamma0=sqrt(pi/8.0_dp)
    if(sd<=0.0_dp.or.horizon<0.0_dp)then;v=0.0_dp;return;end if
    if(abs(mean)<=epsilon(1.0_dp))then
      v=2.0_dp*gamma0*sd*sqrt(horizon)
    else
      x=mean*mean*horizon/(2.0_dp*sd*sd)
      if(mean>0.0_dp)then
        if(x<0.0005_dp)then;q=gamma0*sqrt(2.0_dp*x)
        else if(x<=5000.0_dp)then;q=interp(log(qpx),qpy,log(x))
        else;q=0.25_dp*log(x)+0.49088_dp;end if
        v=2.0_dp*sd*sd/mean*q
      else
        if(x<0.0005_dp)then;q=gamma0*sqrt(2.0_dp*x)
        else if(x<=5.0_dp)then;q=interp(qnx,qny,x)
        else;q=x+0.5_dp;end if
        v=-2.0_dp*sd*sd/mean*q
      end if
    end if
  end function
  real(dp) function dmaxdd(x,sd,horizon,nterms) result(v)
    real(dp),intent(in)::x,sd,horizon;integer,intent(in),optional::nterms
    integer::n,nt;real(dp)::en,pn
    nt=1000;if(present(nterms))nt=nterms
    if(x<=0.0_dp)then;v=0.0_dp;return;end if
    v=0.0_dp
    do n=1,nt
      pn=2.0_dp*sin((real(n,dp)-0.5_dp)*pi)*sd*sd*(real(n,dp)-0.5_dp)*pi*horizon
      en=sd*sd*(real(n,dp)-0.5_dp)**2*pi*pi*horizon/2.0_dp
      v=v+pn*exp(-en/(x*x))/(x**3)
    end do
    v=max(v,0.0_dp)
  end function
  real(dp) function pmaxdd(q,sd,horizon,nterms) result(v)
    real(dp),intent(in)::q,sd,horizon;integer,intent(in),optional::nterms
    integer::n,nt;real(dp)::en,pn
    nt=1000;if(present(nterms))nt=nterms
    if(q<=0.0_dp)then;v=0.0_dp;return;end if
    v=0.0_dp
    do n=1,nt
      pn=2.0_dp*sin((real(n,dp)-0.5_dp)*pi)/((real(n,dp)-0.5_dp)*pi)
      en=sd*sd*(real(n,dp)-0.5_dp)**2*pi*pi*horizon/2.0_dp
      v=v+pn*(1.0_dp-exp(-en/(q*q)))
    end do
    v=min(max(v,0.0_dp),1.0_dp)
  end function
  real(dp) function rmaxdd(mean,sd,horizon) result(v)
    real(dp),intent(in)::mean,sd;integer,intent(in)::horizon
    real(dp)::path,peak,dd;integer::i
    path=0.0_dp;peak=0.0_dp;v=0.0_dp
    do i=1,horizon
      path=path+mean+sd*rnorm_lcg();peak=max(peak,path);dd=peak-path;v=max(v,dd)
    end do
  end function
  real(dp) function realized_max_drawdown(returns) result(v)
    real(dp),intent(in)::returns(:);real(dp)::wealth,peak;integer::i
    wealth=1.0_dp;peak=1.0_dp;v=0.0_dp
    do i=1,size(returns);wealth=wealth*(1.0_dp+returns(i));peak=max(peak,wealth);v=max(v,(peak-wealth)/peak);end do
  end function
  real(dp) function interp(x,y,x0) result(v)
    real(dp),intent(in)::x(:),y(:),x0;integer::i,n
    n=size(x);if(x0<=x(1))then;v=y(1);return;else if(x0>=x(n))then;v=y(n);return;end if
    do i=1,n-1;if(x0<=x(i+1))then;v=y(i)+(y(i+1)-y(i))*(x0-x(i))/(x(i+1)-x(i));return;end if;end do;v=y(n)
  end function
end module fbasics_drawdown
