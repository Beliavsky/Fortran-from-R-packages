! Common links from gamlss.dist R/make-link-gamlss.R.
! SPDX-License-Identifier: GPL-3.0-only
module gamlss_links
   use gamlss_kinds, only : dp, pi
   use gamlss_special, only : normal_cdf, normal_quantile, expm1_v
   implicit none
   private
   integer,parameter,public :: LINK_IDENTITY=1,LINK_LOG=2,LINK_LOGIT=3,LINK_PROBIT=4
   integer,parameter,public :: LINK_CAUCHIT=5,LINK_CLOGLOG=6,LINK_SQRT=7,LINK_INVERSE=8
   integer,parameter,public :: LINK_INV_MU2=9,LINK_MU2=10,LINK_LOGSHIFT0=11,LINK_LOGSHIFT1=12,LINK_LOGSHIFT2=13
   integer,parameter,public :: LINK_M11=14,LINK_02=15,LINK_05=16
   public :: linkfun,linkinv,mu_eta
contains
   elemental real(dp) function logistic(x) result(v)
      real(dp),intent(in)::x
      if(x>=0)then;v=1/(1+exp(-min(x,700.0_dp)));else;v=exp(max(x,-700.0_dp))/(1+exp(max(x,-700.0_dp)));end if
   end function
   elemental real(dp) function linkfun(mu,link) result(v)
      real(dp),intent(in)::mu;integer,intent(in)::link;real(dp)::delta,lo,hi
      select case(link)
      case(LINK_IDENTITY);v=mu
      case(LINK_LOG);v=log(mu)
      case(LINK_LOGIT);v=log(mu/(1-mu))
      case(LINK_PROBIT);v=normal_quantile(mu)
      case(LINK_CAUCHIT);v=tan(pi*(mu-.5_dp))
      case(LINK_CLOGLOG);v=log(-log(1-mu))
      case(LINK_SQRT);v=sqrt(mu)
      case(LINK_INVERSE);v=1/mu
      case(LINK_INV_MU2);v=1/(mu*mu)
      case(LINK_MU2);v=mu*mu
      case(LINK_LOGSHIFT0);v=log(mu-1e-5_dp)
      case(LINK_LOGSHIFT1);v=log(mu-1+1e-5_dp)
      case(LINK_LOGSHIFT2);v=log(mu-2+1e-5_dp)
      case(LINK_M11,LINK_02,LINK_05)
         delta=1e-10_dp
         if(link==LINK_M11)then;lo=-1-delta;hi=1+delta
         else if(link==LINK_02)then;lo=0;hi=2+delta
         else;lo=0;hi=5+delta;end if
         v=log((mu-lo)/(hi-mu))
      case default;v=mu
      end select
   end function
   elemental real(dp) function linkinv(eta,link) result(v)
      real(dp),intent(in)::eta;integer,intent(in)::link;real(dp)::e,delta,lo,hi
      select case(link)
      case(LINK_IDENTITY);v=eta
      case(LINK_LOG);v=max(exp(min(eta,700.0_dp)),tiny(1.0_dp))
      case(LINK_LOGIT);v=logistic(eta)
      case(LINK_PROBIT);v=max(tiny(1.0_dp),min(1-tiny(1.0_dp),normal_cdf(eta)))
      case(LINK_CAUCHIT);v=.5_dp+atan(eta)/pi
      case(LINK_CLOGLOG);v=max(tiny(1.0_dp),min(1-tiny(1.0_dp),-expm1_v(-exp(min(eta,700.0_dp)))))
      case(LINK_SQRT);v=eta*eta
      case(LINK_INVERSE);v=1/eta
      case(LINK_INV_MU2);v=1/sqrt(eta)
      case(LINK_MU2);v=sqrt(eta)
      case(LINK_LOGSHIFT0);v=1e-5_dp+max(tiny(1.0_dp),exp(min(eta,700.0_dp)))
      case(LINK_LOGSHIFT1);v=1+max(tiny(1.0_dp),exp(min(eta,700.0_dp)))
      case(LINK_LOGSHIFT2);v=2+max(tiny(1.0_dp),exp(min(eta,700.0_dp)))
      case(LINK_M11,LINK_02,LINK_05)
         delta=1e-10_dp
         if(link==LINK_M11)then;lo=-1-delta;hi=1+delta
         else if(link==LINK_02)then;lo=0;hi=2+delta
         else;lo=0;hi=5+delta;end if
         e=exp(max(-700.0_dp,min(700.0_dp,eta)));v=(hi*e+lo)/(1+e)
      case default;v=eta
      end select
   end function
   elemental real(dp) function mu_eta(eta,link) result(v)
      real(dp),intent(in)::eta;integer,intent(in)::link;real(dp)::m,delta,lo,hi,e
      select case(link)
      case(LINK_IDENTITY);v=1
      case(LINK_LOG);v=max(exp(min(eta,700.0_dp)),tiny(1.0_dp))
      case(LINK_LOGIT);m=logistic(eta);v=max(m*(1-m),tiny(1.0_dp))
      case(LINK_PROBIT);v=max(exp(-.5_dp*eta*eta)/sqrt(2*pi),tiny(1.0_dp))
      case(LINK_CAUCHIT);v=max(1/(pi*(1+eta*eta)),tiny(1.0_dp))
      case(LINK_CLOGLOG);e=exp(min(eta,700.0_dp));v=max(e*exp(-e),tiny(1.0_dp))
      case(LINK_SQRT);v=2*eta
      case(LINK_INVERSE);v=-1/(eta*eta)
      case(LINK_INV_MU2);v=-1/(2*eta**1.5_dp)
      case(LINK_MU2);v=.5_dp/sqrt(eta)
      case(LINK_LOGSHIFT0,LINK_LOGSHIFT1,LINK_LOGSHIFT2);v=max(exp(min(eta,700.0_dp)),tiny(1.0_dp))
      case(LINK_M11,LINK_02,LINK_05)
         delta=1e-10_dp
         if(link==LINK_M11)then;lo=-1-delta;hi=1+delta
         else if(link==LINK_02)then;lo=0;hi=2+delta
         else;lo=0;hi=5+delta;end if
         m=logistic(eta);v=(hi-lo)*m*(1-m)
      case default;v=1
      end select
   end function
end module gamlss_links
