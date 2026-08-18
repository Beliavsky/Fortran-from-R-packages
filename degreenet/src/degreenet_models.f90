! SPDX-License-Identifier: GPL-3.0-or-later
! Based on 'statnet' project software (statnet.org).
module degreenet_models
  use degreenet_kinds, only : dp, huge_neg
  use degreenet_math, only : geom_pmf, poisson_pmf
  use degreenet_distributions, only : dyule, ddp, dwar, ddqe, dghdi, &
    dnb_degree, dpe_pmf, dcmp_natural, dpln
  use degreenet_compound, only : dgyule, dgeodp, dnbyule, dnbwar, dgwar
  implicit none
  private
  integer, parameter, public :: MODEL_YULE=1, MODEL_DP=2, MODEL_WARING=3
  integer, parameter, public :: MODEL_DQE=4, MODEL_GHDI=5, MODEL_NB=6
  integer, parameter, public :: MODEL_PE=7, MODEL_CMP=8, MODEL_PLN=9
  integer, parameter, public :: MODEL_GYULE=10, MODEL_GEODP=11
  integer, parameter, public :: MODEL_NBYULE=12, MODEL_NBWAR=13
  integer, parameter, public :: MODEL_GWAR=14, MODEL_GEOM=15, MODEL_POIS=16
  public :: model_pmf, loglik_model, model_cdf, grouped_probability
  public :: complete_fit_stats

contains
  real(dp) function model_pmf(model,par,x,cutoff) result(p)
    integer,intent(in)::model,x,cutoff
    real(dp),intent(in)::par(:)
    select case(model)
    case(MODEL_YULE); p=dyule(par(1),x,cutoff)
    case(MODEL_DP); p=ddp(par(1),x,cutoff)
    case(MODEL_WARING); p=dwar(par(1:2),x,cutoff)
    case(MODEL_DQE); p=ddqe(par(1:2),x,cutoff)
    case(MODEL_GHDI); p=dghdi(par(1:4),x,cutoff)
    case(MODEL_NB); p=dnb_degree(par(1:2),x,cutoff,.true.)
    case(MODEL_PE); p=dpe_pmf(par(1:2),x,cutoff)
    case(MODEL_CMP); p=dcmp_natural(par(1),par(2),x,cutoff,1e-10_dp)
    case(MODEL_PLN); p=dpln(par(1:2),x,cutoff,.true.,32)
    case(MODEL_GYULE); p=dgyule(par(1:2),x,cutoff)
    case(MODEL_GEODP); p=dgeodp(par(1:2),x,cutoff)
    case(MODEL_NBYULE); p=dnbyule(par(1:3),x,cutoff)
    case(MODEL_NBWAR); p=dnbwar(par(1:4),x,cutoff)
    case(MODEL_GWAR); p=dgwar(par(1:3),x,cutoff)
    case(MODEL_GEOM)
      if(x<cutoff.or.par(1)<=1.0_dp)then;p=0.0_dp
      else;p=geom_pmf(x-cutoff,1.0_dp/par(1));end if
    case(MODEL_POIS)
      if(x<cutoff.or.par(1)<0.0_dp)then;p=0.0_dp
      else;p=poisson_pmf(x-cutoff,par(1));end if
    case default; p=0.0_dp
    end select
  end function model_pmf

  real(dp) function loglik_model(model,par,x,cutoff,cutabove) result(ll)
    integer,intent(in)::model,cutoff,cutabove
    real(dp),intent(in)::par(:)
    integer,intent(in)::x(:)
    integer::i,k,n
    real(dp)::p,z
    ll=0.0_dp;n=0
    z=0.0_dp
    if(cutabove<1000)then
      do k=cutoff,cutabove;z=z+model_pmf(model,par,k,cutoff);end do
      if(z<=0.0_dp)then;ll=huge_neg;return;end if
    else;z=1.0_dp;end if
    do i=1,size(x)
      if(x(i)<cutoff.or.x(i)>cutabove)cycle
      p=model_pmf(model,par,x(i),cutoff)
      if(p<=0.0_dp.or..not.(p<huge(1.0_dp)))then;ll=huge_neg;return;end if
      ll=ll+log(p);n=n+1
    end do
    if(n==0)ll=huge_neg
    if(z>0.0_dp)ll=ll-real(n,dp)*log(z)
  end function loglik_model

  real(dp) function model_cdf(model,par,q,cutoff) result(cdf)
    integer,intent(in)::model,q,cutoff
    real(dp),intent(in)::par(:)
    integer::k
    cdf=0.0_dp
    if(q<cutoff)return
    do k=cutoff,q;cdf=cdf+model_pmf(model,par,k,cutoff);end do
    cdf=min(1.0_dp,max(0.0_dp,cdf))
  end function model_cdf

  real(dp) function grouped_probability(model,par,code,cutoff) result(p)
    integer,intent(in)::model,code,cutoff
    real(dp),intent(in)::par(:)
    integer::lo,hi,k
    select case(code)
    case(:4);lo=code;hi=code
    case(5);lo=5;hi=10
    case(6);lo=11;hi=20
    case(7);lo=21;hi=100
    case(8)
      p=max(0.0_dp,1.0_dp-model_cdf(model,par,100,cutoff));return
    case(9);lo=5;hi=20
    case(10);lo=5;hi=100
    case default;p=0.0_dp;return
    end select
    lo=max(lo,cutoff);p=0.0_dp
    do k=lo,hi;p=p+model_pmf(model,par,k,cutoff);end do
  end function grouped_probability

  subroutine complete_fit_stats(loglik,n,np,aicc,bic)
    real(dp),intent(in)::loglik
    integer,intent(in)::n,np
    real(dp),intent(out)::aicc,bic
    if(n>np+1)then
      aicc=-2.0_dp*loglik+2.0_dp*np+2.0_dp*np*(np+1.0_dp)/real(n-np-1,dp)
    else;aicc=huge(1.0_dp);end if
    bic=-2.0_dp*loglik+real(np,dp)*log(real(max(n,1),dp))
  end subroutine complete_fit_stats
end module degreenet_models
