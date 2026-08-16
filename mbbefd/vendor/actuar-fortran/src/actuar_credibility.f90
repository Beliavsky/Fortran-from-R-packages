! SPDX-License-Identifier: GPL-2.0-or-later
module actuar_credibility
  use actuar_kinds, only : dp
  use actuar_special, only : nan_dp
  use actuar_types, only : credibility_result
  implicit none
  private
  public :: buhlmann_straub, bayes_poisson_gamma, bayes_bernoulli_beta
  public :: bayes_normal_normal, credibility_premium
contains

  function buhlmann_straub(ratios,weights,iterative,tolerance,max_iterations) result(res)
    real(dp), intent(in) :: ratios(:,:),weights(:,:)
    logical, intent(in), optional :: iterative
    real(dp), intent(in), optional :: tolerance
    integer, intent(in), optional :: max_iterations
    type(credibility_result) :: res
    real(dp), allocatable :: ws(:),rw(:),cred(:)
    real(dp) :: totalw,s2,a,xw,den,tol,old,zmean
    integer :: i,j,ncontracts,ntotal,iter,maxit
    logical :: do_iter
    if(any(shape(ratios)/=shape(weights)) .or. size(ratios,1)<2 .or. size(ratios,2)<2) then
      res%message='invalid ratio or weight dimensions'; return
    end if
    allocate(ws(size(ratios,1)),rw(size(ratios,1)),cred(size(ratios,1)))
    ws=sum(weights,dim=2); rw=0.0_dp
    do i=1,size(ratios,1)
      if(ws(i)>0.0_dp) rw(i)=sum(weights(i,:)*ratios(i,:))/ws(i)
    end do
    ncontracts=count(ws>0.0_dp); ntotal=count(weights>0.0_dp)
    if(ncontracts<2 .or. ntotal<=ncontracts) then
      res%message='insufficient nonzero experience'; return
    end if
    s2=0.0_dp
    do i=1,size(ratios,1)
      do j=1,size(ratios,2)
        s2=s2+weights(i,j)*(ratios(i,j)-rw(i))**2
      end do
    end do
    s2=s2/real(ntotal-ncontracts,dp)
    totalw=sum(ws); xw=sum(ws*rw)/totalw
    den=totalw**2-sum(ws**2)
    if(abs(den)<1.0e-14_dp) then
      res%message='degenerate portfolio weights'; return
    end if
    a=totalw*(sum(ws*(rw-xw)**2)-real(ncontracts-1,dp)*s2)/den
    a=max(0.0_dp,a)
    do_iter=.false.; if(present(iterative)) do_iter=iterative
    tol=sqrt(epsilon(1.0_dp)); if(present(tolerance)) tol=tolerance
    maxit=100; if(present(max_iterations)) maxit=max_iterations
    if(do_iter .and. a>0.0_dp) then
      do iter=1,maxit
        old=a
        cred=1.0_dp/(1.0_dp+s2/(ws*a))
        where(ws<=0.0_dp) cred=0.0_dp
        zmean=sum(cred*rw)/sum(cred)
        a=sum(cred*(rw-zmean)**2)/real(ncontracts-1,dp)
        if(abs(a-old)<=tol*max(1.0_dp,abs(old))) exit
      end do
    end if
    if(a>0.0_dp) then
      cred=1.0_dp/(1.0_dp+s2/(ws*a))
      where(ws<=0.0_dp) cred=0.0_dp
      zmean=sum(cred*rw)/sum(cred)
    else
      cred=0.0_dp; zmean=xw
    end if
    allocate(res%means(size(rw)),res%weights(size(ws)),res%estimates(size(rw)))
    res%means=rw; res%weights=cred
    res%estimates=zmean+cred*(rw-zmean)
    res%collective_mean=zmean; res%process_variance=s2
    res%structural_variance=a
    if(a>0.0_dp) res%k=s2/a
    res%ok=.true.; res%message='Buhlmann-Straub credibility'
  end function buhlmann_straub

  pure function credibility_premium(collective,individual,z) result(p)
    real(dp), intent(in) :: collective,individual,z
    real(dp) :: p
    p=collective+z*(individual-collective)
  end function credibility_premium

  function bayes_poisson_gamma(observations,shape,scale) result(res)
    real(dp), intent(in) :: observations(:),shape,scale
    type(credibility_result) :: res
    real(dp) :: collective,z,individual
    integer :: n
    if(shape<=0.0_dp .or. scale<=0.0_dp) then
      res%message='invalid gamma prior'; return
    end if
    n=size(observations); collective=shape*scale
    individual=merge(sum(observations)/real(n,dp),0.0_dp,n>0)
    z=real(n,dp)/(real(n,dp)+1.0_dp/scale)
    allocate(res%means(1),res%weights(1),res%estimates(1))
    res%means=individual; res%weights=z
    res%estimates=credibility_premium(collective,individual,z)
    res%collective_mean=collective
    res%process_variance=collective
    res%structural_variance=collective*scale
    res%k=1.0_dp/scale; res%ok=.true.; res%message='Poisson-gamma Bayes credibility'
  end function bayes_poisson_gamma

  function bayes_bernoulli_beta(observations,shape1,shape2) result(res)
    real(dp), intent(in) :: observations(:),shape1,shape2
    type(credibility_result) :: res
    real(dp) :: collective,z,individual,k
    integer :: n
    if(shape1<=0.0_dp .or. shape2<=0.0_dp) then
      res%message='invalid beta prior'; return
    end if
    n=size(observations); k=shape1+shape2; collective=shape1/k
    individual=merge(sum(observations)/real(n,dp),0.0_dp,n>0)
    z=real(n,dp)/(real(n,dp)+k)
    allocate(res%means(1),res%weights(1),res%estimates(1))
    res%means=individual; res%weights=z
    res%estimates=credibility_premium(collective,individual,z)
    res%collective_mean=collective; res%k=k
    res%process_variance=shape1*shape2/(k*(k+1.0_dp))
    res%structural_variance=shape1*shape2/(k*k*(k+1.0_dp))
    res%ok=.true.; res%message='Bernoulli-beta Bayes credibility'
  end function bayes_bernoulli_beta

  function bayes_normal_normal(observations,prior_mean,prior_sd,likelihood_sd) result(res)
    real(dp), intent(in) :: observations(:),prior_mean,prior_sd,likelihood_sd
    type(credibility_result) :: res
    real(dp) :: z,individual,k
    integer :: n
    if(prior_sd<=0.0_dp .or. likelihood_sd<=0.0_dp) then
      res%message='invalid normal model standard deviation'; return
    end if
    n=size(observations); k=(likelihood_sd/prior_sd)**2
    individual=merge(sum(observations)/real(n,dp),0.0_dp,n>0)
    z=real(n,dp)/(real(n,dp)+k)
    allocate(res%means(1),res%weights(1),res%estimates(1))
    res%means=individual; res%weights=z
    res%estimates=credibility_premium(prior_mean,individual,z)
    res%collective_mean=prior_mean; res%k=k
    res%process_variance=likelihood_sd**2; res%structural_variance=prior_sd**2
    res%ok=.true.; res%message='normal-normal Bayes credibility'
  end function bayes_normal_normal

end module actuar_credibility
