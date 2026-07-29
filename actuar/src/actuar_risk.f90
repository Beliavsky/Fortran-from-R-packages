! SPDX-License-Identifier: GPL-2.0-or-later
module actuar_risk
  use actuar_kinds, only : dp
  use actuar_special, only : nan_dp
  use actuar_phase_type, only : pphtype
  implicit none
  private

  abstract interface
    function mgf_callback(x) result(y)
      import dp
      real(dp), intent(in) :: x
      real(dp) :: y
    end function mgf_callback
  end interface

  public :: adjustment_coefficient, adjustment_coefficient_poisson
  public :: ruin_exponential, ruin_cramer_lundberg_ph
contains

  function adjustment_coefficient(mgf_claim,mgf_wait,premium_rate,upper_bound) result(r)
    procedure(mgf_callback) :: mgf_claim,mgf_wait
    real(dp), intent(in) :: premium_rate,upper_bound
    real(dp) :: r,lo,hi,mid,flo,fhi,fmid
    integer :: iter
    if(premium_rate<=0.0_dp .or. upper_bound<=0.0_dp) then
      r=nan_dp(); return
    end if
    lo=max(1.0e-12_dp,upper_bound*1.0e-10_dp); hi=upper_bound
    flo=mgf_claim(lo)*mgf_wait(-lo*premium_rate)-1.0_dp
    fhi=mgf_claim(hi)*mgf_wait(-hi*premium_rate)-1.0_dp
    if(flo>0.0_dp) then
      r=0.0_dp; return
    end if
    if(fhi<=0.0_dp) then
      r=nan_dp(); return
    end if
    do iter=1,200
      mid=0.5_dp*(lo+hi)
      fmid=mgf_claim(mid)*mgf_wait(-mid*premium_rate)-1.0_dp
      if(fmid>0.0_dp) then
        hi=mid; fhi=fmid
      else
        lo=mid; flo=fmid
      end if
      if(abs(hi-lo)<=1.0e-12_dp*max(1.0_dp,mid)) exit
    end do
    r=0.5_dp*(lo+hi)
  end function adjustment_coefficient

  function adjustment_coefficient_poisson(mgf_claim,claim_rate,premium_rate,upper_bound) result(r)
    procedure(mgf_callback) :: mgf_claim
    real(dp), intent(in) :: claim_rate,premium_rate,upper_bound
    real(dp) :: r,lo,hi,mid,flo,fhi,fmid
    integer :: iter
    if(claim_rate<=0.0_dp .or. premium_rate<=0.0_dp) then
      r=nan_dp(); return
    end if
    lo=max(1.0e-12_dp,upper_bound*1.0e-10_dp); hi=upper_bound
    flo=claim_rate*(mgf_claim(lo)-1.0_dp)-premium_rate*lo
    fhi=claim_rate*(mgf_claim(hi)-1.0_dp)-premium_rate*hi
    if(flo>0.0_dp .or. fhi<=0.0_dp) then
      r=nan_dp(); return
    end if
    do iter=1,200
      mid=0.5_dp*(lo+hi)
      fmid=claim_rate*(mgf_claim(mid)-1.0_dp)-premium_rate*mid
      if(fmid>0.0_dp) then
        hi=mid; fhi=fmid
      else
        lo=mid; flo=fmid
      end if
      if(abs(hi-lo)<=1.0e-12_dp*max(1.0_dp,mid)) exit
    end do
    r=0.5_dp*(lo+hi)
  end function adjustment_coefficient_poisson

  pure function ruin_exponential(initial_surplus,claim_rate,claim_size_rate,premium_rate) result(psi)
    real(dp), intent(in) :: initial_surplus,claim_rate,claim_size_rate,premium_rate
    real(dp) :: psi,lambda_scaled
    if(initial_surplus<0.0_dp .or. claim_rate<=0.0_dp .or. &
       claim_size_rate<=0.0_dp .or. premium_rate<=0.0_dp) then
      psi=nan_dp(); return
    end if
    lambda_scaled=claim_rate/premium_rate
    if(lambda_scaled>=claim_size_rate) then
      psi=1.0_dp
    else
      psi=lambda_scaled/claim_size_rate* &
        exp(-(claim_size_rate-lambda_scaled)*initial_surplus)
    end if
  end function ruin_exponential

  function ruin_cramer_lundberg_ph(initial_surplus,claim_rate,premium_rate, &
    initial_prob,subgenerator) result(psi)
    real(dp), intent(in) :: initial_surplus,claim_rate,premium_rate
    real(dp), intent(in) :: initial_prob(:),subgenerator(:,:)
    real(dp) :: psi
    real(dp), allocatable :: invt(:,:),rhs(:,:),aug(:,:),row(:),pi_ruin(:),q(:,:)
    real(dp) :: pivot,factor
    integer :: n,i,k,p
    if(claim_rate<=0.0_dp .or. premium_rate<=0.0_dp .or. initial_surplus<0.0_dp) then
      psi=nan_dp(); return
    end if
    n=size(initial_prob)
    if(size(subgenerator,1)/=n .or. size(subgenerator,2)/=n) then
      psi=nan_dp(); return
    end if
    allocate(rhs(n,n)); rhs=0.0_dp
    do i=1,n; rhs(i,i)=1.0_dp; end do
    allocate(aug(n,2*n),row(2*n),invt(n,n))
    aug(:,:n)=subgenerator; aug(:,n+1:)=rhs
    do k=1,n
      p=k
      do i=k+1,n
        if(abs(aug(i,k))>abs(aug(p,k))) p=i
      end do
      if(abs(aug(p,k))<1.0e-14_dp) then
        psi=nan_dp(); return
      end if
      if(p/=k) then; row=aug(k,:); aug(k,:)=aug(p,:); aug(p,:)=row; end if
      pivot=aug(k,k); aug(k,:)=aug(k,:)/pivot
      do i=1,n
        if(i==k) cycle
        factor=aug(i,k); aug(i,:)=aug(i,:)-factor*aug(k,:)
      end do
    end do
    invt=aug(:,n+1:)
    allocate(pi_ruin(n),q(n,n))
    pi_ruin=-claim_rate/premium_rate*matmul(initial_prob,invt)
    q=subgenerator
    do i=1,n
      q(:,i)=q(:,i)-sum(subgenerator,dim=2)*pi_ruin(i)
    end do
    psi=1.0_dp-pphtype(initial_surplus,pi_ruin,q)
    psi=max(0.0_dp,min(1.0_dp,psi))
  end function ruin_cramer_lundberg_ph

end module actuar_risk
