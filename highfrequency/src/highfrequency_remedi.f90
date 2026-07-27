! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from highfrequency 1.0.2 by Kris Boudt, Jonathan Cornelissen,
! Scott Payseur, Onno Kleen, Emil Sjoerup, and contributors.
module highfrequency_remedi
  use highfrequency_kinds, only: dp
  implicit none
  private
  public :: remedi, choose_remedi_kn

contains

  function remedi(prices,kn,lags,make_correlation) result(values)
    real(dp),intent(in)::prices(:)
    integer,intent(in)::kn,lags(:)
    logical,intent(in),optional::make_correlation
    real(dp),allocatable::values(:)
    real(dp),allocatable::raw(:)
    logical::cor
    integer::n,i,j,lag,start_idx,end_idx
    n=size(prices)
    cor=.false.;if(present(make_correlation))cor=make_correlation
    allocate(raw(size(lags)+merge(1,0,cor)))
    raw=0.0_dp
    do j=1,size(raw)
      if(cor .and. j==1)then
        lag=0
      else
        lag=lags(j-merge(1,0,cor))
      end if
      start_idx=2*kn+1
      end_idx=n-kn-lag
      if(start_idx<=end_idx .and. kn>=1 .and. lag>=0)then
        do i=start_idx,end_idx
          raw(j)=raw(j)+(prices(i+lag)-prices(i+lag+kn))*(prices(i)-prices(i-2*kn))
        end do
        raw(j)=raw(j)/real(n,dp)
      end if
    end do
    if(cor)then
      allocate(values(size(lags)))
      if(abs(raw(1))>tiny(1.0_dp))then
        values=raw(2:)/raw(1)
      else
        values=0.0_dp
      end if
    else
      allocate(values(size(lags)))
      values=raw
    end if
  end function remedi

  integer function choose_remedi_kn(prices,kn_max,tolerance,flat_size) result(kn)
    real(dp),intent(in)::prices(:)
    integer,intent(in),optional::kn_max,flat_size
    real(dp),intent(in),optional::tolerance
    integer::kmax,m,i,k,lags(4)
    real(dp)::tol,base,threshold
    real(dp),allocatable::err(:),v(:),base_vec(:)
    kmax=10;if(present(kn_max))kmax=kn_max
    m=3;if(present(flat_size))m=flat_size
    tol=0.05_dp;if(present(tolerance))tol=tolerance
    lags=[0,1,2,3]
    allocate(err(kmax+m+2))
    base_vec=remedi(prices,1,[0])
    base=base_vec(1)
    do k=1,size(err)
      v=remedi(prices,k,lags)
      err(k)=(v(1)-v(2)-v(3)+v(4)-base)**2
    end do
    threshold=tol*maxval(err(1:max(1,kmax/2)))
    kn=1
    do i=1,kmax
      if(i+m-1<=size(err))then
        if(maxval(err(i:i+m-1))<=threshold)then
          kn=i
          return
        end if
      end if
    end do
    kn=minloc(err(1:kmax),dim=1)
  end function choose_remedi_kn

end module highfrequency_remedi
