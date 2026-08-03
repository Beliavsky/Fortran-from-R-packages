! SPDX-License-Identifier: MIT
! PortfolioTesteR modern Fortran translation
module ptr_data
  use ptr_kinds, only : dp, i8
  use ptr_utils, only : nan_dp, is_finite
  implicit none
  private
  public :: generate_sample_prices, forward_fill, align_to_indices
  public :: invert_signal, standardize_panel

contains

  subroutine generate_sample_prices(n_periods,n_assets,prices,seed)
    integer,intent(in)::n_periods,n_assets
    real(dp),allocatable,intent(out)::prices(:,:)
    integer(i8),intent(in),optional::seed
    integer(i8)::state
    real(dp)::u1,u2,common,idio,drift
    integer::t,j
    state=104729_i8;if(present(seed))state=max(1_i8,seed)
    allocate(prices(n_periods,n_assets));prices(1,:)=100.0_dp
    do t=2,n_periods
      common=0.0075_dp*normal_rng()
      do j=1,n_assets
        idio=0.012_dp*normal_rng();drift=0.0008_dp+0.00015_dp*real(j-1,dp)
        prices(t,j)=prices(t-1,j)*exp(drift+0.55_dp*common+idio)
      end do
    end do
  contains
    real(dp) function uniform_rng()
      state=mod(48271_i8*state,2147483647_i8)
      uniform_rng=max(1.0e-15_dp,real(state,dp)/2147483647.0_dp)
    end function uniform_rng
    real(dp) function normal_rng()
      u1=uniform_rng();u2=uniform_rng()
      normal_rng=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*acos(-1.0_dp)*u2)
    end function normal_rng
  end subroutine generate_sample_prices

  subroutine forward_fill(panel,out)
    real(dp),intent(in)::panel(:,:)
    real(dp),allocatable,intent(out)::out(:,:)
    integer::t,j
    allocate(out(size(panel,1),size(panel,2)));out=panel
    do j=1,size(panel,2)
      do t=2,size(panel,1)
        if(.not.is_finite(out(t,j)).and.is_finite(out(t-1,j)))out(t,j)=out(t-1,j)
      end do
    end do
  end subroutine forward_fill

  subroutine align_to_indices(source_index,source,target_index,out,forward_fill_missing)
    integer,intent(in)::source_index(:),target_index(:)
    real(dp),intent(in)::source(:,:)
    real(dp),allocatable,intent(out)::out(:,:)
    logical,intent(in),optional::forward_fill_missing
    logical::ffill
    integer::i,j,k,best
    ffill=.true.;if(present(forward_fill_missing))ffill=forward_fill_missing
    allocate(out(size(target_index),size(source,2)));out=nan_dp()
    do i=1,size(target_index)
      best=0
      do k=1,size(source_index)
        if(source_index(k)==target_index(i))then;best=k;exit;end if
        if(ffill.and.source_index(k)<=target_index(i))best=k
        if(source_index(k)>target_index(i).and.ffill)exit
      end do
      if(best>0)then
        do j=1,size(source,2);out(i,j)=source(best,j);end do
      end if
    end do
  end subroutine align_to_indices

  subroutine invert_signal(signal,out)
    real(dp),intent(in)::signal(:,:)
    real(dp),allocatable,intent(out)::out(:,:)
    allocate(out(size(signal,1),size(signal,2)));out=nan_dp()
    where(is_finite(signal))out=-signal
  end subroutine invert_signal

  subroutine standardize_panel(panel,out,by_row)
    real(dp),intent(in)::panel(:,:)
    real(dp),allocatable,intent(out)::out(:,:)
    logical,intent(in),optional::by_row
    logical::rowwise
    real(dp)::mu,sd
    integer::i,j,n
    rowwise=.true.;if(present(by_row))rowwise=by_row
    allocate(out(size(panel,1),size(panel,2)));out=nan_dp()
    if(rowwise)then
      do i=1,size(panel,1)
        mu=0.0_dp;n=0
        do j=1,size(panel,2);if(is_finite(panel(i,j)))then;mu=mu+panel(i,j);n=n+1;end if;end do
        if(n<2)cycle;mu=mu/real(n,dp);sd=0.0_dp
        do j=1,size(panel,2);if(is_finite(panel(i,j)))sd=sd+(panel(i,j)-mu)**2;end do
        sd=sqrt(sd/real(n-1,dp));if(sd<=0.0_dp)cycle
        do j=1,size(panel,2);if(is_finite(panel(i,j)))out(i,j)=(panel(i,j)-mu)/sd;end do
      end do
    else
      do j=1,size(panel,2)
        mu=0.0_dp;n=0
        do i=1,size(panel,1);if(is_finite(panel(i,j)))then;mu=mu+panel(i,j);n=n+1;end if;end do
        if(n<2)cycle;mu=mu/real(n,dp);sd=0.0_dp
        do i=1,size(panel,1);if(is_finite(panel(i,j)))sd=sd+(panel(i,j)-mu)**2;end do
        sd=sqrt(sd/real(n-1,dp));if(sd<=0.0_dp)cycle
        do i=1,size(panel,1);if(is_finite(panel(i,j)))out(i,j)=(panel(i,j)-mu)/sd;end do
      end do
    end if
  end subroutine standardize_panel

end module ptr_data
