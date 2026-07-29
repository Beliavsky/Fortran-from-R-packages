! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2016 Marius Hofert, Kurt Hornik and Alexander J. McNeil
module qrmtools_brownian
  use qrmtools_kinds, only : dp
  use qrmtools_types, only : brownian_result
  use qrmtools_stats, only : random_normal, normal_quantile, seed_random
  implicit none
  private
  public :: r_brownian, de_browning
contains
  function r_brownian(n_paths,t,d,drift,volatility,process_type,initial,u,seed) result(output)
    integer, intent(in) :: n_paths
    real(dp), intent(in) :: t(:)
    integer, intent(in), optional :: d
    real(dp), intent(in), optional :: drift(:),volatility(:),initial(:),u(:,:)
    character(len=*), intent(in), optional :: process_type
    integer, intent(in), optional :: seed
    type(brownian_result) :: output
    integer :: dimensions,n_steps,i,j,k,row
    real(dp), allocatable :: mu(:),sigma(:),s0(:),w(:,:),terminal(:,:)
    real(dp) :: z,dt
    character(len=8) :: kind
    dimensions=1; if(present(d))dimensions=d
    n_steps=size(t)-1
    if(n_paths<1 .or. dimensions<1 .or. n_steps<0 .or. abs(t(1))>1.0e-14_dp) then
      output%message='Invalid Brownian dimensions or time grid.'; return
    end if
    if(n_steps>0) then
      if(any(t(2:)<=t(:size(t)-1))) then; output%message='Times must be strictly increasing.'; return; end if
    end if
    allocate(mu(dimensions),sigma(dimensions),s0(dimensions))
    mu=0.0_dp; sigma=1.0_dp; s0=1.0_dp
    if(present(drift)) then
      if(size(drift)==1) then; mu=drift(1); else if(size(drift)==dimensions) then; mu=drift
      else; output%message='Drift has wrong length.'; return; end if
    end if
    if(present(volatility)) then
      if(size(volatility)==1) then; sigma=volatility(1)
      else if(size(volatility)==dimensions) then; sigma=volatility
      else; output%message='Volatility has wrong length.'; return; end if
    end if
    if(any(sigma<=0.0_dp)) then; output%message='Volatilities must be positive.'; return; end if
    if(present(initial)) then
      if(size(initial)==1) then; s0=initial(1)
      else if(size(initial)==dimensions) then; s0=initial
      else; output%message='Initial value has wrong length.'; return; end if
    end if
    kind='BM'; if(present(process_type))kind=adjustl(process_type)
    if(present(seed))call seed_random(seed)
    if(present(u)) then
      if(size(u,1)/=n_paths*n_steps .or. size(u,2)/=dimensions) then
        output%message='Uniform input has wrong dimensions.'; return
      end if
    end if
    allocate(output%paths(n_paths,n_steps+1,dimensions),w(n_paths,dimensions))
    output%paths=0.0_dp; w=0.0_dp
    do k=1,n_steps
      dt=t(k+1)-t(k)
      do i=1,n_paths
        row=(i-1)*n_steps+k
        do j=1,dimensions
          if(present(u)) then; z=normal_quantile(u(row,j)); else; z=random_normal(); end if
          w(i,j)=w(i,j)+sqrt(dt)*z
          output%paths(i,k+1,j)=mu(j)*t(k+1)+sigma(j)*w(i,j)
        end do
      end do
    end do
    select case(trim(kind))
    case('GBM','gbm')
      if(any(s0<=0.0_dp)) then; output%message='GBM initial values must be positive.'; return; end if
      do j=1,dimensions; output%paths(:,:,j)=s0(j)*exp(output%paths(:,:,j)); end do
    case('BB','bb')
      if(t(size(t))<=0.0_dp) then; output%message='Brownian bridge requires positive terminal time.'; return; end if
      allocate(terminal(n_paths,dimensions)); terminal=output%paths(:,n_steps+1,:)
      do k=1,n_steps+1
        output%paths(:,k,:)=output%paths(:,k,:)-(t(k)/t(size(t)))*terminal
      end do
    case('BM','bm')
    case default
      output%message='Unknown Brownian process type.'; return
    end select
    output%ok=.true.
  end function r_brownian

  function de_browning(x,t,drift,volatility,process_type) result(increments)
    real(dp), intent(in) :: x(:,:,:),t(:)
    real(dp), intent(in), optional :: drift(:),volatility(:)
    character(len=*), intent(in), optional :: process_type
    real(dp), allocatable :: increments(:,:,:)
    real(dp), allocatable :: mu(:),sigma(:),work(:,:,:)
    integer :: n_paths,n_steps,dimensions,i,j,k
    real(dp) :: dt
    character(len=8) :: kind
    n_paths=size(x,1); n_steps=size(x,2)-1; dimensions=size(x,3)
    allocate(increments(n_paths,n_steps,dimensions),mu(dimensions),sigma(dimensions),work(n_paths,n_steps+1,dimensions))
    mu=0.0_dp; sigma=1.0_dp
    if(present(drift)) then
      if(size(drift)==1) mu=drift(1)
      if(size(drift)==dimensions) mu=drift
    end if
    if(present(volatility)) then
      if(size(volatility)==1) sigma=volatility(1)
      if(size(volatility)==dimensions) sigma=volatility
    end if
    kind='BM'; if(present(process_type))kind=adjustl(process_type)
    work=x
    if(trim(kind)=='GBM' .or. trim(kind)=='gbm') then
      do j=1,dimensions; do i=1,n_paths
        work(i,:,j)=log(x(i,:,j)/x(i,1,j))
      end do; end do
    end if
    do k=1,n_steps
      dt=t(k+1)-t(k)
      do j=1,dimensions
        increments(:,k,j)=((work(:,k+1,j)-mu(j)*t(k+1))/sigma(j)-&
          (work(:,k,j)-mu(j)*t(k))/sigma(j))/sqrt(dt)
      end do
    end do
  end function de_browning
end module qrmtools_brownian
