! SPDX-License-Identifier: GPL-3.0-only
! Native Fortran translation of the computational core of rvinecopulib.
module rvine_cvine
  use rvine_kinds, only : dp, eps_prob, clamp_prob
  use rvine_math, only : random_uniform
  use rvine_bicop, only : bicop_model, make_bicop, bicop_indep
  use rvine_fit, only : select_bicop, fit_bicop
  implicit none
  private
  public :: cvine_model, make_cvine

  type :: cvine_model
    integer :: dimension = 0
    integer, allocatable :: order(:)
    type(bicop_model), allocatable :: pair(:,:)
    integer :: nobs = 0
    real(dp) :: loglik = 0.0_dp
  contains
    procedure :: pdf => cvine_pdf
    procedure :: log_density => cvine_log_density
    procedure :: fit => cvine_fit
    procedure :: fit_fixed => cvine_fit_fixed
    procedure :: rosenblatt => cvine_rosenblatt
    procedure :: inverse_rosenblatt => cvine_inverse_rosenblatt
    procedure :: simulate => cvine_simulate
    procedure :: cdf => cvine_cdf
    procedure :: npars => cvine_npars
    procedure :: aic => cvine_aic
    procedure :: bic => cvine_bic
    procedure :: truncate => cvine_truncate
  end type cvine_model

contains

  function make_cvine(dimension, order) result(model)
    integer, intent(in) :: dimension
    integer, intent(in), optional :: order(:)
    type(cvine_model) :: model
    integer :: i,j
    model%dimension=dimension
    allocate(model%order(dimension),model%pair(dimension,dimension))
    if (present(order)) then
      model%order=order
    else
      model%order=[(i,i=1,dimension)]
    end if
    do i=1,dimension
      do j=1,dimension
        model%pair(i,j)=make_bicop(bicop_indep)
      end do
    end do
  end function make_cvine

  integer function cvine_npars(self) result(n)
    class(cvine_model), intent(in) :: self
    integer :: i,j
    n=0
    do i=1,self%dimension-1
      do j=i+1,self%dimension
        n=n+self%pair(i,j)%npar
      end do
    end do
  end function cvine_npars

  real(dp) function cvine_aic(self) result(value)
    class(cvine_model), intent(in) :: self
    value=-2.0_dp*self%loglik+2.0_dp*real(self%npars(),dp)
  end function cvine_aic

  real(dp) function cvine_bic(self) result(value)
    class(cvine_model), intent(in) :: self
    if (self%nobs>0) then
      value=-2.0_dp*self%loglik+log(real(self%nobs,dp))*real(self%npars(),dp)
    else
      value=huge(1.0_dp)
    end if
  end function cvine_bic

  subroutine cvine_truncate(self,level)
    class(cvine_model), intent(inout) :: self
    integer, intent(in) :: level
    integer :: i,j
    do i=1,self%dimension-1
      if (i>level) then
        do j=i+1,self%dimension
          self%pair(i,j)=make_bicop(bicop_indep)
        end do
      end if
    end do
  end subroutine cvine_truncate

  real(dp) function cvine_pdf(self,u) result(value)
    class(cvine_model), intent(in) :: self
    real(dp), intent(in) :: u(:)
    real(dp) :: ld
    ld=self%log_density(u)
    value=exp(min(log(huge(1.0_dp)),ld))
  end function cvine_pdf

  real(dp) function cvine_log_density(self,u) result(value)
    class(cvine_model), intent(in) :: self
    real(dp), intent(in) :: u(:)
    integer :: d,t,j
    real(dp), allocatable :: q(:,:)
    real(dp) :: x,y,den
    d=self%dimension
    allocate(q(0:d-1,d))
    q=0.0_dp
    do j=1,d
      q(0,j)=clamp_prob(u(self%order(j)))
    end do
    value=0.0_dp
    do t=1,d-1
      x=q(t-1,t)
      do j=t+1,d
        y=q(t-1,j)
        den=self%pair(t,j)%pdf(x,y)
        value=value+log(max(tiny(1.0_dp),den))
        q(t,j)=self%pair(t,j)%hfunc1(x,y)
      end do
    end do
  end function cvine_log_density

  subroutine cvine_fit(self,data,families,criterion,allow_rotations)
    class(cvine_model), intent(inout) :: self
    real(dp), intent(in) :: data(:,:)
    integer, intent(in), optional :: families(:)
    character(len=*), intent(in), optional :: criterion
    logical, intent(in), optional :: allow_rotations
    integer :: d,n,t,j,k
    real(dp), allocatable :: q(:,:,:),pairdata(:,:)
    real(dp) :: x,y
    d=self%dimension; n=size(data,2)
    allocate(q(0:d-1,d,n),pairdata(2,n)); q=0.0_dp
    do j=1,d
      q(0,j,:)=data(self%order(j),:)
    end do
    do t=1,d-1
      pairdata(1,:)=q(t-1,t,:)
      do j=t+1,d
        pairdata(2,:)=q(t-1,j,:)
        call select_dispatch(pairdata,self%pair(t,j),families,criterion,allow_rotations)
        do k=1,n
          x=pairdata(1,k); y=pairdata(2,k)
          q(t,j,k)=self%pair(t,j)%hfunc1(x,y)
        end do
      end do
    end do
    self%nobs=n; self%loglik=0.0_dp
    do k=1,n
      self%loglik=self%loglik+self%log_density(data(:,k))
    end do
  end subroutine cvine_fit

  subroutine select_dispatch(data,model,families,criterion,allow_rotations)
    real(dp), intent(in) :: data(:,:)
    type(bicop_model), intent(out) :: model
    integer, intent(in), optional :: families(:)
    character(len=*), intent(in), optional :: criterion
    logical, intent(in), optional :: allow_rotations
    if (present(families)) then
      if (present(criterion)) then
        if (present(allow_rotations)) then
          call select_bicop(data,model,families,criterion,allow_rotations)
        else
          call select_bicop(data,model,families,criterion)
        end if
      else if (present(allow_rotations)) then
        call select_bicop(data,model,families,allow_rotations=allow_rotations)
      else
        call select_bicop(data,model,families)
      end if
    else if (present(criterion)) then
      if (present(allow_rotations)) then
        call select_bicop(data,model,criterion=criterion,allow_rotations=allow_rotations)
      else
        call select_bicop(data,model,criterion=criterion)
      end if
    else if (present(allow_rotations)) then
      call select_bicop(data,model,allow_rotations=allow_rotations)
    else
      call select_bicop(data,model)
    end if
  end subroutine select_dispatch

  subroutine cvine_fit_fixed(self,data)
    class(cvine_model), intent(inout) :: self
    real(dp), intent(in) :: data(:,:)
    integer :: d,n,t,j,k,fam,rot
    real(dp), allocatable :: q(:,:,:),pairdata(:,:)
    type(bicop_model) :: fitted
    d=self%dimension; n=size(data,2)
    allocate(q(0:d-1,d,n),pairdata(2,n)); q=0.0_dp
    do j=1,d
      q(0,j,:)=data(self%order(j),:)
    end do
    do t=1,d-1
      pairdata(1,:)=q(t-1,t,:)
      do j=t+1,d
        pairdata(2,:)=q(t-1,j,:)
        fam=self%pair(t,j)%family; rot=self%pair(t,j)%rotation
        call fit_bicop(pairdata,fam,rot,fitted)
        self%pair(t,j)=fitted
        do k=1,n
          q(t,j,k)=fitted%hfunc1(pairdata(1,k),pairdata(2,k))
        end do
      end do
    end do
    self%nobs=n; self%loglik=0.0_dp
    do k=1,n
      self%loglik=self%loglik+self%log_density(data(:,k))
    end do
  end subroutine cvine_fit_fixed

  subroutine cvine_rosenblatt(self,data,z)
    class(cvine_model), intent(in) :: self
    real(dp), intent(in) :: data(:,:)
    real(dp), intent(out) :: z(:,:)
    integer :: d,n,t,j,k
    real(dp), allocatable :: q(:,:,:)
    d=self%dimension; n=size(data,2)
    allocate(q(0:d-1,d,n)); q=0.0_dp; z=0.0_dp
    do j=1,d
      q(0,j,:)=data(self%order(j),:)
    end do
    z(self%order(1),:)=q(0,1,:)
    do t=1,d-1
      do j=t+1,d
        do k=1,n
          q(t,j,k)=self%pair(t,j)%hfunc1(q(t-1,t,k),q(t-1,j,k))
        end do
      end do
      z(self%order(t+1),:)=q(t,t+1,:)
    end do
    z=min(1.0_dp-eps_prob,max(eps_prob,z))
  end subroutine cvine_rosenblatt

  subroutine cvine_inverse_rosenblatt(self,z,data)
    class(cvine_model), intent(in) :: self
    real(dp), intent(in) :: z(:,:)
    real(dp), intent(out) :: data(:,:)
    integer :: d,n,k,j,iter
    real(dp), allocatable :: uord(:),word(:)
    real(dp) :: lo,hi,mid,target,cval
    d=self%dimension; n=size(z,2)
    allocate(uord(d),word(d)); data=0.0_dp
    do k=1,n
      do j=1,d
        word(j)=clamp_prob(z(self%order(j),k))
      end do
      uord(1)=word(1)
      do j=2,d
        target=word(j); lo=eps_prob; hi=1.0_dp-eps_prob
        do iter=1,65
          mid=0.5_dp*(lo+hi); uord(j)=mid
          cval=conditional_last(self,uord,j)
          if (cval<target) then
            lo=mid
          else
            hi=mid
          end if
        end do
        uord(j)=0.5_dp*(lo+hi)
      end do
      do j=1,d
        data(self%order(j),k)=uord(j)
      end do
    end do
  end subroutine cvine_inverse_rosenblatt

  real(dp) function conditional_last(self,uord,m) result(value)
    class(cvine_model), intent(in) :: self
    real(dp), intent(in) :: uord(:)
    integer, intent(in) :: m
    integer :: t,j
    real(dp), allocatable :: q(:,:)
    allocate(q(0:m-1,m)); q=0.0_dp
    q(0,1:m)=uord(1:m)
    do t=1,m-1
      do j=t+1,m
        q(t,j)=self%pair(t,j)%hfunc1(q(t-1,t),q(t-1,j))
      end do
    end do
    value=q(m-1,m)
  end function conditional_last

  subroutine cvine_simulate(self,n,data)
    class(cvine_model), intent(in) :: self
    integer, intent(in) :: n
    real(dp), intent(out) :: data(:,:)
    real(dp), allocatable :: z(:,:)
    integer :: i,j
    allocate(z(self%dimension,n))
    do j=1,n
      do i=1,self%dimension
        z(i,j)=random_uniform()
      end do
    end do
    call self%inverse_rosenblatt(z,data)
  end subroutine cvine_simulate

  real(dp) function cvine_cdf(self,u,n_mc) result(value)
    class(cvine_model), intent(in) :: self
    real(dp), intent(in) :: u(:)
    integer, intent(in), optional :: n_mc
    integer :: n,i,j,count
    real(dp), allocatable :: sample(:,:)
    logical :: below
    n=10000; if (present(n_mc)) n=n_mc
    allocate(sample(self%dimension,n)); call self%simulate(n,sample)
    count=0
    do j=1,n
      below=.true.
      do i=1,self%dimension
        if (sample(i,j)>u(i)) then
          below=.false.; exit
        end if
      end do
      if (below) count=count+1
    end do
    value=real(count,dp)/real(n,dp)
  end function cvine_cdf

end module rvine_cvine
