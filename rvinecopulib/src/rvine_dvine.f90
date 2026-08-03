! SPDX-License-Identifier: GPL-3.0-only
! Native Fortran translation of the computational core of rvinecopulib.
module rvine_dvine
  use rvine_kinds, only : dp, eps_prob, clamp_prob
  use rvine_math, only : random_uniform
  use rvine_bicop, only : bicop_model, make_bicop, bicop_indep
  use rvine_fit, only : select_bicop, fit_bicop
  implicit none
  private
  public :: dvine_model, make_dvine

  type :: dvine_model
    integer :: dimension = 0
    integer, allocatable :: order(:)
    type(bicop_model), allocatable :: pair(:,:)
    integer :: nobs = 0
    real(dp) :: loglik = 0.0_dp
  contains
    procedure :: pdf => dvine_pdf
    procedure :: log_density => dvine_log_density
    procedure :: fit => dvine_fit
    procedure :: fit_fixed => dvine_fit_fixed
    procedure :: rosenblatt => dvine_rosenblatt
    procedure :: inverse_rosenblatt => dvine_inverse_rosenblatt
    procedure :: simulate => dvine_simulate
    procedure :: cdf => dvine_cdf
    procedure :: npars => dvine_npars
    procedure :: aic => dvine_aic
    procedure :: bic => dvine_bic
    procedure :: truncate => dvine_truncate
  end type dvine_model

contains

  function make_dvine(dimension, order) result(model)
    integer, intent(in) :: dimension
    integer, intent(in), optional :: order(:)
    type(dvine_model) :: model
    integer :: i, j
    model%dimension = dimension
    allocate(model%order(dimension),model%pair(dimension,dimension))
    if (present(order)) then
      model%order = order
    else
      model%order = [(i,i=1,dimension)]
    end if
    do i=1,dimension
      do j=1,dimension
        model%pair(i,j)=make_bicop(bicop_indep)
      end do
    end do
  end function make_dvine

  integer function dvine_npars(self) result(n)
    class(dvine_model), intent(in) :: self
    integer :: i,j
    n=0
    do i=1,self%dimension-1
      do j=i+1,self%dimension
        n=n+self%pair(i,j)%npar
      end do
    end do
  end function dvine_npars

  real(dp) function dvine_aic(self) result(value)
    class(dvine_model), intent(in) :: self
    value=-2.0_dp*self%loglik+2.0_dp*real(self%npars(),dp)
  end function dvine_aic

  real(dp) function dvine_bic(self) result(value)
    class(dvine_model), intent(in) :: self
    if (self%nobs>0) then
      value=-2.0_dp*self%loglik+log(real(self%nobs,dp))*real(self%npars(),dp)
    else
      value=huge(1.0_dp)
    end if
  end function dvine_bic

  subroutine dvine_truncate(self, level)
    class(dvine_model), intent(inout) :: self
    integer, intent(in) :: level
    integer :: i,j
    do i=1,self%dimension-1
      do j=i+1,self%dimension
        if (j-i>level) self%pair(i,j)=make_bicop(bicop_indep)
      end do
    end do
  end subroutine dvine_truncate

  real(dp) function dvine_pdf(self,u) result(value)
    class(dvine_model), intent(in) :: self
    real(dp), intent(in) :: u(:)
    real(dp) :: logd
    logd=self%log_density(u)
    if (logd>log(huge(1.0_dp))) then
      value=huge(1.0_dp)
    else
      value=exp(logd)
    end if
  end function dvine_pdf

  real(dp) function dvine_log_density(self,u) result(value)
    class(dvine_model), intent(in) :: self
    real(dp), intent(in) :: u(:)
    integer :: d,i,j,dist
    real(dp), allocatable :: left(:,:),right(:,:)
    real(dp) :: x,y,den
    d=self%dimension
    allocate(left(d,d),right(d,d))
    left=0.0_dp; right=0.0_dp
    do i=1,d
      left(i,i)=clamp_prob(u(self%order(i)))
      right(i,i)=left(i,i)
    end do
    value=0.0_dp
    do dist=1,d-1
      do i=1,d-dist
        j=i+dist
        x=left(i,j-1)
        y=right(i+1,j)
        den=self%pair(i,j)%pdf(x,y)
        value=value+log(max(tiny(1.0_dp),den))
        left(i,j)=self%pair(i,j)%hfunc2(x,y)
        right(i,j)=self%pair(i,j)%hfunc1(x,y)
      end do
    end do
  end function dvine_log_density

  subroutine dvine_fit(self,data,families,criterion,allow_rotations)
    class(dvine_model), intent(inout) :: self
    real(dp), intent(in) :: data(:,:)
    integer, intent(in), optional :: families(:)
    character(len=*), intent(in), optional :: criterion
    logical, intent(in), optional :: allow_rotations
    integer :: d,n,i,j,k,dist
    real(dp), allocatable :: left(:,:,:),right(:,:,:),pairdata(:,:)
    real(dp) :: x,y
    d=self%dimension; n=size(data,2)
    allocate(left(d,d,n),right(d,d,n),pairdata(2,n))
    left=0.0_dp; right=0.0_dp
    do i=1,d
      left(i,i,:)=data(self%order(i),:)
      right(i,i,:)=left(i,i,:)
    end do
    do dist=1,d-1
      do i=1,d-dist
        j=i+dist
        pairdata(1,:)=left(i,j-1,:)
        pairdata(2,:)=right(i+1,j,:)
        call select_dispatch(pairdata,self%pair(i,j),families,criterion,allow_rotations)
        do k=1,n
          x=pairdata(1,k); y=pairdata(2,k)
          left(i,j,k)=self%pair(i,j)%hfunc2(x,y)
          right(i,j,k)=self%pair(i,j)%hfunc1(x,y)
        end do
      end do
    end do
    self%nobs=n
    self%loglik=0.0_dp
    do k=1,n
      self%loglik=self%loglik+self%log_density(data(:,k))
    end do
  end subroutine dvine_fit

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
      else
        if (present(allow_rotations)) then
          call select_bicop(data,model,families,allow_rotations=allow_rotations)
        else
          call select_bicop(data,model,families)
        end if
      end if
    else
      if (present(criterion)) then
        if (present(allow_rotations)) then
          call select_bicop(data,model,criterion=criterion,allow_rotations=allow_rotations)
        else
          call select_bicop(data,model,criterion=criterion)
        end if
      else
        if (present(allow_rotations)) then
          call select_bicop(data,model,allow_rotations=allow_rotations)
        else
          call select_bicop(data,model)
        end if
      end if
    end if
  end subroutine select_dispatch

  subroutine dvine_fit_fixed(self,data)
    class(dvine_model), intent(inout) :: self
    real(dp), intent(in) :: data(:,:)
    integer :: d,n,i,j,k,dist,fam,rot
    real(dp), allocatable :: left(:,:,:),right(:,:,:),pairdata(:,:)
    real(dp) :: x,y
    type(bicop_model) :: fitted
    d=self%dimension; n=size(data,2)
    allocate(left(d,d,n),right(d,d,n),pairdata(2,n))
    left=0.0_dp; right=0.0_dp
    do i=1,d
      left(i,i,:)=data(self%order(i),:)
      right(i,i,:)=left(i,i,:)
    end do
    do dist=1,d-1
      do i=1,d-dist
        j=i+dist
        pairdata(1,:)=left(i,j-1,:)
        pairdata(2,:)=right(i+1,j,:)
        fam=self%pair(i,j)%family; rot=self%pair(i,j)%rotation
        call fit_bicop(pairdata,fam,rot,fitted)
        self%pair(i,j)=fitted
        do k=1,n
          x=pairdata(1,k); y=pairdata(2,k)
          left(i,j,k)=fitted%hfunc2(x,y)
          right(i,j,k)=fitted%hfunc1(x,y)
        end do
      end do
    end do
    self%nobs=n
    self%loglik=0.0_dp
    do k=1,n
      self%loglik=self%loglik+self%log_density(data(:,k))
    end do
  end subroutine dvine_fit_fixed

  subroutine dvine_rosenblatt(self,data,z)
    class(dvine_model), intent(in) :: self
    real(dp), intent(in) :: data(:,:)
    real(dp), intent(out) :: z(:,:)
    integer :: d,n,i,j,k,dist
    real(dp), allocatable :: left(:,:,:),right(:,:,:)
    real(dp) :: x,y
    d=self%dimension; n=size(data,2)
    allocate(left(d,d,n),right(d,d,n))
    left=0.0_dp; right=0.0_dp
    do i=1,d
      left(i,i,:)=data(self%order(i),:)
      right(i,i,:)=left(i,i,:)
    end do
    z=0.0_dp
    z(self%order(1),:)=right(1,1,:)
    do dist=1,d-1
      do i=1,d-dist
        j=i+dist
        do k=1,n
          x=left(i,j-1,k); y=right(i+1,j,k)
          left(i,j,k)=self%pair(i,j)%hfunc2(x,y)
          right(i,j,k)=self%pair(i,j)%hfunc1(x,y)
        end do
      end do
      j=dist+1
      z(self%order(j),:)=right(1,j,:)
    end do
    z=min(1.0_dp-eps_prob,max(eps_prob,z))
  end subroutine dvine_rosenblatt

  subroutine dvine_inverse_rosenblatt(self,z,data)
    class(dvine_model), intent(in) :: self
    real(dp), intent(in) :: z(:,:)
    real(dp), intent(out) :: data(:,:)
    integer :: d,n,k,j,iter
    real(dp), allocatable :: uord(:),word(:)
    real(dp) :: lo,hi,mid,target,cval
    d=self%dimension; n=size(z,2)
    allocate(uord(d),word(d))
    data=0.0_dp
    do k=1,n
      do j=1,d
        word(j)=clamp_prob(z(self%order(j),k))
      end do
      uord(1)=word(1)
      do j=2,d
        target=word(j); lo=eps_prob; hi=1.0_dp-eps_prob
        do iter=1,65
          mid=0.5_dp*(lo+hi)
          uord(j)=mid
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
  end subroutine dvine_inverse_rosenblatt

  real(dp) function conditional_last(self,uord,m) result(value)
    class(dvine_model), intent(in) :: self
    real(dp), intent(in) :: uord(:)
    integer, intent(in) :: m
    integer :: i,j,dist
    real(dp), allocatable :: left(:,:),right(:,:)
    real(dp) :: x,y
    allocate(left(m,m),right(m,m))
    left=0.0_dp; right=0.0_dp
    do i=1,m
      left(i,i)=uord(i); right(i,i)=uord(i)
    end do
    do dist=1,m-1
      do i=1,m-dist
        j=i+dist
        x=left(i,j-1); y=right(i+1,j)
        left(i,j)=self%pair(i,j)%hfunc2(x,y)
        right(i,j)=self%pair(i,j)%hfunc1(x,y)
      end do
    end do
    value=right(1,m)
  end function conditional_last

  subroutine dvine_simulate(self,n,data)
    class(dvine_model), intent(in) :: self
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
  end subroutine dvine_simulate

  real(dp) function dvine_cdf(self,u,n_mc) result(value)
    class(dvine_model), intent(in) :: self
    real(dp), intent(in) :: u(:)
    integer, intent(in), optional :: n_mc
    integer :: n,i,j,count
    real(dp), allocatable :: sample(:,:)
    logical :: below
    n=10000
    if (present(n_mc)) n=n_mc
    allocate(sample(self%dimension,n))
    call self%simulate(n,sample)
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
  end function dvine_cdf

end module rvine_dvine
