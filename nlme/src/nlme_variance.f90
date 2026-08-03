! SPDX-License-Identifier: GPL-2.0-or-later
module nlme_variance
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use nlme_kinds, only : dp
  use nlme_status, only : NLME_SUCCESS, NLME_INVALID_ARGUMENT, NLME_DIMENSION_ERROR
  use nlme_types, only : variance_spec, VAR_CONSTANT, VAR_FIXED, VAR_IDENT, VAR_POWER, &
       VAR_EXPONENTIAL, VAR_CONST_POWER, VAR_CONST_PROP
  implicit none
  private
  public :: variance_sd, variance_parameter_count
  public :: variance_to_unconstrained, variance_from_unconstrained
contains
  pure integer function variance_parameter_count(spec, nlevels) result(k)
    type(variance_spec), intent(in) :: spec
    integer, intent(in), optional :: nlevels
    integer :: nl
    nl=0; if (present(nlevels)) nl=nlevels
    select case(spec%kind)
    case(VAR_CONSTANT,VAR_FIXED)
      k=0
    case(VAR_IDENT)
      k=max(0,nl-1)
    case(VAR_POWER,VAR_EXPONENTIAL)
      k=1
    case(VAR_CONST_POWER,VAR_CONST_PROP)
      k=2
    case default
      k=0
    end select
  end function variance_parameter_count

  subroutine variance_sd(spec,covariate,group_index,sd,status)
    type(variance_spec), intent(in) :: spec
    real(dp), intent(in) :: covariate(:)
    integer, intent(in), optional :: group_index(:)
    real(dp), allocatable, intent(out) :: sd(:)
    integer, intent(out) :: status
    integer :: n,i,g
    real(dp) :: v,c,p
    n=size(covariate); allocate(sd(n)); sd=1.0_dp
    select case(spec%kind)
    case(VAR_CONSTANT)
      status=NLME_SUCCESS
    case(VAR_FIXED)
      do i=1,n
        v=abs(covariate(i))
        if (.not.ieee_is_finite(v) .or. v<=0.0_dp) then
          status=NLME_INVALID_ARGUMENT; return
        end if
        sd(i)=sqrt(v)
      end do
      status=NLME_SUCCESS
    case(VAR_IDENT)
      if (.not.present(group_index)) then; status=NLME_INVALID_ARGUMENT; return; end if
      if (size(group_index)/=n) then; status=NLME_DIMENSION_ERROR; return; end if
      do i=1,n
        g=group_index(i)
        if (g<1) then; status=NLME_INVALID_ARGUMENT; return; end if
        if (g==1) then
          sd(i)=1.0_dp
        else
          if (.not.allocated(spec%par) .or. g-1>size(spec%par)) then
            status=NLME_INVALID_ARGUMENT; return
          end if
          sd(i)=exp(spec%par(g-1))
        end if
      end do
      status=NLME_SUCCESS
    case(VAR_POWER)
      if (.not.allocated(spec%par) .or. size(spec%par)<1) then; status=NLME_INVALID_ARGUMENT; return; end if
      p=spec%par(1)
      do i=1,n
        sd(i)=max(abs(covariate(i)),sqrt(tiny(1.0_dp)))**p
      end do
      status=NLME_SUCCESS
    case(VAR_EXPONENTIAL)
      if (.not.allocated(spec%par) .or. size(spec%par)<1) then; status=NLME_INVALID_ARGUMENT; return; end if
      p=spec%par(1)
      sd=exp(max(-300.0_dp,min(300.0_dp,p*covariate)))
      status=NLME_SUCCESS
    case(VAR_CONST_POWER)
      if (.not.allocated(spec%par) .or. size(spec%par)<2) then; status=NLME_INVALID_ARGUMENT; return; end if
      c=exp(spec%par(1)); p=spec%par(2)
      do i=1,n
        sd(i)=c+max(abs(covariate(i)),sqrt(tiny(1.0_dp)))**p
      end do
      status=NLME_SUCCESS
    case(VAR_CONST_PROP)
      if (.not.allocated(spec%par) .or. size(spec%par)<2) then; status=NLME_INVALID_ARGUMENT; return; end if
      c=exp(spec%par(1)); p=exp(spec%par(2))
      do i=1,n
        sd(i)=sqrt(c*c+(p*abs(covariate(i)))**2)
      end do
      status=NLME_SUCCESS
    case default
      status=NLME_INVALID_ARGUMENT
    end select
    if (status==NLME_SUCCESS .and. any(.not.ieee_is_finite(sd))) status=NLME_INVALID_ARGUMENT
  end subroutine variance_sd

  subroutine variance_to_unconstrained(spec,nlevels,x,status)
    type(variance_spec), intent(in) :: spec
    integer, intent(in) :: nlevels
    real(dp), allocatable, intent(out) :: x(:)
    integer, intent(out) :: status
    integer :: k
    k=variance_parameter_count(spec,nlevels); allocate(x(k))
    if (k==0) then; status=NLME_SUCCESS; return; end if
    if (.not.allocated(spec%par) .or. size(spec%par)/=k) then
      status=NLME_INVALID_ARGUMENT; return
    end if
    x=spec%par; status=NLME_SUCCESS
  end subroutine variance_to_unconstrained

  subroutine variance_from_unconstrained(template,nlevels,x,spec,status)
    type(variance_spec), intent(in) :: template
    integer, intent(in) :: nlevels
    real(dp), intent(in) :: x(:)
    type(variance_spec), intent(out) :: spec
    integer, intent(out) :: status
    integer :: k
    spec=template; k=variance_parameter_count(template,nlevels)
    if (size(x)/=k) then; status=NLME_DIMENSION_ERROR; return; end if
    if (allocated(spec%par)) deallocate(spec%par); allocate(spec%par(k)); spec%par=x
    status=NLME_SUCCESS
  end subroutine variance_from_unconstrained
end module nlme_variance
