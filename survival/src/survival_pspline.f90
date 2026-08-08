! survival pspline translation: LGPL-2.0-or-later
! Links to the separately licensed GPL-2.0-or-later splines-fortran dependency.
module survival_pspline
  use survival_kinds, only : dp
  use splines_core, only : spline_design
  implicit none
  private
  public :: pspline_basis
contains
  subroutine pspline_basis(x,nterm,degree,boundary,intercept,basis,penalty,knots,status)
    real(dp), intent(in) :: x(:), boundary(2)
    integer, intent(in) :: nterm, degree
    logical, intent(in) :: intercept
    real(dp), allocatable, intent(out) :: basis(:,:), penalty(:,:), knots(:)
    integer, intent(out), optional :: status
    real(dp), allocatable :: full_basis(:,:), d2(:,:), full_penalty(:,:)
    real(dp) :: dx
    integer :: nvar, i, ierr

    if (present(status)) status=0
    if (nterm<3 .or. degree<1 .or. boundary(1)>=boundary(2)) then
      allocate(basis(0,0),penalty(0,0),knots(0))
      if(present(status)) status=1
      return
    end if
    dx=(boundary(2)-boundary(1))/real(nterm,dp)
    allocate(knots(nterm+2*degree+1))
    do i=0,nterm+degree-1
      knots(i+1)=boundary(1)+dx*real(i-degree,dp)
    end do
    do i=0,degree
      knots(nterm+degree+i+1)=boundary(2)+dx*real(i,dp)
    end do
    call spline_design(knots,x,degree+1,full_basis,status=ierr)
    if(ierr/=0) then
      allocate(basis(0,0),penalty(0,0))
      if(present(status)) status=10+ierr
      return
    end if
    nvar=size(full_basis,2)
    allocate(d2(max(0,nvar-2),nvar))
    d2=0.0_dp
    do i=1,nvar-2
      d2(i,i)=1.0_dp
      d2(i,i+1)=-2.0_dp
      d2(i,i+2)=1.0_dp
    end do
    full_penalty=matmul(transpose(d2),d2)
    if(intercept) then
      allocate(basis(size(x),nvar),penalty(nvar,nvar))
      basis=full_basis
      penalty=full_penalty
    else
      allocate(basis(size(x),nvar-1),penalty(nvar-1,nvar-1))
      basis=full_basis(:,2:nvar)
      penalty=full_penalty(2:nvar,2:nvar)
    end if
  end subroutine pspline_basis
end module survival_pspline
