! SPDX-License-Identifier: GPL-2.0-or-later
module nlme_pdmat
  use nlme_kinds, only : dp
  use nlme_status, only : NLME_SUCCESS, NLME_INVALID_ARGUMENT, NLME_DIMENSION_ERROR
  use nlme_types, only : pd_spec, PD_IDENT, PD_DIAG, PD_LOG_CHOL, PD_COMPOUND_SYMM
  use nlme_linalg, only : symmetrize, cholesky_lower
  implicit none
  private
  public :: pd_matrix, pd_parameter_count, pd_from_matrix, pd_to_unconstrained
contains
  pure integer function pd_parameter_count(spec) result(k)
    type(pd_spec), intent(in) :: spec
    select case(spec%kind)
    case(PD_IDENT)
      k=1
    case(PD_DIAG)
      k=max(0,spec%dim)
    case(PD_LOG_CHOL)
      k=max(0,spec%dim*(spec%dim+1)/2)
    case(PD_COMPOUND_SYMM)
      k=2
    case default
      k=0
    end select
  end function pd_parameter_count

  subroutine pd_matrix(spec,matrix,status)
    type(pd_spec), intent(in) :: spec
    real(dp), allocatable, intent(out) :: matrix(:,:)
    integer, intent(out) :: status
    integer :: n,i,j,k
    real(dp) :: sd,rho,lower
    real(dp), allocatable :: l(:,:)
    n=spec%dim
    if (n<1 .or. .not.allocated(spec%par) .or. size(spec%par)/=pd_parameter_count(spec)) then
      allocate(matrix(0,0)); status=NLME_INVALID_ARGUMENT; return
    end if
    allocate(matrix(n,n)); matrix=0.0_dp
    select case(spec%kind)
    case(PD_IDENT)
      sd=exp(spec%par(1))
      do i=1,n; matrix(i,i)=sd*sd; end do
    case(PD_DIAG)
      do i=1,n; matrix(i,i)=exp(2.0_dp*spec%par(i)); end do
    case(PD_LOG_CHOL)
      allocate(l(n,n)); l=0.0_dp; k=0
      do i=1,n
        k=k+1; l(i,i)=exp(spec%par(k))
        do j=1,i-1
          k=k+1; l(i,j)=spec%par(k)
        end do
      end do
      matrix=matmul(l,transpose(l)); matrix=symmetrize(matrix)
    case(PD_COMPOUND_SYMM)
      sd=exp(spec%par(1)); lower=-1.0_dp/real(max(1,n-1),dp)
      rho=lower+(1.0_dp-lower)/(1.0_dp+exp(-max(-40.0_dp,min(40.0_dp,spec%par(2)))))
      do i=1,n
        matrix(i,i)=sd*sd
        do j=1,i-1
          matrix(i,j)=rho*sd*sd; matrix(j,i)=matrix(i,j)
        end do
      end do
    case default
      status=NLME_INVALID_ARGUMENT; return
    end select
    status=NLME_SUCCESS
  end subroutine pd_matrix

  subroutine pd_from_matrix(a,kind,spec,status)
    real(dp), intent(in) :: a(:,:)
    integer, intent(in) :: kind
    type(pd_spec), intent(out) :: spec
    integer, intent(out) :: status
    integer :: n,i,j,k
    real(dp), allocatable :: l(:,:)
    n=size(a,1)
    if (size(a,2)/=n .or. n<1) then; status=NLME_DIMENSION_ERROR; return; end if
    spec%kind=kind; spec%dim=n
    select case(kind)
    case(PD_IDENT)
      allocate(spec%par(1)); spec%par(1)=0.5_dp*log(max(tiny(1.0_dp),sum([(a(i,i),i=1,n)])/real(n,dp)))
    case(PD_DIAG)
      allocate(spec%par(n))
      do i=1,n; spec%par(i)=0.5_dp*log(max(tiny(1.0_dp),a(i,i))); end do
    case(PD_LOG_CHOL)
      call cholesky_lower(a,l,status)
      if (status/=NLME_SUCCESS) return
      allocate(spec%par(n*(n+1)/2)); k=0
      do i=1,n
        k=k+1; spec%par(k)=log(l(i,i))
        do j=1,i-1
          k=k+1; spec%par(k)=l(i,j)
        end do
      end do
    case(PD_COMPOUND_SYMM)
      allocate(spec%par(2)); spec%par(1)=0.5_dp*log(max(tiny(1.0_dp),sum([(a(i,i),i=1,n)])/real(n,dp)))
      if (n>1) then
        spec%par(2)=0.0_dp
      else
        spec%par(2)=0.0_dp
      end if
    case default
      status=NLME_INVALID_ARGUMENT; return
    end select
    status=NLME_SUCCESS
  end subroutine pd_from_matrix

  subroutine pd_to_unconstrained(spec,x,status)
    type(pd_spec), intent(in) :: spec
    real(dp), allocatable, intent(out) :: x(:)
    integer, intent(out) :: status
    if (.not.allocated(spec%par) .or. size(spec%par)/=pd_parameter_count(spec)) then
      allocate(x(0)); status=NLME_INVALID_ARGUMENT; return
    end if
    allocate(x(size(spec%par))); x=spec%par; status=NLME_SUCCESS
  end subroutine pd_to_unconstrained
end module nlme_pdmat
