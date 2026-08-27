! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2014-2025 Genaro Sucarrat (original R/C package)
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is distributed under the GNU General Public License version 2 only.
module lgarch_linalg
  use lgarch_kinds, only : dp
  use r_linalg, only : shared_cholesky_factor => cholesky_factor
  use r_linalg, only : shared_inverse_matrix => inverse_matrix
  use r_linalg, only : shared_solve_system => solve_system
  use r_linalg, only : shared_spectral_radius => spectral_radius
  use r_linalg, only : shared_symmetric_eigenvalues => symmetric_eigenvalues
  implicit none
  private
  public :: solve_linear, inverse_matrix, cholesky_lower, logdet_spd
  public :: symmetric_eigenvalues, spectral_radius, covariance_matrix, correlation_matrix
contains
  subroutine solve_linear(a,b,x,info)
    real(dp), intent(in) :: a(:,:), b(:)
    real(dp), intent(out) :: x(:)
    integer, intent(out) :: info
    integer :: n
    n=size(a,1)
    if (size(a,2)/=n .or. size(b)/=n .or. size(x)/=n) error stop "solve_linear: size mismatch"
    call shared_solve_system(a,b,x,info)
  end subroutine solve_linear

  subroutine inverse_matrix(a,ainv,info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: ainv(:,:)
    integer, intent(out) :: info
    integer :: n
    real(dp), allocatable :: shared_inverse(:,:)
    n=size(a,1)
    if(size(a,2)/=n .or. any(shape(ainv)/=shape(a))) error stop "inverse_matrix: size mismatch"
    call shared_inverse_matrix(a,shared_inverse,info)
    if (info==0) ainv=shared_inverse
  end subroutine inverse_matrix

  subroutine cholesky_lower(a,l,info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: l(:,:)
    integer, intent(out) :: info
    real(dp), allocatable :: shared_factor(:,:)
    call shared_cholesky_factor(a,shared_factor,info)
    if (info==0) l=shared_factor
  end subroutine cholesky_lower

  subroutine logdet_spd(a,logdet,info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: logdet
    integer, intent(out) :: info
    real(dp), allocatable :: l(:,:)
    integer :: i,n
    n=size(a,1); allocate(l(n,n)); call cholesky_lower(a,l,info)
    if(info/=0) then; logdet=huge(1.0_dp); return; end if
    logdet=0.0_dp
    do i=1,n; logdet=logdet+2.0_dp*log(l(i,i)); end do
  end subroutine logdet_spd

  subroutine symmetric_eigenvalues(a,w,info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: w(:)
    integer, intent(out) :: info
    real(dp), allocatable :: shared_values(:)
    call shared_symmetric_eigenvalues(a,shared_values,info)
    if (info==0) w=shared_values
  end subroutine symmetric_eigenvalues

  real(dp) function spectral_radius(phi) result(rho)
    real(dp), intent(in) :: phi(:)
    real(dp), allocatable :: a(:,:)
    integer :: p,info,i
    p=size(phi)
    if(p==0) then; rho=0.0_dp; return; end if
    allocate(a(p,p)); a=0.0_dp; a(1,:)=phi
    do i=2,p; a(i,i-1)=1.0_dp; end do
    call shared_spectral_radius(a,rho,info)
    if(info/=0) rho=huge(1.0_dp)
  end function spectral_radius

  function covariance_matrix(x) result(cov)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable :: cov(:,:)
    real(dp), allocatable :: xc(:,:)
    real(dp) :: means(size(x,2))
    integer :: j,n
    n=size(x,1); allocate(cov(size(x,2),size(x,2)),xc(n,size(x,2)))
    do j=1,size(x,2); means(j)=sum(x(:,j))/real(n,dp); xc(:,j)=x(:,j)-means(j); end do
    if(n>1) then; cov=matmul(transpose(xc),xc)/real(n-1,dp); else; cov=0.0_dp; end if
  end function covariance_matrix

  function correlation_matrix(x) result(cor)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable :: cor(:,:),cov(:,:)
    real(dp) :: d(size(x,2))
    integer :: i,j
    cov=covariance_matrix(x); allocate(cor(size(cov,1),size(cov,2)))
    do i=1,size(d); d(i)=sqrt(max(cov(i,i),tiny(1.0_dp))); end do
    do j=1,size(d); do i=1,size(d); cor(i,j)=cov(i,j)/(d(i)*d(j)); end do; end do
  end function correlation_matrix
end module lgarch_linalg
