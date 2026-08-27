! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from tsDyn, copyright its original authors.
! This file may be redistributed and/or modified under GPL version 2 or later.
module tsdyn_linalg
  use tsdyn_kinds, only: dp
  use r_linalg, only: shared_cholesky_factor => cholesky_factor
  use r_linalg, only: shared_general_eigen => general_real_eigen
  use r_linalg, only: shared_general_eigenvalues => general_real_eigenvalues
  use r_linalg, only: shared_inverse_matrix => inverse_matrix
  use r_linalg, only: shared_least_squares_svd => least_squares_svd
  use r_linalg, only: shared_symmetric_eigen => symmetric_eigen
  implicit none
  private
  public :: ols_fit, inverse_matrix, cholesky_lower, symmetric_eigen
  public :: general_eigenvalues, general_eigen, matrix_rank, covariance_matrix
contains
  subroutine ols_fit(x, y, beta, fitted, residuals, rank, ssr, info)
    real(dp), intent(in) :: x(:,:), y(:,:)
    real(dp), allocatable, intent(out) :: beta(:,:), fitted(:,:), residuals(:,:)
    integer, intent(out) :: rank, info
    real(dp), intent(out) :: ssr
    integer :: n, p, ny

    n = size(x,1); p = size(x,2); ny = size(y,2)
    if (size(y,1) /= n .or. n < 1 .or. p < 1) then
      info = -1; rank = 0; ssr = huge(1.0_dp)
      allocate(beta(0,0), fitted(0,0), residuals(0,0))
      return
    end if
    allocate(beta(p,ny), fitted(n,ny), residuals(n,ny))
    call shared_least_squares_svd(x,y,beta,rank,info,rcond=1.0e-12_dp)
    if (info /= 0) then
      beta = 0.0_dp; fitted = 0.0_dp; residuals = y
      ssr = huge(1.0_dp)
      return
    end if
    fitted = matmul(x,beta)
    residuals = y - fitted
    ssr = sum(residuals*residuals)
  end subroutine ols_fit

  subroutine inverse_matrix(a, ainv, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: ainv(:,:)
    integer, intent(out) :: info
    integer :: n
    n = size(a,1)
    if (size(a,2) /= n .or. n < 1) then
      info = -1; allocate(ainv(0,0)); return
    end if
    call shared_inverse_matrix(a,ainv,info)
  end subroutine inverse_matrix

  subroutine cholesky_lower(a, l, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: l(:,:)
    integer, intent(out) :: info
    integer :: n
    n=size(a,1)
    if(size(a,2)/=n)then
      info=-1; allocate(l(0,0)); return
    end if
    call shared_cholesky_factor(a,l,info)
  end subroutine cholesky_lower

  subroutine symmetric_eigen(a, values, vectors, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: values(:), vectors(:,:)
    integer, intent(out) :: info
    integer :: n
    n=size(a,1)
    if(size(a,2)/=n)then
      info=-1; allocate(values(0),vectors(0,0)); return
    end if
    call shared_symmetric_eigen(a,values,vectors,info)
  end subroutine symmetric_eigen

  subroutine general_eigenvalues(a, wr, wi, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: wr(:), wi(:)
    integer, intent(out) :: info
    integer :: n
    n=size(a,1)
    if(size(a,2)/=n)then
      info=-1; allocate(wr(0),wi(0)); return
    end if
    call shared_general_eigenvalues(a,wr,wi,info)
  end subroutine general_eigenvalues

  subroutine general_eigen(a, wr, wi, vr, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: wr(:), wi(:), vr(:,:)
    integer, intent(out) :: info
    integer :: n
    n=size(a,1)
    if(size(a,2)/=n)then
      info=-1; allocate(wr(0),wi(0),vr(0,0)); return
    end if
    call shared_general_eigen(a,wr,wi,vr,info)
  end subroutine general_eigen

  integer function matrix_rank(a, tol) result(r)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(in), optional :: tol
    real(dp), allocatable :: beta(:,:),fit(:,:),res(:,:), y(:,:)
    real(dp)::ss
    integer::info
    allocate(y(size(a,1),1)); y=0.0_dp
    call ols_fit(a,y,beta,fit,res,r,ss,info)
    if(present(tol)) r = r + merge(0,0,tol>=0.0_dp)
  end function matrix_rank

  subroutine covariance_matrix(x, cov, center)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable, intent(out) :: cov(:,:)
    logical, intent(in), optional :: center
    real(dp), allocatable :: xc(:,:),mu(:)
    logical :: do_center
    integer :: n,j
    n=size(x,1); do_center=.true.; if(present(center))do_center=center
    allocate(xc(size(x,1),size(x,2))); xc=x
    if(do_center)then
      allocate(mu(size(x,2))); mu=sum(x,dim=1)/real(max(1,n),dp)
      do j=1,size(x,2); xc(:,j)=xc(:,j)-mu(j); end do
    end if
    allocate(cov(size(x,2),size(x,2)))
    cov=matmul(transpose(xc),xc)/real(max(1,n-1),dp)
  end subroutine covariance_matrix
end module tsdyn_linalg
