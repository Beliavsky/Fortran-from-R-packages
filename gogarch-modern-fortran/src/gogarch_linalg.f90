! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of gogarch, copyright (C) 2008-2026 Bernhard Pfaff.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
module gogarch_linalg
   use gogarch_kinds, only : dp
   implicit none
   private

   public :: identity_matrix, covariance_matrix, symmetric_eigen
   public :: inverse_matrix, determinant_matrix, symmetric_sqrt
   public :: symmetric_invsqrt, covariance_to_correlation, outer_product
   public :: is_symmetric, is_orthogonal, spd_logdet

   interface
      subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
         import :: dp
         character(len=1), intent(in) :: jobz, uplo
         integer, intent(in) :: n, lda, lwork
         real(dp), intent(inout) :: a(lda,*)
         real(dp), intent(out) :: w(*)
         real(dp), intent(inout) :: work(*)
         integer, intent(out) :: info
      end subroutine dsyev
      subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: dp
         integer, intent(in) :: n, nrhs, lda, ldb
         real(dp), intent(inout) :: a(lda,*), b(ldb,*)
         integer, intent(out) :: ipiv(*), info
      end subroutine dgesv
      subroutine dgetrf(m, n, a, lda, ipiv, info)
         import :: dp
         integer, intent(in) :: m, n, lda
         real(dp), intent(inout) :: a(lda,*)
         integer, intent(out) :: ipiv(*), info
      end subroutine dgetrf
   end interface

contains

   pure function identity_matrix(n) result(a)
      integer, intent(in) :: n
      real(dp) :: a(n,n)
      integer :: i
      a = 0.0_dp
      do i = 1, n
         a(i,i) = 1.0_dp
      end do
   end function identity_matrix

   pure function outer_product(x, y) result(a)
      real(dp), intent(in) :: x(:), y(:)
      real(dp) :: a(size(x),size(y))
      integer :: i, j
      do j = 1, size(y)
         do i = 1, size(x)
            a(i,j) = x(i)*y(j)
         end do
      end do
   end function outer_product

   function covariance_matrix(x, center) result(cov)
      real(dp), intent(in) :: x(:,:)
      logical, intent(in), optional :: center
      real(dp) :: cov(size(x,2),size(x,2))
      real(dp) :: work(size(x,1),size(x,2)), means(size(x,2))
      logical :: do_center
      integer :: n
      n = size(x,1)
      do_center = .false.
      if (present(center)) do_center = center
      work = x
      if (do_center .and. n > 0) then
         means = sum(x,dim=1)/real(n,dp)
         work = x-spread(means,1,n)
      end if
      if (n > 0) then
         cov = matmul(transpose(work),work)/real(n,dp)
      else
         cov = 0.0_dp
      end if
      cov = 0.5_dp*(cov+transpose(cov))
   end function covariance_matrix

   subroutine symmetric_eigen(a, values, vectors, ok)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: values(size(a,1)), vectors(size(a,1),size(a,2))
      logical, intent(out) :: ok
      real(dp), allocatable :: work(:)
      real(dp) :: query(1), tmp
      integer :: n, info, lwork, i, j
      n = size(a,1)
      if (size(a,2) /= n .or. n < 1) then
         values = 0.0_dp
         vectors = 0.0_dp
         ok = .false.
         return
      end if
      vectors = 0.5_dp*(a+transpose(a))
      call dsyev('V','U',n,vectors,n,values,query,-1,info)
      if (info /= 0) then
         values = 0.0_dp
         vectors = 0.0_dp
         ok = .false.
         return
      end if
      lwork = max(1,int(query(1)))
      allocate(work(lwork))
      call dsyev('V','U',n,vectors,n,values,work,lwork,info)
      ok = info == 0
      if (.not. ok) then
         values = 0.0_dp
         vectors = 0.0_dp
         return
      end if
      do i = 1, n/2
         j = n-i+1
         tmp = values(i); values(i) = values(j); values(j) = tmp
         work(1:n) = vectors(:,i)
         vectors(:,i) = vectors(:,j)
         vectors(:,j) = work(1:n)
      end do
   end subroutine symmetric_eigen

   function inverse_matrix(a, ok) result(ainv)
      real(dp), intent(in) :: a(:,:)
      logical, intent(out) :: ok
      real(dp) :: ainv(size(a,1),size(a,2)), acopy(size(a,1),size(a,2))
      integer :: n, ipiv(size(a,1)), info
      n = size(a,1)
      if (size(a,2) /= n .or. n < 1) then
         ainv = 0.0_dp
         ok = .false.
         return
      end if
      acopy = a
      ainv = identity_matrix(n)
      call dgesv(n,n,acopy,n,ipiv,ainv,n,info)
      ok = info == 0
      if (.not. ok) ainv = 0.0_dp
   end function inverse_matrix

   function determinant_matrix(a, ok) result(det)
      real(dp), intent(in) :: a(:,:)
      logical, intent(out), optional :: ok
      real(dp) :: det, lu(size(a,1),size(a,2))
      integer :: n, ipiv(size(a,1)), info, i, swaps
      n = size(a,1)
      if (size(a,2) /= n .or. n < 1) then
         det = 0.0_dp
         if (present(ok)) ok = .false.
         return
      end if
      lu = a
      call dgetrf(n,n,lu,n,ipiv,info)
      if (info /= 0) then
         det = 0.0_dp
         if (present(ok)) ok = .false.
         return
      end if
      det = 1.0_dp
      swaps = 0
      do i = 1, n
         det = det*lu(i,i)
         if (ipiv(i) /= i) swaps = swaps+1
      end do
      if (mod(swaps,2) == 1) det = -det
      if (present(ok)) ok = .true.
   end function determinant_matrix

   function symmetric_sqrt(a, floor_value, ok) result(root)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(in), optional :: floor_value
      logical, intent(out) :: ok
      real(dp) :: root(size(a,1),size(a,2)), values(size(a,1)), vectors(size(a,1),size(a,2))
      real(dp) :: floorv, scaled(size(a,1),size(a,2))
      integer :: j
      floorv = 0.0_dp
      if (present(floor_value)) floorv = floor_value
      call symmetric_eigen(a,values,vectors,ok)
      if (.not. ok) then
         root = 0.0_dp
         return
      end if
      scaled = vectors
      do j = 1, size(a,1)
         scaled(:,j) = scaled(:,j)*sqrt(max(values(j),floorv))
      end do
      root = matmul(scaled,transpose(vectors))
      root = 0.5_dp*(root+transpose(root))
   end function symmetric_sqrt

   function symmetric_invsqrt(a, floor_value, ok) result(root)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(in), optional :: floor_value
      logical, intent(out) :: ok
      real(dp) :: root(size(a,1),size(a,2)), values(size(a,1)), vectors(size(a,1),size(a,2))
      real(dp) :: floorv, scaled(size(a,1),size(a,2))
      integer :: j
      floorv = 1.0e-12_dp
      if (present(floor_value)) floorv = floor_value
      call symmetric_eigen(a,values,vectors,ok)
      if (.not. ok .or. minval(values) <= 0.0_dp) then
         root = 0.0_dp
         ok = .false.
         return
      end if
      scaled = vectors
      do j = 1, size(a,1)
         scaled(:,j) = scaled(:,j)/sqrt(max(values(j),floorv))
      end do
      root = matmul(scaled,transpose(vectors))
      root = 0.5_dp*(root+transpose(root))
   end function symmetric_invsqrt

   pure function covariance_to_correlation(cov) result(cor)
      real(dp), intent(in) :: cov(:,:)
      real(dp) :: cor(size(cov,1),size(cov,2)), denom
      integer :: i, j, n
      n = size(cov,1)
      cor = 0.0_dp
      do j = 1, n
         do i = 1, n
            denom = sqrt(max(cov(i,i),0.0_dp)*max(cov(j,j),0.0_dp))
            if (denom > tiny(1.0_dp)) cor(i,j) = cov(i,j)/denom
         end do
      end do
      do i = 1, n
         if (cov(i,i) > 0.0_dp) cor(i,i) = 1.0_dp
      end do
   end function covariance_to_correlation

   pure logical function is_symmetric(a, tol)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(in), optional :: tol
      real(dp) :: eps
      eps = 1.0e-10_dp
      if (present(tol)) eps = tol
      is_symmetric = size(a,1) == size(a,2)
      if (is_symmetric) is_symmetric = maxval(abs(a-transpose(a))) <= eps
   end function is_symmetric

   pure logical function is_orthogonal(a, tol)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(in), optional :: tol
      real(dp) :: eps
      eps = 1.0e-9_dp
      if (present(tol)) eps = tol
      is_orthogonal = size(a,1) == size(a,2)
      if (is_orthogonal) then
         is_orthogonal = maxval(abs(matmul(transpose(a),a)-identity_matrix(size(a,1)))) <= eps
      end if
   end function is_orthogonal

   function spd_logdet(a, ok) result(value)
      real(dp), intent(in) :: a(:,:)
      logical, intent(out) :: ok
      real(dp) :: value, eigenvalues(size(a,1)), vectors(size(a,1),size(a,2))
      call symmetric_eigen(a,eigenvalues,vectors,ok)
      if (.not. ok .or. minval(eigenvalues) <= 0.0_dp) then
         value = -huge(1.0_dp)
         ok = .false.
      else
         value = sum(log(eigenvalues))
      end if
   end function spd_logdet

end module gogarch_linalg
