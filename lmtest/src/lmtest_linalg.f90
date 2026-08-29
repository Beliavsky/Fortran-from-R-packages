module lmtest_linalg
   use lmtest_kinds, only : dp
   implicit none
   private
   public :: invert_spd, symmetric_eigenvalues, solve_linear, independent_union
   public :: covariance_matrix, sort_order, reorder_rows, reorder_vector

   interface
      subroutine dpotrf(uplo, n, a, lda, info)
         import dp
         character(len=1), intent(in) :: uplo
         integer, intent(in) :: n, lda
         real(dp), intent(inout) :: a(lda,*)
         integer, intent(out) :: info
      end subroutine dpotrf
      subroutine dpotri(uplo, n, a, lda, info)
         import dp
         character(len=1), intent(in) :: uplo
         integer, intent(in) :: n, lda
         real(dp), intent(inout) :: a(lda,*)
         integer, intent(out) :: info
      end subroutine dpotri
      subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
         import dp
         character(len=1), intent(in) :: jobz, uplo
         integer, intent(in) :: n, lda, lwork
         real(dp), intent(inout) :: a(lda,*)
         real(dp), intent(out) :: w(*)
         real(dp), intent(inout) :: work(*)
         integer, intent(out) :: info
      end subroutine dsyev
      subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
         import dp
         integer, intent(in) :: n, nrhs, lda, ldb
         real(dp), intent(inout) :: a(lda,*), b(ldb,*)
         integer, intent(out) :: ipiv(*), info
      end subroutine dgesv
   end interface

contains

   subroutine invert_spd(a, ainv, info)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: ainv(:,:)
      integer, intent(out) :: info
      integer :: n, i, j

      n = size(a, 1)
      allocate(ainv(n,n))
      if (size(a,2) /= n) then
         info = -1
         ainv = 0.0_dp
         return
      end if
      ainv = a
      call dpotrf('U', n, ainv, n, info)
      if (info /= 0) return
      call dpotri('U', n, ainv, n, info)
      if (info /= 0) return
      do j = 1, n
         do i = j + 1, n
            ainv(i,j) = ainv(j,i)
         end do
      end do
   end subroutine invert_spd

   subroutine symmetric_eigenvalues(a, values, info)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: values(:)
      integer, intent(out) :: info
      real(dp), allocatable :: acopy(:,:), work(:)
      real(dp) :: wq(1)
      integer :: n, lwork

      n = size(a,1)
      allocate(values(n), acopy(n,n))
      acopy = a
      call dsyev('N', 'U', n, acopy, n, values, wq, -1, info)
      if (info /= 0) return
      lwork = max(1, int(wq(1)))
      allocate(work(lwork))
      call dsyev('N', 'U', n, acopy, n, values, work, lwork, info)
   end subroutine symmetric_eigenvalues

   subroutine solve_linear(a, b, x, info)
      real(dp), intent(in) :: a(:,:), b(:)
      real(dp), allocatable, intent(out) :: x(:)
      real(dp), allocatable :: acopy(:,:), rhs(:,:)
      integer, allocatable :: ipiv(:)
      integer, intent(out) :: info
      integer :: n

      n = size(a,1)
      allocate(acopy(n,n), rhs(n,1), ipiv(n), x(n))
      acopy = a
      rhs(:,1) = b
      call dgesv(n, 1, acopy, n, ipiv, rhs, n, info)
      x = rhs(:,1)
   end subroutine solve_linear

   function covariance_matrix(x) result(cov)
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable :: cov(:,:)
      real(dp), allocatable :: xc(:,:)
      real(dp) :: means(size(x,2))
      integer :: n, j

      n = size(x,1)
      allocate(xc(n,size(x,2)), cov(size(x,2),size(x,2)))
      means = sum(x, dim=1) / real(n,dp)
      xc = x
      do j = 1, size(x,2)
         xc(:,j) = xc(:,j) - means(j)
      end do
      if (n > 1) then
         cov = matmul(transpose(xc), xc) / real(n-1,dp)
      else
         cov = 0.0_dp
      end if
   end function covariance_matrix

   subroutine independent_union(x, z, u, tol)
      real(dp), intent(in) :: x(:,:), z(:,:)
      real(dp), allocatable, intent(out) :: u(:,:)
      real(dp), intent(in), optional :: tol
      real(dp), allocatable :: cand(:,:), q(:,:), v(:)
      real(dp) :: eps, nrm, scale
      integer :: n, p, j, r

      n = size(x,1)
      if (size(z,1) /= n) then
         allocate(u(0,0))
         return
      end if
      p = size(x,2) + size(z,2)
      allocate(cand(n,p))
      cand(:,1:size(x,2)) = x
      cand(:,size(x,2)+1:p) = z
      allocate(q(n,p), v(n))
      q = 0.0_dp
      r = 0
      eps = 1.0e-10_dp
      if (present(tol)) eps = tol
      do j = 1, p
         v = cand(:,j)
         scale = max(1.0_dp, sqrt(sum(v*v)))
         if (r > 0) v = v - matmul(q(:,1:r), matmul(transpose(q(:,1:r)), v))
         nrm = sqrt(sum(v*v))
         if (nrm > eps * scale) then
            r = r + 1
            q(:,r) = v / nrm
         end if
      end do
      allocate(u(n,r))
      r = 0
      q = 0.0_dp
      do j = 1, p
         v = cand(:,j)
         scale = max(1.0_dp, sqrt(sum(v*v)))
         if (r > 0) v = v - matmul(q(:,1:r), matmul(transpose(q(:,1:r)), v))
         nrm = sqrt(sum(v*v))
         if (nrm > eps * scale) then
            r = r + 1
            q(:,r) = v / nrm
            u(:,r) = cand(:,j)
         end if
      end do
   end subroutine independent_union

   function sort_order(x) result(idx)
      real(dp), intent(in) :: x(:)
      integer, allocatable :: idx(:)
      integer :: i, j, key
      allocate(idx(size(x)))
      idx = [(i, i=1,size(x))]
      do i = 2, size(x)
         key = idx(i)
         j = i - 1
         do while (j >= 1)
            if (x(idx(j)) <= x(key)) exit
            idx(j+1) = idx(j)
            j = j - 1
         end do
         idx(j+1) = key
      end do
   end function sort_order

   function reorder_rows(x, idx) result(y)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: idx(:)
      real(dp), allocatable :: y(:,:)
      integer :: i
      allocate(y(size(idx),size(x,2)))
      do i = 1, size(idx)
         y(i,:) = x(idx(i),:)
      end do
   end function reorder_rows

   function reorder_vector(x, idx) result(y)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: idx(:)
      real(dp), allocatable :: y(:)
      integer :: i
      allocate(y(size(idx)))
      do i = 1, size(idx)
         y(i) = x(idx(i))
      end do
   end function reorder_vector

end module lmtest_linalg
