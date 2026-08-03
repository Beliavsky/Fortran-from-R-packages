! SPDX-License-Identifier: GPL-2.0-only
module fincov_linalg
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use fincov_kinds, only : dp
   use fincov_status, only : fincov_ok, fincov_invalid_input, fincov_size_mismatch, &
      fincov_singular_matrix, fincov_no_convergence, fincov_allocation_failure
   implicit none
   private

   public :: sample_covariance, column_variances, center_columns
   public :: solve_linear_system, symmetric_eigen_jacobi
   public :: frobenius_norm_squared, spectral_norm_squared, minimum_eigenvalue
   public :: project_simplex, matrix_is_symmetric, add_diagonal

   interface solve_linear_system
      module procedure solve_matrix_rhs
      module procedure solve_vector_rhs
   end interface solve_linear_system
contains
   subroutine sample_covariance(x, covariance, status)
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable, intent(out) :: covariance(:,:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: centered(:,:)
      integer :: n, p, alloc_stat

      n = size(x, 1)
      p = size(x, 2)
      if (n < 2 .or. p < 1) then
         allocate(covariance(0,0))
         if (present(status)) status = fincov_invalid_input
         return
      end if

      allocate(centered(n,p), covariance(p,p), stat=alloc_stat)
      if (alloc_stat /= 0) then
         if (allocated(covariance)) deallocate(covariance)
         allocate(covariance(0,0))
         if (present(status)) status = fincov_allocation_failure
         return
      end if
      call center_columns(x, centered)
      covariance = matmul(transpose(centered), centered) / real(n - 1, dp)
      covariance = 0.5_dp * (covariance + transpose(covariance))
      if (present(status)) status = fincov_ok
   end subroutine sample_covariance

   subroutine center_columns(x, centered)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(out) :: centered(:,:)
      real(dp) :: means(size(x,2))
      integer :: j

      if (size(centered,1) /= size(x,1) .or. size(centered,2) /= size(x,2)) return
      means = sum(x, dim=1) / real(size(x,1), dp)
      do j = 1, size(x,2)
         centered(:,j) = x(:,j) - means(j)
      end do
   end subroutine center_columns

   subroutine column_variances(x, variances, status)
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable, intent(out) :: variances(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: centered(:,:)
      integer :: n, p, j, alloc_stat

      n = size(x,1)
      p = size(x,2)
      if (n < 2 .or. p < 1) then
         allocate(variances(0))
         if (present(status)) status = fincov_invalid_input
         return
      end if
      allocate(centered(n,p), variances(p), stat=alloc_stat)
      if (alloc_stat /= 0) then
         if (allocated(variances)) deallocate(variances)
         allocate(variances(0))
         if (present(status)) status = fincov_allocation_failure
         return
      end if
      call center_columns(x, centered)
      do j = 1, p
         variances(j) = dot_product(centered(:,j), centered(:,j)) / real(n - 1, dp)
      end do
      if (present(status)) status = fincov_ok
   end subroutine column_variances

   subroutine solve_matrix_rhs(a, b, x, status)
      real(dp), intent(in) :: a(:,:), b(:,:)
      real(dp), allocatable, intent(out) :: x(:,:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: aa(:,:), bb(:,:)
      real(dp) :: factor, pivot_value, scale, tmp
      integer :: n, nrhs, i, j, k, pivot, alloc_stat

      n = size(a,1)
      nrhs = size(b,2)
      if (n < 1 .or. size(a,2) /= n .or. size(b,1) /= n .or. nrhs < 1) then
         allocate(x(0,0))
         if (present(status)) status = fincov_size_mismatch
         return
      end if
      allocate(aa(n,n), bb(n,nrhs), x(n,nrhs), stat=alloc_stat)
      if (alloc_stat /= 0) then
         if (allocated(x)) deallocate(x)
         allocate(x(0,0))
         if (present(status)) status = fincov_allocation_failure
         return
      end if
      aa = a
      bb = b
      scale = max(1.0_dp, maxval(abs(aa)))

      do k = 1, n
         pivot = k
         pivot_value = abs(aa(k,k))
         do i = k + 1, n
            if (abs(aa(i,k)) > pivot_value) then
               pivot = i
               pivot_value = abs(aa(i,k))
            end if
         end do
         if (pivot_value <= 100.0_dp * epsilon(1.0_dp) * scale) then
            x = 0.0_dp
            if (present(status)) status = fincov_singular_matrix
            return
         end if
         if (pivot /= k) then
            do j = 1, n
               tmp = aa(k,j)
               aa(k,j) = aa(pivot,j)
               aa(pivot,j) = tmp
            end do
            do j = 1, nrhs
               tmp = bb(k,j)
               bb(k,j) = bb(pivot,j)
               bb(pivot,j) = tmp
            end do
         end if

         do i = k + 1, n
            factor = aa(i,k) / aa(k,k)
            aa(i,k) = 0.0_dp
            if (k < n) aa(i,k+1:n) = aa(i,k+1:n) - factor * aa(k,k+1:n)
            bb(i,:) = bb(i,:) - factor * bb(k,:)
         end do
      end do

      x = 0.0_dp
      do i = n, 1, -1
         if (i < n) then
            x(i,:) = (bb(i,:) - matmul(aa(i,i+1:n), x(i+1:n,:))) / aa(i,i)
         else
            x(i,:) = bb(i,:) / aa(i,i)
         end if
      end do
      if (present(status)) status = fincov_ok
   end subroutine solve_matrix_rhs

   subroutine solve_vector_rhs(a, b, x, status)
      real(dp), intent(in) :: a(:,:), b(:)
      real(dp), allocatable, intent(out) :: x(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: bmat(:,:), xmat(:,:)
      integer :: local_status, alloc_stat

      if (size(a,1) /= size(b)) then
         allocate(x(0))
         if (present(status)) status = fincov_size_mismatch
         return
      end if
      allocate(bmat(size(b),1), stat=alloc_stat)
      if (alloc_stat /= 0) then
         allocate(x(0))
         if (present(status)) status = fincov_allocation_failure
         return
      end if
      bmat(:,1) = b
      call solve_matrix_rhs(a, bmat, xmat, local_status)
      if (local_status == fincov_ok) then
         allocate(x(size(b)), stat=alloc_stat)
         if (alloc_stat == 0) then
            x = xmat(:,1)
         else
            allocate(x(0))
            local_status = fincov_allocation_failure
         end if
      else
         allocate(x(0))
      end if
      if (present(status)) status = local_status
   end subroutine solve_vector_rhs

   subroutine symmetric_eigen_jacobi(a, values, vectors, status, tolerance, max_iterations)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: values(:), vectors(:,:)
      integer, intent(out), optional :: status
      real(dp), intent(in), optional :: tolerance
      integer, intent(in), optional :: max_iterations
      real(dp), allocatable :: work(:,:), column_tmp(:)
      real(dp) :: tol, offmax, app, aqq, apq, tau, t, c, s, aik, aqk, vip, viq, tmp
      integer :: n, i, j, k, p, q, iteration, max_iter, alloc_stat, best

      n = size(a,1)
      if (n < 1 .or. size(a,2) /= n) then
         allocate(values(0), vectors(0,0))
         if (present(status)) status = fincov_size_mismatch
         return
      end if
      if (.not. matrix_is_symmetric(a)) then
         allocate(values(0), vectors(0,0))
         if (present(status)) status = fincov_invalid_input
         return
      end if
      allocate(work(n,n), values(n), vectors(n,n), column_tmp(n), stat=alloc_stat)
      if (alloc_stat /= 0) then
         if (allocated(values)) deallocate(values)
         if (allocated(vectors)) deallocate(vectors)
         allocate(values(0), vectors(0,0))
         if (present(status)) status = fincov_allocation_failure
         return
      end if
      work = 0.5_dp * (a + transpose(a))
      vectors = 0.0_dp
      do i = 1, n
         vectors(i,i) = 1.0_dp
      end do
      tol = 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(work)))
      if (present(tolerance)) tol = max(tolerance, epsilon(1.0_dp))
      max_iter = max(100, 50*n*n)
      if (present(max_iterations)) max_iter = max(1, max_iterations)

      if (n == 1) then
         values(1) = work(1,1)
         if (present(status)) status = fincov_ok
         return
      end if

      do iteration = 1, max_iter
         offmax = 0.0_dp
         p = 1
         q = 2
         do i = 1, n - 1
            do j = i + 1, n
               if (abs(work(i,j)) > offmax) then
                  offmax = abs(work(i,j))
                  p = i
                  q = j
               end if
            end do
         end do
         if (offmax <= tol) exit

         app = work(p,p)
         aqq = work(q,q)
         apq = work(p,q)
         tau = (aqq - app) / (2.0_dp * apq)
         if (tau >= 0.0_dp) then
            t = 1.0_dp / (tau + sqrt(1.0_dp + tau*tau))
         else
            t = -1.0_dp / (-tau + sqrt(1.0_dp + tau*tau))
         end if
         c = 1.0_dp / sqrt(1.0_dp + t*t)
         s = t*c

         do k = 1, n
            if (k /= p .and. k /= q) then
               aik = work(k,p)
               aqk = work(k,q)
               work(k,p) = c*aik - s*aqk
               work(p,k) = work(k,p)
               work(k,q) = s*aik + c*aqk
               work(q,k) = work(k,q)
            end if
         end do
         work(p,p) = app - t*apq
         work(q,q) = aqq + t*apq
         work(p,q) = 0.0_dp
         work(q,p) = 0.0_dp

         do k = 1, n
            vip = vectors(k,p)
            viq = vectors(k,q)
            vectors(k,p) = c*vip - s*viq
            vectors(k,q) = s*vip + c*viq
         end do
      end do

      values = [(work(i,i), i=1,n)]
      do i = 1, n - 1
         best = i
         do j = i + 1, n
            if (values(j) > values(best)) best = j
         end do
         if (best /= i) then
            tmp = values(i)
            values(i) = values(best)
            values(best) = tmp
            column_tmp = vectors(:,i)
            vectors(:,i) = vectors(:,best)
            vectors(:,best) = column_tmp
         end if
      end do

      if (present(status)) then
         if (iteration > max_iter .and. offmax > tol) then
            status = fincov_no_convergence
         else
            status = fincov_ok
         end if
      end if
   end subroutine symmetric_eigen_jacobi

   pure function frobenius_norm_squared(a) result(value)
      real(dp), intent(in) :: a(:,:)
      real(dp) :: value
      value = sum(a*a)
   end function frobenius_norm_squared

   function spectral_norm_squared(a, status) result(value)
      real(dp), intent(in) :: a(:,:)
      integer, intent(out), optional :: status
      real(dp) :: value
      real(dp), allocatable :: gram(:,:), eigenvalues(:), eigenvectors(:,:)
      integer :: local_status

      if (size(a,1) < 1 .or. size(a,2) < 1) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         if (present(status)) status = fincov_invalid_input
         return
      end if
      gram = matmul(transpose(a), a)
      call symmetric_eigen_jacobi(gram, eigenvalues, eigenvectors, local_status)
      if (local_status == fincov_ok .or. local_status == fincov_no_convergence) then
         value = max(0.0_dp, maxval(eigenvalues))
      else
         value = ieee_value(0.0_dp, ieee_quiet_nan)
      end if
      if (present(status)) status = local_status
   end function spectral_norm_squared

   function minimum_eigenvalue(a, status) result(value)
      real(dp), intent(in) :: a(:,:)
      integer, intent(out), optional :: status
      real(dp) :: value
      real(dp), allocatable :: eigenvalues(:), eigenvectors(:,:)
      integer :: local_status

      call symmetric_eigen_jacobi(a, eigenvalues, eigenvectors, local_status)
      if (size(eigenvalues) > 0) then
         value = minval(eigenvalues)
      else
         value = ieee_value(0.0_dp, ieee_quiet_nan)
      end if
      if (present(status)) status = local_status
   end function minimum_eigenvalue

   subroutine project_simplex(v, w, status)
      real(dp), intent(in) :: v(:)
      real(dp), intent(out) :: w(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: sorted(:)
      real(dp) :: cumulative, theta, key, total
      integer :: n, i, j, rho, alloc_stat

      n = size(v)
      if (n < 1 .or. size(w) /= n) then
         if (size(w) > 0) w = 0.0_dp
         if (present(status)) status = fincov_size_mismatch
         return
      end if
      allocate(sorted(n), stat=alloc_stat)
      if (alloc_stat /= 0) then
         w = 0.0_dp
         if (present(status)) status = fincov_allocation_failure
         return
      end if
      sorted = v
      do i = 2, n
         key = sorted(i)
         j = i - 1
         do while (j >= 1)
            if (sorted(j) >= key) exit
            sorted(j+1) = sorted(j)
            j = j - 1
         end do
         sorted(j+1) = key
      end do

      cumulative = 0.0_dp
      rho = 0
      do i = 1, n
         cumulative = cumulative + sorted(i)
         if (sorted(i) + (1.0_dp - cumulative) / real(i,dp) > 0.0_dp) rho = i
      end do
      if (rho == 0) then
         w = 1.0_dp / real(n,dp)
      else
         theta = (sum(sorted(1:rho)) - 1.0_dp) / real(rho,dp)
         w = max(v - theta, 0.0_dp)
         total = sum(w)
         if (total > 0.0_dp) w = w / total
      end if
      if (present(status)) status = fincov_ok
   end subroutine project_simplex

   pure logical function matrix_is_symmetric(a, tolerance)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(in), optional :: tolerance
      real(dp) :: tol

      if (size(a,1) /= size(a,2)) then
         matrix_is_symmetric = .false.
         return
      end if
      tol = 1000.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(a)))
      if (present(tolerance)) tol = tolerance
      matrix_is_symmetric = maxval(abs(a - transpose(a))) <= tol
   end function matrix_is_symmetric

   pure subroutine add_diagonal(a, diagonal)
      real(dp), intent(inout) :: a(:,:)
      real(dp), intent(in) :: diagonal(:)
      integer :: i, n

      n = min(size(diagonal), min(size(a,1), size(a,2)))
      do i = 1, n
         a(i,i) = a(i,i) + diagonal(i)
      end do
   end subroutine add_diagonal
end module fincov_linalg
