! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
!
! Lynch comparative method translated from ape R/compar.lynch.R.
! Upstream copyright/provenance are documented in NOTICE.md.
module ape_compar_lynch
   use r_kinds, only : dp
   use r_linalg, only : inverse_matrix, spd_inverse_logdet
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
   private

   real(dp), parameter :: two_pi = 6.2831853071795864769252867665590058_dp

   type, public :: compar_lynch_result
      real(dp), allocatable :: environmental_covariance(:, :)
      real(dp), allocatable :: phylogenetic_covariance(:, :)
      real(dp), allocatable :: phylogenetic_effect(:, :)
      real(dp), allocatable :: environmental_effect(:, :)
      real(dp), allocatable :: mean(:)
      real(dp) :: log_likelihood = -huge(1.0_dp)
      integer :: iterations = 0
      logical :: converged = .false.
   end type compar_lynch_result

   public :: compar_lynch_fit

contains

   subroutine compar_lynch_fit(x, g, result, info, eps, initial_fraction, max_iter)
      !! Fits Lynch's iterative phylogenetic/environmental variance decomposition.
      real(dp), intent(in) :: x(:, :) !! Trait matrix with species in rows and continuous characters in columns.
      real(dp), intent(in) :: g(:, :) !! Positive-definite phylogenetic covariance matrix with one row/column per species.
      type(compar_lynch_result), intent(out) :: result !! Covariances, fitted effects, means, likelihood, and convergence state.
      integer, intent(out) :: info !! Zero on convergence or nonzero for invalid input, singular matrices, or iteration limit.
      real(dp), intent(in), optional :: eps !! Elementwise effect-convergence tolerance; default `1e-4` as in ape.
      real(dp), intent(in), optional :: initial_fraction !! Deterministic initial phylogenetic variance fraction; default `0.5`.
      integer, intent(in), optional :: max_iter !! Iteration safety limit; default 1000 because upstream has no explicit cap.
      real(dp), allocatable :: a0(:, :)
      real(dp), allocatable :: a1(:, :)
      real(dp), allocatable :: d_inv(:, :)
      real(dp), allocatable :: e0(:, :)
      real(dp), allocatable :: e1(:, :)
      real(dp), allocatable :: g_inv(:, :)
      real(dp), allocatable :: info_matrix(:, :)
      real(dp), allocatable :: inverse_v(:, :)
      real(dp), allocatable :: r_inv(:, :)
      real(dp), allocatable :: residual(:)
      real(dp), allocatable :: trvare(:, :)
      real(dp), allocatable :: u0(:)
      real(dp), allocatable :: v(:, :)
      real(dp), allocatable :: vara(:, :)
      real(dp), allocatable :: vare(:, :)
      real(dp), allocatable :: vcvz(:, :)
      real(dp), allocatable :: work(:, :)
      real(dp), allocatable :: z(:)
      real(dp) :: fraction
      real(dp) :: logdet
      real(dp) :: quadratic
      real(dp) :: tolerance
      integer :: bi
      integer :: bj
      integer :: i
      integer :: iter
      integer :: iter_limit
      integer :: j
      integer :: k
      integer :: n
      integer :: nk
      integer :: status

      result = compar_lynch_result()
      info = 0
      n = size(x, 1)
      k = size(x, 2)
      if (n < 2 .or. k < 1 .or. size(g, 1) /= n .or. size(g, 2) /= n) then
         info = 1
         return
      end if
      if (any(.not. ieee_is_finite(x)) .or. any(.not. ieee_is_finite(g))) then
         info = 2
         return
      end if
      tolerance = 1.0e-4_dp
      if (present(eps)) tolerance = eps
      fraction = 0.5_dp
      if (present(initial_fraction)) fraction = initial_fraction
      iter_limit = 1000
      if (present(max_iter)) iter_limit = max_iter
      if (tolerance <= 0.0_dp .or. fraction <= 0.0_dp .or. fraction >= 1.0_dp .or. iter_limit < 1) then
         info = 3
         return
      end if

      call inverse_matrix(g, g_inv, status)
      if (status /= 0) then
         info = 10 + status
         return
      end if
      nk = n * k
      allocate(vcvz(k, k), vara(k, k), vare(k, k), trvare(k, k), u0(k))
      call sample_covariance(x, vcvz)
      vara = fraction * vcvz
      vare = (1.0_dp - fraction) * vcvz
      u0 = column_means(x)
      allocate(a0(n, k), e0(n, k), a1(n, k), e1(n, k))
      a0 = 0.0_dp
      e0 = 0.0_dp
      a1 = 1.0_dp
      e1 = 1.0_dp
      z = reshape(x, [nk])

      do iter = 1, iter_limit
         result%iterations = iter
         a1 = a0
         e1 = e0
         call inverse_matrix(kronecker_cov(vare, identity_matrix(n)), r_inv, status)
         if (status /= 0) then
            info = 20 + status
            return
         end if
         call inverse_matrix(kronecker_cov(vara, g), d_inv, status)
         if (status /= 0) then
            info = 30 + status
            return
         end if
         call inverse_matrix(r_inv + d_inv, info_matrix, status)
         if (status /= 0) then
            info = 40 + status
            return
         end if
         residual = z - repeated_means(u0, n)
         a0 = reshape(matmul(info_matrix, matmul(r_inv, residual)), [n, k])
         e0 = reshape(residual, [n, k]) - a0

         do i = 1, k
            bi = (i - 1) * n
            do j = 1, k
               bj = (j - 1) * n
               trvare(i, j) = trace_matrix(info_matrix(bi + 1:bi + n, bj + 1:bj + n))
            end do
         end do
         call sample_covariance(e0, work)
         vare = (real(n - 1, dp) * work + trvare) / real(n, dp)
         do i = 1, k
            bi = (i - 1) * n
            do j = 1, k
               bj = (j - 1) * n
               vara(i, j) = (dot_product(a0(:, i), matmul(g_inv, a0(:, j))) + &
                  trace_matrix(matmul(g_inv, info_matrix(bi + 1:bi + n, bj + 1:bj + n)))) / real(n, dp)
            end do
         end do
         u0 = column_means(x - a0)
         if (max(maxval(abs(a0 - a1)), maxval(abs(e0 - e1))) <= tolerance) then
            result%converged = .true.
            exit
         end if
      end do
      if (.not. result%converged) then
         info = 50
         return
      end if

      v = kronecker_cov(vara, g) + kronecker_cov(vare, identity_matrix(n))
      call spd_inverse_logdet(v, inverse_v, logdet, status)
      if (status /= 0) then
         info = 60 + status
         return
      end if
      residual = z - repeated_means(u0, n)
      quadratic = dot_product(residual, matmul(inverse_v, residual))
      result%log_likelihood = -real(n, dp) * log(two_pi) - 0.5_dp * logdet - 0.5_dp * quadratic
      result%environmental_covariance = vare
      result%phylogenetic_covariance = vara
      result%phylogenetic_effect = a0
      result%environmental_effect = e0
      result%mean = u0
   end subroutine compar_lynch_fit

   pure subroutine sample_covariance(x, covariance)
      !! Computes the R `var` sample covariance of matrix columns using denominator `n-1`.
      real(dp), intent(in) :: x(:, :) !! Observation-by-variable matrix.
      real(dp), allocatable, intent(out) :: covariance(:, :) !! Symmetric sample covariance of columns.
      real(dp), allocatable :: centered(:, :)
      real(dp) :: mean(size(x, 2))
      integer :: n

      n = size(x, 1)
      mean = column_means(x)
      allocate(centered(size(x, 1), size(x, 2)))
      centered = x - spread(mean, 1, n)
      allocate(covariance(size(x, 2), size(x, 2)))
      covariance = matmul(transpose(centered), centered) / real(n - 1, dp)
   end subroutine sample_covariance

   pure function column_means(x) result(mean)
      !! Returns arithmetic means of matrix columns.
      real(dp), intent(in) :: x(:, :) !! Observation-by-variable matrix.
      real(dp) :: mean(size(x, 2))

      mean = sum(x, dim=1) / real(size(x, 1), dp)
   end function column_means

   pure function repeated_means(mean, n) result(vector)
      !! Repeats each trait mean `n` times to match R's column-major vectorization of the trait matrix.
      real(dp), intent(in) :: mean(:) !! One mean for each trait column.
      integer, intent(in) :: n !! Number of species rows per trait block.
      real(dp) :: vector(n * size(mean))
      integer :: j

      do j = 1, size(mean)
         vector((j - 1) * n + 1:j * n) = mean(j)
      end do
   end function repeated_means

   pure function identity_matrix(n) result(identity)
      !! Constructs an `n` by `n` identity matrix.
      integer, intent(in) :: n !! Requested positive matrix order.
      real(dp) :: identity(n, n)
      integer :: i

      identity = 0.0_dp
      do i = 1, n
         identity(i, i) = 1.0_dp
      end do
   end function identity_matrix

   pure function kronecker_cov(left, right) result(product)
      !! Forms the Kronecker product using the block ordering of R's `%x%` operator.
      real(dp), intent(in) :: left(:, :) !! Trait covariance matrix defining block multipliers.
      real(dp), intent(in) :: right(:, :) !! Species covariance matrix copied into each block.
      real(dp) :: product(size(left, 1) * size(right, 1), size(left, 2) * size(right, 2))
      integer :: i
      integer :: j
      integer :: nr
      integer :: nc

      nr = size(right, 1)
      nc = size(right, 2)
      do i = 1, size(left, 1)
         do j = 1, size(left, 2)
            product((i - 1) * nr + 1:i * nr, (j - 1) * nc + 1:j * nc) = left(i, j) * right
         end do
      end do
   end function kronecker_cov

   pure real(dp) function trace_matrix(matrix) result(trace)
      !! Returns the trace of a matrix's main diagonal.
      real(dp), intent(in) :: matrix(:, :) !! Matrix whose main diagonal is summed.
      integer :: i

      trace = 0.0_dp
      do i = 1, min(size(matrix, 1), size(matrix, 2))
         trace = trace + matrix(i, i)
      end do
   end function trace_matrix

end module ape_compar_lynch
