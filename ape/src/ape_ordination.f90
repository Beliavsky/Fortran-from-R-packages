! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
!
! Principal-coordinate analysis translated from ape R/pcoa.R.
! Upstream copyright and provenance are documented in NOTICE.md.
module ape_ordination
   use r_kinds, only : dp
   use r_linalg, only : symmetric_eigen, general_real_eigenvalues
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
   private

   integer, parameter, public :: pcoa_none = 0
   integer, parameter, public :: pcoa_lingoes = 1
   integer, parameter, public :: pcoa_cailliez = 2

   type, public :: pcoa_result
      integer :: correction = pcoa_none
      integer :: positive_rank = 0
      real(dp) :: trace = 0.0_dp
      real(dp) :: corrected_trace = 0.0_dp
      real(dp) :: correction_constant = 0.0_dp
      real(dp), allocatable :: eigenvalues(:)
      real(dp), allocatable :: relative_eigenvalues(:)
      real(dp), allocatable :: corrected_eigenvalues(:)
      real(dp), allocatable :: corrected_relative_eigenvalues(:)
      real(dp), allocatable :: cumulative_eigenvalues(:)
      real(dp), allocatable :: broken_stick(:)
      real(dp), allocatable :: cumulative_broken_stick(:)
      real(dp), allocatable :: vectors(:, :)
      real(dp), allocatable :: corrected_vectors(:, :)
   end type pcoa_result

   public :: pcoa

contains

   subroutine pcoa(distance, result, info, correction)
      !! Performs ape-style principal-coordinate analysis of a square distance matrix.
      real(dp), intent(in) :: distance(:, :) !! Symmetric distance matrix with zero diagonal and shape `(n,n)`.
      type(pcoa_result), intent(out) :: result !! Eigenvalue diagnostics, coordinates, and optional corrected coordinates.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid input or a failed eigendecomposition.
      character(len=*), intent(in), optional :: correction !! Negative-eigenvalue correction: `none`, `lingoes`, or `cailliez`.
      real(dp), allocatable :: centered(:, :)
      real(dp), allocatable :: corrected_kernel(:, :)
      real(dp), allocatable :: eig(:)
      real(dp), allocatable :: eig_cor(:)
      real(dp), allocatable :: rel_cor(:)
      real(dp), allocatable :: vectors(:, :)
      real(dp), allocatable :: vectors_cor(:, :)
      real(dp) :: eps
      real(dp) :: minimum_eigenvalue
      integer :: mode
      integer :: n
      integer :: status

      result = pcoa_result()
      info = 0
      n = size(distance, 1)
      if (n < 2 .or. size(distance, 2) /= n) then
         info = 1
         return
      end if
      if (.not. valid_distance_matrix(distance)) then
         info = 2
         return
      end if
      call parse_correction(correction, mode, status)
      if (status /= 0) then
         info = 3
         return
      end if

      eps = sqrt(epsilon(1.0_dp))
      centered = gower_center(-0.5_dp * distance * distance)
      result%trace = trace_matrix(centered)
      call symmetric_eigen(centered, eig, vectors, status, descending=.true.)
      if (status /= 0) then
         info = 10 + status
         return
      end if
      minimum_eigenvalue = minval(eig)
      where (abs(eig) < eps) eig = 0.0_dp
      result%eigenvalues = eig
      result%relative_eigenvalues = eig / result%trace
      result%positive_rank = count(eig > eps)
      call scaled_positive_vectors(eig, vectors, eps, result%vectors)

      if (minimum_eigenvalue > -eps) then
         result%correction = pcoa_none
         call positive_summary(eig, result%trace, eps, result%cumulative_eigenvalues, &
            result%broken_stick, result%cumulative_broken_stick)
         return
      end if

      call negative_eigen_summary(eig, minimum_eigenvalue, result%trace, eps, rel_cor)
      result%corrected_relative_eigenvalues = rel_cor
      if (mode == pcoa_none) then
         result%correction = pcoa_none
         result%cumulative_eigenvalues = cumulative_sum(rel_cor)
         call padded_broken_stick(count(rel_cor > eps), n, result%broken_stick)
         result%cumulative_broken_stick = cumulative_sum(result%broken_stick)
         return
      end if

      result%correction = mode
      select case (mode)
      case (pcoa_lingoes)
         result%correction_constant = -minimum_eigenvalue
         corrected_kernel = -0.5_dp * (distance * distance + 2.0_dp * result%correction_constant)
         call zero_diagonal(corrected_kernel)
      case (pcoa_cailliez)
         call cailliez_constant(distance, centered, result%correction_constant, status)
         if (status /= 0) then
            info = 20 + status
            return
         end if
         corrected_kernel = -0.5_dp * (distance + result%correction_constant)**2
         call zero_diagonal(corrected_kernel)
      end select

      corrected_kernel = gower_center(corrected_kernel)
      result%corrected_trace = trace_matrix(corrected_kernel)
      call symmetric_eigen(corrected_kernel, eig_cor, vectors_cor, status, descending=.true.)
      if (status /= 0) then
         info = 30 + status
         return
      end if
      where (abs(eig_cor) < eps) eig_cor = 0.0_dp
      if (minval(eig_cor) <= -eps) then
         info = 40
         return
      end if
      result%corrected_eigenvalues = eig_cor
      result%positive_rank = count(eig_cor > eps)
      result%corrected_relative_eigenvalues = eig_cor / result%corrected_trace
      result%cumulative_eigenvalues = cumulative_sum(result%corrected_relative_eigenvalues)
      call padded_broken_stick(result%positive_rank, n, result%broken_stick)
      result%cumulative_broken_stick = cumulative_sum(result%broken_stick)
      call scaled_positive_vectors(eig_cor, vectors_cor, eps, result%corrected_vectors)
   end subroutine pcoa

   pure logical function valid_distance_matrix(distance) result(valid)
      !! Checks finite nonnegative symmetry and a numerically zero diagonal.
      real(dp), intent(in) :: distance(:, :) !! Candidate square distance matrix.
      real(dp) :: scale
      real(dp) :: tolerance
      integer :: i
      integer :: j
      integer :: n

      n = size(distance, 1)
      valid = size(distance, 2) == n
      if (.not. valid) return
      if (.not. all(ieee_is_finite(distance))) then
         valid = .false.
         return
      end if
      scale = max(1.0_dp, maxval(abs(distance)))
      tolerance = 100.0_dp * epsilon(1.0_dp) * scale
      do i = 1, n
         if (abs(distance(i, i)) > tolerance) then
            valid = .false.
            return
         end if
         do j = 1, n
            if (distance(i, j) < -tolerance) then
               valid = .false.
               return
            end if
            if (abs(distance(i, j) - distance(j, i)) > tolerance) then
               valid = .false.
               return
            end if
         end do
      end do
   end function valid_distance_matrix

   pure function gower_center(matrix) result(centered)
      !! Double-centers a square matrix with `(I-11'/n) A (I-11'/n)`.
      real(dp), intent(in) :: matrix(:, :) !! Square matrix to double-center.
      real(dp) :: centered(size(matrix, 1), size(matrix, 2))
      real(dp) :: grand_mean
      real(dp), allocatable :: row_mean(:)
      integer :: i
      integer :: j
      integer :: n

      n = size(matrix, 1)
      allocate(row_mean(n))
      do i = 1, n
         row_mean(i) = sum(matrix(i, :)) / real(n, dp)
      end do
      grand_mean = sum(row_mean) / real(n, dp)
      do j = 1, n
         do i = 1, n
            centered(i, j) = matrix(i, j) - row_mean(i) - row_mean(j) + grand_mean
         end do
      end do
      centered = 0.5_dp * (centered + transpose(centered))
   end function gower_center

   pure real(dp) function trace_matrix(matrix) result(value)
      !! Returns the trace of a square matrix.
      real(dp), intent(in) :: matrix(:, :) !! Square matrix whose diagonal is summed.
      integer :: i

      value = 0.0_dp
      do i = 1, size(matrix, 1)
         value = value + matrix(i, i)
      end do
   end function trace_matrix

   subroutine cailliez_constant(distance, delta1, constant, info)
      !! Computes ape's Cailliez additive constant from the Gower-Legendre block matrix.
      real(dp), intent(in) :: distance(:, :) !! Original symmetric distance matrix with shape `(n,n)`.
      real(dp), intent(in) :: delta1(:, :) !! Gower-centered `-0.5*D**2` matrix with shape `(n,n)`.
      real(dp), intent(out) :: constant !! Maximum real part of the block-matrix eigenvalues.
      integer, intent(out) :: info !! Zero on success or a linear-algebra status code.
      real(dp), allocatable :: delta2(:, :)
      real(dp), allocatable :: imaginary(:)
      real(dp), allocatable :: real_part(:)
      real(dp), allocatable :: special(:, :)
      integer :: i
      integer :: n

      n = size(distance, 1)
      allocate(delta2(n, n), special(2 * n, 2 * n))
      delta2 = gower_center(-0.5_dp * distance)
      special = 0.0_dp
      special(1:n, n + 1:2 * n) = 2.0_dp * delta1
      do i = 1, n
         special(n + i, i) = -1.0_dp
      end do
      special(n + 1:2 * n, n + 1:2 * n) = -4.0_dp * delta2
      call general_real_eigenvalues(special, real_part, imaginary, info)
      if (info /= 0) then
         constant = 0.0_dp
         return
      end if
      constant = maxval(real_part)
   end subroutine cailliez_constant

   pure subroutine zero_diagonal(matrix)
      !! Replaces the diagonal of a square matrix by exact zeros.
      real(dp), intent(inout) :: matrix(:, :) !! Square matrix whose diagonal is modified in place.
      integer :: i

      do i = 1, size(matrix, 1)
         matrix(i, i) = 0.0_dp
      end do
   end subroutine zero_diagonal

   subroutine parse_correction(correction, mode, info)
      !! Parses the user-facing PCoA correction name.
      character(len=*), intent(in), optional :: correction !! Optional correction name.
      integer, intent(out) :: mode !! One of `pcoa_none`, `pcoa_lingoes`, or `pcoa_cailliez`.
      integer, intent(out) :: info !! Zero for a recognized correction and one otherwise.
      character(len=:), allocatable :: name

      mode = pcoa_none
      info = 0
      if (.not. present(correction)) return
      name = lower_ascii(trim(adjustl(correction)))
      if (len(name) == 0) then
         info = 1
      else if (index('none', name) == 1) then
         mode = pcoa_none
      else if (index('lingoes', name) == 1) then
         mode = pcoa_lingoes
      else if (index('cailliez', name) == 1) then
         mode = pcoa_cailliez
      else
         info = 1
      end if
   end subroutine parse_correction

   pure function lower_ascii(text) result(lower)
      !! Converts ASCII uppercase letters to lowercase without locale dependence.
      character(len=*), intent(in) :: text !! Input ASCII text.
      character(len=len(text)) :: lower
      integer :: code
      integer :: i

      lower = text
      do i = 1, len(text)
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code + 32)
      end do
   end function lower_ascii

   subroutine scaled_positive_vectors(values, vectors, eps, coordinates)
      !! Scales eigenvectors by square roots of positive eigenvalues as PCoA coordinates.
      real(dp), intent(in) :: values(:) !! Eigenvalues in descending order.
      real(dp), intent(in) :: vectors(:, :) !! Corresponding eigenvectors by column.
      real(dp), intent(in) :: eps !! Positive-eigenvalue threshold.
      real(dp), allocatable, intent(out) :: coordinates(:, :) !! Principal coordinates for positive axes only.
      integer :: j
      integer :: k

      k = count(values > eps)
      allocate(coordinates(size(vectors, 1), k))
      do j = 1, k
         coordinates(:, j) = vectors(:, j) * sqrt(values(j))
      end do
   end subroutine scaled_positive_vectors

   pure subroutine positive_summary(values, trace, eps, cumulative, broken, cumulative_broken)
      !! Builds the positive-axis relative-eigenvalue and broken-stick summaries.
      real(dp), intent(in) :: values(:) !! Eigenvalues in descending order.
      real(dp), intent(in) :: trace !! Sum of all eigenvalues before zero-thresholding.
      real(dp), intent(in) :: eps !! Positive-eigenvalue threshold.
      real(dp), allocatable, intent(out) :: cumulative(:) !! Cumulative relative eigenvalues for positive axes.
      real(dp), allocatable, intent(out) :: broken(:) !! Broken-stick expectations for positive axes.
      real(dp), allocatable, intent(out) :: cumulative_broken(:) !! Cumulative broken-stick expectations.
      real(dp), allocatable :: relative(:)
      integer :: k

      k = count(values > eps)
      allocate(relative(k))
      if (k > 0) relative = values(1:k) / trace
      cumulative = cumulative_sum(relative)
      call padded_broken_stick(k, k, broken)
      cumulative_broken = cumulative_sum(broken)
   end subroutine positive_summary

   pure subroutine negative_eigen_summary(values, minimum_value, trace, eps, relative_corrected)
      !! Reproduces ape's corrected relative eigenvalues when negative axes are retained diagnostically.
      real(dp), intent(in) :: values(:) !! Thresholded eigenvalues in descending order.
      real(dp), intent(in) :: minimum_value !! Most negative eigenvalue before zero-thresholding.
      real(dp), intent(in) :: trace !! Trace of the original Gower-centered matrix.
      real(dp), intent(in) :: eps !! Zero-eigenvalue threshold.
      real(dp), allocatable, intent(out) :: relative_corrected(:) !! ape Eq. 9.27 values, including zero-axis relocation.
      real(dp), allocatable :: temporary(:)
      integer :: first_zero
      integer :: i
      integer :: k
      integer :: n

      n = size(values)
      allocate(temporary(n), relative_corrected(n))
      temporary = (values - minimum_value) / (trace - real(n - 1, dp) * minimum_value)
      first_zero = 0
      do i = 1, n
         if (abs(values(i)) < eps) then
            first_zero = i
            exit
         end if
      end do
      if (first_zero == 0) then
         relative_corrected = temporary
         return
      end if
      k = 0
      do i = 1, n
         if (i == first_zero) cycle
         k = k + 1
         relative_corrected(k) = temporary(i)
      end do
      relative_corrected(n) = 0.0_dp
   end subroutine negative_eigen_summary

   pure subroutine padded_broken_stick(k, n, values)
      !! Returns broken-stick expectations for `k` axes padded with zeros to length `n`.
      integer, intent(in) :: k !! Number of positive axes receiving broken-stick mass.
      integer, intent(in) :: n !! Output length, which must be at least `k`.
      real(dp), allocatable, intent(out) :: values(:) !! Broken-stick expectations followed by zero padding.
      integer :: i
      integer :: j

      allocate(values(n))
      values = 0.0_dp
      if (k <= 0) return
      do i = 1, k
         do j = i, k
            values(i) = values(i) + 1.0_dp / real(j, dp)
         end do
         values(i) = values(i) / real(k, dp)
      end do
   end subroutine padded_broken_stick

   pure function cumulative_sum(values) result(cumulative)
      !! Returns the ordinary cumulative sum of a real vector.
      real(dp), intent(in) :: values(:) !! Input vector.
      real(dp) :: cumulative(size(values))
      integer :: i

      if (size(values) == 0) return
      cumulative(1) = values(1)
      do i = 2, size(values)
         cumulative(i) = cumulative(i - 1) + values(i)
      end do
   end function cumulative_sum

end module ape_ordination
