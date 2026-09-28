! SPDX-License-Identifier: GPL-3.0-only
module mbend_api
   use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
   use mbend_kinds, only : dp
   use r_linalg, only : symmetric_eigen, symmetric_eigenvalues
   implicit none
   private

   integer, parameter, public :: mbend_success = 0
   integer, parameter, public :: mbend_err_shape = 1
   integer, parameter, public :: mbend_err_asymmetric = 2
   integer, parameter, public :: mbend_err_weight_shape = 3
   integer, parameter, public :: mbend_err_weight_asymmetric = 4
   integer, parameter, public :: mbend_err_negative_weight = 5
   integer, parameter, public :: mbend_err_method = 6
   integer, parameter, public :: mbend_err_small_positive = 7
   integer, parameter, public :: mbend_err_already_pd = 8
   integer, parameter, public :: mbend_err_zero_weights = 9
   integer, parameter, public :: mbend_err_eigen = 10
   integer, parameter, public :: mbend_err_lrs_spectrum = 11

   type, public :: bend_result
      real(dp), allocatable :: bent(:, :)
      real(dp), allocatable :: init_ev(:)
      real(dp), allocatable :: final_ev(:)
      real(dp) :: min_dev = 0.0_dp
      real(dp) :: max_dev = 0.0_dp
      integer :: loc_min_dev(2) = 0
      integer :: loc_max_dev(2) = 0
      real(dp) :: ave_dev = 0.0_dp
      real(dp) :: aad = 0.0_dp
      real(dp) :: cor = 0.0_dp
      real(dp) :: rmsd = 0.0_dp
      logical :: has_weighted_stats = .false.
      integer :: w_gt_0 = 0
      real(dp) :: waad = 0.0_dp
      real(dp) :: wcor = 0.0_dp
      real(dp) :: wrmsd = 0.0_dp
      logical :: is_correlation = .false.
      logical :: used_weights = .false.
      logical :: converged = .false.
      integer :: iterations = 0
   end type bend_result

   public :: bend
   public :: dp

contains

   subroutine bend(inmat, result, info, wtmat, reciprocal, max_iter, small_positive, method)
      !! Bend a symmetric matrix to positive-definite form using the mbend algorithms.
      real(dp), intent(in) :: inmat(:, :) !! Symmetric matrix to bend; must be square and contain a negative eigenvalue.
      type(bend_result), intent(out) :: result !! Bent matrix, spectra, diagnostics, convergence state, and weighted statistics.
      integer, intent(out) :: info !! Zero on success; otherwise one of the public `mbend_err_*` status codes.
      real(dp), intent(in), optional :: wtmat(:, :) !! Optional nonnegative symmetric weight matrix with the same shape as `inmat`.
      logical, intent(in), optional :: reciprocal !! If true, replace each nonzero weight by its reciprocal before scaling.
      integer, intent(in), optional :: max_iter !! Maximum bending iterations; defaults to 10000.
      real(dp), intent(in), optional :: small_positive !! HJ eigenvalue floor; default 1e-4 and valid on (0, 0.1].
      character(len=*), intent(in), optional :: method !! Bending method: `hj` for Jorjani et al. or `lrs` for Schaeffer.

      character(len=3) :: method_name
      integer :: eigen_info, i, iter, j, m, max_iterations, n
      logical :: use_reciprocal
      real(dp) :: max_weight, s_value, small_value, w_value
      real(dp), allocatable :: bent(:, :), eigenval(:), eigenvec(:, :), udtu(:, :), weights(:, :)

      call reset_result(result)
      info = mbend_success
      n = size(inmat, 1)

      if (n < 1 .or. size(inmat, 2) /= n) then
         info = mbend_err_shape
         return
      end if
      if (.not. is_symmetric(inmat)) then
         info = mbend_err_asymmetric
         return
      end if

      method_name = 'hj '
      if (present(method)) method_name = trim(method)
      if (trim(method_name) /= 'hj' .and. trim(method_name) /= 'lrs') then
         info = mbend_err_method
         return
      end if

      max_iterations = 10000
      if (present(max_iter)) max_iterations = max_iter
      small_value = 0.0001_dp
      if (present(small_positive)) small_value = small_positive
      use_reciprocal = .false.
      if (present(reciprocal)) use_reciprocal = reciprocal

      allocate(bent(n, n), weights(n, n))
      bent = inmat
      if (present(wtmat)) then
         result%used_weights = .true.
         if (size(wtmat, 1) /= n .or. size(wtmat, 2) /= n) then
            info = mbend_err_weight_shape
            return
         end if
         if (.not. is_symmetric(wtmat)) then
            info = mbend_err_weight_asymmetric
            return
         end if
         if (any(wtmat < 0.0_dp)) then
            info = mbend_err_negative_weight
            return
         end if
         weights = wtmat
      else
         weights = 1.0_dp
      end if

      call symmetric_eigen(bent, eigenval, eigenvec, eigen_info, descending=.true.)
      if (eigen_info /= 0) then
         info = mbend_err_eigen
         return
      end if
      if (minval(eigenval) >= 0.0_dp) then
         info = mbend_err_already_pd
         return
      end if
      result%init_ev = eigenval

      if (use_reciprocal) then
         do j = 1, n
            do i = 1, n
               if (abs(weights(i, j)) > 0.0_dp) weights(i, j) = 1.0_dp / weights(i, j)
            end do
         end do
      end if

      result%is_correlation = has_unit_diagonal(inmat)
      if (result%is_correlation) then
         do i = 1, n
            weights(i, i) = 0.0_dp
         end do
      end if

      max_weight = maxval(weights)
      if (max_weight <= 0.0_dp) then
         info = mbend_err_zero_weights
         return
      end if
      weights = weights / max_weight

      if (trim(method_name) == 'hj') then
         if (small_value <= 0.0_dp .or. small_value > 0.1_dp) then
            info = mbend_err_small_positive
            return
         end if
      end if

      m = count(eigenval < 0.0_dp)
      iter = 0
      allocate(udtu(n, n))
      do while (iter < max_iterations .and. m > 0)
         if (trim(method_name) == 'hj') then
            where (eigenval < small_value) eigenval = small_value
         else
            if (n - m < 1) then
               info = mbend_err_lrs_spectrum
               return
            end if
            s_value = 2.0_dp * sum(pack(eigenval, eigenval < 0.0_dp))
            w_value = 100.0_dp * s_value * s_value + 1.0_dp
            small_value = eigenval(n - m)
            do j = n - m + 1, n
               eigenval(j) = small_value * (s_value - eigenval(j))**2 / w_value
            end do
            call rebuild_from_eigen(eigenvec, eigenval, udtu)
            call symmetric_eigenvalues(udtu, eigenval, eigen_info, descending=.true.)
            if (eigen_info /= 0) then
               info = mbend_err_eigen
               return
            end if
            where (eigenval < 0.0_dp) eigenval = small_value / 10.0_dp
         end if

         call rebuild_from_eigen(eigenvec, eigenval, udtu)
         bent = bent - (bent - udtu) * weights
         call symmetric_eigen(bent, eigenval, eigenvec, eigen_info, descending=.true.)
         if (eigen_info /= 0) then
            info = mbend_err_eigen
            return
         end if
         m = count(eigenval < 0.0_dp)
         iter = iter + 1
      end do

      result%bent = bent
      result%final_ev = eigenval
      result%iterations = iter
      result%converged = m == 0
      call compute_diagnostics(inmat, bent, weights, result)
   end subroutine bend

   pure subroutine rebuild_from_eigen(vectors, values, matrix)
      !! Rebuild `V diag(values) V^T` for orthonormal eigenvectors stored by column.
      real(dp), intent(in) :: vectors(:, :) !! Orthonormal eigenvector matrix with one eigenvector per column.
      real(dp), intent(in) :: values(:) !! Eigenvalues corresponding to the columns of `vectors`.
      real(dp), intent(out) :: matrix(:, :) !! Reconstructed symmetric matrix with shape `(n, n)`.

      integer :: j

      matrix = 0.0_dp
      do j = 1, size(values)
         matrix = matrix + values(j) * spread(vectors(:, j), 2, size(values)) * &
                  spread(vectors(:, j), 1, size(values))
      end do
   end subroutine rebuild_from_eigen

   pure logical function is_symmetric(matrix)
      !! Test symmetry using a tolerance similar to R's matrix symmetry check.
      real(dp), intent(in) :: matrix(:, :) !! Matrix whose symmetry is tested.

      real(dp) :: scale, tolerance

      if (size(matrix, 1) /= size(matrix, 2)) then
         is_symmetric = .false.
         return
      end if
      scale = max(1.0_dp, maxval(abs(matrix)))
      tolerance = 100.0_dp * epsilon(1.0_dp) * scale
      is_symmetric = maxval(abs(matrix - transpose(matrix))) <= tolerance
   end function is_symmetric

   pure logical function has_unit_diagonal(matrix)
      !! Test whether every diagonal entry is exactly one, matching mbend correlation detection.
      real(dp), intent(in) :: matrix(:, :) !! Square matrix whose diagonal is inspected.

      integer :: i

      has_unit_diagonal = .true.
      do i = 1, size(matrix, 1)
         if (abs(matrix(i, i) - 1.0_dp) > 0.0_dp) then
            has_unit_diagonal = .false.
            return
         end if
      end do
   end function has_unit_diagonal

   pure subroutine reset_result(result)
      !! Initialize scalar metadata in a result before a bending attempt.
      type(bend_result), intent(out) :: result !! Result object reset to deterministic defaults.

      real(dp) :: nan_value

      nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
      result%min_dev = nan_value
      result%max_dev = nan_value
      result%loc_min_dev = 0
      result%loc_max_dev = 0
      result%ave_dev = nan_value
      result%aad = nan_value
      result%cor = nan_value
      result%rmsd = nan_value
      result%has_weighted_stats = .false.
      result%w_gt_0 = 0
      result%waad = nan_value
      result%wcor = nan_value
      result%wrmsd = nan_value
      result%is_correlation = .false.
      result%used_weights = .false.
      result%converged = .false.
      result%iterations = 0
   end subroutine reset_result

   pure subroutine compute_diagnostics(inmat, bent, weights, result)
      !! Compute the deviation, ordinary-correlation, and optional weighted diagnostics returned by R `bend`.
      real(dp), intent(in) :: inmat(:, :) !! Original symmetric matrix.
      real(dp), intent(in) :: bent(:, :) !! Bent symmetric matrix with the same shape as `inmat`.
      real(dp), intent(in) :: weights(:, :) !! Final normalized elementwise bending weights.
      type(bend_result), intent(inout) :: result !! Result object receiving diagnostic fields.

      integer :: i, j, k, n, n_stats
      real(dp), allocatable :: dev(:), full_dev(:, :), w(:), x(:), y(:)

      n = size(inmat, 1)
      allocate(full_dev(n, n))
      full_dev = bent - inmat
      call first_extreme_location(full_dev, .true., result%loc_min_dev)
      call first_extreme_location(full_dev, .false., result%loc_max_dev)

      if (result%is_correlation) then
         n_stats = n * (n - 1) / 2
      else
         n_stats = n * (n + 1) / 2
      end if
      allocate(dev(n_stats), w(n_stats), x(n_stats), y(n_stats))

      k = 0
      do j = 1, n
         if (result%is_correlation) then
            do i = 1, j - 1
               k = k + 1
               dev(k) = full_dev(i, j)
               w(k) = weights(i, j)
               x(k) = inmat(i, j)
               y(k) = bent(i, j)
            end do
         else
            do i = 1, j
               k = k + 1
               dev(k) = full_dev(i, j)
               w(k) = weights(i, j)
               x(k) = inmat(i, j)
               y(k) = bent(i, j)
            end do
         end if
      end do

      result%min_dev = minval(dev)
      result%max_dev = maxval(dev)
      result%ave_dev = sum(dev) / real(n_stats, dp)
      result%aad = sum(abs(dev)) / real(n_stats, dp)
      result%cor = pearson_correlation(x, y)
      result%rmsd = sqrt(sum(dev * dev) / real(n_stats, dp))

      result%has_weighted_stats = any(.not. zero_or_one(w))
      if (result%has_weighted_stats) call compute_weighted_diagnostics(dev, w, x, y, result)
   end subroutine compute_diagnostics

   pure subroutine first_extreme_location(matrix, minimum, location)
      !! Return the first row/column location of a global extreme in Fortran/R column-major order.
      real(dp), intent(in) :: matrix(:, :) !! Matrix searched for its global minimum or maximum.
      logical, intent(in) :: minimum !! If true search for the minimum; otherwise search for the maximum.
      integer, intent(out) :: location(2) !! One-based `(row, column)` location of the first matching extreme.

      integer :: i, j
      real(dp) :: extreme

      if (minimum) then
         extreme = minval(matrix)
      else
         extreme = maxval(matrix)
      end if
      location = 0
      do j = 1, size(matrix, 2)
         do i = 1, size(matrix, 1)
            if (abs(matrix(i, j) - extreme) <= 0.0_dp) then
               location = [i, j]
               return
            end if
         end do
      end do
   end subroutine first_extreme_location

   pure elemental logical function zero_or_one(value)
      !! Test exact membership in the set `{0, 1}` without real equality comparisons.
      real(dp), intent(in) :: value !! Normalized weight to classify.

      zero_or_one = abs(value) <= 0.0_dp .or. abs(value - 1.0_dp) <= 0.0_dp
   end function zero_or_one

   pure real(dp) function pearson_correlation(x, y)
      !! Pearson correlation of two vectors, returning NaN when correlation is undefined.
      real(dp), intent(in) :: x(:) !! First observation vector.
      real(dp), intent(in) :: y(:) !! Second observation vector with the same length as `x`.

      real(dp) :: denom, mean_x, mean_y
      real(dp), allocatable :: centered_x(:), centered_y(:)

      if (size(x) /= size(y) .or. size(x) < 2) then
         pearson_correlation = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      mean_x = sum(x) / real(size(x), dp)
      mean_y = sum(y) / real(size(y), dp)
      centered_x = x - mean_x
      centered_y = y - mean_y
      denom = sqrt(sum(centered_x * centered_x) * sum(centered_y * centered_y))
      if (denom <= tiny(1.0_dp)) then
         pearson_correlation = ieee_value(0.0_dp, ieee_quiet_nan)
      else
         pearson_correlation = sum(centered_x * centered_y) / denom
      end if
   end function pearson_correlation

   pure subroutine compute_weighted_diagnostics(dev, w, x, y, result)
      !! Compute mbend's weighted AAD, correlation, and RMSD using positive weights only.
      real(dp), intent(in) :: dev(:) !! Deviations corresponding to upper-triangle matrix elements.
      real(dp), intent(in) :: w(:) !! Normalized bending weights corresponding to `dev`.
      real(dp), intent(in) :: x(:) !! Original matrix elements corresponding to `dev`.
      real(dp), intent(in) :: y(:) !! Bent matrix elements corresponding to `dev`.
      type(bend_result), intent(inout) :: result !! Result object receiving weighted diagnostic fields.

      integer :: i, k, n_positive
      real(dp), allocatable :: dpos(:), q(:), wpos(:), xpos(:), ypos(:)

      n_positive = count(w > 0.0_dp)
      result%w_gt_0 = n_positive
      if (n_positive == 0) return
      allocate(dpos(n_positive), wpos(n_positive), xpos(n_positive), ypos(n_positive), q(n_positive))
      k = 0
      do i = 1, size(w)
         if (w(i) > 0.0_dp) then
            k = k + 1
            dpos(k) = dev(i)
            wpos(k) = w(i)
            xpos(k) = x(i)
            ypos(k) = y(i)
         end if
      end do
      q = 1.0_dp / wpos
      result%waad = sum(abs(dpos / wpos)) / sum(q)
      result%wcor = weighted_correlation(xpos, ypos, q)
      result%wrmsd = sqrt(sum((dpos / wpos)**2) / sum(q * q))
   end subroutine compute_weighted_diagnostics

   pure real(dp) function weighted_correlation(x, y, weight)
      !! Weighted Pearson correlation; normalization of the positive weights is immaterial.
      real(dp), intent(in) :: x(:) !! First observation vector.
      real(dp), intent(in) :: y(:) !! Second observation vector with the same length as `x`.
      real(dp), intent(in) :: weight(:) !! Positive observation weights with the same length as `x`.

      real(dp) :: denom, mean_x, mean_y, sum_weight
      real(dp), allocatable :: centered_x(:), centered_y(:)

      if (size(x) /= size(y) .or. size(x) /= size(weight) .or. size(x) < 2) then
         weighted_correlation = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      sum_weight = sum(weight)
      if (sum_weight <= 0.0_dp) then
         weighted_correlation = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      mean_x = sum(weight * x) / sum_weight
      mean_y = sum(weight * y) / sum_weight
      centered_x = x - mean_x
      centered_y = y - mean_y
      denom = sqrt(sum(weight * centered_x * centered_x) * sum(weight * centered_y * centered_y))
      if (denom <= tiny(1.0_dp)) then
         weighted_correlation = ieee_value(0.0_dp, ieee_quiet_nan)
      else
         weighted_correlation = sum(weight * centered_x * centered_y) / denom
      end if
   end function weighted_correlation

end module mbend_api
