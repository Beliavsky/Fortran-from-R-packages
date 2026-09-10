! SPDX-License-Identifier: GPL-2.0-or-later
! Weighted and oblique SSA kernels translated from Rssa 1.1.
module rssa_oblique
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
   use rssa_kinds, only : dp
   use rssa_matrices, only : hankel_matrix
   use rssa_metrics, only : wcor_default
   use rssa_reconstruction, only : elementary_series_ssa
   use rssa_types, only : rssa_invalid_input, rssa_not_supported, rssa_numerical_failure, rssa_success, ssa_result
   use r_linalg, only : least_squares_svd, symmetric_eigen, thin_svd
   implicit none
   private
   public :: decompose_ossa, decompose_wossa, fossa, fossa_ssa, owcor_ssa, wcor_ossa

   interface fossa
      module procedure fossa_ssa
   end interface fossa
contains

   function fossa_ssa(object, indices, filter, gamma, normalize) result(out)
      type(ssa_result), intent(in) :: object !! Real 1-D SSA decomposition to refine by filter-adjusted O-SSA.
      integer, intent(in), optional :: indices(:) !! One-based eigentriples to rotate; defaults to all stored components.
      real(dp), intent(in), optional :: filter(:) !! Reversed finite impulse-response coefficients; default is [-1, 1].
      real(dp), intent(in), optional :: gamma !! Filter-adjustment weight; omitted or infinite means filter-only FOSSA.
      logical, intent(in), optional :: normalize !! Use normalized right vectors before filtering; default is true.
      type(ssa_result) :: out
      integer, allocatable :: idx(:)
      real(dp), allocatable :: base(:, :), covariance(:, :), eigenvalues(:), filtered(:, :)
      real(dp), allocatable :: impulse(:), new_left(:, :), new_right(:, :), rotation(:, :), scaled_left(:, :)
      real(dp), allocatable :: scaled_right(:, :), stacked(:, :)
      real(dp) :: adjustment
      integer :: i, local_info, n_components, n_rows
      logical :: filtered_only, use_normalized

      out = object
      out%info = rssa_invalid_input
      if (.not. allocated(object%sigma)) return
      if (.not. allocated(object%u) .or. .not. allocated(object%v)) return
      if (.not. allocated(object%series)) return
      if (size(object%u, 2) /= size(object%sigma)) return
      if (size(object%v, 2) /= size(object%sigma)) return
      if (object%window < 1 .or. object%window > size(object%series)) return

      call selected_indices(size(object%sigma), indices, idx, local_info)
      if (local_info /= rssa_success) return
      n_components = size(idx)
      if (n_components < 1) return

      if (present(filter)) then
         if (size(filter) < 1 .or. size(filter) > size(object%v, 1)) return
         allocate(impulse(size(filter)))
         impulse = filter
      else
         allocate(impulse(2))
         impulse = [-1.0_dp, 1.0_dp]
      end if

      use_normalized = .true.
      if (present(normalize)) use_normalized = normalize
      filtered_only = .true.
      adjustment = 0.0_dp
      if (present(gamma)) then
         if (ieee_is_nan(gamma)) return
         if (ieee_is_finite(gamma)) then
            filtered_only = .false.
            adjustment = gamma
         end if
      end if

      allocate(base(size(object%v, 1), n_components))
      if (use_normalized) then
         base = object%v(:, idx)
      else
         base = object%v(:, idx) * spread(object%sigma(idx), 1, size(object%v, 1))
      end if
      call filter_columns(base, impulse, filtered, local_info)
      if (local_info /= rssa_success) return

      if (filtered_only) then
         n_rows = size(filtered, 1)
         allocate(stacked(n_rows, n_components))
         stacked = filtered
      else
         n_rows = size(base, 1) + size(filtered, 1)
         allocate(stacked(n_rows, n_components))
         stacked(1:size(base, 1), :) = base
         stacked(size(base, 1) + 1:n_rows, :) = adjustment * filtered
      end if
      covariance = matmul(transpose(stacked), stacked)
      call symmetric_eigen(covariance, eigenvalues, rotation, local_info, descending=.true.)
      if (local_info /= 0 .or. size(rotation, 1) /= n_components) then
         out%info = rssa_numerical_failure
         return
      end if

      if (use_normalized) then
         allocate(scaled_left(size(object%u, 1), n_components))
         scaled_left = object%u(:, idx) * spread(object%sigma(idx), 1, size(object%u, 1))
         new_left = matmul(scaled_left, rotation)
         new_right = matmul(object%v(:, idx), rotation)
      else
         allocate(scaled_right(size(object%v, 1), n_components))
         scaled_right = object%v(:, idx) * spread(object%sigma(idx), 1, size(object%v, 1))
         new_left = matmul(object%u(:, idx), rotation)
         new_right = matmul(scaled_right, rotation)
      end if

      do i = 1, n_components
         out%u(:, idx(i)) = new_left(:, i)
         out%v(:, idx(i)) = new_right(:, i)
         out%sigma(idx(i)) = 1.0_dp
      end do
      if (minval(idx) <= object%n_special_right + object%n_special_left) then
         out%n_special_right = 0
         out%n_special_left = 0
      end if
      out%info = rssa_success
   end function fossa_ssa

   function decompose_ossa(object) result(out)
      type(ssa_result), intent(in) :: object !! Existing oblique SSA object for which continuation is requested.
      type(ssa_result) :: out

      out = object
      out%info = rssa_not_supported
   end function decompose_ossa

   function decompose_wossa(series, window, neig, column_weights, row_weights) result(out)
      real(dp), intent(in) :: series(:) !! Real source series used to form the Hankel trajectory matrix.
      integer, intent(in) :: window !! Window length L in 1..N.
      integer, intent(in), optional :: neig !! Number of leading weighted-oblique eigentriples to retain.
      real(dp), intent(in), optional :: column_weights(:) !! Nonnegative column-space metric weights of length L.
      real(dp), intent(in), optional :: row_weights(:) !! Nonnegative row-space metric weights of length N-L+1.
      type(ssa_result) :: out
      real(dp), allocatable :: col_inv(:), col_scale(:), h(:, :), row_inv(:), row_scale(:), singular_values(:)
      real(dp), allocatable :: u(:, :), vt(:, :), wh(:, :)
      real(dp) :: nu_col, nv_col
      integer :: i, info, k, rank_max

      if (window < 1 .or. window > size(series)) then
         out%info = rssa_invalid_input
         return
      end if
      k = size(series) - window + 1
      if (present(column_weights)) then
         if (size(column_weights) /= window .or. any(column_weights < 0.0_dp)) then
            out%info = rssa_invalid_input
            return
         end if
      end if
      if (present(row_weights)) then
         if (size(row_weights) /= k .or. any(row_weights < 0.0_dp)) then
            out%info = rssa_invalid_input
            return
         end if
      end if

      allocate(out%series(size(series)))
      out%series = series
      out%window = window
      h = hankel_matrix(series, window)
      allocate(out%trajectory(size(h, 1), size(h, 2)))
      out%trajectory = h
      allocate(col_scale(window), row_scale(k), col_inv(window), row_inv(k))
      col_scale = 1.0_dp
      row_scale = 1.0_dp
      if (present(column_weights)) col_scale = sqrt(column_weights)
      if (present(row_weights)) row_scale = sqrt(row_weights)
      call reciprocal_scale(col_scale, col_inv)
      call reciprocal_scale(row_scale, row_inv)

      allocate(wh(window, k))
      wh = h * spread(col_scale, 2, k) * spread(row_scale, 1, window)
      rank_max = min(window, k)
      i = min(50, rank_max)
      if (present(neig)) i = max(1, min(neig, rank_max))
      call thin_svd(wh, u, singular_values, vt, info)
      if (info /= 0) then
         out%info = rssa_numerical_failure
         return
      end if

      allocate(out%sigma(i))
      allocate(out%u(window, i))
      allocate(out%v(k, i))
      out%sigma = singular_values(:i)
      out%u = u(:, :i) * spread(col_inv, 2, i)
      out%v = transpose(vt(:i, :)) * spread(row_inv, 2, i)
      do i = 1, size(out%sigma)
         nu_col = sqrt(sum(out%u(:, i) * out%u(:, i)))
         nv_col = sqrt(sum(out%v(:, i) * out%v(:, i)))
         if (nu_col <= epsilon(1.0_dp) .or. nv_col <= epsilon(1.0_dp)) then
            out%info = rssa_numerical_failure
            return
         end if
         out%u(:, i) = out%u(:, i) / nu_col
         out%v(:, i) = out%v(:, i) / nv_col
         out%sigma(i) = out%sigma(i) * nu_col * nv_col
      end do
      out%info = rssa_success
   end function decompose_wossa

   function owcor_ssa(object, group_labels, basis_indices, info) result(matrix)
      type(ssa_result), intent(in) :: object !! Oblique SSA decomposition containing U, V, sigma, series, and window.
      integer, intent(in), optional :: group_labels(:) !! Positive group label for each selected basis component.
      integer, intent(in), optional :: basis_indices(:) !! One-based OSSA basis indices; defaults to all stored
      !! components.
      integer, intent(out), optional :: info !! Zero on success or a translation status code on invalid/numerical input.
      real(dp), allocatable :: matrix(:, :)
      integer, allocatable :: idx(:), labels(:)
      real(dp), allocatable :: cov(:, :), grouped(:, :), hf(:, :), lm(:, :), mx(:, :), rm(:, :), transformed(:, :)
      real(dp) :: denom
      integer :: i, j, linear_info, n_groups, n_basis

      call prepare_groups(object, group_labels, basis_indices, idx, labels, grouped, linear_info)
      if (linear_info /= rssa_success) then
         allocate(matrix(0, 0))
         if (present(info)) info = linear_info
         return
      end if
      n_basis = size(idx)
      n_groups = size(grouped, 2)
      call pseudo_inverse(object%u(:, idx), lm, linear_info)
      if (linear_info /= 0) then
         allocate(matrix(0, 0))
         if (present(info)) info = rssa_numerical_failure
         return
      end if
      call pseudo_inverse(object%v(:, idx), rm, linear_info)
      if (linear_info /= 0) then
         allocate(matrix(0, 0))
         if (present(info)) info = rssa_numerical_failure
         return
      end if

      allocate(mx(n_basis * n_basis, n_groups))
      do j = 1, n_groups
         hf = hankel_matrix(grouped(:, j), object%window)
         transformed = matmul(lm, matmul(hf, transpose(rm)))
         mx(:, j) = reshape(transformed, [n_basis * n_basis])
      end do
      cov = matmul(transpose(mx), mx)
      allocate(matrix(n_groups, n_groups))
      matrix = 0.0_dp
      do j = 1, n_groups
         do i = 1, n_groups
            denom = sqrt(max(0.0_dp, cov(i, i) * cov(j, j)))
            if (denom > 0.0_dp) matrix(i, j) = max(-1.0_dp, min(1.0_dp, cov(i, j) / denom))
         end do
      end do
      if (present(info)) info = rssa_success
   end function owcor_ssa

   function wcor_ossa(object, group_labels, basis_indices, info) result(matrix)
      type(ssa_result), intent(in) :: object !! Oblique SSA decomposition whose reconstructed groups are correlated.
      integer, intent(in), optional :: group_labels(:) !! Positive group label for each selected basis component.
      integer, intent(in), optional :: basis_indices(:) !! One-based OSSA basis indices; defaults to all stored
      !! components.
      integer, intent(out), optional :: info !! Zero on success or invalid-input status.
      real(dp), allocatable :: matrix(:, :)
      integer, allocatable :: idx(:), labels(:)
      real(dp), allocatable :: grouped(:, :)
      integer :: i, j, local_info, n_groups

      call prepare_groups(object, group_labels, basis_indices, idx, labels, grouped, local_info)
      if (local_info /= rssa_success) then
         allocate(matrix(0, 0))
         if (present(info)) info = local_info
         return
      end if
      n_groups = size(grouped, 2)
      allocate(matrix(n_groups, n_groups))
      do j = 1, n_groups
         do i = 1, n_groups
            matrix(i, j) = wcor_default(grouped(:, i), grouped(:, j), object%window)
         end do
      end do
      if (present(info)) info = rssa_success
   end function wcor_ossa


   subroutine filter_columns(vectors, impulse, filtered, info)
      real(dp), intent(in) :: vectors(:, :) !! Right-vector columns to filter independently.
      real(dp), intent(in) :: impulse(:) !! Reversed finite impulse-response coefficients.
      real(dp), allocatable, intent(out) :: filtered(:, :) !! Valid filtered columns with K-q+1 rows.
      integer, intent(out) :: info !! Zero on success or invalid-input status.
      integer :: i, j, q, s

      q = size(impulse)
      if (q < 1 .or. q > size(vectors, 1)) then
         allocate(filtered(0, 0))
         info = rssa_invalid_input
         return
      end if
      allocate(filtered(size(vectors, 1) - q + 1, size(vectors, 2)))
      filtered = 0.0_dp
      do j = 1, size(vectors, 2)
         do i = 1, size(filtered, 1)
            do s = 1, q
               filtered(i, j) = filtered(i, j) + impulse(s) * vectors(i + s - 1, j)
            end do
         end do
      end do
      info = rssa_success
   end subroutine filter_columns

   subroutine selected_indices(n_available, requested, idx, info)
      integer, intent(in) :: n_available !! Number of stored eigentriples available for selection.
      integer, intent(in), optional :: requested(:) !! Optional one-based eigentriple indices.
      integer, allocatable, intent(out) :: idx(:) !! Sorted unique selected indices.
      integer, intent(out) :: info !! Zero on success or invalid-input status.
      integer :: i, j, tmp

      if (present(requested)) then
         allocate(idx(size(requested)))
         idx = requested
      else
         allocate(idx(n_available))
         idx = [(i, i = 1, n_available)]
      end if
      if (size(idx) < 1 .or. any(idx < 1) .or. any(idx > n_available)) then
         info = rssa_invalid_input
         return
      end if
      do i = 2, size(idx)
         tmp = idx(i)
         j = i - 1
         do while (j >= 1)
            if (idx(j) <= tmp) exit
            idx(j + 1) = idx(j)
            j = j - 1
         end do
         idx(j + 1) = tmp
      end do
      do i = 2, size(idx)
         if (idx(i) == idx(i - 1)) then
            info = rssa_invalid_input
            return
         end if
      end do
      info = rssa_success
   end subroutine selected_indices

   pure subroutine reciprocal_scale(scale, inverse)
      real(dp), intent(in) :: scale(:) !! Nonnegative square-root metric scale values.
      real(dp), intent(out) :: inverse(:) !! Pseudoinverse scale, zero where the input is numerically zero.
      integer :: i

      do i = 1, size(scale)
         if (abs(scale(i)) >= 1.0e-6_dp) then
            inverse(i) = 1.0_dp / scale(i)
         else
            inverse(i) = 0.0_dp
         end if
      end do
   end subroutine reciprocal_scale

   subroutine pseudo_inverse(a, inverse, info)
      real(dp), intent(in) :: a(:, :) !! Matrix whose Moore-Penrose pseudoinverse is required.
      real(dp), allocatable, intent(out) :: inverse(:, :) !! SVD least-squares pseudoinverse with shape (n,m).
      integer, intent(out) :: info !! Zero on successful SVD least-squares solution.
      real(dp), allocatable :: identity(:, :)
      integer :: i, rank

      allocate(identity(size(a, 1), size(a, 1)))
      identity = 0.0_dp
      do i = 1, size(identity, 1)
         identity(i, i) = 1.0_dp
      end do
      allocate(inverse(size(a, 2), size(a, 1)))
      call least_squares_svd(a, identity, inverse, rank, info)
   end subroutine pseudo_inverse

   subroutine prepare_groups(object, group_labels, basis_indices, idx, labels, grouped, info)
      type(ssa_result), intent(in) :: object !! SSA decomposition from which elementary reconstructed series are grouped.
      integer, intent(in), optional :: group_labels(:) !! Positive group label corresponding to each selected component.
      integer, intent(in), optional :: basis_indices(:) !! Selected one-based component indices; defaults to all
      !! components.
      integer, allocatable, intent(out) :: idx(:) !! Concrete selected component indices used by the calculation.
      integer, allocatable, intent(out) :: labels(:) !! Concrete positive group labels used by the calculation.
      real(dp), allocatable, intent(out) :: grouped(:, :) !! Reconstructed group series, one group per column.
      integer, intent(out) :: info !! Zero on success or invalid-input status.
      real(dp), allocatable :: one(:)
      integer :: g, i, n_groups
      logical, allocatable :: seen(:)

      info = rssa_invalid_input
      if (.not. allocated(object%series) .or. .not. allocated(object%sigma)) return
      if (.not. allocated(object%u) .or. .not. allocated(object%v)) return
      if (object%window < 1 .or. object%window > size(object%series)) return
      if (present(basis_indices)) then
         allocate(idx(size(basis_indices)))
         idx = basis_indices
      else
         allocate(idx(size(object%sigma)))
         idx = [(i, i = 1, size(idx))]
      end if
      if (size(idx) < 1 .or. any(idx < 1) .or. any(idx > size(object%sigma))) return
      if (present(group_labels)) then
         if (size(group_labels) /= size(idx) .or. any(group_labels < 1)) return
         allocate(labels(size(group_labels)))
         labels = group_labels
      else
         allocate(labels(size(idx)))
         labels = [(i, i = 1, size(labels))]
      end if
      n_groups = maxval(labels)
      allocate(seen(n_groups))
      seen = .false.
      do i = 1, size(labels)
         seen(labels(i)) = .true.
      end do
      if (.not. all(seen)) return
      allocate(grouped(size(object%series), n_groups))
      grouped = 0.0_dp
      do i = 1, size(idx)
         g = labels(i)
         one = elementary_series_ssa(object, idx(i))
         if (size(one) /= size(object%series)) return
         grouped(:, g) = grouped(:, g) + one
      end do
      info = rssa_success
   end subroutine prepare_groups
end module rssa_oblique
