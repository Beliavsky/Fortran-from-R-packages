! SPDX-License-Identifier: GPL-2.0-or-later
! Iterative and ESPRIT-based oblique SSA kernels translated from Rssa 1.1.
! Derived from upstream R/ossa.R (Copyright (c) 2014 Alex Shlemov) and
! R/eossa.R (Copyright (c) 2017-2018 Alex Shlemov).
module rssa_iterative_oblique
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
   use rssa_decomposition, only : ssa_1d
   use rssa_kinds, only : dp
   use rssa_matrices, only : hankelize_matrix
   use rssa_reconstruction, only : elementary_series_ssa
   use rssa_types, only : rssa_invalid_input, rssa_not_supported, rssa_numerical_failure, rssa_success
   use rssa_types, only : ssa_result
   use r_linalg, only : general_complex_eigen, least_squares_svd, thin_svd
   use svd, only : propack_svd_dense, propack_svd_result
   implicit none
   private
   public :: eossa, eossa_ssa, iossa, iossa_ssa

   interface iossa
      module procedure iossa_ssa
   end interface iossa

   interface eossa
      module procedure eossa_ssa
   end interface eossa
contains

   function iossa_ssa(object, indices, group_labels, tol, kappa, maxiter, kappa_balance, iterations, converged) result(out)
      type(ssa_result), intent(in) :: object !! Real 1-D SSA decomposition whose selected eigentriples are refined.
      integer, intent(in), optional :: indices(:) !! One-based eigentriples to refine; defaults to all stored components.
      integer, intent(in), optional :: group_labels(:) !! Positive nested-group label for every selected eigentriple.
      real(dp), intent(in), optional :: tol !! RMS reconstruction-change tolerance; default is 1.0e-5.
      real(dp), intent(in), optional :: kappa !! Adjacent-group singular-value separation factor; default is 2.
      integer, intent(in), optional :: maxiter !! Maximum I-OSSA iterations; default is 100.
      real(dp), intent(in), optional :: kappa_balance !! Fraction of separation rescaling applied to left vectors; default 0.5.
      integer, intent(out), optional :: iterations !! Number of iterations performed, including the converged iteration.
      logical, intent(out), optional :: converged !! True when the maximum group RMS change falls below the tolerance.
      type(ssa_result) :: out
      integer, allocatable :: idx(:), labels(:), ranks(:), starts(:)
      real(dp), allocatable :: fs(:, :), fs_new(:, :), original_sigma(:), original_u(:, :), original_v(:, :)
      real(dp), allocatable :: basis_u(:, :), basis_v(:, :), work_sigma(:), work_u(:, :), work_v(:, :)
      real(dp), allocatable :: sigma(:), y(:, :), z(:, :), group_sigma(:), group_u(:, :), group_v(:, :)
      real(dp), allocatable :: one(:)
      real(dp) :: balance_value, delta, max_delta, separation, tolerance_value, ratio, multiplier
      integer :: g, i, info_local, iter, max_iterations, n_groups, n_selected, pos
      logical :: did_converge
      type(ssa_result) :: group_ssa

      out = object
      out%info = rssa_invalid_input
      did_converge = .false.
      iter = 0
      if (present(iterations)) iterations = 0
      if (present(converged)) converged = .false.
      if (.not. valid_real_ssa(object)) return

      call select_components(size(object%sigma), indices, idx, info_local)
      if (info_local /= rssa_success) return
      n_selected = size(idx)
      call concrete_group_labels(n_selected, group_labels, labels, n_groups, info_local)
      if (info_local /= rssa_success) return
      allocate(ranks(n_groups), starts(n_groups))
      do g = 1, n_groups
         ranks(g) = count(labels == g)
      end do
      starts(1) = 1
      do g = 2, n_groups
         starts(g) = starts(g - 1) + ranks(g - 1)
      end do

      tolerance_value = 1.0e-5_dp
      if (present(tol)) tolerance_value = tol
      if (.not. ieee_is_finite(tolerance_value) .or. tolerance_value < 0.0_dp) return
      max_iterations = 100
      if (present(maxiter)) max_iterations = maxiter
      if (max_iterations < 1) return
      separation = 2.0_dp
      if (present(kappa)) separation = kappa
      if (ieee_is_nan(separation)) return
      balance_value = 0.5_dp
      if (present(kappa_balance)) balance_value = kappa_balance
      if (.not. ieee_is_finite(balance_value)) return
      if (balance_value < 0.0_dp .or. balance_value > 1.0_dp) return

      call orthogonal_selected(object, idx, original_sigma, original_u, original_v, info_local)
      if (info_local /= rssa_success) then
         out%info = info_local
         return
      end if

      allocate(fs(size(object%series), n_groups), fs_new(size(object%series), n_groups))
      fs = 0.0_dp
      do i = 1, n_selected
         one = elementary_series_ssa(object, idx(i))
         if (size(one) /= size(object%series)) return
         fs(:, labels(i)) = fs(:, labels(i)) + one
      end do

      allocate(basis_u(size(object%u, 1), n_selected), basis_v(size(object%v, 1), n_selected))
      allocate(work_sigma(n_selected), work_u(size(object%u, 1), n_selected), work_v(size(object%v, 1), n_selected))
      sigma = original_sigma
      y = original_u
      z = original_v

      do iter = 1, max_iterations
         pos = 1
         do g = 1, n_groups
            group_ssa = ssa_1d(fs(:, g), window=object%window, neig=ranks(g))
            if (group_ssa%info /= rssa_success) then
               out%info = group_ssa%info
               if (present(iterations)) iterations = iter
               return
            end if
            if (size(group_ssa%sigma) < ranks(g)) then
               out%info = rssa_numerical_failure
               if (present(iterations)) iterations = iter
               return
            end if
            group_sigma = group_ssa%sigma(1:ranks(g))
            group_u = group_ssa%u(:, 1:ranks(g))
            group_v = group_ssa%v(:, 1:ranks(g))
            work_sigma(pos:pos + ranks(g) - 1) = group_sigma
            work_u(:, pos:pos + ranks(g) - 1) = group_u
            work_v(:, pos:pos + ranks(g) - 1) = group_v
            pos = pos + ranks(g)
         end do

         if (separation > 0.0_dp) then
            do g = 1, n_groups - 1
               i = starts(g) + ranks(g) - 1
               pos = starts(g + 1)
               if (work_sigma(pos) <= epsilon(1.0_dp)) cycle
               ratio = work_sigma(i) / work_sigma(pos)
               if (ratio < separation) then
                  if (ratio <= epsilon(1.0_dp)) then
                     out%info = rssa_numerical_failure
                     if (present(iterations)) iterations = iter
                     return
                  end if
                  multiplier = separation / ratio
                  work_sigma(pos:pos + ranks(g + 1) - 1) = &
                     work_sigma(pos:pos + ranks(g + 1) - 1) / multiplier
                  work_u(:, pos:pos + ranks(g + 1) - 1) = &
                     work_u(:, pos:pos + ranks(g + 1) - 1) * multiplier ** balance_value
                  work_v(:, pos:pos + ranks(g + 1) - 1) = &
                     work_v(:, pos:pos + ranks(g + 1) - 1) * multiplier ** (1.0_dp - balance_value)
               end if
            end do
         end if

         basis_u = work_u
         basis_v = work_v
         call svd_to_low_rank_svd(original_sigma, original_u, original_v, basis_u, basis_v, &
                                  sigma, y, z, info_local)
         if (info_local /= rssa_success) then
            out%info = info_local
            if (present(iterations)) iterations = iter
            return
         end if

         do g = 1, n_groups
            call low_rank_group_series(sigma, y, z, starts(g), ranks(g), fs_new(:, g))
         end do
         max_delta = 0.0_dp
         do g = 1, n_groups
            delta = sqrt(sum((fs_new(:, g) - fs(:, g)) ** 2) / real(size(fs, 1), dp))
            max_delta = max(max_delta, delta)
         end do
         fs = fs_new
         if (max_delta < tolerance_value) then
            did_converge = .true.
            exit
         end if
      end do

      call save_oblique(object, idx, sigma, y, z, out, info_local)
      if (info_local /= rssa_success) then
         out%info = info_local
         if (present(iterations)) iterations = iter
         if (present(converged)) converged = did_converge
         return
      end if
      out%info = rssa_success
      if (present(iterations)) iterations = iter
      if (present(converged)) converged = did_converge
   end function iossa_ssa

   function eossa_ssa(object, indices, k, subspace, solve_method, cluster_labels) result(out)
      type(ssa_result), intent(in) :: object !! Real 1-D SSA decomposition whose selected eigentriples are ESPRIT-rotated.
      integer, intent(in), optional :: indices(:) !! One-based eigentriples used by EOSSA; defaults to all stored components.
      integer, intent(in), optional :: k !! Number of ESPRIT root clusters; default is 2.
      character(len=*), intent(in), optional :: subspace !! ESPRIT subspace; the translated 1-D kernel supports "column".
      character(len=*), intent(in), optional :: solve_method !! Shift solve; the translated kernel supports least squares "ls".
      integer, allocatable, intent(out), optional :: cluster_labels(:) !! Cluster label for each sorted selected eigentriple.
      type(ssa_result) :: out
      integer, allocatable :: idx(:), labels(:)
      real(dp), allocatable :: c(:, :), shifted(:, :), v_scaled(:, :), x(:, :)
      real(dp), allocatable :: y(:, :), z(:, :), sigma(:)
      complex(dp), allocatable :: eigvec(:, :), roots(:)
      integer :: info_local, n_clusters, n_selected, rank_ls
      character(len=:), allocatable :: space, solver

      out = object
      out%info = rssa_invalid_input
      if (present(cluster_labels)) allocate(cluster_labels(0))
      if (.not. valid_real_ssa(object)) return
      call select_components(size(object%sigma), indices, idx, info_local)
      if (info_local /= rssa_success) return
      n_selected = size(idx)
      if (n_selected < 1) return

      n_clusters = 2
      if (present(k)) n_clusters = k
      if (n_clusters < 1 .or. n_clusters > n_selected) return
      space = "column"
      if (present(subspace)) space = trim(adjustl(subspace))
      solver = "ls"
      if (present(solve_method)) solver = trim(adjustl(solve_method))
      if (space /= "column") then
         out%info = rssa_not_supported
         return
      end if
      if (solver /= "ls") then
         out%info = rssa_not_supported
         return
      end if
      if (size(object%u, 1) < 2) return

      call shift_matrix_1d(object%u(:, idx), shifted, info_local)
      if (info_local /= rssa_success) then
         out%info = info_local
         return
      end if
      call general_complex_eigen(cmplx(shifted, 0.0_dp, kind=dp), roots, eigvec, info_local)
      if (info_local /= 0) then
         out%info = rssa_numerical_failure
         return
      end if
      call sort_eigenpairs_frequency(roots, eigvec)
      call cluster_real_basis(eigvec, roots, n_clusters, c, labels, info_local)
      if (info_local /= rssa_success) then
         out%info = info_local
         return
      end if

      y = matmul(object%u(:, idx), c)
      allocate(v_scaled(size(object%v, 1), n_selected))
      v_scaled = object%v(:, idx) * spread(object%sigma(idx), 1, size(object%v, 1))
      allocate(x(n_selected, size(v_scaled, 1)))
      call least_squares_svd(c, transpose(v_scaled), x, rank_ls, info_local)
      if (info_local /= 0 .or. rank_ls < n_selected) then
         out%info = rssa_numerical_failure
         return
      end if
      z = transpose(x)
      allocate(sigma(n_selected))
      sigma = 1.0_dp
      call save_oblique(object, idx, sigma, y, z, out, info_local)
      if (info_local /= rssa_success) then
         out%info = info_local
         return
      end if
      out%info = rssa_success
      if (present(cluster_labels)) then
         if (allocated(cluster_labels)) deallocate(cluster_labels)
         allocate(cluster_labels(size(labels)))
         cluster_labels = labels
      end if
   end function eossa_ssa

   pure logical function valid_real_ssa(object) result(valid)
      type(ssa_result), intent(in) :: object !! Candidate real SSA decomposition to validate.

      valid = allocated(object%series) .and. allocated(object%sigma)
      if (.not. valid) return
      valid = allocated(object%u) .and. allocated(object%v)
      if (.not. valid) return
      valid = object%window >= 1 .and. object%window <= size(object%series)
      if (.not. valid) return
      valid = size(object%u, 2) == size(object%sigma)
      if (.not. valid) return
      valid = size(object%v, 2) == size(object%sigma)
   end function valid_real_ssa

   subroutine select_components(n_available, requested, idx, info)
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
      if (size(idx) < 1) then
         info = rssa_invalid_input
         return
      end if
      if (any(idx < 1) .or. any(idx > n_available)) then
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
   end subroutine select_components

   subroutine concrete_group_labels(n_selected, supplied, labels, n_groups, info)
      integer, intent(in) :: n_selected !! Number of selected eigentriples requiring group membership.
      integer, intent(in), optional :: supplied(:) !! Optional positive group label for each selected eigentriple.
      integer, allocatable, intent(out) :: labels(:) !! Concrete contiguous group labels in 1..n_groups.
      integer, intent(out) :: n_groups !! Number of distinct nested groups.
      integer, intent(out) :: info !! Zero on success or invalid-input status.
      logical, allocatable :: seen(:)
      integer :: i

      info = rssa_invalid_input
      if (n_selected < 1) return
      allocate(labels(n_selected))
      if (present(supplied)) then
         if (size(supplied) /= n_selected) return
         if (any(supplied < 1)) return
         labels = supplied
      else
         labels = [(i, i = 1, n_selected)]
      end if
      n_groups = maxval(labels)
      allocate(seen(n_groups))
      seen = .false.
      do i = 1, n_selected
         seen(labels(i)) = .true.
      end do
      if (.not. all(seen)) return
      info = rssa_success
   end subroutine concrete_group_labels

   subroutine orthogonal_selected(object, idx, sigma, u, v, info)
      type(ssa_result), intent(in) :: object !! SSA decomposition containing the selected low-rank trajectory terms.
      integer, intent(in) :: idx(:) !! One-based selected eigentriples to orthogonalize jointly.
      real(dp), allocatable, intent(out) :: sigma(:) !! Orthogonal singular values of the selected low-rank matrix.
      real(dp), allocatable, intent(out) :: u(:, :) !! Orthogonal left singular vectors.
      real(dp), allocatable, intent(out) :: v(:, :) !! Orthogonal right singular vectors.
      integer, intent(out) :: info !! Zero on success or numerical-failure status.
      type(propack_svd_result) :: dec
      real(dp), allocatable :: matrix(:, :)
      integer :: j

      allocate(matrix(size(object%u, 1), size(object%v, 1)))
      matrix = 0.0_dp
      do j = 1, size(idx)
         matrix = matrix + object%sigma(idx(j)) * spread(object%u(:, idx(j)), 2, size(matrix, 2)) * &
                  spread(object%v(:, idx(j)), 1, size(matrix, 1))
      end do
      dec = propack_svd_dense(matrix, size(idx))
      if (dec%info /= 0) then
         allocate(sigma(0), u(0, 0), v(0, 0))
         info = rssa_numerical_failure
         return
      end if
      sigma = dec%d
      u = dec%u
      v = dec%v
      info = rssa_success
   end subroutine orthogonal_selected

   subroutine svd_to_low_rank_svd(d, u, v, basis_l, basis_r, sigma, y, z, info)
      real(dp), intent(in) :: d(:) !! Singular values of the original orthogonal low-rank decomposition.
      real(dp), intent(in) :: u(:, :) !! Original orthonormal left singular vectors.
      real(dp), intent(in) :: v(:, :) !! Original orthonormal right singular vectors.
      real(dp), intent(in) :: basis_l(:, :) !! Candidate left basis spanning the target component subspaces.
      real(dp), intent(in) :: basis_r(:, :) !! Candidate right basis spanning the target component subspaces.
      real(dp), allocatable, intent(out) :: sigma(:) !! Oblique singular coefficients in decreasing order.
      real(dp), allocatable, intent(out) :: y(:, :) !! Projected oblique left vectors.
      real(dp), allocatable, intent(out) :: z(:, :) !! Projected oblique right vectors.
      integer, intent(out) :: info !! Zero on success or numerical-failure status.
      real(dp), allocatable :: ub(:, :), vb(:, :), core(:, :), left_projected(:, :), right_projected(:, :)
      real(dp), allocatable :: reversed_u(:, :), reversed_v(:, :)
      type(propack_svd_result) :: dec
      integer :: i, n

      n = size(d)
      if (n < 1 .or. size(u, 2) /= n .or. size(v, 2) /= n) then
         allocate(sigma(0), y(0, 0), z(0, 0))
         info = rssa_invalid_input
         return
      end if
      if (size(basis_l, 2) /= n .or. size(basis_r, 2) /= n) then
         allocate(sigma(0), y(0, 0), z(0, 0))
         info = rssa_invalid_input
         return
      end if
      if (any(d <= epsilon(1.0_dp))) then
         allocate(sigma(0), y(0, 0), z(0, 0))
         info = rssa_numerical_failure
         return
      end if
      ub = matmul(transpose(u), basis_l)
      vb = matmul(transpose(v), basis_r)
      core = matmul(transpose(ub), vb / spread(d, 2, n))
      dec = propack_svd_dense(core, n)
      if (dec%info /= 0 .or. any(dec%d <= epsilon(1.0_dp))) then
         allocate(sigma(0), y(0, 0), z(0, 0))
         info = rssa_numerical_failure
         return
      end if
      allocate(sigma(n), reversed_u(n, n), reversed_v(n, n))
      do i = 1, n
         sigma(i) = 1.0_dp / dec%d(n - i + 1)
         reversed_u(:, i) = dec%u(:, n - i + 1)
         reversed_v(:, i) = dec%v(:, n - i + 1)
      end do
      left_projected = matmul(u, ub)
      right_projected = matmul(v, vb)
      y = matmul(left_projected, reversed_u)
      z = matmul(right_projected, reversed_v)
      info = rssa_success
   end subroutine svd_to_low_rank_svd

   subroutine low_rank_group_series(sigma, y, z, first, count, series)
      real(dp), intent(in) :: sigma(:) !! Oblique low-rank coefficients.
      real(dp), intent(in) :: y(:, :) !! Oblique left vectors.
      real(dp), intent(in) :: z(:, :) !! Oblique right vectors.
      integer, intent(in) :: first !! One-based first component belonging to the group.
      integer, intent(in) :: count !! Number of consecutive components belonging to the group.
      real(dp), intent(out) :: series(:) !! Hankelized reconstructed group series.
      real(dp), allocatable :: matrix(:, :), one(:)
      integer :: j

      allocate(matrix(size(y, 1), size(z, 1)))
      matrix = 0.0_dp
      do j = first, first + count - 1
         matrix = matrix + sigma(j) * spread(y(:, j), 2, size(z, 1)) * spread(z(:, j), 1, size(y, 1))
      end do
      one = hankelize_matrix(matrix)
      if (size(one) == size(series)) then
         series = one
      else
         series = 0.0_dp
      end if
   end subroutine low_rank_group_series

   subroutine save_oblique(original, idx, sigma, y, z, out, info)
      type(ssa_result), intent(in) :: original !! Original decomposition to copy and update at selected positions.
      integer, intent(in) :: idx(:) !! One-based positions receiving the oblique eigentriples.
      real(dp), intent(in) :: sigma(:) !! Unnormalized oblique coefficients.
      real(dp), intent(in) :: y(:, :) !! Unnormalized oblique left vectors.
      real(dp), intent(in) :: z(:, :) !! Unnormalized oblique right vectors.
      type(ssa_result), intent(out) :: out !! Updated decomposition with Euclidean-normalized stored vectors.
      integer, intent(out) :: info !! Zero on success or numerical-failure status.
      real(dp) :: ynorm, znorm
      integer :: i

      out = original
      if (size(idx) /= size(sigma) .or. size(y, 2) /= size(idx) .or. size(z, 2) /= size(idx)) then
         out%info = rssa_invalid_input
         info = out%info
         return
      end if
      do i = 1, size(idx)
         ynorm = sqrt(sum(y(:, i) ** 2))
         znorm = sqrt(sum(z(:, i) ** 2))
         if (ynorm <= epsilon(1.0_dp) .or. znorm <= epsilon(1.0_dp)) then
            out%info = rssa_numerical_failure
            info = out%info
            return
         end if
         out%sigma(idx(i)) = sigma(i) * ynorm * znorm
         out%u(:, idx(i)) = y(:, i) / ynorm
         out%v(:, idx(i)) = z(:, i) / znorm
      end do
      if (minval(idx) <= original%n_special_right) then
         out%n_special_right = 0
         out%n_special_left = 0
      else if (minval(idx) <= original%n_special_right + original%n_special_left) then
         out%n_special_left = 0
      end if
      out%info = rssa_success
      info = rssa_success
   end subroutine save_oblique

   subroutine shift_matrix_1d(u, shift, info)
      real(dp), intent(in) :: u(:, :) !! Subspace basis whose one-step shift matrix is estimated.
      real(dp), allocatable, intent(out) :: shift(:, :) !! Least-squares shift operator in the supplied basis.
      integer, intent(out) :: info !! Zero on success or numerical-failure status.
      real(dp), allocatable :: left(:, :), right(:, :)
      integer :: rank_ls

      if (size(u, 1) < 2 .or. size(u, 2) < 1) then
         allocate(shift(0, 0))
         info = rssa_invalid_input
         return
      end if
      left = u(1:size(u, 1) - 1, :)
      right = u(2:size(u, 1), :)
      allocate(shift(size(u, 2), size(u, 2)))
      call least_squares_svd(left, right, shift, rank_ls, info)
      if (info /= 0 .or. rank_ls < size(u, 2)) then
         shift = 0.0_dp
         info = rssa_numerical_failure
      else
         info = rssa_success
      end if
   end subroutine shift_matrix_1d

   subroutine sort_eigenpairs_frequency(values, vectors)
      complex(dp), intent(inout) :: values(:) !! Eigenvalues sorted in ascending absolute argument.
      complex(dp), intent(inout) :: vectors(:, :) !! Right eigenvectors permuted consistently with eigenvalues.
      complex(dp) :: value_tmp
      complex(dp), allocatable :: vector_tmp(:)
      real(dp) :: key_frequency
      integer :: i, j

      allocate(vector_tmp(size(vectors, 1)))
      do i = 2, size(values)
         value_tmp = values(i)
         vector_tmp = vectors(:, i)
         key_frequency = abs(atan2(aimag(value_tmp), real(value_tmp, dp)))
         j = i - 1
         do while (j >= 1)
            if (abs(atan2(aimag(values(j)), real(values(j), dp))) <= key_frequency) exit
            values(j + 1) = values(j)
            vectors(:, j + 1) = vectors(:, j)
            j = j - 1
         end do
         values(j + 1) = value_tmp
         vectors(:, j + 1) = vector_tmp
      end do
   end subroutine sort_eigenpairs_frequency

   subroutine cluster_real_basis(vectors, roots, k, basis, labels, info)
      complex(dp), intent(in) :: vectors(:, :) !! Frequency-sorted complex ESPRIT eigenvectors.
      complex(dp), intent(in) :: roots(:) !! Frequency-sorted ESPRIT roots used for complete-link clustering.
      integer, intent(in) :: k !! Number of requested root clusters.
      real(dp), allocatable, intent(out) :: basis(:, :) !! Real basis spanning each clustered complex eigenspace.
      integer, allocatable, intent(out) :: labels(:) !! Contiguous cluster labels for sorted roots.
      integer, intent(out) :: info !! Zero on success or numerical-failure status.
      integer, allocatable :: members(:, :), counts(:), cluster_order(:), cluster_map(:), group_idx(:)
      logical, allocatable :: active(:)
      real(dp), allocatable :: realified(:, :), singular_values(:), u(:, :), vt(:, :)
      real(dp) :: distance, best_distance
      integer :: a, b, best_a, best_b, c, i, j, n, n_active, m, label

      n = size(roots)
      if (n < 1 .or. k < 1 .or. k > n .or. size(vectors, 2) /= n) then
         allocate(basis(0, 0), labels(0))
         info = rssa_invalid_input
         return
      end if
      allocate(members(n, n), counts(n), active(n))
      members = 0
      counts = 1
      active = .true.
      do i = 1, n
         members(1, i) = i
      end do
      n_active = n
      do while (n_active > k)
         best_distance = huge(1.0_dp)
         best_a = 0
         best_b = 0
         do a = 1, n - 1
            if (.not. active(a)) cycle
            do b = a + 1, n
               if (.not. active(b)) cycle
               distance = complete_link_distance(roots, members(:, a), counts(a), members(:, b), counts(b))
               if (distance < best_distance) then
                  best_distance = distance
                  best_a = a
                  best_b = b
               end if
            end do
         end do
         if (best_a == 0 .or. best_b == 0) then
            allocate(basis(0, 0), labels(0))
            info = rssa_numerical_failure
            return
         end if
         members(counts(best_a) + 1:counts(best_a) + counts(best_b), best_a) = members(1:counts(best_b), best_b)
         counts(best_a) = counts(best_a) + counts(best_b)
         active(best_b) = .false.
         n_active = n_active - 1
      end do

      allocate(cluster_order(k), cluster_map(n), labels(n))
      cluster_order = 0
      c = 0
      do a = 1, n
         if (.not. active(a)) cycle
         c = c + 1
         cluster_order(c) = a
      end do
      call sort_clusters_by_first_member(cluster_order, members, counts)
      cluster_map = 0
      do label = 1, k
         a = cluster_order(label)
         do j = 1, counts(a)
            cluster_map(members(j, a)) = label
         end do
      end do
      labels = cluster_map

      allocate(basis(size(vectors, 1), n))
      basis = 0.0_dp
      do label = 1, k
         m = count(labels == label)
         allocate(group_idx(m))
         c = 0
         do i = 1, n
            if (labels(i) == label) then
               c = c + 1
               group_idx(c) = i
            end if
         end do
         allocate(realified(size(vectors, 1), 2 * m))
         do j = 1, m
            realified(:, j) = real(vectors(:, group_idx(j)), dp)
            realified(:, m + j) = aimag(vectors(:, group_idx(j)))
         end do
         call thin_svd(realified, u, singular_values, vt, info)
         if (info /= 0) then
            deallocate(group_idx, realified)
            info = rssa_numerical_failure
            return
         end if
         do j = 1, m
            basis(:, group_idx(j)) = u(:, j)
         end do
         deallocate(group_idx, realified, singular_values, u, vt)
      end do
      info = rssa_success
   end subroutine cluster_real_basis

   pure real(dp) function complete_link_distance(roots, left_members, n_left, right_members, n_right) result(distance)
      complex(dp), intent(in) :: roots(:) !! ESPRIT roots represented in the clustering plane.
      integer, intent(in) :: left_members(:) !! Root indices belonging to the first cluster.
      integer, intent(in) :: n_left !! Active number of entries in left_members.
      integer, intent(in) :: right_members(:) !! Root indices belonging to the second cluster.
      integer, intent(in) :: n_right !! Active number of entries in right_members.
      real(dp) :: dx, dy
      integer :: i, j

      distance = 0.0_dp
      do i = 1, n_left
         do j = 1, n_right
            dx = real(roots(left_members(i)), dp) - real(roots(right_members(j)), dp)
            dy = abs(aimag(roots(left_members(i)))) - abs(aimag(roots(right_members(j))))
            distance = max(distance, hypot(dx, dy))
         end do
      end do
   end function complete_link_distance

   subroutine sort_clusters_by_first_member(order, members, counts)
      integer, intent(inout) :: order(:) !! Active cluster IDs sorted by their minimum root index.
      integer, intent(in) :: members(:, :) !! Root membership table for all working clusters.
      integer, intent(in) :: counts(:) !! Active member count for every working cluster.
      integer :: i, j, tmp, key, current

      do i = 2, size(order)
         tmp = order(i)
         key = minval(members(1:counts(tmp), tmp))
         j = i - 1
         do while (j >= 1)
            current = minval(members(1:counts(order(j)), order(j)))
            if (current <= key) exit
            order(j + 1) = order(j)
            j = j - 1
         end do
         order(j + 1) = tmp
      end do
   end subroutine sort_clusters_by_first_member

end module rssa_iterative_oblique
