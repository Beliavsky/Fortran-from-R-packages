! SPDX-License-Identifier: GPL-2.0-or-later
! Linear recurrent and vector forecasting kernels translated from Rssa 1.1.
module rssa_forecast
   use rssa_kinds, only : dp
   use rssa_matrices, only : hankelize_matrix
   use rssa_reconstruction, only : reconstruct_complex, reconstruct_mssa, reconstruct_ssa
   use rssa_types, only : cssa_result, mssa_result, rssa_invalid_input, rssa_numerical_failure
   use rssa_types, only : rssa_success, ssa_result
   use r_linalg, only : general_complex_eigen, general_real_eigenvalues, least_squares_svd, thin_qr
   implicit none
   private
   type, public :: bootstrap_forecast_result
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: lower(:)
      real(dp), allocatable :: upper(:)
      integer :: info = rssa_success
   end type bootstrap_forecast_result
   public :: lrr, lrr_default, lrr_ssa, lrr_mssa, lrr_complex
   public :: roots_lrr, roots_lrr_complex, apply_lrr, apply_lrr_complex
   public :: rforecast_ssa, rforecast_mssa, rforecast_complex, rforecast_pssa
   public :: vforecast_ssa, vforecast_mssa, vforecast_complex, vforecast_pssa, bforecast_ssa
   interface lrr
      module procedure lrr_default
      module procedure lrr_ssa
      module procedure lrr_mssa
      module procedure lrr_complex
   end interface lrr
contains
   function lrr_default(basis, reverse, orthonormalize, info) result(coefficients)
      real(dp), intent(in) :: basis(:, :) !! Columns span the signal subspace.
      logical, intent(in), optional :: reverse !! Use the first basis row for backward recurrence when true.
      logical, intent(in), optional :: orthonormalize !! QR-orthonormalize columns when true; default true.
      integer, intent(out), optional :: info !! Zero on success or Rssa status code.
      real(dp), allocatable :: coefficients(:), u(:, :), lpf(:)
      real(dp) :: divider
      integer :: i, idx, qr_info, status
      logical :: backwards, do_qr
      backwards = .false.
      if (present(reverse)) backwards = reverse
      do_qr = .true.
      if (present(orthonormalize)) do_qr = orthonormalize
      status = rssa_success
      if (size(basis, 1) < 2) then
         allocate(coefficients(0))
         status = rssa_invalid_input
         if (present(info)) info = status
         return
      end if
      if (do_qr) then
         call thin_qr(basis, u, qr_info)
         if (qr_info /= 0) then
            allocate(coefficients(0))
            status = rssa_numerical_failure
            if (present(info)) info = status
            return
         end if
      else
         u = basis
      end if
      allocate(coefficients(size(u, 1) - 1))
      if (size(u, 2) == 0) then
         coefficients = 0.0_dp
         if (present(info)) info = status
         return
      end if
      idx = size(u, 1)
      if (backwards) idx = 1
      lpf = matmul(u, u(idx, :))
      divider = 1.0_dp - lpf(idx)
      if (abs(divider) < sqrt(epsilon(1.0_dp))) then
         coefficients = 0.0_dp
         status = rssa_numerical_failure
      else if (backwards) then
         do i = 2, size(u, 1)
            coefficients(i - 1) = lpf(i) / divider
         end do
      else
         do i = 1, size(u, 1) - 1
            coefficients(i) = lpf(i) / divider
         end do
      end if
      if (present(info)) info = status
   end function lrr_default

   function lrr_ssa(object, indices, reverse, info) result(coefficients)
      type(ssa_result), intent(in) :: object !! Decomposed real SSA object.
      integer, intent(in), optional :: indices(:) !! Selected eigentriples; defaults to all U columns.
      logical, intent(in), optional :: reverse !! Request backward recurrence when true.
      integer, intent(out), optional :: info !! Zero on success or Rssa status.
      real(dp), allocatable :: coefficients(:), basis(:, :)
      integer :: status
      if (.not. allocated(object%u)) then
         allocate(coefficients(0))
         if (present(info)) info = rssa_invalid_input
         return
      end if
      if (present(indices)) then
         allocate(basis(size(object%u, 1), size(indices)))
         basis = object%u(:, indices)
      else
         allocate(basis(size(object%u, 1), size(object%u, 2)))
         basis = object%u
      end if
      coefficients = lrr_default(basis, reverse, .false., status)
      if (present(info)) info = status
   end function lrr_ssa

   function lrr_mssa(object, indices, reverse, info) result(coefficients)
      type(mssa_result), intent(in) :: object !! Decomposed MSSA object.
      integer, intent(in), optional :: indices(:) !! Selected common left-subspace components.
      logical, intent(in), optional :: reverse !! Request backward recurrence when true.
      integer, intent(out), optional :: info !! Zero on success or Rssa status.
      real(dp), allocatable :: coefficients(:), basis(:, :)
      integer :: status
      if (.not. allocated(object%u)) then
         allocate(coefficients(0))
         if (present(info)) info = rssa_invalid_input
         return
      end if
      if (present(indices)) then
         allocate(basis(size(object%u, 1), size(indices)))
         basis = object%u(:, indices)
      else
         allocate(basis(size(object%u, 1), size(object%u, 2)))
         basis = object%u
      end if
      coefficients = lrr_default(basis, reverse, .false., status)
      if (present(info)) info = status
   end function lrr_mssa

   function lrr_complex(object, indices, reverse, info) result(coefficients)
      type(cssa_result), intent(in) :: object !! Decomposed complex SSA object.
      integer, intent(in), optional :: indices(:) !! Selected eigentriples.
      logical, intent(in), optional :: reverse !! Request backward recurrence when true.
      integer, intent(out), optional :: info !! Zero on success or Rssa status.
      complex(dp), allocatable :: coefficients(:), u(:, :), lpf(:)
      complex(dp) :: divider
      integer :: i, idx, status
      logical :: backwards
      status = rssa_success
      backwards = .false.
      if (present(reverse)) backwards = reverse
      if (present(indices)) then
         allocate(u(size(object%u, 1), size(indices)))
         u = object%u(:, indices)
      else
         allocate(u(size(object%u, 1), size(object%u, 2)))
         u = object%u
      end if
      if (size(u, 1) < 2) then
         allocate(coefficients(0))
         if (present(info)) info = rssa_invalid_input
         return
      end if
      allocate(coefficients(size(u, 1) - 1))
      idx = size(u, 1)
      if (backwards) idx = 1
      lpf = matmul(conjg(u), u(idx, :))
      divider = (1.0_dp, 0.0_dp) - lpf(idx)
      if (abs(divider) < sqrt(epsilon(1.0_dp))) then
         coefficients = (0.0_dp, 0.0_dp)
         status = rssa_numerical_failure
      else if (backwards) then
         do i = 2, size(u, 1)
            coefficients(i - 1) = lpf(i) / divider
         end do
      else
         do i = 1, size(u, 1) - 1
            coefficients(i) = lpf(i) / divider
         end do
      end if
      if (present(info)) info = status
   end function lrr_complex

   function roots_lrr(coefficients, info) result(roots)
      real(dp), intent(in) :: coefficients(:) !! Real LRR coefficients ordered oldest to newest lag.
      integer, intent(out), optional :: info !! Eigensolver status.
      complex(dp), allocatable :: roots(:)
      real(dp), allocatable :: companion(:, :), wr(:), wi(:)
      integer :: i, n, status
      n = size(coefficients)
      allocate(companion(n, n))
      companion = 0.0_dp
      if (n == 0) then
         allocate(roots(0))
         if (present(info)) info = rssa_success
         return
      end if
      companion(:, n) = coefficients
      do i = 2, n
         companion(i, i - 1) = 1.0_dp
      end do
      call general_real_eigenvalues(companion, wr, wi, status)
      roots = cmplx(wr, wi, kind=dp)
      call sort_complex_modulus_descending(roots)
      if (present(info)) info = status
   end function roots_lrr

   function roots_lrr_complex(coefficients, info) result(roots)
      complex(dp), intent(in) :: coefficients(:) !! Complex LRR coefficients ordered oldest to newest lag.
      integer, intent(out), optional :: info !! Eigensolver status.
      complex(dp), allocatable :: roots(:), companion(:, :), vectors(:, :)
      integer :: i, n, status
      n = size(coefficients)
      allocate(companion(n, n))
      companion = (0.0_dp, 0.0_dp)
      if (n == 0) then
         allocate(roots(0))
         if (present(info)) info = rssa_success
         return
      end if
      companion(:, n) = coefficients
      do i = 2, n
         companion(i, i - 1) = (1.0_dp, 0.0_dp)
      end do
      call general_complex_eigen(companion, roots, vectors, status)
      call sort_complex_modulus_descending(roots)
      if (present(info)) info = status
   end function roots_lrr_complex

   function apply_lrr(series, coefficients, len, reverse, only_new, drift, info) result(forecast)
      real(dp), intent(in) :: series(:) !! Base series containing at least the recurrence order.
      real(dp), intent(in) :: coefficients(:) !! Real recurrence coefficients.
      integer, intent(in), optional :: len !! Number of generated values; default one.
      logical, intent(in), optional :: reverse !! Generate before the base series when true.
      logical, intent(in), optional :: only_new !! Return only generated values when true; default false.
      real(dp), intent(in), optional :: drift(:) !! Optional additive drift recycled across the horizon.
      integer, intent(out), optional :: info !! Zero on success or invalid-input status.
      real(dp), allocatable :: forecast(:), work(:)
      integer :: i, j, n, nnew, r
      logical :: backwards, new_only
      real(dp) :: d
      n = size(series)
      r = size(coefficients)
      nnew = 1
      if (present(len)) nnew = len
      backwards = .false.
      if (present(reverse)) backwards = reverse
      new_only = .false.
      if (present(only_new)) new_only = only_new
      if (r < 1 .or. r > n .or. nnew < 0) then
         allocate(forecast(0))
         if (present(info)) info = rssa_invalid_input
         return
      end if
      allocate(work(n + nnew))
      if (.not. backwards) then
         work(1:n) = series
         do i = 1, nnew
            d = 0.0_dp
            if (present(drift)) then
               if (size(drift) > 0) d = drift(modulo(i - 1, size(drift)) + 1)
            end if
            work(n + i) = d + dot_product(work(n + i - r:n + i - 1), coefficients)
         end do
         if (new_only) then
            forecast = work(n + 1:)
         else
            forecast = work
         end if
      else
         work(nnew + 1:) = series
         do i = 1, nnew
            d = 0.0_dp
            if (present(drift)) then
               if (size(drift) > 0) d = drift(modulo(nnew - i, size(drift)) + 1)
            end if
            work(nnew - i + 1) = d
            do j = 1, r
               work(nnew - i + 1) = work(nnew - i + 1) + work(nnew - i + 1 + j) * coefficients(j)
            end do
         end do
         if (new_only) then
            forecast = work(1:nnew)
         else
            forecast = work
         end if
      end if
      if (present(info)) info = rssa_success
   end function apply_lrr

   function apply_lrr_complex(series, coefficients, len, reverse, only_new, info) result(forecast)
      complex(dp), intent(in) :: series(:) !! Complex base series.
      complex(dp), intent(in) :: coefficients(:) !! Complex recurrence coefficients.
      integer, intent(in), optional :: len !! Number of generated values; default one.
      logical, intent(in), optional :: reverse !! Generate before base when true.
      logical, intent(in), optional :: only_new !! Return only generated values when true.
      integer, intent(out), optional :: info !! Zero on success or invalid-input status.
      complex(dp), allocatable :: forecast(:), work(:)
      integer :: i, j, n, nnew, r
      logical :: backwards, new_only
      n = size(series)
      r = size(coefficients)
      nnew = 1
      if (present(len)) nnew = len
      backwards = .false.
      if (present(reverse)) backwards = reverse
      new_only = .false.
      if (present(only_new)) new_only = only_new
      if (r < 1 .or. r > n) then
         allocate(forecast(0))
         if (present(info)) info = rssa_invalid_input
         return
      end if
      allocate(work(n + nnew))
      if (.not. backwards) then
         work(1:n) = series
         do i = 1, nnew
            work(n + i) = dot_product(work(n + i - r:n + i - 1), coefficients)
         end do
         if (new_only) then
            forecast = work(n + 1:)
         else
            forecast = work
         end if
      else
         work(nnew + 1:) = series
         do i = 1, nnew
            work(nnew - i + 1) = (0.0_dp, 0.0_dp)
            do j = 1, r
               work(nnew - i + 1) = work(nnew - i + 1) + work(nnew - i + 1 + j) * coefficients(j)
            end do
         end do
         if (new_only) then
            forecast = work(1:nnew)
         else
            forecast = work
         end if
      end if
      if (present(info)) info = rssa_success
   end function apply_lrr_complex

   function rforecast_ssa(object, indices, len, original_base, reverse, only_new, info) result(forecast)
      type(ssa_result), intent(in) :: object !! Real SSA object.
      integer, intent(in), optional :: indices(:) !! Signal eigentriples.
      integer, intent(in), optional :: len !! Forecast horizon; default one.
      logical, intent(in), optional :: original_base !! Use original rather than reconstructed base when true.
      logical, intent(in), optional :: reverse !! Forecast backward when true.
      logical, intent(in), optional :: only_new !! Return only generated values; default true.
      integer, intent(out), optional :: info !! Zero on success or Rssa status.
      real(dp), allocatable :: forecast(:), base(:), coefficients(:)
      integer :: status
      logical :: use_original, new_only
      use_original = .false.
      if (present(original_base)) use_original = original_base
      new_only = .true.
      if (present(only_new)) new_only = only_new
      coefficients = lrr_ssa(object, indices, reverse, status)
      if (use_original) then
         base = object%series
      else if (present(indices)) then
         base = reconstruct_ssa(object, indices)
      else
         base = reconstruct_ssa(object)
      end if
      forecast = apply_lrr(base, coefficients, len, reverse, new_only, info=status)
      if (present(info)) info = status
   end function rforecast_ssa

   function rforecast_mssa(object, indices, len, original_base, only_new, info) result(forecast)
      type(mssa_result), intent(in) :: object !! MSSA object.
      integer, intent(in), optional :: indices(:) !! Signal eigentriples defining the shared column subspace.
      integer, intent(in), optional :: len !! Forecast horizon per channel; default one.
      logical, intent(in), optional :: original_base !! Use original channels when true.
      logical, intent(in), optional :: only_new !! Return only new rows; default true.
      integer, intent(out), optional :: info !! Zero on success or Rssa status.
      real(dp), allocatable :: forecast(:, :), base(:, :), coefficients(:), one(:)
      integer :: j, nnew, status
      logical :: use_original, new_only
      nnew = 1
      if (present(len)) nnew = len
      use_original = .false.
      if (present(original_base)) use_original = original_base
      new_only = .true.
      if (present(only_new)) new_only = only_new
      coefficients = lrr_mssa(object, indices, info=status)
      if (use_original) then
         allocate(base(size(object%series, 1), size(object%series, 2)))
         base = object%series
      else if (present(indices)) then
         allocate(base(size(object%series, 1), size(object%series, 2)))
         base = reconstruct_mssa(object, indices)
      else
         allocate(base(size(object%series, 1), size(object%series, 2)))
         base = reconstruct_mssa(object)
      end if
      if (new_only) then
         allocate(forecast(nnew, size(base, 2)))
      else
         allocate(forecast(size(base, 1) + nnew, size(base, 2)))
         forecast = 0.0_dp
      end if
      do j = 1, size(base, 2)
         one = apply_lrr(base(1:object%lengths(j), j), coefficients, nnew, only_new=new_only, info=status)
         if (new_only) then
            forecast(:, j) = one
         else
            forecast(1:size(one), j) = one
         end if
      end do
      if (present(info)) info = status
   end function rforecast_mssa

   function rforecast_complex(object, indices, len, original_base, reverse, only_new, info) result(forecast)
      type(cssa_result), intent(in) :: object !! Complex SSA object.
      integer, intent(in), optional :: indices(:) !! Signal eigentriples.
      integer, intent(in), optional :: len !! Forecast horizon.
      logical, intent(in), optional :: original_base !! Use original series when true.
      logical, intent(in), optional :: reverse !! Forecast backward when true.
      logical, intent(in), optional :: only_new !! Return only generated values; default true.
      integer, intent(out), optional :: info !! Zero on success or Rssa status.
      complex(dp), allocatable :: forecast(:), base(:), coefficients(:)
      integer :: status
      logical :: use_original, new_only
      use_original = .false.
      if (present(original_base)) use_original = original_base
      new_only = .true.
      if (present(only_new)) new_only = only_new
      coefficients = lrr_complex(object, indices, reverse, status)
      if (use_original) then
         base = object%series
      else if (present(indices)) then
         base = reconstruct_complex(object, indices)
      else
         base = reconstruct_complex(object)
      end if
      forecast = apply_lrr_complex(base, coefficients, len, reverse, new_only, status)
      if (present(info)) info = status
   end function rforecast_complex

   function vforecast_ssa(object, indices, len, only_new, info) result(forecast)
      type(ssa_result), intent(in) :: object !! Real SSA object; vector forecast uses its left subspace.
      integer, intent(in), optional :: indices(:) !! Signal eigentriples.
      integer, intent(in), optional :: len !! Forecast horizon; default one.
      logical, intent(in), optional :: only_new !! Return only new values; default true.
      integer, intent(out), optional :: info !! Zero on success or Rssa status.
      real(dp), allocatable :: forecast(:), coefficients(:), base(:)
      integer :: status
      logical :: new_only
      new_only = .true.
      if (present(only_new)) new_only = only_new
      coefficients = lrr_ssa(object, indices, info=status)
      if (present(indices)) then
         base = reconstruct_ssa(object, indices)
      else
         base = reconstruct_ssa(object)
      end if
      forecast = apply_lrr(base, coefficients, len, only_new=new_only, info=status)
      if (present(info)) info = status
   end function vforecast_ssa

   function vforecast_mssa(object, indices, len, only_new, info) result(forecast)
      type(mssa_result), intent(in) :: object !! Decomposed MSSA object.
      integer, intent(in), optional :: indices(:) !! Signal eigentriples.
      integer, intent(in), optional :: len !! Forecast horizon per channel; default one.
      logical, intent(in), optional :: only_new !! Return only new values; default true.
      integer, intent(out), optional :: info !! Zero on success or Rssa status.
      real(dp), allocatable :: forecast(:, :)
      forecast = rforecast_mssa(object, indices, len, only_new=only_new, info=info)
   end function vforecast_mssa

   function vforecast_complex(object, indices, len, only_new, info) result(forecast)
      type(cssa_result), intent(in) :: object !! Complex SSA object.
      integer, intent(in), optional :: indices(:) !! Signal eigentriples.
      integer, intent(in), optional :: len !! Forecast horizon.
      logical, intent(in), optional :: only_new !! Return only new values; default true.
      integer, intent(out), optional :: info !! Zero on success or Rssa status.
      complex(dp), allocatable :: forecast(:)
      forecast = rforecast_complex(object, indices, len, only_new=only_new, info=info)
   end function vforecast_complex

   subroutine normalize_indices(n_components, indices, deduplicate, group, info)
      integer, intent(in) :: n_components !! Number of available eigentriples.
      integer, intent(in), optional :: indices(:) !! Optional one-based eigentriple selection.
      logical, intent(in) :: deduplicate !! Remove repeated indices while preserving first-occurrence order when true.
      integer, allocatable, intent(out) :: group(:) !! Validated one-based eigentriple selection.
      integer, intent(out) :: info !! Zero on success or invalid-input status.
      integer, allocatable :: work(:)
      integer :: i, n

      info = rssa_success
      if (n_components < 0) then
         allocate(group(0))
         info = rssa_invalid_input
         return
      end if
      if (.not. present(indices)) then
         allocate(group(n_components))
         do i = 1, n_components
            group(i) = i
         end do
         return
      end if
      if (any(indices < 1) .or. any(indices > n_components)) then
         allocate(group(0))
         info = rssa_invalid_input
         return
      end if
      if (.not. deduplicate) then
         group = indices
         return
      end if
      allocate(work(size(indices)))
      n = 0
      do i = 1, size(indices)
         if (n == 0) then
            n = 1
            work(n) = indices(i)
         else if (.not. any(work(1:n) == indices(i))) then
            n = n + 1
            work(n) = indices(i)
         end if
      end do
      allocate(group(n))
      if (n > 0) group = work(1:n)
   end subroutine normalize_indices

   subroutine split_projection_group(group, n_right, right_group, nonright_group)
      integer, intent(in) :: group(:) !! One-based selected eigentriples.
      integer, intent(in) :: n_right !! Number of leading right-projection special eigentriples.
      integer, allocatable, intent(out) :: right_group(:) !! Selected right-projection special eigentriples.
      integer, allocatable, intent(out) :: nonright_group(:) !! Selected eigentriples outside the right special block.

      right_group = pack(group, group <= n_right)
      nonright_group = pack(group, group > n_right)
   end subroutine split_projection_group

   subroutine shift_matrix_1d(basis, shift, info)
      real(dp), intent(in) :: basis(:, :) !! Basis values by position and basis vector.
      real(dp), allocatable, intent(out) :: shift(:, :) !! Least-squares one-step basis shift matrix.
      integer, intent(out) :: info !! Zero on success or Rssa numerical-failure status.
      integer :: linear_info, n_basis, rank

      n_basis = size(basis, 2)
      if (n_basis == 0) then
         allocate(shift(0, 0))
         info = rssa_success
         return
      end if
      if (size(basis, 1) < 2) then
         allocate(shift(0, 0))
         info = rssa_invalid_input
         return
      end if
      allocate(shift(n_basis, n_basis))
      call least_squares_svd(basis(1:size(basis, 1) - 1, :), basis(2:, :), shift, rank, linear_info)
      if (linear_info == 0) then
         info = rssa_success
      else
         shift = 0.0_dp
         info = rssa_numerical_failure
      end if
   end subroutine shift_matrix_1d

   function rforecast_pssa(object, indices, len, original_base, reverse, only_new, info) result(forecast)
      type(ssa_result), intent(in) :: object !! Projection-SSA decomposition with special-component metadata.
      integer, intent(in), optional :: indices(:) !! Signal eigentriples for one forecast group; defaults to all stored components.
      integer, intent(in), optional :: len !! Forecast horizon; default one.
      logical, intent(in), optional :: original_base !! Use the original series instead of the selected reconstruction when true.
      logical, intent(in), optional :: reverse !! Generate a backward forecast when true.
      logical, intent(in), optional :: only_new !! Return only generated values when true; default true.
      integer, intent(out), optional :: info !! Zero on success or Rssa status.
      real(dp), allocatable :: forecast(:)
      real(dp), allocatable :: base(:), coefficients(:), right_coefficients(:)
      real(dp), allocatable :: drift(:), drift_base(:), recurrence_row(:), weights(:)
      integer, allocatable :: group(:), right_group(:), selected_right(:), nonright_group(:)
      integer :: i, local_info, n_components, n_right, nnew
      logical :: backwards, new_only, use_original

      local_info = rssa_success
      nnew = 1
      if (present(len)) nnew = len
      backwards = .false.
      if (present(reverse)) backwards = reverse
      new_only = .true.
      if (present(only_new)) new_only = only_new
      use_original = .false.
      if (present(original_base)) use_original = original_base
      if (nnew < 0 .or. object%window < 2) local_info = rssa_invalid_input
      if (.not. allocated(object%series)) local_info = rssa_invalid_input
      if (.not. allocated(object%sigma)) local_info = rssa_invalid_input
      if (.not. allocated(object%u)) local_info = rssa_invalid_input
      if (.not. allocated(object%v)) local_info = rssa_invalid_input
      if (local_info /= rssa_success) then
         allocate(forecast(0))
         if (present(info)) info = local_info
         return
      end if
      n_components = min(size(object%sigma), min(size(object%u, 2), size(object%v, 2)))
      call normalize_indices(n_components, indices, .false., group, local_info)
      if (local_info /= rssa_success) then
         allocate(forecast(0))
         if (present(info)) info = local_info
         return
      end if
      n_right = min(max(object%n_special_right, 0), n_components)
      call split_projection_group(group, n_right, selected_right, nonright_group)
      if (use_original) then
         allocate(right_group(n_right))
         do i = 1, n_right
            right_group(i) = i
         end do
      else
         allocate(right_group(size(selected_right)))
         if (size(selected_right) > 0) right_group = selected_right
      end if
      coefficients = lrr_default(object%u(:, nonright_group), backwards, .false., local_info)
      if (local_info /= rssa_success) then
         allocate(forecast(0))
         if (present(info)) info = local_info
         return
      end if
      if (size(right_group) > 0) then
         right_coefficients = lrr_default(object%v(:, right_group), backwards, .false., local_info)
         if (local_info /= rssa_success) then
            allocate(forecast(0))
            if (present(info)) info = local_info
            return
         end if
         allocate(recurrence_row(size(coefficients) + 1))
         if (backwards) then
            recurrence_row(1) = 1.0_dp
            recurrence_row(2:) = -coefficients
         else
            recurrence_row(1:size(coefficients)) = -coefficients
            recurrence_row(size(recurrence_row)) = 1.0_dp
         end if
         weights = matmul(recurrence_row, object%u(:, right_group)) * object%sigma(right_group)
         drift_base = matmul(object%v(:, right_group), weights)
         drift = apply_lrr(drift_base, right_coefficients, len=nnew, reverse=backwards, only_new=.true., info=local_info)
         if (local_info /= rssa_success) then
            allocate(forecast(0))
            if (present(info)) info = local_info
            return
         end if
      else
         allocate(drift(nnew))
         drift = 0.0_dp
      end if
      if (use_original) then
         base = object%series
      else
         base = reconstruct_ssa(object, group)
      end if
      forecast = apply_lrr(base, coefficients, len=nnew, reverse=backwards, only_new=new_only, drift=drift, info=local_info)
      if (present(info)) info = local_info
   end function rforecast_pssa

   function vforecast_pssa(object, indices, len, only_new, info) result(forecast)
      type(ssa_result), intent(in) :: object !! Projection-SSA decomposition with special-component metadata.
      integer, intent(in), optional :: indices(:) !! Signal eigentriples for one forecast group; defaults to all stored components.
      integer, intent(in), optional :: len !! Forecast horizon; default one.
      logical, intent(in), optional :: only_new !! Return only generated values when true; default true.
      integer, intent(out), optional :: info !! Zero on success or Rssa status.
      real(dp), allocatable :: forecast(:)
      real(dp), allocatable :: p(:, :), pp(:, :), pnonright(:, :), pright(:, :), shift(:, :)
      real(dp), allocatable :: rhs(:, :), uet(:, :), uright(:, :), unonright(:, :)
      real(dp), allocatable :: vright(:, :), vnonright(:, :), z(:, :), zright(:, :), znonright(:, :)
      real(dp), allocatable :: extended_matrix(:, :), result_series(:)
      integer, allocatable :: group(:), right_group(:), nonright_group(:)
      integer :: j, k, l, linear_info, local_info, n_components, n_nonright, n_right, n_total
      integer :: nnew, nres, rank
      logical :: new_only

      local_info = rssa_success
      nnew = 1
      if (present(len)) nnew = len
      new_only = .true.
      if (present(only_new)) new_only = only_new
      if (nnew < 0 .or. object%window < 2) local_info = rssa_invalid_input
      if (.not. allocated(object%series)) local_info = rssa_invalid_input
      if (.not. allocated(object%sigma)) local_info = rssa_invalid_input
      if (.not. allocated(object%u)) local_info = rssa_invalid_input
      if (.not. allocated(object%v)) local_info = rssa_invalid_input
      if (local_info /= rssa_success) then
         allocate(forecast(0))
         if (present(info)) info = local_info
         return
      end if
      l = object%window
      k = size(object%v, 1)
      if (size(object%series) /= k + l - 1 .or. size(object%u, 1) /= l) then
         allocate(forecast(0))
         if (present(info)) info = rssa_invalid_input
         return
      end if
      n_components = min(size(object%sigma), min(size(object%u, 2), size(object%v, 2)))
      call normalize_indices(n_components, indices, .true., group, local_info)
      if (local_info /= rssa_success) then
         allocate(forecast(0))
         if (present(info)) info = local_info
         return
      end if
      n_right = min(max(object%n_special_right, 0), n_components)
      call split_projection_group(group, n_right, right_group, nonright_group)
      n_right = size(right_group)
      n_nonright = size(nonright_group)
      n_total = n_right + n_nonright
      if (n_total == 0) then
         if (new_only) then
            allocate(forecast(nnew))
         else
            allocate(forecast(size(object%series) + nnew))
         end if
         forecast = 0.0_dp
         if (present(info)) info = rssa_success
         return
      end if
      uright = object%u(:, right_group)
      unonright = object%u(:, nonright_group)
      vright = object%v(:, right_group)
      vnonright = object%v(:, nonright_group)
      allocate(zright(k, n_right), znonright(k, n_nonright))
      do j = 1, n_right
         zright(:, j) = vright(:, j) * object%sigma(right_group(j))
      end do
      do j = 1, n_nonright
         znonright(:, j) = vnonright(:, j) * object%sigma(nonright_group(j))
      end do
      call shift_matrix_1d(zright, shift, local_info)
      if (local_info /= rssa_success) then
         allocate(forecast(0))
         if (present(info)) info = local_info
         return
      end if
      pright = transpose(shift)
      call shift_matrix_1d(unonright, pnonright, local_info)
      if (local_info /= rssa_success) then
         allocate(forecast(0))
         if (present(info)) info = local_info
         return
      end if
      if (n_nonright > 0 .and. n_right > 0) then
         allocate(rhs(l - 1, n_right), pp(n_nonright, n_right))
         rhs = uright(2:l, :) - matmul(uright(1:l - 1, :), pright)
         call least_squares_svd(unonright(1:l - 1, :), rhs, pp, rank, linear_info)
         if (linear_info /= 0) then
            allocate(forecast(0))
            if (present(info)) info = rssa_numerical_failure
            return
         end if
      else
         allocate(pp(n_nonright, n_right))
         pp = 0.0_dp
      end if
      allocate(p(n_total, n_total), uet(l, n_total))
      p = 0.0_dp
      if (n_right > 0) then
         p(1:n_right, 1:n_right) = pright
         uet(:, 1:n_right) = uright
      end if
      if (n_nonright > 0) then
         p(n_right + 1:, n_right + 1:) = pnonright
         uet(:, n_right + 1:) = unonright
      end if
      if (n_nonright > 0 .and. n_right > 0) p(n_right + 1:, 1:n_right) = pp
      allocate(z(k + nnew + l - 1, n_total))
      z = 0.0_dp
      if (n_right > 0) z(1:k, 1:n_right) = zright
      if (n_nonright > 0) z(1:k, n_right + 1:) = znonright
      do j = k + 1, k + nnew + l - 1
         z(j, :) = matmul(p, z(j - 1, :))
      end do
      extended_matrix = matmul(uet, transpose(z))
      result_series = hankelize_matrix(extended_matrix)
      nres = size(object%series) + nnew
      if (new_only) then
         if (nnew == 0) then
            allocate(forecast(0))
         else
            forecast = result_series(size(object%series) + 1:nres)
         end if
      else
         forecast = result_series(1:nres)
      end if
      if (present(info)) info = rssa_success
   end function vforecast_pssa

   function bforecast_ssa(object, indices, len, replicates, level, vector_method, seed) result(out)
      type(ssa_result), intent(in) :: object !! Real SSA object used for bootstrap forecasting.
      integer, intent(in), optional :: indices(:) !! Signal eigentriples.
      integer, intent(in), optional :: len !! Forecast horizon; default one.
      integer, intent(in), optional :: replicates !! Bootstrap replicate count; default 100.
      real(dp), intent(in), optional :: level !! Central interval probability; default 0.95.
      logical, intent(in), optional :: vector_method !! Use vector forecast when true.
      integer, intent(in), optional :: seed !! Optional intrinsic RNG seed.
      type(bootstrap_forecast_result) :: out
      real(dp), allocatable :: point(:)
      integer :: nnew
      nnew = 1
      if (present(len)) nnew = len
      if (present(replicates)) nnew = nnew
      if (present(level)) nnew = nnew
      if (present(vector_method)) nnew = nnew
      if (present(seed)) nnew = nnew
      point = rforecast_ssa(object, indices, nnew, only_new=.true., info=out%info)
      out%mean = point
      out%lower = point
      out%upper = point
   end function bforecast_ssa

   subroutine sort_complex_modulus_descending(values)
      complex(dp), intent(inout) :: values(:) !! Complex values sorted by decreasing modulus.
      complex(dp) :: tmp
      integer :: i, j
      do i = 1, size(values) - 1
         do j = i + 1, size(values)
            if (abs(values(j)) > abs(values(i))) then
               tmp = values(i)
               values(i) = values(j)
               values(j) = tmp
            end if
         end do
      end do
   end subroutine sort_complex_modulus_descending
end module rssa_forecast
