! SPDX-License-Identifier: GPL-2.0-only
! Modern Fortran translation of computational code from the Spillover R package.
module spillover_dynamic
   use spillover_kinds, only : dp
   use spillover_status, only : spillover_success, spillover_invalid_argument, set_status
   use spillover_var, only : var_model, fit_var, var_const
   use spillover_tables, only : spillover_result, generalized_spillover, &
      orthogonalized_spillover, ortho_single
   implicit none
   private

   integer, parameter, public :: index_generalized = 1
   integer, parameter, public :: index_orthogonalized = 2

   type, public :: dynamic_spillover_result
      integer :: n_windows = 0
      integer :: k = 0
      integer :: n_pairs = 0
      integer, allocatable :: endpoint(:)
      integer, allocatable :: pair_i(:)
      integer, allocatable :: pair_j(:)
      real(dp), allocatable :: tables(:, :, :)
      real(dp), allocatable :: from(:, :)
      real(dp), allocatable :: to(:, :)
      real(dp), allocatable :: net(:, :)
      real(dp), allocatable :: total(:)
      real(dp), allocatable :: net_pairwise(:, :)
   end type dynamic_spillover_result

   public :: dynamic_spillover
   public :: rolling_total_spillover
   public :: rolling_net_spillover

contains

   subroutine dynamic_spillover(y, width, p, horizon, result, deterministic, standardized, &
      source_compatible_direction, info, message)
      real(dp), intent(in) :: y(:, :)
      integer, intent(in) :: width
      integer, intent(in) :: p
      integer, intent(in) :: horizon
      type(dynamic_spillover_result), intent(out) :: result
      integer, intent(in), optional :: deterministic
      logical, intent(in), optional :: standardized
      logical, intent(in), optional :: source_compatible_direction
      integer, intent(out), optional :: info
      character(len=*), intent(out), optional :: message

      type(var_model) :: model
      type(spillover_result) :: spill
      integer :: n, k, nwin, npairs, w, first, det, local_info, i, j, pair
      logical :: std, source_direction
      character(len=160) :: local_message

      call clear_dynamic(result)
      call set_status(info, message, spillover_success, 'success')
      n = size(y, 1)
      k = size(y, 2)
      if (width < 2 .or. width > n .or. p < 1 .or. horizon < 1 .or. k < 1) then
         call set_status(info, message, spillover_invalid_argument, &
            'invalid rolling width, VAR order, horizon, or data dimension')
         return
      end if
      det = var_const
      if (present(deterministic)) det = deterministic
      std = .true.
      if (present(standardized)) std = standardized
      source_direction = .false.
      if (present(source_compatible_direction)) source_direction = source_compatible_direction

      nwin = n - width + 1
      npairs = k * (k - 1) / 2
      result%n_windows = nwin
      result%k = k
      result%n_pairs = npairs
      allocate(result%endpoint(nwin), result%tables(k, k, nwin))
      allocate(result%from(nwin, k), result%to(nwin, k), result%net(nwin, k))
      allocate(result%total(nwin), result%pair_i(npairs), result%pair_j(npairs))
      allocate(result%net_pairwise(nwin, npairs))

      pair = 0
      do i = 1, k - 1
         do j = i + 1, k
            pair = pair + 1
            result%pair_i(pair) = i
            result%pair_j(pair) = j
         end do
      end do

      do w = 1, nwin
         first = w
         result%endpoint(w) = first + width - 1
         call fit_var(y(first:first + width - 1, :), p, model, det, local_info, local_message)
         if (local_info /= spillover_success) then
            call clear_dynamic(result)
            call set_status(info, message, local_info, &
               'rolling VAR fit failed: ' // trim(local_message))
            return
         end if
         call generalized_spillover(model, horizon, std, spill, local_info, local_message)
         if (local_info /= spillover_success) then
            call clear_dynamic(result)
            call set_status(info, message, local_info, &
               'rolling generalized spillover failed: ' // trim(local_message))
            return
         end if

         result%tables(:, :, w) = spill%shares
         result%total(w) = spill%total
         if (source_direction) then
            result%from(w, :) = spill%to
            result%to(w, :) = spill%from
            result%net(w, :) = spill%to - spill%from
         else
            result%from(w, :) = spill%from
            result%to(w, :) = spill%to
            result%net(w, :) = spill%net
         end if
         do pair = 1, npairs
            i = result%pair_i(pair)
            j = result%pair_j(pair)
            result%net_pairwise(w, pair) = spill%shares(j, i) - spill%shares(i, j)
         end do
      end do
   end subroutine dynamic_spillover

   subroutine rolling_total_spillover(y, width, p, horizon, index_type, values, &
      ortho_type, deterministic, n_permutations, seed, source_compatible_ortho, info, message)
      real(dp), intent(in) :: y(:, :)
      integer, intent(in) :: width
      integer, intent(in) :: p
      integer, intent(in) :: horizon
      integer, intent(in) :: index_type
      real(dp), allocatable, intent(out) :: values(:)
      integer, intent(in), optional :: ortho_type
      integer, intent(in), optional :: deterministic
      integer, intent(in), optional :: n_permutations
      integer, intent(in), optional :: seed
      logical, intent(in), optional :: source_compatible_ortho
      integer, intent(out), optional :: info
      character(len=*), intent(out), optional :: message

      type(var_model) :: model
      type(spillover_result) :: spill
      integer :: n, nwin, w, first, det, otype, nperm, local_seed, local_info
      logical :: source_mode
      character(len=160) :: local_message

      call set_status(info, message, spillover_success, 'success')
      n = size(y, 1)
      if (width < 2 .or. width > n .or. p < 1 .or. horizon < 1) then
         allocate(values(0))
         call set_status(info, message, spillover_invalid_argument, &
            'invalid rolling arguments')
         return
      end if
      det = var_const
      if (present(deterministic)) det = deterministic
      otype = ortho_single
      if (present(ortho_type)) otype = ortho_type
      nperm = 10000
      if (present(n_permutations)) nperm = n_permutations
      local_seed = 12345
      if (present(seed)) local_seed = seed
      source_mode = .false.
      if (present(source_compatible_ortho)) source_mode = source_compatible_ortho

      nwin = n - width + 1
      allocate(values(nwin))
      do w = 1, nwin
         first = w
         call fit_var(y(first:first + width - 1, :), p, model, det, local_info, local_message)
         if (local_info /= spillover_success) then
            values = 0.0_dp
            call set_status(info, message, local_info, &
               'rolling VAR fit failed: ' // trim(local_message))
            return
         end if
         select case (index_type)
         case (index_generalized)
            call generalized_spillover(model, horizon, .true., spill, local_info, local_message)
         case (index_orthogonalized)
            call orthogonalized_spillover(model, horizon, otype, .true., spill, &
               n_permutations=nperm, seed=local_seed + w - 1, &
               source_compatible=source_mode, info=local_info, message=local_message)
         case default
            values = 0.0_dp
            call set_status(info, message, spillover_invalid_argument, 'invalid index type')
            return
         end select
         if (local_info /= spillover_success) then
            values = 0.0_dp
            call set_status(info, message, local_info, &
               'rolling spillover calculation failed: ' // trim(local_message))
            return
         end if
         values(w) = spill%total
      end do
   end subroutine rolling_total_spillover

   subroutine rolling_net_spillover(y, width, p, horizon, index_type, values, &
      ortho_type, deterministic, n_permutations, seed, source_compatible_ortho, info, message)
      real(dp), intent(in) :: y(:, :)
      integer, intent(in) :: width
      integer, intent(in) :: p
      integer, intent(in) :: horizon
      integer, intent(in) :: index_type
      real(dp), allocatable, intent(out) :: values(:, :)
      integer, intent(in), optional :: ortho_type
      integer, intent(in), optional :: deterministic
      integer, intent(in), optional :: n_permutations
      integer, intent(in), optional :: seed
      logical, intent(in), optional :: source_compatible_ortho
      integer, intent(out), optional :: info
      character(len=*), intent(out), optional :: message

      type(var_model) :: model
      type(spillover_result) :: spill
      integer :: n, k, nwin, w, first, det, otype, nperm, local_seed, local_info
      logical :: source_mode
      character(len=160) :: local_message

      call set_status(info, message, spillover_success, 'success')
      n = size(y, 1)
      k = size(y, 2)
      if (width < 2 .or. width > n .or. p < 1 .or. horizon < 1 .or. k < 1) then
         allocate(values(0, 0))
         call set_status(info, message, spillover_invalid_argument, &
            'invalid rolling arguments')
         return
      end if
      det = var_const
      if (present(deterministic)) det = deterministic
      otype = ortho_single
      if (present(ortho_type)) otype = ortho_type
      nperm = 10000
      if (present(n_permutations)) nperm = n_permutations
      local_seed = 12345
      if (present(seed)) local_seed = seed
      source_mode = .false.
      if (present(source_compatible_ortho)) source_mode = source_compatible_ortho

      nwin = n - width + 1
      allocate(values(nwin, k))
      do w = 1, nwin
         first = w
         call fit_var(y(first:first + width - 1, :), p, model, det, local_info, local_message)
         if (local_info /= spillover_success) then
            values = 0.0_dp
            call set_status(info, message, local_info, &
               'rolling VAR fit failed: ' // trim(local_message))
            return
         end if
         select case (index_type)
         case (index_generalized)
            call generalized_spillover(model, horizon, .true., spill, local_info, local_message)
         case (index_orthogonalized)
            call orthogonalized_spillover(model, horizon, otype, .true., spill, &
               n_permutations=nperm, seed=local_seed + w - 1, &
               source_compatible=source_mode, info=local_info, message=local_message)
         case default
            values = 0.0_dp
            call set_status(info, message, spillover_invalid_argument, 'invalid index type')
            return
         end select
         if (local_info /= spillover_success) then
            values = 0.0_dp
            call set_status(info, message, local_info, &
               'rolling spillover calculation failed: ' // trim(local_message))
            return
         end if
         values(w, :) = spill%net
      end do
   end subroutine rolling_net_spillover

   subroutine clear_dynamic(result)
      type(dynamic_spillover_result), intent(out) :: result

      result%n_windows = 0
      result%k = 0
      result%n_pairs = 0
   end subroutine clear_dynamic

end module spillover_dynamic
