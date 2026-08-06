! SPDX-License-Identifier: GPL-2.0-only
! Modern Fortran translation of computational code from the Spillover R package.
module spillover_tables
   use spillover_kinds, only : dp
   use spillover_status, only : spillover_success, spillover_invalid_argument, &
      spillover_iteration_limit, set_status
   use spillover_random, only : spillover_rng
   use spillover_var, only : var_model
   use spillover_fevd, only : generalized_fevd, orthogonalized_fevd
   implicit none
   private

   integer, parameter, public :: ortho_single = 1
   integer, parameter, public :: ortho_partial = 2
   integer, parameter, public :: ortho_total = 3

   type, public :: spillover_result
      real(dp), allocatable :: shares(:, :)
      real(dp), allocatable :: from(:)
      real(dp), allocatable :: to(:)
      real(dp), allocatable :: net(:)
      real(dp) :: total = 0.0_dp
      logical :: standardized = .true.
   end type spillover_result

   type, public :: orthogonal_average_result
      real(dp), allocatable :: average(:, :)
      real(dp), allocatable :: minimum(:, :)
      real(dp), allocatable :: maximum(:, :)
      integer :: n_permutations = 0
   end type orthogonal_average_result

   public :: generalized_spillover
   public :: orthogonalized_spillover
   public :: orthogonal_average_exact
   public :: orthogonal_average_sample
   public :: spillover_from_shares
   public :: compatibility_table

contains

   subroutine generalized_spillover(model, horizon, standardized, result, info, message)
      type(var_model), intent(in) :: model
      integer, intent(in) :: horizon
      logical, intent(in) :: standardized
      type(spillover_result), intent(out) :: result
      integer, intent(out), optional :: info
      character(len=*), intent(out), optional :: message

      real(dp), allocatable :: fevd(:, :, :)
      integer :: local_info
      character(len=160) :: local_message

      call generalized_fevd(model, horizon, .true., fevd, local_info, local_message)
      if (local_info /= spillover_success) then
         call empty_result(result)
         call set_status(info, message, local_info, trim(local_message))
         return
      end if
      call spillover_from_shares(100.0_dp * fevd(:, :, horizon), standardized, result, info, message)
   end subroutine generalized_spillover

   subroutine orthogonalized_spillover(model, horizon, ortho_type, standardized, result, &
      n_permutations, seed, source_compatible, exact_limit, info, message)
      type(var_model), intent(in) :: model
      integer, intent(in) :: horizon
      integer, intent(in) :: ortho_type
      logical, intent(in) :: standardized
      type(spillover_result), intent(out) :: result
      integer, intent(in), optional :: n_permutations
      integer, intent(in), optional :: seed
      logical, intent(in), optional :: source_compatible
      integer, intent(in), optional :: exact_limit
      integer, intent(out), optional :: info
      character(len=*), intent(out), optional :: message

      real(dp), allocatable :: fevd(:, :, :), shares(:, :)
      type(orthogonal_average_result) :: average_result
      integer :: local_info, nperm, local_seed, limit, effective_type
      logical :: source_mode
      character(len=160) :: local_message

      call set_status(info, message, spillover_success, 'success')
      source_mode = .false.
      if (present(source_compatible)) source_mode = source_compatible
      effective_type = ortho_type
      if (source_mode) then
         if (ortho_type == ortho_partial) effective_type = ortho_total
         if (ortho_type == ortho_total) effective_type = ortho_partial
      end if

      select case (effective_type)
      case (ortho_single)
         call orthogonalized_fevd(model, horizon, fevd, info=local_info, message=local_message)
         if (local_info == spillover_success) shares = 100.0_dp * fevd(:, :, horizon)
      case (ortho_partial)
         nperm = 10000
         if (present(n_permutations)) nperm = n_permutations
         local_seed = 12345
         if (present(seed)) local_seed = seed
         call orthogonal_average_sample(model, horizon, nperm, local_seed, average_result, &
            local_info, local_message)
         if (local_info == spillover_success) shares = average_result%average
      case (ortho_total)
         limit = 9
         if (present(exact_limit)) limit = exact_limit
         call orthogonal_average_exact(model, horizon, average_result, limit, &
            local_info, local_message)
         if (local_info == spillover_success) shares = average_result%average
      case default
         call empty_result(result)
         call set_status(info, message, spillover_invalid_argument, 'invalid orthogonalization type')
         return
      end select

      if (local_info /= spillover_success) then
         call empty_result(result)
         call set_status(info, message, local_info, trim(local_message))
         return
      end if
      call spillover_from_shares(shares, standardized, result, info, message)
   end subroutine orthogonalized_spillover

   subroutine orthogonal_average_exact(model, horizon, result, exact_limit, info, message)
      type(var_model), intent(in) :: model
      integer, intent(in) :: horizon
      type(orthogonal_average_result), intent(out) :: result
      integer, intent(in), optional :: exact_limit
      integer, intent(out), optional :: info
      character(len=*), intent(out), optional :: message

      integer, allocatable :: perm(:)
      real(dp), allocatable :: fevd(:, :, :), table(:, :)
      integer :: k, count, local_info, limit
      logical :: more
      character(len=160) :: local_message

      call clear_average(result)
      call set_status(info, message, spillover_success, 'success')
      k = model%k
      limit = 9
      if (present(exact_limit)) limit = exact_limit
      if (k < 1 .or. horizon < 1) then
         call set_status(info, message, spillover_invalid_argument, &
            'invalid model or horizon for exact ordering average')
         return
      end if
      if (k > limit) then
         call set_status(info, message, spillover_iteration_limit, &
            'exact permutation average exceeds the configured dimension limit')
         return
      end if

      allocate(perm(k), result%average(k, k), result%minimum(k, k), result%maximum(k, k))
      perm = [(count, count = 1, k)]
      result%average = 0.0_dp
      result%minimum = huge(1.0_dp)
      result%maximum = -huge(1.0_dp)
      count = 0
      more = .true.
      do while (more)
         call orthogonalized_fevd(model, horizon, fevd, perm, local_info, local_message)
         if (local_info /= spillover_success) then
            call clear_average(result)
            call set_status(info, message, local_info, trim(local_message))
            return
         end if
         table = 100.0_dp * fevd(:, :, horizon)
         result%average = result%average + table
         result%minimum = min(result%minimum, table)
         result%maximum = max(result%maximum, table)
         count = count + 1
         more = next_permutation(perm)
      end do
      result%average = result%average / real(count, dp)
      result%n_permutations = count
   end subroutine orthogonal_average_exact

   subroutine orthogonal_average_sample(model, horizon, n_permutations, seed, result, info, message)
      type(var_model), intent(in) :: model
      integer, intent(in) :: horizon
      integer, intent(in) :: n_permutations
      integer, intent(in) :: seed
      type(orthogonal_average_result), intent(out) :: result
      integer, intent(out), optional :: info
      character(len=*), intent(out), optional :: message

      type(spillover_rng) :: rng
      integer, allocatable :: perm(:)
      real(dp), allocatable :: fevd(:, :, :), table(:, :)
      integer :: k, draw, local_info
      character(len=160) :: local_message

      call clear_average(result)
      call set_status(info, message, spillover_success, 'success')
      k = model%k
      if (k < 1 .or. horizon < 1 .or. n_permutations < 1) then
         call set_status(info, message, spillover_invalid_argument, &
            'invalid model, horizon, or number of sampled permutations')
         return
      end if

      allocate(perm(k), result%average(k, k), result%minimum(k, k), result%maximum(k, k))
      result%average = 0.0_dp
      result%minimum = huge(1.0_dp)
      result%maximum = -huge(1.0_dp)
      call rng%seed(seed)
      do draw = 1, n_permutations
         call rng%permutation(perm)
         call orthogonalized_fevd(model, horizon, fevd, perm, local_info, local_message)
         if (local_info /= spillover_success) then
            call clear_average(result)
            call set_status(info, message, local_info, trim(local_message))
            return
         end if
         table = 100.0_dp * fevd(:, :, horizon)
         result%average = result%average + table
         result%minimum = min(result%minimum, table)
         result%maximum = max(result%maximum, table)
      end do
      result%average = result%average / real(n_permutations, dp)
      result%n_permutations = n_permutations
   end subroutine orthogonal_average_sample

   subroutine spillover_from_shares(shares, standardized, result, info, message)
      real(dp), intent(in) :: shares(:, :)
      logical, intent(in) :: standardized
      type(spillover_result), intent(out) :: result
      integer, intent(out), optional :: info
      character(len=*), intent(out), optional :: message

      integer :: k, i

      call empty_result(result)
      call set_status(info, message, spillover_success, 'success')
      k = size(shares, 1)
      if (k < 1 .or. size(shares, 2) /= k) then
         call set_status(info, message, spillover_invalid_argument, &
            'spillover shares must be a nonempty square matrix')
         return
      end if
      if (any(shares < -100.0_dp * epsilon(1.0_dp))) then
         call set_status(info, message, spillover_invalid_argument, &
            'spillover shares must be nonnegative')
         return
      end if

      allocate(result%shares(k, k), result%from(k), result%to(k), result%net(k))
      result%shares = shares
      result%standardized = standardized
      do i = 1, k
         result%from(i) = sum(shares(i, :)) - shares(i, i)
         result%to(i) = sum(shares(:, i)) - shares(i, i)
      end do
      result%net = result%to - result%from
      result%total = sum(result%from) / real(k, dp)
   end subroutine spillover_from_shares

   subroutine compatibility_table(result, table, info, message)
      type(spillover_result), intent(in) :: result
      real(dp), allocatable, intent(out) :: table(:, :)
      integer, intent(out), optional :: info
      character(len=*), intent(out), optional :: message

      integer :: k, i

      call set_status(info, message, spillover_success, 'success')
      if (.not. allocated(result%shares)) then
         allocate(table(0, 0))
         call set_status(info, message, spillover_invalid_argument, &
            'spillover result is uninitialized')
         return
      end if
      k = size(result%shares, 1)
      allocate(table(k + 2, k + 1))
      table = 0.0_dp
      table(1:k, 1:k) = result%shares
      table(k + 1, 1:k) = result%to
      do i = 1, k
         table(k + 2, i) = sum(result%shares(:, i))
      end do
      if (result%standardized) then
         table(1:k, k + 1) = result%from / real(k, dp)
         table(k + 1, k + 1) = result%total
         table(k + 2, k + 1) = sum(table(k + 2, 1:k)) / real(k, dp)
      else
         table(1:k, k + 1) = result%from
         table(k + 1, k + 1) = result%total
         table(k + 2, k + 1) = sum(table(k + 2, 1:k))
      end if
   end subroutine compatibility_table

   logical function next_permutation(values)
      integer, intent(inout) :: values(:)
      integer :: i, j, left, right, tmp

      i = size(values) - 1
      do while (i >= 1)
         if (values(i) < values(i + 1)) exit
         i = i - 1
      end do
      if (i < 1) then
         next_permutation = .false.
         return
      end if

      j = size(values)
      do while (values(j) <= values(i))
         j = j - 1
      end do
      tmp = values(i)
      values(i) = values(j)
      values(j) = tmp

      left = i + 1
      right = size(values)
      do while (left < right)
         tmp = values(left)
         values(left) = values(right)
         values(right) = tmp
         left = left + 1
         right = right - 1
      end do
      next_permutation = .true.
   end function next_permutation

   subroutine empty_result(result)
      type(spillover_result), intent(out) :: result

      result%total = 0.0_dp
      result%standardized = .true.
   end subroutine empty_result

   subroutine clear_average(result)
      type(orthogonal_average_result), intent(out) :: result

      result%n_permutations = 0
   end subroutine clear_average

end module spillover_tables
