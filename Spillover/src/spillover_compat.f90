! SPDX-License-Identifier: GPL-2.0-only
! Modern Fortran translation of computational code from the Spillover R package.
module spillover_compat
   use spillover_kinds, only : dp
   use spillover_var, only : var_model
   use spillover_fevd, only : generalized_fevd
   use spillover_tables, only : spillover_result, generalized_spillover, &
      orthogonalized_spillover
   use spillover_dynamic, only : rolling_total_spillover, rolling_net_spillover
   implicit none
   private

   public :: g_fevd
   public :: g_spillover
   public :: o_spillover
   public :: roll_spillover
   public :: total_dynamic_spillover
   public :: roll_net

contains

   subroutine g_fevd(model, n_ahead, normalized, fevd, info, message)
      type(var_model), intent(in) :: model
      integer, intent(in) :: n_ahead
      logical, intent(in) :: normalized
      real(dp), allocatable, intent(out) :: fevd(:, :, :)
      integer, intent(out), optional :: info
      character(len=*), intent(out), optional :: message

      call generalized_fevd(model, n_ahead, normalized, fevd, info, message)
   end subroutine g_fevd

   subroutine g_spillover(model, n_ahead, standardized, result, info, message)
      type(var_model), intent(in) :: model
      integer, intent(in) :: n_ahead
      logical, intent(in) :: standardized
      type(spillover_result), intent(out) :: result
      integer, intent(out), optional :: info
      character(len=*), intent(out), optional :: message

      call generalized_spillover(model, n_ahead, standardized, result, info, message)
   end subroutine g_spillover

   subroutine o_spillover(model, n_ahead, ortho_type, standardized, result, &
      n_permutations, seed, source_compatible, info, message)
      type(var_model), intent(in) :: model
      integer, intent(in) :: n_ahead
      integer, intent(in) :: ortho_type
      logical, intent(in) :: standardized
      type(spillover_result), intent(out) :: result
      integer, intent(in), optional :: n_permutations
      integer, intent(in), optional :: seed
      logical, intent(in), optional :: source_compatible
      integer, intent(out), optional :: info
      character(len=*), intent(out), optional :: message

      call orthogonalized_spillover(model, n_ahead, ortho_type, standardized, result, &
         n_permutations=n_permutations, seed=seed, source_compatible=source_compatible, &
         info=info, message=message)
   end subroutine o_spillover

   subroutine roll_spillover(y, width, p, n_ahead, index_type, values, &
      ortho_type, n_permutations, seed, info, message)
      real(dp), intent(in) :: y(:, :)
      integer, intent(in) :: width, p, n_ahead, index_type
      real(dp), allocatable, intent(out) :: values(:)
      integer, intent(in), optional :: ortho_type, n_permutations, seed
      integer, intent(out), optional :: info
      character(len=*), intent(out), optional :: message

      call rolling_total_spillover(y, width, p, n_ahead, index_type, values, &
         ortho_type=ortho_type, n_permutations=n_permutations, seed=seed, &
         info=info, message=message)
   end subroutine roll_spillover

   subroutine total_dynamic_spillover(y, width, p, n_ahead, index_type, values, &
      ortho_type, n_permutations, seed, info, message)
      real(dp), intent(in) :: y(:, :)
      integer, intent(in) :: width, p, n_ahead, index_type
      real(dp), allocatable, intent(out) :: values(:)
      integer, intent(in), optional :: ortho_type, n_permutations, seed
      integer, intent(out), optional :: info
      character(len=*), intent(out), optional :: message

      call rolling_total_spillover(y, width, p, n_ahead, index_type, values, &
         ortho_type=ortho_type, n_permutations=n_permutations, seed=seed, &
         info=info, message=message)
   end subroutine total_dynamic_spillover

   subroutine roll_net(y, width, p, n_ahead, index_type, values, &
      ortho_type, n_permutations, seed, info, message)
      real(dp), intent(in) :: y(:, :)
      integer, intent(in) :: width, p, n_ahead, index_type
      real(dp), allocatable, intent(out) :: values(:, :)
      integer, intent(in), optional :: ortho_type, n_permutations, seed
      integer, intent(out), optional :: info
      character(len=*), intent(out), optional :: message

      call rolling_net_spillover(y, width, p, n_ahead, index_type, values, &
         ortho_type=ortho_type, n_permutations=n_permutations, seed=seed, &
         info=info, message=message)
   end subroutine roll_net

end module spillover_compat
