! SPDX-License-Identifier: GPL-3.0-only
module neighbours_nmof_adapter
   use neighbours, only: numeric_neighbour_config, numeric_neighbour, &
      neighbours_rng_state => rng_state
   use nmof_kinds, only: dp
   use nmof_rng, only: nmof_rng_state => rng_state
   implicit none
   private
   public :: nmof_numeric_neighbour
contains
   subroutine nmof_numeric_neighbour(x, xn, rng, iteration, total_iterations, context)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: xn(:)
      type(nmof_rng_state), intent(inout) :: rng
      integer, intent(in) :: iteration, total_iterations
      class(*), intent(in), optional :: context
      type(neighbours_rng_state) :: local_rng
      integer :: status

      call consume_iterations(iteration, total_iterations)
      local_rng%state = rng%state
      if (present(context)) then
         select type (config => context)
         type is (numeric_neighbour_config)
            call numeric_neighbour(config, x, xn, local_rng, status=status)
         class default
            xn = x
         end select
      else
         xn = x
      end if
      rng%state = local_rng%state
   end subroutine nmof_numeric_neighbour

   subroutine consume_iterations(iteration, total_iterations)
      integer, intent(in) :: iteration, total_iterations
      if (iteration < 0 .and. total_iterations < 0) error stop 1
   end subroutine consume_iterations
end module neighbours_nmof_adapter
