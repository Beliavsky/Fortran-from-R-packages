! SPDX-License-Identifier: Artistic-2.0
module ldhmm_simulation
   use ldhmm_distribution, only : ecld_random
   use ldhmm_kinds, only : dp
   use ldhmm_math, only : categorical_random, correlation
   use ldhmm_modeling, only : ldhmm_state_distribution
   use ldhmm_status, only : LDHMM_SUCCESS, LDHMM_INVALID_ARGUMENT
   use ldhmm_types, only : ecld_type, ldhmm_model
   implicit none
   private

   public :: ldhmm_simulate_state_transition, ldhmm_simulate_abs_acf

contains

   function ldhmm_simulate_state_transition(input_model, init, status) result(model)
      type(ldhmm_model), intent(in) :: input_model
      integer, intent(in), optional :: init
      integer, intent(out), optional :: status
      type(ldhmm_model) :: model
      type(ecld_type) :: distribution
      integer :: i, n, state

      model = input_model
      call clear_simulation_outputs(model)
      if (present(init)) then
         n = init
         if (n < 1 .or. .not. allocated(input_model%delta)) then
            if (present(status)) status = LDHMM_INVALID_ARGUMENT
            return
         end if
         allocate(model%states_local(n), model%observations(n))
         do i = 1, n
            state = categorical_random(input_model%delta)
            model%states_local(i) = state
            distribution = ldhmm_state_distribution(input_model, state)
            model%observations(i) = ecld_random(distribution)
         end do
      else
         if (.not. allocated(input_model%states_local)) then
            if (present(status)) status = LDHMM_INVALID_ARGUMENT
            return
         end if
         n = size(input_model%states_local)
         allocate(model%states_local(n), model%observations(n))
         do i = 1, n
            state = input_model%states_local(i)
            if (state < 1 .or. state > input_model%m) state = 1
            state = categorical_random(input_model%gamma(state,:))
            model%states_local(i) = state
            distribution = ldhmm_state_distribution(input_model, state)
            model%observations(i) = ecld_random(distribution)
         end do
      end if
      if (present(status)) status = LDHMM_SUCCESS
   end function ldhmm_simulate_state_transition

   function ldhmm_simulate_abs_acf(model, n, lag_max, status) result(acf)
      type(ldhmm_model), intent(in) :: model
      integer, intent(in), optional :: n, lag_max
      integer, intent(out), optional :: status
      real(dp), allocatable :: acf(:)
      type(ldhmm_model) :: initial, advanced
      integer :: sample_size, maximum_lag, k, local_status

      sample_size = 10000
      maximum_lag = 5
      if (present(n)) sample_size = n
      if (present(lag_max)) maximum_lag = lag_max
      if (sample_size < 2 .or. maximum_lag < 1) then
         allocate(acf(0))
         if (present(status)) status = LDHMM_INVALID_ARGUMENT
         return
      end if
      initial = ldhmm_simulate_state_transition(model, init=sample_size, status=local_status)
      advanced = initial
      allocate(acf(maximum_lag))
      do k = 1, maximum_lag
         advanced = ldhmm_simulate_state_transition(advanced, status=local_status)
         acf(k) = correlation(abs(initial%observations), abs(advanced%observations))
      end do
      if (present(status)) status = local_status
   end function ldhmm_simulate_abs_acf

   subroutine clear_simulation_outputs(model)
      type(ldhmm_model), intent(inout) :: model
      if (allocated(model%observations)) deallocate(model%observations)
      if (allocated(model%states_prob)) deallocate(model%states_prob)
      if (allocated(model%states_local)) deallocate(model%states_local)
      if (allocated(model%states_global)) deallocate(model%states_global)
      if (allocated(model%states_local_stats)) deallocate(model%states_local_stats)
      if (allocated(model%states_global_stats)) deallocate(model%states_global_stats)
   end subroutine clear_simulation_outputs

end module ldhmm_simulation
