! This file is part of multiAssetOptions-fortran.
! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module mao_pricing
   use mao_kinds, only: dp
   use mao_status, only: status_type, clear_status, set_status, &
      mao_allocation_error, mao_solver_failure, mao_step_failure
   use mao_types, only: pricing_config, pricing_result, exercise_american, &
      timestep_constant, validate_config
   use mao_grid, only: build_grid
   use mao_payoff, only: payoff_values
   use mao_operator, only: build_fdm_operator
   use mao_sparse, only: csr_matrix, csr_matvec, csr_shifted_identity, &
      bicgstab_solve
   implicit none
   private

   public :: price_multi_asset

contains

   subroutine price_multi_asset(config, result, status, solver_tolerance, &
      solver_max_iterations)
      type(pricing_config), intent(in) :: config
      type(pricing_result), intent(out) :: result
      type(status_type), intent(out) :: status
      real(dp), intent(in), optional :: solver_tolerance
      integer, intent(in), optional :: solver_max_iterations

      type(csr_matrix) :: operator, system
      type(status_type) :: local_status
      real(dp), allocatable :: payoff(:), current(:), previous(:), candidate(:)
      real(dp), allocatable :: rhs(:), rhs_penalty(:), operator_value(:)
      real(dp), allocatable :: penalty(:)
      logical, allocatable :: active_old(:), active_new(:)
      real(dp), allocatable :: history(:,:), time_history(:)
      integer, allocatable :: linear_history(:), penalty_history(:)
      real(dp) :: remaining, dt, current_dt, next_dt, theta_step
      real(dp) :: large, relative_change, max_relative_change
      real(dp) :: linear_tolerance, residual_norm
      integer :: nnode, capacity, count, smooth_index, penalty_iter
      integer :: total_linear_iterations, linear_iterations, max_linear_iterations
      integer :: stat, max_steps

      call clear_status(status)
      call validate_config(config, local_status)
      if (.not. local_status%ok()) then
         status = local_status
         return
      end if

      call build_grid(config,result%grid,local_status)
      if (.not. local_status%ok()) then
         status = local_status
         return
      end if
      call payoff_values(config%opt%pay_type,config%opt%pc_flag, &
         config%opt%strike,result%grid,payoff,local_status)
      if (.not. local_status%ok()) then
         status = local_status
         return
      end if
      call build_fdm_operator(result%grid,config%opt%rf,config%opt%q, &
         config%opt%vol,config%opt%rho,operator,local_status)
      if (.not. local_status%ok()) then
         status = local_status
         return
      end if

      nnode = result%grid%n_nodes
      allocate(current(nnode), previous(nnode), candidate(nnode), rhs(nnode), &
         rhs_penalty(nnode), operator_value(nnode), penalty(nnode), &
         active_old(nnode), active_new(nnode), stat=stat)
      if (stat /= 0) then
         call set_status(status, mao_allocation_error, &
            'unable to allocate pricing workspace')
         return
      end if

      if (config%time%ts_type == timestep_constant) then
         capacity = config%time%n_steps + 1
         dt = config%opt%ttm / real(config%time%n_steps,dp)
         max_steps = config%time%n_steps + 2
      else
         capacity = 64
         dt = config%time%dt_init
         max_steps = 1000000
      end if
      allocate(history(nnode,capacity), time_history(capacity), &
         linear_history(capacity), penalty_history(capacity), stat=stat)
      if (stat /= 0) then
         call set_status(status, mao_allocation_error, &
            'unable to allocate pricing history')
         return
      end if

      current = payoff
      remaining = config%opt%ttm
      count = 1
      history(:,1) = current
      time_history(1) = remaining
      linear_history(1) = 0
      penalty_history(1) = 0
      smooth_index = 0
      if (config%opt%exercise_type == exercise_american) then
         large = 1.0_dp / config%fd%tol
      else
         large = 0.0_dp
      end if

      if (present(solver_tolerance)) then
         linear_tolerance = solver_tolerance
      else
         linear_tolerance = max(1.0e-12_dp,min(1.0e-9_dp,0.01_dp*config%fd%tol))
      end if
      if (present(solver_max_iterations)) then
         max_linear_iterations = solver_max_iterations
      else
         max_linear_iterations = max(500,min(20000,5*nnode))
      end if

      do while (remaining > 10.0_dp*epsilon(remaining)*config%opt%ttm)
         if (count >= max_steps) then
            call set_status(status, mao_step_failure, &
               'timestep loop exceeded its safety limit')
            return
         end if
         current_dt = min(dt,remaining)
         if (current_dt <= 100.0_dp*tiny(1.0_dp)) then
            call set_status(status, mao_step_failure, &
               'adaptive timestep became too small')
            return
         end if

         previous = current
         if (smooth_index < config%fd%max_smooth) then
            theta_step = 1.0_dp
            smooth_index = smooth_index + 1
         else
            theta_step = config%fd%theta
         end if

         call csr_matvec(operator,previous,operator_value)
         rhs = previous + (1.0_dp-theta_step)*current_dt*operator_value
         candidate = previous
         total_linear_iterations = 0

         do penalty_iter = 1, config%fd%max_iter
            active_old = payoff > candidate
            penalty = 0.0_dp
            where (active_old) penalty = large

            call csr_shifted_identity(operator,-theta_step*current_dt, &
               penalty,system,local_status)
            if (.not. local_status%ok()) then
               status = local_status
               return
            end if
            rhs_penalty = rhs + penalty*payoff
            call bicgstab_solve(system,rhs_penalty,candidate,linear_tolerance, &
               max_linear_iterations,linear_iterations,residual_norm,local_status)
            total_linear_iterations = total_linear_iterations + linear_iterations
            if (.not. local_status%ok()) then
               call set_status(status, mao_solver_failure, &
                  trim(local_status%message)//'; residual='//real_to_string(residual_norm))
               return
            end if

            active_new = payoff > candidate
            relative_change = maxval(abs(candidate-current) / &
               max(1.0_dp,abs(candidate)))
            if (config%opt%exercise_type /= exercise_american .or. &
                all(active_new .eqv. active_old) .or. &
                relative_change < config%fd%tol .or. &
                penalty_iter >= config%fd%max_iter) exit
            current = candidate
         end do

         current = candidate
         remaining = max(0.0_dp,remaining-current_dt)

         if (config%time%ts_type /= timestep_constant .and. remaining > 0.0_dp) then
            max_relative_change = maxval(abs(current-previous) / &
               max(config%time%scale_d,abs(current),abs(previous)))
            if (max_relative_change <= 100.0_dp*tiny(1.0_dp)) then
               next_dt = remaining
            else
               next_dt = current_dt * config%time%d_norm / max_relative_change
            end if
            dt = min(next_dt,remaining)
         end if

         count = count + 1
         if (count > capacity) then
            call grow_history(history,time_history,linear_history, &
               penalty_history,capacity,local_status)
            if (.not. local_status%ok()) then
               status = local_status
               return
            end if
         end if
         history(:,count) = current
         time_history(count) = remaining
         linear_history(count) = total_linear_iterations
         penalty_history(count) = penalty_iter
      end do

      allocate(result%value(nnode,count),result%time(count), &
         result%linear_iterations(count),result%penalty_iterations(count),stat=stat)
      if (stat /= 0) then
         call set_status(status, mao_allocation_error, &
            'unable to allocate final pricing result')
         return
      end if
      result%value = history(:,1:count)
      result%time = time_history(1:count)
      result%linear_iterations = linear_history(1:count)
      result%penalty_iterations = penalty_history(1:count)
   end subroutine price_multi_asset

   subroutine grow_history(history,times,linear_iterations, &
      penalty_iterations,capacity,status)
      real(dp), allocatable, intent(inout) :: history(:,:), times(:)
      integer, allocatable, intent(inout) :: linear_iterations(:)
      integer, allocatable, intent(inout) :: penalty_iterations(:)
      integer, intent(inout) :: capacity
      type(status_type), intent(out) :: status
      real(dp), allocatable :: new_history(:,:), new_times(:)
      integer, allocatable :: new_linear(:), new_penalty(:)
      integer :: new_capacity, stat

      call clear_status(status)
      new_capacity = 2*capacity
      allocate(new_history(size(history,1),new_capacity), &
         new_times(new_capacity),new_linear(new_capacity), &
         new_penalty(new_capacity),stat=stat)
      if (stat /= 0) then
         call set_status(status, mao_allocation_error, &
            'unable to grow pricing history')
         return
      end if
      new_history(:,1:capacity) = history
      new_times(1:capacity) = times
      new_linear(1:capacity) = linear_iterations
      new_penalty(1:capacity) = penalty_iterations
      call move_alloc(new_history,history)
      call move_alloc(new_times,times)
      call move_alloc(new_linear,linear_iterations)
      call move_alloc(new_penalty,penalty_iterations)
      capacity = new_capacity
   end subroutine grow_history

   function real_to_string(x) result(text)
      real(dp), intent(in) :: x
      character(len=:), allocatable :: text
      character(len=48) :: buffer

      write(buffer,'(es16.8)') x
      text = trim(adjustl(buffer))
   end function real_to_string

end module mao_pricing
