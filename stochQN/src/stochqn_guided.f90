module stochqn_guided
   use stochqn_kinds, only : dp
   use stochqn_core, only : olbfgs_t, sqn_t, adaqn_t, stochqn_request_t, &
      task_calc_grad, task_calc_grad_same_batch, task_calc_grad_big_batch, &
      task_calc_hess_vec, task_calc_fun_value, task_invalid, info_ok
   implicit none
   private

   type, public :: stochqn_run_result_t
      integer :: iterations = 0
      integer :: gradient_evaluations = 0
      integer :: hessian_vector_evaluations = 0
      integer :: function_evaluations = 0
      integer :: last_info = info_ok
      logical :: valid = .true.
   end type stochqn_run_result_t

   abstract interface
      subroutine stochastic_gradient_interface(x, task, iteration, gradient, user_data)
         import dp
         real(dp), intent(in) :: x(:)
         integer, intent(in) :: task, iteration
         real(dp), intent(out) :: gradient(:)
         class(*), intent(inout), optional :: user_data
      end subroutine stochastic_gradient_interface

      subroutine stochastic_hessian_vector_interface(x, vector, iteration, product, user_data)
         import dp
         real(dp), intent(in) :: x(:), vector(:)
         integer, intent(in) :: iteration
         real(dp), intent(out) :: product(:)
         class(*), intent(inout), optional :: user_data
      end subroutine stochastic_hessian_vector_interface

      function stochastic_objective_interface(x, iteration, user_data) result(value)
         import dp
         real(dp), intent(in) :: x(:)
         integer, intent(in) :: iteration
         class(*), intent(inout), optional :: user_data
         real(dp) :: value
      end function stochastic_objective_interface

      function step_schedule_interface(iteration) result(multiplier)
         import dp
         integer, intent(in) :: iteration
         real(dp) :: multiplier
      end function step_schedule_interface
   end interface

   public :: stochastic_gradient_interface, stochastic_hessian_vector_interface
   public :: stochastic_objective_interface, step_schedule_interface
   public :: optimize_olbfgs, optimize_sqn, optimize_adaqn

contains

   subroutine evaluate_gradient(callback, x, task, iteration, gradient, user_data)
      procedure(stochastic_gradient_interface) :: callback
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: task, iteration
      real(dp), intent(out) :: gradient(:)
      class(*), intent(inout), optional :: user_data
      if (present(user_data)) then
         call callback(x, task, iteration, gradient, user_data)
      else
         call callback(x, task, iteration, gradient)
      end if
   end subroutine evaluate_gradient

   subroutine evaluate_hessian_vector(callback, x, vector, iteration, product, user_data)
      procedure(stochastic_hessian_vector_interface) :: callback
      real(dp), intent(in) :: x(:), vector(:)
      integer, intent(in) :: iteration
      real(dp), intent(out) :: product(:)
      class(*), intent(inout), optional :: user_data
      if (present(user_data)) then
         call callback(x, vector, iteration, product, user_data)
      else
         call callback(x, vector, iteration, product)
      end if
   end subroutine evaluate_hessian_vector

   function evaluate_objective(callback, x, iteration, user_data) result(value)
      procedure(stochastic_objective_interface) :: callback
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: iteration
      class(*), intent(inout), optional :: user_data
      real(dp) :: value
      if (present(user_data)) then
         value = callback(x, iteration, user_data)
      else
         value = callback(x, iteration)
      end if
   end function evaluate_objective

   pure real(dp) function default_short_schedule(iteration)
      integer, intent(in) :: iteration
      default_short_schedule = 1.0_dp / sqrt(real(iteration, dp) / 10.0_dp + 1.0_dp)
   end function default_short_schedule

   pure real(dp) function default_long_schedule(iteration)
      integer, intent(in) :: iteration
      default_long_schedule = 1.0_dp / sqrt(real(iteration, dp) / 100.0_dp + 1.0_dp)
   end function default_long_schedule

   real(dp) function choose_step(initial_step, iteration, schedule, long_schedule) result(step)
      real(dp), intent(in) :: initial_step
      integer, intent(in) :: iteration
      procedure(step_schedule_interface), optional :: schedule
      logical, intent(in) :: long_schedule
      if (present(schedule)) then
         step = initial_step * schedule(iteration)
      else if (long_schedule) then
         step = initial_step * default_long_schedule(iteration)
      else
         step = initial_step * default_short_schedule(iteration)
      end if
   end function choose_step

   subroutine optimize_olbfgs(x, n_iterations, initial_step, gradient_callback, result, &
                              user_data, step_schedule, mem_size, hessian_initialization, &
                              y_regularization, min_curvature, check_finite)
      real(dp), intent(inout) :: x(:)
      integer, intent(in) :: n_iterations
      real(dp), intent(in) :: initial_step
      procedure(stochastic_gradient_interface) :: gradient_callback
      type(stochqn_run_result_t), intent(out) :: result
      class(*), intent(inout), optional :: user_data
      procedure(step_schedule_interface), optional :: step_schedule
      integer, intent(in), optional :: mem_size
      real(dp), intent(in), optional :: hessian_initialization, y_regularization, min_curvature
      logical, intent(in), optional :: check_finite
      type(olbfgs_t) :: optimizer
      type(stochqn_request_t) :: request
      real(dp), allocatable :: gradient(:)
      real(dp) :: step, h0, yr, curv
      integer :: m, stat
      logical :: check

      result = stochqn_run_result_t()
      m = 10
      if (present(mem_size)) m = mem_size
      h0 = 0.0_dp
      if (present(hessian_initialization)) h0 = hessian_initialization
      yr = 0.0_dp
      if (present(y_regularization)) yr = y_regularization
      curv = 1.0e-4_dp
      if (present(min_curvature)) curv = min_curvature
      check = .true.
      if (present(check_finite)) check = check_finite
      call optimizer%initialize(size(x), m, h0, yr, curv, check, stat)
      if (stat /= 0 .or. n_iterations < 0 .or. initial_step <= 0.0_dp) then
         result%valid = .false.
         return
      end if

      allocate(gradient(size(x)), source=0.0_dp)
      call optimizer%advance(initial_step, x, gradient, request)
      do
         if (request%task == task_calc_grad .and. request%iteration >= n_iterations) exit
         select case (request%task)
         case (task_calc_grad, task_calc_grad_same_batch)
            call evaluate_gradient(gradient_callback, request%x, request%task, request%iteration, &
                                   gradient, user_data)
            result%gradient_evaluations = result%gradient_evaluations + 1
         case default
            result%valid = .false.
            exit
         end select
         step = choose_step(initial_step, request%iteration, step_schedule, .false.)
         call optimizer%advance(step, x, gradient, request)
         result%last_info = request%info
         if (request%task == task_invalid) then
            result%valid = .false.
            exit
         end if
      end do
      result%iterations = optimizer%get_iteration()
   end subroutine optimize_olbfgs

   subroutine optimize_sqn(x, n_iterations, initial_step, gradient_callback, result, &
                           hessian_vector_callback, user_data, step_schedule, mem_size, &
                           bfgs_update_frequency, min_curvature, use_gradient_difference, &
                           y_regularization, check_finite)
      real(dp), intent(inout) :: x(:)
      integer, intent(in) :: n_iterations
      real(dp), intent(in) :: initial_step
      procedure(stochastic_gradient_interface) :: gradient_callback
      type(stochqn_run_result_t), intent(out) :: result
      procedure(stochastic_hessian_vector_interface), optional :: hessian_vector_callback
      class(*), intent(inout), optional :: user_data
      procedure(step_schedule_interface), optional :: step_schedule
      integer, intent(in), optional :: mem_size, bfgs_update_frequency
      real(dp), intent(in), optional :: min_curvature, y_regularization
      logical, intent(in), optional :: use_gradient_difference, check_finite
      type(sqn_t) :: optimizer
      type(stochqn_request_t) :: request
      real(dp), allocatable :: gradient(:), hessian_vector(:)
      real(dp) :: step, curv, yr
      integer :: m, freq, stat
      logical :: grad_diff, check

      result = stochqn_run_result_t()
      m = 10
      if (present(mem_size)) m = mem_size
      freq = 20
      if (present(bfgs_update_frequency)) freq = bfgs_update_frequency
      curv = 1.0e-4_dp
      if (present(min_curvature)) curv = min_curvature
      yr = 0.0_dp
      if (present(y_regularization)) yr = y_regularization
      grad_diff = .false.
      if (present(use_gradient_difference)) grad_diff = use_gradient_difference
      check = .true.
      if (present(check_finite)) check = check_finite
      if (.not. grad_diff .and. .not. present(hessian_vector_callback)) then
         result%valid = .false.
         return
      end if
      call optimizer%initialize(size(x), m, freq, curv, grad_diff, yr, check, stat)
      if (stat /= 0 .or. n_iterations < 0 .or. initial_step <= 0.0_dp) then
         result%valid = .false.
         return
      end if

      allocate(gradient(size(x)), source=0.0_dp)
      allocate(hessian_vector(size(x)), source=0.0_dp)
      call optimizer%advance(initial_step, x, gradient, hessian_vector, request)
      do
         if (request%task == task_calc_grad .and. request%iteration >= n_iterations) exit
         select case (request%task)
         case (task_calc_grad, task_calc_grad_big_batch)
            call evaluate_gradient(gradient_callback, request%x, request%task, request%iteration, &
                                   gradient, user_data)
            result%gradient_evaluations = result%gradient_evaluations + 1
         case (task_calc_hess_vec)
            call evaluate_hessian_vector(hessian_vector_callback, request%x, request%vec, &
                                         request%iteration, hessian_vector, user_data)
            result%hessian_vector_evaluations = result%hessian_vector_evaluations + 1
         case default
            result%valid = .false.
            exit
         end select
         step = choose_step(initial_step, request%iteration, step_schedule, .false.)
         call optimizer%advance(step, x, gradient, hessian_vector, request)
         result%last_info = request%info
         if (request%task == task_invalid) then
            result%valid = .false.
            exit
         end if
      end do
      result%iterations = optimizer%get_iteration()
   end subroutine optimize_sqn

   subroutine optimize_adaqn(x, n_iterations, initial_step, gradient_callback, result, &
                             objective_callback, user_data, step_schedule, mem_size, fisher_size, &
                             bfgs_update_frequency, maximum_increase, min_curvature, &
                             scaling_regularization, rmsprop_weight, use_gradient_difference, &
                             y_regularization, check_finite)
      real(dp), intent(inout) :: x(:)
      integer, intent(in) :: n_iterations
      real(dp), intent(in) :: initial_step
      procedure(stochastic_gradient_interface) :: gradient_callback
      type(stochqn_run_result_t), intent(out) :: result
      procedure(stochastic_objective_interface), optional :: objective_callback
      class(*), intent(inout), optional :: user_data
      procedure(step_schedule_interface), optional :: step_schedule
      integer, intent(in), optional :: mem_size, fisher_size, bfgs_update_frequency
      real(dp), intent(in), optional :: maximum_increase, min_curvature
      real(dp), intent(in), optional :: scaling_regularization, rmsprop_weight, y_regularization
      logical, intent(in), optional :: use_gradient_difference, check_finite
      type(adaqn_t) :: optimizer
      type(stochqn_request_t) :: request
      real(dp), allocatable :: gradient(:)
      real(dp) :: step, f, maxinc, curv, scal, rms, yr
      integer :: m, fsize, freq, stat
      logical :: grad_diff, check

      result = stochqn_run_result_t()
      m = 10
      if (present(mem_size)) m = mem_size
      fsize = 100
      if (present(fisher_size)) fsize = fisher_size
      freq = 20
      if (present(bfgs_update_frequency)) freq = bfgs_update_frequency
      maxinc = 1.01_dp
      if (present(maximum_increase)) maxinc = maximum_increase
      curv = 1.0e-4_dp
      if (present(min_curvature)) curv = min_curvature
      scal = 1.0e-4_dp
      if (present(scaling_regularization)) scal = scaling_regularization
      rms = 0.9_dp
      if (present(rmsprop_weight)) rms = rmsprop_weight
      yr = 0.0_dp
      if (present(y_regularization)) yr = y_regularization
      grad_diff = .false.
      if (present(use_gradient_difference)) grad_diff = use_gradient_difference
      check = .true.
      if (present(check_finite)) check = check_finite
      if (maxinc > 0.0_dp .and. .not. present(objective_callback)) then
         result%valid = .false.
         return
      end if
      call optimizer%initialize(size(x), m, fsize, freq, maxinc, curv, scal, rms, &
                                grad_diff, yr, check, stat)
      if (stat /= 0 .or. n_iterations < 0 .or. initial_step <= 0.0_dp) then
         result%valid = .false.
         return
      end if

      allocate(gradient(size(x)), source=0.0_dp)
      f = 0.0_dp
      call optimizer%advance(initial_step, x, f, gradient, request)
      do
         if (request%task == task_calc_grad .and. request%iteration >= n_iterations) exit
         select case (request%task)
         case (task_calc_grad, task_calc_grad_big_batch)
            call evaluate_gradient(gradient_callback, request%x, request%task, request%iteration, &
                                   gradient, user_data)
            result%gradient_evaluations = result%gradient_evaluations + 1
         case (task_calc_fun_value)
            f = evaluate_objective(objective_callback, request%x, request%iteration, user_data)
            result%function_evaluations = result%function_evaluations + 1
         case default
            result%valid = .false.
            exit
         end select
         step = choose_step(initial_step, request%iteration, step_schedule, .true.)
         call optimizer%advance(step, x, f, gradient, request)
         result%last_info = request%info
         if (request%task == task_invalid) then
            result%valid = .false.
            exit
         end if
      end do
      result%iterations = optimizer%get_iteration()
   end subroutine optimize_adaqn

end module stochqn_guided
