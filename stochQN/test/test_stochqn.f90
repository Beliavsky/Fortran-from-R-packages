program test_stochqn
   use stochqn_kinds, only : dp
   use stochqn_core
   use stochqn_guided
   use stochqn_logistic
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none

   type :: quad_data_t
      real(dp), allocatable :: a(:,:), b(:)
   end type quad_data_t

   integer :: failures

   failures = 0
   call test_low_level_olbfgs(failures)
   call test_guided_optimizers(failures)
   call test_adaqn_fisher_memory(failures)
   call test_adaqn_validation_rollback(failures)
   call test_logistic_derivatives(failures)
   call test_logistic_model(failures)

   if (failures /= 0) then
      print '(a,i0)', 'FAILED TESTS: ', failures
      error stop 1
   end if
   print '(a)', 'All stochQN tests passed.'

contains

   subroutine assert_true(condition, name, failures)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: name
      integer, intent(inout) :: failures
      if (.not. condition) then
         print '(a)', 'FAIL: '//trim(name)
         failures = failures + 1
      end if
   end subroutine assert_true


   subroutine test_low_level_olbfgs(failures)
      integer, intent(inout) :: failures
      type(olbfgs_t) :: optimizer
      type(stochqn_request_t) :: request
      real(dp) :: x(2), gradient(2)
      integer :: stat

      x = [2.0_dp, -1.0_dp]
      gradient = 0.0_dp
      call optimizer%initialize(2, mem_size=4, hessian_initialization=1.0_dp, &
           min_curvature=0.0_dp, stat=stat)
      call assert_true(stat == 0, 'oLBFGS initialization', failures)
      call optimizer%advance(0.1_dp, x, gradient, request)
      call assert_true(request%task == task_calc_grad, 'oLBFGS first gradient request', failures)

      gradient = request%x
      call optimizer%advance(0.1_dp, x, gradient, request)
      call assert_true(request%task == task_calc_grad_same_batch, &
                       'oLBFGS same-batch request', failures)
      call assert_true(request%status == status_updated, 'oLBFGS update status', failures)

      gradient = request%x
      call optimizer%advance(0.1_dp, x, gradient, request)
      call assert_true(request%task == task_calc_grad, 'oLBFGS loop request', failures)
      call assert_true(optimizer%correction_pairs() == 1, 'oLBFGS correction pair', failures)

      x = [2.0_dp, -1.0_dp]
      gradient = 0.0_dp
      call optimizer%initialize(2, mem_size=4, hessian_initialization=1.0_dp, &
           min_curvature=0.0_dp, stat=stat)
      call optimizer%advance(0.1_dp, x, gradient, request)
      do while (optimizer%get_iteration() < 5 .or. request%task /= task_calc_grad)
         gradient = request%x
         call optimizer%advance(0.1_dp, x, gradient, request)
      end do
      call assert_true(maxval(abs(x - [1.18098_dp, -0.59049_dp])) < 1.0e-14_dp, &
                       'oLBFGS original C reference sequence', failures)
      call assert_true(optimizer%correction_pairs() == 4, &
                       'oLBFGS original C reference memory', failures)

      call optimizer%initialize(2, mem_size=0, stat=stat)
      call assert_true(stat /= 0, 'reject zero memory size', failures)
   end subroutine test_low_level_olbfgs

   subroutine test_guided_optimizers(failures)
      integer, intent(inout) :: failures
      type(stochqn_run_result_t) :: result
      type(quad_data_t) :: data
      real(dp) :: x(3), solution(3), initial_norm, final_norm

      allocate(data%a(3, 3), data%b(3))
      data%a = 0.0_dp
      data%a(1, 1) = 1.0_dp
      data%a(2, 2) = 2.0_dp
      data%a(3, 3) = 4.0_dp
      data%b = [1.0_dp, -2.0_dp, 0.5_dp]
      solution = [1.0_dp, -1.0_dp, 0.125_dp]

      x = [4.0_dp, -3.0_dp, 2.0_dp]
      call optimize_olbfgs(x, 300, 0.2_dp, quadratic_gradient, result, data, &
                           mem_size=6, hessian_initialization=0.0_dp)
      call assert_true(result%valid, 'guided oLBFGS valid', failures)
      call assert_true(maxval(abs(x - solution)) < 2.0e-8_dp, 'guided oLBFGS solution', failures)

      x = [4.0_dp, -3.0_dp, 2.0_dp]
      call optimize_sqn(x, 400, 0.12_dp, quadratic_gradient, result, quadratic_hess_vec, data, &
                        mem_size=6, bfgs_update_frequency=5)
      call assert_true(result%valid, 'guided SQN Hessian-vector valid', failures)
      call assert_true(maxval(abs(x - solution)) < 1.0e-5_dp, &
                       'guided SQN Hessian-vector solution', failures)
      call assert_true(result%hessian_vector_evaluations > 0, 'SQN Hessian-vector used', failures)

      x = [4.0_dp, -3.0_dp, 2.0_dp]
      call optimize_sqn(x, 400, 0.12_dp, quadratic_gradient, result, user_data=data, &
                        mem_size=6, bfgs_update_frequency=5, use_gradient_difference=.true.)
      call assert_true(result%valid, 'guided SQN gradient-difference valid', failures)
      call assert_true(maxval(abs(x - solution)) < 1.0e-5_dp, &
                       'guided SQN gradient-difference solution', failures)

      x = [4.0_dp, -3.0_dp, 2.0_dp]
      initial_norm = norm2(matmul(data%a, x) - data%b)
      call optimize_adaqn(x, 500, 0.05_dp, quadratic_gradient, result, user_data=data, &
                          mem_size=6, fisher_size=20, bfgs_update_frequency=5, &
                          maximum_increase=0.0_dp, rmsprop_weight=0.9_dp, use_gradient_difference=.true.)
      final_norm = norm2(matmul(data%a, x) - data%b)
      call assert_true(result%valid, 'guided adaQN Fisher valid', failures)
      call assert_true(final_norm < 0.01_dp * initial_norm, 'guided adaQN gradient-difference progress', failures)
   end subroutine test_guided_optimizers


   subroutine test_adaqn_fisher_memory(failures)
      integer, intent(inout) :: failures
      type(adaqn_t) :: optimizer
      type(stochqn_request_t) :: request
      real(dp) :: x(1), gradient(1), f
      integer :: stat, k

      x = 2.0_dp
      gradient = 0.0_dp
      f = 0.0_dp
      call optimizer%initialize(1, mem_size=4, fisher_size=4, bfgs_update_frequency=2, &
           maximum_increase=0.0_dp, min_curvature=0.0_dp, &
           scaling_regularization=1.0e-4_dp, rmsprop_weight=0.9_dp, &
           use_gradient_difference=.false., stat=stat)
      call assert_true(stat == 0, 'adaQN Fisher initialization', failures)
      call optimizer%advance(0.05_dp, x, f, gradient, request)
      do k = 1, 6
         call assert_true(request%task == task_calc_grad, 'adaQN Fisher gradient request', failures)
         gradient = request%x
         call optimizer%advance(0.05_dp, x, f, gradient, request)
      end do
      call assert_true(optimizer%correction_pairs() > 0, 'adaQN Fisher correction pair', failures)
      call assert_true(all(ieee_is_finite(x)), 'adaQN Fisher finite iterate', failures)
   end subroutine test_adaqn_fisher_memory


   subroutine test_adaqn_validation_rollback(failures)
      integer, intent(inout) :: failures
      type(adaqn_t) :: optimizer
      type(stochqn_request_t) :: request
      real(dp) :: x(1), gradient(1), f
      integer :: stat, k

      x = 2.0_dp
      gradient = 0.0_dp
      f = 0.0_dp
      call optimizer%initialize(1, mem_size=4, fisher_size=4, bfgs_update_frequency=2, &
           maximum_increase=1.01_dp, min_curvature=0.0_dp, &
           scaling_regularization=1.0e-4_dp, rmsprop_weight=0.9_dp, &
           use_gradient_difference=.false., stat=stat)
      call assert_true(stat == 0, 'adaQN rollback initialization', failures)
      call optimizer%advance(0.05_dp, x, f, gradient, request)

      do k = 1, 2
         gradient = request%x
         call optimizer%advance(0.05_dp, x, f, gradient, request)
      end do
      call assert_true(request%task == task_calc_fun_value, &
                       'adaQN first validation request', failures)
      f = 0.5_dp * request%x(1) * request%x(1)
      call optimizer%advance(0.05_dp, x, f, gradient, request)

      do k = 1, 2
         gradient = request%x
         call optimizer%advance(0.05_dp, x, f, gradient, request)
      end do
      call assert_true(request%task == task_calc_fun_value, &
                       'adaQN later validation request', failures)
      f = 1.0e6_dp
      call optimizer%advance(0.05_dp, x, f, gradient, request)
      call assert_true(request%info == info_function_increased, &
                       'adaQN validation rollback info', failures)
      call assert_true(request%task == task_calc_grad, &
                       'adaQN rollback resumes gradients', failures)
   end subroutine test_adaqn_validation_rollback

   subroutine test_logistic_derivatives(failures)
      integer, intent(inout) :: failures
      real(dp) :: x(6, 3), y(6), beta(3), gradient(3), finite_difference(3)
      real(dp) :: vector(3), hessian_vector(3), fd_hessian_vector(3)
      real(dp) :: beta_plus(3), beta_minus(3), grad_plus(3), grad_minus(3)
      real(dp) :: h, fplus, fminus
      integer :: j

      x = reshape([ &
         1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
        -2.0_dp,-1.0_dp, 0.0_dp, 0.5_dp, 1.0_dp, 2.0_dp, &
         0.5_dp,-1.0_dp, 1.0_dp,-0.5_dp, 2.0_dp, 1.5_dp], shape(x))
      y = [0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
      beta = [0.1_dp, 0.4_dp, -0.2_dp]
      vector = [0.3_dp, -0.7_dp, 0.2_dp]
      h = 1.0e-6_dp

      call logistic_gradient(beta, x, y, gradient, lambda=1.0e-3_dp)
      do j = 1, 3
         beta_plus = beta
         beta_minus = beta
         beta_plus(j) = beta_plus(j) + h
         beta_minus(j) = beta_minus(j) - h
         fplus = logistic_loss(beta_plus, x, y, lambda=1.0e-3_dp)
         fminus = logistic_loss(beta_minus, x, y, lambda=1.0e-3_dp)
         finite_difference(j) = (fplus - fminus) / (2.0_dp * h)
      end do
      call assert_true(maxval(abs(gradient - finite_difference)) < 2.0e-9_dp, &
                       'logistic gradient finite difference', failures)

      call logistic_hessian_vector(beta, vector, x, y, hessian_vector, lambda=1.0e-3_dp)
      beta_plus = beta + h * vector
      beta_minus = beta - h * vector
      call logistic_gradient(beta_plus, x, y, grad_plus, lambda=1.0e-3_dp)
      call logistic_gradient(beta_minus, x, y, grad_minus, lambda=1.0e-3_dp)
      fd_hessian_vector = (grad_plus - grad_minus) / (2.0_dp * h)
      call assert_true(maxval(abs(hessian_vector - fd_hessian_vector)) < 2.0e-9_dp, &
                       'logistic Hessian-vector finite difference', failures)
   end subroutine test_logistic_derivatives

   subroutine test_logistic_model(failures)
      integer, intent(inout) :: failures
      type(stochastic_logistic_t) :: model
      real(dp) :: x(12, 2), y(12)
      real(dp), allocatable :: probability(:), coefficients(:)
      integer, allocatable :: predicted(:)
      integer :: i, stat, correct

      do i = 1, 12
         x(i, 1) = -2.0_dp + 4.0_dp * real(i - 1, dp) / 11.0_dp
         x(i, 2) = sin(real(i, dp))
         y(i) = merge(1.0_dp, 0.0_dp, 1.5_dp * x(i, 1) - 0.3_dp * x(i, 2) > 0.0_dp)
      end do
      call model%initialize(2, optimizer_kind=logistic_olbfgs, fit_intercept=.true., &
                            lambda=1.0e-4_dp, initial_step=0.2_dp, stat=stat)
      call assert_true(stat == 0, 'logistic model initialize', failures)
      do i = 1, 150
         call model%partial_fit(x, y, stat=stat)
         if (stat /= 0) exit
      end do
      call assert_true(stat == 0, 'logistic model partial fit', failures)
      probability = model%predict_probability(x)
      predicted = model%predict_class(x)
      correct = count(predicted == int(y))
      call assert_true(correct >= 11, 'logistic model classification accuracy', failures)
      coefficients = model%get_coefficients()
      call assert_true(size(coefficients) == 3, 'logistic coefficient size', failures)
      call assert_true(all(probability >= 0.0_dp .and. probability <= 1.0_dp), &
                       'logistic probabilities bounded', failures)
   end subroutine test_logistic_model

   subroutine quadratic_gradient(x, task, iteration, gradient, user_data)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: task, iteration
      real(dp), intent(out) :: gradient(:)
      class(*), intent(inout), optional :: user_data
      select type (data => user_data)
      type is (quad_data_t)
         gradient = matmul(data%a, x) - data%b
      class default
         gradient = 0.0_dp
      end select
      if (task + iteration < -huge(1)) gradient = 0.0_dp
   end subroutine quadratic_gradient

   subroutine quadratic_hess_vec(x, vector, iteration, product, user_data)
      real(dp), intent(in) :: x(:), vector(:)
      integer, intent(in) :: iteration
      real(dp), intent(out) :: product(:)
      class(*), intent(inout), optional :: user_data
      select type (data => user_data)
      type is (quad_data_t)
         product = matmul(data%a, vector)
      class default
         product = 0.0_dp
      end select
      if (sum(x) + real(iteration, dp) < -huge(1.0_dp)) product = 0.0_dp
   end subroutine quadratic_hess_vec

end program test_stochqn
