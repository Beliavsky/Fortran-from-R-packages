module stochqn_logistic
   use stochqn_kinds, only : dp
   use stochqn_core, only : olbfgs_t, sqn_t, adaqn_t, stochqn_request_t, &
      task_calc_grad, task_calc_grad_same_batch, task_calc_grad_big_batch, &
      task_calc_hess_vec, task_calc_fun_value, task_invalid, info_ok
   implicit none
   private

   integer, parameter, public :: logistic_olbfgs = 1
   integer, parameter, public :: logistic_sqn = 2
   integer, parameter, public :: logistic_adaqn = 3

   type, public :: stochastic_logistic_t
      private
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: history_x(:,:)
      real(dp), allocatable :: history_y(:)
      real(dp), allocatable :: history_w(:)
      type(olbfgs_t) :: olbfgs
      type(sqn_t) :: sqn
      type(adaqn_t) :: adaqn
      type(stochqn_request_t) :: request
      integer :: optimizer_kind = logistic_adaqn
      integer :: n_features = 0
      logical :: fit_intercept = .true.
      logical :: initialized = .false.
      real(dp) :: lambda = 1.0e-5_dp
      real(dp) :: initial_step = 1.0e-2_dp
   contains
      procedure, public :: initialize => logistic_model_initialize
      procedure, public :: partial_fit => logistic_model_partial_fit
      procedure, public :: predict_probability => logistic_model_predict_probability
      procedure, public :: predict_class => logistic_model_predict_class
      procedure, public :: get_coefficients => logistic_model_get_coefficients
      procedure, public :: get_iteration => logistic_model_get_iteration
      procedure, public :: clear_history => logistic_model_clear_history
   end type stochastic_logistic_t

   public :: logistic_loss, logistic_gradient, logistic_hessian_vector
   public :: logistic_predict_probability

contains

   elemental real(dp) function stable_sigmoid(z)
      real(dp), intent(in) :: z
      if (z >= 0.0_dp) then
         stable_sigmoid = 1.0_dp / (1.0_dp + exp(-z))
      else
         stable_sigmoid = exp(z) / (1.0_dp + exp(z))
      end if
   end function stable_sigmoid

   function logistic_loss(coefficients, x, y, lambda, weights) result(loss)
      real(dp), intent(in) :: coefficients(:)
      real(dp), intent(in) :: x(:,:), y(:)
      real(dp), intent(in), optional :: lambda, weights(:)
      real(dp) :: loss
      real(dp), allocatable :: eta(:), terms(:), w(:)
      real(dp) :: reg, weight_sum

      if (size(x, 1) /= size(y) .or. size(x, 2) /= size(coefficients) .or. size(y) == 0) then
         loss = huge(1.0_dp)
         return
      end if
      allocate(eta(size(y)), terms(size(y)), w(size(y)))
      eta = matmul(x, coefficients)
      terms = max(eta, 0.0_dp) - y * eta + log(1.0_dp + exp(-abs(eta)))
      if (present(weights)) then
         if (size(weights) /= size(y) .or. any(weights < 0.0_dp)) then
            loss = huge(1.0_dp)
            return
         end if
         w = weights
      else
         w = 1.0_dp
      end if
      weight_sum = sum(w)
      if (weight_sum <= 0.0_dp) then
         loss = huge(1.0_dp)
         return
      end if
      reg = 1.0e-5_dp
      if (present(lambda)) reg = lambda
      loss = sum(w * terms) / weight_sum + reg * dot_product(coefficients, coefficients)
   end function logistic_loss

   subroutine logistic_gradient(coefficients, x, y, gradient, lambda, weights)
      real(dp), intent(in) :: coefficients(:)
      real(dp), intent(in) :: x(:,:), y(:)
      real(dp), intent(out) :: gradient(:)
      real(dp), intent(in), optional :: lambda, weights(:)
      real(dp), allocatable :: probabilities(:), residual(:), w(:)
      real(dp) :: reg, weight_sum

      if (size(x, 1) /= size(y) .or. size(x, 2) /= size(coefficients) .or. &
          size(gradient) /= size(coefficients) .or. size(y) == 0) then
         gradient = 0.0_dp
         return
      end if
      allocate(probabilities(size(y)), residual(size(y)), w(size(y)))
      probabilities = stable_sigmoid(matmul(x, coefficients))
      if (present(weights)) then
         if (size(weights) /= size(y) .or. any(weights < 0.0_dp)) then
            gradient = 0.0_dp
            return
         end if
         w = weights
      else
         w = 1.0_dp
      end if
      weight_sum = sum(w)
      if (weight_sum <= 0.0_dp) then
         gradient = 0.0_dp
         return
      end if
      residual = w * (probabilities - y)
      reg = 1.0e-5_dp
      if (present(lambda)) reg = lambda
      gradient = matmul(transpose(x), residual) / weight_sum + 2.0_dp * reg * coefficients
   end subroutine logistic_gradient

   subroutine logistic_hessian_vector(coefficients, vector, x, y, product, lambda, weights)
      real(dp), intent(in) :: coefficients(:), vector(:)
      real(dp), intent(in) :: x(:,:), y(:)
      real(dp), intent(out) :: product(:)
      real(dp), intent(in), optional :: lambda, weights(:)
      real(dp), allocatable :: probabilities(:), diagonal(:), w(:), xv(:)
      real(dp) :: reg, weight_sum

      if (size(x, 1) /= size(y) .or. size(x, 2) /= size(coefficients) .or. &
          size(vector) /= size(coefficients) .or. size(product) /= size(coefficients)) then
         product = 0.0_dp
         return
      end if
      allocate(probabilities(size(y)), diagonal(size(y)), w(size(y)), xv(size(y)))
      probabilities = stable_sigmoid(matmul(x, coefficients))
      if (present(weights)) then
         if (size(weights) /= size(y) .or. any(weights < 0.0_dp)) then
            product = 0.0_dp
            return
         end if
         w = weights
      else
         w = 1.0_dp
      end if
      weight_sum = sum(w)
      if (weight_sum <= 0.0_dp) then
         product = 0.0_dp
         return
      end if
      diagonal = w * probabilities * (1.0_dp - probabilities)
      xv = matmul(x, vector)
      reg = 1.0e-5_dp
      if (present(lambda)) reg = lambda
      product = matmul(transpose(x), diagonal * xv) / weight_sum + 2.0_dp * reg * vector
   end subroutine logistic_hessian_vector

   function logistic_predict_probability(coefficients, x) result(probability)
      real(dp), intent(in) :: coefficients(:), x(:,:)
      real(dp), allocatable :: probability(:)
      allocate(probability(size(x, 1)))
      if (size(x, 2) /= size(coefficients)) then
         probability = 0.0_dp
      else
         probability = stable_sigmoid(matmul(x, coefficients))
      end if
   end function logistic_predict_probability

   subroutine make_design_matrix(raw_x, fit_intercept, design)
      real(dp), intent(in) :: raw_x(:,:)
      logical, intent(in) :: fit_intercept
      real(dp), allocatable, intent(out) :: design(:,:)
      if (fit_intercept) then
         allocate(design(size(raw_x, 1), size(raw_x, 2) + 1))
         design(:, 1) = 1.0_dp
         design(:, 2:) = raw_x
      else
         allocate(design(size(raw_x, 1), size(raw_x, 2)), source=raw_x)
      end if
   end subroutine make_design_matrix

   subroutine logistic_model_initialize(self, n_features, optimizer_kind, fit_intercept, &
                                        lambda, initial_step, initial_coefficients, stat)
      class(stochastic_logistic_t), intent(inout) :: self
      integer, intent(in) :: n_features
      integer, intent(in), optional :: optimizer_kind
      logical, intent(in), optional :: fit_intercept
      real(dp), intent(in), optional :: lambda, initial_step, initial_coefficients(:)
      integer, intent(out), optional :: stat
      integer :: kind, npar, ierr
      real(dp), allocatable :: zero_gradient(:), zero_hessian(:)
      real(dp) :: zero_function

      ierr = 0
      kind = logistic_adaqn
      if (present(optimizer_kind)) kind = optimizer_kind
      self%fit_intercept = .true.
      if (present(fit_intercept)) self%fit_intercept = fit_intercept
      if (n_features <= 0 .or. kind < logistic_olbfgs .or. kind > logistic_adaqn) ierr = 1
      self%lambda = 1.0e-5_dp
      if (present(lambda)) self%lambda = lambda
      self%initial_step = merge(1.0e-3_dp, 1.0e-2_dp, kind == logistic_sqn)
      if (present(initial_step)) self%initial_step = initial_step
      if (self%lambda < 0.0_dp .or. self%initial_step <= 0.0_dp) ierr = 1
      npar = n_features + merge(1, 0, self%fit_intercept)
      if (present(initial_coefficients)) then
         if (size(initial_coefficients) /= npar) ierr = 1
      end if
      if (ierr /= 0) then
         self%initialized = .false.
         if (present(stat)) stat = ierr
         return
      end if

      if (allocated(self%coefficients)) deallocate(self%coefficients)
      allocate(self%coefficients(npar), source=0.0_dp)
      if (present(initial_coefficients)) self%coefficients = initial_coefficients
      self%n_features = n_features
      self%optimizer_kind = kind
      call self%clear_history()

      select case (kind)
      case (logistic_olbfgs)
         call self%olbfgs%initialize(npar, mem_size=10, hessian_initialization=0.0_dp, &
              y_regularization=0.0_dp, min_curvature=1.0e-4_dp, check_finite=.true., stat=ierr)
      case (logistic_sqn)
         call self%sqn%initialize(npar, mem_size=10, bfgs_update_frequency=20, &
              min_curvature=1.0e-4_dp, use_gradient_difference=.false., &
              y_regularization=0.0_dp, check_finite=.true., stat=ierr)
      case (logistic_adaqn)
         call self%adaqn%initialize(npar, mem_size=10, fisher_size=100, &
              bfgs_update_frequency=20, maximum_increase=0.0_dp, min_curvature=1.0e-4_dp, &
              scaling_regularization=1.0e-4_dp, rmsprop_weight=0.9_dp, &
              use_gradient_difference=.false., y_regularization=0.0_dp, &
              check_finite=.true., stat=ierr)
      end select
      if (ierr == 0) then
         allocate(zero_gradient(npar), source=0.0_dp)
         allocate(zero_hessian(npar), source=0.0_dp)
         zero_function = 0.0_dp
         select case (kind)
         case (logistic_olbfgs)
            call self%olbfgs%advance(self%initial_step, self%coefficients, zero_gradient, self%request)
         case (logistic_sqn)
            call self%sqn%advance(self%initial_step, self%coefficients, zero_gradient, &
                                  zero_hessian, self%request)
         case (logistic_adaqn)
            call self%adaqn%advance(self%initial_step, self%coefficients, zero_function, &
                                    zero_gradient, self%request)
         end select
      end if
      self%initialized = ierr == 0
      if (present(stat)) stat = ierr
   end subroutine logistic_model_initialize

   subroutine logistic_model_clear_history(self)
      class(stochastic_logistic_t), intent(inout) :: self
      if (allocated(self%history_x)) deallocate(self%history_x)
      if (allocated(self%history_y)) deallocate(self%history_y)
      if (allocated(self%history_w)) deallocate(self%history_w)
   end subroutine logistic_model_clear_history

   subroutine append_history(self, x, y, weights)
      class(stochastic_logistic_t), intent(inout) :: self
      real(dp), intent(in) :: x(:,:), y(:), weights(:)
      real(dp), allocatable :: new_x(:,:), new_y(:), new_w(:)
      integer :: old_n, new_n

      if (.not. allocated(self%history_y)) then
         allocate(self%history_x(size(x, 1), size(x, 2)), source=x)
         allocate(self%history_y(size(y)), source=y)
         allocate(self%history_w(size(weights)), source=weights)
         return
      end if
      old_n = size(self%history_y)
      new_n = old_n + size(y)
      allocate(new_x(new_n, size(x, 2)))
      allocate(new_y(new_n), new_w(new_n))
      new_x(1:old_n, :) = self%history_x
      new_x(old_n + 1:, :) = x
      new_y(1:old_n) = self%history_y
      new_y(old_n + 1:) = y
      new_w(1:old_n) = self%history_w
      new_w(old_n + 1:) = weights
      call move_alloc(new_x, self%history_x)
      call move_alloc(new_y, self%history_y)
      call move_alloc(new_w, self%history_w)
   end subroutine append_history

   integer function logistic_model_get_iteration(self)
      class(stochastic_logistic_t), intent(in) :: self
      select case (self%optimizer_kind)
      case (logistic_olbfgs)
         logistic_model_get_iteration = self%olbfgs%get_iteration()
      case (logistic_sqn)
         logistic_model_get_iteration = self%sqn%get_iteration()
      case default
         logistic_model_get_iteration = self%adaqn%get_iteration()
      end select
   end function logistic_model_get_iteration

   subroutine logistic_model_partial_fit(self, x, y, weights, validation_x, validation_y, &
                                         validation_weights, info, stat)
      class(stochastic_logistic_t), intent(inout) :: self
      real(dp), intent(in) :: x(:,:), y(:)
      real(dp), intent(in), optional :: weights(:)
      real(dp), intent(in), optional :: validation_x(:,:), validation_y(:), validation_weights(:)
      integer, intent(out), optional :: info, stat
      real(dp), allocatable :: batch_w(:), eval_x(:,:), eval_y(:), eval_w(:), design(:,:)
      real(dp), allocatable :: gradient(:), hessian_vector(:)
      real(dp) :: function_value, step
      integer :: start_iteration, ierr
      logical :: use_validation, special_history_used, updated

      ierr = 0
      if (.not. self%initialized .or. size(x, 1) /= size(y) .or. &
          size(x, 2) /= self%n_features .or. size(y) == 0) ierr = 1
      if (present(weights)) then
         if (size(weights) /= size(y) .or. any(weights < 0.0_dp)) ierr = 1
      end if
      use_validation = present(validation_x) .or. present(validation_y) .or. present(validation_weights)
      if (use_validation) then
         if (.not. present(validation_x) .or. .not. present(validation_y)) ierr = 1
         if (ierr == 0) then
            if (size(validation_x, 1) /= size(validation_y) .or. &
                size(validation_x, 2) /= self%n_features) ierr = 1
            if (present(validation_weights)) then
               if (size(validation_weights) /= size(validation_y)) ierr = 1
            end if
         end if
      end if
      if (ierr /= 0) then
         if (present(stat)) stat = ierr
         if (present(info)) info = info_ok
         return
      end if

      allocate(batch_w(size(y)))
      if (present(weights)) then
         batch_w = weights
      else
         batch_w = 1.0_dp
      end if
      allocate(gradient(size(self%coefficients)), source=0.0_dp)
      allocate(hessian_vector(size(self%coefficients)), source=0.0_dp)
      function_value = 0.0_dp
      start_iteration = self%get_iteration()
      special_history_used = .false.
      updated = .false.

      do
         select case (self%request%task)
         case (task_calc_grad, task_calc_grad_same_batch)
            allocate(eval_x(size(x, 1), size(x, 2)), source=x)
            allocate(eval_y(size(y)), source=y)
            allocate(eval_w(size(batch_w)), source=batch_w)
         case (task_calc_grad_big_batch, task_calc_hess_vec, task_calc_fun_value)
            special_history_used = .true.
            if (use_validation) then
               allocate(eval_x(size(validation_x, 1), size(validation_x, 2)), source=validation_x)
               allocate(eval_y(size(validation_y)), source=validation_y)
               allocate(eval_w(size(validation_y)), source=1.0_dp)
               if (present(validation_weights)) eval_w = validation_weights
            else if (allocated(self%history_y)) then
               allocate(eval_x(size(self%history_x, 1), size(self%history_x, 2)), source=self%history_x)
               allocate(eval_y(size(self%history_y)), source=self%history_y)
               allocate(eval_w(size(self%history_w)), source=self%history_w)
            else
               allocate(eval_x(size(x, 1), size(x, 2)), source=x)
               allocate(eval_y(size(y)), source=y)
               allocate(eval_w(size(batch_w)), source=batch_w)
            end if
         case default
            ierr = 1
            exit
         end select

         call make_design_matrix(eval_x, self%fit_intercept, design)
         select case (self%request%task)
         case (task_calc_grad, task_calc_grad_same_batch, task_calc_grad_big_batch)
            call logistic_gradient(self%request%x, design, eval_y, gradient, self%lambda, eval_w)
         case (task_calc_hess_vec)
            call logistic_hessian_vector(self%request%x, self%request%vec, design, eval_y, &
                                         hessian_vector, self%lambda, eval_w)
         case (task_calc_fun_value)
            function_value = logistic_loss(self%request%x, design, eval_y, self%lambda, eval_w)
         end select
         deallocate(eval_x, eval_y, eval_w, design)

         select case (self%optimizer_kind)
         case (logistic_olbfgs)
            step = self%initial_step / sqrt(real(self%olbfgs%get_iteration(), dp) / 10.0_dp + 1.0_dp)
            call self%olbfgs%advance(step, self%coefficients, gradient, self%request)
         case (logistic_sqn)
            step = self%initial_step / sqrt(real(self%sqn%get_iteration(), dp) / 10.0_dp + 1.0_dp)
            call self%sqn%advance(step, self%coefficients, gradient, hessian_vector, self%request)
         case (logistic_adaqn)
            step = self%initial_step / sqrt(real(self%adaqn%get_iteration(), dp) / 100.0_dp + 1.0_dp)
            call self%adaqn%advance(step, self%coefficients, function_value, gradient, self%request)
         end select
         if (self%request%task == task_invalid) then
            ierr = 1
            exit
         end if
         updated = self%get_iteration() > start_iteration
         if (updated) then
            if (self%optimizer_kind /= logistic_olbfgs) exit
            if (self%request%task == task_calc_grad) exit
         end if
      end do

      if (special_history_used .and. .not. use_validation) call self%clear_history()
      if (self%optimizer_kind /= logistic_olbfgs) call append_history(self, x, y, batch_w)
      if (present(info)) info = self%request%info
      if (present(stat)) stat = ierr
   end subroutine logistic_model_partial_fit

   function logistic_model_get_coefficients(self) result(coefficients)
      class(stochastic_logistic_t), intent(in) :: self
      real(dp), allocatable :: coefficients(:)
      if (allocated(self%coefficients)) then
         allocate(coefficients(size(self%coefficients)), source=self%coefficients)
      else
         allocate(coefficients(0))
      end if
   end function logistic_model_get_coefficients

   function logistic_model_predict_probability(self, x) result(probability)
      class(stochastic_logistic_t), intent(in) :: self
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable :: probability(:)
      real(dp), allocatable :: design(:,:)
      if (.not. self%initialized .or. size(x, 2) /= self%n_features) then
         allocate(probability(size(x, 1)), source=0.0_dp)
         return
      end if
      call make_design_matrix(x, self%fit_intercept, design)
      probability = logistic_predict_probability(self%coefficients, design)
   end function logistic_model_predict_probability

   function logistic_model_predict_class(self, x, threshold) result(classification)
      class(stochastic_logistic_t), intent(in) :: self
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(in), optional :: threshold
      integer, allocatable :: classification(:)
      real(dp), allocatable :: probability(:)
      real(dp) :: cutoff
      cutoff = 0.5_dp
      if (present(threshold)) cutoff = threshold
      probability = self%predict_probability(x)
      allocate(classification(size(probability)))
      classification = merge(1, 0, probability >= cutoff)
   end function logistic_model_predict_class

end module stochqn_logistic
