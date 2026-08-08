module stochqn_core
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use stochqn_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: task_invalid = 100
   integer, parameter, public :: task_calc_grad = 101
   integer, parameter, public :: task_calc_grad_same_batch = 102
   integer, parameter, public :: task_calc_grad_big_batch = 103
   integer, parameter, public :: task_calc_hess_vec = 104
   integer, parameter, public :: task_calc_fun_value = 105

   integer, parameter, public :: info_ok = 200
   integer, parameter, public :: info_function_increased = 201
   integer, parameter, public :: info_curvature_too_small = 202
   integer, parameter, public :: info_invalid_direction = 203

   integer, parameter, public :: status_no_update = 0
   integer, parameter, public :: status_updated = 1
   integer, parameter, public :: status_invalid_input = -1000

   type, public :: stochqn_request_t
      integer :: task = task_invalid
      integer :: info = info_ok
      integer :: status = status_no_update
      integer :: iteration = 0
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: vec(:)
   contains
      procedure :: clear => request_clear
   end type stochqn_request_t

   type :: bfgs_memory_t
      real(dp), allocatable :: s(:,:)
      real(dp), allocatable :: y(:,:)
      real(dp), allocatable :: rho(:)
      real(dp), allocatable :: alpha(:)
      real(dp), allocatable :: s_backup(:)
      real(dp), allocatable :: y_backup(:)
      integer :: mem_size = 0
      integer :: mem_used = 0
      integer :: next_slot = 1
      integer :: update_frequency = 1
      real(dp) :: y_regularization = 0.0_dp
      real(dp) :: min_curvature = 0.0_dp
   contains
      procedure :: initialize => bfgs_initialize
      procedure :: reset => bfgs_reset
   end type bfgs_memory_t

   type :: fisher_memory_t
      real(dp), allocatable :: gradients(:,:)
      integer :: mem_size = 0
      integer :: mem_used = 0
      integer :: next_slot = 1
   contains
      procedure :: initialize => fisher_initialize
      procedure :: reset => fisher_reset
      procedure :: add => fisher_add
   end type fisher_memory_t

   type, public :: olbfgs_t
      private
      type(bfgs_memory_t) :: memory
      real(dp), allocatable :: previous_gradient(:)
      real(dp) :: hessian_initialization = 0.0_dp
      integer :: iteration = 0
      integer :: section = 0
      logical :: check_finite = .true.
      integer :: n = 0
      logical :: initialized = .false.
   contains
      procedure, public :: initialize => olbfgs_initialize
      procedure, public :: advance => olbfgs_advance
      procedure, public :: reset => olbfgs_reset
      procedure, public :: get_iteration => olbfgs_get_iteration
      procedure, public :: correction_pairs => olbfgs_correction_pairs
   end type olbfgs_t

   type, public :: sqn_t
      private
      type(bfgs_memory_t) :: memory
      real(dp), allocatable :: previous_gradient(:)
      real(dp), allocatable :: x_sum(:)
      real(dp), allocatable :: previous_average(:)
      logical :: use_gradient_difference = .false.
      integer :: iteration = 0
      integer :: section = 0
      logical :: check_finite = .true.
      integer :: n = 0
      logical :: initialized = .false.
   contains
      procedure, public :: initialize => sqn_initialize
      procedure, public :: advance => sqn_advance
      procedure, public :: reset => sqn_reset
      procedure, public :: get_iteration => sqn_get_iteration
      procedure, public :: correction_pairs => sqn_correction_pairs
   end type sqn_t

   type, public :: adaqn_t
      private
      type(bfgs_memory_t) :: memory
      type(fisher_memory_t) :: fisher
      real(dp), allocatable :: h0(:)
      real(dp), allocatable :: previous_gradient(:)
      real(dp), allocatable :: x_sum(:)
      real(dp), allocatable :: previous_average(:)
      real(dp), allocatable :: gradient_sum_squares(:)
      real(dp) :: previous_function = 0.0_dp
      real(dp) :: maximum_increase = 1.01_dp
      real(dp) :: scaling_regularization = 1.0e-4_dp
      real(dp) :: rmsprop_weight = 0.9_dp
      logical :: use_gradient_difference = .false.
      integer :: iteration = 0
      integer :: section = 0
      logical :: check_finite = .true.
      integer :: n = 0
      logical :: initialized = .false.
   contains
      procedure, public :: initialize => adaqn_initialize
      procedure, public :: advance => adaqn_advance
      procedure, public :: reset => adaqn_reset
      procedure, public :: get_iteration => adaqn_get_iteration
      procedure, public :: correction_pairs => adaqn_correction_pairs
      procedure, public :: set_previous_function => adaqn_set_previous_function
   end type adaqn_t

   public :: stochqn_task_name, stochqn_info_name, stochqn_status_name

contains

   subroutine request_clear(self)
      class(stochqn_request_t), intent(inout) :: self
      self%task = task_invalid
      self%info = info_ok
      self%status = status_no_update
      self%iteration = 0
      if (allocated(self%x)) deallocate(self%x)
      if (allocated(self%vec)) deallocate(self%vec)
   end subroutine request_clear

   subroutine set_request(request, task, x, iteration, info, status, vec)
      type(stochqn_request_t), intent(inout) :: request
      integer, intent(in) :: task, iteration, info, status
      real(dp), intent(in) :: x(:)
      real(dp), intent(in), optional :: vec(:)

      call request%clear()
      request%task = task
      request%iteration = iteration
      request%info = info
      request%status = status
      allocate(request%x(size(x)), source=x)
      if (present(vec)) allocate(request%vec(size(vec)), source=vec)
   end subroutine set_request

   subroutine set_invalid_request(request, x, iteration)
      type(stochqn_request_t), intent(inout) :: request
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: iteration
      call set_request(request, task_invalid, x, iteration, info_ok, status_invalid_input)
   end subroutine set_invalid_request

   subroutine bfgs_initialize(self, n, mem_size, update_frequency, min_curvature, y_regularization, stat)
      class(bfgs_memory_t), intent(inout) :: self
      integer, intent(in) :: n, mem_size, update_frequency
      real(dp), intent(in) :: min_curvature, y_regularization
      integer, intent(out) :: stat

      stat = 0
      if (n <= 0 .or. mem_size <= 0 .or. update_frequency <= 0) then
         stat = 1
         return
      end if
      if (min_curvature < 0.0_dp .or. y_regularization < 0.0_dp) then
         stat = 1
         return
      end if

      if (allocated(self%s)) deallocate(self%s, self%y, self%rho, self%alpha, &
                                        self%s_backup, self%y_backup)
      allocate(self%s(n, mem_size), self%y(n, mem_size))
      allocate(self%rho(mem_size), self%alpha(mem_size))
      allocate(self%s_backup(n), self%y_backup(n))
      self%s = 0.0_dp
      self%y = 0.0_dp
      self%rho = 0.0_dp
      self%alpha = 0.0_dp
      self%s_backup = 0.0_dp
      self%y_backup = 0.0_dp
      self%mem_size = mem_size
      self%update_frequency = update_frequency
      self%min_curvature = min_curvature
      self%y_regularization = y_regularization
      call self%reset()
   end subroutine bfgs_initialize

   subroutine bfgs_reset(self)
      class(bfgs_memory_t), intent(inout) :: self
      self%mem_used = 0
      self%next_slot = 1
   end subroutine bfgs_reset

   subroutine fisher_initialize(self, n, mem_size, stat)
      class(fisher_memory_t), intent(inout) :: self
      integer, intent(in) :: n, mem_size
      integer, intent(out) :: stat

      stat = 0
      if (n <= 0 .or. mem_size <= 0) then
         stat = 1
         return
      end if
      if (allocated(self%gradients)) deallocate(self%gradients)
      allocate(self%gradients(n, mem_size))
      self%gradients = 0.0_dp
      self%mem_size = mem_size
      call self%reset()
   end subroutine fisher_initialize

   subroutine fisher_reset(self)
      class(fisher_memory_t), intent(inout) :: self
      self%mem_used = 0
      self%next_slot = 1
   end subroutine fisher_reset

   subroutine fisher_add(self, gradient)
      class(fisher_memory_t), intent(inout) :: self
      real(dp), intent(in) :: gradient(:)
      if (.not. allocated(self%gradients)) return
      if (size(gradient) /= size(self%gradients, 1)) return
      self%gradients(:, self%next_slot) = gradient
      self%next_slot = modulo(self%next_slot, self%mem_size) + 1
      self%mem_used = min(self%mem_used + 1, self%mem_size)
   end subroutine fisher_add

   pure logical function vector_is_finite(x)
      real(dp), intent(in) :: x(:)
      vector_is_finite = all(ieee_is_finite(x))
   end function vector_is_finite

   pure integer function earliest_slot(memory)
      type(bfgs_memory_t), intent(in) :: memory
      if (memory%mem_used == memory%mem_size) then
         earliest_slot = memory%next_slot
      else
         earliest_slot = 1
      end if
   end function earliest_slot

   subroutine inverse_hessian_product(memory, gradient, direction, h0_scalar, h0_diagonal)
      type(bfgs_memory_t), intent(inout) :: memory
      real(dp), intent(in) :: gradient(:)
      real(dp), intent(out) :: direction(:)
      real(dp), intent(in) :: h0_scalar
      real(dp), intent(in), optional :: h0_diagonal(:)
      integer :: k, slot, first, newest
      real(dp) :: beta, yy, sy, scale

      direction = gradient
      first = earliest_slot(memory)

      do k = memory%mem_used, 1, -1
         slot = modulo(first - 1 + k - 1, memory%mem_size) + 1
         sy = dot_product(memory%y(:, slot), memory%s(:, slot))
         memory%rho(k) = 1.0_dp / sy
         memory%alpha(k) = memory%rho(k) * dot_product(memory%s(:, slot), direction)
         direction = direction - memory%alpha(k) * memory%y(:, slot)
      end do

      if (present(h0_diagonal)) then
         direction = h0_diagonal * direction
      else if (h0_scalar > 0.0_dp) then
         direction = h0_scalar * direction
      else
         newest = modulo(first - 1 + memory%mem_used - 1, memory%mem_size) + 1
         sy = dot_product(memory%s(:, newest), memory%y(:, newest))
         yy = dot_product(memory%y(:, newest), memory%y(:, newest))
         scale = sy / yy
         direction = scale * direction
      end if

      do k = 1, memory%mem_used
         slot = modulo(first - 1 + k - 1, memory%mem_size) + 1
         beta = memory%rho(k) * dot_product(memory%y(:, slot), direction)
         direction = direction + (memory%alpha(k) - beta) * memory%s(:, slot)
      end do
   end subroutine inverse_hessian_product

   subroutine update_gradient_squares(gradient, gradient_sum_squares, rmsprop_weight)
      real(dp), intent(in) :: gradient(:)
      real(dp), intent(inout) :: gradient_sum_squares(:)
      real(dp), intent(in) :: rmsprop_weight
      if (rmsprop_weight > 0.0_dp .and. rmsprop_weight < 1.0_dp) then
         gradient_sum_squares = rmsprop_weight * gradient_sum_squares + &
              (1.0_dp - rmsprop_weight) * gradient * gradient
      else
         gradient_sum_squares = gradient_sum_squares + gradient * gradient
      end if
   end subroutine update_gradient_squares

   subroutine take_step(memory, x, gradient, step_size, check_finite, direction, ok, &
                        h0_scalar, gradient_sum_squares, scaling_regularization, &
                        rmsprop_weight, h0_diagonal)
      type(bfgs_memory_t), intent(inout) :: memory
      real(dp), intent(inout) :: x(:)
      real(dp), intent(in) :: gradient(:)
      real(dp), intent(in) :: step_size
      logical, intent(in) :: check_finite
      real(dp), intent(out) :: direction(:)
      logical, intent(out) :: ok
      real(dp), intent(in), optional :: h0_scalar
      real(dp), intent(inout), optional :: gradient_sum_squares(:)
      real(dp), intent(in), optional :: scaling_regularization, rmsprop_weight
      real(dp), intent(out), optional :: h0_diagonal(:)
      real(dp) :: h0, eps_scal, rms

      h0 = 0.0_dp
      if (present(h0_scalar)) h0 = h0_scalar
      eps_scal = 0.0_dp
      if (present(scaling_regularization)) eps_scal = scaling_regularization
      rms = 0.0_dp
      if (present(rmsprop_weight)) rms = rmsprop_weight

      if (present(gradient_sum_squares)) then
         call update_gradient_squares(gradient, gradient_sum_squares, rms)
      end if

      if (memory%mem_used == 0) then
         direction = gradient
         if (present(gradient_sum_squares)) then
            direction = direction / sqrt(gradient_sum_squares + eps_scal)
         end if
      else
         if (present(gradient_sum_squares)) then
            if (.not. present(h0_diagonal)) then
               ok = .false.
               return
            end if
            h0_diagonal = 1.0_dp / sqrt(gradient_sum_squares + eps_scal)
            call inverse_hessian_product(memory, gradient, direction, h0, h0_diagonal)
         else
            call inverse_hessian_product(memory, gradient, direction, h0)
         end if
      end if

      ok = .true.
      if (check_finite) then
         if (.not. vector_is_finite(direction)) ok = .false.
         if (norm2(direction) > 1.0e3_dp * real(size(direction), dp)) ok = .false.
      end if
      if (.not. ok) then
         call memory%reset()
         return
      end if
      x = x - step_size * direction
   end subroutine take_step

   subroutine backup_slot(memory)
      type(bfgs_memory_t), intent(inout) :: memory
      memory%s_backup = memory%s(:, memory%next_slot)
      memory%y_backup = memory%y(:, memory%next_slot)
   end subroutine backup_slot

   subroutine restore_slot(memory)
      type(bfgs_memory_t), intent(inout) :: memory
      memory%s(:, memory%next_slot) = memory%s_backup
      memory%y(:, memory%next_slot) = memory%y_backup
   end subroutine restore_slot

   subroutine accept_or_reject_pair(memory, info)
      type(bfgs_memory_t), intent(inout) :: memory
      integer, intent(out) :: info
      real(dp) :: ss, sy, curvature, threshold

      ss = dot_product(memory%s(:, memory%next_slot), memory%s(:, memory%next_slot))
      sy = dot_product(memory%s(:, memory%next_slot), memory%y(:, memory%next_slot))
      threshold = max(memory%min_curvature, 10.0_dp * epsilon(1.0_dp))
      info = info_ok
      if (.not. ieee_is_finite(ss) .or. .not. ieee_is_finite(sy) .or. ss <= tiny(1.0_dp)) then
         call restore_slot(memory)
         info = info_curvature_too_small
         return
      end if
      curvature = sy / ss
      if (curvature <= threshold) then
         call restore_slot(memory)
         info = info_curvature_too_small
         return
      end if
      memory%next_slot = modulo(memory%next_slot, memory%mem_size) + 1
      memory%mem_used = min(memory%mem_used + 1, memory%mem_size)
   end subroutine accept_or_reject_pair

   subroutine set_s_from_average(memory, average, previous_average)
      type(bfgs_memory_t), intent(inout) :: memory
      real(dp), intent(in) :: average(:), previous_average(:)
      call backup_slot(memory)
      memory%s(:, memory%next_slot) = average - previous_average
   end subroutine set_s_from_average

   subroutine set_y_gradient_difference(memory, gradient, previous_gradient, info)
      type(bfgs_memory_t), intent(inout) :: memory
      real(dp), intent(in) :: gradient(:), previous_gradient(:)
      integer, intent(out) :: info
      memory%y(:, memory%next_slot) = gradient - previous_gradient + &
           memory%y_regularization * memory%s(:, memory%next_slot)
      call accept_or_reject_pair(memory, info)
   end subroutine set_y_gradient_difference

   subroutine set_y_hessian_vector(memory, hessian_vector, info)
      type(bfgs_memory_t), intent(inout) :: memory
      real(dp), intent(in) :: hessian_vector(:)
      integer, intent(out) :: info
      memory%y(:, memory%next_slot) = hessian_vector
      call accept_or_reject_pair(memory, info)
   end subroutine set_y_hessian_vector

   subroutine set_y_fisher(memory, fisher, info)
      type(bfgs_memory_t), intent(inout) :: memory
      type(fisher_memory_t), intent(in) :: fisher
      integer, intent(out) :: info
      integer :: j
      real(dp), allocatable :: y(:)
      real(dp) :: projection

      allocate(y(size(memory%s, 1)), source=0.0_dp)
      do j = 1, fisher%mem_used
         projection = dot_product(fisher%gradients(:, j), memory%s(:, memory%next_slot))
         y = y + projection * fisher%gradients(:, j)
      end do
      if (fisher%mem_used > 0) y = y / real(fisher%mem_used, dp)
      memory%y(:, memory%next_slot) = y
      call accept_or_reject_pair(memory, info)
   end subroutine set_y_fisher

   subroutine average_and_clear(sum_vector, count, average)
      real(dp), intent(inout) :: sum_vector(:)
      integer, intent(in) :: count
      real(dp), intent(out) :: average(:)
      average = sum_vector / real(count, dp)
      sum_vector = 0.0_dp
   end subroutine average_and_clear

   subroutine olbfgs_initialize(self, n, mem_size, hessian_initialization, y_regularization, &
                                min_curvature, check_finite, stat)
      class(olbfgs_t), intent(inout) :: self
      integer, intent(in) :: n
      integer, intent(in), optional :: mem_size
      real(dp), intent(in), optional :: hessian_initialization, y_regularization, min_curvature
      logical, intent(in), optional :: check_finite
      integer, intent(out), optional :: stat
      integer :: m, ierr
      real(dp) :: h0, yr, curv

      m = 10
      if (present(mem_size)) m = mem_size
      h0 = 0.0_dp
      if (present(hessian_initialization)) h0 = hessian_initialization
      yr = 0.0_dp
      if (present(y_regularization)) yr = y_regularization
      curv = 0.0_dp
      if (present(min_curvature)) curv = min_curvature
      if (h0 < 0.0_dp) then
         ierr = 1
      else
         call self%memory%initialize(n, m, 1, curv, yr, ierr)
      end if
      if (ierr == 0) then
         if (allocated(self%previous_gradient)) deallocate(self%previous_gradient)
         allocate(self%previous_gradient(n), source=0.0_dp)
         self%n = n
         self%hessian_initialization = h0
         self%check_finite = .true.
         if (present(check_finite)) self%check_finite = check_finite
         self%initialized = .true.
         call self%reset()
      else
         self%initialized = .false.
      end if
      if (present(stat)) stat = ierr
   end subroutine olbfgs_initialize

   subroutine olbfgs_reset(self)
      class(olbfgs_t), intent(inout) :: self
      self%iteration = 0
      self%section = 0
      if (allocated(self%previous_gradient)) self%previous_gradient = 0.0_dp
      call self%memory%reset()
   end subroutine olbfgs_reset

   subroutine olbfgs_advance(self, step_size, x, gradient, request)
      class(olbfgs_t), intent(inout) :: self
      real(dp), intent(in) :: step_size
      real(dp), intent(inout) :: x(:)
      real(dp), intent(in) :: gradient(:)
      type(stochqn_request_t), intent(inout) :: request
      real(dp), allocatable :: direction(:)
      logical :: ok
      integer :: info

      if (.not. self%initialized .or. size(x) /= self%n .or. size(gradient) /= self%n .or. &
          step_size <= 0.0_dp) then
         call set_invalid_request(request, x, self%iteration)
         return
      end if

      select case (self%section)
      case (0)
         self%section = 1
         call set_request(request, task_calc_grad, x, self%iteration, info_ok, status_no_update)
      case (1)
         self%previous_gradient = gradient
         allocate(direction(self%n))
         call take_step(self%memory, x, gradient, step_size, self%check_finite, direction, ok, &
                        h0_scalar=self%hessian_initialization)
         self%iteration = self%iteration + 1
         if (ok) then
            call backup_slot(self%memory)
            self%memory%s(:, self%memory%next_slot) = -step_size * direction
            self%section = 2
            call set_request(request, task_calc_grad_same_batch, x, self%iteration, info_ok, status_updated)
         else
            self%section = 1
            call set_request(request, task_calc_grad, x, self%iteration, info_invalid_direction, status_no_update)
         end if
      case (2)
         call set_y_gradient_difference(self%memory, gradient, self%previous_gradient, info)
         self%section = 1
         call set_request(request, task_calc_grad, x, self%iteration, info, status_no_update)
      case default
         call set_invalid_request(request, x, self%iteration)
      end select
   end subroutine olbfgs_advance

   integer function olbfgs_get_iteration(self)
      class(olbfgs_t), intent(in) :: self
      olbfgs_get_iteration = self%iteration
   end function olbfgs_get_iteration

   integer function olbfgs_correction_pairs(self)
      class(olbfgs_t), intent(in) :: self
      olbfgs_correction_pairs = self%memory%mem_used
   end function olbfgs_correction_pairs

   subroutine sqn_initialize(self, n, mem_size, bfgs_update_frequency, min_curvature, &
                             use_gradient_difference, y_regularization, check_finite, stat)
      class(sqn_t), intent(inout) :: self
      integer, intent(in) :: n
      integer, intent(in), optional :: mem_size, bfgs_update_frequency
      real(dp), intent(in), optional :: min_curvature, y_regularization
      logical, intent(in), optional :: use_gradient_difference, check_finite
      integer, intent(out), optional :: stat
      integer :: m, freq, ierr
      real(dp) :: curv, yr

      m = 10
      if (present(mem_size)) m = mem_size
      freq = 20
      if (present(bfgs_update_frequency)) freq = bfgs_update_frequency
      curv = 1.0e-4_dp
      if (present(min_curvature)) curv = min_curvature
      yr = 0.0_dp
      if (present(y_regularization)) yr = y_regularization
      call self%memory%initialize(n, m, freq, curv, yr, ierr)
      if (ierr == 0) then
         if (allocated(self%previous_gradient)) then
            deallocate(self%previous_gradient, self%x_sum, self%previous_average)
         end if
         allocate(self%previous_gradient(n), source=0.0_dp)
         allocate(self%x_sum(n), source=0.0_dp)
         allocate(self%previous_average(n), source=0.0_dp)
         self%n = n
         self%use_gradient_difference = .false.
         if (present(use_gradient_difference)) self%use_gradient_difference = use_gradient_difference
         self%check_finite = .true.
         if (present(check_finite)) self%check_finite = check_finite
         self%initialized = .true.
         call self%reset()
      else
         self%initialized = .false.
      end if
      if (present(stat)) stat = ierr
   end subroutine sqn_initialize

   subroutine sqn_reset(self)
      class(sqn_t), intent(inout) :: self
      self%iteration = 0
      self%section = 0
      if (allocated(self%previous_gradient)) self%previous_gradient = 0.0_dp
      if (allocated(self%x_sum)) self%x_sum = 0.0_dp
      if (allocated(self%previous_average)) self%previous_average = 0.0_dp
      call self%memory%reset()
   end subroutine sqn_reset

   subroutine sqn_regular_request(self, x, request, info, status)
      class(sqn_t), intent(inout) :: self
      real(dp), intent(in) :: x(:)
      type(stochqn_request_t), intent(inout) :: request
      integer, intent(in) :: info, status
      self%section = 1
      call set_request(request, task_calc_grad, x, self%iteration, info, status)
   end subroutine sqn_regular_request

   subroutine sqn_advance(self, step_size, x, gradient, hessian_vector, request)
      class(sqn_t), intent(inout) :: self
      real(dp), intent(in) :: step_size
      real(dp), intent(inout) :: x(:)
      real(dp), intent(in) :: gradient(:)
      real(dp), intent(in) :: hessian_vector(:)
      type(stochqn_request_t), intent(inout) :: request
      real(dp), allocatable :: direction(:), average(:)
      logical :: ok
      integer :: info, status

      if (.not. self%initialized .or. size(x) /= self%n .or. size(gradient) /= self%n .or. &
          size(hessian_vector) /= self%n .or. step_size <= 0.0_dp) then
         call set_invalid_request(request, x, self%iteration)
         return
      end if
      allocate(average(self%n))

      select case (self%section)
      case (0)
         call sqn_regular_request(self, x, request, info_ok, status_no_update)
      case (1)
         allocate(direction(self%n))
         call take_step(self%memory, x, gradient, step_size, self%check_finite, direction, ok)
         self%iteration = self%iteration + 1
         status = merge(status_updated, status_no_update, ok)
         info = merge(info_ok, info_invalid_direction, ok)
         self%x_sum = self%x_sum + x

         if (modulo(self%iteration, self%memory%update_frequency) /= 0) then
            call sqn_regular_request(self, x, request, info, status)
            return
         end if

         if (self%iteration == self%memory%update_frequency) then
            call average_and_clear(self%x_sum, self%memory%update_frequency, average)
            self%previous_average = average
            if (self%use_gradient_difference) then
               self%section = 2
               call set_request(request, task_calc_grad_big_batch, self%previous_average, &
                                self%iteration, info, status)
            else
               call sqn_regular_request(self, x, request, info, status)
            end if
            return
         end if

         average = self%x_sum / real(self%memory%update_frequency, dp)
         call set_s_from_average(self%memory, average, self%previous_average)
         if (self%use_gradient_difference) then
            self%section = 3
            call set_request(request, task_calc_grad_big_batch, average, self%iteration, info, status)
         else
            self%section = 4
            call set_request(request, task_calc_hess_vec, average, self%iteration, info, status, &
                             self%memory%s(:, self%memory%next_slot))
         end if
      case (2)
         self%previous_gradient = gradient
         call sqn_regular_request(self, x, request, info_ok, status_no_update)
      case (3)
         average = self%x_sum / real(self%memory%update_frequency, dp)
         call set_y_gradient_difference(self%memory, gradient, self%previous_gradient, info)
         if (info == info_ok) then
            self%previous_gradient = gradient
            self%previous_average = average
         end if
         self%x_sum = 0.0_dp
         call sqn_regular_request(self, x, request, info, status_no_update)
      case (4)
         average = self%x_sum / real(self%memory%update_frequency, dp)
         self%previous_average = average
         self%x_sum = 0.0_dp
         call set_y_hessian_vector(self%memory, hessian_vector, info)
         call sqn_regular_request(self, x, request, info, status_no_update)
      case default
         call set_invalid_request(request, x, self%iteration)
      end select
   end subroutine sqn_advance

   integer function sqn_get_iteration(self)
      class(sqn_t), intent(in) :: self
      sqn_get_iteration = self%iteration
   end function sqn_get_iteration

   integer function sqn_correction_pairs(self)
      class(sqn_t), intent(in) :: self
      sqn_correction_pairs = self%memory%mem_used
   end function sqn_correction_pairs

   subroutine adaqn_initialize(self, n, mem_size, fisher_size, bfgs_update_frequency, &
                               maximum_increase, min_curvature, scaling_regularization, &
                               rmsprop_weight, use_gradient_difference, y_regularization, &
                               check_finite, stat)
      class(adaqn_t), intent(inout) :: self
      integer, intent(in) :: n
      integer, intent(in), optional :: mem_size, fisher_size, bfgs_update_frequency
      real(dp), intent(in), optional :: maximum_increase, min_curvature
      real(dp), intent(in), optional :: scaling_regularization, rmsprop_weight, y_regularization
      logical, intent(in), optional :: use_gradient_difference, check_finite
      integer, intent(out), optional :: stat
      integer :: m, fsize, freq, ierr, ierr2
      real(dp) :: maxinc, curv, scal, rms, yr
      logical :: grad_diff

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

      ierr = 0
      if (maxinc < 0.0_dp .or. scal <= 0.0_dp .or. rms < 0.0_dp .or. rms >= 1.0_dp) ierr = 1
      if (.not. grad_diff .and. fsize <= 0) ierr = 1
      if (ierr == 0) call self%memory%initialize(n, m, freq, curv, yr, ierr)
      if (ierr == 0 .and. .not. grad_diff) then
         call self%fisher%initialize(n, fsize, ierr2)
         if (ierr2 /= 0) ierr = ierr2
      end if

      if (ierr == 0) then
         if (allocated(self%h0)) then
            deallocate(self%h0, self%previous_gradient, self%x_sum, self%previous_average, &
                       self%gradient_sum_squares)
         end if
         allocate(self%h0(n), source=0.0_dp)
         allocate(self%previous_gradient(n), source=0.0_dp)
         allocate(self%x_sum(n), source=0.0_dp)
         allocate(self%previous_average(n), source=0.0_dp)
         allocate(self%gradient_sum_squares(n), source=0.0_dp)
         self%n = n
         self%maximum_increase = maxinc
         self%scaling_regularization = scal
         self%rmsprop_weight = rms
         self%use_gradient_difference = grad_diff
         self%check_finite = .true.
         if (present(check_finite)) self%check_finite = check_finite
         self%initialized = .true.
         call self%reset()
      else
         self%initialized = .false.
      end if
      if (present(stat)) stat = ierr
   end subroutine adaqn_initialize

   subroutine adaqn_reset(self)
      class(adaqn_t), intent(inout) :: self
      self%iteration = 0
      self%section = 0
      self%previous_function = 0.0_dp
      if (allocated(self%h0)) self%h0 = 0.0_dp
      if (allocated(self%previous_gradient)) self%previous_gradient = 0.0_dp
      if (allocated(self%x_sum)) self%x_sum = 0.0_dp
      if (allocated(self%previous_average)) self%previous_average = 0.0_dp
      if (allocated(self%gradient_sum_squares)) self%gradient_sum_squares = 0.0_dp
      call self%memory%reset()
      call self%fisher%reset()
   end subroutine adaqn_reset

   subroutine adaqn_regular_request(self, x, request, info, status)
      class(adaqn_t), intent(inout) :: self
      real(dp), intent(in) :: x(:)
      type(stochqn_request_t), intent(inout) :: request
      integer, intent(in) :: info, status
      self%section = 1
      call set_request(request, task_calc_grad, x, self%iteration, info, status)
   end subroutine adaqn_regular_request

   subroutine adaqn_update_y(self, x, average, request, info, status)
      class(adaqn_t), intent(inout) :: self
      real(dp), intent(in) :: x(:), average(:)
      type(stochqn_request_t), intent(inout) :: request
      integer, intent(in) :: info, status
      integer :: pair_info

      if (self%use_gradient_difference) then
         self%section = 4
         call set_request(request, task_calc_grad_big_batch, average, self%iteration, info, status)
      else
         call set_y_fisher(self%memory, self%fisher, pair_info)
         if (pair_info == info_ok) self%previous_average = average
         self%x_sum = 0.0_dp
         call adaqn_regular_request(self, x, request, pair_info, status_no_update)
      end if
   end subroutine adaqn_update_y

   subroutine adaqn_advance(self, step_size, x, function_value, gradient, request)
      class(adaqn_t), intent(inout) :: self
      real(dp), intent(in) :: step_size, function_value
      real(dp), intent(inout) :: x(:)
      real(dp), intent(in) :: gradient(:)
      type(stochqn_request_t), intent(inout) :: request
      real(dp), allocatable :: direction(:), average(:)
      logical :: ok
      integer :: info, status

      if (.not. self%initialized .or. size(x) /= self%n .or. size(gradient) /= self%n .or. &
          step_size <= 0.0_dp) then
         call set_invalid_request(request, x, self%iteration)
         return
      end if
      allocate(average(self%n))

      select case (self%section)
      case (0)
         call adaqn_regular_request(self, x, request, info_ok, status_no_update)
      case (1)
         if (.not. self%use_gradient_difference) call self%fisher%add(gradient)
         allocate(direction(self%n))
         call take_step(self%memory, x, gradient, step_size, self%check_finite, direction, ok, &
                        gradient_sum_squares=self%gradient_sum_squares, &
                        scaling_regularization=self%scaling_regularization, &
                        rmsprop_weight=self%rmsprop_weight, h0_diagonal=self%h0)
         self%iteration = self%iteration + 1
         status = merge(status_updated, status_no_update, ok)
         info = merge(info_ok, info_invalid_direction, ok)
         self%x_sum = self%x_sum + x

         if (modulo(self%iteration, self%memory%update_frequency) /= 0) then
            call adaqn_regular_request(self, x, request, info, status)
            return
         end if

         if (self%iteration == self%memory%update_frequency) then
            call average_and_clear(self%x_sum, self%memory%update_frequency, average)
            self%previous_average = average
            if (self%use_gradient_difference) then
               self%section = 2
               call set_request(request, task_calc_grad_big_batch, average, self%iteration, info, status)
            else if (self%maximum_increase > 0.0_dp) then
               self%section = 3
               call set_request(request, task_calc_fun_value, average, self%iteration, info, status)
            else
               call adaqn_regular_request(self, x, request, info, status)
            end if
            return
         end if

         average = self%x_sum / real(self%memory%update_frequency, dp)
         if (self%maximum_increase > 0.0_dp) then
            self%section = 5
            call set_request(request, task_calc_fun_value, average, self%iteration, info, status)
            return
         end if
         call set_s_from_average(self%memory, average, self%previous_average)
         call adaqn_update_y(self, x, average, request, info, status)
      case (2)
         self%previous_gradient = gradient
         if (self%maximum_increase > 0.0_dp) then
            self%section = 3
            call set_request(request, task_calc_fun_value, self%previous_average, &
                             self%iteration, info_ok, status_no_update)
         else
            call adaqn_regular_request(self, x, request, info_ok, status_no_update)
         end if
      case (3)
         self%previous_function = function_value
         call adaqn_regular_request(self, x, request, info_ok, status_no_update)
      case (4)
         average = self%x_sum / real(self%memory%update_frequency, dp)
         call set_y_gradient_difference(self%memory, gradient, self%previous_gradient, info)
         if (info == info_ok) then
            self%previous_gradient = gradient
            self%previous_average = average
         end if
         self%x_sum = 0.0_dp
         call adaqn_regular_request(self, x, request, info, status_no_update)
      case (5)
         average = self%x_sum / real(self%memory%update_frequency, dp)
         if (.not. ieee_is_finite(function_value) .or. &
             function_value > self%maximum_increase * self%previous_function) then
            call self%memory%reset()
            call self%fisher%reset()
            x = self%previous_average
            self%x_sum = 0.0_dp
            call adaqn_regular_request(self, x, request, info_function_increased, status_updated)
         else
            self%previous_function = function_value
            call set_s_from_average(self%memory, average, self%previous_average)
            call adaqn_update_y(self, x, average, request, info_ok, status_no_update)
         end if
      case default
         call set_invalid_request(request, x, self%iteration)
      end select
   end subroutine adaqn_advance

   integer function adaqn_get_iteration(self)
      class(adaqn_t), intent(in) :: self
      adaqn_get_iteration = self%iteration
   end function adaqn_get_iteration

   integer function adaqn_correction_pairs(self)
      class(adaqn_t), intent(in) :: self
      adaqn_correction_pairs = self%memory%mem_used
   end function adaqn_correction_pairs

   subroutine adaqn_set_previous_function(self, value)
      class(adaqn_t), intent(inout) :: self
      real(dp), intent(in) :: value
      self%previous_function = value
   end subroutine adaqn_set_previous_function

   pure function stochqn_task_name(task) result(name)
      integer, intent(in) :: task
      character(len=:), allocatable :: name
      select case (task)
      case (task_calc_grad)
         name = 'calc_grad'
      case (task_calc_grad_same_batch)
         name = 'calc_grad_same_batch'
      case (task_calc_grad_big_batch)
         name = 'calc_grad_big_batch'
      case (task_calc_hess_vec)
         name = 'calc_hess_vec'
      case (task_calc_fun_value)
         name = 'calc_fun_value_batch'
      case default
         name = 'invalid_input'
      end select
   end function stochqn_task_name

   pure function stochqn_info_name(info) result(name)
      integer, intent(in) :: info
      character(len=:), allocatable :: name
      select case (info)
      case (info_ok)
         name = 'no_problems_encountered'
      case (info_function_increased)
         name = 'function_increased'
      case (info_curvature_too_small)
         name = 'curvature_too_small'
      case (info_invalid_direction)
         name = 'search_direction_invalid'
      case default
         name = 'unknown_info'
      end select
   end function stochqn_info_name

   pure function stochqn_status_name(status) result(name)
      integer, intent(in) :: status
      character(len=:), allocatable :: name
      select case (status)
      case (status_no_update)
         name = 'did_not_update_x'
      case (status_updated)
         name = 'updated_x'
      case default
         name = 'received_invalid_input'
      end select
   end function stochqn_status_name

end module stochqn_core
