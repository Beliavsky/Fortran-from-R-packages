! SPDX-License-Identifier: Artistic-2.0
module ldhmm_parameters
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use ldhmm_kinds, only : dp
   use ldhmm_math, only : solve_linear_system
   use ldhmm_status, only : LDHMM_SUCCESS, LDHMM_INVALID_ARGUMENT, LDHMM_NUMERICAL_ERROR
   use ldhmm_types, only : ldhmm_model
   implicit none
   private

   public :: ldhmm_create, ldhmm_validate, ldhmm_gamma_init
   public :: ldhmm_natural_to_working, ldhmm_working_to_natural
   public :: ldhmm_stationary_distribution, ldhmm_parameter_count

contains

   function ldhmm_create(m, param, gamma_matrix, delta, stationary, optimizer, status) result(model)
      integer, intent(in) :: m
      real(dp), intent(in) :: param(:, :), gamma_matrix(:, :)
      real(dp), intent(in), optional :: delta(:)
      logical, intent(in), optional :: stationary
      character(len=*), intent(in), optional :: optimizer
      integer, intent(out), optional :: status
      type(ldhmm_model) :: model
      integer :: local_status

      model%m = m
      model%param_nbr = size(param,2)
      model%stationary = .true.
      if (present(stationary)) model%stationary = stationary
      if (present(optimizer)) model%mle_optimizer = optimizer
      allocate(model%param(size(param,1),size(param,2)))
      allocate(model%gamma(size(gamma_matrix,1),size(gamma_matrix,2)))
      model%param = param
      model%gamma = gamma_matrix

      local_status = ldhmm_validate(model, require_delta=.false.)
      if (local_status /= LDHMM_SUCCESS) then
         if (present(status)) status = local_status
         return
      end if

      allocate(model%delta(m))
      if (present(delta)) then
         if (size(delta) /= m) then
            model%delta = 0.0_dp
            if (present(status)) status = LDHMM_INVALID_ARGUMENT
            return
         end if
         model%delta = max(delta, 0.0_dp)
         if (sum(model%delta) <= 0.0_dp) then
            if (present(status)) status = LDHMM_INVALID_ARGUMENT
            return
         end if
         model%delta = model%delta / sum(model%delta)
      else if (model%stationary) then
         call ldhmm_stationary_distribution(model%gamma, model%delta, local_status)
         if (local_status /= LDHMM_SUCCESS) then
            if (present(status)) status = local_status
            return
         end if
      else
         model%delta = 1.0_dp / real(m, dp)
      end if
      if (present(status)) status = LDHMM_SUCCESS
   end function ldhmm_create

   integer function ldhmm_validate(model, require_delta) result(status)
      type(ldhmm_model), intent(in) :: model
      logical, intent(in), optional :: require_delta
      logical :: need_delta
      integer :: i
      real(dp), parameter :: tolerance = 1.0e-8_dp

      need_delta = .true.
      if (present(require_delta)) need_delta = require_delta
      status = LDHMM_INVALID_ARGUMENT
      if (model%m < 1) return
      if (.not. allocated(model%param) .or. .not. allocated(model%gamma)) return
      if (size(model%param,1) /= model%m) return
      if (size(model%param,2) /= 2 .and. size(model%param,2) /= 3) return
      if (size(model%gamma,1) /= model%m .or. size(model%gamma,2) /= model%m) return
      if (any(.not. ieee_is_finite(model%param))) return
      if (any(.not. ieee_is_finite(model%gamma))) return
      if (any(model%param(:,2) <= 0.0_dp)) return
      if (model%param_nbr == 3) then
         if (any(model%param(:,3) <= 0.0_dp)) return
      end if
      if (any(model%gamma < 0.0_dp)) return
      do i = 1, model%m
         if (abs(sum(model%gamma(i,:))-1.0_dp) > tolerance) return
      end do
      if (need_delta) then
         if (.not. allocated(model%delta)) return
         if (size(model%delta) /= model%m) return
         if (any(model%delta < 0.0_dp)) return
         if (abs(sum(model%delta)-1.0_dp) > tolerance) return
      end if
      status = LDHMM_SUCCESS
   end function ldhmm_validate

   function ldhmm_gamma_init(m, p1, p2, probability, min_gamma) result(gamma_matrix)
      integer, intent(in) :: m
      real(dp), intent(in), optional :: p1, p2
      real(dp), intent(in), optional :: probability(:, :)
      real(dp), intent(in), optional :: min_gamma
      real(dp), allocatable :: gamma_matrix(:, :)
      real(dp) :: first, second, floor_value
      integer :: i, j, distance

      first = 0.04_dp
      second = 0.01_dp
      floor_value = 0.0_dp
      if (present(p1)) first = p1
      if (present(p2)) second = p2
      if (present(min_gamma)) floor_value = max(0.0_dp, min_gamma)
      allocate(gamma_matrix(m,m))

      if (present(probability)) then
         if (size(probability,1) /= m .or. size(probability,2) /= m) then
            gamma_matrix = 0.0_dp
            return
         end if
         gamma_matrix = probability
         where (abs(gamma_matrix) <= tiny(1.0_dp))
            gamma_matrix = gamma_matrix + floor_value
         end where
         do i = 1, m
            if (sum(gamma_matrix(i,:)) <= 0.0_dp) then
               gamma_matrix(i,:) = 1.0_dp / real(m, dp)
            else
               gamma_matrix(i,:) = gamma_matrix(i,:) / sum(gamma_matrix(i,:))
            end if
         end do
         return
      end if

      if (m == 1) then
         gamma_matrix(1,1) = 1.0_dp
         return
      end if
      if (m == 2) then
         gamma_matrix = reshape([1.0_dp-first, first, first, 1.0_dp-first], [m,m], order=[2,1])
         return
      end if

      gamma_matrix = 0.0_dp
      do i = 1, m
         do j = 1, m
            if (i == j) cycle
            distance = abs(i-j)
            if (distance == 1) gamma_matrix(i,j) = first
            if (distance == 2) gamma_matrix(i,j) = second
         end do
         gamma_matrix(i,i) = 1.0_dp - sum(gamma_matrix(i,:))
      end do
   end function ldhmm_gamma_init

   subroutine ldhmm_stationary_distribution(gamma_matrix, delta, status)
      real(dp), intent(in) :: gamma_matrix(:, :)
      real(dp), allocatable, intent(out) :: delta(:)
      integer, intent(out) :: status
      real(dp), allocatable :: a(:, :), b(:), solution(:)
      integer :: n, solve_status

      n = size(gamma_matrix,1)
      if (n < 1 .or. size(gamma_matrix,2) /= n) then
         allocate(delta(0))
         status = LDHMM_INVALID_ARGUMENT
         return
      end if
      allocate(a(n,n), b(n))
      a = transpose(gamma_matrix) - identity_matrix(n)
      b = 0.0_dp
      a(n,:) = 1.0_dp
      b(n) = 1.0_dp
      call solve_linear_system(a, b, solution, solve_status)
      if (solve_status /= 0 .or. any(.not. ieee_is_finite(solution))) then
         allocate(delta(n))
         delta = 1.0_dp / real(n, dp)
         status = LDHMM_NUMERICAL_ERROR
         return
      end if
      allocate(delta(n))
      delta = max(solution, 0.0_dp)
      if (sum(delta) <= 0.0_dp) then
         delta = 1.0_dp / real(n, dp)
         status = LDHMM_NUMERICAL_ERROR
      else
         delta = delta / sum(delta)
         status = LDHMM_SUCCESS
      end if
   end subroutine ldhmm_stationary_distribution

   integer function ldhmm_parameter_count(model) result(n)
      type(ldhmm_model), intent(in) :: model
      n = model%m*model%param_nbr
      if (model%m > 1) n = n + model%m*(model%m-1)
      if (.not. model%stationary .and. model%m > 1) n = n + model%m-1
   end function ldhmm_parameter_count

   function ldhmm_natural_to_working(model, mu_scale) result(parameters)
      type(ldhmm_model), intent(in) :: model
      real(dp), intent(in), optional :: mu_scale
      real(dp), allocatable :: parameters(:)
      real(dp) :: scale
      integer :: i, j, k, n

      scale = 1.0_dp
      if (present(mu_scale)) scale = max(abs(mu_scale), sqrt(tiny(1.0_dp)))
      n = ldhmm_parameter_count(model)
      allocate(parameters(n))
      k = 0
      do i = 1, model%m
         do j = 1, model%param_nbr
            k = k + 1
            if (j == 1) then
               parameters(k) = model%param(i,j) / scale
            else
               parameters(k) = log(max(abs(model%param(i,j)), tiny(1.0_dp)))
            end if
         end do
      end do
      if (model%m == 1) return
      do j = 1, model%m
         do i = 1, model%m
            if (i == j) cycle
            k = k + 1
            parameters(k) = log(max(model%gamma(i,j), tiny(1.0_dp)) / &
               max(model%gamma(i,i), tiny(1.0_dp)))
         end do
      end do
      if (.not. model%stationary) then
         do i = 2, model%m
            k = k + 1
            parameters(k) = log(max(model%delta(i), tiny(1.0_dp)) / &
               max(model%delta(1), tiny(1.0_dp)))
         end do
      end if
   end function ldhmm_natural_to_working

   function ldhmm_working_to_natural(template, parameters, mu_scale, status) result(model)
      type(ldhmm_model), intent(in) :: template
      real(dp), intent(in) :: parameters(:)
      real(dp), intent(in), optional :: mu_scale
      integer, intent(out), optional :: status
      type(ldhmm_model) :: model
      real(dp) :: scale
      real(dp), allocatable :: delta_values(:)
      integer :: i, j, k, local_status

      model = template
      scale = 1.0_dp
      if (present(mu_scale)) scale = max(abs(mu_scale), sqrt(tiny(1.0_dp)))
      if (size(parameters) /= ldhmm_parameter_count(template)) then
         if (present(status)) status = LDHMM_INVALID_ARGUMENT
         return
      end if
      if (.not. allocated(model%param)) allocate(model%param(template%m,template%param_nbr))
      if (.not. allocated(model%gamma)) allocate(model%gamma(template%m,template%m))
      if (.not. allocated(model%delta)) allocate(model%delta(template%m))

      k = 0
      do i = 1, template%m
         do j = 1, template%param_nbr
            k = k + 1
            if (j == 1) then
               model%param(i,j) = parameters(k) * scale
            else
               model%param(i,j) = exp(min(parameters(k), log(1.0e100_dp)))
            end if
         end do
      end do
      model%gamma = identity_matrix(template%m)
      if (template%m == 1) then
         model%delta = 1.0_dp
         if (present(status)) status = LDHMM_SUCCESS
         return
      end if
      do j = 1, template%m
         do i = 1, template%m
            if (i == j) cycle
            k = k + 1
            model%gamma(i,j) = exp(min(parameters(k), log(1.0e100_dp)))
         end do
      end do
      do i = 1, template%m
         model%gamma(i,:) = model%gamma(i,:) / sum(model%gamma(i,:))
      end do

      if (template%stationary) then
         call ldhmm_stationary_distribution(model%gamma, delta_values, local_status)
         model%delta = delta_values
      else
         allocate(delta_values(template%m))
         delta_values(1) = 1.0_dp
         do i = 2, template%m
            k = k + 1
            delta_values(i) = exp(min(parameters(k), log(1.0e100_dp)))
         end do
         model%delta = delta_values / sum(delta_values)
         local_status = LDHMM_SUCCESS
      end if
      if (present(status)) status = local_status
   end function ldhmm_working_to_natural

   pure function identity_matrix(n) result(identity)
      integer, intent(in) :: n
      real(dp) :: identity(n,n)
      integer :: i
      identity = 0.0_dp
      do i = 1, n
         identity(i,i) = 1.0_dp
      end do
   end function identity_matrix

end module ldhmm_parameters
