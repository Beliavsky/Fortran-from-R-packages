! SPDX-License-Identifier: Artistic-2.0
module ldhmm_optimization
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use ldhmm_kinds, only : dp
   use ldhmm_math, only : finite_count, finite_mean_abs
   use ldhmm_modeling, only : ldhmm_decode, ldhmm_mllk
   use ldhmm_parameters, only : ldhmm_gamma_init, ldhmm_natural_to_working, &
      ldhmm_parameter_count, ldhmm_stationary_distribution, ldhmm_working_to_natural
   use ldhmm_status, only : LDHMM_SUCCESS, LDHMM_INVALID_ARGUMENT, &
      LDHMM_MAX_ITERATIONS, LDHMM_LINE_SEARCH_FAILED, LDHMM_NUMERICAL_ERROR
   use ldhmm_types, only : ldhmm_fit_control, ldhmm_model
   implicit none
   private

   public :: ldhmm_fit

contains

   function ldhmm_fit(input_model, x, control, status) result(model)
      type(ldhmm_model), intent(in) :: input_model
      real(dp), intent(in) :: x(:)
      type(ldhmm_fit_control), intent(in), optional :: control
      integer, intent(out), optional :: status
      type(ldhmm_model) :: model, base_model
      type(ldhmm_fit_control) :: settings
      real(dp), allocatable :: initial(:), solution(:), normalized_gamma(:, :)
      real(dp), allocatable :: stationary_delta(:)
      real(dp) :: objective, mu_scale
      integer :: local_status, iterations, nobs, npar
      character(len=24) :: optimizer

      settings = ldhmm_fit_control()
      settings%optimizer = input_model%mle_optimizer
      if (present(control)) settings = control
      optimizer = lower_string(trim(settings%optimizer))
      if (size(x) < 1) then
         model = input_model
         if (present(status)) status = LDHMM_INVALID_ARGUMENT
         return
      end if

      base_model = input_model
      normalized_gamma = ldhmm_gamma_init(base_model%m, probability=base_model%gamma, &
         min_gamma=settings%min_gamma)
      base_model%gamma = normalized_gamma
      if (base_model%stationary) then
         call ldhmm_stationary_distribution(base_model%gamma, stationary_delta, local_status)
         if (local_status == LDHMM_SUCCESS) base_model%delta = stationary_delta
      end if
      mu_scale = finite_mean_abs(x)
      if (.not. ieee_is_finite(mu_scale) .or. mu_scale <= sqrt(tiny(1.0_dp))) mu_scale = 1.0_dp
      initial = ldhmm_natural_to_working(base_model, mu_scale)

      select case (optimizer)
      case ('nelder-mead', 'nelder_mead', 'nm', 'nlm')
         call nelder_mead_minimize(base_model, x, mu_scale, initial, settings, &
            solution, objective, iterations, local_status)
      case ('bfgs', 'rvmmin', 'optimx')
         call bfgs_minimize(base_model, x, mu_scale, initial, settings, &
            solution, objective, iterations, local_status)
      case default
         call bfgs_minimize(base_model, x, mu_scale, initial, settings, &
            solution, objective, iterations, local_status)
      end select

      model = ldhmm_working_to_natural(base_model, solution, mu_scale)
      model%mllk = objective
      model%iterations = iterations
      model%return_code = local_status
      model%fitted = local_status == LDHMM_SUCCESS .or. local_status == LDHMM_MAX_ITERATIONS
      npar = ldhmm_parameter_count(model)
      nobs = finite_count(x)
      model%aic = 2.0_dp*(model%mllk+real(npar,dp))
      if (nobs > 0) then
         model%bic = 2.0_dp*model%mllk + real(npar,dp)*log(real(nobs,dp))
      else
         model%bic = huge(1.0_dp)
      end if
      if (settings%decode) model = ldhmm_decode(model, x)
      if (present(status)) status = local_status
   end function ldhmm_fit

   real(dp) function evaluate_objective(base_model, observations, mu_scale, &
      parameters) result(value)
      type(ldhmm_model), intent(in) :: base_model
      real(dp), intent(in) :: observations(:), mu_scale, parameters(:)
      type(ldhmm_model) :: candidate
      integer :: objective_status

      candidate = ldhmm_working_to_natural(base_model, parameters, mu_scale, &
         objective_status)
      if (objective_status /= LDHMM_SUCCESS) then
         value = huge(1.0_dp)
         return
      end if
      value = ldhmm_mllk(candidate, observations, objective_status)
      if (.not. ieee_is_finite(value)) value = huge(1.0_dp)
   end function evaluate_objective

   subroutine bfgs_minimize(base_model, observations, mu_scale, x0, control, &
      solution, value, iterations, status)
      type(ldhmm_model), intent(in) :: base_model
      real(dp), intent(in) :: observations(:), mu_scale, x0(:)
      type(ldhmm_fit_control), intent(in) :: control
      real(dp), allocatable, intent(out) :: solution(:)
      real(dp), intent(out) :: value
      integer, intent(out) :: iterations, status
      real(dp), allocatable :: x(:), x_new(:), gradient(:), gradient_new(:)
      real(dp), allocatable :: direction(:), hessian_inverse(:, :), identity(:, :)
      real(dp), allocatable :: s(:), y(:), hy(:), outer_ss(:, :), update(:, :)
      real(dp) :: alpha, directional_derivative, candidate_value, ys, yhy
      real(dp) :: gradient_norm, step_norm
      integer :: i, line_iteration, n

      n = size(x0)
      allocate(x(n), x_new(n), gradient(n), gradient_new(n), direction(n))
      allocate(hessian_inverse(n,n), identity(n,n), s(n), y(n), hy(n))
      allocate(outer_ss(n,n), update(n,n), solution(n))
      identity = 0.0_dp
      do i = 1, n
         identity(i,i) = 1.0_dp
      end do
      hessian_inverse = identity
      x = x0
      value = evaluate_objective(base_model, observations, mu_scale, x)
      call numerical_gradient(base_model, observations, mu_scale, x, &
         control%gradient_step, gradient)
      status = LDHMM_MAX_ITERATIONS
      iterations = 0

      do iterations = 1, control%max_iterations
         gradient_norm = maxval(abs(gradient))
         if (gradient_norm <= control%tolerance*(1.0_dp+abs(value))) then
            status = LDHMM_SUCCESS
            exit
         end if
         direction = -matmul(hessian_inverse, gradient)
         directional_derivative = dot_product(gradient, direction)
         if (.not. ieee_is_finite(directional_derivative) .or. &
             directional_derivative >= 0.0_dp) then
            direction = -gradient
            directional_derivative = -dot_product(gradient, gradient)
            hessian_inverse = identity
         end if

         alpha = 1.0_dp
         candidate_value = huge(1.0_dp)
         do line_iteration = 1, 40
            x_new = x + alpha*direction
            candidate_value = evaluate_objective(base_model, observations, &
               mu_scale, x_new)
            if (ieee_is_finite(candidate_value)) then
               if (candidate_value <= value + &
                   1.0e-4_dp*alpha*directional_derivative) exit
            end if
            alpha = 0.5_dp*alpha
         end do
         if (line_iteration > 40) then
            status = LDHMM_LINE_SEARCH_FAILED
            exit
         end if

         call numerical_gradient(base_model, observations, mu_scale, x_new, &
            control%gradient_step, gradient_new)
         s = x_new - x
         y = gradient_new - gradient
         ys = dot_product(y, s)
         step_norm = maxval(abs(s))
         if (ys > sqrt(epsilon(1.0_dp))*sqrt(dot_product(y,y)*dot_product(s,s))) then
            hy = matmul(hessian_inverse, y)
            yhy = dot_product(y, hy)
            call outer_product(s, s, outer_ss)
            call outer_product(hy, s, update)
            hessian_inverse = hessian_inverse + &
               ((ys+yhy)/(ys*ys))*outer_ss - (update+transpose(update))/ys
         else
            hessian_inverse = identity
         end if
         x = x_new
         gradient = gradient_new
         value = candidate_value
         if (control%print_level >= 2) then
            write(*,'(a,i0,a,es16.8,a,es12.4)') 'iteration ', iterations, &
               ' objective=', value, ' gradient=', maxval(abs(gradient))
         end if
         if (step_norm <= control%tolerance*(1.0_dp+maxval(abs(x)))) then
            status = LDHMM_SUCCESS
            exit
         end if
      end do
      if (iterations > control%max_iterations) iterations = control%max_iterations
      solution = x
   end subroutine bfgs_minimize

   subroutine nelder_mead_minimize(base_model, observations, mu_scale, x0, &
      control, solution, value, iterations, status)
      type(ldhmm_model), intent(in) :: base_model
      real(dp), intent(in) :: observations(:), mu_scale, x0(:)
      type(ldhmm_fit_control), intent(in) :: control
      real(dp), allocatable, intent(out) :: solution(:)
      real(dp), intent(out) :: value
      integer, intent(out) :: iterations, status
      real(dp), allocatable :: simplex(:, :), values(:), centroid(:), reflected(:)
      real(dp), allocatable :: expanded(:), contracted(:)
      real(dp) :: reflected_value, expanded_value, contracted_value
      real(dp) :: spread_value, spread_x, step
      integer :: j, n

      n = size(x0)
      allocate(simplex(n,n+1), values(n+1), centroid(n), reflected(n))
      allocate(expanded(n), contracted(n), solution(n))
      simplex(:,1) = x0
      do j = 1, n
         simplex(:,j+1) = x0
         step = control%initial_simplex_step*max(1.0_dp,abs(x0(j)))
         simplex(j,j+1) = simplex(j,j+1) + step
      end do
      do j = 1, n+1
         values(j) = evaluate_objective(base_model, observations, mu_scale, &
            simplex(:,j))
      end do
      status = LDHMM_MAX_ITERATIONS
      iterations = 0

      do iterations = 1, control%max_iterations
         call sort_simplex(simplex, values)
         spread_value = maxval(abs(values-values(1)))
         spread_x = 0.0_dp
         do j = 2, n+1
            spread_x = max(spread_x, maxval(abs(simplex(:,j)-simplex(:,1))))
         end do
         if (spread_value <= control%tolerance*(1.0_dp+abs(values(1))) .and. &
             spread_x <= sqrt(control%tolerance)*(1.0_dp+maxval(abs(simplex(:,1))))) then
            status = LDHMM_SUCCESS
            exit
         end if

         centroid = 0.0_dp
         do j = 1, n
            centroid = centroid + simplex(:,j)
         end do
         centroid = centroid/real(n,dp)
         reflected = centroid + (centroid-simplex(:,n+1))
         reflected_value = evaluate_objective(base_model, observations, mu_scale, &
            reflected)

         if (reflected_value < values(1)) then
            expanded = centroid + 2.0_dp*(reflected-centroid)
            expanded_value = evaluate_objective(base_model, observations, &
               mu_scale, expanded)
            if (expanded_value < reflected_value) then
               simplex(:,n+1) = expanded
               values(n+1) = expanded_value
            else
               simplex(:,n+1) = reflected
               values(n+1) = reflected_value
            end if
         else if (reflected_value < values(n)) then
            simplex(:,n+1) = reflected
            values(n+1) = reflected_value
         else
            if (reflected_value < values(n+1)) then
               contracted = centroid + 0.5_dp*(reflected-centroid)
            else
               contracted = centroid + 0.5_dp*(simplex(:,n+1)-centroid)
            end if
            contracted_value = evaluate_objective(base_model, observations, &
               mu_scale, contracted)
            if (contracted_value < min(reflected_value,values(n+1))) then
               simplex(:,n+1) = contracted
               values(n+1) = contracted_value
            else
               do j = 2, n+1
                  simplex(:,j) = simplex(:,1) + 0.5_dp*(simplex(:,j)-simplex(:,1))
                  values(j) = evaluate_objective(base_model, observations, &
                     mu_scale, simplex(:,j))
               end do
            end if
         end if
         if (control%print_level >= 2) then
            write(*,'(a,i0,a,es16.8)') 'iteration ', iterations, &
               ' objective=', minval(values)
         end if
      end do
      call sort_simplex(simplex, values)
      solution = simplex(:,1)
      value = values(1)
      if (.not. ieee_is_finite(value)) status = LDHMM_NUMERICAL_ERROR
      if (iterations > control%max_iterations) iterations = control%max_iterations
   end subroutine nelder_mead_minimize

   subroutine numerical_gradient(base_model, observations, mu_scale, x, &
      relative_step, gradient)
      type(ldhmm_model), intent(in) :: base_model
      real(dp), intent(in) :: observations(:), mu_scale, x(:), relative_step
      real(dp), intent(out) :: gradient(:)
      real(dp), allocatable :: plus(:), minus(:)
      real(dp) :: step, fplus, fminus
      integer :: i

      allocate(plus(size(x)), minus(size(x)))
      do i = 1, size(x)
         plus = x
         minus = x
         step = max(relative_step, sqrt(epsilon(1.0_dp)))*max(1.0_dp,abs(x(i)))
         plus(i) = plus(i) + step
         minus(i) = minus(i) - step
         fplus = evaluate_objective(base_model, observations, mu_scale, plus)
         fminus = evaluate_objective(base_model, observations, mu_scale, minus)
         if (ieee_is_finite(fplus) .and. ieee_is_finite(fminus)) then
            gradient(i) = (fplus-fminus)/(2.0_dp*step)
         else
            gradient(i) = 0.0_dp
         end if
      end do
   end subroutine numerical_gradient

   subroutine outer_product(x, y, matrix)
      real(dp), intent(in) :: x(:), y(:)
      real(dp), intent(out) :: matrix(:, :)
      integer :: i, j
      do j = 1, size(y)
         do i = 1, size(x)
            matrix(i,j) = x(i)*y(j)
         end do
      end do
   end subroutine outer_product

   subroutine sort_simplex(simplex, values)
      real(dp), intent(inout) :: simplex(:, :), values(:)
      real(dp), allocatable :: column(:)
      real(dp) :: value_temp
      integer :: i, j, best

      allocate(column(size(simplex,1)))
      do i = 1, size(values)-1
         best = i
         do j = i+1, size(values)
            if (values(j) < values(best)) best = j
         end do
         if (best /= i) then
            value_temp = values(i)
            values(i) = values(best)
            values(best) = value_temp
            column = simplex(:,i)
            simplex(:,i) = simplex(:,best)
            simplex(:,best) = column
         end if
      end do
   end subroutine sort_simplex

   pure function lower_string(text) result(lower)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: lower
      integer :: i, code
      do i = 1, len(text)
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            lower(i:i) = achar(code+32)
         else
            lower(i:i) = text(i:i)
         end if
      end do
   end function lower_string

end module ldhmm_optimization
