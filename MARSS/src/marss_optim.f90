! SPDX-License-Identifier: GPL-2.0-only
module marss_optim_mod
   use marss_kinds, only : dp
   use marss_types, only : marss_model, marss_fit_result, marss_kf_result
   use marss_kalman, only : marss_kf
   use marss_utils, only : identity_matrix
   use r_linalg, only : cholesky_factor
   implicit none
   private
   public :: marss_optim

   type :: optim_flags
      logical :: b = .true.
      logical :: u = .true.
      logical :: q = .true.
      logical :: z = .true.
      logical :: a = .true.
      logical :: r = .true.
      logical :: x0 = .true.
      logical :: v0 = .false.
   end type optim_flags

contains

   subroutine marss_optim(start_model, max_iter, tol, fit, estimate_b, estimate_u, estimate_q, &
      estimate_z, estimate_a, estimate_r, estimate_x0, estimate_v0)
      type(marss_model), intent(in) :: start_model !! Starting time-invariant MARSS model, y may contain NaNs.
      integer, intent(in) :: max_iter !! Maximum BFGS iterations, which must be positive.
      real(dp), intent(in) :: tol !! Relative objective and gradient convergence tolerance, which must be positive.
      type(marss_fit_result), intent(out) :: fit !! Optimized model, Kalman result, status, iteration count, and log likelihood.
      logical, intent(in), optional :: estimate_b !! Optimize the state-transition matrix B, defaulting to true.
      logical, intent(in), optional :: estimate_u !! Optimize the state intercept U, defaulting to true.
      logical, intent(in), optional :: estimate_q !! Optimize Q with lower-Cholesky coordinates, defaults to true.
      logical, intent(in), optional :: estimate_z !! Optimize observation loading matrix Z, defaulting to true.
      logical, intent(in), optional :: estimate_a !! Optimize observation intercept A, defaulting to true.
      logical, intent(in), optional :: estimate_r !! Optimize R with lower-Cholesky coordinates, defaults to true.
      logical, intent(in), optional :: estimate_x0 !! Optimize initial-state mean x0, defaulting to true.
      logical, intent(in), optional :: estimate_v0 !! Optimize V0 with lower-Cholesky coordinates, defaults to false.
      type(optim_flags) :: flags
      type(marss_model) :: final_model
      real(dp), allocatable :: theta(:)
      real(dp), allocatable :: theta_new(:)
      real(dp), allocatable :: gradient(:)
      real(dp), allocatable :: gradient_new(:)
      real(dp), allocatable :: h_inv(:, :)
      real(dp), allocatable :: direction(:)
      real(dp), allocatable :: s(:)
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: left(:, :)
      real(dp), allocatable :: right(:, :)
      real(dp) :: alpha
      real(dp) :: f
      real(dp) :: f_new
      real(dp) :: gtd
      real(dp) :: rho
      real(dp) :: ys
      integer :: info
      integer :: iter
      integer :: p

      fit%converged = .false.
      fit%iterations = 0
      fit%info = 0
      if (max_iter < 1 .or. tol <= 0.0_dp) then
         fit%info = 1
         return
      end if
      call set_flags(flags, estimate_b, estimate_u, estimate_q, estimate_z, estimate_a, estimate_r, &
         estimate_x0, estimate_v0)
      if ((flags%b .and. allocated(start_model%b_t)) .or. &
         (flags%u .and. allocated(start_model%u_t)) .or. &
         (flags%q .and. allocated(start_model%q_t)) .or. &
         (flags%z .and. allocated(start_model%z_t)) .or. &
         (flags%a .and. allocated(start_model%a_t)) .or. &
         (flags%r .and. allocated(start_model%r_t))) then
         fit%info = 2
         return
      end if
      if ((flags%q .and. allocated(start_model%q_noise)) .or. &
         (flags%r .and. allocated(start_model%r_noise)) .or. &
         (flags%v0 .and. allocated(start_model%v0_noise))) then
         fit%info = 3
         return
      end if
      call pack_optim(start_model, flags, theta, info)
      if (info /= 0) then
         fit%info = 10 + info
         return
      end if
      p = size(theta)
      if (p == 0) then
         fit%model = start_model
         call marss_kf(fit%model, fit%kf)
         if (.not. fit%kf%ok) then
            fit%info = 20 + fit%kf%info
            return
         end if
         fit%loglik = fit%kf%loglik
         fit%converged = .true.
         return
      end if
      allocate(theta_new(p), gradient(p), gradient_new(p), h_inv(p, p), direction(p), s(p), y(p))
      allocate(left(p, p), right(p, p))
      call objective(start_model, flags, theta, f, info)
      if (info /= 0) then
         fit%info = 30 + info
         return
      end if
      call objective_gradient(start_model, flags, theta, f, gradient, info)
      if (info /= 0) then
         fit%info = 40 + info
         return
      end if
      h_inv = identity_matrix(p)

      do iter = 1, max_iter
         if (maxval(abs(gradient)) <= tol) then
            fit%converged = .true.
            exit
         end if
         direction = -matmul(h_inv, gradient)
         gtd = dot_product(gradient, direction)
         if (gtd >= -epsilon(1.0_dp)) then
            h_inv = identity_matrix(p)
            direction = -gradient
            gtd = -dot_product(gradient, gradient)
         end if
         call armijo_step(start_model, flags, theta, f, gradient, direction, theta_new, f_new, alpha, info)
         if (info /= 0) then
            fit%info = 50 + info
            return
         end if
         call objective_gradient(start_model, flags, theta_new, f_new, gradient_new, info)
         if (info /= 0) then
            fit%info = 60 + info
            return
         end if
         s = theta_new - theta
         y = gradient_new - gradient
         ys = dot_product(y, s)
         if (ys > sqrt(epsilon(1.0_dp)) * max(1.0_dp, sqrt(dot_product(y, y) * dot_product(s, s)))) then
            rho = 1.0_dp / ys
            left = identity_matrix(p) - rho * outer_product(s, y)
            right = identity_matrix(p) - rho * outer_product(y, s)
            h_inv = matmul(matmul(left, h_inv), right) + rho * outer_product(s, s)
         else
            h_inv = identity_matrix(p)
         end if
         if (abs(f_new - f) <= tol * (1.0_dp + abs(f))) fit%converged = .true.
         theta = theta_new
         gradient = gradient_new
         f = f_new
         if (fit%converged) exit
      end do
      fit%iterations = min(iter, max_iter)
      call unpack_optim(start_model, flags, theta, final_model, info)
      if (info /= 0) then
         fit%info = 70 + info
         return
      end if
      fit%model = final_model
      call marss_kf(fit%model, fit%kf)
      if (.not. fit%kf%ok) then
         fit%info = 80 + fit%kf%info
         return
      end if
      fit%loglik = fit%kf%loglik
      fit%info = 0
   end subroutine marss_optim

   pure subroutine set_flags(flags, estimate_b, estimate_u, estimate_q, estimate_z, estimate_a, &
      estimate_r, estimate_x0, estimate_v0)
      type(optim_flags), intent(out) :: flags !! Resolved logical flags controlling inclusion of each model block.
      logical, intent(in), optional :: estimate_b !! Optional override for optimizing B.
      logical, intent(in), optional :: estimate_u !! Optional override for optimizing U.
      logical, intent(in), optional :: estimate_q !! Optional override for optimizing Q.
      logical, intent(in), optional :: estimate_z !! Optional override for optimizing Z.
      logical, intent(in), optional :: estimate_a !! Optional override for optimizing A.
      logical, intent(in), optional :: estimate_r !! Optional override for optimizing R.
      logical, intent(in), optional :: estimate_x0 !! Optional override for optimizing x0.
      logical, intent(in), optional :: estimate_v0 !! Optional override for optimizing V0.

      flags = optim_flags()
      if (present(estimate_b)) flags%b = estimate_b
      if (present(estimate_u)) flags%u = estimate_u
      if (present(estimate_q)) flags%q = estimate_q
      if (present(estimate_z)) flags%z = estimate_z
      if (present(estimate_a)) flags%a = estimate_a
      if (present(estimate_r)) flags%r = estimate_r
      if (present(estimate_x0)) flags%x0 = estimate_x0
      if (present(estimate_v0)) flags%v0 = estimate_v0
   end subroutine set_flags

   pure subroutine pack_optim(model, flags, theta, info)
      type(marss_model), intent(in) :: model !! Model whose selected blocks are transformed into the optimizer vector.
      type(optim_flags), intent(in) :: flags !! Block-selection flags for the optimizer vector.
      real(dp), allocatable, intent(out) :: theta(:) !! Packed optimizer coordinates with log-Cholesky covariance diagonals.
      integer, intent(out) :: info !! Zero on success, or a positive code if a selected covariance is not positive definite.
      real(dp), allocatable :: factor(:, :)
      integer :: k
      integer :: m
      integer :: n
      integer :: p

      m = size(model%b, 1)
      n = size(model%z, 1)
      p = 0
      if (flags%b) p = p + m * m
      if (flags%u) p = p + m
      if (flags%q) p = p + m * (m + 1) / 2
      if (flags%z) p = p + n * m
      if (flags%a) p = p + n
      if (flags%r) p = p + n * (n + 1) / 2
      if (flags%x0) p = p + m
      if (flags%v0) p = p + m * (m + 1) / 2
      allocate(theta(p))
      k = 0
      if (flags%b) call pack_full(model%b, theta, k)
      if (flags%u) call pack_vector(model%u, theta, k)
      if (flags%q) then
         call cholesky_factor(model%q, factor, info)
         if (info /= 0) return
         call pack_cholesky(factor, theta, k)
         deallocate(factor)
      end if
      if (flags%z) call pack_full(model%z, theta, k)
      if (flags%a) call pack_vector(model%a, theta, k)
      if (flags%r) then
         call cholesky_factor(model%r, factor, info)
         if (info /= 0) return
         call pack_cholesky(factor, theta, k)
         deallocate(factor)
      end if
      if (flags%x0) call pack_vector(model%x0, theta, k)
      if (flags%v0) then
         call cholesky_factor(model%v0, factor, info)
         if (info /= 0) return
         call pack_cholesky(factor, theta, k)
      end if
      info = 0
   end subroutine pack_optim

   pure subroutine unpack_optim(template, flags, theta, model, info)
      type(marss_model), intent(in) :: template !! Template supplying observations, tinitx, and all unoptimized blocks.
      type(optim_flags), intent(in) :: flags !! Block-selection flags defining the optimizer vector layout.
      real(dp), intent(in) :: theta(:) !! Optimizer coordinates produced by pack_optim or subsequent BFGS steps.
      type(marss_model), intent(out) :: model !! Reconstructed model with positive-definite optimized covariance blocks.
      integer, intent(out) :: info !! Zero on success, or nonzero when theta has the wrong length.
      real(dp), allocatable :: factor(:, :)
      integer :: expected
      integer :: k
      integer :: m
      integer :: n

      model = template
      m = size(model%b, 1)
      n = size(model%z, 1)
      expected = optim_count(model, flags)
      if (size(theta) /= expected) then
         info = 1
         return
      end if
      k = 0
      if (flags%b) call unpack_full(theta, k, model%b)
      if (flags%u) call unpack_vector(theta, k, model%u)
      if (flags%q) then
         call unpack_cholesky(theta, k, m, factor)
         model%q = matmul(factor, transpose(factor))
         deallocate(factor)
      end if
      if (flags%z) call unpack_full(theta, k, model%z)
      if (flags%a) call unpack_vector(theta, k, model%a)
      if (flags%r) then
         call unpack_cholesky(theta, k, n, factor)
         model%r = matmul(factor, transpose(factor))
         deallocate(factor)
      end if
      if (flags%x0) call unpack_vector(theta, k, model%x0)
      if (flags%v0) then
         call unpack_cholesky(theta, k, m, factor)
         model%v0 = matmul(factor, transpose(factor))
      end if
      info = 0
   end subroutine unpack_optim

   pure integer function optim_count(model, flags) result(p)
      type(marss_model), intent(in) :: model !! Model supplying state and observation dimensions.
      type(optim_flags), intent(in) :: flags !! Block-selection flags whose active coordinate count is requested.
      integer :: m
      integer :: n

      m = size(model%b, 1)
      n = size(model%z, 1)
      p = 0
      if (flags%b) p = p + m * m
      if (flags%u) p = p + m
      if (flags%q) p = p + m * (m + 1) / 2
      if (flags%z) p = p + n * m
      if (flags%a) p = p + n
      if (flags%r) p = p + n * (n + 1) / 2
      if (flags%x0) p = p + m
      if (flags%v0) p = p + m * (m + 1) / 2
   end function optim_count

   subroutine objective(template, flags, theta, value, info)
      type(marss_model), intent(in) :: template !! Template model for reconstructing optimizer coordinates.
      type(optim_flags), intent(in) :: flags !! Block-selection flags defining theta.
      real(dp), intent(in) :: theta(:) !! Optimizer coordinates at which the negative log likelihood is evaluated.
      real(dp), intent(out) :: value !! Negative Kalman innovations log likelihood, or huge on failure.
      integer, intent(out) :: info !! Zero on success, or nonzero if reconstruction/filtering fails.
      type(marss_model) :: model
      type(marss_kf_result) :: kf

      call unpack_optim(template, flags, theta, model, info)
      if (info /= 0) then
         value = huge(1.0_dp)
         return
      end if
      call marss_kf(model, kf, smoother=.false.)
      if (.not. kf%ok) then
         info = 100 + kf%info
         value = huge(1.0_dp)
         return
      end if
      value = -kf%loglik
      info = 0
   end subroutine objective

   subroutine objective_gradient(template, flags, theta, f0, gradient, info)
      type(marss_model), intent(in) :: template !! Template model for finite-difference objective evaluations.
      type(optim_flags), intent(in) :: flags !! Block-selection flags defining theta.
      real(dp), intent(in) :: theta(:) !! Current optimizer coordinates.
      real(dp), intent(in) :: f0 !! Objective value at theta, used for one-sided fallback differences.
      real(dp), intent(out) :: gradient(:) !! Central finite-difference gradient of the negative log likelihood.
      integer, intent(out) :: info !! Zero on success, or the first failed objective-evaluation code.
      real(dp), allocatable :: work(:)
      real(dp) :: fm
      real(dp) :: fp
      real(dp) :: h
      integer :: im
      integer :: ip
      integer :: j

      if (size(gradient) /= size(theta)) then
         info = 1
         return
      end if
      allocate(work(size(theta)))
      do j = 1, size(theta)
         h = 1.0e-5_dp * (1.0_dp + abs(theta(j)))
         work = theta
         work(j) = theta(j) + h
         call objective(template, flags, work, fp, ip)
         work(j) = theta(j) - h
         call objective(template, flags, work, fm, im)
         if (ip == 0 .and. im == 0) then
            gradient(j) = (fp - fm) / (2.0_dp * h)
         else if (ip == 0) then
            gradient(j) = (fp - f0) / h
         else if (im == 0) then
            gradient(j) = (f0 - fm) / h
         else
            info = 10 + j
            return
         end if
      end do
      info = 0
   end subroutine objective_gradient

   subroutine armijo_step(template, flags, theta, f, gradient, direction, theta_new, f_new, alpha, info)
      type(marss_model), intent(in) :: template !! Template model for line-search objective evaluations.
      type(optim_flags), intent(in) :: flags !! Block-selection flags defining theta.
      real(dp), intent(in) :: theta(:) !! Current optimizer coordinates.
      real(dp), intent(in) :: f !! Current negative log likelihood.
      real(dp), intent(in) :: gradient(:) !! Current objective gradient.
      real(dp), intent(in) :: direction(:) !! Descent direction proposed by BFGS.
      real(dp), intent(out) :: theta_new(:) !! Accepted optimizer coordinates.
      real(dp), intent(out) :: f_new !! Objective value at theta_new.
      real(dp), intent(out) :: alpha !! Accepted line-search step length.
      integer, intent(out) :: info !! Zero on success, or nonzero if no acceptable step is found.
      real(dp), parameter :: c1 = 1.0e-4_dp
      real(dp) :: slope
      integer :: eval_info
      integer :: trial

      slope = dot_product(gradient, direction)
      alpha = 1.0_dp
      do trial = 1, 32
         theta_new = theta + alpha * direction
         call objective(template, flags, theta_new, f_new, eval_info)
         if (eval_info == 0) then
            if (f_new <= f + c1 * alpha * slope) then
               info = 0
               return
            end if
         end if
         alpha = 0.5_dp * alpha
      end do
      theta_new = theta
      f_new = f
      info = 1
   end subroutine armijo_step

   pure subroutine pack_full(a, theta, k)
      real(dp), intent(in) :: a(:, :) !! Full matrix appended in Fortran column-major order.
      real(dp), intent(inout) :: theta(:) !! Destination optimizer vector.
      integer, intent(inout) :: k !! Number of optimizer coordinates already filled on entry and after appending on return.
      integer :: i
      integer :: j

      do j = 1, size(a, 2)
         do i = 1, size(a, 1)
            k = k + 1
            theta(k) = a(i, j)
         end do
      end do
   end subroutine pack_full

   pure subroutine unpack_full(theta, k, a)
      real(dp), intent(in) :: theta(:) !! Source optimizer vector.
      integer, intent(inout) :: k !! Number of optimizer coordinates already consumed on entry and after reading on return.
      real(dp), intent(out) :: a(:, :) !! Reconstructed full matrix in column-major order.
      integer :: i
      integer :: j

      do j = 1, size(a, 2)
         do i = 1, size(a, 1)
            k = k + 1
            a(i, j) = theta(k)
         end do
      end do
   end subroutine unpack_full

   pure subroutine pack_vector(a, theta, k)
      real(dp), intent(in) :: a(:) !! Vector appended to the optimizer coordinates.
      real(dp), intent(inout) :: theta(:) !! Destination optimizer vector.
      integer, intent(inout) :: k !! Number of optimizer coordinates already filled on entry and after appending on return.

      theta(k + 1:k + size(a)) = a
      k = k + size(a)
   end subroutine pack_vector

   pure subroutine unpack_vector(theta, k, a)
      real(dp), intent(in) :: theta(:) !! Source optimizer vector.
      integer, intent(inout) :: k !! Number of optimizer coordinates already consumed on entry and after reading on return.
      real(dp), intent(out) :: a(:) !! Reconstructed vector.

      a = theta(k + 1:k + size(a))
      k = k + size(a)
   end subroutine unpack_vector

   pure subroutine pack_cholesky(factor, theta, k)
      real(dp), intent(in) :: factor(:, :) !! Lower-triangular Cholesky factor with strictly positive diagonal.
      real(dp), intent(inout) :: theta(:) !! Destination optimizer vector receiving log diagonals and raw lower entries.
      integer, intent(inout) :: k !! Number of optimizer coordinates already filled on entry and after appending on return.
      integer :: i
      integer :: j

      do j = 1, size(factor, 1)
         do i = j, size(factor, 1)
            k = k + 1
            if (i == j) then
               theta(k) = log(factor(i, j))
            else
               theta(k) = factor(i, j)
            end if
         end do
      end do
   end subroutine pack_cholesky

   pure subroutine unpack_cholesky(theta, k, n, factor)
      real(dp), intent(in) :: theta(:) !! Source optimizer vector containing log diagonals and raw lower entries.
      integer, intent(inout) :: k !! Number of optimizer coordinates already consumed on entry and after reading on return.
      integer, intent(in) :: n !! Order of the square lower-Cholesky factor to reconstruct.
      real(dp), allocatable, intent(out) :: factor(:, :) !! Reconstructed lower-triangular factor with positive diagonal.
      integer :: i
      integer :: j

      allocate(factor(n, n))
      factor = 0.0_dp
      do j = 1, n
         do i = j, n
            k = k + 1
            if (i == j) then
               factor(i, j) = exp(theta(k))
            else
               factor(i, j) = theta(k)
            end if
         end do
      end do
   end subroutine unpack_cholesky

   pure function outer_product(x, y) result(a)
      real(dp), intent(in) :: x(:) !! Left vector in the outer product.
      real(dp), intent(in) :: y(:) !! Right vector in the outer product.
      real(dp) :: a(size(x), size(y))

      a = spread(x, 2, size(y)) * spread(y, 1, size(x))
   end function outer_product

end module marss_optim_mod
