! SPDX-License-Identifier: GPL-2.0-or-later
!
! Modern free-form Fortran translation of the computational core of the
! R package ucminf 1.2.3. The optimization algorithm and original Fortran
! code are by Hans Bruun Nielsen, IMM/DTU (2000). The R implementation is
! by Stig Bousgaard Mortensen, with later maintenance/modifications noted
! in upstream/DESCRIPTION and upstream/ucminf.R.
!
! This file preserves the UCMINF algorithm: inverse-Hessian BFGS updates,
! soft line search, trust-region-style step monitoring, analytic or finite-
! difference gradients, and the original convergence codes.
module ucminf
   use ucminf_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: UCMINF_GRAD_FORWARD = 1
   integer, parameter, public :: UCMINF_GRAD_CENTRAL = 2

   type, public :: ucminf_options
      logical :: trace = .false.
      real(dp) :: grtol = 1.0e-6_dp
      real(dp) :: xtol = 1.0e-12_dp
      real(dp) :: stepmax = 1.0_dp
      integer :: maxeval = 500
      integer :: grad_method = UCMINF_GRAD_FORWARD
      real(dp) :: gradstep(2) = [1.0e-6_dp, 1.0e-8_dp]
      real(dp), allocatable :: invhessian_lt(:)
   end type ucminf_options

   type, public :: ucminf_result
      real(dp), allocatable :: par(:)
      real(dp) :: value = 0.0_dp
      integer :: convergence = 0
      character(len=:), allocatable :: message
      real(dp), allocatable :: invhessian(:,:)
      real(dp), allocatable :: invhessian_lt(:)
      real(dp) :: maxgradient = 0.0_dp
      real(dp) :: laststep = 0.0_dp
      real(dp) :: stepmax = 0.0_dp
      integer :: neval = 0
   end type ucminf_result

   type, public :: ucminf_gradient_check
      real(dp) :: max_gradient = 0.0_dp
      real(dp) :: max_forward_error = 0.0_dp
      real(dp) :: max_backward_error = 0.0_dp
      real(dp) :: max_extrapolated_error = 0.0_dp
      integer :: forward_index = 0
      integer :: backward_index = 0
      integer :: extrapolated_index = 0
      logical :: success = .false.
   end type ucminf_gradient_check

   abstract interface
      function ucminf_objective(x) result(value)
         import :: dp
         real(dp), intent(in) :: x(:)
         real(dp) :: value
      end function ucminf_objective

      subroutine ucminf_gradient(x, g)
         import :: dp
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: g(:)
      end subroutine ucminf_gradient

      function ucminf_objective_context(x, context) result(value)
         import :: dp
         real(dp), intent(in) :: x(:)
         class(*), intent(inout) :: context
         real(dp) :: value
      end function ucminf_objective_context

      subroutine ucminf_gradient_context(x, g, context)
         import :: dp
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: g(:)
         class(*), intent(inout) :: context
      end subroutine ucminf_gradient_context
   end interface

   public :: dp
   public :: ucminf_minimize
   public :: ucminf_minimize_context
   public :: ucminf_check_gradient

contains

   subroutine ucminf_minimize(par, fn, result, gr, options)
      real(dp), intent(in) :: par(:)
      procedure(ucminf_objective) :: fn
      type(ucminf_result), intent(out) :: result
      procedure(ucminf_gradient), optional :: gr
      type(ucminf_options), intent(in), optional :: options

      type(ucminf_options) :: opt
      real(dp), allocatable :: x(:), g(:), g_old(:), step(:), s(:), y(:)
      real(dp), allocatable :: invh(:), chol_work(:), dy(:)
      real(dp) :: delta, fx, fx_new, nmg, nmh, nmx, alpha
      real(dp) :: slopes(2), threshold, ys, ydy, coef
      integer :: n, nn, icontr, neval, line_evals, fail
      logical :: given_hessian, reduce_step, use_delta

      opt = ucminf_options()
      if (present(options)) opt = options

      n = size(par)
      nn = n * (n + 1) / 2
      allocate(result%par(n), result%invhessian(n,n), result%invhessian_lt(nn))
      result%par = par
      result%invhessian = 0.0_dp
      result%invhessian_lt = 0.0_dp
      result%stepmax = opt%stepmax
      result%neval = 0

      if (n <= 0) then
         call set_failure(result, -2)
         return
      end if
      if (opt%stepmax <= 0.0_dp) then
         call set_failure(result, -4)
         return
      end if
      if (opt%grtol <= 0.0_dp .or. opt%xtol <= 0.0_dp) then
         call set_failure(result, -5)
         return
      end if
      if (opt%maxeval <= 0) then
         call set_failure(result, -6)
         return
      end if
      if (opt%grad_method /= UCMINF_GRAD_FORWARD .and. &
          opt%grad_method /= UCMINF_GRAD_CENTRAL) then
         call set_failure(result, -9, "Invalid finite-difference gradient method.")
         return
      end if
      if (any(opt%gradstep <= 0.0_dp)) then
         call set_failure(result, -10, "gradstep entries must be positive.")
         return
      end if

      allocate(x(n), g(n), g_old(n), step(n), s(n), y(n), dy(n), invh(nn))
      x = par
      delta = opt%stepmax
      given_hessian = allocated(opt%invhessian_lt)

      if (given_hessian) then
         if (size(opt%invhessian_lt) /= nn) then
            call set_failure(result, -11, "invhessian_lt has the wrong length.")
            return
         end if
         invh = opt%invhessian_lt
         allocate(chol_work(nn))
         chol_work = invh
         call packed_cholesky(n, chol_work, fail)
         deallocate(chol_work)
         if (fail /= 0) then
            call set_failure(result, -7)
            return
         end if
         use_delta = .false.
      else
         call packed_identity(n, invh)
         use_delta = .true.
      end if

      call evaluate_fg(x, fn, g, fx, gr, opt%grad_method, opt%gradstep)
      neval = 1
      nmh = 0.0_dp
      nmx = norm2(x)
      nmg = maxval(abs(g))
      icontr = 0

      if (nmg <= opt%grtol) then
         icontr = 1
      end if

      do while (icontr == 0)
         if (opt%trace) call print_trace(neval, fx, nmg, x)

         g_old = g
         call packed_sym_matvec(n, invh, g, step)
         step = -step

         reduce_step = .false.
         nmh = norm2(step)
         if (nmh <= opt%xtol * (opt%xtol + nmx)) then
            icontr = 2
            exit
         end if
         if (nmh > delta .or. use_delta) then
            reduce_step = .true.
            step = step * (delta / nmh)
            nmh = delta
            use_delta = .false.
         end if

         line_evals = 5
         call soft_line_search(x, fx, g, step, alpha, fx_new, slopes, line_evals, &
            fn, gr, opt%grad_method, opt%gradstep)
         if (opt%trace) call print_line(alpha, slopes)

         if (alpha <= 0.0_dp) then
            icontr = 4
            nmh = 0.0_dp
            exit
         end if

         neval = neval + line_evals
         nmg = maxval(abs(g))
         fx = fx_new
         s = alpha * step
         x = x + s
         nmx = norm2(x)
         nmh = norm2(s)

         if (alpha < 1.0_dp) then
            delta = 0.35_dp * delta
         else if (reduce_step) then
            if (slopes(2) < 0.7_dp * slopes(1)) delta = 3.0_dp * delta
         end if

         y = g - g_old
         ys = dot_product(y, s)
         if (ys > 1.0e-8_dp * nmh * norm2(y)) then
            call packed_sym_matvec(n, invh, y, dy)
            ydy = dot_product(y, dy)
            coef = (1.0_dp + ydy / ys) / ys
            call packed_rank2_update(n, invh, s, s, 0.5_dp * coef)
            call packed_rank2_update(n, invh, s, dy, -1.0_dp / ys)
         end if

         threshold = opt%xtol * (opt%xtol + nmx)
         delta = max(delta, threshold)
         if (neval >= opt%maxeval) icontr = 3
         if (nmh <= threshold) icontr = 2
         if (nmg <= opt%grtol) icontr = 1
      end do

      result%par = x
      result%value = fx
      result%convergence = icontr
      result%message = convergence_message(icontr)
      result%maxgradient = nmg
      result%laststep = nmh
      result%stepmax = delta
      result%neval = neval
      result%invhessian_lt = invh
      call unpack_lower_symmetric(n, invh, result%invhessian)

      if (opt%trace) then
         if (icontr == 1 .or. icontr == 2 .or. icontr == 4) then
            write(*,'(a)') " Optimization has converged"
         else
            write(*,'(a,i0,a)') " Optimization stopped after ", neval, " function evaluations"
         end if
      end if
   end subroutine ucminf_minimize

   subroutine evaluate_fg(x, fn, g, f, gr, grad_method, gradstep)
      real(dp), intent(in) :: x(:)
      procedure(ucminf_objective) :: fn
      real(dp), intent(out) :: g(:)
      real(dp), intent(out) :: f
      procedure(ucminf_gradient), optional :: gr
      integer, intent(in) :: grad_method
      real(dp), intent(in) :: gradstep(2)

      f = fn(x)
      if (present(gr)) then
         call gr(x, g)
      else
         call numerical_gradient(x, f, fn, g, grad_method, gradstep)
      end if
   end subroutine evaluate_fg

   subroutine numerical_gradient(x, f, fn, g, grad_method, gradstep)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: f
      procedure(ucminf_objective) :: fn
      real(dp), intent(out) :: g(:)
      integer, intent(in) :: grad_method
      real(dp), intent(in) :: gradstep(2)

      real(dp) :: x2(size(x)), dx, f2, f3
      integer :: i

      do i = 1, size(x)
         x2 = x
         dx = abs(x2(i)) * gradstep(1) + gradstep(2)
         x2(i) = x2(i) + dx
         f2 = fn(x2)
         if (grad_method == UCMINF_GRAD_FORWARD) then
            g(i) = (f2 - f) / dx
         else
            x2(i) = x2(i) - 2.0_dp * dx
            f3 = fn(x2)
            g(i) = (f2 - f3) / (2.0_dp * dx)
         end if
      end do
   end subroutine numerical_gradient

   subroutine soft_line_search(x, f, g, h, alpha, fn_new, slopes, neval, fn, gr, grad_method, gradstep)
      real(dp), intent(in) :: x(:), f, h(:)
      real(dp), intent(inout) :: g(:)
      real(dp), intent(out) :: alpha, fn_new, slopes(2)
      integer, intent(inout) :: neval
      procedure(ucminf_objective) :: fn
      procedure(ucminf_gradient), optional :: gr
      integer, intent(in) :: grad_method
      real(dp), intent(in) :: gradstep(2)

      logical :: ok
      integer :: maxeval
      real(dp) :: xfd(3,3), trial_x(size(x)), trial_g(size(x))
      real(dp) :: a, b, c, d, fi0, sl0, slthr, trial_f

      alpha = 0.0_dp
      fn_new = f
      maxeval = neval
      neval = 0
      slopes(1) = dot_product(g, h)
      slopes(2) = slopes(1)
      if (slopes(1) >= 0.0_dp) return

      fi0 = f
      sl0 = 0.05_dp * slopes(1)
      slthr = 0.995_dp * slopes(1)
      ok = .false.
      xfd = 0.0_dp
      xfd(1,1) = 0.0_dp
      xfd(2,1) = f
      xfd(3,1) = slopes(1)
      b = 1.0_dp

      do
         xfd(1,2) = b
         trial_x = x + b * h
         call evaluate_fg(trial_x, fn, trial_g, xfd(2,2), gr, grad_method, gradstep)
         neval = neval + 1
         xfd(3,2) = dot_product(trial_g, h)
         if (b < 1.5_dp) slopes(2) = xfd(3,2)

         if (xfd(2,2) <= fi0 + sl0 * xfd(1,2)) then
            if (xfd(3,2) <= abs(slthr)) then
               ok = .true.
               alpha = xfd(1,2)
               fn_new = xfd(2,2)
               slopes(2) = xfd(3,2)
               g = trial_g
               if (b < 2.0_dp .and. xfd(3,2) < slthr) then
                  xfd(:,1) = xfd(:,2)
                  b = 2.0_dp
                  cycle
               end if
            end if
         end if

         d = xfd(1,2) - xfd(1,1)
         exit
      end do

      do
         if (ok .or. neval == maxeval) return

         c = xfd(2,2) - xfd(2,1) - d * xfd(3,1)
         if (c > 1.0e-15_dp * real(size(x), dp) * xfd(1,2)) then
            a = xfd(1,1) - 0.5_dp * xfd(3,1) * (d**2 / c)
            d = 0.1_dp * d
            xfd(1,3) = min(max(xfd(1,1) + d, a), xfd(1,2) - d)
         else
            xfd(1,3) = 0.5_dp * (xfd(1,1) + xfd(1,2))
         end if

         trial_x = x + xfd(1,3) * h
         call evaluate_fg(trial_x, fn, trial_g, trial_f, gr, grad_method, gradstep)
         xfd(2,3) = trial_f
         neval = neval + 1
         xfd(3,3) = dot_product(trial_g, h)

         if (xfd(2,3) < fi0 + sl0 * xfd(1,3)) then
            ok = .true.
            alpha = xfd(1,3)
            fn_new = xfd(2,3)
            slopes(2) = xfd(3,3)
            g = trial_g
            xfd(:,1) = xfd(:,3)
         else
            xfd(:,2) = xfd(:,3)
         end if

         d = xfd(1,2) - xfd(1,1)
         ok = ok .and. (abs(xfd(3,3)) <= abs(slthr))
         ok = ok .or. (d <= 0.0_dp)
      end do
   end subroutine soft_line_search

   subroutine ucminf_minimize_context(par, fn, context, result, gr, options)
      real(dp), intent(in) :: par(:)
      procedure(ucminf_objective_context) :: fn
      type(ucminf_result), intent(out) :: result
      procedure(ucminf_gradient_context), optional :: gr
      type(ucminf_options), intent(in), optional :: options

      class(*), intent(inout) :: context

      type(ucminf_options) :: opt
      real(dp), allocatable :: x(:), g(:), g_old(:), step(:), s(:), y(:)
      real(dp), allocatable :: invh(:), chol_work(:), dy(:)
      real(dp) :: delta, fx, fx_new, nmg, nmh, nmx, alpha
      real(dp) :: slopes(2), threshold, ys, ydy, coef
      integer :: n, nn, icontr, neval, line_evals, fail
      logical :: given_hessian, reduce_step, use_delta

      opt = ucminf_options()
      if (present(options)) opt = options

      n = size(par)
      nn = n * (n + 1) / 2
      allocate(result%par(n), result%invhessian(n,n), result%invhessian_lt(nn))
      result%par = par
      result%invhessian = 0.0_dp
      result%invhessian_lt = 0.0_dp
      result%stepmax = opt%stepmax
      result%neval = 0

      if (n <= 0) then
         call set_failure(result, -2)
         return
      end if
      if (opt%stepmax <= 0.0_dp) then
         call set_failure(result, -4)
         return
      end if
      if (opt%grtol <= 0.0_dp .or. opt%xtol <= 0.0_dp) then
         call set_failure(result, -5)
         return
      end if
      if (opt%maxeval <= 0) then
         call set_failure(result, -6)
         return
      end if
      if (opt%grad_method /= UCMINF_GRAD_FORWARD .and. &
          opt%grad_method /= UCMINF_GRAD_CENTRAL) then
         call set_failure(result, -9, "Invalid finite-difference gradient method.")
         return
      end if
      if (any(opt%gradstep <= 0.0_dp)) then
         call set_failure(result, -10, "gradstep entries must be positive.")
         return
      end if

      allocate(x(n), g(n), g_old(n), step(n), s(n), y(n), dy(n), invh(nn))
      x = par
      delta = opt%stepmax
      given_hessian = allocated(opt%invhessian_lt)

      if (given_hessian) then
         if (size(opt%invhessian_lt) /= nn) then
            call set_failure(result, -11, "invhessian_lt has the wrong length.")
            return
         end if
         invh = opt%invhessian_lt
         allocate(chol_work(nn))
         chol_work = invh
         call packed_cholesky(n, chol_work, fail)
         deallocate(chol_work)
         if (fail /= 0) then
            call set_failure(result, -7)
            return
         end if
         use_delta = .false.
      else
         call packed_identity(n, invh)
         use_delta = .true.
      end if

      call evaluate_fg_context(x, fn, context, g, fx, gr, opt%grad_method, opt%gradstep)
      neval = 1
      nmh = 0.0_dp
      nmx = norm2(x)
      nmg = maxval(abs(g))
      icontr = 0

      if (nmg <= opt%grtol) then
         icontr = 1
      end if

      do while (icontr == 0)
         if (opt%trace) call print_trace(neval, fx, nmg, x)

         g_old = g
         call packed_sym_matvec(n, invh, g, step)
         step = -step

         reduce_step = .false.
         nmh = norm2(step)
         if (nmh <= opt%xtol * (opt%xtol + nmx)) then
            icontr = 2
            exit
         end if
         if (nmh > delta .or. use_delta) then
            reduce_step = .true.
            step = step * (delta / nmh)
            nmh = delta
            use_delta = .false.
         end if

         line_evals = 5
         call soft_line_search_context(x, fx, g, step, alpha, fx_new, slopes, line_evals, &
            fn, context, gr, opt%grad_method, opt%gradstep)
         if (opt%trace) call print_line(alpha, slopes)

         if (alpha <= 0.0_dp) then
            icontr = 4
            nmh = 0.0_dp
            exit
         end if

         neval = neval + line_evals
         nmg = maxval(abs(g))
         fx = fx_new
         s = alpha * step
         x = x + s
         nmx = norm2(x)
         nmh = norm2(s)

         if (alpha < 1.0_dp) then
            delta = 0.35_dp * delta
         else if (reduce_step) then
            if (slopes(2) < 0.7_dp * slopes(1)) delta = 3.0_dp * delta
         end if

         y = g - g_old
         ys = dot_product(y, s)
         if (ys > 1.0e-8_dp * nmh * norm2(y)) then
            call packed_sym_matvec(n, invh, y, dy)
            ydy = dot_product(y, dy)
            coef = (1.0_dp + ydy / ys) / ys
            call packed_rank2_update(n, invh, s, s, 0.5_dp * coef)
            call packed_rank2_update(n, invh, s, dy, -1.0_dp / ys)
         end if

         threshold = opt%xtol * (opt%xtol + nmx)
         delta = max(delta, threshold)
         if (neval >= opt%maxeval) icontr = 3
         if (nmh <= threshold) icontr = 2
         if (nmg <= opt%grtol) icontr = 1
      end do

      result%par = x
      result%value = fx
      result%convergence = icontr
      result%message = convergence_message(icontr)
      result%maxgradient = nmg
      result%laststep = nmh
      result%stepmax = delta
      result%neval = neval
      result%invhessian_lt = invh
      call unpack_lower_symmetric(n, invh, result%invhessian)

      if (opt%trace) then
         if (icontr == 1 .or. icontr == 2 .or. icontr == 4) then
            write(*,'(a)') " Optimization has converged"
         else
            write(*,'(a,i0,a)') " Optimization stopped after ", neval, " function evaluations"
         end if
      end if
   end subroutine ucminf_minimize_context

   subroutine evaluate_fg_context(x, fn, context, g, f, gr, grad_method, gradstep)
      real(dp), intent(in) :: x(:)
      procedure(ucminf_objective_context) :: fn
      class(*), intent(inout) :: context
      real(dp), intent(out) :: g(:)
      real(dp), intent(out) :: f
      procedure(ucminf_gradient_context), optional :: gr
      integer, intent(in) :: grad_method
      real(dp), intent(in) :: gradstep(2)

      f = fn(x, context)
      if (present(gr)) then
         call gr(x, g, context)
      else
         call numerical_gradient_context(x, f, fn, context, g, grad_method, gradstep)
      end if
   end subroutine evaluate_fg_context

   subroutine numerical_gradient_context(x, f, fn, context, g, grad_method, gradstep)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: f
      procedure(ucminf_objective_context) :: fn
      class(*), intent(inout) :: context
      real(dp), intent(out) :: g(:)
      integer, intent(in) :: grad_method
      real(dp), intent(in) :: gradstep(2)

      real(dp) :: x2(size(x)), dx, f2, f3
      integer :: i

      do i = 1, size(x)
         x2 = x
         dx = abs(x2(i)) * gradstep(1) + gradstep(2)
         x2(i) = x2(i) + dx
         f2 = fn(x2, context)
         if (grad_method == UCMINF_GRAD_FORWARD) then
            g(i) = (f2 - f) / dx
         else
            x2(i) = x2(i) - 2.0_dp * dx
            f3 = fn(x2, context)
            g(i) = (f2 - f3) / (2.0_dp * dx)
         end if
      end do
   end subroutine numerical_gradient_context

   subroutine soft_line_search_context(x, f, g, h, alpha, fn_new, slopes, neval, fn, context, gr, grad_method, gradstep)
      real(dp), intent(in) :: x(:), f, h(:)
      real(dp), intent(inout) :: g(:)
      real(dp), intent(out) :: alpha, fn_new, slopes(2)
      integer, intent(inout) :: neval
      procedure(ucminf_objective_context) :: fn
      class(*), intent(inout) :: context
      procedure(ucminf_gradient_context), optional :: gr
      integer, intent(in) :: grad_method
      real(dp), intent(in) :: gradstep(2)

      logical :: ok
      integer :: maxeval
      real(dp) :: xfd(3,3), trial_x(size(x)), trial_g(size(x))
      real(dp) :: a, b, c, d, fi0, sl0, slthr, trial_f

      alpha = 0.0_dp
      fn_new = f
      maxeval = neval
      neval = 0
      slopes(1) = dot_product(g, h)
      slopes(2) = slopes(1)
      if (slopes(1) >= 0.0_dp) return

      fi0 = f
      sl0 = 0.05_dp * slopes(1)
      slthr = 0.995_dp * slopes(1)
      ok = .false.
      xfd = 0.0_dp
      xfd(1,1) = 0.0_dp
      xfd(2,1) = f
      xfd(3,1) = slopes(1)
      b = 1.0_dp

      do
         xfd(1,2) = b
         trial_x = x + b * h
         call evaluate_fg_context(trial_x, fn, context, trial_g, xfd(2,2), gr, grad_method, gradstep)
         neval = neval + 1
         xfd(3,2) = dot_product(trial_g, h)
         if (b < 1.5_dp) slopes(2) = xfd(3,2)

         if (xfd(2,2) <= fi0 + sl0 * xfd(1,2)) then
            if (xfd(3,2) <= abs(slthr)) then
               ok = .true.
               alpha = xfd(1,2)
               fn_new = xfd(2,2)
               slopes(2) = xfd(3,2)
               g = trial_g
               if (b < 2.0_dp .and. xfd(3,2) < slthr) then
                  xfd(:,1) = xfd(:,2)
                  b = 2.0_dp
                  cycle
               end if
            end if
         end if

         d = xfd(1,2) - xfd(1,1)
         exit
      end do

      do
         if (ok .or. neval == maxeval) return

         c = xfd(2,2) - xfd(2,1) - d * xfd(3,1)
         if (c > 1.0e-15_dp * real(size(x), dp) * xfd(1,2)) then
            a = xfd(1,1) - 0.5_dp * xfd(3,1) * (d**2 / c)
            d = 0.1_dp * d
            xfd(1,3) = min(max(xfd(1,1) + d, a), xfd(1,2) - d)
         else
            xfd(1,3) = 0.5_dp * (xfd(1,1) + xfd(1,2))
         end if

         trial_x = x + xfd(1,3) * h
         call evaluate_fg_context(trial_x, fn, context, trial_g, trial_f, gr, grad_method, gradstep)
         xfd(2,3) = trial_f
         neval = neval + 1
         xfd(3,3) = dot_product(trial_g, h)

         if (xfd(2,3) < fi0 + sl0 * xfd(1,3)) then
            ok = .true.
            alpha = xfd(1,3)
            fn_new = xfd(2,3)
            slopes(2) = xfd(3,3)
            g = trial_g
            xfd(:,1) = xfd(:,3)
         else
            xfd(:,2) = xfd(:,3)
         end if

         d = xfd(1,2) - xfd(1,1)
         ok = ok .and. (abs(xfd(3,3)) <= abs(slthr))
         ok = ok .or. (d <= 0.0_dp)
      end do
   end subroutine soft_line_search_context

   subroutine ucminf_check_gradient(x, fn, gr, step, report)
      real(dp), intent(inout) :: x(:)
      procedure(ucminf_objective) :: fn
      procedure(ucminf_gradient) :: gr
      real(dp), intent(in) :: step
      type(ucminf_gradient_check), intent(out) :: report

      real(dp) :: g(size(x)), g1(size(x)), f, f1, xi, h, af, ab, ae, er
      integer :: i

      report = ucminf_gradient_check()
      if (size(x) <= 0 .or. abs(step) <= tiny(1.0_dp)) return

      f = fn(x)
      call gr(x, g)
      report%max_gradient = maxval(abs(g))

      do i = 1, size(x)
         xi = x(i)

         x(i) = xi + step
         h = x(i) - xi
         if (abs(h) <= tiny(1.0_dp)) then
            x(i) = xi
            return
         end if
         f1 = fn(x)
         call gr(x, g1)
         af = (f1 - f) / h
         er = af - g(i)
         if (abs(er) > abs(report%max_forward_error)) then
            report%max_forward_error = er
            report%forward_index = i
         end if

         x(i) = xi - 0.5_dp * step
         h = x(i) - xi
         if (abs(h) <= tiny(1.0_dp)) then
            x(i) = xi
            return
         end if
         f1 = fn(x)
         call gr(x, g1)
         ab = (f1 - f) / h
         er = ab - g(i)
         if (abs(er) > abs(report%max_backward_error)) then
            report%max_backward_error = er
            report%backward_index = i
         end if

         ae = (2.0_dp * ab + af) / 3.0_dp
         er = ae - g(i)
         if (abs(er) > abs(report%max_extrapolated_error)) then
            report%max_extrapolated_error = er
            report%extrapolated_index = i
         end if

         x(i) = xi
      end do
      report%success = .true.
   end subroutine ucminf_check_gradient

   subroutine packed_identity(n, a)
      integer, intent(in) :: n
      real(dp), intent(out) :: a(:)
      integer :: i, j, k

      k = 0
      do j = 1, n
         do i = j, n
            k = k + 1
            if (i == j) then
               a(k) = 1.0_dp
            else
               a(k) = 0.0_dp
            end if
         end do
      end do
   end subroutine packed_identity

   subroutine packed_sym_matvec(n, a, x, y)
      integer, intent(in) :: n
      real(dp), intent(in) :: a(:), x(:)
      real(dp), intent(out) :: y(:)
      integer :: i, j, k

      y = 0.0_dp
      k = 0
      do j = 1, n
         do i = j, n
            k = k + 1
            y(i) = y(i) + a(k) * x(j)
            if (i /= j) y(j) = y(j) + a(k) * x(i)
         end do
      end do
   end subroutine packed_sym_matvec

   subroutine packed_rank2_update(n, a, x, y, alpha)
      integer, intent(in) :: n
      real(dp), intent(inout) :: a(:)
      real(dp), intent(in) :: x(:), y(:), alpha
      integer :: i, j, k

      k = 0
      do j = 1, n
         do i = j, n
            k = k + 1
            a(k) = a(k) + alpha * (x(i) * y(j) + y(i) * x(j))
         end do
      end do
   end subroutine packed_rank2_update

   subroutine packed_cholesky(n, a, fail)
      integer, intent(in) :: n
      real(dp), intent(inout) :: a(:)
      integer, intent(out) :: fail

      real(dp) :: full(n,n), l(n,n), s
      integer :: i, j, k

      call unpack_lower_symmetric(n, a, full)
      l = 0.0_dp
      fail = 0
      do j = 1, n
         do i = j, n
            s = full(i,j)
            if (j > 1) s = s - dot_product(l(i,1:j-1), l(j,1:j-1))
            if (i == j) then
               if (s <= 0.0_dp) then
                  fail = j
                  return
               end if
               l(j,j) = sqrt(s)
            else
               l(i,j) = s / l(j,j)
            end if
         end do
      end do

      k = 0
      do j = 1, n
         do i = j, n
            k = k + 1
            a(k) = l(i,j)
         end do
      end do
   end subroutine packed_cholesky

   subroutine unpack_lower_symmetric(n, packed, full)
      integer, intent(in) :: n
      real(dp), intent(in) :: packed(:)
      real(dp), intent(out) :: full(:,:)
      integer :: i, j, k

      full = 0.0_dp
      k = 0
      do j = 1, n
         do i = j, n
            k = k + 1
            full(i,j) = packed(k)
            full(j,i) = packed(k)
         end do
      end do
   end subroutine unpack_lower_symmetric

   subroutine set_failure(result, code, custom_message)
      type(ucminf_result), intent(inout) :: result
      integer, intent(in) :: code
      character(len=*), intent(in), optional :: custom_message

      result%convergence = code
      if (present(custom_message)) then
         result%message = custom_message
      else
         result%message = convergence_message(code)
      end if
   end subroutine set_failure

   function convergence_message(code) result(message)
      integer, intent(in) :: code
      character(len=:), allocatable :: message

      select case (code)
      case (1)
         message = "Stopped by small gradient (grtol)."
      case (2)
         message = "Stopped by small step (xtol)."
      case (3)
         message = "Stopped by function evaluation limit (maxeval)."
      case (4)
         message = "Stopped by zero step from line search."
      case (-2)
         message = "Computation did not start: length(par) = 0."
      case (-4)
         message = "Computation did not start: stepmax is too small."
      case (-5)
         message = "Computation did not start: grtol or xtol <= 0."
      case (-6)
         message = "Computation did not start: maxeval <= 0."
      case (-7)
         message = "Computation did not start: given Hessian not positive definite."
      case default
         message = "Unknown UCMINF status."
      end select
   end function convergence_message

   subroutine print_trace(neval, fx, nmg, x)
      integer, intent(in) :: neval
      real(dp), intent(in) :: fx, nmg, x(:)
      integer :: i

      write(*,'(a,i0,a,es12.4,a,es12.4)') " neval = ", neval, ", F(x) = ", fx, ", max|g(x)| = ", nmg
      write(*,'(a)',advance='no') " x = "
      do i = 1, size(x)
         if (i > 1) write(*,'(a)',advance='no') ", "
         write(*,'(es12.4)',advance='no') x(i)
      end do
      write(*,*)
   end subroutine print_trace

   subroutine print_line(alpha, slopes)
      real(dp), intent(in) :: alpha, slopes(2)
      write(*,'(a,es12.4,a,es12.4,a,es12.4)') &
         " Line search: alpha = ", alpha, ", dphi(0) = ", slopes(1), ", dphi(alpha) = ", slopes(2)
   end subroutine print_line

end module ucminf
