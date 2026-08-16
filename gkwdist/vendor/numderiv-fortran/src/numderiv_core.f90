! SPDX-License-Identifier: GPL-2.0-or-later
module numderiv_core
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
   use numderiv_kinds, only : dp
   use numderiv_types, only : deriv_options, gend_result, first_deriv_options, hessian_options, &
      nd_success, nd_invalid_argument, nd_nonfinite_value, nd_shape_mismatch
   use numderiv_callbacks, only : scalar_real_function, vector_real_function, &
      scalar_complex_function, vector_complex_function
   implicit none
   private

   public :: grad, grad_elementwise, grad_complex, grad_elementwise_complex
   public :: jacobian, jacobian_complex, hessian, hessian_complex, gend

contains

   subroutine grad(func, x, gradient, method, side, options, status, message)
      procedure(scalar_real_function) :: func
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: gradient(:)
      character(len=*), intent(in), optional :: method
      integer, intent(in), optional :: side(:)
      type(deriv_options), intent(in), optional :: options
      integer, intent(out), optional :: status
      character(:), allocatable, intent(out), optional :: message

      type(deriv_options) :: opts
      character(:), allocatable :: meth
      integer, allocatable :: direction(:)
      real(dp), allocatable :: a(:, :), h(:), xp(:), xm(:)
      real(dp) :: f0, fp, fm, factor
      integer :: i, k, m, n, stat
      character(:), allocatable :: msg

      n = size(x)
      gradient = quiet_nan()
      call begin_status(status, message)
      if (size(gradient) /= n) then
         call finish_status(nd_shape_mismatch, 'gradient must have size(x)', status, message)
         return
      end if

      opts = first_deriv_options()
      if (present(options)) opts = options
      call validate_options(opts, stat, msg)
      if (stat /= nd_success) then
         call finish_status(stat, msg, status, message)
         return
      end if
      call make_side(n, side, direction, stat, msg)
      if (stat /= nd_success) then
         call finish_status(stat, msg, status, message)
         return
      end if

      meth = 'richardson'
      if (present(method)) meth = lower(trim(method))
      select case (meth)
      case ('simple')
         f0 = func(x)
         if (.not. ieee_is_finite(f0)) then
            call finish_status(nd_nonfinite_value, 'function returned a nonfinite value at x', &
               status, message)
            return
         end if
         allocate(xp(n))
         do i = 1, n
            xp = x
            if (direction(i) == 0) direction(i) = 1
            xp(i) = xp(i) + opts%eps * real(direction(i), dp)
            fp = func(xp)
            if (.not. ieee_is_finite(fp)) then
               call finish_status(nd_nonfinite_value, &
                  'function returned a nonfinite value in simple difference', status, message)
               return
            end if
            gradient(i) = (fp - f0) / (opts%eps * real(direction(i), dp))
         end do

      case ('richardson')
         allocate(a(opts%r, n), h(n), xp(n), xm(n))
         h = initial_steps(x, opts)
         do k = 1, opts%r
            do i = 1, n
               xp = x
               xm = x
               select case (direction(i))
               case (1)
                  xp(i) = x(i) + 2.0_dp * h(i)
               case (-1)
                  xm(i) = x(i) - 2.0_dp * h(i)
               case default
                  xp(i) = x(i) + h(i)
                  xm(i) = x(i) - h(i)
               end select
               if (k > 1) then
                  if (abs(a(k - 1, i)) < 1.0e-20_dp) then
                     a(k, i) = 0.0_dp
                     cycle
                  end if
               end if
               fp = func(xp)
               fm = func(xm)
               if (.not. finite_pair(fp, fm)) then
                  call finish_status(nd_nonfinite_value, &
                     'function returned a nonfinite value in Richardson difference', &
                     status, message)
                  return
               end if
               a(k, i) = (fp - fm) / (2.0_dp * h(i))
            end do
            if (opts%show_details) then
               write(*, '(a,i0)') 'first-order approximation level ', k
               write(*, '(*(es20.12,1x))') a(k, :)
            end if
            h = h / opts%v
         end do
         do m = 1, opts%r - 1
            factor = opts%v ** (2 * m)
            do k = 1, opts%r - m
               a(k, :) = (factor * a(k + 1, :) - a(k, :)) / (factor - 1.0_dp)
            end do
            if (opts%show_details) then
               write(*, '(a,i0)') 'Richardson improvement group ', m
               do k = 1, opts%r - m
                  write(*, '(*(es20.12,1x))') a(k, :)
               end do
            end if
         end do
         gradient = a(1, :)

      case default
         call finish_status(nd_invalid_argument, &
            'method must be "simple" or "Richardson"; use grad_complex for complex step', &
            status, message)
         return
      end select
      call finish_status(nd_success, '', status, message)
   end subroutine grad


   subroutine grad_elementwise(func, x, gradient, method, side, options, status, message)
      procedure(vector_real_function) :: func
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: gradient(:)
      character(len=*), intent(in), optional :: method
      integer, intent(in), optional :: side(:)
      type(deriv_options), intent(in), optional :: options
      integer, intent(out), optional :: status
      character(:), allocatable, intent(out), optional :: message

      type(deriv_options) :: opts
      character(:), allocatable :: meth
      integer, allocatable :: direction(:)
      real(dp), allocatable :: f0(:), fp(:), fm(:), a(:, :), h(:), ph(:), mh(:)
      real(dp) :: factor
      integer :: k, m, n, stat
      character(:), allocatable :: msg

      n = size(x)
      gradient = quiet_nan()
      call begin_status(status, message)
      if (size(gradient) /= n) then
         call finish_status(nd_shape_mismatch, 'gradient must have size(x)', status, message)
         return
      end if
      opts = first_deriv_options()
      if (present(options)) opts = options
      call validate_options(opts, stat, msg)
      if (stat /= nd_success) then
         call finish_status(stat, msg, status, message)
         return
      end if
      call make_side(n, side, direction, stat, msg)
      if (stat /= nd_success) then
         call finish_status(stat, msg, status, message)
         return
      end if
      meth = 'richardson'
      if (present(method)) meth = lower(trim(method))

      f0 = func(x)
      if (size(f0) /= n) then
         call finish_status(nd_shape_mismatch, &
            'grad_elementwise requires a result with size(x)', status, message)
         return
      end if
      if (.not. all(ieee_is_finite(f0))) then
         call finish_status(nd_nonfinite_value, 'function returned a nonfinite value at x', &
            status, message)
         return
      end if

      select case (meth)
      case ('simple')
         where (direction == 0) direction = 1
         fp = func(x + opts%eps * real(direction, dp))
         if (size(fp) /= n) then
            call finish_status(nd_shape_mismatch, 'function result size changed', status, message)
            return
         end if
         if (.not. all(ieee_is_finite(fp))) then
            call finish_status(nd_nonfinite_value, &
               'function returned a nonfinite value in simple difference', status, message)
            return
         end if
         gradient = (fp - f0) / (opts%eps * real(direction, dp))

      case ('richardson')
         allocate(a(opts%r, n), h(n), ph(n), mh(n))
         h = initial_steps(x, opts)
         do k = 1, opts%r
            ph = h
            mh = h
            where (direction == 1)
               ph = 2.0_dp * h
               mh = 0.0_dp
            elsewhere (direction == -1)
               ph = 0.0_dp
               mh = 2.0_dp * h
            end where
            fp = func(x + ph)
            fm = func(x - mh)
            if (size(fp) /= n .or. size(fm) /= n) then
               call finish_status(nd_shape_mismatch, 'function result size changed', status, message)
               return
            end if
            if (.not. all(ieee_is_finite(fp)) .or. .not. all(ieee_is_finite(fm))) then
               call finish_status(nd_nonfinite_value, &
                  'function returned a nonfinite value in Richardson difference', status, message)
               return
            end if
            a(k, :) = (fp - fm) / (2.0_dp * h)
            h = h / opts%v
         end do
         do m = 1, opts%r - 1
            factor = opts%v ** (2 * m)
            do k = 1, opts%r - m
               a(k, :) = (factor * a(k + 1, :) - a(k, :)) / (factor - 1.0_dp)
            end do
         end do
         gradient = a(1, :)

      case default
         call finish_status(nd_invalid_argument, &
            'method must be "simple" or "Richardson"; use grad_elementwise_complex for complex step', &
            status, message)
         return
      end select
      call finish_status(nd_success, '', status, message)
   end subroutine grad_elementwise


   subroutine grad_complex(func, x, gradient, status, message)
      procedure(scalar_complex_function) :: func
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: gradient(:)
      integer, intent(out), optional :: status
      character(:), allocatable, intent(out), optional :: message

      complex(dp), allocatable :: z(:)
      complex(dp) :: value
      real(dp), parameter :: step = epsilon(1.0_dp)
      integer :: i, n

      n = size(x)
      gradient = quiet_nan()
      call begin_status(status, message)
      if (size(gradient) /= n) then
         call finish_status(nd_shape_mismatch, 'gradient must have size(x)', status, message)
         return
      end if
      allocate(z(n))
      z = cmplx(x, 0.0_dp, dp)
      do i = 1, n
         z(i) = cmplx(x(i), step, dp)
         value = func(z)
         z(i) = cmplx(x(i), 0.0_dp, dp)
         if (.not. finite_complex(value)) then
            call finish_status(nd_nonfinite_value, &
               'complex function returned a nonfinite value', status, message)
            return
         end if
         gradient(i) = aimag(value) / step
      end do
      call finish_status(nd_success, '', status, message)
   end subroutine grad_complex


   subroutine grad_elementwise_complex(func, x, gradient, status, message)
      procedure(vector_complex_function) :: func
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: gradient(:)
      integer, intent(out), optional :: status
      character(:), allocatable, intent(out), optional :: message

      complex(dp), allocatable :: z(:), value(:)
      real(dp), parameter :: step = epsilon(1.0_dp)
      integer :: n

      n = size(x)
      gradient = quiet_nan()
      call begin_status(status, message)
      if (size(gradient) /= n) then
         call finish_status(nd_shape_mismatch, 'gradient must have size(x)', status, message)
         return
      end if
      allocate(z(n))
      z = cmplx(x, step, dp)
      value = func(z)
      if (size(value) /= n) then
         call finish_status(nd_shape_mismatch, &
            'grad_elementwise_complex requires a result with size(x)', status, message)
         return
      end if
      if (.not. all_finite_complex(value)) then
         call finish_status(nd_nonfinite_value, &
            'complex function returned a nonfinite value', status, message)
         return
      end if
      gradient = aimag(value) / step
      call finish_status(nd_success, '', status, message)
   end subroutine grad_elementwise_complex


   subroutine jacobian(func, x, jac, method, side, options, status, message)
      procedure(vector_real_function) :: func
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: jac(:, :)
      character(len=*), intent(in), optional :: method
      integer, intent(in), optional :: side(:)
      type(deriv_options), intent(in), optional :: options
      integer, intent(out), optional :: status
      character(:), allocatable, intent(out), optional :: message

      type(deriv_options) :: opts
      character(:), allocatable :: meth
      integer, allocatable :: direction(:)
      real(dp), allocatable :: f0(:), fp(:), fm(:), a(:, :, :), h(:), xp(:), xm(:)
      real(dp) :: factor
      integer :: i, k, m, n, nout, stat
      character(:), allocatable :: msg

      call begin_status(status, message)
      n = size(x)
      f0 = func(x)
      nout = size(f0)
      allocate(jac(nout, n))
      jac = quiet_nan()
      if (.not. all(ieee_is_finite(f0))) then
         call finish_status(nd_nonfinite_value, 'function returned a nonfinite value at x', &
            status, message)
         return
      end if
      opts = first_deriv_options()
      if (present(options)) opts = options
      call validate_options(opts, stat, msg)
      if (stat /= nd_success) then
         call finish_status(stat, msg, status, message)
         return
      end if
      call make_side(n, side, direction, stat, msg)
      if (stat /= nd_success) then
         call finish_status(stat, msg, status, message)
         return
      end if
      meth = 'richardson'
      if (present(method)) meth = lower(trim(method))

      select case (meth)
      case ('simple')
         allocate(xp(n))
         do i = 1, n
            xp = x
            if (direction(i) == 0) direction(i) = 1
            xp(i) = xp(i) + opts%eps * real(direction(i), dp)
            fp = func(xp)
            if (size(fp) /= nout) then
               call finish_status(nd_shape_mismatch, 'function result size changed', status, message)
               return
            end if
            if (.not. all(ieee_is_finite(fp))) then
               call finish_status(nd_nonfinite_value, &
                  'function returned a nonfinite value in simple difference', status, message)
               return
            end if
            jac(:, i) = (fp - f0) / (opts%eps * real(direction(i), dp))
         end do

      case ('richardson')
         allocate(a(nout, opts%r, n), h(n), xp(n), xm(n))
         h = initial_steps(x, opts)
         do k = 1, opts%r
            do i = 1, n
               xp = x
               xm = x
               select case (direction(i))
               case (1)
                  xp(i) = x(i) + 2.0_dp * h(i)
               case (-1)
                  xm(i) = x(i) - 2.0_dp * h(i)
               case default
                  xp(i) = x(i) + h(i)
                  xm(i) = x(i) - h(i)
               end select
               fp = func(xp)
               fm = func(xm)
               if (size(fp) /= nout .or. size(fm) /= nout) then
                  call finish_status(nd_shape_mismatch, 'function result size changed', status, message)
                  return
               end if
               if (.not. all(ieee_is_finite(fp)) .or. .not. all(ieee_is_finite(fm))) then
                  call finish_status(nd_nonfinite_value, &
                     'function returned a nonfinite value in Richardson difference', &
                     status, message)
                  return
               end if
               a(:, k, i) = (fp - fm) / (2.0_dp * h(i))
            end do
            h = h / opts%v
         end do
         do m = 1, opts%r - 1
            factor = opts%v ** (2 * m)
            do k = 1, opts%r - m
               a(:, k, :) = (factor * a(:, k + 1, :) - a(:, k, :)) / (factor - 1.0_dp)
            end do
         end do
         jac = a(:, 1, :)

      case default
         call finish_status(nd_invalid_argument, &
            'method must be "simple" or "Richardson"; use jacobian_complex for complex step', &
            status, message)
         return
      end select
      call finish_status(nd_success, '', status, message)
   end subroutine jacobian


   subroutine jacobian_complex(func, x, jac, status, message)
      procedure(vector_complex_function) :: func
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: jac(:, :)
      integer, intent(out), optional :: status
      character(:), allocatable, intent(out), optional :: message

      complex(dp), allocatable :: z(:), value(:)
      real(dp), parameter :: step = epsilon(1.0_dp)
      integer :: i, n, nout

      call begin_status(status, message)
      n = size(x)
      allocate(z(n))
      z = cmplx(x, 0.0_dp, dp)
      if (n > 0) z(1) = cmplx(x(1), step, dp)
      value = func(z)
      nout = size(value)
      allocate(jac(nout, n))
      jac = quiet_nan()
      if (.not. all_finite_complex(value)) then
         call finish_status(nd_nonfinite_value, &
            'complex function returned a nonfinite value', status, message)
         return
      end if
      if (n == 0) then
         call finish_status(nd_success, '', status, message)
         return
      end if
      jac(:, 1) = aimag(value) / step
      z(1) = cmplx(x(1), 0.0_dp, dp)
      do i = 2, n
         z(i) = cmplx(x(i), step, dp)
         value = func(z)
         z(i) = cmplx(x(i), 0.0_dp, dp)
         if (size(value) /= nout) then
            call finish_status(nd_shape_mismatch, 'function result size changed', status, message)
            return
         end if
         if (.not. all_finite_complex(value)) then
            call finish_status(nd_nonfinite_value, &
               'complex function returned a nonfinite value', status, message)
            return
         end if
         jac(:, i) = aimag(value) / step
      end do
      call finish_status(nd_success, '', status, message)
   end subroutine jacobian_complex


   subroutine gend(func, x, result, options)
      procedure(vector_real_function) :: func
      real(dp), intent(in) :: x(:)
      type(gend_result), intent(out) :: result
      type(deriv_options), intent(in), optional :: options

      type(deriv_options) :: opts
      real(dp), allocatable :: f1(:), f2(:), d_approx(:, :), h_approx(:, :)
      real(dp), allocatable :: hdiag(:, :), xp(:), xm(:), h0(:)
      real(dp) :: h
      integer :: i, j, k, n, nout, u, stat
      character(:), allocatable :: msg

      opts = first_deriv_options()
      if (present(options)) opts = options
      call validate_options(opts, stat, msg)
      result%options = opts
      result%d = opts%d
      result%x = x
      result%p = size(x)
      result%status = stat
      result%message = msg
      if (stat /= nd_success) return

      result%f0 = func(x)
      nout = size(result%f0)
      n = size(x)
      if (.not. all(ieee_is_finite(result%f0))) then
         result%status = nd_nonfinite_value
         result%message = 'function returned a nonfinite value at x'
         return
      end if
      allocate(result%dmat(nout, n * (n + 3) / 2))
      allocate(d_approx(nout, opts%r), h_approx(nout, opts%r))
      allocate(hdiag(nout, n), xp(n), xm(n), h0(n))
      result%dmat = 0.0_dp
      hdiag = 0.0_dp
      h0 = initial_steps(x, opts)

      do i = 1, n
         h = h0(i)
         do k = 1, opts%r
            xp = x
            xm = x
            xp(i) = x(i) + h
            xm(i) = x(i) - h
            f1 = func(xp)
            f2 = func(xm)
            if (size(f1) /= nout .or. size(f2) /= nout) then
               result%status = nd_shape_mismatch
               result%message = 'function result size changed'
               return
            end if
            if (.not. all(ieee_is_finite(f1)) .or. .not. all(ieee_is_finite(f2))) then
               result%status = nd_nonfinite_value
               result%message = 'function returned a nonfinite value in genD'
               return
            end if
            d_approx(:, k) = (f1 - f2) / (2.0_dp * h)
            h_approx(:, k) = (f1 - 2.0_dp * result%f0 + f2) / (h * h)
            h = h / opts%v
         end do
         call extrapolate_columns(d_approx, opts)
         call extrapolate_columns(h_approx, opts)
         result%dmat(:, i) = d_approx(:, 1)
         hdiag(:, i) = h_approx(:, 1)
      end do

      u = n
      do i = 1, n
         do j = 1, i
            u = u + 1
            if (i == j) then
               result%dmat(:, u) = hdiag(:, i)
            else
               h = 1.0_dp
               do k = 1, opts%r
                  xp = x
                  xm = x
                  xp(i) = x(i) + h0(i) / opts%v ** (k - 1)
                  xp(j) = x(j) + h0(j) / opts%v ** (k - 1)
                  xm(i) = x(i) - h0(i) / opts%v ** (k - 1)
                  xm(j) = x(j) - h0(j) / opts%v ** (k - 1)
                  f1 = func(xp)
                  f2 = func(xm)
                  if (size(f1) /= nout .or. size(f2) /= nout) then
                     result%status = nd_shape_mismatch
                     result%message = 'function result size changed'
                     return
                  end if
                  if (.not. all(ieee_is_finite(f1)) .or. .not. all(ieee_is_finite(f2))) then
                     result%status = nd_nonfinite_value
                     result%message = 'function returned a nonfinite value in genD'
                     return
                  end if
                  h = h0(i) / opts%v ** (k - 1)
                  d_approx(:, k) = (f1 - 2.0_dp * result%f0 + f2 - &
                     hdiag(:, i) * h * h - &
                     hdiag(:, j) * (h0(j) / opts%v ** (k - 1)) ** 2) / &
                     (2.0_dp * h * (h0(j) / opts%v ** (k - 1)))
               end do
               call extrapolate_columns(d_approx, opts)
               result%dmat(:, u) = d_approx(:, 1)
            end if
         end do
      end do
      result%status = nd_success
      result%message = ''
   end subroutine gend


   subroutine hessian(func, x, hess, options, status, message)
      procedure(scalar_real_function) :: func
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: hess(:, :)
      type(deriv_options), intent(in), optional :: options
      integer, intent(out), optional :: status
      character(:), allocatable, intent(out), optional :: message

      type(deriv_options) :: opts
      real(dp), allocatable :: h0(:), xp(:), xm(:), diag_approx(:, :)
      real(dp), allocatable :: mixed_approx(:, :), hdiag(:)
      real(dp) :: f0, fp, fm, hi, hj
      integer :: i, j, k, n, stat
      character(:), allocatable :: msg

      call begin_status(status, message)
      opts = hessian_options()
      if (present(options)) opts = options
      call validate_options(opts, stat, msg)
      n = size(x)
      allocate(hess(n, n))
      hess = quiet_nan()
      if (stat /= nd_success) then
         call finish_status(stat, msg, status, message)
         return
      end if

      f0 = func(x)
      if (.not. ieee_is_finite(f0)) then
         call finish_status(nd_nonfinite_value, 'function returned a nonfinite value at x', &
            status, message)
         return
      end if

      allocate(h0(n), xp(n), xm(n), hdiag(n))
      allocate(diag_approx(1, opts%r), mixed_approx(1, opts%r))
      h0 = initial_steps(x, opts)

      do i = 1, n
         do k = 1, opts%r
            hi = h0(i) / opts%v ** (k - 1)
            xp = x
            xm = x
            xp(i) = x(i) + hi
            xm(i) = x(i) - hi
            fp = func(xp)
            fm = func(xm)
            if (.not. finite_pair(fp, fm)) then
               call finish_status(nd_nonfinite_value, &
                  'function returned a nonfinite value in Hessian calculation', &
                  status, message)
               return
            end if
            diag_approx(1, k) = (fp - 2.0_dp * f0 + fm) / (hi * hi)
         end do
         call extrapolate_columns(diag_approx, opts)
         hdiag(i) = diag_approx(1, 1)
         hess(i, i) = hdiag(i)
      end do

      do i = 1, n
         do j = 1, i - 1
            do k = 1, opts%r
               hi = h0(i) / opts%v ** (k - 1)
               hj = h0(j) / opts%v ** (k - 1)
               xp = x
               xm = x
               xp(i) = x(i) + hi
               xp(j) = x(j) + hj
               xm(i) = x(i) - hi
               xm(j) = x(j) - hj
               fp = func(xp)
               fm = func(xm)
               if (.not. finite_pair(fp, fm)) then
                  call finish_status(nd_nonfinite_value, &
                     'function returned a nonfinite value in Hessian calculation', &
                     status, message)
                  return
               end if
               mixed_approx(1, k) = (fp - 2.0_dp * f0 + fm - &
                  hdiag(i) * hi * hi - hdiag(j) * hj * hj) / (2.0_dp * hi * hj)
            end do
            call extrapolate_columns(mixed_approx, opts)
            hess(i, j) = mixed_approx(1, 1)
            hess(j, i) = hess(i, j)
         end do
      end do
      call finish_status(nd_success, '', status, message)
   end subroutine hessian


   subroutine hessian_complex(func, x, hess, options, status, message)
      procedure(scalar_complex_function) :: func
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: hess(:, :)
      type(deriv_options), intent(in), optional :: options
      integer, intent(out), optional :: status
      character(:), allocatable, intent(out), optional :: message

      type(deriv_options) :: opts
      real(dp), allocatable :: a(:, :, :), h(:), xp(:), xm(:), gp(:), gm(:)
      real(dp) :: factor
      integer :: i, k, m, n, stat
      character(:), allocatable :: msg

      call begin_status(status, message)
      opts = hessian_options()
      if (present(options)) opts = options
      call validate_options(opts, stat, msg)
      n = size(x)
      allocate(hess(n, n))
      hess = quiet_nan()
      if (stat /= nd_success) then
         call finish_status(stat, msg, status, message)
         return
      end if
      allocate(a(n, opts%r, n), h(n), xp(n), xm(n), gp(n), gm(n))
      h = initial_steps(x, opts)
      do k = 1, opts%r
         do i = 1, n
            xp = x
            xm = x
            xp(i) = x(i) + h(i)
            xm(i) = x(i) - h(i)
            call complex_gradient_at(func, xp, gp, stat)
            if (stat /= nd_success) then
               call finish_status(stat, 'complex function returned a nonfinite value', &
                  status, message)
               return
            end if
            call complex_gradient_at(func, xm, gm, stat)
            if (stat /= nd_success) then
               call finish_status(stat, 'complex function returned a nonfinite value', &
                  status, message)
               return
            end if
            a(:, k, i) = (gp - gm) / (2.0_dp * h(i))
         end do
         h = h / opts%v
      end do
      do m = 1, opts%r - 1
         factor = opts%v ** (2 * m)
         do k = 1, opts%r - m
            a(:, k, :) = (factor * a(:, k + 1, :) - a(:, k, :)) / (factor - 1.0_dp)
         end do
      end do
      hess = a(:, 1, :)
      hess = 0.5_dp * (hess + transpose(hess))
      call finish_status(nd_success, '', status, message)
   end subroutine hessian_complex


   subroutine complex_gradient_at(func, x, gradient, status)
      procedure(scalar_complex_function) :: func
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: gradient(:)
      integer, intent(out) :: status

      complex(dp), allocatable :: z(:)
      complex(dp) :: value
      real(dp), parameter :: step = epsilon(1.0_dp)
      integer :: i

      allocate(z(size(x)))
      z = cmplx(x, 0.0_dp, dp)
      do i = 1, size(x)
         z(i) = cmplx(x(i), step, dp)
         value = func(z)
         z(i) = cmplx(x(i), 0.0_dp, dp)
         if (.not. finite_complex(value)) then
            status = nd_nonfinite_value
            return
         end if
         gradient(i) = aimag(value) / step
      end do
      status = nd_success
   end subroutine complex_gradient_at


   subroutine extrapolate_columns(a, options)
      real(dp), intent(inout) :: a(:, :)
      type(deriv_options), intent(in) :: options
      real(dp) :: factor
      integer :: k, m

      do m = 1, options%r - 1
         factor = options%v ** (2 * m)
         do k = 1, options%r - m
            a(:, k) = (factor * a(:, k + 1) - a(:, k)) / (factor - 1.0_dp)
         end do
      end do
   end subroutine extrapolate_columns


   pure function initial_steps(x, options) result(h)
      real(dp), intent(in) :: x(:)
      type(deriv_options), intent(in) :: options
      real(dp) :: h(size(x))

      h = abs(options%d * x)
      where (abs(x) < options%zero_tol) h = options%eps
   end function initial_steps


   subroutine validate_options(options, status, message)
      type(deriv_options), intent(in) :: options
      integer, intent(out) :: status
      character(:), allocatable, intent(out) :: message

      if (options%eps <= 0.0_dp) then
         status = nd_invalid_argument
         message = 'eps must be positive'
      else if (options%d <= 0.0_dp) then
         status = nd_invalid_argument
         message = 'd must be positive'
      else if (options%zero_tol < 0.0_dp) then
         status = nd_invalid_argument
         message = 'zero_tol must be nonnegative'
      else if (options%r < 1) then
         status = nd_invalid_argument
         message = 'r must be at least 1'
      else if (options%v <= 1.0_dp) then
         status = nd_invalid_argument
         message = 'v must be greater than 1'
      else
         status = nd_success
         message = ''
      end if
   end subroutine validate_options


   subroutine make_side(n, supplied, direction, status, message)
      integer, intent(in) :: n
      integer, intent(in), optional :: supplied(:)
      integer, allocatable, intent(out) :: direction(:)
      integer, intent(out) :: status
      character(:), allocatable, intent(out) :: message

      allocate(direction(n))
      direction = 0
      if (present(supplied)) then
         if (size(supplied) /= n) then
            status = nd_shape_mismatch
            message = 'side must have size(x)'
            return
         end if
         if (any(abs(supplied) > 1)) then
            status = nd_invalid_argument
            message = 'side values must be -1, 0, or +1'
            return
         end if
         direction = supplied
      end if
      status = nd_success
      message = ''
   end subroutine make_side


   subroutine begin_status(status, message)
      integer, intent(out), optional :: status
      character(:), allocatable, intent(out), optional :: message
      if (present(status)) status = nd_success
      if (present(message)) message = ''
   end subroutine begin_status


   subroutine finish_status(code, text, status, message)
      integer, intent(in) :: code
      character(len=*), intent(in) :: text
      integer, intent(out), optional :: status
      character(:), allocatable, intent(out), optional :: message
      if (present(status)) status = code
      if (present(message)) message = text
   end subroutine finish_status


   pure real(dp) function quiet_nan()
      quiet_nan = ieee_value(0.0_dp, ieee_quiet_nan)
   end function quiet_nan


   pure logical function finite_pair(a, b)
      real(dp), intent(in) :: a, b
      finite_pair = ieee_is_finite(a) .and. ieee_is_finite(b)
   end function finite_pair


   pure logical function finite_complex(z)
      complex(dp), intent(in) :: z
      finite_complex = ieee_is_finite(real(z, dp)) .and. ieee_is_finite(aimag(z))
   end function finite_complex


   pure logical function all_finite_complex(z)
      complex(dp), intent(in) :: z(:)
      all_finite_complex = all(ieee_is_finite(real(z, dp))) .and. &
         all(ieee_is_finite(aimag(z)))
   end function all_finite_complex


   pure function lower(text) result(result)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: result
      integer :: i, code

      do i = 1, len(text)
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            result(i:i) = achar(code + iachar('a') - iachar('A'))
         else
            result(i:i) = text(i:i)
         end if
      end do
   end function lower

end module numderiv_core
