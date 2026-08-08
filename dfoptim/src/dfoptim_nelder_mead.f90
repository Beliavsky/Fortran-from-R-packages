! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of computational code from the R package dfoptim.
! Original authors: Ravi Varadhan, Hans W. Borchers, and Vincent Bechard.

module dfoptim_nelder_mead
   use dfoptim_kinds, only : dp
   use dfoptim_interfaces, only : nmk_control_t, dfoptim_result_t, dfoptim_objective, &
      dfoptim_monitor, dfoptim_success, dfoptim_max_evaluations, dfoptim_stagnation, &
      dfoptim_invalid_input, dfoptim_nonfinite_objective, dfoptim_cancelled
   use dfoptim_utils, only : evaluate_objective, call_monitor, sort_simplex, status_message
   implicit none
   private
   public :: nmk, nmkb

contains

   function nmk(par, fn, control, user_data, monitor) result(result)
      real(dp), intent(in) :: par(:)
      procedure(dfoptim_objective) :: fn
      type(nmk_control_t), intent(in), optional :: control
      class(*), intent(inout), optional :: user_data
      procedure(dfoptim_monitor), optional :: monitor
      type(dfoptim_result_t) :: result
      real(dp), allocatable :: lower(:), upper(:)
      integer, allocatable :: transform(:)

      allocate(lower(size(par)), upper(size(par)), transform(size(par)))
      lower = -huge(1.0_dp)
      upper = huge(1.0_dp)
      transform = 0
      result = nmk_impl(par, fn, lower, upper, transform, .false., control, user_data, monitor)
   end function nmk

   function nmkb(par, fn, lower, upper, control, user_data, monitor) result(result)
      real(dp), intent(in) :: par(:), lower(:), upper(:)
      procedure(dfoptim_objective) :: fn
      type(nmk_control_t), intent(in), optional :: control
      class(*), intent(inout), optional :: user_data
      procedure(dfoptim_monitor), optional :: monitor
      type(dfoptim_result_t) :: result
      integer, allocatable :: transform(:)
      integer :: i, n
      logical :: lo_finite, up_finite

      n = size(par)
      allocate(transform(n))
      transform = 0
      if (size(lower) == n .and. size(upper) == n) then
         do i = 1, n
            lo_finite = abs(lower(i)) < huge(1.0_dp) / 2.0_dp
            up_finite = abs(upper(i)) < huge(1.0_dp) / 2.0_dp
            if (lo_finite .and. up_finite) then
               transform(i) = 1
            else if (lo_finite) then
               transform(i) = 3
            else if (up_finite) then
               transform(i) = 4
            else
               transform(i) = 2
            end if
         end do
      end if
      result = nmk_impl(par, fn, lower, upper, transform, .true., control, user_data, monitor)
   end function nmkb

   function nmk_impl(par, fn, lower, upper, transform, bounded, control_in, user_data, monitor) result(result)
      real(dp), intent(in) :: par(:), lower(:), upper(:)
      integer, intent(in) :: transform(:)
      logical, intent(in) :: bounded
      procedure(dfoptim_objective) :: fn
      type(nmk_control_t), intent(in), optional :: control_in
      class(*), intent(inout), optional :: user_data
      procedure(dfoptim_monitor), optional :: monitor
      type(dfoptim_result_t) :: result
      type(nmk_control_t) :: control
      real(dp), allocatable :: v(:, :), f(:), x0(:), x_external(:), xbar(:), xr(:), xe(:), xc(:)
      real(dp), allocatable :: differences(:, :), delf(:), diameter(:), sgrad(:), xnew(:), sx(:)
      real(dp) :: scale, alpha1, alpha2, rho, gamma, chi, sigma
      real(dp) :: distance, simplex_size, fbc, fr, fe, fc, fnew, fbt, delfb, armtst
      real(dp) :: alpha_restart, diams, sgrad_norm
      integer :: n, nf, iteration, restarts, maxfeval, i, j, happy
      logical :: finite_ok, cancelled, valid

      control = nmk_control_t()
      if (present(control_in)) control = control_in
      n = size(par)
      allocate(result%x(n))
      result%x = par

      valid = n >= 2 .and. size(lower) == n .and. size(upper) == n .and. &
         size(transform) == n .and. control%tol > 0.0_dp .and. control%max_restarts >= 0
      if (bounded) then
         valid = valid .and. any(transform /= 2) .and. all(lower <= par) .and. all(par <= upper)
         do i = 1, n
            select case (transform(i))
            case (1)
               valid = valid .and. lower(i) < par(i) .and. par(i) < upper(i)
            case (3)
               valid = valid .and. lower(i) < par(i)
            case (4)
               valid = valid .and. par(i) < upper(i)
            end select
         end do
      end if
      if (.not. valid) then
         result%convergence = dfoptim_invalid_input
         result%message = status_message(result%convergence)
         return
      end if

      maxfeval = control%maxfeval
      if (maxfeval <= 0) maxfeval = min(5000, max(1500, 20 * n * n))
      if (maxfeval < n + 1) then
         result%convergence = dfoptim_invalid_input
         result%message = status_message(result%convergence)
         return
      end if

      allocate(v(n, n + 1), f(n + 1), x0(n), x_external(n), xbar(n), xr(n), xe(n), xc(n))
      allocate(differences(n, n), delf(n), diameter(n), sgrad(n), xnew(n), sx(n))
      call to_internal(par, lower, upper, transform, x0)
      v = 0.0_dp
      v(:, 1) = x0
      scale = max(1.0_dp, sqrt(sum(x0 * x0)))

      if (control%regular_simplex) then
         alpha1 = scale / (real(n, dp) * sqrt(2.0_dp)) * (sqrt(real(n + 1, dp)) + real(n - 1, dp))
         alpha2 = scale / (real(n, dp) * sqrt(2.0_dp)) * (sqrt(real(n + 1, dp)) - 1.0_dp)
         do j = 2, n + 1
            v(:, j) = x0 + alpha2
            v(j - 1, j) = x0(j - 1) + alpha1
         end do
      else
         do j = 2, n + 1
            v(:, j) = x0
            v(j - 1, j) = x0(j - 1) + scale
         end do
      end if

      nf = 0
      do j = 1, n + 1
         call evaluate_internal(v(:, j), f(j), finite_ok)
         nf = nf + 1
         if (j == 1 .and. .not. finite_ok) then
            result%convergence = dfoptim_nonfinite_objective
            result%feval = nf
            result%message = status_message(result%convergence)
            return
         end if
      end do
      call sort_simplex(v, f)

      rho = 1.0_dp
      gamma = 0.5_dp
      chi = 2.0_dp
      sigma = 0.5_dp
      restarts = 0
      iteration = 0
      cancelled = .false.
      call update_geometry()

      do while (nf < maxfeval .and. restarts < control%max_restarts .and. &
                distance > control%tol .and. simplex_size > 1.0e-6_dp)
         fbc = sum(f) / real(n + 1, dp)
         happy = 0
         iteration = iteration + 1
         xbar = sum(v(:, 1:n), dim=2) / real(n, dp)
         xr = (1.0_dp + rho) * xbar - rho * v(:, n + 1)
         call evaluate_internal(xr, fr, finite_ok)
         nf = nf + 1

         if (fr >= f(1) .and. fr < f(n)) then
            happy = 1; xnew = xr; fnew = fr
         else if (fr < f(1)) then
            xe = (1.0_dp + rho * chi) * xbar - rho * chi * v(:, n + 1)
            if (nf < maxfeval) then
               call evaluate_internal(xe, fe, finite_ok)
               nf = nf + 1
            else
               fe = huge(1.0_dp)
            end if
            if (fe < fr) then
               xnew = xe; fnew = fe
            else
               xnew = xr; fnew = fr
            end if
            happy = 1
         else if (fr >= f(n) .and. fr < f(n + 1)) then
            xc = (1.0_dp + rho * gamma) * xbar - rho * gamma * v(:, n + 1)
            if (nf < maxfeval) then
               call evaluate_internal(xc, fc, finite_ok)
               nf = nf + 1
            else
               fc = huge(1.0_dp)
            end if
            if (fc <= fr) then
               xnew = xc; fnew = fc; happy = 1
            end if
         else
            xc = (1.0_dp - gamma) * xbar + gamma * v(:, n + 1)
            if (nf < maxfeval) then
               call evaluate_internal(xc, fc, finite_ok)
               nf = nf + 1
            else
               fc = huge(1.0_dp)
            end if
            if (fc < f(n + 1)) then
               xnew = xc; fnew = fc; happy = 1
            end if
         end if

         if (happy == 1) then
            fbt = (sum(f(1:n)) + fnew) / real(n + 1, dp)
            delfb = fbt - fbc
            armtst = alpha_restart * sum(sgrad * sgrad)
            if (delfb > -armtst / real(n, dp)) then
               restarts = restarts + 1
               diams = minval(diameter)
               sx = sign(1.0_dp, sgrad)
               where (abs(sgrad) <= tiny(1.0_dp)) sx = 0.0_dp
               do j = 2, n + 1
                  v(:, j) = v(:, 1)
               end do
               do j = 1, n
                  v(j, j + 1) = v(j, j + 1) - diams * sx(j)
               end do
               happy = 0
               if (control%trace) write(*, '(a)') 'Trouble - restarting'
            end if
         end if

         if (happy == 1) then
            v(:, n + 1) = xnew
            f(n + 1) = fnew
            call sort_simplex(v, f)
         else if (restarts < control%max_restarts) then
            do j = 2, n + 1
               v(:, j) = v(:, 1) - sigma * (v(:, j) - v(:, 1))
               if (nf < maxfeval) then
                  call evaluate_internal(v(:, j), f(j), finite_ok)
                  nf = nf + 1
               else
                  f(j) = huge(1.0_dp)
               end if
            end do
            call sort_simplex(v, f)
         end if

         call update_geometry()
         if (control%trace .and. mod(iteration, 2) == 0) then
            write(*, '(a,i0,2x,a,es22.13)') 'iter: ', iteration, 'value: ', &
               merge(-f(1), f(1), control%maximize)
         end if
         call from_internal(v(:, 1), lower, upper, transform, x_external)
         call call_monitor(monitor, x_external, merge(-f(1), f(1), control%maximize), &
            iteration, nf, cancelled, user_data)
         if (cancelled) exit
      end do

      call from_internal(v(:, 1), lower, upper, transform, result%x)
      result%value = merge(-f(1), f(1), control%maximize)
      result%feval = nf
      result%niter = iteration
      result%restarts = restarts
      if (cancelled) then
         result%convergence = dfoptim_cancelled
      else if (distance <= control%tol .or. simplex_size <= 1.0e-6_dp) then
         result%convergence = dfoptim_success
      else if (nf >= maxfeval) then
         result%convergence = dfoptim_max_evaluations
      else if (restarts >= control%max_restarts) then
         result%convergence = dfoptim_stagnation
      else
         result%convergence = dfoptim_stagnation
      end if
      result%message = status_message(result%convergence)

   contains

      subroutine evaluate_internal(z, value, ok)
         real(dp), intent(in) :: z(:)
         real(dp), intent(out) :: value
         logical, intent(out) :: ok
         real(dp) :: original(size(z))
         call from_internal(z, lower, upper, transform, original)
         value = evaluate_objective(fn, original, control%maximize, user_data, ok)
      end subroutine evaluate_internal

      subroutine update_geometry()
         integer :: k
         do k = 1, n
            differences(:, k) = v(:, k + 1) - v(:, 1)
            delf(k) = f(k + 1) - f(1)
            diameter(k) = sqrt(sum(differences(:, k) ** 2))
         end do
         simplex_size = sum(abs(differences)) / max(1.0_dp, sum(abs(v(:, 1))))
         distance = f(n + 1) - f(1)
         sgrad = matmul(differences, delf)
         sgrad_norm = sqrt(sum(sgrad * sgrad))
         if (sgrad_norm > tiny(1.0_dp)) then
            alpha_restart = 1.0e-4_dp * maxval(diameter) / sgrad_norm
         else
            alpha_restart = 0.0_dp
         end if
      end subroutine update_geometry

   end function nmk_impl

   pure subroutine to_internal(x, lower, upper, transform, z)
      real(dp), intent(in) :: x(:), lower(:), upper(:)
      integer, intent(in) :: transform(:)
      real(dp), intent(out) :: z(:)
      integer :: i
      do i = 1, size(x)
         select case (transform(i))
         case (1)
            z(i) = atanh(2.0_dp * (x(i) - lower(i)) / (upper(i) - lower(i)) - 1.0_dp)
         case (3)
            z(i) = log(x(i) - lower(i))
         case (4)
            z(i) = log(upper(i) - x(i))
         case default
            z(i) = x(i)
         end select
      end do
   end subroutine to_internal

   pure subroutine from_internal(z, lower, upper, transform, x)
      real(dp), intent(in) :: z(:), lower(:), upper(:)
      integer, intent(in) :: transform(:)
      real(dp), intent(out) :: x(:)
      integer :: i
      real(dp) :: ez
      do i = 1, size(z)
         select case (transform(i))
         case (1)
            x(i) = lower(i) + 0.5_dp * (upper(i) - lower(i)) * (1.0_dp + tanh(z(i)))
         case (3)
            ez = exp(min(z(i), log(huge(1.0_dp)) - 2.0_dp))
            x(i) = lower(i) + ez
         case (4)
            ez = exp(min(z(i), log(huge(1.0_dp)) - 2.0_dp))
            x(i) = upper(i) - ez
         case default
            x(i) = z(i)
         end select
      end do
   end subroutine from_internal

end module dfoptim_nelder_mead
