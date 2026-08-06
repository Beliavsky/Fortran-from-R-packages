! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran computational translation of tvgarch 2.4.3.
module tvgarch_model
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use garchx_kinds, only : dp
   use garchx_math, only : mean_value, empirical_quantile
   use garchx_linalg, only : invert_matrix
   use garchx_optimize, only : bounded_nelder_mead, numerical_hessian
   use garchx_model, only : garchx_spec, garchx_fit, make_garchx_spec, fit_garchx, &
                            garchx_filter, garchx_simulate, garchx_forecast
   use tvgarch_transition, only : tv_spec, make_tv_spec, tv_component, tv_parameter_count, &
                                  pack_tv_parameters, unpack_tv_parameters, combinations_binary
   implicit none
   private
   real(dp), parameter :: log_two_pi = log(2.0_dp*acos(-1.0_dp))

   type, public :: tvgarch_spec
      type(tv_spec) :: tv
      type(garchx_spec) :: garch
   end type tvgarch_spec

   type, public :: tvgarch_fit
      type(tvgarch_spec) :: spec
      type(garchx_fit) :: hfit
      real(dp), allocatable :: y(:), xtv(:), xreg(:, :)
      real(dp), allocatable :: par_g(:), se_g(:), vcov_g(:, :)
      real(dp), allocatable :: g(:), h(:), sigma2(:), residuals(:)
      real(dp) :: loglik = -huge(1.0_dp)
      integer :: status = 1
      integer :: iterations = 0
      logical :: joint_refined = .false.
   end type tvgarch_fit

   type, public :: tvgarch_simulation
      real(dp), allocatable :: y(:), sigma2(:), g(:), h(:), innovations(:)
      integer :: status = 0
   end type tvgarch_simulation

   public :: make_tvgarch_spec, make_garch_order_spec, tv_objective_value
   public :: tvgarch_objective_value, fit_tvgarch, tvgarch_simulate
   public :: tvgarch_forecast, tvgarch_quantile_path, fitted_tvgarch
contains
   subroutine make_garch_order_spec(order_h, xreg_count, spec, status)
      integer, intent(in) :: order_h(3)
      integer, intent(in), optional :: xreg_count
      type(garchx_spec), intent(out) :: spec
      integer, intent(out) :: status
      integer :: nx, i
      integer, allocatable :: arch(:), garch(:), asym(:)
      status = 0
      if (any(order_h < 0)) then
         status = 1
         call make_garchx_spec(spec)
         return
      end if
      allocate(garch(order_h(1)), arch(order_h(2)), asym(order_h(3)))
      do i = 1, order_h(1); garch(i) = i; end do
      do i = 1, order_h(2); arch(i) = i; end do
      do i = 1, order_h(3); asym(i) = i; end do
      nx = 0
      if (present(xreg_count)) nx = xreg_count
      call make_garchx_spec(spec, arch_lags=arch, garch_lags=garch, &
                            asym_lags=asym, xreg_count=nx)
   end subroutine make_garch_order_spec

   subroutine make_tvgarch_spec(spec, order_g, order_h, speed_option, xreg_count, status)
      type(tvgarch_spec), intent(out) :: spec
      integer, intent(in), optional :: order_g(:), order_h(3), speed_option, xreg_count
      integer, intent(out) :: status
      integer :: option, nx, horder(3)
      status = 0
      option = 2
      nx = 0
      horder = [1, 1, 0]
      if (present(speed_option)) option = speed_option
      if (present(xreg_count)) nx = xreg_count
      if (present(order_h)) horder = order_h
      if (present(order_g)) then
         call make_tv_spec(spec%tv, orders=order_g, speed_option=option)
      else
         call make_tv_spec(spec%tv)
      end if
      call make_garch_order_spec(horder, nx, spec%garch, status)
   end subroutine make_tvgarch_spec

   pure real(dp) function normal_negloglik(y, variance) result(value)
      real(dp), intent(in) :: y(:), variance(:)
      if (size(y) /= size(variance) .or. any(variance <= 0.0_dp) .or. &
          .not. all(ieee_is_finite(variance))) then
         value = huge(1.0_dp)*0.01_dp
      else
         value = 0.5_dp*sum(log_two_pi + log(variance) + y*y/variance)
      end if
   end function normal_negloglik

   logical function valid_tv_parameters(orders, par, xtv, speed_option) result(ok)
      integer, intent(in) :: orders(:), speed_option
      real(dp), intent(in) :: par(:), xtv(:)
      type(tv_spec) :: spec
      integer :: status, s, i, first, last
      integer, allocatable :: binary(:, :)
      real(dp), allocatable :: g(:)
      ok = .false.
      if (size(par) /= tv_parameter_count(orders)) return
      call unpack_tv_parameters(orders, par, speed_option, spec, status)
      if (status /= 0 .or. spec%intercept <= 1.0e-10_dp) return
      if (size(orders) > 0) then
         call combinations_binary(size(orders), binary, status)
         if (status /= 0) return
         do i = 1, size(binary, 1)
            if (spec%intercept + dot_product(real(binary(i, :), dp), spec%sizes) <= 1.0e-10_dp) return
         end do
         first = 1
         do s = 1, size(orders)
            last = first+orders(s)-1
            if (any(spec%locations(first:last) < minval(xtv)+1.0e-8_dp) .or. &
                any(spec%locations(first:last) > maxval(xtv)-1.0e-8_dp)) return
            if (orders(s) > 1) then
               do i = first, last-1
                  if (spec%locations(i+1) < spec%locations(i)) return
               end do
            end if
            first = last+1
         end do
      end if
      call tv_component(spec, xtv, g, status)
      ok = status == 0 .and. all(g > 1.0e-10_dp)
   end function valid_tv_parameters

   subroutine tv_objective_value(par_g, orders, xtv, fixed_h, y, speed_option, value, status, per_observation)
      real(dp), intent(in) :: par_g(:), xtv(:), fixed_h(:), y(:)
      integer, intent(in) :: orders(:), speed_option
      real(dp), intent(out) :: value
      integer, intent(out) :: status
      real(dp), allocatable, intent(out), optional :: per_observation(:)
      type(tv_spec) :: spec
      real(dp), allocatable :: g(:), phi(:)

      status = 0
      value = huge(1.0_dp)*0.01_dp
      if (size(xtv) /= size(y) .or. size(fixed_h) /= size(y) .or. any(fixed_h <= 0.0_dp)) then
         status = 1
         if (present(per_observation)) allocate(per_observation(0))
         return
      end if
      if (.not. valid_tv_parameters(orders, par_g, xtv, speed_option)) then
         status = 2
         if (present(per_observation)) allocate(per_observation(0))
         return
      end if
      call unpack_tv_parameters(orders, par_g, speed_option, spec, status)
      call tv_component(spec, xtv, g, status)
      if (status /= 0) return
      allocate(phi(size(y)))
      phi = y/sqrt(fixed_h)
      value = normal_negloglik(phi, g)
      if (present(per_observation)) then
         allocate(per_observation(size(y)))
         per_observation = 0.5_dp*(log_two_pi+log(g)+phi*phi/g)
      end if
   end subroutine tv_objective_value

   subroutine tvgarch_objective_value(par, fixed_intercept_g, spec, xtv, y, value, status, xreg)
      real(dp), intent(in) :: par(:), fixed_intercept_g, xtv(:), y(:)
      type(tvgarch_spec), intent(in) :: spec
      real(dp), intent(out) :: value
      integer, intent(out) :: status
      real(dp), intent(in), optional :: xreg(:, :)
      integer :: ng_free, nh, hstatus
      real(dp), allocatable :: pg(:), ph(:), g(:), h(:)
      type(tv_spec) :: local_tv

      status = 0
      ng_free = tv_parameter_count(spec%tv%orders)-1
      nh = 1 + size(spec%garch%arch_lags) + size(spec%garch%garch_lags) + &
           size(spec%garch%asym_lags) + spec%garch%xreg_count
      if (size(par) /= ng_free+nh) then
         status = 1; value = huge(1.0_dp)*0.01_dp; return
      end if
      allocate(pg(ng_free+1), ph(nh))
      pg(1) = fixed_intercept_g
      pg(2:) = par(1:ng_free)
      ph = par(ng_free+1:)
      if (.not. valid_tv_parameters(spec%tv%orders, pg, xtv, spec%tv%speed_option)) then
         status = 2; value = huge(1.0_dp)*0.01_dp; return
      end if
      call unpack_tv_parameters(spec%tv%orders, pg, spec%tv%speed_option, local_tv, status)
      call tv_component(local_tv, xtv, g, status)
      if (status /= 0) then
         value = huge(1.0_dp)*0.01_dp; return
      end if
      call garchx_filter(y/sqrt(g), spec%garch, ph, h, hstatus, xreg)
      if (hstatus /= 0 .or. any(h <= 0.0_dp)) then
         status = 3; value = huge(1.0_dp)*0.01_dp; return
      end if
      value = normal_negloglik(y, g*h)
   end subroutine tvgarch_objective_value

   subroutine initialize_tv(spec, xtv, par)
      type(tvgarch_spec), intent(in) :: spec
      real(dp), intent(in) :: xtv(:)
      real(dp), allocatable, intent(out) :: par(:)
      integer :: s, i, j, first, last
      real(dp) :: xmin, xmax, xmean, xsd
      s = size(spec%tv%orders)
      allocate(par(tv_parameter_count(spec%tv%orders)))
      par = 0.1_dp
      par(1) = 1.0_dp
      if (s == 0) return
      par(2+s:1+2*s) = merge(log(10.0_dp), 10.0_dp, spec%tv%speed_option == 2)
      xmin = minval(xtv); xmax = maxval(xtv); xmean = sum(xtv)/real(size(xtv), dp)
      if (size(xtv) > 1) then
         xsd = sqrt(sum((xtv-xmean)**2)/real(size(xtv)-1, dp))
      else
         xsd = 0.0_dp
      end if
      first = 2+2*s
      do i = 1, s
         last = first+spec%tv%orders(i)-1
         if (spec%tv%orders(i) == 1) then
            par(first) = xmean
         else
            do j = 1, spec%tv%orders(i)
               par(first+j-1) = xmin+0.5_dp*xsd + real(j-1, dp)* &
                  max(xmax-xmin-xsd, epsilon(1.0_dp))/real(spec%tv%orders(i)-1, dp)
            end do
         end if
         first = last+1
      end do
   end subroutine initialize_tv

   subroutine tv_bounds(spec, xtv, y, lower, upper)
      type(tvgarch_spec), intent(in) :: spec
      real(dp), intent(in) :: xtv(:), y(:)
      real(dp), allocatable, intent(out) :: lower(:), upper(:)
      integer :: s
      real(dp) :: scale
      s = size(spec%tv%orders)
      allocate(lower(tv_parameter_count(spec%tv%orders)), upper(tv_parameter_count(spec%tv%orders)))
      scale = max(mean_value(y*y), 1.0e-3_dp)
      lower = -20.0_dp*scale
      upper = 20.0_dp*scale
      lower(1) = 1.0e-6_dp
      upper(1) = 20.0_dp*scale
      if (s > 0) then
         if (spec%tv%speed_option == 2) then
            lower(2+s:1+2*s) = log(1.0e-4_dp)
            upper(2+s:1+2*s) = log(250.0_dp)
         else
            lower(2+s:1+2*s) = 1.0e-5_dp
            upper(2+s:1+2*s) = 250.0_dp
         end if
         lower(2+2*s:) = minval(xtv)+1.0e-6_dp
         upper(2+2*s:) = maxval(xtv)-1.0e-6_dp
      end if
   end subroutine tv_bounds

   subroutine fit_tvgarch(y, spec, fit, xtv, xreg, initial_g, initial_h, max_outer, &
                          rel_tol, joint_refine, turbo)
      real(dp), intent(in) :: y(:)
      type(tvgarch_spec), intent(in) :: spec
      type(tvgarch_fit), intent(out) :: fit
      real(dp), intent(in), optional :: xtv(:), xreg(:, :), initial_g(:), initial_h(:), rel_tol
      integer, intent(in), optional :: max_outer
      logical, intent(in), optional :: joint_refine, turbo
      integer :: n, s, ng, nh, status, it, maxit, opt_status, opt_iter, inv_status, i
      real(dp) :: tol, fbest, conv_g, conv_h, fixed_intercept
      logical :: do_joint, fast
      real(dp), allocatable :: xvar(:), h(:), g(:), phi(:), pg(:), pg_new(:), pg_free(:)
      real(dp), allocatable :: lo(:), hi(:), lo_free(:), hi_free(:), old_h(:), old_pg(:)
      real(dp), allocatable :: hess(:, :), hinv(:, :), scores(:, :), meat(:, :), per0(:), perp(:), perm(:)
      real(dp), allocatable :: pjoint(:), ljoint(:), ujoint(:), pjoint_best(:)
      type(tv_spec) :: final_tv

      fit%status = 0
      n = size(y)
      if (n < 20 .or. .not. all(ieee_is_finite(y))) then
         fit%status = 1; return
      end if
      allocate(fit%y(n)); fit%y = y
      allocate(xvar(n))
      if (present(xtv)) then
         if (size(xtv) /= n) then; fit%status = 2; return; end if
         xvar = xtv
      else
         do i = 1, n; xvar(i) = real(i, dp)/real(n, dp); end do
      end if
      allocate(fit%xtv(n)); fit%xtv = xvar
      if (spec%garch%xreg_count > 0) then
         if (.not. present(xreg)) then; fit%status = 3; return; end if
         if (size(xreg,1) /= n .or. size(xreg,2) /= spec%garch%xreg_count) then
            fit%status = 4; return
         end if
         allocate(fit%xreg(n, spec%garch%xreg_count)); fit%xreg = xreg
      else
         allocate(fit%xreg(n,0))
      end if
      fit%spec = spec
      s = size(spec%tv%orders)
      ng = tv_parameter_count(spec%tv%orders)
      nh = 1+size(spec%garch%arch_lags)+size(spec%garch%garch_lags)+ &
           size(spec%garch%asym_lags)+spec%garch%xreg_count
      if (s == 0) then
         if (present(initial_h)) then
            call fit_garchx(y, spec%garch, fit%hfit, xreg=xreg, initial=initial_h, vcov_type='robust')
         else
            call fit_garchx(y, spec%garch, fit%hfit, xreg=xreg, vcov_type='robust')
         end if
         if (.not. allocated(fit%hfit%sigma2)) then; fit%status = 5; return; end if
         allocate(fit%par_g(0), fit%se_g(0), fit%vcov_g(0,0), fit%g(n), fit%h(n), &
                  fit%sigma2(n), fit%residuals(n))
         fit%g = 1.0_dp; fit%h = fit%hfit%sigma2; fit%sigma2 = fit%h
         fit%residuals = y/sqrt(fit%sigma2); fit%loglik = -normal_negloglik(y, fit%sigma2)
         fit%status = fit%hfit%status
         return
      end if
      call initialize_tv(spec, xvar, pg)
      if (present(initial_g)) then
         if (size(initial_g) /= ng) then; fit%status = 6; return; end if
         pg = initial_g
      end if
      call tv_bounds(spec, xvar, y, lo, hi)
      allocate(h(n)); h = 1.0_dp
      call bounded_nelder_mead(tv_full_obj, pg, lo, hi, pg_new, fbest, opt_status, opt_iter, 2000, 1.0e-7_dp)
      pg = pg_new
      fixed_intercept = pg(1)
      allocate(pg_free(ng-1), lo_free(ng-1), hi_free(ng-1))
      pg_free = pg(2:); lo_free = lo(2:); hi_free = hi(2:)
      maxit = 50; if (present(max_outer)) maxit = max_outer
      tol = 1.0e-5_dp; if (present(rel_tol)) tol = rel_tol
      fast = .false.; if (present(turbo)) fast = turbo
      allocate(old_pg(ng-1), old_h(n), phi(n))
      do it = 1, maxit
         call unpack_tv_parameters(spec%tv%orders, [fixed_intercept, pg_free], spec%tv%speed_option, final_tv, status)
         call tv_component(final_tv, xvar, g, status)
         if (status /= 0) then; fit%status = 7; return; end if
         phi = y/sqrt(g)
         old_h = h
         if (it == 1 .and. present(initial_h)) then
            call fit_garchx(phi, spec%garch, fit%hfit, xreg=xreg, initial=initial_h, vcov_type='robust')
         else if (it > 1 .and. allocated(fit%hfit%par)) then
            call fit_garchx(phi, spec%garch, fit%hfit, xreg=xreg, initial=fit%hfit%par, vcov_type='robust')
         else
            call fit_garchx(phi, spec%garch, fit%hfit, xreg=xreg, vcov_type='robust')
         end if
         if (.not. allocated(fit%hfit%sigma2)) then; fit%status = 8; return; end if
         h = fit%hfit%sigma2
         old_pg = pg_free
         call bounded_nelder_mead(tv_free_obj, pg_free, lo_free, hi_free, pg_new, fbest, &
                                  opt_status, opt_iter, 2000, 1.0e-7_dp)
         pg_free = pg_new
         conv_g = maxval(abs(pg_free-old_pg))
         conv_h = maxval(abs(h-old_h))/max(1.0_dp, maxval(abs(old_h)))
         if (max(conv_g, conv_h) < tol .and. it >= 2) exit
      end do
      fit%iterations = min(it, maxit)
      allocate(fit%par_g(ng)); fit%par_g = [fixed_intercept, pg_free]
      call unpack_tv_parameters(spec%tv%orders, fit%par_g, spec%tv%speed_option, final_tv, status)
      fit%spec%tv = final_tv
      call tv_component(final_tv, xvar, fit%g, status)
      phi = y/sqrt(fit%g)
      call fit_garchx(phi, spec%garch, fit%hfit, xreg=xreg, initial=fit%hfit%par, vcov_type='robust')
      allocate(fit%h(n), fit%sigma2(n), fit%residuals(n))
      fit%h = fit%hfit%sigma2
      do_joint = .false.; if (present(joint_refine)) do_joint = joint_refine
      if (do_joint) then
         allocate(pjoint(ng-1+nh), ljoint(ng-1+nh), ujoint(ng-1+nh))
         pjoint = [fit%par_g(2:), fit%hfit%par]
         ljoint(1:ng-1) = lo(2:); ujoint(1:ng-1) = hi(2:)
         ljoint(ng:) = 0.0_dp; ujoint(ng:) = huge(1.0_dp)**0.20_dp
         call bounded_nelder_mead(joint_obj, pjoint, ljoint, ujoint, pjoint_best, fbest, &
                                  opt_status, opt_iter, 3000, 1.0e-7_dp)
         fit%par_g(2:) = pjoint_best(1:ng-1)
         fit%hfit%par = pjoint_best(ng:)
         call unpack_tv_parameters(spec%tv%orders, fit%par_g, spec%tv%speed_option, final_tv, status)
         fit%spec%tv = final_tv
         call tv_component(final_tv, xvar, fit%g, status)
         call garchx_filter(y/sqrt(fit%g), spec%garch, fit%hfit%par, fit%h, status, xreg)
         fit%joint_refined = opt_status == 0
      end if
      fit%sigma2 = fit%g*fit%h
      fit%residuals = y/sqrt(fit%sigma2)
      fit%loglik = -normal_negloglik(y, fit%sigma2)
      allocate(fit%se_g(ng)); fit%se_g = 0.0_dp
      allocate(fit%vcov_g(max(ng-1,0), max(ng-1,0)))
      if (.not. fast .and. ng > 1) then
         call numerical_hessian(tv_free_obj, fit%par_g(2:), hess)
         call invert_matrix(hess, hinv, inv_status)
         if (inv_status == 0) then
            allocate(scores(n,ng-1), meat(ng-1,ng-1))
            call tv_perobs(fit%par_g(2:), per0)
            do i = 1, ng-1
               pg_free = fit%par_g(2:)
               pg_free(i) = pg_free(i)+epsilon(1.0_dp)**0.25_dp*max(1.0_dp,abs(pg_free(i)))
               call tv_perobs(pg_free, perp)
               pg_free = fit%par_g(2:)
               pg_free(i) = pg_free(i)-epsilon(1.0_dp)**0.25_dp*max(1.0_dp,abs(pg_free(i)))
               call tv_perobs(pg_free, perm)
               scores(:,i) = (perp-perm)/(2.0_dp*epsilon(1.0_dp)**0.25_dp* &
                              max(1.0_dp,abs(fit%par_g(i+1))))
            end do
            meat = matmul(transpose(scores), scores)
            fit%vcov_g = matmul(hinv, matmul(meat, hinv))
            fit%se_g(2:) = sqrt(max([(fit%vcov_g(i,i), i=1,ng-1)], 0.0_dp))
         else
            fit%vcov_g = 0.0_dp
         end if
      else
         fit%vcov_g = 0.0_dp
      end if
      fit%status = 0
   contains
      function tv_full_obj(p) result(value)
         real(dp), intent(in) :: p(:)
         real(dp) :: value
         integer :: st
         call tv_objective_value(p, spec%tv%orders, xvar, h, y, spec%tv%speed_option, value, st)
         if (st /= 0) value = huge(1.0_dp)*0.01_dp
      end function tv_full_obj
      function tv_free_obj(p) result(value)
         real(dp), intent(in) :: p(:)
         real(dp) :: value
         integer :: st
         call tv_objective_value([fixed_intercept,p], spec%tv%orders, xvar, h, y, &
                                 spec%tv%speed_option, value, st)
         if (st /= 0) value = huge(1.0_dp)*0.01_dp
      end function tv_free_obj
      function joint_obj(p) result(value)
         real(dp), intent(in) :: p(:)
         real(dp) :: value
         integer :: st
         call tvgarch_objective_value(p, fixed_intercept, spec, xvar, y, value, st, xreg)
         if (st /= 0) value = huge(1.0_dp)*0.01_dp
      end function joint_obj
      subroutine tv_perobs(p, values)
         real(dp), intent(in) :: p(:)
         real(dp), allocatable, intent(out) :: values(:)
         real(dp) :: val
         integer :: st
         call tv_objective_value([fixed_intercept,p], spec%tv%orders, xvar, h, y, &
                                 spec%tv%speed_option, val, st, values)
         if (st /= 0) then
            if (allocated(values)) deallocate(values)
            allocate(values(n)); values = huge(1.0_dp)*0.01_dp
         end if
      end subroutine tv_perobs
   end subroutine fit_tvgarch

   subroutine tvgarch_simulate(n, spec, par_h, result, xtv, xreg, innovations)
      integer, intent(in) :: n
      type(tvgarch_spec), intent(in) :: spec
      real(dp), intent(in) :: par_h(:)
      type(tvgarch_simulation), intent(out) :: result
      real(dp), intent(in), optional :: xtv(:), xreg(:, :), innovations(:)
      real(dp), allocatable :: xv(:), yh(:)
      integer :: i, status
      allocate(xv(n))
      if (present(xtv)) then
         if (size(xtv) /= n) then; result%status=1; return; end if
         xv = xtv
      else
         do i=1,n; xv(i)=real(i,dp)/real(n,dp); end do
      end if
      if (size(spec%tv%orders) > 0) then
         call tv_component(spec%tv, xv, result%g, status)
      else
         allocate(result%g(n)); result%g=1.0_dp; status=0
      end if
      if (status /= 0) then; result%status=2; return; end if
      call garchx_simulate(n, spec%garch, par_h, yh, result%h, result%innovations, status, &
                           xreg=xreg, supplied_innovations=innovations)
      if (status /= 0) then; result%status=3; return; end if
      allocate(result%sigma2(n), result%y(n))
      result%sigma2 = result%g*result%h
      result%y = sqrt(result%g)*yh
      result%status = 0
   end subroutine tvgarch_simulate

   subroutine tvgarch_forecast(fit, n_ahead, forecast, status, new_xtv, new_xreg, n_sim, h_paths)
      type(tvgarch_fit), intent(in) :: fit
      integer, intent(in) :: n_ahead
      real(dp), allocatable, intent(out) :: forecast(:)
      integer, intent(out) :: status
      real(dp), intent(in), optional :: new_xtv(:), new_xreg(:, :)
      integer, intent(in), optional :: n_sim
      real(dp), allocatable, intent(out), optional :: h_paths(:, :)
      real(dp), allocatable :: hf(:), gf(:), xv(:), paths(:, :)
      integer :: i, st
      call garchx_forecast(fit%hfit, n_ahead, hf, st, future_xreg=new_xreg, n_sim=n_sim, paths=paths)
      if (st /= 0) then; status=st; allocate(forecast(0)); return; end if
      if (size(fit%spec%tv%orders) > 0) then
         if (present(new_xtv)) then
            if (size(new_xtv) /= n_ahead) then; status=10; allocate(forecast(0)); return; end if
            xv = new_xtv
            call tv_component(fit%spec%tv, xv, gf, st)
         else
            allocate(gf(n_ahead)); gf=fit%par_g(1); st=0
         end if
      else
         allocate(gf(n_ahead)); gf=1.0_dp; st=0
      end if
      if (st /= 0) then; status=11; allocate(forecast(0)); return; end if
      allocate(forecast(n_ahead)); forecast=hf*gf
      if (present(h_paths)) then
         allocate(h_paths(size(paths,1),size(paths,2)))
         do i=1,size(paths,2); h_paths(:,i)=paths(:,i)*gf; end do
      end if
      status=0
   end subroutine tvgarch_forecast

   subroutine tvgarch_quantile_path(fit, probs, paths, status)
      type(tvgarch_fit), intent(in) :: fit
      real(dp), intent(in) :: probs(:)
      real(dp), allocatable, intent(out) :: paths(:, :)
      integer, intent(out) :: status
      integer :: j
      if (any(probs < 0.0_dp) .or. any(probs > 1.0_dp)) then
         status=1; allocate(paths(0,0)); return
      end if
      allocate(paths(size(fit%y),size(probs)))
      do j=1,size(probs)
         paths(:,j)=sqrt(fit%sigma2)*empirical_quantile(fit%residuals,probs(j))
      end do
      status=0
   end subroutine tvgarch_quantile_path

   subroutine fitted_tvgarch(fit, component, values, status)
      type(tvgarch_fit), intent(in) :: fit
      character(len=*), intent(in) :: component
      real(dp), allocatable, intent(out) :: values(:)
      integer, intent(out) :: status
      select case(trim(adjustl(component)))
      case('tvgarch','sigma2')
         values=fit%sigma2; status=0
      case('garch','h')
         values=fit%h; status=0
      case('tv','g')
         values=fit%g; status=0
      case default
         allocate(values(0)); status=1
      end select
   end subroutine fitted_tvgarch
end module tvgarch_model
