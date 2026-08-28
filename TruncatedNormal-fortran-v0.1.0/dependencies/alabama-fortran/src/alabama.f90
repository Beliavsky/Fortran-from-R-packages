! Modern Fortran computational port of the R package alabama.
! Original package copyright Ravi Varadhan and contributors.
! SPDX-License-Identifier: GPL-2.0-or-later
module alabama
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use numderiv, only : nd_dp => dp, grad, jacobian, deriv_options
   use roptim_mod, only : roptim_control_t, roptim_result_t, roptim_minimize, &
      method_bfgs, method_nelder_mead, method_cg, method_lbfgsb, method_sann
   implicit none
   private

   integer, parameter, public :: dp = kind(1.0d0)

   integer, parameter, public :: al_success = 0
   integer, parameter, public :: al_outer_limit = 7
   integer, parameter, public :: al_lack_progress = 9
   integer, parameter, public :: al_objective_reversal = 11
   integer, parameter, public :: al_invalid_input = 20
   integer, parameter, public :: al_inner_failure = 21

   type, public :: alabama_outer_control_t
      real(dp) :: lam0 = 10.0_dp
      real(dp) :: sig0 = 100.0_dp
      real(dp) :: mu0 = 0.01_dp
      real(dp) :: eps = 1.0e-7_dp
      integer :: itmax = 50
      integer :: ilack_max = 6
      character(len=16) :: method = 'BFGS'
      logical :: trace = .false.
      logical :: nm_init = .false.
      logical :: maximize = .false.
      logical :: kkt2_check = .true.
      real(dp), allocatable :: i_scale(:)
      real(dp), allocatable :: e_scale(:)
   end type alabama_outer_control_t

   type, public :: alabama_inner_control_t
      integer :: max_iterations = 100
      real(dp) :: reltol = sqrt(epsilon(1.0_dp))
      real(dp) :: ndeps = 1.0e-3_dp
      integer :: trace = 0
   end type alabama_inner_control_t

   type, public :: alabama_result_t
      real(dp), allocatable :: par(:)
      real(dp) :: value = huge(1.0_dp)
      integer :: convergence = al_invalid_input
      character(len=:), allocatable :: message
      integer :: outer_iterations = 0
      integer :: function_evaluations = 0
      integer :: gradient_evaluations = 0
      real(dp), allocatable :: lambda(:)
      real(dp) :: sigma = 0.0_dp
      real(dp) :: barrier_value = 0.0_dp
      real(dp) :: penalty = 0.0_dp
      real(dp) :: constraint_norm = 0.0_dp
      real(dp), allocatable :: gradient(:)
      real(dp), allocatable :: hessian(:, :)
      real(dp), allocatable :: ineq(:)
      real(dp), allocatable :: equal(:)
      logical :: kkt1 = .false.
      logical :: kkt2 = .false.
      logical :: kkt2_available = .false.
   end type alabama_result_t

   abstract interface
      function objective_function(x) result(value)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp) :: value
      end function objective_function

      subroutine gradient_function(x, g)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: g(:)
      end subroutine gradient_function

      function constraint_function(x) result(value)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp), allocatable :: value(:)
      end function constraint_function

      function constraint_jacobian_function(x) result(value)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp), allocatable :: value(:, :)
      end function constraint_jacobian_function
   end interface

   public :: auglag, auglag1, auglag2, auglag3
   public :: constr_optim_nl, adpbar, augpen, alabama_legacy

contains

   subroutine auglag(par, fn, result, gr, hin, hin_jac, heq, heq_jac, &
                     control_outer, control_inner)
      real(dp), intent(inout) :: par(:)
      procedure(objective_function) :: fn
      type(alabama_result_t), intent(out) :: result
      procedure(gradient_function), optional :: gr
      procedure(constraint_function), optional :: hin, heq
      procedure(constraint_jacobian_function), optional :: hin_jac, heq_jac
      type(alabama_outer_control_t), intent(in), optional :: control_outer
      type(alabama_inner_control_t), intent(in), optional :: control_inner

      if (.not. present(hin) .and. .not. present(heq)) then
         call fail_result(result, par, al_invalid_input, &
            'auglag requires at least one equality or inequality constraint')
         return
      end if
      call auglag_core(par, fn, result, gr, hin, hin_jac, heq, heq_jac, &
                       control_outer, control_inner)
   end subroutine auglag

   subroutine auglag1(par, fn, heq, result, gr, heq_jac, control_outer, control_inner)
      real(dp), intent(inout) :: par(:)
      procedure(objective_function) :: fn
      procedure(constraint_function) :: heq
      type(alabama_result_t), intent(out) :: result
      procedure(gradient_function), optional :: gr
      procedure(constraint_jacobian_function), optional :: heq_jac
      type(alabama_outer_control_t), intent(in), optional :: control_outer
      type(alabama_inner_control_t), intent(in), optional :: control_inner
      call auglag_core(par, fn, result, gr, heq=heq, heq_jac=heq_jac, &
                       control_outer=control_outer, control_inner=control_inner)
   end subroutine auglag1

   subroutine auglag2(par, fn, hin, result, gr, hin_jac, control_outer, control_inner)
      real(dp), intent(inout) :: par(:)
      procedure(objective_function) :: fn
      procedure(constraint_function) :: hin
      type(alabama_result_t), intent(out) :: result
      procedure(gradient_function), optional :: gr
      procedure(constraint_jacobian_function), optional :: hin_jac
      type(alabama_outer_control_t), intent(in), optional :: control_outer
      type(alabama_inner_control_t), intent(in), optional :: control_inner
      call auglag_core(par, fn, result, gr, hin=hin, hin_jac=hin_jac, &
                       control_outer=control_outer, control_inner=control_inner)
   end subroutine auglag2

   subroutine auglag3(par, fn, hin, heq, result, gr, hin_jac, heq_jac, &
                      control_outer, control_inner)
      real(dp), intent(inout) :: par(:)
      procedure(objective_function) :: fn
      procedure(constraint_function) :: hin, heq
      type(alabama_result_t), intent(out) :: result
      procedure(gradient_function), optional :: gr
      procedure(constraint_jacobian_function), optional :: hin_jac, heq_jac
      type(alabama_outer_control_t), intent(in), optional :: control_outer
      type(alabama_inner_control_t), intent(in), optional :: control_inner
      call auglag_core(par, fn, result, gr, hin, hin_jac, heq, heq_jac, &
                       control_outer, control_inner)
   end subroutine auglag3

   subroutine auglag_core(par, fn, result, gr, hin, hin_jac, heq, heq_jac, &
                          control_outer, control_inner)
      real(dp), intent(inout) :: par(:)
      procedure(objective_function) :: fn
      type(alabama_result_t), intent(out) :: result
      procedure(gradient_function), optional :: gr
      procedure(constraint_function), optional :: hin, heq
      procedure(constraint_jacobian_function), optional :: hin_jac, heq_jac
      type(alabama_outer_control_t), intent(in), optional :: control_outer
      type(alabama_inner_control_t), intent(in), optional :: control_inner

      type(alabama_outer_control_t) :: oc
      type(alabama_inner_control_t) :: ic
      type(roptim_control_t) :: rc
      type(roptim_result_t) :: ir
      real(dp), allocatable :: h(:), e(:), d(:), lam(:), iscale(:), escale(:)
      real(dp), allocatable :: oldpar(:), g(:), hs(:, :)
      real(dp) :: sig, kprev, kval, r, rold, obj, pconv, pfact
      integer :: mi, me, m, i, ilack, stat
      character(:), allocatable :: meth
      logical :: converged

      call init_result(result, par)
      if (size(par) == 0 .or. .not. all(ieee_is_finite(par))) then
         call fail_result(result, par, al_invalid_input, 'invalid starting parameters')
         return
      end if
      oc = alabama_outer_control_t()
      if (present(control_outer)) oc = control_outer
      ic = alabama_inner_control_t()
      if (present(control_inner)) ic = control_inner
      if (oc%itmax <= 0 .or. oc%eps <= 0.0_dp .or. oc%sig0 <= 0.0_dp) then
         call fail_result(result, par, al_invalid_input, 'invalid outer control values')
         return
      end if

      mi = 0
      me = 0
      if (present(hin)) then
         h = invoke_constraint(hin, par)
         mi = size(h)
      else
         allocate(h(0))
      end if
      if (present(heq)) then
         e = invoke_constraint(heq, par)
         me = size(e)
      else
         allocate(e(0))
      end if
      if (mi + me == 0) then
         call fail_result(result, par, al_invalid_input, 'constraint vectors must not both be empty')
         return
      end if
      if (.not. all(ieee_is_finite(h)) .or. .not. all(ieee_is_finite(e))) then
         call fail_result(result, par, al_invalid_input, 'constraints are nonfinite at the start')
         return
      end if
      call make_scale(mi, oc%i_scale, iscale, stat)
      if (stat /= 0) then
         call fail_result(result, par, al_invalid_input, 'i_scale has wrong size or invalid values')
         return
      end if
      call make_scale(me, oc%e_scale, escale, stat)
      if (stat /= 0) then
         call fail_result(result, par, al_invalid_input, 'e_scale has wrong size or invalid values')
         return
      end if
      if (mi > 0) h = h / iscale
      if (me > 0) e = e / escale

      m = mi + me
      allocate(d(m), lam(m), oldpar(size(par)), g(size(par)))
      lam = oc%lam0
      sig = oc%sig0
      call make_d(h, e, lam, sig, d)
      kprev = maxval(abs(d))
      if (kprev > 0.0_dp) then
         sig = oc%sig0 / kprev
      else
         sig = 1.0_dp
      end if
      call make_d(h, e, lam, sig, d)
      obj = invoke_objective(fn, par)
      r = obj
      ilack = 0
      kval = huge(1.0_dp)
      pfact = merge(-1.0_dp, 1.0_dp, oc%maximize)
      converged = .false.

      do i = 1, oc%itmax
         if (oc%trace) call trace_aug(i, par, obj, h, e)
         oldpar = par
         rold = r

         call map_inner_control(ic, oc%maximize, size(par), rc)
         meth = canonical_inner_method(oc%method)
         if (oc%nm_init .and. i == 1) meth = method_nelder_mead
         call roptim_minimize(par, augmented_objective, ir, method=meth, &
                              gradient=augmented_gradient, control=rc)
         result%function_evaluations = result%function_evaluations + ir%function_evaluations
         result%gradient_evaluations = result%gradient_evaluations + ir%gradient_evaluations
         if (ir%convergence /= 0 .and. ir%convergence /= 1) then
            call fail_result(result, par, al_inner_failure, 'inner optimizer failed: '//ir%message)
            return
         end if
         r = ir%value
         if (mi > 0) h = invoke_constraint(hin, par) / iscale
         if (me > 0) e = invoke_constraint(heq, par) / escale
         call make_d(h, e, lam, sig, d)
         kval = maxval(abs(d))
         if (kval <= kprev / 4.0_dp) then
            lam = lam - d * sig
            kprev = kval
         else
            sig = 10.0_dp * sig
         end if
         obj = invoke_objective(fn, par)
         pconv = maxval(abs(par - oldpar))
         if (pconv < oc%eps) then
            ilack = ilack + 1
         else
            ilack = 0
         end if
         if ((ieee_is_finite(r) .and. ieee_is_finite(rold) .and. &
              abs(r - rold) < oc%eps .and. kval < oc%eps) .or. &
              ilack >= oc%ilack_max) then
            converged = .true.
            exit
         end if
      end do

      result%outer_iterations = min(i, oc%itmax)
      result%par = par
      result%value = invoke_objective(fn, par)
      result%sigma = sig
      result%lambda = lam
      result%constraint_norm = kval
      allocate(result%ineq(mi), result%equal(me), result%gradient(size(par)))
      if (mi > 0) result%ineq = invoke_constraint(hin, par)
      if (me > 0) result%equal = invoke_constraint(heq, par)
      call augmented_gradient_plain(par, result%gradient)
      result%kkt1 = maxval(abs(result%gradient)) <= 0.01_dp * (1.0_dp + abs(result%value))

      if (.not. converged .and. result%outer_iterations >= oc%itmax) then
         result%convergence = al_outer_limit
         result%message = 'ALABaMA ran out of outer iterations'
      else if (kval > oc%eps) then
         result%convergence = al_lack_progress
         result%message = 'convergence due to lack of progress in parameter updates'
      else
         result%convergence = al_success
         result%message = 'successful convergence'
      end if

      if (oc%kkt2_check) then
         allocate(result%hessian(size(par), size(par)))
         call jacobian(augmented_gradient_vector, par, result%hessian, method='simple', status=stat)
         if (stat == 0) then
            result%hessian = 0.5_dp * (result%hessian + transpose(result%hessian))
            allocate(hs(size(par), size(par)))
            hs = pfact * result%hessian
            result%kkt2 = is_positive_definite(hs)
            result%kkt2_available = .true.
         end if
      end if

   contains

      function augmented_objective(x, user_data) result(value)
         real(dp), intent(in) :: x(:)
         class(*), intent(inout), optional :: user_data
         real(dp) :: value
         real(dp), allocatable :: hh(:), ee(:), dd(:)
         if (present(user_data)) continue
         allocate(dd(m))
         if (mi > 0) then
            hh = invoke_constraint(hin, x) / iscale
         else
            allocate(hh(0))
         end if
         if (me > 0) then
            ee = invoke_constraint(heq, x) / escale
         else
            allocate(ee(0))
         end if
         call make_d(hh, ee, lam, sig, dd)
         value = invoke_objective(fn, x) - pfact * dot_product(lam, dd) + &
                 pfact * 0.5_dp * sig * dot_product(dd, dd)
      end function augmented_objective

      subroutine augmented_gradient(x, gout, user_data)
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: gout(:)
         class(*), intent(inout), optional :: user_data
         if (present(user_data)) continue
         call augmented_gradient_plain(x, gout)
      end subroutine augmented_gradient

      subroutine augmented_gradient_plain(x, gout)
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: gout(:)
         real(dp), allocatable :: hh(:), ee(:), jh(:, :), je(:, :), gg(:)
         integer :: j
         allocate(gg(size(x)))
         call objective_gradient(x, gg)
         gout = gg
         if (mi > 0) then
            hh = invoke_constraint(hin, x) / iscale
            call inequality_jacobian(x, jh)
            do j = 1, mi
               if (hh(j) <= lam(j) / sig) then
                  gout = gout + (-pfact * lam(j) + pfact * sig * hh(j)) * jh(j, :)
               end if
            end do
         end if
         if (me > 0) then
            ee = invoke_constraint(heq, x) / escale
            call equality_jacobian(x, je)
            do j = 1, me
               gout = gout + (-pfact * lam(mi+j) + pfact * sig * ee(j)) * je(j, :)
            end do
         end if
      end subroutine augmented_gradient_plain

      function augmented_gradient_vector(x) result(value)
         real(nd_dp), intent(in) :: x(:)
         real(nd_dp), allocatable :: value(:)
         allocate(value(size(x)))
         call augmented_gradient_plain(x, value)
      end function augmented_gradient_vector

      subroutine objective_gradient(x, gout)
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: gout(:)
         integer :: s
         type(deriv_options) :: opts
         if (present(gr)) then
            call invoke_gradient(gr, x, gout)
         else
            opts = deriv_options(eps=ic%ndeps)
            call grad(fn, x, gout, method='simple', options=opts, status=s)
            if (s /= 0) gout = huge(1.0_dp)
         end if
      end subroutine objective_gradient

      subroutine inequality_jacobian(x, jout)
         real(dp), intent(in) :: x(:)
         real(dp), allocatable, intent(out) :: jout(:, :)
         integer :: s
         type(deriv_options) :: opts
         allocate(jout(mi, size(x)))
         if (present(hin_jac)) then
            jout = invoke_constraint_jacobian(hin_jac, x)
         else
            opts = deriv_options(eps=ic%ndeps)
            call jacobian(hin, x, jout, method='simple', options=opts, status=s)
            if (s /= 0) jout = huge(1.0_dp)
         end if
         jout = jout / spread(iscale, 2, size(x))
      end subroutine inequality_jacobian

      subroutine equality_jacobian(x, jout)
         real(dp), intent(in) :: x(:)
         real(dp), allocatable, intent(out) :: jout(:, :)
         integer :: s
         type(deriv_options) :: opts
         allocate(jout(me, size(x)))
         if (present(heq_jac)) then
            jout = invoke_constraint_jacobian(heq_jac, x)
         else
            opts = deriv_options(eps=ic%ndeps)
            call jacobian(heq, x, jout, method='simple', options=opts, status=s)
            if (s /= 0) jout = huge(1.0_dp)
         end if
         jout = jout / spread(escale, 2, size(x))
      end subroutine equality_jacobian

   end subroutine auglag_core

   subroutine constr_optim_nl(par, fn, result, gr, hin, hin_jac, heq, heq_jac, &
                              control_outer, control_inner)
      real(dp), intent(inout) :: par(:)
      procedure(objective_function) :: fn
      type(alabama_result_t), intent(out) :: result
      procedure(gradient_function), optional :: gr
      procedure(constraint_function), optional :: hin, heq
      procedure(constraint_jacobian_function), optional :: hin_jac, heq_jac
      type(alabama_outer_control_t), intent(in), optional :: control_outer
      type(alabama_inner_control_t), intent(in), optional :: control_inner
      if (.not. present(hin) .and. .not. present(heq)) then
         call fail_result(result, par, al_invalid_input, &
            'constr_optim_nl requires at least one constraint')
      else if (.not. present(hin)) then
         call augpen(par, fn, heq, result, gr, heq_jac, control_outer, control_inner)
      else if (.not. present(heq)) then
         call adpbar(par, fn, hin, result, gr, hin_jac, control_outer, control_inner)
      else
         call alabama_legacy(par, fn, hin, heq, result, gr, hin_jac, heq_jac, &
                             control_outer, control_inner)
      end if
   end subroutine constr_optim_nl

   subroutine augpen(par, fn, heq, result, gr, heq_jac, control_outer, control_inner)
      real(dp), intent(inout) :: par(:)
      procedure(objective_function) :: fn
      procedure(constraint_function) :: heq
      type(alabama_result_t), intent(out) :: result
      procedure(gradient_function), optional :: gr
      procedure(constraint_jacobian_function), optional :: heq_jac
      type(alabama_outer_control_t), intent(in), optional :: control_outer
      type(alabama_inner_control_t), intent(in), optional :: control_inner
      type(alabama_outer_control_t) :: oc
      type(alabama_inner_control_t) :: ic
      type(roptim_control_t) :: rc
      type(roptim_result_t) :: ir
      real(dp), allocatable :: eq(:), lam(:)
      real(dp) :: sig, kprev, kval, pfact, rval
      integer :: k
      character(:), allocatable :: meth

      call init_result(result, par)
      oc = alabama_outer_control_t(mu0=0.01_dp, sig0=10.0_dp, eps=1.0e-7_dp)
      if (present(control_outer)) oc = control_outer
      ic = alabama_inner_control_t()
      if (present(control_inner)) ic = control_inner
      eq = invoke_constraint(heq, par)
      if (size(eq) == 0) then
         call fail_result(result, par, al_invalid_input, 'empty equality constraint vector')
         return
      end if
      allocate(lam(size(eq)))
      lam = 0.0_dp
      kprev = maxval(abs(eq))
      if (kprev > 0.0_dp) then
         sig = oc%sig0 / kprev
      else
         sig = 1.0_dp
      end if
      kval = huge(1.0_dp)
      rval = invoke_objective(fn, par)
      pfact = merge(-1.0_dp, 1.0_dp, oc%maximize)
      do k = 1, oc%itmax
         call map_inner_control(ic, oc%maximize, size(par), rc)
         meth = canonical_inner_method(oc%method)
         if (oc%nm_init .and. k == 1) meth = method_nelder_mead
         call roptim_minimize(par, penalty_objective, ir, method=meth, &
                              gradient=penalty_gradient, control=rc)
         result%function_evaluations = result%function_evaluations + ir%function_evaluations
         result%gradient_evaluations = result%gradient_evaluations + ir%gradient_evaluations
         rval = ir%value
         eq = invoke_constraint(heq, par)
         kval = maxval(abs(eq))
         if (kval <= kprev / 4.0_dp) then
            lam = lam - eq * sig
            kprev = kval
         else
            sig = 10.0_dp * sig
         end if
         if (kval <= oc%eps) exit
      end do
      result%par = par
      result%value = invoke_objective(fn, par)
      result%outer_iterations = min(k, oc%itmax)
      result%sigma = sig
      result%lambda = lam
      result%constraint_norm = kval
      result%penalty = rval - result%value
      allocate(result%equal(size(eq)), result%ineq(0))
      result%equal = eq
      if (k >= oc%itmax .and. kval > oc%eps) then
         result%convergence = al_outer_limit
         result%message = 'augmented Lagrangian algorithm ran out of iterations'
      else
         result%convergence = al_success
         result%message = 'successful convergence'
      end if

   contains
      function penalty_objective(x, user_data) result(v)
         real(dp), intent(in) :: x(:)
         class(*), intent(inout), optional :: user_data
         real(dp) :: v
         real(dp), allocatable :: ee(:)
         if (present(user_data)) continue
         ee = invoke_constraint(heq, x)
         v = invoke_objective(fn, x) - pfact * dot_product(lam, ee) + &
             pfact * 0.5_dp * sig * dot_product(ee, ee)
      end function penalty_objective
      subroutine penalty_gradient(x, gout, user_data)
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: gout(:)
         class(*), intent(inout), optional :: user_data
         real(dp), allocatable :: ee(:), jj(:, :), gg(:)
         integer :: s
         type(deriv_options) :: opts
         if (present(user_data)) continue
         allocate(gg(size(x)))
         if (present(gr)) then
            call invoke_gradient(gr, x, gg)
         else
            opts = deriv_options(eps=ic%ndeps)
            call grad(fn, x, gg, method='simple', options=opts, status=s)
         end if
         ee = invoke_constraint(heq, x)
         allocate(jj(size(ee), size(x)))
         if (present(heq_jac)) then
            jj = invoke_constraint_jacobian(heq_jac, x)
         else
            opts = deriv_options(eps=ic%ndeps)
            call jacobian(heq, x, jj, method='simple', options=opts, status=s)
         end if
         gout = gg - pfact * matmul(transpose(jj), lam) + &
                pfact * sig * matmul(transpose(jj), ee)
      end subroutine penalty_gradient
   end subroutine augpen

   subroutine adpbar(par, fn, hin, result, gr, hin_jac, control_outer, control_inner)
      real(dp), intent(inout) :: par(:)
      procedure(objective_function) :: fn
      procedure(constraint_function) :: hin
      type(alabama_result_t), intent(out) :: result
      procedure(gradient_function), optional :: gr
      procedure(constraint_jacobian_function), optional :: hin_jac
      type(alabama_outer_control_t), intent(in), optional :: control_outer
      type(alabama_inner_control_t), intent(in), optional :: control_inner
      type(alabama_outer_control_t) :: oc
      type(alabama_inner_control_t) :: ic
      type(roptim_control_t) :: rc
      type(roptim_result_t) :: ir
      real(dp), allocatable :: oldpar(:), h(:)
      real(dp) :: mu, obj, objold, r, rold
      integer :: i
      character(:), allocatable :: meth

      call init_result(result, par)
      oc = alabama_outer_control_t(mu0=0.01_dp, sig0=10.0_dp, eps=1.0e-7_dp)
      if (present(control_outer)) oc = control_outer
      ic = alabama_inner_control_t()
      if (present(control_inner)) ic = control_inner
      h = invoke_constraint(hin, par)
      if (size(h) == 0 .or. any(h <= 0.0_dp)) then
         call fail_result(result, par, al_invalid_input, 'initial value is not strictly inequality-feasible')
         return
      end if
      allocate(oldpar(size(par)))
      mu = oc%mu0
      if (oc%maximize) mu = -mu
      obj = invoke_objective(fn, par)
      ! Preserve the R routine's first r value, evaluated before mu is scaled.
      r = barrier_at(par, par)
      mu = mu * minval(h)
      if (abs(mu) < 1.0e-10_dp) mu = 1.0e-4_dp * sign(1.0_dp, mu)

      do i = 1, oc%itmax
         oldpar = par
         objold = obj
         rold = r
         call map_inner_control(ic, oc%maximize, size(par), rc)
         meth = canonical_inner_method(oc%method)
         if (oc%nm_init .and. i == 1) meth = method_nelder_mead
         call roptim_minimize(par, barrier_objective, ir, method=meth, &
                              gradient=barrier_gradient, control=rc)
         result%function_evaluations = result%function_evaluations + ir%function_evaluations
         result%gradient_evaluations = result%gradient_evaluations + ir%gradient_evaluations
         r = ir%value
         obj = invoke_objective(fn, par)
         if (ieee_is_finite(r) .and. ieee_is_finite(rold) .and. abs(r-rold) < oc%eps) exit
         if ((.not. oc%maximize .and. obj > objold) .or. &
             (oc%maximize .and. obj < objold)) exit
      end do
      result%par = par
      result%value = invoke_objective(fn, par)
      result%outer_iterations = min(i, oc%itmax)
      result%barrier_value = r - result%value
      allocate(result%ineq(size(h)), result%equal(0))
      result%ineq = invoke_constraint(hin, par)
      if (i >= oc%itmax) then
         result%convergence = al_outer_limit
         result%message = 'barrier algorithm ran out of iterations'
      else if ((.not. oc%maximize .and. obj > objold) .or. &
               (oc%maximize .and. obj < objold)) then
         result%convergence = al_objective_reversal
         result%message = 'objective reversed direction during barrier outer iteration'
      else
         result%convergence = al_success
         result%message = 'successful convergence'
      end if

   contains
      function barrier_at(x, xold) result(v)
         real(dp), intent(in) :: x(:), xold(:)
         real(dp) :: v, bar
         real(dp), allocatable :: gx(:), go(:), jo(:, :)
         gx = invoke_constraint(hin, x)
         if (any(gx <= 0.0_dp)) then
            v = huge(1.0_dp)
            return
         end if
         go = invoke_constraint(hin, xold)
         call get_hin_jac(xold, jo)
         bar = sum(go * log(gx) - matmul(jo, x))
         v = invoke_objective(fn, x) - mu * bar
      end function barrier_at
      function barrier_objective(x, user_data) result(v)
         real(dp), intent(in) :: x(:)
         class(*), intent(inout), optional :: user_data
         real(dp) :: v
         if (present(user_data)) continue
         v = barrier_at(x, oldpar)
      end function barrier_objective
      subroutine barrier_gradient(x, gout, user_data)
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: gout(:)
         class(*), intent(inout), optional :: user_data
         real(dp), allocatable :: gx(:), go(:), jo(:, :), gg(:), w(:)
         if (present(user_data)) continue
         allocate(gg(size(x)))
         call get_obj_grad(x, gg)
         gx = invoke_constraint(hin, x)
         go = invoke_constraint(hin, oldpar)
         call get_hin_jac(oldpar, jo)
         allocate(w(size(gx)))
         w = go / gx - 1.0_dp
         gout = gg - mu * matmul(transpose(jo), w)
      end subroutine barrier_gradient
      subroutine get_obj_grad(x, gout)
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: gout(:)
         integer :: s
         type(deriv_options) :: opts
         if (present(gr)) then
            call invoke_gradient(gr, x, gout)
         else
            opts = deriv_options(eps=ic%ndeps)
            call grad(fn, x, gout, method='simple', options=opts, status=s)
         end if
      end subroutine get_obj_grad
      subroutine get_hin_jac(x, jj)
         real(dp), intent(in) :: x(:)
         real(dp), allocatable, intent(out) :: jj(:, :)
         integer :: s
         type(deriv_options) :: opts
         real(dp), allocatable :: gv(:)
         gv = invoke_constraint(hin, x)
         allocate(jj(size(gv), size(x)))
         if (present(hin_jac)) then
            jj = invoke_constraint_jacobian(hin_jac, x)
         else
            opts = deriv_options(eps=ic%ndeps)
            call jacobian(hin, x, jj, method='simple', options=opts, status=s)
         end if
      end subroutine get_hin_jac
   end subroutine adpbar

   subroutine alabama_legacy(par, fn, hin, heq, result, gr, hin_jac, heq_jac, &
                             control_outer, control_inner)
      real(dp), intent(inout) :: par(:)
      procedure(objective_function) :: fn
      procedure(constraint_function) :: hin, heq
      type(alabama_result_t), intent(out) :: result
      procedure(gradient_function), optional :: gr
      procedure(constraint_jacobian_function), optional :: hin_jac, heq_jac
      type(alabama_outer_control_t), intent(in), optional :: control_outer
      type(alabama_inner_control_t), intent(in), optional :: control_inner
      type(alabama_outer_control_t) :: oc
      type(alabama_inner_control_t) :: ic
      type(roptim_control_t) :: rc
      type(roptim_result_t) :: ir
      real(dp), allocatable :: oldpar(:), h(:), e(:), lam(:)
      real(dp) :: mu, sig, kprev, kval, obj, r, rold, pfact, pconv
      integer :: i, stat
      character(:), allocatable :: meth

      call init_result(result, par)
      oc = alabama_outer_control_t(mu0=0.01_dp, sig0=10.0_dp, eps=1.0e-7_dp)
      if (present(control_outer)) oc = control_outer
      ic = alabama_inner_control_t()
      if (present(control_inner)) ic = control_inner
      h = invoke_constraint(hin, par)
      e = invoke_constraint(heq, par)
      if (size(h) == 0 .or. any(h <= 0.0_dp)) then
         call fail_result(result, par, al_invalid_input, 'initial value violates inequality constraints')
         return
      end if
      if (size(e) == 0) then
         call fail_result(result, par, al_invalid_input, 'empty equality constraint vector')
         return
      end if
      allocate(oldpar(size(par)), lam(size(e)))
      lam = 0.0_dp
      pfact = merge(-1.0_dp, 1.0_dp, oc%maximize)
      mu = merge(-oc%mu0, oc%mu0, oc%maximize)
      kprev = maxval(abs(e))
      if (kprev > 0.0_dp) then
         sig = oc%sig0 / kprev
      else
         sig = 1.0_dp
      end if
      kval = huge(1.0_dp)
      obj = invoke_objective(fn, par)
      ! R initializes r with the barrier term alone, before scaling mu.
      r = invoke_objective(fn, par) - mu * barrier_component(par, par)
      mu = mu * minval(h)
      if (abs(mu) < 1.0e-10_dp) mu = 1.0e-4_dp * sign(1.0_dp, mu)

      do i = 1, oc%itmax
         oldpar = par
         rold = r
         call map_inner_control(ic, oc%maximize, size(par), rc)
         meth = canonical_inner_method(oc%method)
         if (oc%nm_init .and. i == 1) meth = method_nelder_mead
         call roptim_minimize(par, combined_objective, ir, method=meth, &
                              gradient=combined_gradient, control=rc)
         result%function_evaluations = result%function_evaluations + ir%function_evaluations
         result%gradient_evaluations = result%gradient_evaluations + ir%gradient_evaluations
         r = ir%value
         h = invoke_constraint(hin, par)
         e = invoke_constraint(heq, par)
         kval = maxval(abs(e))
         if (kval <= kprev / 4.0_dp) then
            lam = lam - e * sig
            kprev = kval
         else
            sig = 10.0_dp * sig
         end if
         obj = invoke_objective(fn, par)
         pconv = maxval(abs(par-oldpar))
         if ((ieee_is_finite(r) .and. ieee_is_finite(rold) .and. &
              abs(r-rold) < oc%eps .and. kval < oc%eps) .or. pconv < 1.0e-12_dp) exit
      end do

      result%par = par
      result%value = invoke_objective(fn, par)
      result%outer_iterations = min(i, oc%itmax)
      result%lambda = lam
      result%sigma = sig
      result%barrier_value = r - result%value
      result%constraint_norm = kval
      allocate(result%ineq(size(h)), result%equal(size(e)))
      result%ineq = h
      result%equal = e
      if (i >= oc%itmax .and. kval > oc%eps) then
         result%convergence = al_outer_limit
         result%message = 'ALABaMA barrier algorithm ran out of iterations'
      else
         result%convergence = al_success
         result%message = 'successful convergence'
      end if
      allocate(result%hessian(size(par), size(par)))
      call jacobian(combined_gradient_vector, par, result%hessian, method='simple', status=stat)

   contains
      function barrier_component(x, xold) result(v)
         real(dp), intent(in) :: x(:), xold(:)
         real(dp) :: v
         real(dp), allocatable :: gx(:), go(:), jo(:, :)
         gx = invoke_constraint(hin, x)
         if (any(gx <= 0.0_dp)) then
            v = huge(1.0_dp)
            return
         end if
         go = invoke_constraint(hin, xold)
         call get_hin_jac2(xold, jo)
         v = sum(go * log(gx) - matmul(jo, x))
      end function barrier_component
      function combined_at(x, xold) result(v)
         real(dp), intent(in) :: x(:), xold(:)
         real(dp) :: v
         real(dp), allocatable :: ee(:)
         ee = invoke_constraint(heq, x)
         v = invoke_objective(fn, x) - mu * barrier_component(x, xold) - pfact * dot_product(lam, ee) + &
             pfact * 0.5_dp * sig * dot_product(ee, ee)
      end function combined_at
      function combined_objective(x, user_data) result(v)
         real(dp), intent(in) :: x(:)
         class(*), intent(inout), optional :: user_data
         real(dp) :: v
         if (present(user_data)) continue
         v = combined_at(x, oldpar)
      end function combined_objective
      subroutine combined_gradient(x, gout, user_data)
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: gout(:)
         class(*), intent(inout), optional :: user_data
         if (present(user_data)) continue
         call combined_gradient_plain(x, gout)
      end subroutine combined_gradient
      subroutine combined_gradient_plain(x, gout)
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: gout(:)
         real(dp), allocatable :: gx(:), go(:), ji(:, :), je(:, :), gg(:), ee(:), w(:)
         allocate(gg(size(x)))
         call get_obj_grad2(x, gg)
         gx = invoke_constraint(hin, x)
         go = invoke_constraint(hin, oldpar)
         call get_hin_jac2(oldpar, ji)
         allocate(w(size(gx)))
         w = go/gx - 1.0_dp
         ee = invoke_constraint(heq, x)
         call get_heq_jac2(x, je)
         gout = gg - mu * matmul(transpose(ji), w) - pfact * matmul(transpose(je), lam) + &
                pfact * sig * matmul(transpose(je), ee)
      end subroutine combined_gradient_plain
      function combined_gradient_vector(x) result(v)
         real(nd_dp), intent(in) :: x(:)
         real(nd_dp), allocatable :: v(:)
         allocate(v(size(x)))
         call combined_gradient_plain(x, v)
      end function combined_gradient_vector
      subroutine get_obj_grad2(x, gout)
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: gout(:)
         integer :: s
         type(deriv_options) :: opts
         if (present(gr)) then
            call invoke_gradient(gr, x, gout)
         else
            opts = deriv_options(eps=ic%ndeps)
            call grad(fn, x, gout, method='simple', options=opts, status=s)
         end if
      end subroutine get_obj_grad2
      subroutine get_hin_jac2(x, jj)
         real(dp), intent(in) :: x(:)
         real(dp), allocatable, intent(out) :: jj(:, :)
         real(dp), allocatable :: tmp(:)
         integer :: s
         type(deriv_options) :: opts
         tmp = invoke_constraint(hin, x)
         allocate(jj(size(tmp), size(x)))
         if (present(hin_jac)) then
            jj = invoke_constraint_jacobian(hin_jac, x)
         else
            opts = deriv_options(eps=ic%ndeps)
            call jacobian(hin, x, jj, method='simple', options=opts, status=s)
         end if
      end subroutine get_hin_jac2
      subroutine get_heq_jac2(x, jj)
         real(dp), intent(in) :: x(:)
         real(dp), allocatable, intent(out) :: jj(:, :)
         integer :: s
         type(deriv_options) :: opts
         allocate(jj(size(e), size(x)))
         if (present(heq_jac)) then
            jj = invoke_constraint_jacobian(heq_jac, x)
         else
            opts = deriv_options(eps=ic%ndeps)
            call jacobian(heq, x, jj, method='simple', options=opts, status=s)
         end if
      end subroutine get_heq_jac2
   end subroutine alabama_legacy

   function invoke_objective(callback, x) result(value)
      procedure(objective_function) :: callback
      real(dp), intent(in) :: x(:)
      real(dp) :: value
      value = callback(x)
   end function invoke_objective

   subroutine invoke_gradient(callback, x, g)
      procedure(gradient_function) :: callback
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: g(:)
      call callback(x, g)
   end subroutine invoke_gradient

   function invoke_constraint(callback, x) result(value)
      procedure(constraint_function) :: callback
      real(dp), intent(in) :: x(:)
      real(dp), allocatable :: value(:)
      value = callback(x)
   end function invoke_constraint

   function invoke_constraint_jacobian(callback, x) result(value)
      procedure(constraint_jacobian_function) :: callback
      real(dp), intent(in) :: x(:)
      real(dp), allocatable :: value(:, :)
      value = callback(x)
   end function invoke_constraint_jacobian

   subroutine init_result(result, par)
      type(alabama_result_t), intent(out) :: result
      real(dp), intent(in) :: par(:)
      allocate(result%par(size(par)))
      result%par = par
      result%message = ''
   end subroutine init_result

   subroutine fail_result(result, par, code, msg)
      type(alabama_result_t), intent(inout) :: result
      real(dp), intent(in) :: par(:)
      integer, intent(in) :: code
      character(len=*), intent(in) :: msg
      if (.not. allocated(result%par)) allocate(result%par(size(par)))
      result%par = par
      result%convergence = code
      result%message = msg
   end subroutine fail_result

   subroutine make_scale(n, requested, scale, status)
      integer, intent(in) :: n
      real(dp), allocatable, intent(in) :: requested(:)
      real(dp), allocatable, intent(out) :: scale(:)
      integer, intent(out) :: status
      allocate(scale(n))
      scale = 1.0_dp
      status = 0
      if (allocated(requested)) then
         if (size(requested) == 1) then
            scale = requested(1)
         else if (size(requested) == n) then
            scale = requested
         else
            status = 1
            return
         end if
      end if
      if (n > 0) then
         if (any(.not. ieee_is_finite(scale)) .or. any(scale == 0.0_dp)) status = 1
      end if
   end subroutine make_scale

   subroutine make_d(h, e, lam, sig, d)
      real(dp), intent(in) :: h(:), e(:), lam(:), sig
      real(dp), intent(out) :: d(:)
      integer :: j, mi
      mi = size(h)
      if (mi > 0) then
         do j = 1, mi
            if (h(j) > lam(j)/sig) then
               d(j) = lam(j)/sig
            else
               d(j) = h(j)
            end if
         end do
      end if
      if (size(e) > 0) d(mi+1:) = e
   end subroutine make_d

   subroutine map_inner_control(ic, maximize, n, rc)
      type(alabama_inner_control_t), intent(in) :: ic
      logical, intent(in) :: maximize
      integer, intent(in) :: n
      type(roptim_control_t), intent(out) :: rc
      rc = roptim_control_t()
      rc%max_iterations = ic%max_iterations
      rc%reltol = ic%reltol
      rc%trace = ic%trace
      rc%fnscale = merge(-1.0_dp, 1.0_dp, maximize)
      allocate(rc%ndeps(n))
      rc%ndeps = ic%ndeps
   end subroutine map_inner_control

   function canonical_inner_method(method) result(out)
      character(len=*), intent(in) :: method
      character(len=:), allocatable :: out
      character(len=:), allocatable :: m
      integer :: i, c
      m = trim(adjustl(method))
      do i = 1, len(m)
         c = iachar(m(i:i))
         if (c >= iachar('A') .and. c <= iachar('Z')) m(i:i) = achar(c + 32)
      end do
      select case (m)
      case ('bfgs', 'nlminb')
         out = method_bfgs
      case ('nelder-mead', 'nelder_mead', 'nm')
         out = method_nelder_mead
      case ('cg')
         out = method_cg
      case ('l-bfgs-b', 'lbfgsb')
         out = method_lbfgsb
      case ('sann')
         out = method_sann
      case default
         out = method_bfgs
      end select
   end function canonical_inner_method

   logical function is_positive_definite(a) result(ok)
      real(dp), intent(in) :: a(:, :)
      real(dp), allocatable :: l(:, :)
      real(dp) :: s
      integer :: n, i, j, k
      n = size(a,1)
      if (size(a,2) /= n) then
         ok = .false.
         return
      end if
      allocate(l(n,n))
      l = 0.0_dp
      do i = 1, n
         do j = 1, i
            s = a(i,j)
            do k = 1, j-1
               s = s - l(i,k)*l(j,k)
            end do
            if (i == j) then
               if (.not. ieee_is_finite(s) .or. s <= 1.0e-12_dp) then
                  ok = .false.
                  return
               end if
               l(i,j) = sqrt(s)
            else
               l(i,j) = s/l(j,j)
            end if
         end do
      end do
      ok = .true.
   end function is_positive_definite

   subroutine trace_aug(i, par, obj, h, e)
      integer, intent(in) :: i
      real(dp), intent(in) :: par(:), obj, h(:), e(:)
      write(*,'(a,i0)') 'Outer iteration: ', i
      if (size(h) > 0) write(*,'(a,es13.5)') 'Min(hin): ', minval(h)
      if (size(e) > 0) write(*,'(a,es13.5)') 'Max(abs(heq)): ', maxval(abs(e))
      write(*,'(a,*(es13.5,1x))') 'par: ', par
      write(*,'(a,es13.5)') 'fval: ', obj
   end subroutine trace_aug

end module alabama
