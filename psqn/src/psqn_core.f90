! SPDX-License-Identifier: Apache-2.0
module psqn_core
  use psqn_types, only : dp, psqn_element_spec, psqn_element_eval, psqn_objective_eval, &
                         psqn_options, psqn_bfgs_options, psqn_auglag_options, psqn_info, psqn_auglag_info, &
                         psqn_converged, psqn_max_it_reached, psqn_conjugate_gradient_failed, &
                         psqn_line_search_failed, psqn_pre_none, psqn_pre_diag, &
                         psqn_pre_cholesky, psqn_pre_block
  use psqn_linalg, only : packed_size, packed_index, sym_packed_matvec, &
                          sym_packed_matvec_indexed, rank_one_update_packed, &
                          packed_to_dense, cholesky_factor_regularized, cholesky_solve, kahan_add
  use psqn_linesearch, only : wolfe_line_search
  use psqn_bfgs, only : psqn_bfgs_optimize
  use psqn_richardson, only : richardson_vector_derivative
  implicit none
  private

  public :: psqn_make_structured_specs
  public :: psqn_optimize_generic, psqn_optimize_structured, psqn_optimize_private_structured
  public :: psqn_aug_lagrang_generic, psqn_aug_lagrang_structured
  public :: psqn_generic_hess, psqn_structured_hess

  type :: element_state
    integer, allocatable :: idx(:)
    real(dp), allocatable :: b(:)
    real(dp), allocatable :: gr(:), gr_old(:), x_old(:), x_new(:)
    logical :: first_call = .true.
  end type element_state

contains

  function psqn_make_structured_specs(global_dim, private_dims) result(specs)
    integer, intent(in) :: global_dim
    integer, intent(in) :: private_dims(:)
    type(psqn_element_spec), allocatable :: specs(:)
    integer :: i, j, pstart, nloc

    if (global_dim < 0 .or. any(private_dims < 0)) error stop "psqn: negative dimension"
    allocate(specs(size(private_dims)))
    pstart = global_dim + 1
    do i = 1, size(private_dims)
      nloc = global_dim + private_dims(i)
      allocate(specs(i)%idx(nloc))
      if (global_dim > 0) specs(i)%idx(1:global_dim) = [(j, j=1,global_dim)]
      if (private_dims(i) > 0) then
        specs(i)%idx(global_dim+1:nloc) = [(j, j=pstart,pstart+private_dims(i)-1)]
      end if
      pstart = pstart + private_dims(i)
    end do
  end function psqn_make_structured_specs

  subroutine psqn_optimize_generic(x, specs, evaluator, result, options, masked, hess)
    real(dp), intent(inout) :: x(:)
    type(psqn_element_spec), intent(in) :: specs(:)
    procedure(psqn_element_eval) :: evaluator
    type(psqn_info), intent(out) :: result
    type(psqn_options), intent(in), optional :: options
    integer, intent(in), optional :: masked(:)
    real(dp), intent(out), optional :: hess(:,:)

    type(psqn_options) :: opt
    type(element_state), allocatable :: states(:)
    logical, allocatable :: mask(:)

    opt = psqn_options()
    if (present(options)) opt = options
    call validate_specs(size(x), specs)
    call init_states(specs, states)
    call make_mask(size(x), masked, mask)
    call optimize_core(x, states, evaluator, opt, mask, result, 0, &
                       constraint_evaluator=null_constraint_evaluator)
    if (present(hess)) call assemble_hess(size(x), states, hess)
  end subroutine psqn_optimize_generic

  subroutine psqn_optimize_structured(x, global_dim, private_dims, evaluator, result, options, masked, hess)
    real(dp), intent(inout) :: x(:)
    integer, intent(in) :: global_dim
    integer, intent(in) :: private_dims(:)
    procedure(psqn_element_eval) :: evaluator
    type(psqn_info), intent(out) :: result
    type(psqn_options), intent(in), optional :: options
    integer, intent(in), optional :: masked(:)
    real(dp), intent(out), optional :: hess(:,:)

    type(psqn_element_spec), allocatable :: specs(:)
    type(element_state), allocatable :: states(:)
    type(psqn_options) :: opt
    logical, allocatable :: mask(:)

    if (size(x) /= global_dim + sum(private_dims)) error stop "psqn structured: wrong x size"
    specs = psqn_make_structured_specs(global_dim, private_dims)
    opt = psqn_options()
    if (present(options)) opt = options
    call init_states(specs, states)
    call make_mask(size(x), masked, mask)
    call optimize_core(x, states, evaluator, opt, mask, result, global_dim, &
                       constraint_evaluator=null_constraint_evaluator)
    if (present(hess)) call assemble_hess(size(x), states, hess)
  end subroutine psqn_optimize_structured

  subroutine psqn_optimize_private_structured(x, global_dim, private_dims, evaluator, value, &
                                              rel_eps, max_it, c1, c2, gr_tol)
    real(dp), intent(inout) :: x(:)
    integer, intent(in) :: global_dim
    integer, intent(in) :: private_dims(:)
    procedure(psqn_element_eval) :: evaluator
    real(dp), intent(out) :: value
    real(dp), intent(in), optional :: rel_eps, c1, c2, gr_tol
    integer, intent(in), optional :: max_it

    type(psqn_bfgs_options) :: bopt
    type(psqn_info) :: bres
    real(dp), allocatable :: private_x(:), local_x(:)
    integer :: ie, pstart, npriv

    if (size(x) /= global_dim + sum(private_dims)) error stop "psqn private: wrong x size"
    bopt = psqn_bfgs_options()
    if (present(rel_eps)) bopt%rel_eps = rel_eps
    if (present(max_it)) bopt%max_it = max_it
    if (present(c1)) bopt%c1 = c1
    if (present(c2)) bopt%c2 = c2
    if (present(gr_tol)) bopt%gr_tol = gr_tol

    value = 0.0_dp
    pstart = global_dim + 1
    do ie = 1, size(private_dims)
      npriv = private_dims(ie)
      if (npriv > 0) then
        private_x = x(pstart:pstart+npriv-1)
        allocate(local_x(global_dim+npriv))
        if (global_dim > 0) local_x(1:global_dim) = x(1:global_dim)
        call psqn_bfgs_optimize(private_x, private_objective, bres, bopt)
        x(pstart:pstart+npriv-1) = private_x
        value = value + bres%value
        deallocate(local_x, private_x)
      end if
      pstart = pstart + npriv
    end do

  contains
    subroutine private_objective(v, f, g, comp_grad)
      real(dp), intent(in) :: v(:)
      real(dp), intent(out) :: f, g(:)
      logical, intent(in) :: comp_grad
      real(dp), allocatable :: glocal(:)
      if (npriv < 1) then
        f = 0.0_dp
        if (comp_grad) g = 0.0_dp
        return
      end if
      local_x(global_dim+1:) = v
      allocate(glocal(global_dim+npriv))
      glocal = 0.0_dp
      call evaluator(ie, local_x, f, glocal, comp_grad)
      if (comp_grad) g = glocal(global_dim+1:)
    end subroutine private_objective
  end subroutine psqn_optimize_private_structured

  subroutine psqn_aug_lagrang_generic(x, specs, evaluator, constraint_specs, constraint_evaluator, &
                                      multipliers, result, options, aug_options, masked, hess)
    real(dp), intent(inout) :: x(:)
    type(psqn_element_spec), intent(in) :: specs(:), constraint_specs(:)
    procedure(psqn_element_eval) :: evaluator, constraint_evaluator
    real(dp), intent(inout) :: multipliers(:)
    type(psqn_auglag_info), intent(out) :: result
    type(psqn_options), intent(in), optional :: options
    type(psqn_auglag_options), intent(in), optional :: aug_options
    integer, intent(in), optional :: masked(:)
    real(dp), intent(out), optional :: hess(:,:)

    call auglag_driver(x, specs, evaluator, constraint_specs, constraint_evaluator, multipliers, &
                       result, options, aug_options, masked, hess, 0)
  end subroutine psqn_aug_lagrang_generic

  subroutine psqn_aug_lagrang_structured(x, global_dim, private_dims, evaluator, constraint_specs, &
                                         constraint_evaluator, multipliers, result, options, &
                                         aug_options, masked, hess)
    real(dp), intent(inout) :: x(:)
    integer, intent(in) :: global_dim
    integer, intent(in) :: private_dims(:)
    procedure(psqn_element_eval) :: evaluator, constraint_evaluator
    type(psqn_element_spec), intent(in) :: constraint_specs(:)
    real(dp), intent(inout) :: multipliers(:)
    type(psqn_auglag_info), intent(out) :: result
    type(psqn_options), intent(in), optional :: options
    type(psqn_auglag_options), intent(in), optional :: aug_options
    integer, intent(in), optional :: masked(:)
    real(dp), intent(out), optional :: hess(:,:)
    type(psqn_element_spec), allocatable :: specs(:)

    if (size(x) /= global_dim + sum(private_dims)) error stop "psqn structured: wrong x size"
    specs = psqn_make_structured_specs(global_dim, private_dims)
    call auglag_driver(x, specs, evaluator, constraint_specs, constraint_evaluator, multipliers, &
                       result, options, aug_options, masked, hess, global_dim)
  end subroutine psqn_aug_lagrang_structured

  subroutine auglag_driver(x, specs, evaluator, constraint_specs, constraint_evaluator, multipliers, &
                           result, options, aug_options, masked, hess, global_dim)
    real(dp), intent(inout) :: x(:)
    type(psqn_element_spec), intent(in) :: specs(:), constraint_specs(:)
    procedure(psqn_element_eval) :: evaluator, constraint_evaluator
    real(dp), intent(inout) :: multipliers(:)
    type(psqn_auglag_info), intent(out) :: result
    type(psqn_options), intent(in), optional :: options
    type(psqn_auglag_options), intent(in), optional :: aug_options
    integer, intent(in), optional :: masked(:)
    real(dp), intent(out), optional :: hess(:,:)
    integer, intent(in) :: global_dim

    type(psqn_options) :: opt
    type(psqn_auglag_options) :: aopt
    type(psqn_info) :: inner
    type(element_state), allocatable :: states(:), cstates(:)
    logical, allocatable :: mask(:)
    real(dp) :: penalty, norm_v
    real(dp), allocatable :: cvals(:), gdummy(:)
    integer :: outer, i

    opt = psqn_options()
    if (present(options)) opt = options
    aopt = psqn_auglag_options()
    if (present(aug_options)) aopt = aug_options
    if (aopt%tau < 1.0_dp) error stop "psqn auglag: tau < 1"
    if (aopt%penalty_start <= 0.0_dp) error stop "psqn auglag: penalty_start <= 0"
    if (size(multipliers) /= size(constraint_specs)) error stop "psqn auglag: wrong multiplier count"

    call validate_specs(size(x), specs)
    call validate_specs(size(x), constraint_specs)
    call init_states(specs, states)
    call init_states(constraint_specs, cstates)
    call make_mask(size(x), masked, mask)
    allocate(cvals(size(cstates)))

    result = psqn_auglag_info()
    penalty = aopt%penalty_start
    do outer = 1, aopt%max_it_outer
      call optimize_core(x, states, evaluator, opt, mask, inner, global_dim, cstates, &
                         constraint_evaluator, multipliers, penalty)
      result%n_eval = result%n_eval + inner%n_eval
      result%n_grad = result%n_grad + inner%n_grad
      result%n_cg = result%n_cg + inner%n_cg
      result%value = inner%value
      result%info = inner%info
      result%n_aug_lagrang = outer
      result%penalty = penalty
      if (inner%info /= psqn_converged) exit

      do i = 1, size(cstates)
        allocate(gdummy(size(cstates(i)%idx)))
        call constraint_evaluator(i, x(cstates(i)%idx), cvals(i), gdummy, .false.)
        deallocate(gdummy)
      end do
      norm_v = sqrt(dot_product(cvals, cvals))
      if (norm_v < aopt%violations_norm_thresh) then
        result%info = psqn_converged
        exit
      end if
      multipliers = multipliers - penalty * cvals
      penalty = penalty * aopt%tau
      result%penalty = penalty
    end do

    if (present(hess)) call assemble_hess(size(x), states, hess)
  end subroutine auglag_driver

  subroutine psqn_generic_hess(x, specs, evaluator, hess, eps, scale, tol, order)
    real(dp), intent(in) :: x(:)
    type(psqn_element_spec), intent(in) :: specs(:)
    procedure(psqn_element_eval) :: evaluator
    real(dp), intent(out) :: hess(:,:)
    real(dp), intent(in), optional :: eps, scale, tol
    integer, intent(in), optional :: order
    type(element_state), allocatable :: states(:)
    real(dp) :: eps_use, scale_use, tol_use
    integer :: order_use

    call validate_specs(size(x), specs)
    call init_states(specs, states)
    eps_use = 1.0e-3_dp; if (present(eps)) eps_use = eps
    scale_use = 2.0_dp; if (present(scale)) scale_use = scale
    tol_use = 1.0e-9_dp; if (present(tol)) tol_use = tol
    order_use = 6; if (present(order)) order_use = order
    call compute_true_hess_states(x, states, evaluator, eps_use, scale_use, tol_use, order_use)
    call assemble_hess(size(x), states, hess)
  end subroutine psqn_generic_hess

  subroutine psqn_structured_hess(x, global_dim, private_dims, evaluator, hess, eps, scale, tol, order)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: global_dim
    integer, intent(in) :: private_dims(:)
    procedure(psqn_element_eval) :: evaluator
    real(dp), intent(out) :: hess(:,:)
    real(dp), intent(in), optional :: eps, scale, tol
    integer, intent(in), optional :: order
    type(psqn_element_spec), allocatable :: specs(:)

    if (size(x) /= global_dim + sum(private_dims)) error stop "psqn structured hess: wrong x size"
    specs = psqn_make_structured_specs(global_dim, private_dims)
    call psqn_generic_hess(x, specs, evaluator, hess, eps, scale, tol, order)
  end subroutine psqn_structured_hess

  subroutine optimize_core(x, states, evaluator, opt, mask, result, global_dim, cstates, &
                           constraint_evaluator, multipliers, penalty)
    real(dp), intent(inout) :: x(:)
    type(element_state), intent(inout) :: states(:)
    procedure(psqn_element_eval) :: evaluator
    type(psqn_options), intent(in) :: opt
    logical, intent(in) :: mask(:)
    type(psqn_info), intent(out) :: result
    integer, intent(in) :: global_dim
    type(element_state), intent(inout), optional :: cstates(:)
    procedure(psqn_element_eval) :: constraint_evaluator
    real(dp), intent(in), optional :: multipliers(:), penalty

    real(dp), allocatable :: gr(:), direction(:)
    real(dp) :: fval, fval_old
    integer :: iter, n_line_search_fail, n_cg_this
    logical :: ls_ok, cg_ok, has_constraints

    if (opt%c1 < 0.0_dp .or. opt%c1 >= opt%c2 .or. opt%c2 >= 1.0_dp) &
      error stop "psqn: invalid Wolfe constants"
    has_constraints = present(cstates)
    if (has_constraints .neqv. present(multipliers)) error stop "psqn: incomplete constraint arguments"
    if (has_constraints .neqv. present(penalty)) error stop "psqn: incomplete constraint arguments"

    call reset_states(states)
    if (has_constraints) call reset_states(cstates)
    allocate(gr(size(x)), direction(size(x)))
    result = psqn_info()
    call eval_total(x, fval, gr, .true.)
    result%n_grad = 1
    call record_states(states)
    if (has_constraints) call record_states(cstates)
    result%info = psqn_max_it_reached
    n_line_search_fail = 0

    do iter = 1, opt%max_it
      fval_old = fval
      call conjugate_gradient(gr, direction, cg_ok, n_cg_this)
      result%n_cg = result%n_cg + n_cg_this
      if (.not. cg_ok) then
        result%info = psqn_conjugate_gradient_failed
        exit
      end if
      direction = -direction

      ls_ok = wolfe_line_search(eval_for_linesearch, fval_old, x, gr, direction, fval, &
                                opt%c1, opt%c2, opt%strong_wolfe, result%n_eval, result%n_grad, opt%trace)
      if (.not. ls_ok) then
        result%info = psqn_line_search_failed
        n_line_search_fail = n_line_search_fail + 1
        if (n_line_search_fail > 2) exit
      else
        n_line_search_fail = 0
      end if

      if (opt%trace > 0) then
        write(*,'(a,i0,2(a,es16.8),a,i0,a,l1)') 'psqn iter=', iter, ' f=', fval, &
          ' |g|=', sqrt(dot_product(gr,gr)), ' cg=', n_cg_this, ' line_search=', ls_ok
      end if

      if (n_line_search_fail < 1 .and. &
          abs(fval-fval_old) < opt%rel_eps*(abs(fval_old)+opt%rel_eps) .and. &
          (opt%gr_tol <= 0.0_dp .or. dot_product(gr,gr) < opt%gr_tol**2)) then
        result%info = psqn_converged
        exit
      end if

      if (n_line_search_fail < 1) then
        call update_states(states, opt%use_bfgs)
        if (has_constraints) call update_states(cstates, opt%use_bfgs)
      else
        call reset_states(states)
        call record_states(states)
        if (has_constraints) then
          call reset_states(cstates)
          call record_states(cstates)
        end if
      end if
    end do
    result%value = fval

  contains

    subroutine eval_for_linesearch(v, f, g, comp_grad)
      real(dp), intent(in) :: v(:)
      real(dp), intent(out) :: f, g(:)
      logical, intent(in) :: comp_grad
      call eval_total(v, f, g, comp_grad)
    end subroutine eval_for_linesearch

    subroutine eval_total(v, f, g, comp_grad)
      real(dp), intent(in) :: v(:)
      real(dp), intent(out) :: f
      real(dp), intent(out) :: g(:)
      logical, intent(in) :: comp_grad
      real(dp), allocatable :: gloc(:), comp(:)
      real(dp) :: fi, fcomp, cval, factor
      integer :: i, j, idxj

      f = 0.0_dp
      fcomp = 0.0_dp
      if (comp_grad) then
        g = 0.0_dp
        allocate(comp(size(g)))
        comp = 0.0_dp
      end if

      do i = 1, size(states)
        states(i)%x_new = v(states(i)%idx)
        allocate(gloc(size(states(i)%idx)))
        gloc = 0.0_dp
        call evaluator(i, states(i)%x_new, fi, gloc, comp_grad)
        call kahan_add(f, fcomp, fi)
        if (comp_grad) then
          states(i)%gr = gloc
          do j = 1, size(states(i)%idx)
            idxj = states(i)%idx(j)
            if (mask(idxj)) states(i)%gr(j) = 0.0_dp
            call kahan_add(g(idxj), comp(idxj), states(i)%gr(j))
          end do
        end if
        deallocate(gloc)
      end do
      f = f - fcomp

      if (has_constraints) then
        do i = 1, size(cstates)
          cstates(i)%x_new = v(cstates(i)%idx)
          allocate(gloc(size(cstates(i)%idx)))
          gloc = 0.0_dp
          call constraint_evaluator(i, cstates(i)%x_new, cval, gloc, comp_grad)
          f = f - multipliers(i) * cval + 0.5_dp * penalty * cval * cval
          if (comp_grad) then
            factor = penalty * cval - multipliers(i)
            cstates(i)%gr = gloc * factor
            do j = 1, size(cstates(i)%idx)
              idxj = cstates(i)%idx(j)
              if (mask(idxj)) cstates(i)%gr(j) = 0.0_dp
              call kahan_add(g(idxj), comp(idxj), cstates(i)%gr(j))
            end do
          end if
          deallocate(gloc)
        end do
      end if
    end subroutine eval_total

    subroutine conjugate_gradient(rhs, sol, ok, niter)
      real(dp), intent(in) :: rhs(:)
      real(dp), intent(out) :: sol(:)
      logical, intent(out) :: ok
      integer, intent(out) :: niter
      real(dp), allocatable :: r(:), p(:), bp(:), z(:), diag(:)
      real(dp) :: rhs_norm, tol_use, eps_cg, old_rz, rz, pbp, alpha, beta, test_norm
      integer :: i, max_cg_use

      allocate(r(size(rhs)), p(size(rhs)), bp(size(rhs)), z(size(rhs)), diag(size(rhs)))
      rhs_norm = sqrt(abs(dot_product(rhs, rhs)))
      tol_use = min(opt%cg_tol, sqrt(rhs_norm))
      eps_cg = tol_use * rhs_norm
      max_cg_use = opt%max_cg
      if (max_cg_use < 1) max_cg_use = size(rhs)

      sol = 0.0_dp
      r = -rhs
      select case (opt%pre_method)
      case (psqn_pre_none)
        p = rhs
        z = r
      case (psqn_pre_diag)
        call total_diag(diag)
        where (abs(diag) < sqrt(tiny(1.0_dp))) diag = sign(sqrt(tiny(1.0_dp)), diag + tiny(1.0_dp))
        z = r / diag
        p = -z
      case (psqn_pre_cholesky)
        call apply_full_preconditioner(r, z)
        p = -z
      case (psqn_pre_block)
        if (global_dim > 0) then
          call apply_block_preconditioner(r, z)
        else
          call total_diag(diag)
          where (abs(diag) < sqrt(tiny(1.0_dp))) diag = sign(sqrt(tiny(1.0_dp)), diag + tiny(1.0_dp))
          z = r / diag
        end if
        p = -z
      case default
        error stop "psqn: invalid preconditioner"
      end select

      if (opt%pre_method == psqn_pre_none) then
        old_rz = dot_product(r, r)
      else
        old_rz = dot_product(r, z)
      end if
      ok = .true.
      niter = 0

      do i = 1, max_cg_use
        niter = niter + 1
        bp = 0.0_dp
        call total_hvp(p, bp)
        pbp = dot_product(p, bp)
        if (pbp <= 0.0_dp) then
          if (i == 1) sol = rhs
          exit
        end if
        alpha = old_rz / pbp
        sol = sol + alpha * p
        r = r + alpha * bp

        select case (opt%pre_method)
        case (psqn_pre_none)
          rz = dot_product(r, r)
          test_norm = sqrt(abs(rz))
        case (psqn_pre_diag)
          z = r / diag
          rz = dot_product(r, z)
          test_norm = sqrt(abs(dot_product(r,r)))
        case (psqn_pre_cholesky)
          call apply_full_preconditioner(r, z)
          rz = dot_product(r, z)
          test_norm = sqrt(abs(dot_product(r,r)))
        case (psqn_pre_block)
          if (global_dim > 0) then
            call apply_block_preconditioner(r, z)
          else
            z = r / diag
          end if
          rz = dot_product(r, z)
          test_norm = sqrt(abs(dot_product(r,r)))
        end select

        if (opt%trace >= 3) then
          write(*,'(a,i0,2(a,es14.6))') '  cg iter=', i, ' residual=', test_norm, ' threshold=', eps_cg
        end if
        if (test_norm < eps_cg) exit
        if (abs(old_rz) <= tiny(1.0_dp)) exit
        beta = rz / old_rz
        old_rz = rz
        if (opt%pre_method == psqn_pre_none) then
          p = beta * p - r
        else
          p = beta * p - z
        end if
      end do
    end subroutine conjugate_gradient

    subroutine total_hvp(v, y)
      real(dp), intent(in) :: v(:)
      real(dp), intent(inout) :: y(:)
      real(dp), allocatable :: yl(:)
      integer :: i, j

      do i = 1, size(states)
        allocate(yl(size(states(i)%idx)))
        call sym_packed_matvec_indexed(states(i)%b, v, states(i)%idx, yl)
        do j = 1, size(states(i)%idx)
          y(states(i)%idx(j)) = y(states(i)%idx(j)) + yl(j)
        end do
        deallocate(yl)
      end do
      if (has_constraints) then
        do i = 1, size(cstates)
          allocate(yl(size(cstates(i)%idx)))
          call sym_packed_matvec_indexed(cstates(i)%b, v, cstates(i)%idx, yl)
          do j = 1, size(cstates(i)%idx)
            y(cstates(i)%idx(j)) = y(cstates(i)%idx(j)) + yl(j)
          end do
          deallocate(yl)
        end do
      end if
    end subroutine total_hvp

    subroutine total_diag(d)
      real(dp), intent(out) :: d(:)
      integer :: i, j
      d = 0.0_dp
      do i = 1, size(states)
        do j = 1, size(states(i)%idx)
          d(states(i)%idx(j)) = d(states(i)%idx(j)) + states(i)%b(packed_index(j,j))
        end do
      end do
      if (has_constraints) then
        do i = 1, size(cstates)
          do j = 1, size(cstates(i)%idx)
            d(cstates(i)%idx(j)) = d(cstates(i)%idx(j)) + cstates(i)%b(packed_index(j,j))
          end do
        end do
      end if
    end subroutine total_diag

    subroutine apply_full_preconditioner(rhs_in, out)
      real(dp), intent(in) :: rhs_in(:)
      real(dp), intent(out) :: out(:)
      real(dp), allocatable :: hfull(:,:), l(:,:)
      logical :: chol_ok
      allocate(hfull(size(rhs_in),size(rhs_in)), l(size(rhs_in),size(rhs_in)))
      call assemble_total_hess(hfull)
      call cholesky_factor_regularized(hfull, l, chol_ok)
      if (.not. chol_ok) then
        out = rhs_in
      else
        call cholesky_solve(l, rhs_in, out)
      end if
    end subroutine apply_full_preconditioner

    subroutine apply_block_preconditioner(rhs_in, out)
      real(dp), intent(in) :: rhs_in(:)
      real(dp), intent(out) :: out(:)
      real(dp), allocatable :: block(:,:), l(:,:), bloc(:), sol(:)
      integer :: i, j, k, npriv
      logical :: chol_ok

      out = rhs_in
      allocate(block(global_dim,global_dim), l(global_dim,global_dim), bloc(global_dim), sol(global_dim))
      block = 0.0_dp
      do i = 1, size(states)
        do j = 1, global_dim
          do k = 1, global_dim
            block(j,k) = block(j,k) + states(i)%b(packed_index(j,k))
          end do
        end do
      end do
      call cholesky_factor_regularized(block, l, chol_ok)
      bloc = rhs_in(1:global_dim)
      call cholesky_solve(l, bloc, sol)
      out(1:global_dim) = sol
      deallocate(block, l, bloc, sol)

      do i = 1, size(states)
        npriv = size(states(i)%idx) - global_dim
        if (npriv <= 0) cycle
        allocate(block(npriv,npriv), l(npriv,npriv), bloc(npriv), sol(npriv))
        do j = 1, npriv
          do k = 1, npriv
            block(j,k) = states(i)%b(packed_index(global_dim+j, global_dim+k))
          end do
        end do
        call cholesky_factor_regularized(block, l, chol_ok)
        bloc = rhs_in(states(i)%idx(global_dim+1:))
        call cholesky_solve(l, bloc, sol)
        out(states(i)%idx(global_dim+1:)) = sol
        deallocate(block, l, bloc, sol)
      end do
    end subroutine apply_block_preconditioner

    subroutine assemble_total_hess(h)
      real(dp), intent(out) :: h(:,:)
      integer :: i, j, k, ii, jj
      h = 0.0_dp
      do i = 1, size(states)
        do j = 1, size(states(i)%idx)
          jj = states(i)%idx(j)
          do k = 1, j
            ii = states(i)%idx(k)
            h(ii,jj) = h(ii,jj) + states(i)%b(packed_index(k,j))
            if (ii /= jj) h(jj,ii) = h(jj,ii) + states(i)%b(packed_index(k,j))
          end do
        end do
      end do
      if (has_constraints) then
        do i = 1, size(cstates)
          do j = 1, size(cstates(i)%idx)
            jj = cstates(i)%idx(j)
            do k = 1, j
              ii = cstates(i)%idx(k)
              h(ii,jj) = h(ii,jj) + cstates(i)%b(packed_index(k,j))
              if (ii /= jj) h(jj,ii) = h(jj,ii) + cstates(i)%b(packed_index(k,j))
            end do
          end do
        end do
      end if
    end subroutine assemble_total_hess

  end subroutine optimize_core

  subroutine null_constraint_evaluator(i, x, f, g, comp_grad)
    integer, intent(in) :: i
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f
    real(dp), intent(out) :: g(:)
    logical, intent(in) :: comp_grad

    ! Required only so optimize_core always receives a procedure with an
    ! explicit interface. It is never called when constraints are absent.
    f = 0.0_dp + 0.0_dp * real(i + size(x), dp)
    if (comp_grad) g = 0.0_dp
  end subroutine null_constraint_evaluator

  subroutine init_states(specs, states)
    type(psqn_element_spec), intent(in) :: specs(:)
    type(element_state), allocatable, intent(out) :: states(:)
    integer :: i, n

    allocate(states(size(specs)))
    do i = 1, size(specs)
      n = size(specs(i)%idx)
      states(i)%idx = specs(i)%idx
      allocate(states(i)%b(packed_size(n)), states(i)%gr(n), states(i)%gr_old(n), &
               states(i)%x_old(n), states(i)%x_new(n))
    end do
    call reset_states(states)
  end subroutine init_states

  subroutine validate_specs(n_par, specs)
    integer, intent(in) :: n_par
    type(psqn_element_spec), intent(in) :: specs(:)
    integer :: i
    if (size(specs) < 1) error stop "psqn: no element functions"
    do i = 1, size(specs)
      if (.not. allocated(specs(i)%idx)) error stop "psqn: unallocated element indices"
      if (size(specs(i)%idx) < 1) error stop "psqn: empty element"
      if (any(specs(i)%idx < 1) .or. any(specs(i)%idx > n_par)) error stop "psqn: index out of range"
    end do
  end subroutine validate_specs

  subroutine make_mask(n_par, masked, mask)
    integer, intent(in) :: n_par
    integer, intent(in), optional :: masked(:)
    logical, allocatable, intent(out) :: mask(:)
    integer :: i
    allocate(mask(n_par))
    mask = .false.
    if (present(masked)) then
      do i = 1, size(masked)
        if (masked(i) < 1 .or. masked(i) > n_par) error stop "psqn: masked index out of range"
        mask(masked(i)) = .true.
      end do
    end if
  end subroutine make_mask

  subroutine reset_states(states)
    type(element_state), intent(inout) :: states(:)
    integer :: i, j, n
    do i = 1, size(states)
      n = size(states(i)%idx)
      states(i)%b = 0.0_dp
      do j = 1, n
        states(i)%b(packed_index(j,j)) = 1.0_dp
      end do
      states(i)%first_call = .true.
    end do
  end subroutine reset_states

  subroutine record_states(states)
    type(element_state), intent(inout) :: states(:)
    integer :: i
    do i = 1, size(states)
      states(i)%x_old = states(i)%x_new
      states(i)%gr_old = states(i)%gr
    end do
  end subroutine record_states

  subroutine update_states(states, use_bfgs)
    type(element_state), intent(inout) :: states(:)
    logical, intent(in) :: use_bfgs
    integer :: i
    do i = 1, size(states)
      call update_one(states(i), use_bfgs)
    end do
  end subroutine update_states

  subroutine update_one(st, use_bfgs)
    type(element_state), intent(inout) :: st
    logical, intent(in) :: use_bfgs
    real(dp), allocatable :: s(:), y(:), bs(:), w(:)
    real(dp) :: sy, sbs, theta, sr, sw, snorm, wnorm, scale
    integer :: j, n
    logical :: all_unchanged

    n = size(st%idx)
    allocate(s(n), y(n), bs(n), w(n))
    s = st%x_new - st%x_old
    all_unchanged = .true.
    do j = 1, n
      if (abs(s(j)) > abs(st%x_new(j)) * epsilon(1.0_dp) * 100.0_dp) then
        all_unchanged = .false.
        exit
      end if
    end do
    if (all_unchanged) then
      call reset_one(st)
      call record_one(st)
      return
    end if

    y = st%gr - st%gr_old
    sy = dot_product(y, s)
    if (use_bfgs) then
      if (abs(sy) <= tiny(1.0_dp)) then
        call reset_one(st)
        call record_one(st)
        return
      end if
      if (st%first_call) then
        st%first_call = .false.
        scale = dot_product(y,y) / sy
        st%b = 0.0_dp
        do j = 1, n
          st%b(packed_index(j,j)) = scale
        end do
      end if

      bs = 0.0_dp
      call sym_packed_matvec(st%b, s, bs)
      sbs = dot_product(s, bs)
      if (sbs <= tiny(1.0_dp)) then
        call reset_one(st)
        call record_one(st)
        return
      end if
      call rank_one_update_packed(st%b, bs, -1.0_dp/sbs)
      if (sy < 0.2_dp * sbs) then
        theta = 0.8_dp * sbs / (sbs - sy)
        y = theta * y + (1.0_dp - theta) * bs
        sr = dot_product(y, s)
        if (abs(sr) <= tiny(1.0_dp)) then
          call reset_one(st)
        else
          call rank_one_update_packed(st%b, y, 1.0_dp/sr)
        end if
      else
        call rank_one_update_packed(st%b, y, 1.0_dp/sy)
      end if
    else
      if (abs(sy) <= tiny(1.0_dp)) then
        call reset_one(st)
        call record_one(st)
        return
      end if
      if (st%first_call) then
        st%first_call = .false.
        scale = dot_product(y,y) / sy
        st%b = 0.0_dp
        do j = 1, n
          st%b(packed_index(j,j)) = scale
        end do
      end if
      w = 0.0_dp
      call sym_packed_matvec(st%b, s, w)
      w = y - w
      sw = dot_product(s, w)
      snorm = sqrt(abs(dot_product(s,s)))
      wnorm = sqrt(abs(dot_product(w,w)))
      if (abs(sw) > 1.0e-8_dp * snorm * wnorm) call rank_one_update_packed(st%b, w, 1.0_dp/sw)
    end if
    call record_one(st)
  end subroutine update_one

  subroutine reset_one(st)
    type(element_state), intent(inout) :: st
    integer :: j
    st%b = 0.0_dp
    do j = 1, size(st%idx)
      st%b(packed_index(j,j)) = 1.0_dp
    end do
    st%first_call = .true.
  end subroutine reset_one

  subroutine record_one(st)
    type(element_state), intent(inout) :: st
    st%x_old = st%x_new
    st%gr_old = st%gr
  end subroutine record_one

  subroutine assemble_hess(n_par, states, hess)
    integer, intent(in) :: n_par
    type(element_state), intent(in) :: states(:)
    real(dp), intent(out) :: hess(:,:)
    integer :: i, j, k, ii, jj

    if (size(hess,1) /= n_par .or. size(hess,2) /= n_par) error stop "psqn: wrong Hessian shape"
    hess = 0.0_dp
    do i = 1, size(states)
      do j = 1, size(states(i)%idx)
        jj = states(i)%idx(j)
        do k = 1, j
          ii = states(i)%idx(k)
          hess(ii,jj) = hess(ii,jj) + states(i)%b(packed_index(k,j))
          if (ii /= jj) hess(jj,ii) = hess(jj,ii) + states(i)%b(packed_index(k,j))
        end do
      end do
    end do
  end subroutine assemble_hess

  subroutine compute_true_hess_states(x, states, evaluator, eps, scale, tol, order)
    real(dp), intent(in) :: x(:)
    type(element_state), intent(inout) :: states(:)
    procedure(psqn_element_eval) :: evaluator
    real(dp), intent(in) :: eps, scale, tol
    integer, intent(in) :: order
    integer :: ie, j, n
    real(dp), allocatable :: xloc(:), deriv(:)

    do ie = 1, size(states)
      n = size(states(ie)%idx)
      xloc = x(states(ie)%idx)
      do j = 1, n
        allocate(deriv(j))
        call richardson_vector_derivative(local_gradient_slice, xloc(j), deriv, eps, scale, tol, order)
        states(ie)%b(packed_index(1,j):packed_index(j,j)) = deriv
        deallocate(deriv)
      end do
    end do

  contains
    subroutine local_gradient_slice(xj, out)
      real(dp), intent(in) :: xj
      real(dp), intent(out) :: out(:)
      real(dp), allocatable :: xtmp(:), gtmp(:)
      real(dp) :: ftmp
      allocate(xtmp(size(xloc)), gtmp(size(xloc)))
      xtmp = xloc
      xtmp(j) = xj
      call evaluator(ie, xtmp, ftmp, gtmp, .true.)
      out = gtmp(1:size(out))
    end subroutine local_gradient_slice
  end subroutine compute_true_hess_states

end module psqn_core
