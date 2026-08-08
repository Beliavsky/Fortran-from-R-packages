! SPDX-License-Identifier: GPL-3.0-only
module bb_drivers
  use bb_kinds, only: dp
  use bb_interfaces, only: bb_scalar_fn, bb_gradient_fn, bb_vector_fn, bb_projection_fn
  use bb_types, only: spg_control, spg_result, sane_control, sane_result, &
    bboptim_control, bbsolve_control, multistart_result
  use bb_spg, only: spg, spg_box, spg_projected, spg_linear
  use bb_nonlinear, only: dfsane
  implicit none
  private

  public :: bboptim, bboptim_box, bboptim_projected, bboptim_linear
  public :: bbsolve, multistart_solve, multistart_optimize, multistart_optimize_box
  public :: multistart_optimize_projected, multistart_optimize_linear

contains

  function bboptim(par, fn, control, gr) result(best)
    real(dp), intent(in) :: par(:)
    procedure(bb_scalar_fn) :: fn
    type(bboptim_control), intent(in), optional :: control
    procedure(bb_gradient_fn), optional :: gr
    type(spg_result) :: best
    type(bboptim_control) :: ctrl
    type(spg_control) :: sc
    type(spg_result) :: cur
    integer :: i, j, total_feval, total_iter
    real(dp) :: best_metric, metric
    logical :: have_best

    ctrl = bboptim_control()
    if (present(control)) ctrl = control
    total_feval = 0
    total_iter = 0
    best_metric = huge(1.0_dp)
    have_best = .false.

    do j = 1, size(ctrl%m_values)
      do i = 1, size(ctrl%methods)
        call fill_spg_control(sc, ctrl, ctrl%methods(i), ctrl%m_values(j))
        if (present(gr)) then
          cur = spg(par, fn, sc, gr)
        else
          cur = spg(par, fn, sc)
        end if
        total_feval = total_feval + cur%feval
        total_iter = total_iter + cur%iter
        cur%method = sc%method
        cur%m = sc%m
        metric = merge(-cur%value, cur%value, ctrl%maximize)
        if (.not. have_best .or. metric < best_metric) then
          best = cur
          best_metric = metric
          have_best = .true.
        end if
        if (cur%succeeded()) exit
      end do
      if (have_best) then
        if (best%succeeded()) exit
      end if
    end do
    best%feval = total_feval
    best%iter = total_iter
  end function bboptim

  function bboptim_box(par, fn, lower, upper, control, gr) result(best)
    real(dp), intent(in) :: par(:), lower(:), upper(:)
    procedure(bb_scalar_fn) :: fn
    type(bboptim_control), intent(in), optional :: control
    procedure(bb_gradient_fn), optional :: gr
    type(spg_result) :: best
    type(bboptim_control) :: ctrl
    type(spg_control) :: sc
    type(spg_result) :: cur
    integer :: i, j, total_feval, total_iter
    real(dp) :: best_metric, metric
    logical :: have_best

    ctrl = bboptim_control()
    if (present(control)) ctrl = control
    total_feval = 0
    total_iter = 0
    best_metric = huge(1.0_dp)
    have_best = .false.
    do j = 1, size(ctrl%m_values)
      do i = 1, size(ctrl%methods)
        call fill_spg_control(sc, ctrl, ctrl%methods(i), ctrl%m_values(j))
        if (present(gr)) then
          cur = spg_box(par, fn, lower, upper, sc, gr)
        else
          cur = spg_box(par, fn, lower, upper, sc)
        end if
        total_feval = total_feval + cur%feval
        total_iter = total_iter + cur%iter
        cur%method = sc%method
        cur%m = sc%m
        metric = merge(-cur%value, cur%value, ctrl%maximize)
        if (.not. have_best .or. metric < best_metric) then
          best = cur
          best_metric = metric
          have_best = .true.
        end if
        if (cur%succeeded()) exit
      end do
      if (have_best) then
        if (best%succeeded()) exit
      end if
    end do
    best%feval = total_feval
    best%iter = total_iter
  end function bboptim_box

  function bboptim_projected(par, fn, project, control, gr) result(best)
    real(dp), intent(in) :: par(:)
    procedure(bb_scalar_fn) :: fn
    procedure(bb_projection_fn) :: project
    type(bboptim_control), intent(in), optional :: control
    procedure(bb_gradient_fn), optional :: gr
    type(spg_result) :: best
    type(bboptim_control) :: ctrl
    type(spg_control) :: sc
    type(spg_result) :: cur
    integer :: i, j, total_feval, total_iter
    real(dp) :: best_metric, metric
    logical :: have_best

    ctrl = bboptim_control()
    if (present(control)) ctrl = control
    total_feval = 0
    total_iter = 0
    best_metric = huge(1.0_dp)
    have_best = .false.
    do j = 1, size(ctrl%m_values)
      do i = 1, size(ctrl%methods)
        call fill_spg_control(sc, ctrl, ctrl%methods(i), ctrl%m_values(j))
        if (present(gr)) then
          cur = spg_projected(par, fn, project, sc, gr)
        else
          cur = spg_projected(par, fn, project, sc)
        end if
        total_feval = total_feval + cur%feval
        total_iter = total_iter + cur%iter
        cur%method = sc%method
        cur%m = sc%m
        metric = merge(-cur%value, cur%value, ctrl%maximize)
        if (.not. have_best .or. metric < best_metric) then
          best = cur
          best_metric = metric
          have_best = .true.
        end if
        if (cur%succeeded()) exit
      end do
      if (have_best) then
        if (best%succeeded()) exit
      end if
    end do
    best%feval = total_feval
    best%iter = total_iter
  end function bboptim_projected

  function bboptim_linear(par, fn, a, b, meq, control, gr, lower, upper) result(best)
    real(dp), intent(in) :: par(:), a(:, :), b(:)
    integer, intent(in) :: meq
    procedure(bb_scalar_fn) :: fn
    type(bboptim_control), intent(in), optional :: control
    procedure(bb_gradient_fn), optional :: gr
    real(dp), intent(in), optional :: lower(:), upper(:)
    type(spg_result) :: best
    type(bboptim_control) :: ctrl
    type(spg_control) :: sc
    type(spg_result) :: cur
    integer :: i, j, total_feval, total_iter
    real(dp) :: best_metric, metric
    logical :: have_best

    ctrl = bboptim_control()
    if (present(control)) ctrl = control
    total_feval = 0
    total_iter = 0
    best_metric = huge(1.0_dp)
    have_best = .false.
    do j = 1, size(ctrl%m_values)
      do i = 1, size(ctrl%methods)
        call fill_spg_control(sc, ctrl, ctrl%methods(i), ctrl%m_values(j))
        if (present(lower) .and. present(upper)) then
          if (present(gr)) then
            cur = spg_linear(par, fn, a, b, meq, sc, gr, lower, upper)
          else
            cur = spg_linear(par, fn, a, b, meq, sc, lower=lower, upper=upper)
          end if
        else if (present(lower)) then
          if (present(gr)) then
            cur = spg_linear(par, fn, a, b, meq, sc, gr=gr, lower=lower)
          else
            cur = spg_linear(par, fn, a, b, meq, sc, lower=lower)
          end if
        else if (present(upper)) then
          if (present(gr)) then
            cur = spg_linear(par, fn, a, b, meq, sc, gr=gr, upper=upper)
          else
            cur = spg_linear(par, fn, a, b, meq, sc, upper=upper)
          end if
        else
          if (present(gr)) then
            cur = spg_linear(par, fn, a, b, meq, sc, gr)
          else
            cur = spg_linear(par, fn, a, b, meq, sc)
          end if
        end if
        total_feval = total_feval + cur%feval
        total_iter = total_iter + cur%iter
        cur%method = sc%method
        cur%m = sc%m
        metric = merge(-cur%value, cur%value, ctrl%maximize)
        if (.not. have_best .or. metric < best_metric) then
          best = cur
          best_metric = metric
          have_best = .true.
        end if
        if (cur%succeeded()) exit
      end do
      if (have_best) then
        if (best%succeeded()) exit
      end if
    end do
    best%feval = total_feval
    best%iter = total_iter
  end function bboptim_linear

  function bbsolve(par, fn, control) result(best)
    real(dp), intent(in) :: par(:)
    procedure(bb_vector_fn) :: fn
    type(bbsolve_control), intent(in), optional :: control
    type(sane_result) :: best
    type(bbsolve_control) :: ctrl
    type(sane_control) :: sc
    type(sane_result) :: cur
    integer :: i, j, k, total_feval, total_iter, nnm
    logical :: nm_values(2), have_best
    real(dp) :: best_value

    ctrl = bbsolve_control()
    if (present(control)) ctrl = control
    if (ctrl%try_nm .and. size(par) > 1 .and. size(par) <= 20) then
      nm_values = [.true., .false.]
      nnm = 2
    else
      nm_values = [.false., .false.]
      nnm = 1
    end if
    total_feval = 0
    total_iter = 0
    best_value = huge(1.0_dp)
    have_best = .false.

    do k = 1, nnm
      do j = 1, size(ctrl%m_values)
        do i = 1, size(ctrl%methods)
          sc = sane_control(maxit=ctrl%maxit, m=ctrl%m_values(j), method=ctrl%methods(i), &
            tol=ctrl%tol, trace=ctrl%trace, triter=ctrl%triter, noimp=ctrl%noimp, &
            nm=nm_values(k), bfgs=.false.)
          cur = dfsane(par, fn, sc)
          total_feval = total_feval + cur%feval
          total_iter = total_iter + cur%iter
          cur%method = sc%method
          cur%m = sc%m
          cur%nm_used = sc%nm
          if (.not. have_best .or. cur%residual < best_value) then
            best = cur
            best_value = cur%residual
            have_best = .true.
          end if
          if (cur%succeeded()) exit
        end do
        if (have_best) then
          if (best%succeeded()) exit
        end if
      end do
      if (have_best) then
        if (best%succeeded()) exit
      end if
    end do
    best%feval = total_feval
    best%iter = total_iter
  end function bbsolve

  function multistart_solve(starts, fn, control) result(ans)
    real(dp), intent(in) :: starts(:, :)
    procedure(bb_vector_fn) :: fn
    type(bbsolve_control), intent(in), optional :: control
    type(multistart_result) :: ans
    type(sane_result) :: fit
    integer :: k

    allocate(ans%par(size(starts,1), size(starts,2)))
    allocate(ans%fvalue(size(starts,2)), ans%converged(size(starts,2)))
    ans%par = 0.0_dp
    ans%fvalue = huge(1.0_dp)
    ans%converged = .false.
    do k = 1, size(starts, 2)
      if (present(control)) then
        fit = bbsolve(starts(:,k), fn, control)
      else
        fit = bbsolve(starts(:,k), fn)
      end if
      ans%par(:,k) = fit%par
      ans%fvalue(k) = fit%residual
      ans%converged(k) = fit%succeeded()
    end do
  end function multistart_solve

  function multistart_optimize(starts, fn, control, gr) result(ans)
    real(dp), intent(in) :: starts(:, :)
    procedure(bb_scalar_fn) :: fn
    type(bboptim_control), intent(in), optional :: control
    procedure(bb_gradient_fn), optional :: gr
    type(multistart_result) :: ans
    type(spg_result) :: fit
    integer :: k

    allocate(ans%par(size(starts,1), size(starts,2)))
    allocate(ans%fvalue(size(starts,2)), ans%converged(size(starts,2)))
    do k = 1, size(starts, 2)
      if (present(control)) then
        if (present(gr)) then
          fit = bboptim(starts(:,k), fn, control, gr)
        else
          fit = bboptim(starts(:,k), fn, control)
        end if
      else
        if (present(gr)) then
          fit = bboptim(starts(:,k), fn, gr=gr)
        else
          fit = bboptim(starts(:,k), fn)
        end if
      end if
      ans%par(:,k) = fit%par
      ans%fvalue(k) = fit%value
      ans%converged(k) = fit%succeeded()
    end do
  end function multistart_optimize

  function multistart_optimize_box(starts, fn, lower, upper, control, gr) result(ans)
    real(dp), intent(in) :: starts(:, :), lower(:), upper(:)
    procedure(bb_scalar_fn) :: fn
    type(bboptim_control), intent(in), optional :: control
    procedure(bb_gradient_fn), optional :: gr
    type(multistart_result) :: ans
    type(spg_result) :: fit
    integer :: k

    allocate(ans%par(size(starts,1), size(starts,2)))
    allocate(ans%fvalue(size(starts,2)), ans%converged(size(starts,2)))
    do k = 1, size(starts, 2)
      if (present(control)) then
        if (present(gr)) then
          fit = bboptim_box(starts(:,k), fn, lower, upper, control, gr)
        else
          fit = bboptim_box(starts(:,k), fn, lower, upper, control)
        end if
      else
        if (present(gr)) then
          fit = bboptim_box(starts(:,k), fn, lower, upper, gr=gr)
        else
          fit = bboptim_box(starts(:,k), fn, lower, upper)
        end if
      end if
      ans%par(:,k) = fit%par
      ans%fvalue(k) = fit%value
      ans%converged(k) = fit%succeeded()
    end do
  end function multistart_optimize_box


  function multistart_optimize_projected(starts, fn, project, control, gr) result(ans)
    real(dp), intent(in) :: starts(:, :)
    procedure(bb_scalar_fn) :: fn
    procedure(bb_projection_fn) :: project
    type(bboptim_control), intent(in), optional :: control
    procedure(bb_gradient_fn), optional :: gr
    type(multistart_result) :: ans
    type(spg_result) :: fit
    integer :: k

    allocate(ans%par(size(starts,1), size(starts,2)))
    allocate(ans%fvalue(size(starts,2)), ans%converged(size(starts,2)))
    do k = 1, size(starts, 2)
      if (present(control)) then
        if (present(gr)) then
          fit = bboptim_projected(starts(:,k), fn, project, control, gr)
        else
          fit = bboptim_projected(starts(:,k), fn, project, control)
        end if
      else
        if (present(gr)) then
          fit = bboptim_projected(starts(:,k), fn, project, gr=gr)
        else
          fit = bboptim_projected(starts(:,k), fn, project)
        end if
      end if
      ans%par(:,k) = fit%par
      ans%fvalue(k) = fit%value
      ans%converged(k) = fit%succeeded()
    end do
  end function multistart_optimize_projected

  function multistart_optimize_linear(starts, fn, a, b, meq, control, gr) result(ans)
    real(dp), intent(in) :: starts(:, :), a(:, :), b(:)
    integer, intent(in) :: meq
    procedure(bb_scalar_fn) :: fn
    type(bboptim_control), intent(in), optional :: control
    procedure(bb_gradient_fn), optional :: gr
    type(multistart_result) :: ans
    type(spg_result) :: fit
    integer :: k

    allocate(ans%par(size(starts,1), size(starts,2)))
    allocate(ans%fvalue(size(starts,2)), ans%converged(size(starts,2)))
    do k = 1, size(starts, 2)
      if (present(control)) then
        if (present(gr)) then
          fit = bboptim_linear(starts(:,k), fn, a, b, meq, control, gr)
        else
          fit = bboptim_linear(starts(:,k), fn, a, b, meq, control)
        end if
      else
        if (present(gr)) then
          fit = bboptim_linear(starts(:,k), fn, a, b, meq, gr=gr)
        else
          fit = bboptim_linear(starts(:,k), fn, a, b, meq)
        end if
      end if
      ans%par(:,k) = fit%par
      ans%fvalue(k) = fit%value
      ans%converged(k) = fit%succeeded()
    end do
  end function multistart_optimize_linear

  subroutine fill_spg_control(sc, ctrl, method, m)
    type(spg_control), intent(out) :: sc
    type(bboptim_control), intent(in) :: ctrl
    integer, intent(in) :: method, m
    sc = spg_control(m=m, maxit=ctrl%maxit, maxfeval=ctrl%maxfeval, method=method, &
      ftol=ctrl%ftol, gtol=ctrl%gtol, eps=ctrl%eps, maximize=ctrl%maximize, &
      trace=ctrl%trace, triter=ctrl%triter)
  end subroutine fill_spg_control

end module bb_drivers
