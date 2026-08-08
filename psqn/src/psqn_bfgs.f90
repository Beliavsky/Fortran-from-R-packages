! SPDX-License-Identifier: Apache-2.0
module psqn_bfgs
  use psqn_types, only : dp, psqn_objective_eval, psqn_info, psqn_bfgs_options, &
                         psqn_converged, psqn_max_it_reached, psqn_line_search_failed
  use psqn_linalg, only : packed_size, sym_packed_matvec, bfgs_inverse_update_packed
  use psqn_linesearch, only : wolfe_line_search
  implicit none
  private
  public :: psqn_bfgs_optimize

contains

  subroutine psqn_bfgs_optimize(x, evaluator, result, options)
    real(dp), intent(inout) :: x(:)
    procedure(psqn_objective_eval) :: evaluator
    type(psqn_info), intent(out) :: result
    type(psqn_bfgs_options), intent(in), optional :: options

    type(psqn_bfgs_options) :: opt
    real(dp), allocatable :: x_old(:), gr(:), gr_old(:), s(:), y(:), wrk(:), dir(:), h(:)
    real(dp) :: fval, fval_old, sy, yhy, scale, err
    integer :: i, j, n, n_line_search_fail
    logical :: first_call, ls_ok, all_unchanged, has_converged

    opt = psqn_bfgs_options()
    if (present(options)) opt = options
    if (opt%c1 < 0.0_dp .or. opt%c1 >= opt%c2 .or. opt%c2 >= 1.0_dp) &
      error stop "psqn_bfgs: invalid Wolfe constants"

    n = size(x)
    allocate(x_old(n), gr(n), gr_old(n), s(n), y(n), wrk(n), dir(n), h(packed_size(n)))
    result = psqn_info()
    call reset_h()
    call evaluator(x, fval, gr, .true.)
    result%n_grad = 1
    call record_state()
    result%info = psqn_max_it_reached
    n_line_search_fail = 0

    do i = 1, opt%max_it
      fval_old = fval
      dir = 0.0_dp
      call sym_packed_matvec(h, gr, dir)
      dir = -dir

      ls_ok = wolfe_line_search(evaluator, fval_old, x, gr, dir, fval, opt%c1, opt%c2, &
                                .true., result%n_eval, result%n_grad, opt%trace)
      if (.not. ls_ok) then
        result%info = psqn_line_search_failed
        n_line_search_fail = n_line_search_fail + 1
        if (n_line_search_fail > 2) exit
      else
        n_line_search_fail = 0
      end if

      if (opt%trace > 0) then
        write(*,'(a,i0,2(a,es16.8),a,l1)') 'bfgs iter=', i, ' f=', fval, &
          ' |g|=', sqrt(dot_product(gr,gr)), ' line_search=', ls_ok
      end if

      err = abs(fval - fval_old)
      has_converged = n_line_search_fail < 1 .and. &
        err < opt%rel_eps * (abs(fval_old) + opt%rel_eps) .and. &
        (opt%abs_eps < 0.0_dp .or. err < opt%abs_eps) .and. &
        (opt%gr_tol < 0.0_dp .or. dot_product(gr,gr) < opt%gr_tol**2)
      if (has_converged) then
        result%info = psqn_converged
        exit
      end if

      if (n_line_search_fail < 2) then
        s = x - x_old
        all_unchanged = .true.
        do j = 1, n
          if (abs(s(j)) > abs(x(j)) * epsilon(1.0_dp) * 100.0_dp) then
            all_unchanged = .false.
            exit
          end if
        end do

        if (.not. all_unchanged) then
          y = gr - gr_old
          sy = dot_product(y, s)
          if (sy > 0.0_dp) then
            if (first_call) then
              first_call = .false.
              scale = sy / dot_product(y, y)
              h = 0.0_dp
              do j = 1, n
                h(j + j*(j-1)/2) = scale
              end do
            end if
            wrk = 0.0_dp
            call sym_packed_matvec(h, y, wrk)
            yhy = dot_product(y, wrk)
            call bfgs_inverse_update_packed(h, s, wrk, yhy, 1.0_dp / sy)
          else
            call reset_h()
          end if
        else
          call reset_h()
        end if
        call record_state()
      else
        call reset_h()
        call record_state()
      end if
    end do

    result%value = fval

  contains
    subroutine reset_h()
      integer :: ii
      h = 0.0_dp
      do ii = 1, n
        h(ii + ii*(ii-1)/2) = 1.0_dp
      end do
      first_call = .true.
    end subroutine reset_h

    subroutine record_state()
      x_old = x
      gr_old = gr
    end subroutine record_state
  end subroutine psqn_bfgs_optimize

end module psqn_bfgs
