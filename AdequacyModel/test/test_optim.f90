! SPDX-License-Identifier: GPL-2.0-or-later
program test_optim
    use adequacy_model, only: dp, optimize_result, bfgs_optimize, nelder_mead_optimize, cg_optimize
    implicit none
    real(dp) :: start(2), data(1)
    type(optimize_result) :: res

    start = [5.0_dp, 5.0_dp]
    data = 0.0_dp

    call bfgs_optimize(quadratic, data, start, res)
    call check(res, 1.0e-5_dp, 'BFGS')

    call nelder_mead_optimize(quadratic, data, start, res)
    call check(res, 1.0e-4_dp, 'Nelder-Mead')

    call cg_optimize(quadratic, data, start, res)
    call check(res, 1.0e-5_dp, 'CG')

    print '(a)', 'test_optim: PASS'
contains
    function quadratic(par, x) result(v)
        real(dp), intent(in) :: par(:), x(:)
        real(dp) :: v
        v = (par(1)-1.0_dp)**2 + (par(2)+2.0_dp)**2 + 0.0_dp*sum(x)
    end function quadratic

    subroutine check(r, tol, label)
        type(optimize_result), intent(in) :: r
        real(dp), intent(in) :: tol
        character(len=*), intent(in) :: label
        if (maxval(abs(r%par - [1.0_dp, -2.0_dp])) > tol) then
            print *, trim(label), r%par, r%value
            error stop 1
        end if
    end subroutine check
end program test_optim
