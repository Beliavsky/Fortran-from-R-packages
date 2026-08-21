! SPDX-License-Identifier: GPL-2.0-or-later
program test_pso
    use adequacy_model, only: dp, optimize_result, pso_optimize
    implicit none
    real(dp) :: lower(2), upper(2), data(1)
    integer, allocatable :: seed(:)
    type(optimize_result) :: res
    integer :: i

    call random_seed(size=i)
    allocate(seed(i))
    seed = 12345
    call random_seed(put=seed)

    lower = [-5.0_dp, -5.0_dp]
    upper = [5.0_dp, 5.0_dp]
    data = 0.0_dp
    call pso_optimize(quadratic, data, lower, upper, res, swarm_size=80, &
                      min_history=80, tol_var=1.0e-12_dp, max_iter=1000)
    if (maxval(abs(res%par - [1.0_dp, -2.0_dp])) > 2.0e-2_dp) then
        print *, res%par, res%value
        error stop 1
    end if
    print '(a)', 'test_pso: PASS'
contains
    function quadratic(par, x) result(v)
        real(dp), intent(in) :: par(:), x(:)
        real(dp) :: v
        v = (par(1)-1.0_dp)**2 + (par(2)+2.0_dp)**2 + 0.0_dp*sum(x)
    end function quadratic
end program test_pso
