! SPDX-License-Identifier: GPL-2.0-or-later
program test_sann
    use adequacy_model, only: dp, optimize_result, sann_optimize
    implicit none
    real(dp) :: start(2), data(1)
    integer :: nseed
    integer, allocatable :: seed(:)
    type(optimize_result) :: res

    call random_seed(size=nseed)
    allocate(seed(nseed))
    seed = 24680
    call random_seed(put=seed)
    start = [4.0_dp, 4.0_dp]
    data = 0.0_dp
    call sann_optimize(quadratic, data, start, res, temp0=2.0_dp, cooling=0.998_dp, max_iter=12000)
    if (res%value > 2.0e-2_dp) then
        print *, res%par, res%value
        error stop 1
    end if
    print '(a)', 'test_sann: PASS'
contains
    function quadratic(par, x) result(v)
        real(dp), intent(in) :: par(:), x(:)
        real(dp) :: v
        v = (par(1)-1.0_dp)**2 + (par(2)+2.0_dp)**2 + 0.0_dp*sum(x)
    end function quadratic
end program test_sann
