module test_minimize_objective
    use rceim, only : dp
    implicit none
contains
    function quadratic(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f
        f = (x(1)+1.0_dp)**2 + (x(2)+2.0_dp)**2
    end function quadratic
end module test_minimize_objective

program test_minimize
    use rceim, only : dp, rceim_options, rceim_result, ceim_optimize
    use test_minimize_objective, only : quadratic
    implicit none
    type(rceim_options) :: opt
    type(rceim_result) :: res
    real(dp) :: lo(2), hi(2)

    lo = [-10.0_dp,-10.0_dp]
    hi = [ 10.0_dp, 10.0_dp]
    opt%n_total = 500
    opt%n_elite = 100
    opt%n_super = 2
    opt%epsilon = 0.005_dp
    opt%max_iter = 80
    opt%seed = 1234
    call ceim_optimize(quadratic, lo, hi, res, opt)
    if (res%value > 2.0e-3_dp) error stop 1
    if (maxval(abs(res%x-[-1.0_dp,-2.0_dp])) > 0.08_dp) error stop 2
end program test_minimize
