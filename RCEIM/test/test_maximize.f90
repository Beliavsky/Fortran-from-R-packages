module test_maximize_objective
    use rceim, only : dp
    implicit none
contains
    function concave(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f
        f = 3.0_dp - (x(1)-1.25_dp)**2
    end function concave
end module test_maximize_objective

program test_maximize
    use rceim, only : dp, rceim_options, rceim_result, ceim_optimize
    use test_maximize_objective, only : concave
    implicit none
    type(rceim_options) :: opt
    type(rceim_result) :: res
    real(dp) :: lo(1), hi(1)

    lo = [-5.0_dp]
    hi = [ 5.0_dp]
    opt%n_total = 300
    opt%n_elite = 60
    opt%minimize = .false.
    opt%epsilon = 0.002_dp
    opt%max_iter = 70
    opt%seed = 991
    call ceim_optimize(concave, lo, hi, res, opt)
    if (abs(res%x(1)-1.25_dp) > 0.05_dp) error stop 1
    if (res%value < 2.99_dp) error stop 2
    if (res%score > -2.99_dp) error stop 3
end program test_maximize
