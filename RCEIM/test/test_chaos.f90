module test_chaos_objective
    use rceim, only : dp
    implicit none
contains
    function flat(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f
        f = 1.0_dp + 0.0_dp*x(1)
    end function flat
end module test_chaos_objective

program test_chaos
    use rceim, only : dp, rceim_options, rceim_result, ceim_optimize
    use test_chaos_objective, only : flat
    implicit none
    type(rceim_options) :: opt
    type(rceim_result) :: res
    real(dp) :: lo(1), hi(1)
    lo = [-2.0_dp]
    hi = [ 2.0_dp]
    opt%n_total = 80
    opt%n_elite = 20
    opt%n_super = 3
    opt%max_iter = 8
    opt%wait_gen = 20
    opt%chaos_gen = 1
    opt%epsilon = 0.01_dp
    opt%alpha = 0.0_dp
    opt%seed = 77
    call ceim_optimize(flat, lo, hi, res, opt)
    if (any(res%elite_x < -2.0_dp) .or. any(res%elite_x > 2.0_dp)) error stop 1
    if (res%history_n < 1) error stop 2
end program test_chaos
