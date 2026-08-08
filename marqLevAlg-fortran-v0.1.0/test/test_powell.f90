program test_powell
    use marqlevalg, only : dp, mla_control, mla_result, marqlev_optimize
    implicit none
    type(mla_control) :: ctl
    type(mla_result) :: r
    real(dp) :: x0(4)

    ctl = mla_control()
    ctl%maxiter = 100
    ctl%epsa = 1.0e-3_dp
    ctl%epsb = 1.0e-3_dp
    ctl%epsd = 1.0e-3_dp
    x0 = [3.0_dp, -1.0_dp, 0.0_dp, 1.0_dp]
    call marqlev_optimize(x0, f, r, ctl)
    if (r%istop /= 1) error stop 1
    if (r%fn_value > 1.0e-3_dp) error stop 2
    if (maxval(abs(r%par)) > 0.1_dp) error stop 3
contains
    function f(x) result(v)
        real(dp), intent(in) :: x(:)
        real(dp) :: v
        v = (x(1) + 10.0_dp * x(2))**2 + 5.0_dp * (x(3) - x(4))**2 &
            + (x(2) - 2.0_dp * x(3))**4 + 10.0_dp * (x(1) - x(4))**4
    end function f
end program test_powell
