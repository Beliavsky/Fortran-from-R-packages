program test_maximize
    use marqlevalg, only : dp, mla_control, mla_result, marqlev_optimize
    implicit none
    type(mla_control) :: ctl
    type(mla_result) :: r
    real(dp) :: x0(2)

    ctl = mla_control()
    ctl%minimize = .false.
    ctl%maxiter = 100
    ctl%epsa = 1.0e-3_dp
    ctl%epsb = 1.0e-3_dp
    ctl%epsd = 1.0e-3_dp
    x0 = [8.0_dp, 9.0_dp]
    call marqlev_optimize(x0, f, g, r, ctl)
    if (r%istop /= 1) error stop 1
    if (maxval(abs(r%par - [5.0_dp, 6.0_dp])) > 2.0e-3_dp) error stop 2
    if (abs(r%fn_value) > 1.0e-6_dp) error stop 3
contains
    function f(x) result(v)
        real(dp), intent(in) :: x(:)
        real(dp) :: v
        v = -4.0_dp * (x(1) - 5.0_dp)**2 - (x(2) - 6.0_dp)**2
    end function f

    subroutine g(x, v)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: v(:)
        v = [-8.0_dp * (x(1) - 5.0_dp), -2.0_dp * (x(2) - 6.0_dp)]
    end subroutine g
end program test_maximize
