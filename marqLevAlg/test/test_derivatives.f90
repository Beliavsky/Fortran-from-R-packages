program test_derivatives
    use marqlevalg, only : dp, deriva, deriva_grad
    implicit none
    real(dp) :: x(2), f0, inf(2, 2), ing(2, 2), g(2)

    x = [1.2_dp, -0.7_dp]
    call deriva(x, f, f0, inf, g)
    if (abs(f0 - f(x)) > 1.0e-12_dp) error stop 1
    if (maxval(abs(g - [2.0_dp * x(1) + 3.0_dp, 6.0_dp * x(2)])) > 1.0e-7_dp) error stop 2
    call deriva_grad(x, gradf, ing)
    if (maxval(abs(ing + reshape([2.0_dp, 0.0_dp, 0.0_dp, 6.0_dp], [2, 2]))) > 1.0e-7_dp) &
        error stop 3
contains
    function f(z) result(v)
        real(dp), intent(in) :: z(:)
        real(dp) :: v
        v = z(1)**2 + 3.0_dp * z(1) + 3.0_dp * z(2)**2
    end function f

    subroutine gradf(z, v)
        real(dp), intent(in) :: z(:)
        real(dp), intent(out) :: v(:)
        v = [2.0_dp * z(1) + 3.0_dp, 6.0_dp * z(2)]
    end subroutine gradf
end program test_derivatives
