module banana_functions
    use n1qn1_module, only : dp
    implicit none
contains
    function banana_value(x, user_data) result(f)
        real(dp), intent(in) :: x(:)
        class(*), intent(inout), optional :: user_data
        real(dp) :: f
        integer :: i
        if (present(user_data)) continue
        f = 1.0_dp
        do i = 2, size(x)
            f = f + 100.0_dp * (x(i) - x(i - 1) ** 2) ** 2 + (1.0_dp - x(i)) ** 2
        end do
    end function banana_value

    subroutine banana_gradient(x, g, user_data)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: g(:)
        class(*), intent(inout), optional :: user_data
        integer :: i, n
        if (present(user_data)) continue
        n = size(x)
        g = 0.0_dp
        g(1) = -400.0_dp * (x(2) - x(1) ** 2) * x(1)
        do i = 2, n - 1
            g(i) = 200.0_dp * (x(i) - x(i - 1) ** 2) &
                 - 400.0_dp * (x(i + 1) - x(i) ** 2) * x(i) &
                 - 2.0_dp * (1.0_dp - x(i))
        end do
        g(n) = 200.0_dp * (x(n) - x(n - 1) ** 2) - 2.0_dp * (1.0_dp - x(n))
    end subroutine banana_gradient
end module banana_functions

program banana_example
    use n1qn1_module, only : dp, n1qn1_control_t, n1qn1_result_t, n1qn1_minimize
    use banana_functions
    implicit none

    type(n1qn1_control_t) :: control
    type(n1qn1_result_t) :: result
    real(dp) :: x0(3)

    x0 = [1.02_dp, 1.02_dp, 1.02_dp]
    control%max_iterations = 100
    control%max_evaluations = 100

    call n1qn1_minimize(banana_value, banana_gradient, x0, result, control)

    write(*, '(a,*(f18.12,1x))') 'x: ', result%x
    write(*, '(a,es20.12)') 'f: ', result%value
    write(*, '(a,i0)') 'evaluations: ', result%function_evaluations
    write(*, '(a,a)') 'status: ', result%message
    write(*, '(a)') 'Hessian:'
    write(*, '(*(f16.8,1x))') result%hessian
end program banana_example
