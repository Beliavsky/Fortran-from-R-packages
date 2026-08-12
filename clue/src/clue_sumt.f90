! SPDX-License-Identifier: GPL-2.0-only
module clue_sumt
    use clue_kinds, only: dp
    implicit none
    private
    public :: sumt_result, sumt_optimize

    abstract interface
        function scalar_fun(x) result(v)
            import dp
            real(dp), intent(in) :: x(:)
            real(dp) :: v
        end function scalar_fun

        subroutine grad_fun(x, g)
            import dp
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: g(:)
        end subroutine grad_fun
    end interface

    type :: sumt_result
        real(dp), allocatable :: x(:)
        real(dp) :: objective = huge(1.0_dp)
        real(dp) :: penalty = huge(1.0_dp)
        real(dp) :: rho = 0.0_dp
        integer :: outer_iterations = 0
    end type sumt_result

contains

    function sumt_optimize(x0, L, P, grad_L, grad_P, eps, q, max_outer, max_inner) result(out)
        real(dp), intent(in) :: x0(:)
        procedure(scalar_fun) :: L, P
        procedure(grad_fun), optional :: grad_L
        procedure(grad_fun), optional :: grad_P
        real(dp), intent(in), optional :: eps, q
        integer, intent(in), optional :: max_outer, max_inner
        type(sumt_result) :: out
        real(dp), allocatable :: x(:), xold(:), g(:), d(:), gnew(:)
        real(dp) :: tol, qq, rho, alpha, f, fn, beta, gg
        integer :: mo, mi, it, j
        logical :: have_analytic_grad

        tol = sqrt(epsilon(1.0_dp))
        if (present(eps)) tol = eps
        qq = 10.0_dp
        if (present(q)) qq = q
        mo = 100
        if (present(max_outer)) mo = max_outer
        mi = 500
        if (present(max_inner)) mi = max_inner

        x = x0
        rho = max(L(x), 1.0e-5_dp) / max(P(x), 1.0e-5_dp)
        allocate(g(size(x)), d(size(x)), gnew(size(x)))
        have_analytic_grad = present(grad_L) .and. present(grad_P)

        do it = 1, mo
            xold = x
            if (have_analytic_grad) then
                call phi_grad_analytic(x, rho, grad_L, grad_P, g)
            else
                call phi_grad_fd(x, rho, L, P, g)
            end if
            d = -g

            do j = 1, mi
                gg = dot_product(g, g)
                if (sqrt(gg) < 1.0e-10_dp) exit
                f = phi_value(x, rho, L, P)
                alpha = 1.0_dp
                do while (alpha > 1.0e-12_dp)
                    fn = phi_value(x + alpha * d, rho, L, P)
                    if (fn <= f + 1.0e-4_dp * alpha * dot_product(g, d)) exit
                    alpha = 0.5_dp * alpha
                end do
                x = x + alpha * d
                if (have_analytic_grad) then
                    call phi_grad_analytic(x, rho, grad_L, grad_P, gnew)
                else
                    call phi_grad_fd(x, rho, L, P, gnew)
                end if
                beta = max(0.0_dp, dot_product(gnew, gnew - g) / max(gg, tiny(1.0_dp)))
                d = -gnew + beta * d
                g = gnew
            end do

            if (maxval(abs(x - xold)) < tol) exit
            rho = qq * rho
        end do

        allocate(out%x(size(x)))
        out%x = x
        out%objective = L(x)
        out%penalty = P(x)
        out%rho = rho
        out%outer_iterations = min(it, mo)
    end function sumt_optimize

    function phi_value(x, rho, L, P) result(v)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: rho
        procedure(scalar_fun) :: L, P
        real(dp) :: v

        v = L(x) + rho * P(x)
    end function phi_value

    subroutine phi_grad_analytic(x, rho, grad_L, grad_P, g)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: rho
        procedure(grad_fun) :: grad_L, grad_P
        real(dp), intent(out) :: g(:)
        real(dp), allocatable :: a(:), b(:)

        allocate(a(size(x)), b(size(x)))
        call grad_L(x, a)
        call grad_P(x, b)
        g = a + rho * b
    end subroutine phi_grad_analytic

    subroutine phi_grad_fd(x, rho, L, P, g)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: rho
        procedure(scalar_fun) :: L, P
        real(dp), intent(out) :: g(:)
        real(dp), allocatable :: xp(:), xm(:)
        real(dp) :: h
        integer :: i

        allocate(xp(size(x)), xm(size(x)))
        do i = 1, size(x)
            h = sqrt(epsilon(1.0_dp)) * (1.0_dp + abs(x(i)))
            xp = x
            xm = x
            xp(i) = xp(i) + h
            xm(i) = xm(i) - h
            g(i) = (phi_value(xp, rho, L, P) - phi_value(xm, rho, L, P)) / (2.0_dp * h)
        end do
    end subroutine phi_grad_fd

end module clue_sumt
