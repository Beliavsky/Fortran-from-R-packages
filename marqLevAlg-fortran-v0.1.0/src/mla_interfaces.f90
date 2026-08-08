! SPDX-License-Identifier: GPL-2.0-or-later
module mla_interfaces
    use mla_kinds, only : dp
    implicit none
    private

    public :: objective_fn, gradient_fn, hessian_fn

    abstract interface
        function objective_fn(x) result(f)
            import :: dp
            real(dp), intent(in) :: x(:)
            real(dp) :: f
        end function objective_fn

        subroutine gradient_fn(x, g)
            import :: dp
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: g(:)
        end subroutine gradient_fn

        subroutine hessian_fn(x, h)
            import :: dp
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: h(:, :)
        end subroutine hessian_fn
    end interface
end module mla_interfaces
