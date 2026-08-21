! SPDX-License-Identifier: GPL-3.0-or-later
module cubintegrate_mod
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use cubature_kinds, only : dp, i8
    use cubature_types, only : cubature_result, cubature_integrand, cuhre_options, divonne_options, &
        suave_options, vegas_options, CUBATURE_BADARG
    use hcubature_mod, only : hcubature
    use pcubature_mod, only : pcubature
    use cuba_mod, only : cuhre, divonne, suave, vegas
    implicit none
    private
    public :: cubintegrate

    procedure(cubature_integrand), pointer, save :: active_integrand => null()

contains

    subroutine cubintegrate(f, lower, upper, fdim, method, result, rel_tol, abs_tol, max_eval)
        procedure(cubature_integrand) :: f
        real(dp), intent(in) :: lower(:), upper(:)
        integer, intent(in) :: fdim
        character(len=*), intent(in) :: method
        type(cubature_result), intent(out) :: result
        real(dp), intent(in), optional :: rel_tol, abs_tol
        integer(i8), intent(in), optional :: max_eval
        real(dp) :: rt, at
        integer(i8) :: me
        real(dp), allocatable :: alower(:), aupper(:)
        logical :: transformed

        rt = 1.0e-5_dp
        at = 1.0e-12_dp
        me = 1000000_i8
        if (present(rel_tol)) rt = rel_tol
        if (present(abs_tol)) at = abs_tol
        if (present(max_eval)) me = max_eval
        if (size(lower) /= size(upper) .or. size(lower) == 0) then
            call bad_result(result, fdim)
            return
        end if

        transformed = any(.not. ieee_is_finite(lower)) .or. any(.not. ieee_is_finite(upper))
        if (transformed) then
            allocate(alower(size(lower)), aupper(size(upper)))
            alower = atan_bound(lower)
            aupper = atan_bound(upper)
            active_integrand => f
            call dispatch_method(transformed_callback, alower, aupper, fdim, method, result, rt, at, me)
            nullify(active_integrand)
        else
            call dispatch_method(f, lower, upper, fdim, method, result, rt, at, me)
        end if
    end subroutine cubintegrate

    subroutine dispatch_method(fun, lower, upper, fdim, method, result, rel_tol, abs_tol, max_eval)
        procedure(cubature_integrand) :: fun
        real(dp), intent(in) :: lower(:), upper(:), rel_tol, abs_tol
        integer, intent(in) :: fdim
        character(len=*), intent(in) :: method
        type(cubature_result), intent(out) :: result
        integer(i8), intent(in) :: max_eval
        type(cuhre_options) :: co
        type(divonne_options) :: dio
        type(suave_options) :: so
        type(vegas_options) :: vo
        character(len=:), allocatable :: name

        name = lowercase(trim(method))
        select case (name)
        case ('hcubature', 'adaptintegrate')
            call hcubature(fun, lower, upper, fdim, result, rel_tol, abs_tol, max_eval)
        case ('pcubature')
            call pcubature(fun, lower, upper, fdim, result, rel_tol, abs_tol, max_eval)
        case ('cuhre')
            co%max_eval = max_eval
            call cuhre(fun, lower, upper, fdim, result, rel_tol, abs_tol, co)
        case ('divonne')
            dio%max_eval = max_eval
            call divonne(fun, lower, upper, fdim, result, rel_tol, abs_tol, dio)
        case ('suave')
            so%max_eval = max_eval
            call suave(fun, lower, upper, fdim, result, rel_tol, abs_tol, so)
        case ('vegas')
            vo%max_eval = max_eval
            call vegas(fun, lower, upper, fdim, result, rel_tol, abs_tol, vo)
        case default
            call bad_result(result, fdim)
        end select
    end subroutine dispatch_method

    subroutine transformed_callback(t, value)
        real(dp), intent(in) :: t(:)
        real(dp), intent(out) :: value(:)
        real(dp) :: x(size(t)), jac
        integer :: d
        if (.not. associated(active_integrand)) then
            value = 0.0_dp
            return
        end if
        x = tan(t)
        jac = 1.0_dp
        do d = 1, size(t)
            jac = jac / (cos(t(d)) * cos(t(d)))
        end do
        call active_integrand(x, value)
        value = jac * value
    end subroutine transformed_callback

    pure function atan_bound(x) result(y)
        real(dp), intent(in) :: x(:)
        real(dp) :: y(size(x))
        real(dp), parameter :: halfpi = 0.5_dp * acos(-1.0_dp)
        integer :: i
        do i = 1, size(x)
            if (ieee_is_finite(x(i))) then
                y(i) = atan(x(i))
            else if (x(i) < 0.0_dp) then
                y(i) = -halfpi
            else
                y(i) = halfpi
            end if
        end do
    end function atan_bound

    pure function lowercase(s) result(out)
        character(len=*), intent(in) :: s
        character(len=len(s)) :: out
        integer :: i, c
        do i = 1, len(s)
            c = iachar(s(i:i))
            if (c >= iachar('A') .and. c <= iachar('Z')) then
                out(i:i) = achar(c + iachar('a') - iachar('A'))
            else
                out(i:i) = s(i:i)
            end if
        end do
    end function lowercase

    subroutine bad_result(result, fdim)
        type(cubature_result), intent(out) :: result
        integer, intent(in) :: fdim
        allocate(result%integral(max(0, fdim)), result%error(max(0, fdim)), result%prob(max(0, fdim)))
        result%integral = 0.0_dp
        result%error = huge(1.0_dp)
        result%prob = 1.0_dp
        result%return_code = CUBATURE_BADARG
        result%evaluations = 0_i8
        result%nregions = 0
    end subroutine bad_result

end module cubintegrate_mod
