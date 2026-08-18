! SPDX-License-Identifier: GPL-3.0-or-later
!
! P-adaptive tensor-product Clenshaw-Curtis cubature.  This preserves the
! mathematical method of Steven G. Johnson's pcubature implementation while
! using a simpler recomputation strategy instead of the C value cache.
module pcubature_mod
    use cubature_kinds, only : dp, i8
    use cubature_types, only : cubature_result, cubature_integrand, cubature_integrand_v, &
        CUBATURE_SUCCESS, CUBATURE_MAXEVAL, CUBATURE_BADARG, ERROR_INDIVIDUAL, errors_converged
    use cubature_utils, only : clenshaw_curtis_rule, safe_pow_int
    implicit none
    private
    public :: pcubature, pcubature_v

contains

    subroutine pcubature_core(f, fv, lower, upper, fdim, result, rel_tol, abs_tol, max_eval, norm)
        procedure(cubature_integrand), optional :: f
        procedure(cubature_integrand_v), optional :: fv
        real(dp), intent(in) :: lower(:), upper(:)
        integer, intent(in) :: fdim
        type(cubature_result), intent(out) :: result
        real(dp), intent(in), optional :: rel_tol, abs_tol
        integer(i8), intent(in), optional :: max_eval
        integer, intent(in), optional :: norm
        real(dp) :: rtol, atol
        integer(i8) :: maxev, need, neval, used
        integer :: nrm, dim, level, n, status
        real(dp), allocatable :: current(:), previous(:), err(:)
        logical :: overflow, have_previous

        rtol = 1.0e-5_dp
        atol = 100.0_dp * epsilon(1.0_dp)
        maxev = 0_i8
        nrm = ERROR_INDIVIDUAL
        if (present(rel_tol)) rtol = rel_tol
        if (present(abs_tol)) atol = abs_tol
        if (present(max_eval)) maxev = max_eval
        if (present(norm)) nrm = norm
        dim = size(lower)
        call init_result(result, max(fdim, 0))
        if (fdim <= 0 .or. dim <= 0 .or. size(upper) /= dim .or. any(upper < lower)) then
            result%return_code = CUBATURE_BADARG
            return
        end if
        if (dim > 20) then
            result%return_code = CUBATURE_BADARG
            return
        end if
        allocate(current(fdim), previous(fdim), err(fdim))
        previous = 0.0_dp
        err = huge(1.0_dp)
        neval = 0_i8
        have_previous = .false.

        do level = 0, 12
            n = 2 ** (level + 1) + 1
            need = safe_pow_int(n, dim, overflow)
            if (overflow) exit
            if (maxev > 0_i8 .and. neval + need > maxev .and. have_previous) exit
            call tensor_cc(f, fv, lower, upper, n, fdim, current, used, status)
            if (status /= 0) then
                result%return_code = CUBATURE_BADARG
                return
            end if
            neval = neval + used
            if (have_previous) then
                err = abs(current - previous)
                if (errors_converged(current, err, atol, rtol, nrm)) exit
            end if
            previous = current
            have_previous = .true.
            if (maxev > 0_i8 .and. neval >= maxev) exit
        end do

        result%integral = current
        result%error = err
        result%evaluations = neval
        result%nregions = 1
        if (have_previous .and. errors_converged(current, err, atol, rtol, nrm)) then
            result%return_code = CUBATURE_SUCCESS
        else
            result%return_code = CUBATURE_MAXEVAL
        end if
        result%prob = 1.0_dp
    end subroutine pcubature_core

    subroutine pcubature(f, lower, upper, fdim, result, rel_tol, abs_tol, max_eval, norm)
        procedure(cubature_integrand) :: f
        real(dp), intent(in) :: lower(:), upper(:)
        integer, intent(in) :: fdim
        type(cubature_result), intent(out) :: result
        real(dp), intent(in), optional :: rel_tol, abs_tol
        integer(i8), intent(in), optional :: max_eval
        integer, intent(in), optional :: norm
        call pcubature_core(f=f, lower=lower, upper=upper, fdim=fdim, result=result, &
            rel_tol=rel_tol, abs_tol=abs_tol, max_eval=max_eval, norm=norm)
    end subroutine pcubature

    subroutine pcubature_v(f, lower, upper, fdim, result, rel_tol, abs_tol, max_eval, norm)
        procedure(cubature_integrand_v) :: f
        real(dp), intent(in) :: lower(:), upper(:)
        integer, intent(in) :: fdim
        type(cubature_result), intent(out) :: result
        real(dp), intent(in), optional :: rel_tol, abs_tol
        integer(i8), intent(in), optional :: max_eval
        integer, intent(in), optional :: norm
        call pcubature_core(fv=f, lower=lower, upper=upper, fdim=fdim, result=result, &
            rel_tol=rel_tol, abs_tol=abs_tol, max_eval=max_eval, norm=norm)
    end subroutine pcubature_v

    subroutine tensor_cc(f, fv, lower, upper, n, fdim, val, neval, status)
        procedure(cubature_integrand), optional :: f
        procedure(cubature_integrand_v), optional :: fv
        real(dp), intent(in) :: lower(:), upper(:)
        integer, intent(in) :: n, fdim
        real(dp), intent(out) :: val(fdim)
        integer(i8), intent(out) :: neval
        integer, intent(out) :: status
        real(dp), allocatable :: nodes(:), weights(:), x(:), fx(:)
        integer, allocatable :: idx(:)
        integer :: dim, d
        real(dp) :: wt, scale
        logical :: done

        dim = size(lower)
        call clenshaw_curtis_rule(n, nodes, weights)
        allocate(x(dim), fx(fdim), idx(dim))
        idx = 1
        val = 0.0_dp
        neval = 0_i8
        scale = product(0.5_dp * (upper - lower))
        done = .false.
        do while (.not. done)
            wt = scale
            do d = 1, dim
                x(d) = 0.5_dp * (lower(d) + upper(d)) + 0.5_dp * (upper(d) - lower(d)) * nodes(idx(d))
                wt = wt * weights(idx(d))
            end do
            call eval_point(f, fv, x, fx)
            val = val + wt * fx
            neval = neval + 1_i8
            call increment_index(idx, n, done)
        end do
        status = 0
    end subroutine tensor_cc

    subroutine increment_index(idx, n, done)
        integer, intent(inout) :: idx(:)
        integer, intent(in) :: n
        logical, intent(out) :: done
        integer :: d
        done = .false.
        do d = 1, size(idx)
            idx(d) = idx(d) + 1
            if (idx(d) <= n) return
            idx(d) = 1
        end do
        done = .true.
    end subroutine increment_index

    subroutine init_result(result, fdim)
        type(cubature_result), intent(out) :: result
        integer, intent(in) :: fdim
        allocate(result%integral(fdim), result%error(fdim), result%prob(fdim))
        result%integral = 0.0_dp
        result%error = huge(1.0_dp)
        result%prob = 1.0_dp
        result%evaluations = 0_i8
        result%return_code = CUBATURE_SUCCESS
        result%nregions = 0
    end subroutine init_result

    subroutine eval_point(f, fv, x, value)
        procedure(cubature_integrand), optional :: f
        procedure(cubature_integrand_v), optional :: fv
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: value(:)
        real(dp) :: xx(size(x), 1), vv(size(value), 1)
        if (present(f)) then
            call f(x, value)
        else if (present(fv)) then
            xx(:, 1) = x
            call fv(xx, vv)
            value = vv(:, 1)
        else
            value = 0.0_dp
        end if
    end subroutine eval_point

end module pcubature_mod
