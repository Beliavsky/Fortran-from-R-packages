! SPDX-License-Identifier: GPL-2.0-or-later
!
! Native Fortran replacements for the small Rmath density bridge used by KFAS.

subroutine dnormf(x, mu, sigma, res)
    use kfas_kinds, only: dp
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    real(real64), intent(in) :: x, mu, sigma
    real(real64), intent(inout) :: res
    real(real64), parameter :: pi = acos(-1.0_real64)
    real(real64) :: log_density

    if (sigma <= 0.0_real64) then
        res = huge(res)
        return
    end if
    log_density = -0.5_real64 * log(2.0_real64 * pi) - log(sigma) &
        - 0.5_real64 * ((x - mu) / sigma)**2
    res = res - log_density
end subroutine dnormf

subroutine dpoisf(x, lambda, res)
    use kfas_kinds, only: dp
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    real(real64), intent(in) :: x, lambda
    real(real64), intent(inout) :: res

    if (lambda < 0.0_real64 .or. x < 0.0_real64) then
        res = -huge(res)
    else if (lambda == 0.0_real64) then
        if (x /= 0.0_real64) res = -huge(res)
    else
        res = res + x * log(lambda) - lambda - log_gamma(x + 1.0_real64)
    end if
end subroutine dpoisf

subroutine dbinomf(x, n, prob, res)
    use kfas_kinds, only: dp
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    real(real64), intent(in) :: x, n, prob
    real(real64), intent(inout) :: res
    real(real64) :: log_density

    if (x < 0.0_real64 .or. x > n .or. prob < 0.0_real64 .or. prob > 1.0_real64) then
        res = -huge(res)
        return
    end if
    log_density = log_gamma(n + 1.0_real64) - log_gamma(x + 1.0_real64) &
        - log_gamma(n - x + 1.0_real64)
    if (prob == 0.0_real64) then
        if (x == 0.0_real64) then
            res = res + log_density
        else
            res = -huge(res)
        end if
    else if (prob == 1.0_real64) then
        if (x == n) then
            res = res + log_density
        else
            res = -huge(res)
        end if
    else
        res = res + log_density + x * log(prob) + (n - x) * log(1.0_real64 - prob)
    end if
end subroutine dbinomf

subroutine dgammaf(x, shape, scale, res)
    use kfas_kinds, only: dp
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    real(real64), intent(in) :: x, shape, scale
    real(real64), intent(inout) :: res
    real(real64) :: log_density

    if (x < 0.0_real64 .or. shape <= 0.0_real64 .or. scale <= 0.0_real64) then
        res = -huge(res)
        return
    end if
    if (x == 0.0_real64) then
        if (shape < 1.0_real64) then
            res = huge(res)
        else if (shape == 1.0_real64) then
            res = res - log(scale)
        else
            res = -huge(res)
        end if
        return
    end if
    log_density = (shape - 1.0_real64) * log(x) - x / scale &
        - log_gamma(shape) - shape * log(scale)
    res = res + log_density
end subroutine dgammaf

subroutine dnbinomf(x, size, mu, res)
    use kfas_kinds, only: dp
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    real(real64), intent(in) :: x, size, mu
    real(real64), intent(inout) :: res
    real(real64) :: log_density

    if (x < 0.0_real64 .or. size <= 0.0_real64 .or. mu < 0.0_real64) then
        res = -huge(res)
        return
    end if
    if (mu == 0.0_real64) then
        if (x /= 0.0_real64) res = -huge(res)
        return
    end if
    log_density = log_gamma(x + size) - log_gamma(size) - log_gamma(x + 1.0_real64) &
        + size * log(size / (size + mu)) + x * log(mu / (size + mu))
    res = res + log_density
end subroutine dnbinomf

integer function finitex(x) result(is_finite)
    use kfas_kinds, only: dp
    use, intrinsic :: iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(real64), intent(in) :: x

    is_finite = merge(1, 0, ieee_is_finite(x))
end function finitex
