! SPDX-License-Identifier: GPL-2.0-or-later
module hyper2_integration
    use hyper2_kinds, only : dp, i8
    use hyper2_types, only : hyper2_model, operator(+)
    use hyper2_likelihood, only : loglik_h2
    use hyper2_models, only : dirichlet
    use cubature, only : cubature_result, adapt_integrate, ERROR_INDIVIDUAL
    implicit none
    private

    type(hyper2_model), allocatable, save :: active_h

    public :: p_to_e, e_to_p, simplex_jacobian, hyper2_B, dhyper2, mgf, mean_hyper2

contains

    function p_to_e(p) result(e)
        real(dp), intent(in) :: p(:)
        real(dp), allocatable :: e(:)
        real(dp) :: tail
        integer :: n, i
        n = size(p)
        allocate(e(n))
        if (n == 0) return
        e(1) = sum(p)
        do i = 2, n
            tail = sum(p(i-1:n))
            if (tail > 0.0_dp) then
                e(i) = p(i-1) / tail
            else
                e(i) = 0.0_dp
            end if
        end do
    end function p_to_e

    function e_to_p(e) result(p)
        real(dp), intent(in) :: e(:)
        real(dp), allocatable :: p(:)
        real(dp) :: prod
        integer :: n, i
        n = size(e)
        allocate(p(n))
        if (n == 0) return
        prod = e(1)
        do i = 1, n - 1
            p(i) = prod * e(i + 1)
            prod = prod * (1.0_dp - e(i + 1))
        end do
        p(n) = prod
    end function e_to_p

    real(dp) function simplex_jacobian(e_full) result(jac)
        real(dp), intent(in) :: e_full(:)
        real(dp) :: prod
        integer :: n, i
        n = size(e_full)
        if (n <= 1) then
            jac = 1.0_dp
            return
        end if
        jac = e_full(1)**n
        prod = 1.0_dp
        do i = 2, n - 1
            prod = prod * (1.0_dp - e_full(i))
            jac = jac * prod
        end do
    end function simplex_jacobian

    subroutine b_integrand(x, value)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: value(:)
        real(dp), allocatable :: e(:), p(:)
        integer :: n
        n = size(x) + 1
        allocate(e(n))
        e(1) = 1.0_dp
        if (size(x) > 0) e(2:n) = x
        p = e_to_p(e)
        value(1) = exp(loglik_h2(p, active_h)) * simplex_jacobian(e)
    end subroutine b_integrand

    real(dp) function hyper2_B(h, rel_tol, abs_tol, max_eval, error_estimate, evaluations) result(val)
        type(hyper2_model), intent(in) :: h
        real(dp), intent(in), optional :: rel_tol, abs_tol
        integer(i8), intent(in), optional :: max_eval
        real(dp), intent(out), optional :: error_estimate
        integer(i8), intent(out), optional :: evaluations
        type(cubature_result) :: res
        real(dp), allocatable :: lo(:), hi(:)
        real(dp) :: rt, at
        integer(i8) :: me
        integer :: dim
        dim = h%size() - 1
        if (dim < 0) then
            val = 0.0_dp
            return
        else if (dim == 0) then
            val = exp(loglik_h2([1.0_dp], h))
            if (present(error_estimate)) error_estimate = 0.0_dp
            if (present(evaluations)) evaluations = 1_i8
            return
        end if
        rt = 1.0e-8_dp
        if (present(rel_tol)) rt = rel_tol
        at = 1.0e-12_dp
        if (present(abs_tol)) at = abs_tol
        me = 1000000_i8
        if (present(max_eval)) me = max_eval
        allocate(lo(dim), hi(dim))
        lo = 0.0_dp
        hi = 1.0_dp
        active_h = h
        call adapt_integrate(b_integrand, lo, hi, 1, res, tol=rt, abs_error=at, max_eval=me, norm=ERROR_INDIVIDUAL)
        val = res%integral(1)
        if (present(error_estimate)) error_estimate = res%error(1)
        if (present(evaluations)) evaluations = res%evaluations
        if (allocated(active_h)) deallocate(active_h)
    end function hyper2_B

    real(dp) function dhyper2(p, h, normalizer) result(d)
        real(dp), intent(in) :: p(:)
        type(hyper2_model), intent(in) :: h
        real(dp), intent(in), optional :: normalizer
        real(dp) :: z
        if (present(normalizer)) then
            z = normalizer
        else
            z = hyper2_B(h)
        end if
        if (z <= 0.0_dp) then
            d = 0.0_dp
        else
            d = exp(loglik_h2(p,h)) / z
        end if
    end function dhyper2

    real(dp) function mgf(h, powers, rel_tol, abs_tol) result(v)
        type(hyper2_model), intent(in) :: h
        real(dp), intent(in) :: powers(:)
        real(dp), intent(in), optional :: rel_tol, abs_tol
        type(hyper2_model) :: shifted, d
        real(dp) :: z0, z1
        if (size(powers) /= h%size()) error stop "mgf powers size mismatch"
        d = dirichlet(h%pnames, powers=powers)
        shifted = h + d
        z0 = hyper2_B(h, rel_tol=rel_tol, abs_tol=abs_tol)
        z1 = hyper2_B(shifted, rel_tol=rel_tol, abs_tol=abs_tol)
        v = z1 / z0
    end function mgf

    function mean_hyper2(h, normalize, rel_tol, abs_tol) result(mu)
        type(hyper2_model), intent(in) :: h
        logical, intent(in), optional :: normalize
        real(dp), intent(in), optional :: rel_tol, abs_tol
        real(dp), allocatable :: mu(:), pw(:)
        logical :: norm
        integer :: i, n
        n = h%size()
        allocate(mu(n), pw(n))
        mu = 0.0_dp
        norm = .true.
        if (present(normalize)) norm = normalize
        do i = 1, n
            pw = 0.0_dp
            pw(i) = 1.0_dp
            mu(i) = mgf(h, pw, rel_tol=rel_tol, abs_tol=abs_tol)
        end do
        if (norm .and. sum(mu) > 0.0_dp) mu = mu / sum(mu)
    end function mean_hyper2

end module hyper2_integration
