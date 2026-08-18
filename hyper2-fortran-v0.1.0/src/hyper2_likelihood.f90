! SPDX-License-Identifier: GPL-2.0-or-later
module hyper2_likelihood
    use hyper2_kinds, only : dp
    use hyper2_types, only : hyper2_model, hyper3_model
    implicit none
    private

    interface loglik
        module procedure loglik_h2
        module procedure loglik_h3
    end interface
    interface gradient
        module procedure gradient_h2
        module procedure gradient_h3
    end interface
    interface gradient_full
        module procedure gradient_full_h2
        module procedure gradient_full_h3
    end interface
    interface hessian_independent
        module procedure hessian_independent_h2
        module procedure hessian_independent_h3
    end interface

    public :: loglik, loglik_h2, loglik_h3, gradient, gradient_full
    public :: gradient_h2, gradient_h3, gradient_full_h2, gradient_full_h3
    public :: hessian_independent, hessian_independent_h2, hessian_independent_h3
    public :: fillup, indep, equalp, power_sum_h2, power_sum_h3

contains

    function fillup(x, total) result(p)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in), optional :: total
        real(dp), allocatable :: p(:)
        real(dp) :: t
        t = 1.0_dp
        if (present(total)) t = total
        allocate(p(size(x) + 1))
        if (size(x) > 0) p(1:size(x)) = x
        p(size(p)) = t - sum(x)
    end function fillup

    function indep(p) result(x)
        real(dp), intent(in) :: p(:)
        real(dp), allocatable :: x(:)
        if (size(p) < 1) then
            allocate(x(0))
        else
            allocate(x(size(p) - 1))
            if (size(x) > 0) x = p(1:size(x))
        end if
    end function indep

    function equalp(n) result(p)
        integer, intent(in) :: n
        real(dp), allocatable :: p(:)
        if (n <= 0) then
            allocate(p(0))
            return
        end if
        allocate(p(n))
        p = 1.0_dp / real(n, dp)
    end function equalp

    real(dp) function power_sum_h2(h)
        type(hyper2_model), intent(in) :: h
        integer :: k
        power_sum_h2 = 0.0_dp
        do k = 1, size(h%terms)
            power_sum_h2 = power_sum_h2 + h%terms(k)%power
        end do
    end function power_sum_h2

    real(dp) function power_sum_h3(h)
        type(hyper3_model), intent(in) :: h
        integer :: k
        power_sum_h3 = 0.0_dp
        do k = 1, size(h%terms)
            power_sum_h3 = power_sum_h3 + h%terms(k)%power
        end do
    end function power_sum_h3

    subroutine normalize_input(p_in, n, p, ok)
        real(dp), intent(in) :: p_in(:)
        integer, intent(in) :: n
        real(dp), allocatable, intent(out) :: p(:)
        logical, intent(out) :: ok
        real(dp), parameter :: tol = 1.0e-6_dp
        ok = .false.
        if (size(p_in) == n - 1) then
            p = fillup(p_in)
        else if (size(p_in) == n) then
            p = p_in
        else
            allocate(p(0))
            return
        end if
        if (any(p < 0.0_dp)) return
        if (abs(sum(p) - 1.0_dp) > tol) return
        ok = .true.
    end subroutine normalize_input

    real(dp) function term_sum_h2(h, k, p) result(s)
        type(hyper2_model), intent(in) :: h
        integer, intent(in) :: k
        real(dp), intent(in) :: p(:)
        integer :: j
        s = 0.0_dp
        do j = 1, size(h%terms(k)%ids)
            s = s + p(h%terms(k)%ids(j))
        end do
    end function term_sum_h2

    real(dp) function term_sum_h3(h, k, p) result(s)
        type(hyper3_model), intent(in) :: h
        integer, intent(in) :: k
        real(dp), intent(in) :: p(:)
        integer :: j
        s = 0.0_dp
        do j = 1, size(h%terms(k)%ids)
            s = s + h%terms(k)%weights(j) * p(h%terms(k)%ids(j))
        end do
    end function term_sum_h3

    real(dp) function loglik_h2(p_in, h, log_scale) result(out)
        real(dp), intent(in) :: p_in(:)
        type(hyper2_model), intent(in) :: h
        logical, intent(in), optional :: log_scale
        real(dp), allocatable :: p(:)
        real(dp) :: s
        logical :: ok, lg
        integer :: k
        lg = .true.
        if (present(log_scale)) lg = log_scale
        call normalize_input(p_in, h%size(), p, ok)
        if (.not. ok) then
            out = -huge(1.0_dp)
            if (.not. lg) out = 0.0_dp
            return
        end if
        out = 0.0_dp
        do k = 1, size(h%terms)
            s = term_sum_h2(h, k, p)
            if (s <= 0.0_dp) then
                out = -huge(1.0_dp)
                if (.not. lg) out = 0.0_dp
                return
            end if
            out = out + h%terms(k)%power * log(s)
        end do
        if (.not. lg) out = exp(out)
    end function loglik_h2

    real(dp) function loglik_h3(p_in, h, log_scale) result(out)
        real(dp), intent(in) :: p_in(:)
        type(hyper3_model), intent(in) :: h
        logical, intent(in), optional :: log_scale
        real(dp), allocatable :: p(:)
        real(dp) :: s
        logical :: ok, lg
        integer :: k
        lg = .true.
        if (present(log_scale)) lg = log_scale
        call normalize_input(p_in, h%size(), p, ok)
        if (.not. ok) then
            out = -huge(1.0_dp)
            if (.not. lg) out = 0.0_dp
            return
        end if
        out = 0.0_dp
        do k = 1, size(h%terms)
            s = term_sum_h3(h, k, p)
            if (s <= 0.0_dp) then
                out = -huge(1.0_dp)
                if (.not. lg) out = 0.0_dp
                return
            end if
            out = out + h%terms(k)%power * log(s)
        end do
        if (.not. lg) out = exp(out)
    end function loglik_h3

    function gradient_full_h2(h, p_in) result(g)
        type(hyper2_model), intent(in) :: h
        real(dp), intent(in) :: p_in(:)
        real(dp), allocatable :: g(:), p(:)
        real(dp) :: s
        integer :: i, j, k
        logical :: ok
        call normalize_input(p_in, h%size(), p, ok)
        allocate(g(h%size()))
        g = 0.0_dp
        if (.not. ok) then
            g = huge(1.0_dp)
            return
        end if
        do k = 1, size(h%terms)
            s = term_sum_h2(h, k, p)
            if (s <= 0.0_dp) then
                g = huge(1.0_dp)
                return
            end if
            do j = 1, size(h%terms(k)%ids)
                i = h%terms(k)%ids(j)
                g(i) = g(i) + h%terms(k)%power / s
            end do
        end do
    end function gradient_full_h2

    function gradient_full_h3(h, p_in) result(g)
        type(hyper3_model), intent(in) :: h
        real(dp), intent(in) :: p_in(:)
        real(dp), allocatable :: g(:), p(:)
        real(dp) :: s
        integer :: i, j, k
        logical :: ok
        call normalize_input(p_in, h%size(), p, ok)
        allocate(g(h%size()))
        g = 0.0_dp
        if (.not. ok) then
            g = huge(1.0_dp)
            return
        end if
        do k = 1, size(h%terms)
            s = term_sum_h3(h, k, p)
            if (s <= 0.0_dp) then
                g = huge(1.0_dp)
                return
            end if
            do j = 1, size(h%terms(k)%ids)
                i = h%terms(k)%ids(j)
                g(i) = g(i) + h%terms(k)%power * h%terms(k)%weights(j) / s
            end do
        end do
    end function gradient_full_h3

    function gradient_h2(h, p_in) result(g)
        type(hyper2_model), intent(in) :: h
        real(dp), intent(in) :: p_in(:)
        real(dp), allocatable :: g(:), gf(:), p(:)
        logical :: ok
        integer :: n
        n = h%size()
        call normalize_input(p_in, n, p, ok)
        allocate(g(max(0, n - 1)))
        g = 0.0_dp
        if (.not. ok .or. n < 2) return
        gf = gradient_full_h2(h, p)
        g = gf(1:n - 1) - gf(n)
    end function gradient_h2

    function gradient_h3(h, p_in) result(g)
        type(hyper3_model), intent(in) :: h
        real(dp), intent(in) :: p_in(:)
        real(dp), allocatable :: g(:), gf(:), p(:)
        logical :: ok
        integer :: n
        n = h%size()
        call normalize_input(p_in, n, p, ok)
        allocate(g(max(0, n - 1)))
        g = 0.0_dp
        if (.not. ok .or. n < 2) return
        gf = gradient_full_h3(h, p)
        g = gf(1:n - 1) - gf(n)
    end function gradient_h3

    function weight_vector_h2(h, k) result(w)
        type(hyper2_model), intent(in) :: h
        integer, intent(in) :: k
        real(dp), allocatable :: w(:)
        integer :: j
        allocate(w(h%size()))
        w = 0.0_dp
        do j = 1, size(h%terms(k)%ids)
            w(h%terms(k)%ids(j)) = 1.0_dp
        end do
    end function weight_vector_h2

    function weight_vector_h3(h, k) result(w)
        type(hyper3_model), intent(in) :: h
        integer, intent(in) :: k
        real(dp), allocatable :: w(:)
        integer :: j
        allocate(w(h%size()))
        w = 0.0_dp
        do j = 1, size(h%terms(k)%ids)
            w(h%terms(k)%ids(j)) = h%terms(k)%weights(j)
        end do
    end function weight_vector_h3

    function hessian_independent_h2(h, p_in) result(hess)
        type(hyper2_model), intent(in) :: h
        real(dp), intent(in) :: p_in(:)
        real(dp), allocatable :: hess(:,:), p(:), w(:), d(:)
        real(dp) :: s
        integer :: n, k, i, j
        logical :: ok
        n = h%size()
        allocate(hess(max(0,n-1), max(0,n-1)))
        hess = 0.0_dp
        call normalize_input(p_in, n, p, ok)
        if (.not. ok .or. n < 2) return
        do k = 1, size(h%terms)
            w = weight_vector_h2(h, k)
            s = dot_product(w, p)
            if (s <= 0.0_dp) cycle
            allocate(d(n - 1))
            d = w(1:n - 1) - w(n)
            do j = 1, n - 1
                do i = 1, n - 1
                    hess(i,j) = hess(i,j) - h%terms(k)%power * d(i) * d(j) / (s*s)
                end do
            end do
            deallocate(d)
        end do
    end function hessian_independent_h2

    function hessian_independent_h3(h, p_in) result(hess)
        type(hyper3_model), intent(in) :: h
        real(dp), intent(in) :: p_in(:)
        real(dp), allocatable :: hess(:,:), p(:), w(:), d(:)
        real(dp) :: s
        integer :: n, k, i, j
        logical :: ok
        n = h%size()
        allocate(hess(max(0,n-1), max(0,n-1)))
        hess = 0.0_dp
        call normalize_input(p_in, n, p, ok)
        if (.not. ok .or. n < 2) return
        do k = 1, size(h%terms)
            w = weight_vector_h3(h, k)
            s = dot_product(w, p)
            if (s <= 0.0_dp) cycle
            allocate(d(n - 1))
            d = w(1:n - 1) - w(n)
            do j = 1, n - 1
                do i = 1, n - 1
                    hess(i,j) = hess(i,j) - h%terms(k)%power * d(i) * d(j) / (s*s)
                end do
            end do
            deallocate(d)
        end do
    end function hessian_independent_h3

end module hyper2_likelihood
