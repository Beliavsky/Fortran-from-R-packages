! SPDX-License-Identifier: GPL-2.0-or-later
module ltsa_arma
    use ltsa_kinds, only : dp
    use ltsa_status, only : ltsa_error, ltsa_success, ltsa_invalid_input, ltsa_nonstationary, set_error
    use ltsa_linalg, only : solve_linear
    implicit none
    private

    public :: tacvf_arma, ar_to_ma, ar_is_stationary

contains

    logical function ar_is_stationary(phi) result(stationary)
        real(dp), intent(in) :: phi(:)
        real(dp), allocatable :: work(:), reduced(:), reflection(:)
        real(dp) :: a, denominator
        integer :: k, l, p
        p = size(phi)
        stationary = .true.
        if (p == 0) return
        allocate(work(p), reflection(p))
        work = phi
        do k = 1, p
            l = p+1-k
            a = work(l)
            reflection(k) = a
            if (abs(a) >= 1.0_dp .or. a /= a) then
                stationary = .false.
                return
            end if
            if (l > 1) then
                denominator = 1.0_dp-a*a
                allocate(reduced(l-1))
                reduced = (work(1:l-1)+a*work(l-1:1:-1))/denominator
                work(1:l-1) = reduced
                deallocate(reduced)
            end if
        end do
        stationary = all(abs(reflection) < 1.0_dp)
    end function ar_is_stationary

    subroutine tacvf_arma(phi, theta, max_lag, sigma2, acvf, error)
        real(dp), intent(in) :: phi(:), theta(:)
        integer, intent(in) :: max_lag
        real(dp), intent(in) :: sigma2
        real(dp), allocatable, intent(out) :: acvf(:)
        type(ltsa_error), intent(out) :: error
        integer :: i, j, k, p, q, order, max_lag_p1
        real(dp), allocatable :: b(:), c(:), theta2(:), phi2(:), a(:,:), rhs(:), solution(:), g(:)
        error%code = ltsa_success
        error%message = ''
        if (max_lag < 0 .or. sigma2 <= 0.0_dp) then
            allocate(acvf(0))
            call set_error(error, ltsa_invalid_input, 'max_lag must be nonnegative and sigma2 positive')
            return
        end if
        if (.not. ar_is_stationary(phi)) then
            allocate(acvf(0))
            call set_error(error, ltsa_nonstationary, 'AR polynomial is not stationary-causal')
            return
        end if
        p = size(phi)
        q = size(theta)
        max_lag_p1 = max_lag+1
        allocate(acvf(max_lag_p1), source=0.0_dp)
        if (max(p,q) == 0) then
            acvf(1) = sigma2
            return
        end if
        order = max(p,q)+1
        allocate(b(order), source=0.0_dp)
        allocate(c(q+1), source=0.0_dp)
        allocate(theta2(q+1), source=0.0_dp)
        allocate(phi2(3*order), source=0.0_dp)
        c(1) = 1.0_dp
        theta2(1) = -1.0_dp
        if (q > 0) theta2(2:) = theta
        phi2(order) = -1.0_dp
        if (p > 0) phi2(order+1:order+p) = phi
        do k = 1, q
            c(k+1) = -theta(k)
            do i = 1, min(p,k)
                c(k+1) = c(k+1)+phi(i)*c(k+1-i)
            end do
        end do
        do k = 0, q
            do i = k, q
                b(k+1) = b(k+1)-theta2(i+1)*c(i-k+1)
            end do
        end do
        if (p == 0) then
            acvf(1:min(max_lag_p1,order)) = sigma2*b(1:min(max_lag_p1,order))
            return
        end if
        allocate(a(order,order), source=0.0_dp)
        allocate(rhs(order))
        do i = 1, order
            do j = 1, order
                if (j == 1) then
                    a(i,j) = phi2(order+i-1)
                else
                    a(i,j) = phi2(order+i-j)+phi2(order+i+j-2)
                end if
            end do
        end do
        rhs = -b
        call solve_linear(a, rhs, solution, error)
        if (.not. error%ok()) then
            deallocate(acvf)
            allocate(acvf(0))
            return
        end if
        allocate(g(max(max_lag_p1,order)), source=0.0_dp)
        g(1:order) = solution
        do i = order+1, size(g)
            g(i) = 0.0_dp
            do j = 1, p
                g(i) = g(i)+phi(j)*g(i-j)
            end do
        end do
        acvf = sigma2*g(1:max_lag_p1)
    end subroutine tacvf_arma

    function ar_to_ma(phi, max_lag) result(psi)
        real(dp), intent(in) :: phi(:)
        integer, intent(in) :: max_lag
        real(dp), allocatable :: psi(:)
        integer :: j, k, p
        p = size(phi)
        allocate(psi(max_lag+1), source=0.0_dp)
        psi(1) = 1.0_dp
        do k = 1, max_lag
            do j = 1, min(k,p)
                psi(k+1) = psi(k+1)+phi(j)*psi(k-j+1)
            end do
        end do
    end function ar_to_ma

end module ltsa_arma
