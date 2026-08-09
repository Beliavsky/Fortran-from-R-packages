! SPDX-License-Identifier: GPL-2.0-only
module nlsr_linalg
    use nlsr_kinds, only : dp
    implicit none
    private
    public :: qr_least_squares, invert_symmetric_positive, solve_linear, covariance_from_jacobian

contains

    subroutine qr_least_squares(a, b, x, qtb, ok)
        real(dp), intent(in) :: a(:,:), b(:)
        real(dp), intent(out) :: x(:), qtb(:)
        logical, intent(out) :: ok
        integer :: m, n, i, j
        real(dp), allocatable :: q(:,:), r(:,:), v(:)
        real(dp) :: t, rn, tol, scale

        m = size(a,1)
        n = size(a,2)
        x = 0.0_dp
        qtb = 0.0_dp
        ok = .false.
        if (size(b) /= m .or. size(x) /= n .or. size(qtb) /= n .or. m < n) return
        allocate(q(m,n), r(n,n), v(m))
        q = 0.0_dp
        r = 0.0_dp
        scale = max(1.0_dp, maxval(abs(a)))
        tol = epsilon(1.0_dp) * scale * real(max(m,n),dp) * 100.0_dp

        do j = 1, n
            v = a(:,j)
            do i = 1, j - 1
                t = dot_product(q(:,i), v)
                r(i,j) = t
                v = v - t*q(:,i)
            end do
            do i = 1, j - 1
                t = dot_product(q(:,i), v)
                r(i,j) = r(i,j) + t
                v = v - t*q(:,i)
            end do
            rn = sqrt(max(0.0_dp, dot_product(v,v)))
            if (rn <= tol) return
            r(j,j) = rn
            q(:,j) = v/rn
        end do
        qtb = matmul(transpose(q), b)
        do i = n, 1, -1
            t = qtb(i)
            if (i < n) t = t - dot_product(r(i,i+1:n), x(i+1:n))
            if (abs(r(i,i)) <= tol) return
            x(i) = t/r(i,i)
        end do
        ok = .true.
    end subroutine qr_least_squares

    subroutine solve_linear(a, b, x, ok)
        real(dp), intent(in) :: a(:,:), b(:)
        real(dp), intent(out) :: x(:)
        logical, intent(out) :: ok
        real(dp), allocatable :: aa(:,:), bb(:), row(:)
        real(dp) :: piv, fac, tol
        integer :: n, i, k, ip

        n = size(a,1)
        x = 0.0_dp
        ok = .false.
        if (size(a,2) /= n .or. size(b) /= n .or. size(x) /= n) return
        allocate(aa(n,n), bb(n), row(n))
        aa = a
        bb = b
        tol = epsilon(1.0_dp)*max(1.0_dp,maxval(abs(a)))*real(n,dp)*100.0_dp
        do k = 1, n
            ip = k
            do i = k + 1, n
                if (abs(aa(i,k)) > abs(aa(ip,k))) ip = i
            end do
            if (abs(aa(ip,k)) <= tol) return
            if (ip /= k) then
                row = aa(k,:); aa(k,:) = aa(ip,:); aa(ip,:) = row
                piv = bb(k); bb(k) = bb(ip); bb(ip) = piv
            end if
            do i = k + 1, n
                fac = aa(i,k)/aa(k,k)
                aa(i,k:n) = aa(i,k:n) - fac*aa(k,k:n)
                bb(i) = bb(i) - fac*bb(k)
            end do
        end do
        do i = n, 1, -1
            piv = bb(i)
            if (i < n) piv = piv - dot_product(aa(i,i+1:n),x(i+1:n))
            x(i) = piv/aa(i,i)
        end do
        ok = .true.
    end subroutine solve_linear

    subroutine invert_symmetric_positive(a, ainv, ok)
        real(dp), intent(in) :: a(:,:)
        real(dp), intent(out) :: ainv(:,:)
        logical, intent(out) :: ok
        real(dp), allocatable :: rhs(:), col(:)
        integer :: n, j
        logical :: one_ok

        n = size(a,1)
        ainv = 0.0_dp
        ok = .false.
        if (size(a,2) /= n .or. size(ainv,1) /= n .or. size(ainv,2) /= n) return
        allocate(rhs(n), col(n))
        do j = 1, n
            rhs = 0.0_dp
            rhs(j) = 1.0_dp
            call solve_linear(a, rhs, col, one_ok)
            if (.not. one_ok) return
            ainv(:,j) = col
        end do
        ainv = 0.5_dp*(ainv + transpose(ainv))
        ok = .true.
    end subroutine invert_symmetric_positive

    subroutine covariance_from_jacobian(jac, ss, nobs, nfree, cov, sigma, ok)
        real(dp), intent(in) :: jac(:,:), ss
        integer, intent(in) :: nobs, nfree
        real(dp), intent(out) :: cov(:,:), sigma
        logical, intent(out) :: ok
        real(dp), allocatable :: xtx(:,:), inv(:,:)
        integer :: p

        p = size(jac,2)
        cov = 0.0_dp
        sigma = huge(1.0_dp)
        ok = .false.
        if (size(cov,1) /= p .or. size(cov,2) /= p .or. nfree < 1) return
        if (nobs <= nfree) return
        allocate(xtx(p,p), inv(p,p))
        xtx = matmul(transpose(jac),jac)
        call invert_symmetric_positive(xtx,inv,ok)
        if (.not. ok) return
        sigma = sqrt(ss/real(nobs-nfree,dp))
        cov = inv*sigma*sigma
    end subroutine covariance_from_jacobian

end module nlsr_linalg
