! SPDX-License-Identifier: GPL-3.0-or-later
module ecos_linalg
    use ecos_types, only : dp
    implicit none
    private
    public :: solve_linear, solve_kkt, least_norm_equalities, vecnorm2

contains

    pure real(dp) function vecnorm2(x) result(v)
        real(dp), intent(in) :: x(:)
        v = sqrt(max(0.0_dp, dot_product(x,x)))
    end function vecnorm2

    subroutine solve_linear(a, b, x, info)
        real(dp), intent(in) :: a(:,:), b(:)
        real(dp), intent(out) :: x(:)
        integer, intent(out) :: info
        real(dp), allocatable :: aa(:,:), bb(:)
        real(dp) :: piv, fac, tmp, scale
        integer :: n, i, j, k, ip
        n = size(b)
        info = 0
        if (size(a,1) /= n .or. size(a,2) /= n .or. size(x) /= n) then
            info = -1
            x = 0.0_dp
            return
        end if
        if (n == 0) return
        allocate(aa(n,n), bb(n))
        aa = a
        bb = b
        scale = max(1.0_dp, maxval(abs(aa)))
        do k = 1, n-1
            ip = k - 1 + maxloc(abs(aa(k:n,k)), dim=1)
            if (abs(aa(ip,k)) <= 100.0_dp*epsilon(1.0_dp)*scale) then
                info = k
                x = 0.0_dp
                return
            end if
            if (ip /= k) then
                do j = k, n
                    tmp = aa(k,j)
                    aa(k,j) = aa(ip,j)
                    aa(ip,j) = tmp
                end do
                tmp = bb(k)
                bb(k) = bb(ip)
                bb(ip) = tmp
            end if
            piv = aa(k,k)
            do i = k+1, n
                fac = aa(i,k)/piv
                aa(i,k) = 0.0_dp
                aa(i,k+1:n) = aa(i,k+1:n) - fac*aa(k,k+1:n)
                bb(i) = bb(i) - fac*bb(k)
            end do
        end do
        if (abs(aa(n,n)) <= 100.0_dp*epsilon(1.0_dp)*scale) then
            info = n
            x = 0.0_dp
            return
        end if
        x(n) = bb(n)/aa(n,n)
        do i = n-1, 1, -1
            x(i) = (bb(i) - dot_product(aa(i,i+1:n),x(i+1:n)))/aa(i,i)
        end do
    end subroutine solve_linear

    subroutine solve_kkt(h, a, rhsx, rhsy, dx, dy, reg, info)
        real(dp), intent(in) :: h(:,:), a(:,:), rhsx(:), rhsy(:), reg
        real(dp), intent(out) :: dx(:), dy(:)
        integer, intent(out) :: info
        real(dp), allocatable :: kkt(:,:), rhs(:), sol(:)
        integer :: n, p, i
        n = size(rhsx)
        p = size(rhsy)
        if (size(h,1) /= n .or. size(h,2) /= n .or. size(a,2) /= n .or. &
            size(a,1) /= p .or. size(dx) /= n .or. size(dy) /= p) then
            info = -1
            dx = 0.0_dp
            dy = 0.0_dp
            return
        end if
        if (p == 0) then
            allocate(kkt(n,n), rhs(n), sol(n))
            kkt = h
            do i = 1, n
                kkt(i,i) = kkt(i,i) + reg
            end do
            rhs = rhsx
            call solve_linear(kkt,rhs,sol,info)
            dx = sol
            return
        end if
        allocate(kkt(n+p,n+p), rhs(n+p), sol(n+p))
        kkt = 0.0_dp
        kkt(1:n,1:n) = h
        do i = 1, n
            kkt(i,i) = kkt(i,i) + reg
        end do
        kkt(1:n,n+1:n+p) = transpose(a)
        kkt(n+1:n+p,1:n) = a
        rhs(1:n) = rhsx
        rhs(n+1:n+p) = rhsy
        call solve_linear(kkt,rhs,sol,info)
        if (info /= 0) then
            ! Increase primal regularization and add tiny dual regularization.
            kkt = 0.0_dp
            kkt(1:n,1:n) = h
            do i = 1, n
                kkt(i,i) = kkt(i,i) + max(1.0e-7_dp,100.0_dp*reg)
            end do
            kkt(1:n,n+1:n+p) = transpose(a)
            kkt(n+1:n+p,1:n) = a
            do i = 1, p
                kkt(n+i,n+i) = -max(1.0e-10_dp,reg)
            end do
            call solve_linear(kkt,rhs,sol,info)
        end if
        if (info == 0) then
            dx = sol(1:n)
            dy = sol(n+1:n+p)
        else
            dx = 0.0_dp
            dy = 0.0_dp
        end if
    end subroutine solve_kkt

    subroutine least_norm_equalities(a, b, x, residual, info)
        real(dp), intent(in) :: a(:,:), b(:)
        real(dp), intent(out) :: x(:)
        real(dp), intent(out) :: residual
        integer, intent(out) :: info
        real(dp), allocatable :: gram(:,:), y(:), rhs(:)
        integer :: p, n, i
        p = size(a,1)
        n = size(a,2)
        info = 0
        x = 0.0_dp
        residual = 0.0_dp
        if (size(b) /= p .or. size(x) /= n) then
            info = -1
            return
        end if
        if (p == 0) return
        allocate(gram(p,p), y(p), rhs(p))
        gram = matmul(a,transpose(a))
        rhs = b
        do i = 1, p
            gram(i,i) = gram(i,i) + 1.0e-12_dp*max(1.0_dp,maxval(abs(gram)))
        end do
        call solve_linear(gram,rhs,y,info)
        if (info /= 0) return
        x = matmul(transpose(a),y)
        residual = vecnorm2(matmul(a,x)-b)
    end subroutine least_norm_equalities

end module ecos_linalg
