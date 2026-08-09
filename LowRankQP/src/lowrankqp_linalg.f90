! SPDX-License-Identifier: GPL-2.0-or-later
module lowrankqp_linalg
    use lowrankqp_kinds, only : dp
    implicit none
    private
    public :: chol_factor, chol_solve, lu_factor, lu_solve

contains

subroutine chol_factor(a, info)
    real(dp), intent(inout) :: a(:,:)
    integer, intent(out) :: info
    integer :: i, j, k, n
    real(dp) :: s

    n = size(a,1)
    info = 0
    if (size(a,2) /= n) then
        info = -1
        return
    end if
    do j = 1, n
        s = a(j,j)
        do k = 1, j - 1
            s = s - a(j,k)*a(j,k)
        end do
        if (s <= 0.0_dp .or. .not. (s < huge(s))) then
            info = j
            return
        end if
        a(j,j) = sqrt(s)
        do i = j + 1, n
            s = a(i,j)
            do k = 1, j - 1
                s = s - a(i,k)*a(j,k)
            end do
            a(i,j) = s/a(j,j)
        end do
        if (j < n) a(j,j+1:n) = 0.0_dp
    end do
end subroutine chol_factor

subroutine chol_solve(l, b, info)
    real(dp), intent(in) :: l(:,:)
    real(dp), intent(inout) :: b(:,:)
    integer, intent(out) :: info
    integer :: i, k, n, nrhs
    real(dp) :: s

    n = size(l,1)
    nrhs = size(b,2)
    info = 0
    if (size(l,2) /= n .or. size(b,1) /= n) then
        info = -1
        return
    end if
    do k = 1, nrhs
        do i = 1, n
            s = b(i,k)
            if (i > 1) s = s - dot_product(l(i,1:i-1), b(1:i-1,k))
            if (abs(l(i,i)) <= tiny(1.0_dp)) then
                info = i
                return
            end if
            b(i,k) = s/l(i,i)
        end do
        do i = n, 1, -1
            s = b(i,k)
            if (i < n) s = s - dot_product(l(i+1:n,i), b(i+1:n,k))
            if (abs(l(i,i)) <= tiny(1.0_dp)) then
                info = i
                return
            end if
            b(i,k) = s/l(i,i)
        end do
    end do
end subroutine chol_solve

subroutine lu_factor(a, piv, info)
    real(dp), intent(inout) :: a(:,:)
    integer, intent(out) :: piv(:)
    integer, intent(out) :: info
    integer :: i, j, k, n, p
    real(dp) :: vmax
    real(dp), allocatable :: row(:)

    n = size(a,1)
    info = 0
    if (size(a,2) /= n .or. size(piv) /= n) then
        info = -1
        return
    end if
    allocate(row(n))
    do k = 1, n
        p = k
        vmax = abs(a(k,k))
        do i = k + 1, n
            if (abs(a(i,k)) > vmax) then
                vmax = abs(a(i,k))
                p = i
            end if
        end do
        piv(k) = p
        if (vmax <= tiny(1.0_dp)) then
            info = k
            return
        end if
        if (p /= k) then
            row = a(k,:)
            a(k,:) = a(p,:)
            a(p,:) = row
        end if
        do i = k + 1, n
            a(i,k) = a(i,k)/a(k,k)
            do j = k + 1, n
                a(i,j) = a(i,j) - a(i,k)*a(k,j)
            end do
        end do
    end do
end subroutine lu_factor

subroutine lu_solve(lu, piv, b, info)
    real(dp), intent(in) :: lu(:,:)
    integer, intent(in) :: piv(:)
    real(dp), intent(inout) :: b(:,:)
    integer, intent(out) :: info
    integer :: i, k, n, nrhs, p
    real(dp) :: s, tmp

    n = size(lu,1)
    nrhs = size(b,2)
    info = 0
    if (size(lu,2) /= n .or. size(piv) /= n .or. size(b,1) /= n) then
        info = -1
        return
    end if
    do k = 1, n
        p = piv(k)
        if (p /= k) then
            do i = 1, nrhs
                tmp = b(k,i)
                b(k,i) = b(p,i)
                b(p,i) = tmp
            end do
        end if
    end do
    do k = 1, nrhs
        do i = 2, n
            b(i,k) = b(i,k) - dot_product(lu(i,1:i-1), b(1:i-1,k))
        end do
        do i = n, 1, -1
            s = b(i,k)
            if (i < n) s = s - dot_product(lu(i,i+1:n), b(i+1:n,k))
            if (abs(lu(i,i)) <= tiny(1.0_dp)) then
                info = i
                return
            end if
            b(i,k) = s/lu(i,i)
        end do
    end do
end subroutine lu_solve

end module lowrankqp_linalg
