! SPDX-License-Identifier: GPL-2.0-or-later
module dykstra_linalg
    use dykstra_kinds, only : dp
    implicit none
    private
    public :: symmetric_eigen_jacobi

contains

    subroutine symmetric_eigen_jacobi(a, values, vectors, info)
        real(dp), intent(in) :: a(:,:)
        real(dp), allocatable, intent(out) :: values(:)
        real(dp), allocatable, intent(out) :: vectors(:,:)
        integer, intent(out) :: info

        real(dp), allocatable :: work(:,:)
        real(dp) :: app, aqq, apq, tau, t, c, s
        real(dp) :: wip, wiq, vip, viq, off, scale
        integer :: n, i, j, p, q, sweep, max_sweeps

        n = size(a, 1)
        info = 0
        allocate(values(n), vectors(n,n), work(n,n))
        if (size(a, 2) /= n) then
            info = -1
            values = 0.0_dp
            vectors = 0.0_dp
            return
        end if

        work = 0.5_dp * (a + transpose(a))
        vectors = 0.0_dp
        do i = 1, n
            vectors(i,i) = 1.0_dp
        end do

        if (n <= 1) then
            if (n == 1) values(1) = work(1,1)
            return
        end if

        max_sweeps = max(30, 8 * n * n)
        do sweep = 1, max_sweeps
            off = 0.0_dp
            p = 1
            q = 2
            do j = 2, n
                do i = 1, j - 1
                    if (abs(work(i,j)) > off) then
                        off = abs(work(i,j))
                        p = i
                        q = j
                    end if
                end do
            end do
            scale = max(1.0_dp, maxval(abs(work)))
            if (off <= 16.0_dp * epsilon(1.0_dp) * scale) exit

            app = work(p,p)
            aqq = work(q,q)
            apq = work(p,q)
            tau = (aqq - app) / (2.0_dp * apq)
            if (tau >= 0.0_dp) then
                t = 1.0_dp / (tau + sqrt(1.0_dp + tau * tau))
            else
                t = -1.0_dp / (-tau + sqrt(1.0_dp + tau * tau))
            end if
            c = 1.0_dp / sqrt(1.0_dp + t * t)
            s = t * c

            do i = 1, n
                if (i /= p .and. i /= q) then
                    wip = work(i,p)
                    wiq = work(i,q)
                    work(i,p) = c * wip - s * wiq
                    work(p,i) = work(i,p)
                    work(i,q) = s * wip + c * wiq
                    work(q,i) = work(i,q)
                end if
            end do
            work(p,p) = app - t * apq
            work(q,q) = aqq + t * apq
            work(p,q) = 0.0_dp
            work(q,p) = 0.0_dp

            do i = 1, n
                vip = vectors(i,p)
                viq = vectors(i,q)
                vectors(i,p) = c * vip - s * viq
                vectors(i,q) = s * vip + c * viq
            end do
        end do

        if (sweep > max_sweeps) info = 1
        do i = 1, n
            values(i) = work(i,i)
        end do
        call sort_descending(values, vectors)
    end subroutine symmetric_eigen_jacobi

    subroutine sort_descending(values, vectors)
        real(dp), intent(inout) :: values(:)
        real(dp), intent(inout) :: vectors(:,:)

        real(dp), allocatable :: tmp(:)
        real(dp) :: tv
        integer :: i, j, k, n

        n = size(values)
        allocate(tmp(size(vectors,1)))
        do i = 1, n - 1
            k = i
            do j = i + 1, n
                if (values(j) > values(k)) k = j
            end do
            if (k /= i) then
                tv = values(i)
                values(i) = values(k)
                values(k) = tv
                tmp = vectors(:,i)
                vectors(:,i) = vectors(:,k)
                vectors(:,k) = tmp
            end if
        end do
    end subroutine sort_descending

end module dykstra_linalg
