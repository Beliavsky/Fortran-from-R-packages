! SPDX-License-Identifier: MIT
! Translated from the computational core of the R package trust 0.1-9.
module trust_linalg
    use trust_kinds, only : dp
    implicit none
    private
    public :: symmetric_eigen_jacobi, vector_norm
contains

    pure real(dp) function vector_norm(x) result(v)
        real(dp), intent(in) :: x(:)
        v = sqrt(sum(x * x))
    end function vector_norm

    subroutine symmetric_eigen_jacobi(a, values, vectors, info)
        real(dp), intent(in) :: a(:, :)
        real(dp), intent(out) :: values(:)
        real(dp), intent(out) :: vectors(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: work(:, :)
        real(dp) :: app, aqq, apq, tau, t, c, s, aik, akq, offmax, scale, threshold
        integer :: n, i, k, p, q, sweep, max_sweeps
        logical :: converged

        n = size(a, 1)
        info = 0
        if (size(a, 2) /= n .or. size(values) /= n .or. &
            size(vectors, 1) /= n .or. size(vectors, 2) /= n) then
            info = 1
            return
        end if

        allocate(work(n, n))
        work = 0.5_dp * (a + transpose(a))
        vectors = 0.0_dp
        do i = 1, n
            vectors(i, i) = 1.0_dp
        end do
        if (n == 1) then
            values(1) = work(1, 1)
            return
        end if

        max_sweeps = 60
        converged = .false.
        do sweep = 1, max_sweeps
            offmax = 0.0_dp
            do q = 2, n
                do p = 1, q - 1
                    offmax = max(offmax, abs(work(p, q)))
                end do
            end do
            scale = max(1.0_dp, maxval(abs(work)))
            threshold = 32.0_dp * epsilon(1.0_dp) * scale
            if (offmax <= threshold) then
                converged = .true.
                exit
            end if

            do q = 2, n
                do p = 1, q - 1
                    apq = work(p, q)
                    if (abs(apq) <= threshold) cycle
                    app = work(p, p)
                    aqq = work(q, q)
                    tau = (aqq - app) / (2.0_dp * apq)
                    if (tau >= 0.0_dp) then
                        t = 1.0_dp / (tau + sqrt(1.0_dp + tau * tau))
                    else
                        t = -1.0_dp / (-tau + sqrt(1.0_dp + tau * tau))
                    end if
                    c = 1.0_dp / sqrt(1.0_dp + t * t)
                    s = t * c

                    do k = 1, n
                        if (k /= p .and. k /= q) then
                            aik = work(k, p)
                            akq = work(k, q)
                            work(k, p) = c * aik - s * akq
                            work(p, k) = work(k, p)
                            work(k, q) = s * aik + c * akq
                            work(q, k) = work(k, q)
                        end if
                    end do
                    work(p, p) = c*c*app - 2.0_dp*s*c*apq + s*s*aqq
                    work(q, q) = s*s*app + 2.0_dp*s*c*apq + c*c*aqq
                    work(p, q) = 0.0_dp
                    work(q, p) = 0.0_dp

                    do k = 1, n
                        aik = vectors(k, p)
                        akq = vectors(k, q)
                        vectors(k, p) = c * aik - s * akq
                        vectors(k, q) = s * aik + c * akq
                    end do
                end do
            end do
        end do

        if (.not. converged) then
            offmax = 0.0_dp
            do q = 2, n
                do p = 1, q - 1
                    offmax = max(offmax, abs(work(p, q)))
                end do
            end do
            scale = max(1.0_dp, maxval(abs(work)))
            if (offmax > 64.0_dp * epsilon(1.0_dp) * scale) then
                info = 2
                return
            end if
        end if

        do i = 1, n
            values(i) = work(i, i)
        end do
        call sort_eigenpairs(values, vectors)
    end subroutine symmetric_eigen_jacobi

    subroutine sort_eigenpairs(values, vectors)
        real(dp), intent(inout) :: values(:)
        real(dp), intent(inout) :: vectors(:, :)
        real(dp) :: tmp
        real(dp), allocatable :: col(:)
        integer :: i, j, k, n

        n = size(values)
        allocate(col(n))
        do i = 1, n - 1
            k = i
            do j = i + 1, n
                if (values(j) < values(k)) k = j
            end do
            if (k /= i) then
                tmp = values(i)
                values(i) = values(k)
                values(k) = tmp
                col = vectors(:, i)
                vectors(:, i) = vectors(:, k)
                vectors(:, k) = col
            end if
        end do
    end subroutine sort_eigenpairs

end module trust_linalg
