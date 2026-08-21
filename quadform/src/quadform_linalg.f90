module quadform_linalg
    use quadform_kinds, only : dp
    implicit none
    private
    public :: solve_linear

    interface solve_linear
        module procedure solve_linear_real
        module procedure solve_linear_complex
    end interface solve_linear

contains

    function solve_linear_real(a, b, info) result(x)
        real(dp), intent(in) :: a(:, :), b(:, :)
        integer, intent(out), optional :: info
        real(dp), allocatable :: x(:, :)
        real(dp), allocatable :: lu(:, :), rhs(:, :), tmp_row(:)
        real(dp) :: pivabs, scale, factor
        integer :: n, nrhs, i, j, k, p, istat

        n = size(a, 1)
        nrhs = size(b, 2)
        istat = 0
        if (size(a, 2) /= n .or. size(b, 1) /= n) then
            allocate(x(0, 0))
            istat = -1
            if (present(info)) info = istat
            return
        end if

        allocate(lu(n, n), rhs(n, nrhs), x(n, nrhs), tmp_row(max(n, nrhs)))
        lu = a
        rhs = b
        scale = max(1.0_dp, maxval(abs(lu)))

        do k = 1, n
            p = k
            pivabs = abs(lu(k, k))
            do i = k + 1, n
                if (abs(lu(i, k)) > pivabs) then
                    p = i
                    pivabs = abs(lu(i, k))
                end if
            end do
            if (pivabs <= epsilon(1.0_dp) * scale * real(max(1, n), dp)) then
                x = 0.0_dp
                istat = k
                if (present(info)) info = istat
                return
            end if
            if (p /= k) then
                tmp_row(1:n) = lu(k, :)
                lu(k, :) = lu(p, :)
                lu(p, :) = tmp_row(1:n)
                if (nrhs > 0) then
                    tmp_row(1:nrhs) = rhs(k, :)
                    rhs(k, :) = rhs(p, :)
                    rhs(p, :) = tmp_row(1:nrhs)
                end if
            end if

            do i = k + 1, n
                factor = lu(i, k) / lu(k, k)
                lu(i, k) = factor
                lu(i, k + 1:n) = lu(i, k + 1:n) - factor * lu(k, k + 1:n)
                rhs(i, :) = rhs(i, :) - factor * rhs(k, :)
            end do
        end do

        x = rhs
        do j = 1, nrhs
            do i = n, 1, -1
                if (i < n) x(i, j) = x(i, j) - sum(lu(i, i + 1:n) * x(i + 1:n, j))
                x(i, j) = x(i, j) / lu(i, i)
            end do
        end do
        if (present(info)) info = istat
    end function solve_linear_real

    function solve_linear_complex(a, b, info) result(x)
        complex(dp), intent(in) :: a(:, :), b(:, :)
        integer, intent(out), optional :: info
        complex(dp), allocatable :: x(:, :)
        complex(dp), allocatable :: lu(:, :), rhs(:, :), tmp_row(:)
        complex(dp) :: factor
        real(dp) :: pivabs, scale
        integer :: n, nrhs, i, j, k, p, istat

        n = size(a, 1)
        nrhs = size(b, 2)
        istat = 0
        if (size(a, 2) /= n .or. size(b, 1) /= n) then
            allocate(x(0, 0))
            istat = -1
            if (present(info)) info = istat
            return
        end if

        allocate(lu(n, n), rhs(n, nrhs), x(n, nrhs), tmp_row(max(n, nrhs)))
        lu = a
        rhs = b
        scale = max(1.0_dp, maxval(abs(lu)))

        do k = 1, n
            p = k
            pivabs = abs(lu(k, k))
            do i = k + 1, n
                if (abs(lu(i, k)) > pivabs) then
                    p = i
                    pivabs = abs(lu(i, k))
                end if
            end do
            if (pivabs <= epsilon(1.0_dp) * scale * real(max(1, n), dp)) then
                x = cmplx(0.0_dp, 0.0_dp, dp)
                istat = k
                if (present(info)) info = istat
                return
            end if
            if (p /= k) then
                tmp_row(1:n) = lu(k, :)
                lu(k, :) = lu(p, :)
                lu(p, :) = tmp_row(1:n)
                if (nrhs > 0) then
                    tmp_row(1:nrhs) = rhs(k, :)
                    rhs(k, :) = rhs(p, :)
                    rhs(p, :) = tmp_row(1:nrhs)
                end if
            end if

            do i = k + 1, n
                factor = lu(i, k) / lu(k, k)
                lu(i, k) = factor
                lu(i, k + 1:n) = lu(i, k + 1:n) - factor * lu(k, k + 1:n)
                rhs(i, :) = rhs(i, :) - factor * rhs(k, :)
            end do
        end do

        x = rhs
        do j = 1, nrhs
            do i = n, 1, -1
                if (i < n) x(i, j) = x(i, j) - sum(lu(i, i + 1:n) * x(i + 1:n, j))
                x(i, j) = x(i, j) / lu(i, i)
            end do
        end do
        if (present(info)) info = istat
    end function solve_linear_complex

end module quadform_linalg
