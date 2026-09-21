module fdapace_basis
    use fdapace_kinds, only : dp
    use fdapace_math, only : pi_dp
    implicit none
    private

    public :: create_basis

contains

    function create_basis(k, pts, basis_type) result(res)
        integer, intent(in) :: k !! Positive number of basis functions to evaluate.
        real(dp), intent(in) :: pts(:) !! Coordinates, normally on [0,1], at which the basis is evaluated.
        character(len=*), intent(in) :: basis_type !! Basis family: cos, sin, fourier, legendre01, or poly.
        real(dp), allocatable :: res(:,:)
        real(dp) :: coeff
        real(dp) :: root2
        integer :: i
        integer :: j
        integer :: n

        if (k < 1) error stop "create_basis: k must be positive"
        allocate(res(size(pts), k))
        root2 = sqrt(2.0_dp)
        select case (trim(basis_type))
        case ("cos")
            res(:, 1) = 1.0_dp
            do j = 2, k
                res(:, j) = root2 * cos(real(j - 1, dp) * pi_dp * pts)
            end do
        case ("sin")
            do j = 1, k
                res(:, j) = root2 * sin(real(j, dp) * pi_dp * pts)
            end do
        case ("fourier")
            res(:, 1) = 1.0_dp
            do j = 2, k
                if (mod(j, 2) == 0) then
                    res(:, j) = root2 * sin(real(j, dp) * pi_dp * pts)
                else
                    res(:, j) = root2 * cos(real(j - 1, dp) * pi_dp * pts)
                end if
            end do
        case ("legendre01")
            res = 0.0_dp
            do n = 1, k
                do j = 1, n
                    coeff = (-1.0_dp)**(n - j) * real(binomial_integer(n - 1, j - 1), dp) &
                            * real(binomial_integer(n + j - 2, j - 1), dp) * sqrt(real(2 * n - 1, dp))
                    if (j == 1) then
                        res(:, n) = res(:, n) + coeff
                    else
                        do i = 1, size(pts)
                            res(i, n) = res(i, n) + coeff * pts(i)**(j - 1)
                        end do
                    end if
                end do
            end do
        case ("poly")
            res(:, 1) = 1.0_dp
            do j = 2, k
                res(:, j) = pts**(j - 1)
            end do
        case default
            error stop "create_basis: unknown basis type"
        end select
    end function create_basis

    pure integer function binomial_integer(n, r) result(value)
        integer, intent(in) :: n !! Nonnegative upper argument of the binomial coefficient.
        integer, intent(in) :: r !! Lower argument of the binomial coefficient, between zero and n.
        integer :: i
        integer :: rr

        if (r < 0 .or. r > n) then
            value = 0
            return
        end if
        rr = min(r, n - r)
        value = 1
        do i = 1, rr
            value = value * (n - rr + i) / i
        end do
    end function binomial_integer

end module fdapace_basis
