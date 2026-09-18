module tclust_rng
    use tclust_kinds, only : dp
    implicit none
    private

    real(dp), parameter :: pi_dp = acos(-1.0_dp)

    public :: set_tclust_seed
    public :: rand_normal
    public :: rand_uniform
    public :: sample_index
    public :: random_orthonormal

contains

    subroutine set_tclust_seed(seed)
        integer, intent(in) :: seed !! Deterministic scalar seed used to initialize the intrinsic RNG.
        integer :: n
        integer :: i
        integer, allocatable :: put(:)

        call random_seed(size=n)
        allocate(put(n))
        do i = 1, n
            put(i) = modulo(seed + 104729 * i, huge(1) - 1)
            if (put(i) <= 0) put(i) = i
        end do
        call random_seed(put=put)
    end subroutine set_tclust_seed

    real(dp) function rand_uniform(a, b) result(x)
        real(dp), intent(in) :: a !! Lower endpoint of the requested uniform interval.
        real(dp), intent(in) :: b !! Upper endpoint of the requested uniform interval.
        real(dp) :: u

        call random_number(u)
        x = a + (b - a) * u
    end function rand_uniform

    real(dp) function rand_normal() result(x)
        real(dp) :: u1
        real(dp) :: u2

        call random_number(u1)
        call random_number(u2)
        u1 = max(u1, tiny(1.0_dp))
        x = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * pi_dp * u2)
    end function rand_normal

    integer function sample_index(n) result(idx)
        integer, intent(in) :: n !! Number of equally likely indices; valid range is 1 through n.
        real(dp) :: u

        call random_number(u)
        idx = min(n, 1 + int(u * real(n, dp)))
    end function sample_index

    subroutine random_orthonormal(q)
        real(dp), intent(out) :: q(:, :) !! Square output matrix whose columns are random orthonormal vectors.
        integer :: i
        integer :: j
        integer :: n
        real(dp) :: nrm

        n = size(q, 1)
        if (size(q, 2) /= n) error stop 'random_orthonormal: q must be square'
        do j = 1, n
            do i = 1, n
                q(i, j) = rand_normal()
            end do
            do i = 1, j - 1
                q(:, j) = q(:, j) - dot_product(q(:, i), q(:, j)) * q(:, i)
            end do
            nrm = sqrt(sum(q(:, j) ** 2))
            if (nrm <= sqrt(tiny(1.0_dp))) then
                q(:, j) = 0.0_dp
                q(j, j) = 1.0_dp
            else
                q(:, j) = q(:, j) / nrm
            end if
        end do
    end subroutine random_orthonormal

end module tclust_rng
