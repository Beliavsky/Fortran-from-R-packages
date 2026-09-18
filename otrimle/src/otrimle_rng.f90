module otrimle_rng
    use otrimle_kinds, only : dp
    implicit none
    private

    real(dp), parameter :: pi_dp = acos(-1.0_dp)

    public :: set_otrimle_seed
    public :: rand_normal
    public :: sample_weighted_index

contains

    subroutine set_otrimle_seed(seed)
        integer, intent(in) :: seed !! Scalar seed used to initialize the intrinsic Fortran random-number generator.
        integer :: i
        integer :: n
        integer, allocatable :: put(:)

        call random_seed(size=n)
        allocate(put(n))
        do i = 1, n
            put(i) = modulo(seed + 104729 * i, huge(1) - 1)
            if (put(i) <= 0) put(i) = i
        end do
        call random_seed(put=put)
    end subroutine set_otrimle_seed

    real(dp) function rand_normal() result(x)
        real(dp) :: u1
        real(dp) :: u2

        call random_number(u1)
        call random_number(u2)
        u1 = max(u1, tiny(1.0_dp))
        x = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * pi_dp * u2)
    end function rand_normal

    integer function sample_weighted_index(weights) result(idx)
        real(dp), intent(in) :: weights(:) !! Nonnegative sampling weights; they need not be normalized.
        real(dp) :: c
        real(dp) :: sw
        real(dp) :: u
        integer :: i

        sw = sum(max(weights, 0.0_dp))
        if (sw <= 0.0_dp) then
            idx = 1
            return
        end if
        call random_number(u)
        u = u * sw
        c = 0.0_dp
        idx = size(weights)
        do i = 1, size(weights)
            c = c + max(weights(i), 0.0_dp)
            if (u <= c) then
                idx = i
                return
            end if
        end do
    end function sample_weighted_index

end module otrimle_rng
