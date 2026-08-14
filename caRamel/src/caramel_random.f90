module caramel_random
    use caramel_kinds, only: dp
    implicit none
    private
    public :: seed_random, random_uniform, random_normal, random_permutation, random_index

contains

    subroutine seed_random(seed)
        integer, intent(in) :: seed
        integer :: n, i
        integer, allocatable :: put(:)

        call random_seed(size=n)
        allocate(put(n))
        do i = 1, n
            put(i) = modulo(seed + 104729 * i + 8191 * i * i, huge(1) - 1)
            if (put(i) <= 0) put(i) = i
        end do
        call random_seed(put=put)
    end subroutine seed_random

    real(dp) function random_uniform() result(u)
        call random_number(u)
    end function random_uniform

    real(dp) function random_normal() result(z)
        real(dp) :: u1, u2
        real(dp), parameter :: twopi = 2.0_dp * acos(-1.0_dp)

        call random_number(u1)
        call random_number(u2)
        u1 = max(u1, tiny(1.0_dp))
        z = sqrt(-2.0_dp * log(u1)) * cos(twopi * u2)
    end function random_normal

    integer function random_index(n) result(idx)
        integer, intent(in) :: n
        real(dp) :: u

        if (n <= 0) then
            idx = 0
            return
        end if
        call random_number(u)
        idx = min(n, 1 + int(u * real(n, dp)))
    end function random_index

    subroutine random_permutation(n, p)
        integer, intent(in) :: n
        integer, intent(out) :: p(n)
        integer :: i, j, tmp

        p = [(i, i=1,n)]
        do i = n, 2, -1
            j = random_index(i)
            tmp = p(i)
            p(i) = p(j)
            p(j) = tmp
        end do
    end subroutine random_permutation

end module caramel_random
