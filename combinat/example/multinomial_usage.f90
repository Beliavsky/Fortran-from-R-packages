program multinomial_usage
    use, intrinsic :: iso_fortran_env, only : int64
    use combinat
    implicit none

    integer, allocatable :: draws(:, :)
    type(combinat_rng_state) :: rng
    integer :: info
    integer :: i

    call rng_seed(rng, 123_int64)
    call rmultz2([5, 10, 20], [0.2_dp, 0.3_dp, 0.5_dp], rng, draws, info=info)
    if (info /= combinat_success) error stop 'rmultz2 failed'

    do i = 1, size(draws, 2)
        print '(a, i0, a, *(1x, i0))', 'draw ', i, ':', draws(:, i)
    end do
end program multinomial_usage
