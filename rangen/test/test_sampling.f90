program test_sampling
    use rangen
    implicit none
    integer, allocatable :: ix(:)
    real(dp), allocatable :: x(:), y(:), m(:,:)
    logical :: seen(10)
    integer :: i

    call seed_all(2026_i8)
    ix = sample_int(10, 7, .false.)
    if (size(ix) /= 7) error stop "sample_int size"
    if (minval(ix) < 1 .or. maxval(ix) > 10) error stop "sample_int range"
    seen = .false.
    do i = 1, size(ix)
        if (seen(ix(i))) error stop "sample_int duplicate without replacement"
        seen(ix(i)) = .true.
    end do

    ix = sample_int(10, 400, .true.)
    if (.not. any(ix == 10)) error stop "replacement sampler never reached upper endpoint"

    x = [(real(i, dp), i=1,10)]
    y = sample_real(x, 5, .false.)
    if (size(y) /= 5) error stop "sample_real size"
    if (minval(y) < 1.0_dp .or. maxval(y) > 10.0_dp) error stop "sample_real range"

    m = reshape([(real(i,dp), i=1,12)], [3,4])
    y = pack(col_sample(m, [2,2,2,2], [.false.,.false.,.false.,.false.]), .true.)
    if (size(y) /= 8) error stop "col_sample dimensions"

    if (nano_time() <= 0.0_dp) error stop "nano_time"
    print *, "test_sampling: PASS"
end program test_sampling
