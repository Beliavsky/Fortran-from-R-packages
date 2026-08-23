program test_pcg
    use rangen_kinds, only : i8
    use rangen_pcg32, only : pcg32_state
    implicit none
    type(pcg32_state) :: rng
    integer(i8) :: got(5)
    integer(i8), parameter :: expected(5) = [0_i8, 1971522493_i8, 242089394_i8, 3457789919_i8, 3637502659_i8]
    integer :: i

    call rng%seed(42_i8, 1442695040888963407_i8)
    do i = 1, 5
        got(i) = rng%next_uint32()
    end do
    if (any(got /= expected)) error stop "PCG32 reference sequence mismatch"
    print *, "test_pcg: PASS"
end program test_pcg
