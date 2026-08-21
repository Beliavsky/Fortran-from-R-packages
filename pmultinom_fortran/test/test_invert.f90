program test_invert
    use pmultinom_module, only : dp, invert_pmultinom
    implicit none

    integer :: failures

    failures = 0
    call check_int('one lower 0', invert_pmultinom([1.0_dp], 0.95_dp, lower=[0.0_dp]), 1)
    call check_int('one lower 1', invert_pmultinom([1.0_dp], 0.95_dp, lower=[1.0_dp]), 2)
    call check_int('one upper 0', invert_pmultinom([1.0_dp], 0.95_dp, upper=[0.0_dp]), 1)
    call check_int('three equal lower', invert_pmultinom([1.0_dp/3.0_dp, 1.0_dp/3.0_dp, 1.0_dp/3.0_dp], &
                   0.95_dp, lower=[0.0_dp, 0.0_dp, 0.0_dp]), 11)

    if (failures /= 0) then
        print '(a,i0)', 'test_invert: FAIL ', failures
        error stop 1
    end if
    print '(a)', 'test_invert: PASS'

contains

    subroutine check_int(label, got, expected)
        character(len=*), intent(in) :: label
        integer, intent(in) :: got, expected
        if (got /= expected) then
            failures = failures + 1
            print '(a,2(1x,i0))', trim(label)//' got/expected:', got, expected
        end if
    end subroutine check_int

end program test_invert
