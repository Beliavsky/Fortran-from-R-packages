program test_shape
    use ecpdist, only: dp, ecp_shape
    implicit none
    integer :: fails

    fails = 0
    call check_close(ecp_shape(2.0_dp, 0.3_dp, 30.0_dp, 'bowley'), &
        0.19715028103014534_dp, 5.0e-12_dp, 'bowley', fails)
    call check_close(ecp_shape(2.0_dp, 0.3_dp, 30.0_dp, 'moors'), &
        0.54159658342916764_dp, 5.0e-12_dp, 'moors', fails)

    if (fails /= 0) error stop 1
    print '(a)', 'test_shape: PASS'

contains

    subroutine check_close(got, expected, atol, label, nfail)
        real(dp), intent(in) :: got, expected, atol
        character(len=*), intent(in) :: label
        integer, intent(inout) :: nfail
        if (abs(got - expected) > atol) then
            print '(a,2es24.15)', trim(label)//': ', got, expected
            nfail = nfail + 1
        end if
    end subroutine check_close

end program test_shape
