program test_special
    use goftest_kinds, only : dp
    use goftest_special, only : bessel_k_frac
    implicit none

    call check_close(bessel_k_frac(0.25_dp, 0.1_dp), 2.685156871876059_dp, 3.0e-12_dp, 'K1/4(.1)')
    call check_close(bessel_k_frac(0.75_dp, 0.1_dp), 5.596702511268127_dp, 5.0e-12_dp, 'K3/4(.1)')
    call check_close(bessel_k_frac(0.25_dp, 1.0_dp), 0.4307397744485855_dp, 2.0e-13_dp, 'K1/4(1)')
    call check_close(bessel_k_frac(0.75_dp, 5.0_dp), 0.003886159254974276_dp, 3.0e-14_dp, 'K3/4(5)')

    print '(a)', 'test_special: PASS'

contains

    subroutine check_close(actual, expected, atol, label)
        real(dp), intent(in) :: actual, expected, atol
        character(len=*), intent(in) :: label
        if (abs(actual - expected) > atol) then
            write(*,'(a,2(1x,es24.16))') trim(label)//' mismatch:', actual, expected
            error stop 1
        end if
    end subroutine check_close

end program test_special
