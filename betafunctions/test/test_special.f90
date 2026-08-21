program test_special
  use betafunctions
  implicit none
  call check_close(reg_incomplete_beta(0.3_dp, 2.5_dp, 4.0_dp), &
                   0.3521975859068914_dp, 3.0e-13_dp, 'incomplete beta')
  call check_close(inv_reg_incomplete_beta(0.3521975859068914_dp, 2.5_dp, 4.0_dp), &
                   0.3_dp, 2.0e-12_dp, 'inverse beta')
  call check_close(binomial_cdf_le(4, 10, 0.3_dp), 0.8497316674_dp, 3.0e-13_dp, 'binomial cdf')
  call check_close(chi_square_sf(10.0_dp, 5), 0.07523524614651217_dp, 3.0e-12_dp, 'chi-square sf')
  print '(a)', 'test_special: PASS'
contains
  subroutine check_close(got, want, tol, label)
    real(dp), intent(in) :: got, want, tol
    character(*), intent(in) :: label
    if (abs(got - want) > tol) then
      write(*,'(a,2es24.15)') trim(label)//' mismatch: ', got, want
      error stop 1
    end if
  end subroutine check_close
end program test_special
