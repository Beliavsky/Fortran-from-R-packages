program test_fitting
  use betafunctions
  implicit none
  type(beta_params) :: p
  real(dp) :: empty(0)

  p = beta4_fit(empty, 0.5625_dp, 0.006510416666666667_dp, &
                -0.30983866769659335_dp, 2.5854545454545454_dp, .true.)
  call check_close(p%alpha, 5.0_dp, 2.0e-12_dp, 'beta4 alpha')
  call check_close(p%beta, 3.0_dp, 2.0e-12_dp, 'beta4 beta')
  call check_close(p%l, 0.25_dp, 2.0e-12_dp, 'beta4 l')
  call check_close(p%u, 0.75_dp, 2.0e-12_dp, 'beta4 u')

  print '(a)', 'test_fitting: PASS'
contains
  subroutine check_close(got, want, tol, label)
    real(dp), intent(in) :: got, want, tol
    character(*), intent(in) :: label
    if (abs(got - want) > tol) then
      write(*,'(a,2es24.15)') trim(label)//' mismatch: ', got, want
      error stop 1
    end if
  end subroutine check_close
end program test_fitting
