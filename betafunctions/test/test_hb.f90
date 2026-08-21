program test_hb
  use betafunctions
  implicit none
  type(beta_params) :: p
  type(classification_result) :: r
  real(dp) :: cuts(1)

  p%alpha = 5.0_dp
  p%beta = 3.0_dp
  p%l = 0.1_dp
  p%u = 0.9_dp
  p%n = 10.0_dp
  p%k = 0.8_dp
  cuts = [5.0_dp]
  call hb_classify_params(p, cuts, r)
  call check_close(r%accuracy_matrix(1,1), 0.1315246348537801_dp, 3.0e-8_dp, 'HB accuracy 11')
  call check_close(r%accuracy_matrix(1,2), 0.08565757356002823_dp, 3.0e-8_dp, 'HB accuracy 12')
  call check_close(r%accuracy_matrix(2,1), 0.09503786514622248_dp, 3.0e-8_dp, 'HB accuracy 21')
  call check_close(r%accuracy_matrix(2,2), 0.6877799264399896_dp, 3.0e-8_dp, 'HB accuracy 22')
  call check_close(sum(r%consistency_matrix), 1.0_dp, 2.0e-10_dp, 'HB consistency normalization')

  print '(a)', 'test_hb: PASS'
contains
  subroutine check_close(got, want, tol, label)
    real(dp), intent(in) :: got, want, tol
    character(*), intent(in) :: label
    if (abs(got - want) > tol) then
      write(*,'(a,2es24.15)') trim(label)//' mismatch: ', got, want
      error stop 1
    end if
  end subroutine check_close
end program test_hb
