program test_classification
  use betafunctions
  implicit none
  type(beta_params) :: p
  type(classification_result) :: r
  type(accuracy_stats) :: a
  real(dp) :: cuts(1), truecuts(1)

  p%alpha = 5.0_dp
  p%beta = 3.0_dp
  p%l = 0.1_dp
  p%u = 0.9_dp
  p%etl = 20.0_dp
  cuts = [0.5_dp]
  truecuts = [0.5_dp]
  call ll_classify_params(p, cuts, r, truecuts)

  call check_close(r%accuracy_matrix(1,1), 0.15132528760375122_dp, 2.0e-8_dp, 'LL accuracy 11')
  call check_close(r%accuracy_matrix(1,2), 0.08005989025735502_dp, 2.0e-8_dp, 'LL accuracy 12')
  call check_close(r%accuracy_matrix(2,1), 0.07523721239626660_dp, 2.0e-8_dp, 'LL accuracy 21')
  call check_close(r%accuracy_matrix(2,2), 0.6933776097426485_dp, 2.0e-8_dp, 'LL accuracy 22')
  call check_close(r%overall_accuracy, 0.8447028973463817_dp, 3.0e-8_dp, 'LL overall accuracy')
  call check_close(r%overall_consistency%p, 0.7889732242630928_dp, 3.0e-8_dp, 'LL consistency p')
  call check_close(r%overall_consistency%p_c, 0.6443078453454085_dp, 3.0e-8_dp, 'LL consistency pc')
  call check_close(r%overall_consistency%kappa, 0.40671512436974366_dp, 4.0e-8_dp, 'LL kappa')

  a = ca_stats(10.0_dp, 20.0_dp, 5.0_dp, 4.0_dp)
  call check_close(a%sensitivity, 10.0_dp/14.0_dp, 1.0e-15_dp, 'sensitivity')
  call check_close(a%specificity, 0.8_dp, 1.0e-15_dp, 'specificity')

  call check_close(etl(60.0_dp, 100.0_dp, 0.0_dp, 100.0_dp, 0.8_dp), 116.0_dp, 1.0e-13_dp, 'ETL')
  call check_close(reliability_from_etl(60.0_dp, 100.0_dp, 0.0_dp, 100.0_dp, 116.0_dp), 0.8_dp, 1.0e-13_dp, 'R.ETL')

  print '(a)', 'test_classification: PASS'
contains
  subroutine check_close(got, want, tol, label)
    real(dp), intent(in) :: got, want, tol
    character(*), intent(in) :: label
    if (abs(got - want) > tol) then
      write(*,'(a,2es24.15)') trim(label)//' mismatch: ', got, want
      error stop 1
    end if
  end subroutine check_close
end program test_classification
