program test_pcapp_unit
  use pcapp_api, only: dp, qn, cor_fk, scale_adv
  use pcapp_api, only: l1median, l1median_bfgs, l1median_cg, l1median_hocr
  use pcapp_api, only: l1median_nlm, l1median_nm, l1median_vazh, scale_result, median_result
  implicit none

  real(dp) :: x1(5), a1(5), a2(5), b1(7), b2(7), c1(8), c2(8)
  real(dp) :: x(4, 2), medx(5, 2)
  type(scale_result) :: scl
  type(median_result) :: med

  x1 = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
  call assert_close(qn(x1), 1.875177073757389_dp, 1.0e-12_dp, 'Qn finite-sample value')

  a1 = [1.0_dp, 2.0_dp, 3.0_dp, 5.0_dp, 4.0_dp]
  a2 = [1.0_dp, 2.0_dp, 3.0_dp, 3.0_dp, 5.0_dp]
  call assert_close(cor_fk(a1, a2), 0.7378647873726218_dp, 1.0e-12_dp, 'Kendall example A')
  b1 = [8.0_dp, 6.0_dp, 7.0_dp, 5.0_dp, 3.0_dp, 0.0_dp, 9.0_dp]
  b2 = [3.0_dp, 1.0_dp, 4.0_dp, 1.0_dp, 5.0_dp, 9.0_dp, 2.0_dp]
  call assert_close(cor_fk(b1, b2), -0.3903600291794133_dp, 1.0e-12_dp, 'Kendall example B')
  c1 = [1.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 3.0_dp, 4.0_dp, 4.0_dp]
  c2 = [1.0_dp, 2.0_dp, 1.0_dp, 3.0_dp, 3.0_dp, 5.0_dp, 5.0_dp, 5.0_dp]
  call assert_close(cor_fk(c1, c2), 0.8695652173913043_dp, 1.0e-12_dp, 'Kendall example C')

  x(:, 1) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
  x(:, 2) = [10.0_dp, 12.0_dp, 14.0_dp, 16.0_dp]
  scl = scale_adv(x)
  call assert_close(scl%center(1), 2.5_dp, 1.0e-14_dp, 'ScaleAdv center 1')
  call assert_close(scl%center(2), 13.0_dp, 1.0e-14_dp, 'ScaleAdv center 2')
  call assert_close(sum(scl%x(:, 1)), 0.0_dp, 1.0e-13_dp, 'ScaleAdv standardized mean')
  call assert_close(sum(scl%x(:, 1)**2), 3.0_dp, 1.0e-12_dp, 'ScaleAdv standardized sumsq')

  medx(:, 1) = [-1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
  medx(:, 2) = [0.0_dp, 0.0_dp, -1.0_dp, 1.0_dp, 0.0_dp]
  med = l1median(medx)
  call assert_vec_zero(med%par, 1.0e-7_dp, 'l1median')
  med = l1median_bfgs(medx)
  call assert_vec_zero(med%par, 1.0e-7_dp, 'l1median_BFGS')
  med = l1median_cg(medx)
  call assert_vec_zero(med%par, 1.0e-7_dp, 'l1median_CG')
  med = l1median_hocr(medx)
  call assert_vec_zero(med%par, 1.0e-7_dp, 'l1median_HoCr')
  med = l1median_nlm(medx)
  call assert_vec_zero(med%par, 1.0e-7_dp, 'l1median_NLM')
  med = l1median_nm(medx, maxit=500)
  call assert_vec_zero(med%par, 1.0e-5_dp, 'l1median_NM')
  med = l1median_vazh(medx)
  call assert_vec_zero(med%par, 1.0e-7_dp, 'l1median_VaZh')

  print '(a)', 'All pcaPP unit tests passed.'

contains

  subroutine assert_close(actual, expected, tol, label)
    real(dp), intent(in) :: actual !! Value produced by the translated routine.
    real(dp), intent(in) :: expected !! Reference value expected by the deterministic test.
    real(dp), intent(in) :: tol !! Maximum allowed absolute error.
    character(len=*), intent(in) :: label !! Human-readable test label printed on failure.

    if (abs(actual - expected) > tol) then
      write (*, '(a,2es24.14)') 'FAILED '//trim(label)//': ', actual, expected
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_vec_zero(actual, tol, label)
    real(dp), intent(in) :: actual(:) !! Vector expected to be numerically zero.
    real(dp), intent(in) :: tol !! Maximum allowed Euclidean norm.
    character(len=*), intent(in) :: label !! Human-readable test label printed on failure.

    if (sqrt(sum(actual**2)) > tol) then
      write (*, '(a,*(es16.7,1x))') 'FAILED '//trim(label)//': ', actual
      error stop 1
    end if
  end subroutine assert_vec_zero

end program test_pcapp_unit
