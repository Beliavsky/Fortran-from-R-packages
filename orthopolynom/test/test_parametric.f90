program test_parametric
  use orthopolynom
  implicit none
  integer :: failures, j
  type(polylist_t) :: p
  real(dp), parameter :: x = 0.23_dp
  real(dp), parameter :: gegen_ref(0:4) = [ &
    1.0_dp, 0.368_dp, -0.647648_dp, -0.596990208_dp, 0.321998479104_dp ]
  real(dp), parameter :: jacp_ref(0:4) = [ &
    1.0_dp, 0.1065_dp, -0.585857625_dp, -0.2991931119375_dp, 0.34185519953552324_dp ]
  real(dp), parameter :: glag_ref(0:4) = [ &
    1.0_dp, 1.47_dp, 1.70045_dp, 1.777487166666667_dp, 1.74654103375_dp ]
  real(dp), parameter :: jacg_ref(0:4) = [ &
    1.0_dp, -0.0124242424242424_dp, -0.04014080737165426_dp, &
    0.010370463427242174_dp, -0.0002897976356686171_dp ]
  real(dp), parameter :: gherm_ref(0:4) = [ &
    1.0_dp, 0.46_dp, -4.5884_dp, -3.950664_dp, 38.56061456_dp ]

  failures = 0

  p = gegenbauer_polynomials(4, 0.8_dp)
  do j = 0, 4
    call check_close(p%item(j + 1)%evaluate(x), gegen_ref(j), 2.0e-12_dp, &
      'Gegenbauer independent reference', failures)
  end do

  p = jacobi_p_polynomials(4, 0.3_dp, 0.8_dp)
  do j = 0, 4
    call check_close(p%item(j + 1)%evaluate(x), jacp_ref(j), 2.0e-12_dp, &
      'Jacobi P independent reference', failures)
  end do

  p = glaguerre_polynomials(4, 0.7_dp)
  do j = 0, 4
    call check_close(p%item(j + 1)%evaluate(x), glag_ref(j), 2.0e-12_dp, &
      'generalized Laguerre independent reference', failures)
  end do

  p = jacobi_g_polynomials(4, 2.3_dp, 0.8_dp)
  do j = 0, 4
    call check_close(p%item(j + 1)%evaluate(x), jacg_ref(j), 2.0e-12_dp, &
      'Jacobi G independent reference', failures)
  end do

  p = ghermite_h_polynomials(4, 0.7_dp)
  do j = 0, 4
    call check_close(p%item(j + 1)%evaluate(x), gherm_ref(j), 2.0e-11_dp, &
      'generalized Hermite independent reference', failures)
  end do

  call check_normalization(gegenbauer_polynomials(5, 0.8_dp), &
    gegenbauer_polynomials(5, 0.8_dp, .true.), gegenbauer_inner_products(5, 0.8_dp), &
    'Gegenbauer normalized', failures)
  call check_normalization(ghermite_h_polynomials(5, 0.7_dp), &
    ghermite_h_polynomials(5, 0.7_dp, .true.), ghermite_h_inner_products(5, 0.7_dp), &
    'generalized Hermite normalized', failures)
  call check_normalization(glaguerre_polynomials(5, 0.7_dp), &
    glaguerre_polynomials(5, 0.7_dp, .true.), glaguerre_inner_products(5, 0.7_dp), &
    'generalized Laguerre normalized', failures)
  call check_normalization(jacobi_p_polynomials(5, 0.3_dp, 0.8_dp), &
    jacobi_p_polynomials(5, 0.3_dp, 0.8_dp, .true.), jacobi_p_inner_products(5, 0.3_dp, 0.8_dp), &
    'Jacobi P normalized', failures)
  call check_normalization(jacobi_g_polynomials(5, 2.3_dp, 0.8_dp), &
    jacobi_g_polynomials(5, 2.3_dp, 0.8_dp, .true.), jacobi_g_inner_products(5, 2.3_dp, 0.8_dp), &
    'Jacobi G normalized', failures)
  call check_normalization(jacobi_g_polynomials(5, 0.0_dp, 0.3_dp), &
    jacobi_g_polynomials(5, 0.0_dp, 0.3_dp, .true.), jacobi_g_inner_products(5, 0.0_dp, 0.3_dp), &
    'Jacobi G p=0 normalized branch', failures)

  ! Exercise special parameter branches and n=0 boundaries.
  p = gegenbauer_polynomials(0, 0.0_dp)
  call check_close(p%item(1)%evaluate(x), 1.0_dp, 1.0e-13_dp, 'Gegenbauer alpha=0 n=0', failures)
  p = jacobi_g_polynomials(0, 0.0_dp, 0.3_dp)
  call check_close(p%item(1)%evaluate(x), 1.0_dp, 1.0e-13_dp, 'Jacobi G p=0 n=0', failures)

  if (failures /= 0) then
    print '(a,i0)', 'test_parametric: FAIL ', failures
    error stop 1
  end if
  print '(a)', 'test_parametric: PASS'

contains

  subroutine check_normalization(un, norm, h, label, failures)
    type(polylist_t), intent(in) :: un, norm
    real(dp), intent(in) :: h(:)
    character(*), intent(in) :: label
    integer, intent(inout) :: failures
    integer :: k
    real(dp), allocatable :: a(:), b(:)

    do k = 1, un%size()
      a = un%item(k)%coefficients()
      b = norm%item(k)%coefficients()
      if (size(a) /= size(b)) then
        failures = failures + 1
        print '(a,a)', 'normalization size mismatch: ', label
      else if (maxval(abs(b - a / sqrt(h(k)))) > 2.0e-10_dp) then
        failures = failures + 1
        print '(a,a,1x,i0)', 'normalization mismatch: ', label, k - 1
      end if
    end do
  end subroutine check_normalization

  subroutine check_close(got, expected, tol, label, failures)
    real(dp), intent(in) :: got, expected, tol
    character(*), intent(in) :: label
    integer, intent(inout) :: failures
    if (abs(got - expected) > tol) then
      failures = failures + 1
      print '(a,a,2(1x,es16.8))', 'value mismatch: ', label, got, expected
    end if
  end subroutine check_close

end program test_parametric
