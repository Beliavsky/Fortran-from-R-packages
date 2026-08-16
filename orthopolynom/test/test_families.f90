program test_families
  use orthopolynom
  implicit none
  integer :: failures
  type(polylist_t) :: p
  type(recurrence_t) :: r
  real(dp), allocatable :: h(:)

  failures = 0

  p = chebyshev_t_polynomials(3)
  call check_coef(p%item(1), [1.0_dp], 1.0e-13_dp, 'T0', failures)
  call check_coef(p%item(2), [0.0_dp, 1.0_dp], 1.0e-13_dp, 'T1', failures)
  call check_coef(p%item(3), [-1.0_dp, 0.0_dp, 2.0_dp], 1.0e-13_dp, 'T2', failures)
  call check_coef(p%item(4), [0.0_dp, -3.0_dp, 0.0_dp, 4.0_dp], 1.0e-13_dp, 'T3', failures)

  p = chebyshev_u_polynomials(3)
  call check_coef(p%item(4), [0.0_dp, -4.0_dp, 0.0_dp, 8.0_dp], 1.0e-13_dp, 'U3', failures)

  p = chebyshev_c_polynomials(3)
  call check_coef(p%item(2), [0.0_dp, 0.5_dp], 1.0e-13_dp, 'C1', failures)
  call check_coef(p%item(3), [-1.0_dp, 0.0_dp, 0.5_dp], 1.0e-13_dp, 'C2', failures)

  p = chebyshev_s_polynomials(3)
  call check_coef(p%item(3), [-1.0_dp, 0.0_dp, 1.0_dp], 1.0e-13_dp, 'S2', failures)

  p = schebyshev_t_polynomials(3)
  call check_coef(p%item(2), [-1.0_dp, 2.0_dp], 1.0e-13_dp, 'shifted T1', failures)
  call check_coef(p%item(3), [1.0_dp, -8.0_dp, 8.0_dp], 1.0e-13_dp, 'shifted T2', failures)

  p = schebyshev_u_polynomials(2)
  call check_coef(p%item(2), [-2.0_dp, 4.0_dp], 1.0e-13_dp, 'shifted U1', failures)
  call check_coef(p%item(3), [3.0_dp, -16.0_dp, 16.0_dp], 1.0e-13_dp, 'shifted U2', failures)

  p = legendre_polynomials(3)
  call check_coef(p%item(3), [-0.5_dp, 0.0_dp, 1.5_dp], 1.0e-13_dp, 'Legendre P2', failures)
  call check_coef(p%item(4), [0.0_dp, -1.5_dp, 0.0_dp, 2.5_dp], 1.0e-13_dp, 'Legendre P3', failures)

  p = slegendre_polynomials(2)
  call check_coef(p%item(2), [-1.0_dp, 2.0_dp], 1.0e-13_dp, 'shifted Legendre P1', failures)
  call check_coef(p%item(3), [1.0_dp, -6.0_dp, 6.0_dp], 1.0e-13_dp, 'shifted Legendre P2', failures)

  p = hermite_h_polynomials(3)
  call check_coef(p%item(3), [-2.0_dp, 0.0_dp, 4.0_dp], 1.0e-13_dp, 'Hermite H2', failures)
  call check_coef(p%item(4), [0.0_dp, -12.0_dp, 0.0_dp, 8.0_dp], 1.0e-13_dp, 'Hermite H3', failures)

  p = hermite_he_polynomials(3)
  call check_coef(p%item(3), [-1.0_dp, 0.0_dp, 1.0_dp], 1.0e-13_dp, 'Hermite He2', failures)
  call check_coef(p%item(4), [0.0_dp, -3.0_dp, 0.0_dp, 1.0_dp], 1.0e-13_dp, 'Hermite He3', failures)

  p = laguerre_polynomials(2)
  call check_coef(p%item(2), [1.0_dp, -1.0_dp], 1.0e-13_dp, 'Laguerre L1', failures)
  call check_coef(p%item(3), [1.0_dp, -2.0_dp, 0.5_dp], 1.0e-13_dp, 'Laguerre L2', failures)

  p = glaguerre_polynomials(2, 1.0_dp)
  call check_coef(p%item(2), [2.0_dp, -1.0_dp], 1.0e-13_dp, 'gen Laguerre L1', failures)
  call check_coef(p%item(3), [3.0_dp, -3.0_dp, 0.5_dp], 1.0e-13_dp, 'gen Laguerre L2', failures)

  p = gegenbauer_polynomials(3, 1.0_dp)
  call check_coef(p%item(4), [0.0_dp, -4.0_dp, 0.0_dp, 8.0_dp], 1.0e-13_dp, 'Gegenbauer alpha=1', failures)

  p = jacobi_p_polynomials(1, 1.0_dp, 2.0_dp)
  call check_coef(p%item(2), [-0.5_dp, 2.5_dp], 1.0e-12_dp, 'Jacobi P1(1,2)', failures)

  p = ghermite_h_polynomials(3, 0.0_dp)
  call check_coef(p%item(4), [0.0_dp, -12.0_dp, 0.0_dp, 8.0_dp], 1.0e-13_dp, 'gen Hermite mu=0', failures)

  p = spherical_polynomials(3)
  call check_coef(p%item(4), [0.0_dp, -1.5_dp, 0.0_dp, 2.5_dp], 1.0e-13_dp, 'spherical alias', failures)

  p = ultraspherical_polynomials(3, 1.0_dp)
  call check_coef(p%item(4), [0.0_dp, -4.0_dp, 0.0_dp, 8.0_dp], 1.0e-13_dp, 'ultraspherical alias', failures)

  h = legendre_inner_products(3)
  call check_vector(h, [2.0_dp, 2.0_dp/3.0_dp, 2.0_dp/5.0_dp, 2.0_dp/7.0_dp], &
    1.0e-13_dp, 'Legendre inner products', failures)
  h = chebyshev_t_inner_products(2)
  call check_vector(h, [acos(-1.0_dp), 0.5_dp*acos(-1.0_dp), 0.5_dp*acos(-1.0_dp)], &
    1.0e-13_dp, 'Chebyshev T inner products', failures)
  h = glaguerre_inner_products(2, 1.0_dp)
  call check_vector(h, [1.0_dp, 2.0_dp, 3.0_dp], 1.0e-12_dp, 'gen Laguerre norms', failures)

  call check_close(chebyshev_t_weight(0.5_dp), 2.0_dp/sqrt(3.0_dp), 1.0e-13_dp, 'T weight', failures)
  call check_close(schebyshev_t_weight(0.25_dp), 4.0_dp/sqrt(3.0_dp), 1.0e-13_dp, 'shifted T weight', failures)
  call check_close(jacobi_p_weight(0.25_dp, 1.0_dp, 2.0_dp), 0.75_dp*1.25_dp**2, &
    1.0e-13_dp, 'Jacobi P weight', failures)
  call check_close(jacobi_g_weight(0.25_dp, 2.0_dp, 1.0_dp), 0.75_dp, 1.0e-13_dp, 'Jacobi G weight', failures)
  call check_close(ghermite_h_weight(2.0_dp, 0.5_dp), 2.0_dp*exp(-4.0_dp), 1.0e-13_dp, &
    'generalized Hermite weight', failures)

  r = legendre_recurrences(2, normalized=.false.)
  call check_vector(r%c, [1.0_dp, 2.0_dp, 3.0_dp], 1.0e-13_dp, 'Legendre recurrence c', failures)
  call check_vector(r%e, [1.0_dp, 3.0_dp, 5.0_dp], 1.0e-13_dp, 'Legendre recurrence e', failures)
  call check_vector(r%f, [0.0_dp, 1.0_dp, 2.0_dp], 1.0e-13_dp, 'Legendre recurrence f', failures)

  p = legendre_polynomials(2, normalized=.true.)
  call check_close(p%item(1)%coef(1), 1.0_dp/sqrt(2.0_dp), 1.0e-13_dp, 'normalized P0', failures)
  call check_close(p%item(2)%coef(2), sqrt(3.0_dp/2.0_dp), 1.0e-13_dp, 'normalized P1', failures)

  if (failures /= 0) then
    print '(a,i0)', 'test_families: FAIL ', failures
    error stop 1
  end if
  print '(a)', 'test_families: PASS'

contains

  subroutine check_coef(poly, expected, tol, label, failures)
    type(polynomial_t), intent(in) :: poly
    real(dp), intent(in) :: expected(:), tol
    character(*), intent(in) :: label
    integer, intent(inout) :: failures
    real(dp), allocatable :: got(:)
    got = poly%coefficients()
    if (size(got) /= size(expected)) then
      failures = failures + 1
      print '(a,a)', 'size mismatch: ', label
    else if (maxval(abs(got - expected)) > tol) then
      failures = failures + 1
      print '(a,a,1x,es12.4)', 'coefficient mismatch: ', label, maxval(abs(got - expected))
    end if
  end subroutine check_coef

  subroutine check_vector(got, expected, tol, label, failures)
    real(dp), intent(in) :: got(:), expected(:), tol
    character(*), intent(in) :: label
    integer, intent(inout) :: failures
    if (size(got) /= size(expected) .or. maxval(abs(got - expected)) > tol) then
      failures = failures + 1
      print '(a,a)', 'vector mismatch: ', label
    end if
  end subroutine check_vector

  subroutine check_close(got, expected, tol, label, failures)
    real(dp), intent(in) :: got, expected, tol
    character(*), intent(in) :: label
    integer, intent(inout) :: failures
    if (abs(got - expected) > tol) then
      failures = failures + 1
      print '(a,a,2(1x,es16.8))', 'value mismatch: ', label, got, expected
    end if
  end subroutine check_close

end program test_families
