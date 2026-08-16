program test_helpers
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_negative_inf, ieee_positive_inf
  use orthopolynom
  implicit none
  integer :: failures
  type(recurrence_t) :: r
  type(monic_recurrence_t) :: mr
  type(polylist_t) :: p, dpdx, ip
  type(real_matrix_list_t) :: matrices
  type(real_vector_list_t) :: vectors, roots, powers
  type(polynomial_function_list_t) :: funcs
  integer, allocatable :: orders(:)
  real(dp) :: x(3), ninf, pinf

  failures = 0
  ninf = ieee_value(0.0_dp, ieee_negative_inf)
  pinf = ieee_value(0.0_dp, ieee_positive_inf)

  call check_close(pochhammer(2.5_dp, 3), 2.5_dp*3.5_dp*4.5_dp, 1.0e-13_dp, 'pochhammer', failures)
  call check_close(lpochhammer(2.5_dp, 3), log(2.5_dp*3.5_dp*4.5_dp), 1.0e-13_dp, 'lpochhammer', failures)
  call check_close(lpochhammer(2.5_dp, 0), 1.0_dp, 0.0_dp, 'lpochhammer upstream n=0', failures)

  r = legendre_recurrences(3)
  mr = monic_polynomial_recurrences(r)
  call check_vector(mr%a, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], 1.0e-13_dp, 'monic a', failures)
  call check_vector(mr%b, [0.0_dp, 1.0_dp/3.0_dp, 4.0_dp/15.0_dp, 9.0_dp/35.0_dp], &
    1.0e-13_dp, 'monic b', failures)

  p = monic_polynomials(mr)
  call check_coef(p%item(3), [-1.0_dp/3.0_dp, 0.0_dp, 1.0_dp], 1.0e-13_dp, 'monic P2', failures)

  matrices = jacobi_matrices(mr)
  call check_int(matrices%size(), 3, 'jacobi matrix count', failures)
  call check_close(matrices%item(2)%value(1,2), sqrt(1.0_dp/3.0_dp), 1.0e-13_dp, &
    'jacobi offdiag', failures)

  roots = polynomial_roots(monic_polynomial_recurrences(chebyshev_t_recurrences(3)))
  call check_int(roots%size(), 4, 'root list size', failures)
  call check_vector(roots%item(4)%value, [-sqrt(3.0_dp)/2.0_dp, 0.0_dp, sqrt(3.0_dp)/2.0_dp], &
    1.0e-9_dp, 'Chebyshev T3 roots', failures)

  p = orthogonal_polynomials(legendre_recurrences(2))
  vectors = polynomial_coefficients(p)
  call check_vector(vectors%item(3)%value, [-0.5_dp, 0.0_dp, 1.5_dp], 1.0e-13_dp, &
    'polynomial coefficients', failures)

  dpdx = polynomial_derivatives(p)
  call check_coef(dpdx%item(3), [0.0_dp, 3.0_dp], 1.0e-13_dp, 'polynomial derivative', failures)

  ip = polynomial_integrals(p)
  call check_coef(ip%item(2), [0.0_dp, 0.0_dp, 0.5_dp], 1.0e-13_dp, 'polynomial integral', failures)

  orders = polynomial_orders(p)
  call check_int_vector(orders, [0, 1, 2], 'polynomial orders', failures)

  powers = polynomial_powers(p)
  call check_vector(powers%item(1)%value, [1.0_dp], 1.0e-13_dp, 'power 0', failures)
  call check_vector(powers%item(2)%value, [0.0_dp, 1.0_dp], 1.0e-13_dp, 'power 1', failures)
  call check_vector(powers%item(3)%value, [1.0_dp/3.0_dp, 0.0_dp, 2.0_dp/3.0_dp], &
    1.0e-13_dp, 'power 2', failures)

  x = [-0.5_dp, 0.0_dp, 0.5_dp]
  vectors = polynomial_values(p, x)
  call check_vector(vectors%item(3)%value, [-0.125_dp, -0.5_dp, -0.125_dp], 1.0e-13_dp, &
    'polynomial values', failures)

  funcs = polynomial_functions(p)
  call check_close(funcs%item(3)%evaluate(0.5_dp), -0.125_dp, 1.0e-13_dp, 'function wrapper', failures)

  call check_vector(scale_x([0.0_dp, 5.0_dp, 10.0_dp], -1.0_dp, 1.0_dp), [-1.0_dp, 0.0_dp, 1.0_dp], &
    1.0e-13_dp, 'scale finite', failures)
  call check_vector(scale_x([0.0_dp, 5.0_dp, 10.0_dp], ninf, pinf), [0.0_dp, 5.0_dp, 10.0_dp], &
    1.0e-13_dp, 'scale both infinite', failures)
  call check_vector(scale_x([0.0_dp, 5.0_dp, 10.0_dp], ninf, 2.0_dp), [-8.0_dp, -3.0_dp, 2.0_dp], &
    1.0e-13_dp, 'scale upper finite', failures)
  call check_vector(scale_x([0.0_dp, 5.0_dp, 10.0_dp], -2.0_dp, pinf), [-2.0_dp, 3.0_dp, 8.0_dp], &
    1.0e-13_dp, 'scale lower finite', failures)

  if (failures /= 0) then
    print '(a,i0)', 'test_helpers: FAIL ', failures
    error stop 1
  end if
  print '(a)', 'test_helpers: PASS'

contains

  subroutine check_coef(poly, expected, tol, label, failures)
    type(polynomial_t), intent(in) :: poly
    real(dp), intent(in) :: expected(:), tol
    character(*), intent(in) :: label
    integer, intent(inout) :: failures
    real(dp), allocatable :: got(:)
    got = poly%coefficients()
    if (size(got) /= size(expected) .or. maxval(abs(got - expected)) > tol) then
      failures = failures + 1
      print '(a,a)', 'coefficient mismatch: ', label
    end if
  end subroutine check_coef

  subroutine check_vector(got, expected, tol, label, failures)
    real(dp), intent(in) :: got(:), expected(:), tol
    character(*), intent(in) :: label
    integer, intent(inout) :: failures
    if (size(got) /= size(expected)) then
      failures = failures + 1
      print '(a,a)', 'vector size mismatch: ', label
    else if (size(got) > 0 .and. maxval(abs(got - expected)) > tol) then
      failures = failures + 1
      print '(a,a,1x,es12.4)', 'vector mismatch: ', label, maxval(abs(got - expected))
    end if
  end subroutine check_vector

  subroutine check_int_vector(got, expected, label, failures)
    integer, intent(in) :: got(:), expected(:)
    character(*), intent(in) :: label
    integer, intent(inout) :: failures
    if (size(got) /= size(expected) .or. any(got /= expected)) then
      failures = failures + 1
      print '(a,a)', 'integer vector mismatch: ', label
    end if
  end subroutine check_int_vector

  subroutine check_int(got, expected, label, failures)
    integer, intent(in) :: got, expected
    character(*), intent(in) :: label
    integer, intent(inout) :: failures
    if (got /= expected) then
      failures = failures + 1
      print '(a,a,2(1x,i0))', 'integer mismatch: ', label, got, expected
    end if
  end subroutine check_int

  subroutine check_close(got, expected, tol, label, failures)
    real(dp), intent(in) :: got, expected, tol
    character(*), intent(in) :: label
    integer, intent(inout) :: failures
    if (abs(got - expected) > tol) then
      failures = failures + 1
      print '(a,a,2(1x,es16.8))', 'value mismatch: ', label, got, expected
    end if
  end subroutine check_close

end program test_helpers
