module polynom_compat
  use polynom_kinds, only : dp
  use polynom_status, only : poly_status_t
  use polynom_core, only : polynomial_t, polylist_t
  use polynom_algorithms, only : derivative, integral_polynomial, definite_integral, &
    poly_from_roots, poly_from_values, polynomial_gcd, polynomial_lcm, &
    polylist_gcd, polylist_lcm
  use polynom_roots, only : polynomial_roots
  use polynom_orthogonal, only : orthogonal_polynomials
  implicit none
  private

  interface predict_polynomial
    module procedure predict_polynomial_scalar
    module procedure predict_polynomial_vector
  end interface predict_polynomial

  interface integral
    module procedure integral_indefinite
    module procedure integral_definite_limits
  end interface integral

  interface gcd
    module procedure gcd_pair
    module procedure gcd_list
  end interface gcd

  interface lcm
    module procedure lcm_pair
    module procedure lcm_list
  end interface lcm

  public :: poly_calc, poly_from_zeros, solve_polynomial, poly_orth
  public :: deriv_polynomial, predict_polynomial, coef_polynomial
  public :: integral, gcd, lcm

contains

  function poly_calc(x, y, tol, status) result(p)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: y(:)
    real(dp), intent(in), optional :: tol
    type(poly_status_t), intent(out), optional :: status
    type(polynomial_t) :: p
    if (present(y)) then
      p = poly_from_values(x, y, tol, status)
    else
      p = poly_from_roots(x)
    end if
  end function poly_calc

  function poly_from_zeros(zeros) result(p)
    real(dp), intent(in) :: zeros(:)
    type(polynomial_t) :: p
    p = poly_from_roots(zeros)
  end function poly_from_zeros

  function solve_polynomial(p, intercept, status, tol, max_iter) result(roots)
    type(polynomial_t), intent(in) :: p
    real(dp), intent(in), optional :: intercept
    type(poly_status_t), intent(out), optional :: status
    real(dp), intent(in), optional :: tol
    integer, intent(in), optional :: max_iter
    complex(dp), allocatable :: roots(:)
    roots = polynomial_roots(p, intercept, status, tol, max_iter)
  end function solve_polynomial

  function poly_orth(x, degree, norm, status) result(list)
    real(dp), intent(in) :: x(:)
    integer, intent(in), optional :: degree
    logical, intent(in), optional :: norm
    type(poly_status_t), intent(out), optional :: status
    type(polylist_t) :: list
    list = orthogonal_polynomials(x, degree, norm, status)
  end function poly_orth

  function deriv_polynomial(p, order) result(d)
    type(polynomial_t), intent(in) :: p
    integer, intent(in), optional :: order
    type(polynomial_t) :: d
    d = derivative(p, order)
  end function deriv_polynomial

  pure real(dp) function predict_polynomial_scalar(p, x)
    type(polynomial_t), intent(in) :: p
    real(dp), intent(in) :: x
    predict_polynomial_scalar = p%evaluate(x)
  end function predict_polynomial_scalar

  pure function predict_polynomial_vector(p, x) result(y)
    type(polynomial_t), intent(in) :: p
    real(dp), intent(in) :: x(:)
    real(dp) :: y(size(x))
    y = p%evaluate(x)
  end function predict_polynomial_vector

  function coef_polynomial(p) result(coef)
    type(polynomial_t), intent(in) :: p
    real(dp), allocatable :: coef(:)
    coef = p%coefficients()
  end function coef_polynomial

  function integral_indefinite(p, constant) result(q)
    type(polynomial_t), intent(in) :: p
    real(dp), intent(in), optional :: constant
    type(polynomial_t) :: q
    q = integral_polynomial(p, constant)
  end function integral_indefinite

  real(dp) function integral_definite_limits(p, limits)
    type(polynomial_t), intent(in) :: p
    real(dp), intent(in) :: limits(2)
    integral_definite_limits = definite_integral(p, limits(1), limits(2))
  end function integral_definite_limits

  function gcd_pair(a, b, tol) result(g)
    type(polynomial_t), intent(in) :: a, b
    real(dp), intent(in), optional :: tol
    type(polynomial_t) :: g
    g = polynomial_gcd(a, b, tol)
  end function gcd_pair

  function gcd_list(list, tol, status) result(g)
    type(polylist_t), intent(in) :: list
    real(dp), intent(in), optional :: tol
    type(poly_status_t), intent(out), optional :: status
    type(polynomial_t) :: g
    g = polylist_gcd(list, tol, status)
  end function gcd_list

  function lcm_pair(a, b, tol) result(value)
    type(polynomial_t), intent(in) :: a, b
    real(dp), intent(in), optional :: tol
    type(polynomial_t) :: value
    value = polynomial_lcm(a, b, tol)
  end function lcm_pair

  function lcm_list(list, tol, status) result(value)
    type(polylist_t), intent(in) :: list
    real(dp), intent(in), optional :: tol
    type(poly_status_t), intent(out), optional :: status
    type(polynomial_t) :: value
    value = polylist_lcm(list, tol, status)
  end function lcm_list

end module polynom_compat
