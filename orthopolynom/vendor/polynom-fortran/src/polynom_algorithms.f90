module polynom_algorithms
  use polynom_kinds, only : dp
  use polynom_status, only : poly_status_t, poly_ok, poly_invalid_argument, &
    poly_duplicate_abscissa, set_status
  use polynom_core, only : polynomial_t, polylist_t, polynomial, operator(+), &
    operator(-), operator(*), operator(/), operator(**), poly_rem, poly_divmod
  implicit none
  private

  public :: derivative, derivative_polylist, integral_polynomial, integral_polylist, definite_integral
  public :: monic, change_origin
  public :: poly_from_roots, poly_from_values, poly_from_values_matrix
  public :: polynomial_gcd, polynomial_lcm, polylist_gcd, polylist_lcm
  public :: round_coefficients, floor_coefficients, ceiling_coefficients
  public :: truncate_coefficients, significant_coefficients

contains

  function derivative(p, order) result(d)
    type(polynomial_t), intent(in) :: p
    integer, intent(in), optional :: order
    type(polynomial_t) :: d
    real(dp), allocatable :: c(:)
    integer :: k, n, i

    n = 1
    if (present(order)) n = max(0, order)
    d = p
    do k = 1, n
      if (d%degree() == 0) then
        d = polynomial(0.0_dp)
        exit
      end if
      allocate(c(d%degree()))
      do i = 1, size(c)
        c(i) = real(i, dp) * d%coef(i + 1)
      end do
      d = polynomial(c)
      deallocate(c)
    end do
  end function derivative

  function derivative_polylist(list, order) result(result)
    type(polylist_t), intent(in) :: list
    integer, intent(in), optional :: order
    type(polylist_t) :: result
    integer :: i
    if (.not. allocated(list%item)) then
      allocate(result%item(0))
      return
    end if
    allocate(result%item(size(list%item)))
    do i = 1, size(list%item)
      result%item(i) = derivative(list%item(i), order)
    end do
  end function derivative_polylist

  function integral_polynomial(p, constant) result(q)
    type(polynomial_t), intent(in) :: p
    real(dp), intent(in), optional :: constant
    type(polynomial_t) :: q
    real(dp), allocatable :: c(:)
    integer :: i

    allocate(c(size(p%coef) + 1), source=0.0_dp)
    if (present(constant)) c(1) = constant
    do i = 1, size(p%coef)
      c(i + 1) = p%coef(i) / real(i, dp)
    end do
    q = polynomial(c)
  end function integral_polynomial

  function integral_polylist(list, constant) result(result)
    type(polylist_t), intent(in) :: list
    real(dp), intent(in), optional :: constant
    type(polylist_t) :: result
    integer :: i
    if (.not. allocated(list%item)) then
      allocate(result%item(0))
      return
    end if
    allocate(result%item(size(list%item)))
    do i = 1, size(list%item)
      result%item(i) = integral_polynomial(list%item(i), constant)
    end do
  end function integral_polylist

  real(dp) function definite_integral(p, lower, upper)
    type(polynomial_t), intent(in) :: p
    real(dp), intent(in) :: lower, upper
    type(polynomial_t) :: q
    q = integral_polynomial(p)
    definite_integral = q%evaluate(upper) - q%evaluate(lower)
  end function definite_integral

  function monic(p, status, tol) result(q)
    type(polynomial_t), intent(in) :: p
    type(poly_status_t), intent(out), optional :: status
    real(dp), intent(in), optional :: tol
    type(polynomial_t) :: q
    real(dp) :: eps

    eps = 100.0_dp * epsilon(1.0_dp)
    if (present(tol)) eps = max(0.0_dp, tol)
    if (p%is_zero(eps)) then
      q = polynomial(0.0_dp)
      call set_status(status, poly_invalid_argument, 'the zero polynomial has no monic form')
    else
      q = p / p%coef(size(p%coef))
      call set_status(status, poly_ok, '')
    end if
  end function monic

  function change_origin(p, origin) result(q)
    type(polynomial_t), intent(in) :: p
    real(dp), intent(in) :: origin
    type(polynomial_t) :: q
    real(dp), allocatable :: c(:)
    integer :: n, j, k

    n = p%degree()
    allocate(c(n + 1), source=0.0_dp)
    do k = 0, n
      do j = 0, k
        c(j + 1) = c(j + 1) + p%coef(k + 1) * binomial_real(k, j) * origin**(k - j)
      end do
    end do
    q = polynomial(c)
  end function change_origin

  pure real(dp) function binomial_real(n, k)
    integer, intent(in) :: n, k
    integer :: i, kk
    if (k < 0 .or. k > n) then
      binomial_real = 0.0_dp
      return
    end if
    kk = min(k, n - k)
    binomial_real = 1.0_dp
    do i = 1, kk
      binomial_real = binomial_real * real(n - kk + i, dp) / real(i, dp)
    end do
  end function binomial_real

  function poly_from_roots(roots) result(p)
    real(dp), intent(in) :: roots(:)
    type(polynomial_t) :: p
    integer :: i
    p = polynomial(1.0_dp)
    do i = 1, size(roots)
      p = p * polynomial([-roots(i), 1.0_dp])
    end do
  end function poly_from_roots

  function poly_from_values(x, y, tol, status) result(p)
    real(dp), intent(in) :: x(:), y(:)
    real(dp), intent(in), optional :: tol
    type(poly_status_t), intent(out), optional :: status
    type(polynomial_t) :: p, basis
    real(dp), allocatable :: xu(:), yu(:)
    logical, allocatable :: used(:)
    real(dp) :: eps, denom, spread
    integer :: i, j, m, count

    eps = sqrt(epsilon(1.0_dp))
    if (present(tol)) eps = max(0.0_dp, tol)
    if (size(x) /= size(y) .or. size(x) == 0) then
      p = polynomial(0.0_dp)
      call set_status(status, poly_invalid_argument, 'x and y must have equal positive size')
      return
    end if

    allocate(used(size(x)), source=.false.)
    allocate(xu(size(x)), yu(size(y)))
    count = 0
    do i = 1, size(x)
      if (used(i)) cycle
      count = count + 1
      xu(count) = x(i)
      yu(count) = y(i)
      used(i) = .true.
      spread = 0.0_dp
      do j = i + 1, size(x)
        if (.not. (x(j) < x(i) .or. x(j) > x(i))) then
          used(j) = .true.
          spread = max(spread, abs(y(j) - y(i)))
        end if
      end do
      if (spread > eps) then
        call set_status(status, poly_duplicate_abscissa, &
          'duplicated x values have inconsistent y values; first value used')
      end if
    end do
    m = count
    if (m == 1) then
      p = polynomial(yu(1))
      if (.not. present(status)) return
      if (status%code == 0) call set_status(status, poly_ok, '')
      return
    end if

    p = polynomial(0.0_dp)
    do i = 1, m
      basis = polynomial(1.0_dp)
      denom = 1.0_dp
      do j = 1, m
        if (j == i) cycle
        basis = basis * polynomial([-xu(j), 1.0_dp])
        denom = denom * (xu(i) - xu(j))
      end do
      p = p + basis * (yu(i) / denom)
    end do
    where (abs(p%coef) < eps) p%coef = 0.0_dp
    p = polynomial(p%coef, eps)
    if (present(status)) then
      if (status%code == 0) call set_status(status, poly_ok, '')
    end if
  end function poly_from_values

  function poly_from_values_matrix(x, y, tol, status) result(list)
    real(dp), intent(in) :: x(:), y(:,:)
    real(dp), intent(in), optional :: tol
    type(poly_status_t), intent(out), optional :: status
    type(polylist_t) :: list
    type(poly_status_t) :: local_status
    integer :: j

    if (size(y, 1) /= size(x)) then
      allocate(list%item(0))
      call set_status(status, poly_invalid_argument, 'rows of y must equal size of x')
      return
    end if
    allocate(list%item(size(y, 2)))
    do j = 1, size(y, 2)
      list%item(j) = poly_from_values(x, y(:, j), tol, local_status)
      if (.not. local_status%succeeded()) then
        call set_status(status, local_status%code, local_status%message)
      end if
    end do
    if (present(status)) then
      if (status%code == 0) call set_status(status, poly_ok, '')
    end if
  end function poly_from_values_matrix

  recursive function polynomial_gcd(a, b, tol) result(g)
    type(polynomial_t), intent(in) :: a, b
    real(dp), intent(in), optional :: tol
    type(polynomial_t) :: g, r
    real(dp) :: eps

    eps = sqrt(epsilon(1.0_dp))
    if (present(tol)) eps = max(0.0_dp, tol)
    if (b%is_zero(eps)) then
      g = a
    else if (b%degree() == 0) then
      g = polynomial(1.0_dp)
    else
      r = poly_rem(a, b, tol=eps)
      g = polynomial_gcd(b, r, eps)
    end if
  end function polynomial_gcd

  function polynomial_lcm(a, b, tol) result(l)
    type(polynomial_t), intent(in) :: a, b
    real(dp), intent(in), optional :: tol
    type(polynomial_t) :: l, g
    real(dp) :: eps

    eps = sqrt(epsilon(1.0_dp))
    if (present(tol)) eps = max(0.0_dp, tol)
    if (a%is_zero(eps) .or. b%is_zero(eps)) then
      l = polynomial(0.0_dp)
    else
      g = polynomial_gcd(a, b, eps)
      l = (a / g) * b
    end if
  end function polynomial_lcm

  function polylist_gcd(list, tol, status) result(g)
    type(polylist_t), intent(in) :: list
    real(dp), intent(in), optional :: tol
    type(poly_status_t), intent(out), optional :: status
    type(polynomial_t) :: g
    integer :: i

    if (.not. allocated(list%item) .or. size(list%item) < 2) then
      g = polynomial(0.0_dp)
      call set_status(status, poly_invalid_argument, 'at least two polynomials are required')
      return
    end if
    g = list%item(1)
    do i = 2, size(list%item)
      g = polynomial_gcd(g, list%item(i), tol)
    end do
    call set_status(status, poly_ok, '')
  end function polylist_gcd

  function polylist_lcm(list, tol, status) result(l)
    type(polylist_t), intent(in) :: list
    real(dp), intent(in), optional :: tol
    type(poly_status_t), intent(out), optional :: status
    type(polynomial_t) :: l
    integer :: i

    if (.not. allocated(list%item) .or. size(list%item) < 2) then
      l = polynomial(0.0_dp)
      call set_status(status, poly_invalid_argument, 'at least two polynomials are required')
      return
    end if
    l = list%item(1)
    do i = 2, size(list%item)
      l = polynomial_lcm(l, list%item(i), tol)
    end do
    call set_status(status, poly_ok, '')
  end function polylist_lcm

  function round_coefficients(p, digits) result(q)
    type(polynomial_t), intent(in) :: p
    integer, intent(in), optional :: digits
    type(polynomial_t) :: q
    real(dp) :: factor
    integer :: d
    d = 0
    if (present(digits)) d = digits
    factor = 10.0_dp**d
    q = polynomial(anint(p%coef * factor) / factor)
  end function round_coefficients

  function floor_coefficients(p) result(q)
    type(polynomial_t), intent(in) :: p
    type(polynomial_t) :: q
    q = polynomial(real(floor(p%coef), dp))
  end function floor_coefficients

  function ceiling_coefficients(p) result(q)
    type(polynomial_t), intent(in) :: p
    type(polynomial_t) :: q
    q = polynomial(ceiling(p%coef, kind=kind(1)) * 1.0_dp)
  end function ceiling_coefficients

  function truncate_coefficients(p) result(q)
    type(polynomial_t), intent(in) :: p
    type(polynomial_t) :: q
    q = polynomial(aint(p%coef))
  end function truncate_coefficients

  function significant_coefficients(p, digits) result(q)
    type(polynomial_t), intent(in) :: p
    integer, intent(in) :: digits
    type(polynomial_t) :: q
    real(dp), allocatable :: c(:)
    real(dp) :: scale
    integer :: i
    allocate(c(size(p%coef)))
    do i = 1, size(c)
      if (abs(p%coef(i)) <= tiny(1.0_dp)) then
        c(i) = 0.0_dp
      else
        scale = 10.0_dp**(digits - 1 - floor(log10(abs(p%coef(i)))))
        c(i) = anint(p%coef(i) * scale) / scale
      end if
    end do
    q = polynomial(c)
  end function significant_coefficients

end module polynom_algorithms
