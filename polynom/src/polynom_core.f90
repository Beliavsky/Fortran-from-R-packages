module polynom_core
  use polynom_kinds, only : dp
  use polynom_status, only : poly_status_t, poly_ok, poly_divide_by_zero, set_status
  implicit none
  private

  real(dp), parameter :: default_tol = 0.0_dp

  type, public :: polynomial_t
    real(dp), allocatable :: coef(:)
  contains
    procedure :: degree => polynomial_degree
    procedure :: is_zero => polynomial_is_zero
    procedure :: evaluate_scalar
    procedure :: evaluate_vector
    generic :: evaluate => evaluate_scalar, evaluate_vector
    procedure :: coefficients => polynomial_coefficients
    procedure :: to_string => polynomial_to_string
  end type polynomial_t

  type, public :: polylist_t
    type(polynomial_t), allocatable :: item(:)
  contains
    procedure :: size => polylist_size
  end type polylist_t

  interface polynomial
    module procedure make_polynomial
    module procedure make_polynomial_scalar
    module procedure make_polynomial_integer
    module procedure make_polynomial_integer_scalar
    module procedure make_polynomial_default
  end interface polynomial

  interface operator(+)
    module procedure add_pp, add_pr, add_rp
  end interface
  interface operator(-)
    module procedure subtract_pp, subtract_pr, subtract_rp, negate_p
  end interface
  interface operator(*)
    module procedure multiply_pp, multiply_pr, multiply_rp
  end interface
  interface operator(/)
    module procedure divide_pp, divide_pr
  end interface
  interface operator(**)
    module procedure power_pi
  end interface
  interface operator(==)
    module procedure equal_pp
  end interface
  interface operator(/=)
    module procedure not_equal_pp
  end interface

  public :: polynomial
  public :: operator(+), operator(-), operator(*), operator(/), operator(**)
  public :: operator(==), operator(/=)
  public :: poly_divmod, poly_rem, sum_polynomials, product_polynomials
  public :: trim_polynomial

contains

  function make_polynomial(coef, tol) result(p)
    real(dp), intent(in) :: coef(:)
    real(dp), intent(in), optional :: tol
    type(polynomial_t) :: p
    real(dp) :: eps
    integer :: n

    eps = default_tol
    if (present(tol)) eps = max(0.0_dp, tol)
    n = size(coef)
    if (n == 0) then
      allocate(p%coef(1), source=0.0_dp)
      return
    end if
    do while (n > 1 .and. abs(coef(n)) <= eps)
      n = n - 1
    end do
    allocate(p%coef(n))
    p%coef = coef(:n)
    where (abs(p%coef) <= eps) p%coef = 0.0_dp
  end function make_polynomial

  function make_polynomial_scalar(coef) result(p)
    real(dp), intent(in) :: coef
    type(polynomial_t) :: p
    allocate(p%coef(1), source=coef)
  end function make_polynomial_scalar


  function make_polynomial_integer(coef) result(p)
    integer, intent(in) :: coef(:)
    type(polynomial_t) :: p
    p = make_polynomial(real(coef, dp))
  end function make_polynomial_integer

  function make_polynomial_integer_scalar(coef) result(p)
    integer, intent(in) :: coef
    type(polynomial_t) :: p
    p = make_polynomial_scalar(real(coef, dp))
  end function make_polynomial_integer_scalar

  function make_polynomial_default() result(p)
    type(polynomial_t) :: p
    p = make_polynomial([0.0_dp, 1.0_dp])
  end function make_polynomial_default

  function trim_polynomial(p, tol) result(q)
    type(polynomial_t), intent(in) :: p
    real(dp), intent(in), optional :: tol
    type(polynomial_t) :: q
    q = polynomial(p%coef, tol)
  end function trim_polynomial

  integer function polynomial_degree(self)
    class(polynomial_t), intent(in) :: self
    if (.not. allocated(self%coef)) then
      polynomial_degree = 0
    else
      polynomial_degree = size(self%coef) - 1
    end if
  end function polynomial_degree

  logical function polynomial_is_zero(self, tol)
    class(polynomial_t), intent(in) :: self
    real(dp), intent(in), optional :: tol
    real(dp) :: eps
    eps = default_tol
    if (present(tol)) eps = max(0.0_dp, tol)
    polynomial_is_zero = .not. allocated(self%coef) .or. all(abs(self%coef) <= eps)
  end function polynomial_is_zero

  function polynomial_coefficients(self) result(coef)
    class(polynomial_t), intent(in) :: self
    real(dp), allocatable :: coef(:)
    if (allocated(self%coef)) then
      coef = self%coef
    else
      allocate(coef(1), source=0.0_dp)
    end if
  end function polynomial_coefficients

  pure real(dp) function evaluate_scalar(self, x)
    class(polynomial_t), intent(in) :: self
    real(dp), intent(in) :: x
    integer :: i
    evaluate_scalar = 0.0_dp
    if (.not. allocated(self%coef)) return
    do i = size(self%coef), 1, -1
      evaluate_scalar = self%coef(i) + x * evaluate_scalar
    end do
  end function evaluate_scalar

  pure function evaluate_vector(self, x) result(y)
    class(polynomial_t), intent(in) :: self
    real(dp), intent(in) :: x(:)
    real(dp) :: y(size(x))
    integer :: i
    y = 0.0_dp
    if (.not. allocated(self%coef)) return
    do i = size(self%coef), 1, -1
      y = self%coef(i) + x * y
    end do
  end function evaluate_vector

  integer function polylist_size(self)
    class(polylist_t), intent(in) :: self
    if (allocated(self%item)) then
      polylist_size = size(self%item)
    else
      polylist_size = 0
    end if
  end function polylist_size

  function add_pp(a, b) result(c)
    type(polynomial_t), intent(in) :: a, b
    type(polynomial_t) :: c
    real(dp), allocatable :: x(:)
    integer :: n
    n = max(size(a%coef), size(b%coef))
    allocate(x(n), source=0.0_dp)
    x(:size(a%coef)) = x(:size(a%coef)) + a%coef
    x(:size(b%coef)) = x(:size(b%coef)) + b%coef
    c = polynomial(x)
  end function add_pp

  function add_pr(a, b) result(c)
    type(polynomial_t), intent(in) :: a
    real(dp), intent(in) :: b
    type(polynomial_t) :: c
    c = a
    c%coef(1) = c%coef(1) + b
    c = trim_polynomial(c)
  end function add_pr

  function add_rp(a, b) result(c)
    real(dp), intent(in) :: a
    type(polynomial_t), intent(in) :: b
    type(polynomial_t) :: c
    c = b + a
  end function add_rp

  function subtract_pp(a, b) result(c)
    type(polynomial_t), intent(in) :: a, b
    type(polynomial_t) :: c
    real(dp), allocatable :: x(:)
    integer :: n
    n = max(size(a%coef), size(b%coef))
    allocate(x(n), source=0.0_dp)
    x(:size(a%coef)) = x(:size(a%coef)) + a%coef
    x(:size(b%coef)) = x(:size(b%coef)) - b%coef
    c = polynomial(x)
  end function subtract_pp

  function subtract_pr(a, b) result(c)
    type(polynomial_t), intent(in) :: a
    real(dp), intent(in) :: b
    type(polynomial_t) :: c
    c = a + (-b)
  end function subtract_pr

  function subtract_rp(a, b) result(c)
    real(dp), intent(in) :: a
    type(polynomial_t), intent(in) :: b
    type(polynomial_t) :: c
    c = polynomial(a) - b
  end function subtract_rp

  function negate_p(a) result(c)
    type(polynomial_t), intent(in) :: a
    type(polynomial_t) :: c
    c = polynomial(-a%coef)
  end function negate_p

  function multiply_pp(a, b) result(c)
    type(polynomial_t), intent(in) :: a, b
    type(polynomial_t) :: c
    real(dp), allocatable :: x(:)
    integer :: i, j
    allocate(x(size(a%coef) + size(b%coef) - 1), source=0.0_dp)
    do j = 1, size(b%coef)
      do i = 1, size(a%coef)
        x(i + j - 1) = x(i + j - 1) + a%coef(i) * b%coef(j)
      end do
    end do
    c = polynomial(x)
  end function multiply_pp

  function multiply_pr(a, b) result(c)
    type(polynomial_t), intent(in) :: a
    real(dp), intent(in) :: b
    type(polynomial_t) :: c
    c = polynomial(a%coef * b)
  end function multiply_pr

  function multiply_rp(a, b) result(c)
    real(dp), intent(in) :: a
    type(polynomial_t), intent(in) :: b
    type(polynomial_t) :: c
    c = b * a
  end function multiply_rp

  subroutine poly_divmod(numerator, denominator, quotient, remainder, status, tol)
    type(polynomial_t), intent(in) :: numerator, denominator
    type(polynomial_t), intent(out) :: quotient, remainder
    type(poly_status_t), intent(out), optional :: status
    real(dp), intent(in), optional :: tol
    real(dp), allocatable :: work(:), q(:)
    real(dp) :: eps, factor, scale
    integer :: dn, dd, k, j

    eps = 100.0_dp * epsilon(1.0_dp)
    if (present(tol)) eps = max(0.0_dp, tol)
    if (denominator%is_zero(eps)) then
      quotient = polynomial(0.0_dp)
      remainder = numerator
      call set_status(status, poly_divide_by_zero, 'polynomial division by zero')
      return
    end if
    dn = numerator%degree()
    dd = denominator%degree()
    if (dn < dd) then
      quotient = polynomial(0.0_dp)
      remainder = numerator
      call set_status(status, poly_ok, '')
      return
    end if
    work = numerator%coef
    allocate(q(dn - dd + 1), source=0.0_dp)
    scale = max(1.0_dp, maxval(abs(work)), maxval(abs(denominator%coef)))
    do k = dn - dd, 0, -1
      factor = work(dd + k + 1) / denominator%coef(dd + 1)
      q(k + 1) = factor
      do j = 0, dd
        work(j + k + 1) = work(j + k + 1) - factor * denominator%coef(j + 1)
      end do
    end do
    where (abs(work) <= eps * scale) work = 0.0_dp
    quotient = polynomial(q, eps * scale)
    if (dd == 0) then
      remainder = polynomial(0.0_dp)
    else
      remainder = polynomial(work(:dd), eps * scale)
    end if
    call set_status(status, poly_ok, '')
  end subroutine poly_divmod

  function divide_pp(a, b) result(c)
    type(polynomial_t), intent(in) :: a, b
    type(polynomial_t) :: c, r
    call poly_divmod(a, b, c, r)
  end function divide_pp

  function divide_pr(a, b) result(c)
    type(polynomial_t), intent(in) :: a
    real(dp), intent(in) :: b
    type(polynomial_t) :: c
    if (abs(b) <= tiny(1.0_dp)) then
      c = polynomial(0.0_dp)
    else
      c = polynomial(a%coef / b)
    end if
  end function divide_pr

  function poly_rem(a, b, status, tol) result(r)
    type(polynomial_t), intent(in) :: a, b
    type(poly_status_t), intent(out), optional :: status
    real(dp), intent(in), optional :: tol
    type(polynomial_t) :: r, q
    call poly_divmod(a, b, q, r, status, tol)
  end function poly_rem

  function power_pi(a, exponent) result(c)
    type(polynomial_t), intent(in) :: a
    integer, intent(in) :: exponent
    type(polynomial_t) :: c, base
    integer :: n
    if (exponent < 0) then
      c = polynomial(0.0_dp)
      return
    end if
    c = polynomial(1.0_dp)
    base = a
    n = exponent
    do while (n > 0)
      if (mod(n, 2) == 1) c = c * base
      n = n / 2
      if (n > 0) base = base * base
    end do
  end function power_pi

  logical function equal_pp(a, b)
    type(polynomial_t), intent(in) :: a, b
    equal_pp = size(a%coef) == size(b%coef)
    if (equal_pp) equal_pp = all(.not. (a%coef < b%coef .or. a%coef > b%coef))
  end function equal_pp

  logical function not_equal_pp(a, b)
    type(polynomial_t), intent(in) :: a, b
    not_equal_pp = .not. (a == b)
  end function not_equal_pp

  function sum_polynomials(list) result(total)
    type(polylist_t), intent(in) :: list
    type(polynomial_t) :: total
    integer :: i
    total = polynomial(0.0_dp)
    if (.not. allocated(list%item)) return
    do i = 1, size(list%item)
      total = total + list%item(i)
    end do
  end function sum_polynomials

  function product_polynomials(list) result(total)
    type(polylist_t), intent(in) :: list
    type(polynomial_t) :: total
    integer :: i
    total = polynomial(1.0_dp)
    if (.not. allocated(list%item)) return
    do i = 1, size(list%item)
      total = total * list%item(i)
    end do
  end function product_polynomials

  function polynomial_to_string(self, decreasing, digits) result(text)
    class(polynomial_t), intent(in) :: self
    logical, intent(in), optional :: decreasing
    integer, intent(in), optional :: digits
    character(len=:), allocatable :: text
    character(len=128) :: term, num, fmt
    logical :: dec
    integer :: i, first, last, step, d, power
    real(dp) :: value

    dec = .false.
    if (present(decreasing)) dec = decreasing
    d = 7
    if (present(digits)) d = max(1, digits)
    write(fmt, '(a,i0,a)') '(g0.', d, ')'
    text = ''
    if (all(abs(self%coef) <= default_tol)) then
      text = '0'
      return
    end if
    if (dec) then
      first = size(self%coef)
      last = 1
      step = -1
    else
      first = 1
      last = size(self%coef)
      step = 1
    end if
    i = first
    do
      value = self%coef(i)
      power = i - 1
      if (abs(value) > default_tol) then
        write(num, fmt) abs(value)
        num = adjustl(num)
        term = ''
        if (len(text) == 0) then
          if (value < 0.0_dp) term = '-'
        else
          if (value < 0.0_dp) then
            term = '- '
          else
            term = '+ '
          end if
        end if
        if (power == 0) then
          term = trim(term) // trim(num)
        else if (power == 1) then
          if (abs(abs(value) - 1.0_dp) <= default_tol) then
            term = trim(term) // 'x'
          else
            term = trim(term) // trim(num) // '*x'
          end if
        else
          if (abs(abs(value) - 1.0_dp) <= default_tol) then
            write(term(len_trim(term)+1:), '(a,i0)') 'x^', power
          else
            term = trim(term) // trim(num)
            write(term(len_trim(term)+1:), '(a,i0)') '*x^', power
          end if
        end if
        if (len(text) == 0) then
          text = trim(term)
        else
          text = text // ' ' // trim(term)
        end if
      end if
      if (i == last) exit
      i = i + step
    end do
  end function polynomial_to_string

end module polynom_core
