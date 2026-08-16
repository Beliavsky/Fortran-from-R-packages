module orthopolynom_core
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use polynom, only : dp, polynomial_t, polylist_t, polynomial, derivative, integral_polynomial, &
    polynom_roots => polynomial_roots, operator(+), operator(-), operator(*), operator(/)
  use orthopolynom_types
  implicit none
  private

  public :: pochhammer, lpochhammer
  public :: orthogonal_polynomials_from_recurrence, orthonormal_polynomials_from_recurrence
  public :: orthogonal_polynomials, orthonormal_polynomials
  public :: monic_polynomial_recurrences, monic_polynomials
  public :: jacobi_matrices
  public :: polynomial_coefficients, polynomial_derivatives, polynomial_integrals
  public :: polynomial_orders, polynomial_powers, polynomial_roots_from_recurrence, polynomial_roots
  public :: polynomial_values, polynomial_functions, scale_x

interface orthogonal_polynomials
    module procedure orthogonal_polynomials_from_recurrence
  end interface orthogonal_polynomials

  interface orthonormal_polynomials
    module procedure orthonormal_polynomials_from_recurrence
  end interface orthonormal_polynomials

  interface polynomial_roots
    module procedure polynomial_roots_from_recurrence
  end interface polynomial_roots

contains

  pure real(dp) function pochhammer(z, n) result(value)
    real(dp), intent(in) :: z
    integer, intent(in) :: n
    integer :: i

    if (n < 0) error stop 'pochhammer: n is negative'
    value = 1.0_dp
    do i = 1, n
      value = value * (z + real(i - 1, dp))
    end do
  end function pochhammer

  pure real(dp) function lpochhammer(z, n) result(value)
    real(dp), intent(in) :: z
    integer, intent(in) :: n

    if (n < 0) error stop 'lpochhammer: n is negative'
    if (n == 0) then
      ! Upstream orthopolynom returns 1 here, although log((z)_0) is 0.
      value = 1.0_dp
    else
      value = log_gamma(z + real(n, dp)) - log_gamma(z)
    end if
  end function lpochhammer

  function orthogonal_polynomials_from_recurrence(recurrence) result(polynomials)
    type(recurrence_t), intent(in) :: recurrence
    type(polylist_t) :: polynomials
    integer :: n, j
    type(polynomial_t) :: p0, pjm1, pj, pjp1, monomial
    real(dp) :: cj, dj, ej, fj

    call validate_recurrence(recurrence)
    n = recurrence%size() - 1
    allocate(polynomials%item(n + 1))
    p0 = polynomial(1.0_dp)
    polynomials%item(1) = p0
    do j = 0, n - 1
      cj = recurrence%c(j + 1)
      dj = recurrence%d(j + 1)
      ej = recurrence%e(j + 1)
      fj = recurrence%f(j + 1)
      if (abs(cj) <= tiny(1.0_dp)) error stop 'orthogonal_polynomials: zero recurrence c coefficient'
      monomial = polynomial([dj, ej])
      if (j == 0) then
        pjp1 = (monomial * p0) / cj
      else
        pjm1 = polynomials%item(j)
        pj = polynomials%item(j + 1)
        pjp1 = (monomial * pj - fj * pjm1) / cj
      end if
      polynomials%item(j + 2) = pjp1
    end do
  end function orthogonal_polynomials_from_recurrence

  function orthonormal_polynomials_from_recurrence(recurrence, p0) result(polynomials)
    type(recurrence_t), intent(in) :: recurrence
    type(polynomial_t), intent(in) :: p0
    type(polylist_t) :: polynomials
    integer :: n, j
    type(polynomial_t) :: pjm1, pj, pjp1, monomial
    real(dp) :: cj, dj, ej, fj

    call validate_recurrence(recurrence)
    n = recurrence%size() - 1
    allocate(polynomials%item(n + 1))
    polynomials%item(1) = p0
    do j = 0, n - 1
      cj = recurrence%c(j + 1)
      dj = recurrence%d(j + 1)
      ej = recurrence%e(j + 1)
      fj = recurrence%f(j + 1)
      if (abs(cj) <= tiny(1.0_dp)) error stop 'orthonormal_polynomials: zero recurrence c coefficient'
      monomial = polynomial([dj, ej])
      if (j == 0) then
        pjp1 = (monomial * p0) / cj
      else
        pjm1 = polynomials%item(j)
        pj = polynomials%item(j + 1)
        pjp1 = (monomial * pj - fj * pjm1) / cj
      end if
      polynomials%item(j + 2) = pjp1
    end do
  end function orthonormal_polynomials_from_recurrence

  function monic_polynomial_recurrences(recurrence) result(monic)
    type(recurrence_t), intent(in) :: recurrence
    type(monic_recurrence_t) :: monic
    integer :: n, j

    call validate_recurrence(recurrence)
    n = recurrence%size()
    allocate(monic%a(n), monic%b(n))
    do j = 1, n
      if (abs(recurrence%e(j)) <= tiny(1.0_dp)) then
        monic%a(j) = 0.0_dp
      else
        monic%a(j) = -recurrence%d(j) / recurrence%e(j)
      end if
      if (j == 1) then
        monic%b(j) = 0.0_dp
      else if (abs(recurrence%e(j - 1)) <= tiny(1.0_dp) .or. &
          abs(recurrence%e(j)) <= tiny(1.0_dp)) then
        monic%b(j) = 0.0_dp
      else
        monic%b(j) = recurrence%c(j - 1) * recurrence%f(j) / &
          (recurrence%e(j - 1) * recurrence%e(j))
      end if
    end do
  end function monic_polynomial_recurrences

  function monic_polynomials(monic) result(polynomials)
    type(monic_recurrence_t), intent(in) :: monic
    type(polylist_t) :: polynomials
    integer :: n, j
    type(polynomial_t) :: p0, pjm1, pj, pjp1, monomial

    call validate_monic_recurrence(monic)
    n = monic%size() - 1
    allocate(polynomials%item(n + 1))
    p0 = polynomial(1.0_dp)
    polynomials%item(1) = p0
    do j = 0, n - 1
      monomial = polynomial([-monic%a(j + 1), 1.0_dp])
      if (j == 0) then
        pjp1 = monomial
      else
        pjm1 = polynomials%item(j)
        pj = polynomials%item(j + 1)
        pjp1 = monomial * pj - monic%b(j + 1) * pjm1
      end if
      polynomials%item(j + 2) = pjp1
    end do
  end function monic_polynomials

  function jacobi_matrices(monic) result(matrices)
    type(monic_recurrence_t), intent(in) :: monic
    type(real_matrix_list_t) :: matrices
    real(dp), allocatable :: full(:,:)
    integer :: n, i, j

    call validate_monic_recurrence(monic)
    n = monic%size() - 1
    if (n <= 0) then
      allocate(matrices%item(0))
      return
    end if
    allocate(full(n, n), source=0.0_dp)
    do i = 1, n
      full(i, i) = monic%a(i)
    end do
    do i = 2, n
      if (monic%b(i) < 0.0_dp) error stop 'jacobi_matrices: negative monic b coefficient'
      full(i, i - 1) = sqrt(monic%b(i))
      full(i - 1, i) = full(i, i - 1)
    end do
    allocate(matrices%item(n))
    do j = 1, n
      allocate(matrices%item(j)%value(j, j))
      matrices%item(j)%value = full(:j, :j)
    end do
  end function jacobi_matrices

  function polynomial_coefficients(polynomials) result(coefficients)
    type(polylist_t), intent(in) :: polynomials
    type(real_vector_list_t) :: coefficients
    integer :: j, n

    n = polynomials%size()
    allocate(coefficients%item(n))
    do j = 1, n
      coefficients%item(j)%value = polynomials%item(j)%coefficients()
    end do
  end function polynomial_coefficients

  function polynomial_derivatives(polynomials) result(derivatives)
    type(polylist_t), intent(in) :: polynomials
    type(polylist_t) :: derivatives
    integer :: j, n

    n = polynomials%size()
    allocate(derivatives%item(n))
    do j = 1, n
      derivatives%item(j) = derivative(polynomials%item(j))
    end do
  end function polynomial_derivatives

  function polynomial_integrals(polynomials) result(integrals)
    type(polylist_t), intent(in) :: polynomials
    type(polylist_t) :: integrals
    integer :: j, n

    n = polynomials%size()
    allocate(integrals%item(n))
    do j = 1, n
      integrals%item(j) = integral_polynomial(polynomials%item(j))
    end do
  end function polynomial_integrals

  function polynomial_orders(polynomials) result(orders)
    type(polylist_t), intent(in) :: polynomials
    integer, allocatable :: orders(:)
    integer :: j, n

    n = polynomials%size()
    allocate(orders(n))
    do j = 1, n
      orders(j) = polynomials%item(j)%degree()
    end do
  end function polynomial_orders

  function polynomial_powers(polynomials) result(powers)
    type(polylist_t), intent(in) :: polynomials
    type(real_vector_list_t) :: powers
    real(dp), allocatable :: lambda(:,:), inverse(:,:)
    integer :: n, j, m

    n = polynomials%size()
    allocate(lambda(n, n), source=0.0_dp)
    do j = 1, n
      m = min(j, size(polynomials%item(j)%coef))
      lambda(j, 1:m) = polynomials%item(j)%coef(:m)
    end do
    call invert_lower_triangular(lambda, inverse)
    allocate(powers%item(n))
    do j = 1, n
      allocate(powers%item(j)%value(j))
      powers%item(j)%value = inverse(j, :j)
    end do
  end function polynomial_powers

  function polynomial_roots_from_recurrence(monic) result(roots)
    type(monic_recurrence_t), intent(in) :: monic
    type(real_vector_list_t) :: roots
    type(polylist_t) :: polys
    complex(dp), allocatable :: z(:)
    integer :: j, n

    polys = monic_polynomials(monic)
    n = polys%size() - 1
    allocate(roots%item(n + 1))
    allocate(roots%item(1)%value(0))
    do j = 1, n
      z = polynom_roots(polys%item(j + 1))
      allocate(roots%item(j + 1)%value(size(z)))
      roots%item(j + 1)%value = real(z, dp)
      call sort_real(roots%item(j + 1)%value)
    end do
  end function polynomial_roots_from_recurrence

  function polynomial_values(polynomials, x) result(values)
    type(polylist_t), intent(in) :: polynomials
    real(dp), intent(in) :: x(:)
    type(real_vector_list_t) :: values
    integer :: j, n

    n = polynomials%size()
    allocate(values%item(n))
    do j = 1, n
      values%item(j)%value = polynomials%item(j)%evaluate(x)
    end do
  end function polynomial_values

  function polynomial_functions(polynomials) result(functions)
    type(polylist_t), intent(in) :: polynomials
    type(polynomial_function_list_t) :: functions
    integer :: j, n

    n = polynomials%size()
    allocate(functions%item(n))
    do j = 1, n
      functions%item(j)%polynomial = polynomials%item(j)
    end do
  end function polynomial_functions

  function scale_x(x, u, v, a, b) result(y)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: u, v
    real(dp), intent(in), optional :: a, b
    real(dp) :: y(size(x))
    real(dp) :: aa, bb

    if (size(x) == 0) return
    aa = minval(x)
    bb = maxval(x)
    if (present(a)) aa = a
    if (present(b)) bb = b
    if (aa > bb) error stop 'scale_x: domain lower bound exceeds upper bound'
    if (u > v) error stop 'scale_x: target lower bound exceeds upper bound'

    if (.not. ieee_is_finite(u) .and. .not. ieee_is_finite(v)) then
      y = x
    else if (.not. ieee_is_finite(u)) then
      y = x + v - bb
    else if (.not. ieee_is_finite(v)) then
      y = x + u - aa
    else
      if (abs(bb - aa) <= tiny(1.0_dp)) error stop 'scale_x: finite rescaling requires a nonzero domain width'
      y = u + (v - u) * (x - aa) / (bb - aa)
    end if
  end function scale_x

  subroutine validate_recurrence(recurrence)
    type(recurrence_t), intent(in) :: recurrence
    integer :: n

    if (.not. allocated(recurrence%c) .or. .not. allocated(recurrence%d) .or. &
        .not. allocated(recurrence%e) .or. .not. allocated(recurrence%f)) then
      error stop 'recurrence arrays are not allocated'
    end if
    n = size(recurrence%c)
    if (size(recurrence%d) /= n .or. size(recurrence%e) /= n .or. size(recurrence%f) /= n) then
      error stop 'recurrence arrays have inconsistent lengths'
    end if
  end subroutine validate_recurrence

  subroutine validate_monic_recurrence(monic)
    type(monic_recurrence_t), intent(in) :: monic
    if (.not. allocated(monic%a) .or. .not. allocated(monic%b)) then
      error stop 'monic recurrence arrays are not allocated'
    end if
    if (size(monic%a) /= size(monic%b)) error stop 'monic recurrence arrays have inconsistent lengths'
  end subroutine validate_monic_recurrence

  subroutine invert_lower_triangular(a, inverse)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: inverse(:,:)
    integer :: n, i, j, k
    real(dp) :: s

    n = size(a, 1)
    if (size(a, 2) /= n) error stop 'invert_lower_triangular: matrix is not square'
    allocate(inverse(n, n), source=0.0_dp)
    do j = 1, n
      do i = 1, n
        if (i < j) cycle
        if (abs(a(i, i)) <= tiny(1.0_dp)) error stop 'polynomial_powers: singular coefficient matrix'
        s = merge(1.0_dp, 0.0_dp, i == j)
        do k = j, i - 1
          s = s - a(i, k) * inverse(k, j)
        end do
        inverse(i, j) = s / a(i, i)
      end do
    end do
  end subroutine invert_lower_triangular

  subroutine sort_real(x)
    real(dp), intent(inout) :: x(:)
    integer :: i, j
    real(dp) :: key
    do i = 2, size(x)
      key = x(i)
      j = i - 1
      do while (j >= 1)
        if (x(j) <= key) exit
        x(j + 1) = x(j)
        j = j - 1
      end do
      x(j + 1) = key
    end do
  end subroutine sort_real

end module orthopolynom_core
