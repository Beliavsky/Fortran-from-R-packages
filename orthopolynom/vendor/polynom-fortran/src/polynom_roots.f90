module polynom_roots
  use polynom_kinds, only : dp
  use polynom_status, only : poly_status_t, poly_ok, poly_root_failure, set_status
  use polynom_core, only : polynomial_t, polynomial, operator(-)
  use polynom_algorithms, only : derivative
  implicit none
  private

  type, public :: polynomial_summary_t
    complex(dp), allocatable :: zeros(:)
    complex(dp), allocatable :: stationary_points(:)
    complex(dp), allocatable :: inflexion_points(:)
  end type polynomial_summary_t

  public :: polynomial_roots, summarize_polynomial

contains

  function polynomial_roots(p, intercept, status, tol, max_iter) result(roots)
    type(polynomial_t), intent(in) :: p
    real(dp), intent(in), optional :: intercept
    type(poly_status_t), intent(out), optional :: status
    real(dp), intent(in), optional :: tol
    integer, intent(in), optional :: max_iter
    complex(dp), allocatable :: roots(:)
    type(polynomial_t) :: q
    real(dp) :: eps, radius, pi
    complex(dp), allocatable :: z(:), znew(:)
    complex(dp) :: f, fp, corr, denom, s
    integer :: n, nzero, i, j, iter, maxit, active_degree
    logical :: converged

    q = p
    if (present(intercept)) q = p - intercept
    eps = 1.0e-12_dp
    if (present(tol)) eps = max(tol, 10.0_dp * epsilon(1.0_dp))
    maxit = 5000
    if (present(max_iter)) maxit = max(10, max_iter)

    nzero = 0
    do while (nzero < q%degree() .and. abs(q%coef(nzero + 1)) <= eps)
      nzero = nzero + 1
    end do
    active_degree = q%degree() - nzero
    if (active_degree <= 0) then
      allocate(roots(nzero), source=cmplx(0.0_dp, 0.0_dp, dp))
      call set_status(status, poly_ok, '')
      return
    end if

    if (active_degree == 1) then
      allocate(roots(nzero + 1), source=cmplx(0.0_dp, 0.0_dp, dp))
      roots(nzero + 1) = cmplx(-q%coef(nzero + 1) / q%coef(nzero + 2), 0.0_dp, dp)
      call sort_roots(roots)
      call set_status(status, poly_ok, '')
      return
    end if

    n = active_degree
    allocate(z(n), znew(n))
    radius = 1.0_dp + maxval(abs(q%coef(nzero + 1:nzero + n))) / abs(q%coef(nzero + n + 1))
    pi = acos(-1.0_dp)
    do i = 1, n
      z(i) = radius * exp(cmplx(0.0_dp, 2.0_dp * pi * (real(i - 1, dp) + 0.25_dp) / real(n, dp), dp))
    end do

    converged = .false.
    do iter = 1, maxit
      do i = 1, n
        call evaluate_complex_with_derivative(q%coef(nzero + 1:), z(i), f, fp)
        s = cmplx(0.0_dp, 0.0_dp, dp)
        do j = 1, n
          if (j == i) cycle
          if (abs(z(i) - z(j)) > tiny(1.0_dp)) s = s + 1.0_dp / (z(i) - z(j))
        end do
        denom = fp - f * s
        if (abs(denom) <= tiny(1.0_dp)) then
          corr = f / merge(fp, cmplx(1.0_dp, 0.0_dp, dp), abs(fp) > tiny(1.0_dp))
        else
          corr = f / denom
        end if
        znew(i) = z(i) - corr
      end do
      if (maxval(abs(znew - z)) <= eps * (1.0_dp + maxval(abs(znew)))) then
        converged = .true.
        z = znew
        exit
      end if
      z = znew
    end do

    allocate(roots(nzero + n), source=cmplx(0.0_dp, 0.0_dp, dp))
    roots(nzero + 1:) = z
    do i = 1, size(roots)
      if (abs(aimag(roots(i))) <= 100.0_dp * eps * (1.0_dp + abs(real(roots(i), dp)))) then
        roots(i) = cmplx(real(roots(i), dp), 0.0_dp, dp)
      end if
    end do
    call sort_roots(roots)
    if (converged) then
      call set_status(status, poly_ok, '')
    else
      call set_status(status, poly_root_failure, 'root iteration reached the maximum iteration count')
    end if
  end function polynomial_roots

  subroutine evaluate_complex_with_derivative(coef, z, value, derivative_value)
    real(dp), intent(in) :: coef(:)
    complex(dp), intent(in) :: z
    complex(dp), intent(out) :: value, derivative_value
    integer :: i
    value = cmplx(coef(size(coef)), 0.0_dp, dp)
    derivative_value = cmplx(0.0_dp, 0.0_dp, dp)
    do i = size(coef) - 1, 1, -1
      derivative_value = value + z * derivative_value
      value = cmplx(coef(i), 0.0_dp, dp) + z * value
    end do
  end subroutine evaluate_complex_with_derivative

  subroutine sort_roots(roots)
    complex(dp), intent(inout) :: roots(:)
    complex(dp) :: key
    integer :: i, j
    do i = 2, size(roots)
      key = roots(i)
      j = i - 1
      do while (j >= 1)
        if (.not. root_after(roots(j), key)) exit
        roots(j + 1) = roots(j)
        j = j - 1
      end do
      roots(j + 1) = key
    end do
  end subroutine sort_roots

  pure logical function root_after(a, b)
    complex(dp), intent(in) :: a, b
    if (real(a, dp) > real(b, dp)) then
      root_after = .true.
    else if (real(a, dp) < real(b, dp)) then
      root_after = .false.
    else
      root_after = aimag(a) > aimag(b)
    end if
  end function root_after

  function summarize_polynomial(p, status) result(summary)
    type(polynomial_t), intent(in) :: p
    type(poly_status_t), intent(out), optional :: status
    type(polynomial_summary_t) :: summary
    type(polynomial_t) :: d1, d2
    type(poly_status_t) :: local_status

    allocate(summary%zeros(0), summary%stationary_points(0), summary%inflexion_points(0))
    summary%zeros = polynomial_roots(p, status=local_status)
    d1 = derivative(p)
    summary%stationary_points = polynomial_roots(d1)
    d2 = derivative(d1)
    summary%inflexion_points = polynomial_roots(d2)
    call set_status(status, local_status%code, local_status%message)
  end function summarize_polynomial

end module polynom_roots
