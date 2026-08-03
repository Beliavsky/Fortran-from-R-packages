! SPDX-License-Identifier: GPL-2.0-or-later
module maxlik_linalg
  use maxlik_kinds, only: dp
  implicit none
  private

  public :: vector_norm, outer_product, identity_matrix
  public :: solve_linear, invert_matrix, symmetric_eigenvalues
  public :: symmetric_condition_number, rectangular_condition_number, is_negative_definite

contains

  pure real(dp) function vector_norm(x) result(value)
    real(dp), intent(in) :: x(:)
    value = sqrt(max(0.0_dp, dot_product(x, x)))
  end function vector_norm

  pure function outer_product(x, y) result(a)
    real(dp), intent(in) :: x(:), y(:)
    real(dp) :: a(size(x), size(y))
    integer :: j
    do j = 1, size(y)
      a(:, j) = x * y(j)
    end do
  end function outer_product

  pure function identity_matrix(n) result(a)
    integer, intent(in) :: n
    real(dp) :: a(n, n)
    integer :: i
    a = 0.0_dp
    do i = 1, n
      a(i, i) = 1.0_dp
    end do
  end function identity_matrix

  subroutine solve_linear(a, b, x, status, tolerance)
    real(dp), intent(in) :: a(:, :), b(:)
    real(dp), intent(out) :: x(:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: tolerance

    real(dp), allocatable :: work(:, :), rhs(:), row(:)
    real(dp) :: pivot, factor, tol
    integer :: i, k, p, n

    n = size(b)
    status = 1
    x = 0.0_dp
    if (size(a, 1) /= n .or. size(a, 2) /= n .or. size(x) /= n) return

    tol = 100.0_dp * epsilon(1.0_dp)
    if (present(tolerance)) tol = max(tolerance, epsilon(1.0_dp))
    allocate(work(n, n), rhs(n), row(n))
    work = a
    rhs = b

    do k = 1, n - 1
      p = k - 1 + maxloc(abs(work(k:n, k)), dim=1)
      pivot = abs(work(p, k))
      if (pivot <= tol * max(1.0_dp, maxval(abs(work)))) return
      if (p /= k) then
        row = work(k, :)
        work(k, :) = work(p, :)
        work(p, :) = row
        factor = rhs(k)
        rhs(k) = rhs(p)
        rhs(p) = factor
      end if
      do i = k + 1, n
        factor = work(i, k) / work(k, k)
        work(i, k) = 0.0_dp
        work(i, k + 1:n) = work(i, k + 1:n) - factor * work(k, k + 1:n)
        rhs(i) = rhs(i) - factor * rhs(k)
      end do
    end do

    if (abs(work(n, n)) <= tol * max(1.0_dp, maxval(abs(work)))) return
    do i = n, 1, -1
      if (i < n) then
        x(i) = (rhs(i) - dot_product(work(i, i + 1:n), x(i + 1:n))) / work(i, i)
      else
        x(i) = rhs(i) / work(i, i)
      end if
    end do
    status = 0
  end subroutine solve_linear

  subroutine invert_matrix(a, inverse, status, tolerance)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(out) :: inverse(:, :)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: tolerance

    real(dp), allocatable :: e(:), column(:)
    integer :: j, n, solve_status

    n = size(a, 1)
    status = 1
    inverse = 0.0_dp
    if (size(a, 2) /= n .or. size(inverse, 1) /= n .or. size(inverse, 2) /= n) return
    allocate(e(n), column(n))
    do j = 1, n
      e = 0.0_dp
      e(j) = 1.0_dp
      call solve_linear(a, e, column, solve_status, tolerance)
      if (solve_status /= 0) return
      inverse(:, j) = column
    end do
    status = 0
  end subroutine invert_matrix

  subroutine symmetric_eigenvalues(a, eigenvalues, status, tolerance, max_sweeps)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(out) :: eigenvalues(:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: tolerance
    integer, intent(in), optional :: max_sweeps

    real(dp), allocatable :: work(:, :)
    real(dp) :: app, aqq, apq, tau, t, c, s, wip, wiq, tol
    integer :: i, p, q, sweep, n, limit

    n = size(a, 1)
    status = 1
    eigenvalues = 0.0_dp
    if (size(a, 2) /= n .or. size(eigenvalues) /= n) return
    tol = 100.0_dp * epsilon(1.0_dp)
    if (present(tolerance)) tol = max(tolerance, epsilon(1.0_dp))
    limit = max(30, 10 * n * n)
    if (present(max_sweeps)) limit = max_sweeps

    allocate(work(n, n))
    work = 0.5_dp * (a + transpose(a))
    do sweep = 1, limit
      apq = 0.0_dp
      p = 1
      q = min(2, n)
      do i = 1, n - 1
        if (maxval(abs(work(i, i + 1:n))) > apq) then
          q = i + maxloc(abs(work(i, i + 1:n)), dim=1)
          p = i
          apq = abs(work(p, q))
        end if
      end do
      if (n <= 1 .or. apq <= tol * max(1.0_dp, maxval(abs(work)))) then
        do i = 1, n
          eigenvalues(i) = work(i, i)
        end do
        call sort_ascending(eigenvalues)
        status = 0
        return
      end if

      app = work(p, p)
      aqq = work(q, q)
      tau = (aqq - app) / (2.0_dp * work(p, q))
      if (tau >= 0.0_dp) then
        t = 1.0_dp / (tau + sqrt(1.0_dp + tau * tau))
      else
        t = -1.0_dp / (-tau + sqrt(1.0_dp + tau * tau))
      end if
      c = 1.0_dp / sqrt(1.0_dp + t * t)
      s = t * c

      do i = 1, n
        if (i == p .or. i == q) cycle
        wip = work(i, p)
        wiq = work(i, q)
        work(i, p) = c * wip - s * wiq
        work(p, i) = work(i, p)
        work(i, q) = s * wip + c * wiq
        work(q, i) = work(i, q)
      end do
      work(p, p) = c * c * app - 2.0_dp * s * c * work(p, q) + s * s * aqq
      work(q, q) = s * s * app + 2.0_dp * s * c * work(p, q) + c * c * aqq
      work(p, q) = 0.0_dp
      work(q, p) = 0.0_dp
    end do

    do i = 1, n
      eigenvalues(i) = work(i, i)
    end do
  end subroutine symmetric_eigenvalues

  pure subroutine sort_ascending(x)
    real(dp), intent(inout) :: x(:)
    real(dp) :: key
    integer :: i, j
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
  end subroutine sort_ascending

  real(dp) function symmetric_condition_number(a, status, tolerance) result(value)
    real(dp), intent(in) :: a(:, :)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: tolerance

    real(dp), allocatable :: eigenvalues(:)
    real(dp) :: largest, smallest, tol

    allocate(eigenvalues(size(a, 1)))
    call symmetric_eigenvalues(a, eigenvalues, status)
    if (status /= 0) then
      value = huge(1.0_dp)
      return
    end if
    tol = sqrt(epsilon(1.0_dp))
    if (present(tolerance)) tol = tolerance
    largest = maxval(abs(eigenvalues))
    if (count(abs(eigenvalues) > tol * max(1.0_dp, largest)) == 0) then
      value = huge(1.0_dp)
      status = 1
      return
    end if
    smallest = minval(abs(eigenvalues), mask=abs(eigenvalues) > tol * max(1.0_dp, largest))
    if (largest <= 0.0_dp .or. smallest <= 0.0_dp .or. &
        count(abs(eigenvalues) > tol * max(1.0_dp, largest)) < size(eigenvalues)) then
      value = huge(1.0_dp)
      status = 1
    else
      value = largest / smallest
      status = 0
    end if
  end function symmetric_condition_number

  real(dp) function rectangular_condition_number(a, status, normalize) result(value)
    real(dp), intent(in) :: a(:, :)
    integer, intent(out) :: status
    logical, intent(in), optional :: normalize

    real(dp), allocatable :: work(:, :), gram(:, :), norms(:)
    logical :: do_normalize
    integer :: j

    do_normalize = .false.
    if (present(normalize)) do_normalize = normalize
    allocate(work(size(a, 1), size(a, 2)))
    work = a
    if (do_normalize) then
      allocate(norms(size(a, 2)))
      do j = 1, size(a, 2)
        norms(j) = vector_norm(work(:, j))
        if (norms(j) <= tiny(1.0_dp)) then
          value = huge(1.0_dp)
          status = 1
          return
        end if
        work(:, j) = work(:, j) / norms(j)
      end do
    end if
    allocate(gram(size(a, 2), size(a, 2)))
    gram = matmul(transpose(work), work)
    value = symmetric_condition_number(gram, status)
    if (status == 0) value = sqrt(value)
  end function rectangular_condition_number

  logical function is_negative_definite(a, tolerance) result(ok)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(in), optional :: tolerance

    real(dp), allocatable :: eigenvalues(:)
    real(dp) :: tol
    integer :: status

    allocate(eigenvalues(size(a, 1)))
    call symmetric_eigenvalues(a, eigenvalues, status)
    tol = 100.0_dp * epsilon(1.0_dp)
    if (present(tolerance)) tol = tolerance
    ok = status == 0 .and. maxval(eigenvalues) < tol
  end function is_negative_definite

end module maxlik_linalg
