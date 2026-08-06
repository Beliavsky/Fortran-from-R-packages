! SPDX-License-Identifier: MIT
module uncorbets_linalg
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use uncorbets_kinds, only : dp
  use uncorbets_types, only : status_type, set_status, uncorbets_ok, &
      uncorbets_invalid_input, uncorbets_not_pos_semidefinite, &
      uncorbets_singular_matrix, uncorbets_no_convergence
  implicit none
  private

  public :: symmetric_eigen, symmetric_sqrt, solve_linear, inverse_matrix
  public :: frobenius_norm, project_simplex, is_finite_matrix, is_symmetric

contains

  logical function is_finite_matrix(a)
    real(dp), intent(in) :: a(:, :)
    is_finite_matrix = all(ieee_is_finite(a))
  end function is_finite_matrix

  logical function is_symmetric(a, tolerance)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(in), optional :: tolerance
    real(dp) :: tol
    if (size(a, 1) /= size(a, 2)) then
      is_symmetric = .false.
      return
    end if
    tol = 100.0_dp * epsilon(1.0_dp)
    if (present(tolerance)) tol = tolerance
    is_symmetric = maxval(abs(a - transpose(a))) <= tol * &
        max(1.0_dp, maxval(abs(a)))
  end function is_symmetric

  real(dp) function frobenius_norm(a)
    real(dp), intent(in) :: a(:, :)
    frobenius_norm = sqrt(sum(a * a))
  end function frobenius_norm

  subroutine symmetric_eigen(a, values, vectors, status, descending)
    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable, intent(out) :: values(:)
    real(dp), allocatable, intent(out) :: vectors(:, :)
    type(status_type), intent(out) :: status
    logical, intent(in), optional :: descending

    real(dp), allocatable :: work(:, :)
    real(dp) :: app, aqq, apq, tau, t, c, s, akp, akq, vkp, vkq
    real(dp) :: off, scale, tol, tmp_value
    real(dp), allocatable :: tmp_vector(:)
    integer :: n, i, j, k, p, q, sweep, max_sweeps
    logical :: desc

    call set_status(status, uncorbets_ok, 'ok')
    n = size(a, 1)
    if (n < 1 .or. size(a, 2) /= n) then
      call set_status(status, uncorbets_invalid_input, &
          'eigen input must be a nonempty square matrix')
      allocate(values(0), vectors(0, 0))
      return
    end if
    if (.not. is_finite_matrix(a) .or. .not. is_symmetric(a)) then
      call set_status(status, uncorbets_invalid_input, &
          'eigen input must be finite and symmetric')
      allocate(values(0), vectors(0, 0))
      return
    end if

    allocate(work(n, n), values(n), vectors(n, n), tmp_vector(n))
    work = 0.5_dp * (a + transpose(a))
    vectors = 0.0_dp
    do i = 1, n
      vectors(i, i) = 1.0_dp
    end do

    scale = max(1.0_dp, maxval(abs(work)))
    tol = 50.0_dp * epsilon(1.0_dp) * scale
    max_sweeps = max(50, 20 * n * n)

    do sweep = 1, max_sweeps
      off = 0.0_dp
      p = 1
      q = 1
      do j = 2, n
        do i = 1, j - 1
          if (abs(work(i, j)) > off) then
            off = abs(work(i, j))
            p = i
            q = j
          end if
        end do
      end do
      if (off <= tol) exit

      app = work(p, p)
      aqq = work(q, q)
      apq = work(p, q)
      tau = (aqq - app) / (2.0_dp * apq)
      if (tau >= 0.0_dp) then
        t = 1.0_dp / (tau + sqrt(1.0_dp + tau * tau))
      else
        t = -1.0_dp / (-tau + sqrt(1.0_dp + tau * tau))
      end if
      c = 1.0_dp / sqrt(1.0_dp + t * t)
      s = t * c

      do k = 1, n
        if (k /= p .and. k /= q) then
          akp = work(k, p)
          akq = work(k, q)
          work(k, p) = c * akp - s * akq
          work(p, k) = work(k, p)
          work(k, q) = s * akp + c * akq
          work(q, k) = work(k, q)
        end if
      end do
      work(p, p) = app - t * apq
      work(q, q) = aqq + t * apq
      work(p, q) = 0.0_dp
      work(q, p) = 0.0_dp

      do k = 1, n
        vkp = vectors(k, p)
        vkq = vectors(k, q)
        vectors(k, p) = c * vkp - s * vkq
        vectors(k, q) = s * vkp + c * vkq
      end do
    end do

    if (sweep > max_sweeps) then
      call set_status(status, uncorbets_no_convergence, &
          'Jacobi eigensolver did not converge')
      return
    end if

    do i = 1, n
      values(i) = work(i, i)
    end do

    desc = .false.
    if (present(descending)) desc = descending
    do i = 1, n - 1
      k = i
      do j = i + 1, n
        if ((desc .and. values(j) > values(k)) .or. &
            (.not. desc .and. values(j) < values(k))) k = j
      end do
      if (k /= i) then
        tmp_value = values(i)
        values(i) = values(k)
        values(k) = tmp_value
        tmp_vector = vectors(:, i)
        vectors(:, i) = vectors(:, k)
        vectors(:, k) = tmp_vector
      end if
    end do
  end subroutine symmetric_eigen

  subroutine symmetric_sqrt(a, root, status)
    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable, intent(out) :: root(:, :)
    type(status_type), intent(out) :: status

    real(dp), allocatable :: values(:), vectors(:, :), scaled(:, :)
    real(dp) :: tol, scale
    integer :: n, j

    call symmetric_eigen(a, values, vectors, status, descending=.true.)
    if (.not. status%ok()) then
      allocate(root(0, 0))
      return
    end if
    n = size(values)
    scale = max(1.0_dp, maxval(abs(values)))
    tol = 1000.0_dp * epsilon(1.0_dp) * scale
    if (minval(values) < -tol) then
      call set_status(status, uncorbets_not_pos_semidefinite, &
          'matrix has a negative eigenvalue')
      allocate(root(0, 0))
      return
    end if
    values = sqrt(max(values, 0.0_dp))
    allocate(scaled(n, n), root(n, n))
    scaled = vectors
    do j = 1, n
      scaled(:, j) = scaled(:, j) * values(j)
    end do
    root = matmul(scaled, transpose(vectors))
    root = 0.5_dp * (root + transpose(root))
    call set_status(status, uncorbets_ok, 'ok')
  end subroutine symmetric_sqrt

  subroutine solve_linear(a, b, x, status)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(in) :: b(:, :)
    real(dp), allocatable, intent(out) :: x(:, :)
    type(status_type), intent(out) :: status

    real(dp), allocatable :: aug(:, :), temp_row(:)
    real(dp) :: pivot_abs, factor, scale
    integer :: n, nrhs, i, j, k, pivot

    call set_status(status, uncorbets_ok, 'ok')
    n = size(a, 1)
    nrhs = size(b, 2)
    if (n < 1 .or. size(a, 2) /= n .or. size(b, 1) /= n) then
      call set_status(status, uncorbets_invalid_input, &
          'incompatible dimensions in linear solve')
      allocate(x(0, 0))
      return
    end if
    if (.not. all(ieee_is_finite(a)) .or. .not. all(ieee_is_finite(b))) then
      call set_status(status, uncorbets_invalid_input, &
          'linear solve inputs must be finite')
      allocate(x(0, 0))
      return
    end if

    allocate(aug(n, n + nrhs), temp_row(n + nrhs), x(n, nrhs))
    aug(:, 1:n) = a
    aug(:, n + 1:n + nrhs) = b
    scale = max(1.0_dp, maxval(abs(a)))

    do k = 1, n
      pivot = k
      pivot_abs = abs(aug(k, k))
      do i = k + 1, n
        if (abs(aug(i, k)) > pivot_abs) then
          pivot = i
          pivot_abs = abs(aug(i, k))
        end if
      end do
      if (pivot_abs <= 100.0_dp * epsilon(1.0_dp) * scale) then
        call set_status(status, uncorbets_singular_matrix, &
            'matrix is singular to working precision')
        deallocate(x)
        allocate(x(0, 0))
        return
      end if
      if (pivot /= k) then
        temp_row = aug(k, :)
        aug(k, :) = aug(pivot, :)
        aug(pivot, :) = temp_row
      end if

      aug(k, :) = aug(k, :) / aug(k, k)
      do i = 1, n
        if (i /= k) then
          factor = aug(i, k)
          aug(i, :) = aug(i, :) - factor * aug(k, :)
        end if
      end do
    end do

    do j = 1, nrhs
      x(:, j) = aug(:, n + j)
    end do
  end subroutine solve_linear

  subroutine inverse_matrix(a, inverse, status)
    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable, intent(out) :: inverse(:, :)
    type(status_type), intent(out) :: status
    real(dp), allocatable :: identity(:, :)
    integer :: n, i

    n = size(a, 1)
    if (n < 1 .or. size(a, 2) /= n) then
      call set_status(status, uncorbets_invalid_input, &
          'inverse input must be a nonempty square matrix')
      allocate(inverse(0, 0))
      return
    end if
    allocate(identity(n, n))
    identity = 0.0_dp
    do i = 1, n
      identity(i, i) = 1.0_dp
    end do
    call solve_linear(a, identity, inverse, status)
  end subroutine inverse_matrix

  subroutine project_simplex(y, x)
    real(dp), intent(in) :: y(:)
    real(dp), intent(out) :: x(:)

    real(dp), allocatable :: u(:)
    real(dp) :: cumulative, theta, tmp
    integer :: n, i, j, rho

    n = size(y)
    if (size(x) /= n) error stop 'project_simplex: size mismatch'
    allocate(u(n))
    u = y
    do i = 2, n
      tmp = u(i)
      j = i - 1
      do while (j >= 1)
        if (u(j) >= tmp) exit
        u(j + 1) = u(j)
        j = j - 1
      end do
      u(j + 1) = tmp
    end do

    cumulative = 0.0_dp
    rho = 1
    do i = 1, n
      cumulative = cumulative + u(i)
      if (u(i) - (cumulative - 1.0_dp) / real(i, dp) > 0.0_dp) rho = i
    end do
    theta = (sum(u(1:rho)) - 1.0_dp) / real(rho, dp)
    x = max(y - theta, 0.0_dp)
    if (sum(x) > 0.0_dp) x = x / sum(x)
  end subroutine project_simplex

end module uncorbets_linalg
