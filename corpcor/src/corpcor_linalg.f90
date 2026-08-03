! SPDX-License-Identifier: GPL-3.0-or-later
module corpcor_linalg
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use corpcor_kinds, only : dp
  use corpcor_types, only : svd_result, rank_condition_result, corpcor_success, &
    corpcor_invalid_argument, corpcor_dimension_error, corpcor_numerical_error
  implicit none
  private

  public :: identity_matrix, symmetric_eigen, fast_svd, pseudoinverse
  public :: matrix_power, rank_condition, is_positive_definite
  public :: make_positive_definite, covariance_to_correlation

contains

  pure function identity_matrix(n) result(a)
    integer, intent(in) :: n
    real(dp) :: a(n, n)
    integer :: i
    a = 0.0_dp
    do i = 1, n
      a(i, i) = 1.0_dp
    end do
  end function identity_matrix

  subroutine symmetric_eigen(a, values, vectors, status, tolerance, max_sweeps)
    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable, intent(out) :: values(:)
    real(dp), allocatable, intent(out) :: vectors(:, :)
    integer, intent(out), optional :: status
    real(dp), intent(in), optional :: tolerance
    integer, intent(in), optional :: max_sweeps
    real(dp), allocatable :: b(:, :)
    real(dp) :: app, aqq, apq, tau, t, c, s, bpj, bqj, vjp, vjq
    real(dp) :: tol, max_off
    integer :: n, p, q, j, sweep, nsweep, istat

    istat = corpcor_success
    n = size(a, 1)
    if (n < 1 .or. size(a, 2) /= n) then
      allocate(values(0), vectors(0, 0))
      istat = corpcor_dimension_error
      if (present(status)) status = istat
      return
    end if
    if (maxval(abs(a - transpose(a))) > 1000.0_dp * epsilon(1.0_dp) * &
        max(1.0_dp, maxval(abs(a)))) then
      allocate(values(0), vectors(0, 0))
      istat = corpcor_invalid_argument
      if (present(status)) status = istat
      return
    end if

    allocate(b(n, n), values(n), vectors(n, n))
    b = 0.5_dp * (a + transpose(a))
    vectors = identity_matrix(n)
    tol = 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(b)))
    if (present(tolerance)) tol = max(tolerance, 0.0_dp)
    nsweep = max(50, 20 * n * n)
    if (present(max_sweeps)) nsweep = max(1, max_sweeps)

    do sweep = 1, nsweep
      max_off = 0.0_dp
      do p = 1, n - 1
        do q = p + 1, n
          apq = b(p, q)
          max_off = max(max_off, abs(apq))
          if (abs(apq) <= tol) cycle
          app = b(p, p)
          aqq = b(q, q)
          tau = (aqq - app) / (2.0_dp * apq)
          if (tau >= 0.0_dp) then
            t = 1.0_dp / (tau + sqrt(1.0_dp + tau * tau))
          else
            t = -1.0_dp / (-tau + sqrt(1.0_dp + tau * tau))
          end if
          c = 1.0_dp / sqrt(1.0_dp + t * t)
          s = t * c

          do j = 1, n
            if (j /= p .and. j /= q) then
              bpj = b(p, j)
              bqj = b(q, j)
              b(p, j) = c * bpj - s * bqj
              b(j, p) = b(p, j)
              b(q, j) = s * bpj + c * bqj
              b(j, q) = b(q, j)
            end if
          end do
          b(p, p) = c*c*app - 2.0_dp*s*c*apq + s*s*aqq
          b(q, q) = s*s*app + 2.0_dp*s*c*apq + c*c*aqq
          b(p, q) = 0.0_dp
          b(q, p) = 0.0_dp
          do j = 1, n
            vjp = vectors(j, p)
            vjq = vectors(j, q)
            vectors(j, p) = c * vjp - s * vjq
            vectors(j, q) = s * vjp + c * vjq
          end do
        end do
      end do
      if (max_off <= tol) exit
    end do
    if (max_off > 10.0_dp * tol) istat = corpcor_numerical_error
    do j = 1, n
      values(j) = b(j, j)
    end do
    call sort_eigen_descending(values, vectors)
    if (present(status)) status = istat
  end subroutine symmetric_eigen

  subroutine sort_eigen_descending(values, vectors)
    real(dp), intent(inout) :: values(:)
    real(dp), intent(inout) :: vectors(:, :)
    real(dp) :: tv
    real(dp), allocatable :: col(:)
    integer :: i, j, k, n
    n = size(values)
    allocate(col(size(vectors, 1)))
    do i = 1, n - 1
      k = i
      do j = i + 1, n
        if (values(j) > values(k)) k = j
      end do
      if (k /= i) then
        tv = values(i)
        values(i) = values(k)
        values(k) = tv
        col = vectors(:, i)
        vectors(:, i) = vectors(:, k)
        vectors(:, k) = col
      end if
    end do
  end subroutine sort_eigen_descending

  function fast_svd(a, tolerance) result(res)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(in), optional :: tolerance
    type(svd_result) :: res
    real(dp), allocatable :: eval(:), evec(:, :), gram(:, :)
    real(dp) :: tol, maxd
    integer :: m, n, rmax, i, k, istat
    logical, allocatable :: keep(:)

    m = size(a, 1)
    n = size(a, 2)
    if (m < 1 .or. n < 1) then
      allocate(res%d(0), res%u(m, 0), res%v(n, 0))
      res%status = corpcor_dimension_error
      return
    end if

    if (m >= n) then
      gram = matmul(transpose(a), a)
      call symmetric_eigen(gram, eval, evec, istat)
      rmax = n
      maxd = sqrt(max(0.0_dp, eval(1)))
      tol = real(max(m, n), dp) * maxd * epsilon(1.0_dp)
      if (present(tolerance)) tol = max(tolerance, 0.0_dp)
      allocate(keep(rmax))
      keep = sqrt(max(0.0_dp, eval)) > tol
      k = count(keep)
      allocate(res%d(k), res%u(m, k), res%v(n, k))
      k = 0
      do i = 1, rmax
        if (.not. keep(i)) cycle
        k = k + 1
        res%d(k) = sqrt(max(0.0_dp, eval(i)))
        res%v(:, k) = evec(:, i)
        res%u(:, k) = matmul(a, res%v(:, k)) / res%d(k)
      end do
    else
      gram = matmul(a, transpose(a))
      call symmetric_eigen(gram, eval, evec, istat)
      rmax = m
      maxd = sqrt(max(0.0_dp, eval(1)))
      tol = real(max(m, n), dp) * maxd * epsilon(1.0_dp)
      if (present(tolerance)) tol = max(tolerance, 0.0_dp)
      allocate(keep(rmax))
      keep = sqrt(max(0.0_dp, eval)) > tol
      k = count(keep)
      allocate(res%d(k), res%u(m, k), res%v(n, k))
      k = 0
      do i = 1, rmax
        if (.not. keep(i)) cycle
        k = k + 1
        res%d(k) = sqrt(max(0.0_dp, eval(i)))
        res%u(:, k) = evec(:, i)
        res%v(:, k) = matmul(transpose(a), res%u(:, k)) / res%d(k)
      end do
    end if
    res%tol = tol
    res%rank = size(res%d)
    res%status = istat
  end function fast_svd

  function pseudoinverse(a, tolerance, status) result(pinv)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(in), optional :: tolerance
    integer, intent(out), optional :: status
    real(dp), allocatable :: pinv(:, :)
    type(svd_result) :: s
    integer :: i

    s = fast_svd(a, tolerance)
    allocate(pinv(size(a, 2), size(a, 1)))
    pinv = 0.0_dp
    do i = 1, s%rank
      pinv = pinv + spread(s%v(:, i) / s%d(i), 2, size(a, 1)) * &
        spread(s%u(:, i), 1, size(a, 2))
    end do
    if (present(status)) status = s%status
  end function pseudoinverse

  function matrix_power(a, alpha, pseudo, tolerance, status) result(power)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(in) :: alpha
    logical, intent(in), optional :: pseudo
    real(dp), intent(in), optional :: tolerance
    integer, intent(out), optional :: status
    real(dp), allocatable :: power(:, :)
    real(dp), allocatable :: eval(:), evec(:, :)
    real(dp) :: tol, evpow, nearest
    integer :: n, i, istat
    logical :: use_pseudo

    n = size(a, 1)
    allocate(power(n, size(a, 2)))
    power = 0.0_dp
    if (size(a, 2) /= n) then
      istat = corpcor_dimension_error
      if (present(status)) status = istat
      return
    end if
    call symmetric_eigen(a, eval, evec, istat)
    if (istat == corpcor_invalid_argument .or. size(eval) == 0) then
      if (present(status)) status = istat
      return
    end if
    tol = real(n, dp) * maxval(abs(eval)) * epsilon(1.0_dp)
    if (present(tolerance)) tol = max(tolerance, 0.0_dp)
    where (abs(eval) <= tol) eval = 0.0_dp
    use_pseudo = .false.
    if (present(pseudo)) use_pseudo = pseudo

    do i = 1, n
      if (use_pseudo .and. abs(eval(i)) <= tol) cycle
      if (abs(eval(i)) <= tol .and. alpha < 0.0_dp) then
        power = ieee_value(0.0_dp, ieee_quiet_nan)
        istat = corpcor_numerical_error
        exit
      end if
      if (eval(i) < 0.0_dp) then
        nearest = anint(alpha)
        if (abs(alpha - nearest) > 100.0_dp * epsilon(1.0_dp)) then
          power = ieee_value(0.0_dp, ieee_quiet_nan)
          istat = corpcor_numerical_error
          exit
        end if
      end if
      evpow = eval(i) ** alpha
      power = power + evpow * spread(evec(:, i), 2, n) * spread(evec(:, i), 1, n)
    end do
    if (present(status)) status = istat
  end function matrix_power

  function rank_condition(a, tolerance) result(res)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(in), optional :: tolerance
    type(rank_condition_result) :: res
    type(svd_result) :: s
    real(dp) :: maxd, mind

    s = fast_svd(a, tolerance)
    res%status = s%status
    res%rank = s%rank
    res%tol = s%tol
    if (s%rank == 0) then
      res%condition = huge(1.0_dp)
      return
    end if
    maxd = s%d(1)
    if (s%rank < min(size(a, 1), size(a, 2))) then
      res%condition = huge(1.0_dp)
    else
      mind = s%d(s%rank)
      if (mind <= 0.0_dp) then
        res%condition = huge(1.0_dp)
      else
        res%condition = maxd / mind
      end if
    end if
  end function rank_condition

  function is_positive_definite(a, tolerance, status) result(ok)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(in), optional :: tolerance
    integer, intent(out), optional :: status
    logical :: ok
    real(dp), allocatable :: eval(:), evec(:, :)
    real(dp) :: tol
    integer :: istat, n

    n = size(a, 1)
    ok = .false.
    if (size(a, 2) /= n) then
      if (present(status)) status = corpcor_dimension_error
      return
    end if
    call symmetric_eigen(a, eval, evec, istat)
    if (size(eval) == 0) then
      if (present(status)) status = istat
      return
    end if
    tol = real(n, dp) * maxval(abs(eval)) * epsilon(1.0_dp)
    if (present(tolerance)) tol = max(tolerance, 0.0_dp)
    ok = all(eval > tol)
    if (present(status)) status = istat
  end function is_positive_definite

  function make_positive_definite(a, tolerance, status) result(out)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(in), optional :: tolerance
    integer, intent(out), optional :: status
    real(dp), allocatable :: out(:, :)
    real(dp), allocatable :: eval(:), evec(:, :)
    real(dp) :: tol, delta, tau
    integer :: i, n, istat

    n = size(a, 1)
    allocate(out(n, size(a, 2)))
    out = a
    if (size(a, 2) /= n) then
      if (present(status)) status = corpcor_dimension_error
      return
    end if
    call symmetric_eigen(a, eval, evec, istat)
    tol = real(n, dp) * max(1.0_dp, maxval(abs(eval))) * epsilon(1.0_dp)
    if (present(tolerance)) tol = max(tolerance, 0.0_dp)
    delta = 2.0_dp * tol
    do i = 1, n
      tau = max(0.0_dp, delta - eval(i))
      out = out + tau * spread(evec(:, i), 2, n) * spread(evec(:, i), 1, n)
    end do
    out = 0.5_dp * (out + transpose(out))
    if (present(status)) status = istat
  end function make_positive_definite

  function covariance_to_correlation(a, status) result(r)
    real(dp), intent(in) :: a(:, :)
    integer, intent(out), optional :: status
    real(dp), allocatable :: r(:, :)
    real(dp), allocatable :: sd(:)
    integer :: i, j, n, istat

    n = size(a, 1)
    allocate(r(n, size(a, 2)))
    r = 0.0_dp
    istat = corpcor_success
    if (size(a, 2) /= n) then
      istat = corpcor_dimension_error
      if (present(status)) status = istat
      return
    end if
    allocate(sd(n))
    do i = 1, n
      sd(i) = sqrt(max(0.0_dp, a(i, i)))
    end do
    do j = 1, n
      do i = 1, n
        if (sd(i) > 0.0_dp .and. sd(j) > 0.0_dp) then
          r(i, j) = a(i, j) / (sd(i) * sd(j))
        else
          r(i, j) = 0.0_dp
        end if
      end do
      if (sd(j) > 0.0_dp) r(j, j) = 1.0_dp
    end do
    if (present(status)) status = istat
  end function covariance_to_correlation

end module corpcor_linalg
