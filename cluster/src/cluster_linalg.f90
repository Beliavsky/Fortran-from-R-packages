! SPDX-License-Identifier: GPL-2.0-or-later
module cluster_linalg
  use fastcluster_kinds, only: dp
  implicit none
  private
  public :: inverse_matrix, determinant_spd, symmetric_eigen, sort_eigenpairs

contains

  subroutine inverse_matrix(a, ainv, status)
    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable, intent(out) :: ainv(:, :)
    integer, intent(out) :: status

    real(dp), allocatable :: aug(:, :), rowtmp(:)
    real(dp) :: pivot, factor, scale
    integer :: i, k, n, p

    status = 0
    n = size(a, 1)
    if (n < 1 .or. size(a, 2) /= n) then
      allocate(ainv(0, 0))
      status = 1
      return
    end if
    allocate(aug(n, 2*n), rowtmp(2*n), ainv(n, n))
    aug = 0.0_dp
    aug(:, 1:n) = a
    do i = 1, n
      aug(i, n+i) = 1.0_dp
    end do
    do k = 1, n
      p = k
      scale = abs(aug(k, k))
      do i = k + 1, n
        if (abs(aug(i, k)) > scale) then
          scale = abs(aug(i, k))
          p = i
        end if
      end do
      if (scale <= 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(a)))) then
        status = 2
        ainv = 0.0_dp
        return
      end if
      if (p /= k) then
        rowtmp = aug(k, :)
        aug(k, :) = aug(p, :)
        aug(p, :) = rowtmp
      end if
      pivot = aug(k, k)
      aug(k, :) = aug(k, :) / pivot
      do i = 1, n
        if (i == k) cycle
        factor = aug(i, k)
        if (abs(factor) > tiny(1.0_dp)) aug(i, :) = aug(i, :) - factor * aug(k, :)
      end do
    end do
    ainv = aug(:, n+1:2*n)
  end subroutine inverse_matrix

  subroutine determinant_spd(a, determinant, status)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(out) :: determinant
    integer, intent(out) :: status

    real(dp), allocatable :: l(:, :)
    real(dp) :: value
    integer :: i, j, k, n

    status = 0
    determinant = 0.0_dp
    n = size(a, 1)
    if (n < 1 .or. size(a, 2) /= n) then
      status = 1
      return
    end if
    allocate(l(n, n))
    l = 0.0_dp
    do i = 1, n
      do j = 1, i
        value = a(i, j)
        do k = 1, j - 1
          value = value - l(i, k) * l(j, k)
        end do
        if (i == j) then
          if (value <= 0.0_dp) then
            status = 2
            return
          end if
          l(i, j) = sqrt(value)
        else
          l(i, j) = value / l(j, j)
        end if
      end do
    end do
    determinant = 1.0_dp
    do i = 1, n
      determinant = determinant * l(i, i) * l(i, i)
    end do
  end subroutine determinant_spd

  subroutine symmetric_eigen(a, values, vectors, status, tolerance, max_sweeps)
    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable, intent(out) :: values(:), vectors(:, :)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: tolerance
    integer, intent(in), optional :: max_sweeps

    real(dp), allocatable :: work(:, :)
    real(dp) :: app, aqq, apq, c, s, tau, t, tol, off, vip, viq, wip, wiq
    integer :: i, j, n, p, q, sweep, nsweep

    status = 0
    n = size(a, 1)
    if (n < 1 .or. size(a, 2) /= n) then
      allocate(values(0), vectors(0, 0))
      status = 1
      return
    end if
    tol = 100.0_dp * epsilon(1.0_dp)
    if (present(tolerance)) tol = tolerance
    nsweep = max(50, 10*n*n)
    if (present(max_sweeps)) nsweep = max_sweeps
    allocate(work(n, n), values(n), vectors(n, n))
    work = 0.5_dp * (a + transpose(a))
    vectors = 0.0_dp
    do i = 1, n
      vectors(i, i) = 1.0_dp
    end do
    do sweep = 1, nsweep
      off = 0.0_dp
      p = 1
      q = min(2, n)
      do j = 2, n
        do i = 1, j - 1
          if (abs(work(i, j)) > off) then
            off = abs(work(i, j))
            p = i
            q = j
          end if
        end do
      end do
      if (off <= tol * max(1.0_dp, maxval(abs(work)))) exit
      app = work(p, p)
      aqq = work(q, q)
      apq = work(p, q)
      tau = (aqq - app) / (2.0_dp * apq)
      if (tau >= 0.0_dp) then
        t = 1.0_dp / (tau + sqrt(1.0_dp + tau*tau))
      else
        t = -1.0_dp / (-tau + sqrt(1.0_dp + tau*tau))
      end if
      c = 1.0_dp / sqrt(1.0_dp + t*t)
      s = t*c
      do i = 1, n
        if (i /= p .and. i /= q) then
          wip = work(i, p)
          wiq = work(i, q)
          work(i, p) = c*wip - s*wiq
          work(p, i) = work(i, p)
          work(i, q) = s*wip + c*wiq
          work(q, i) = work(i, q)
        end if
      end do
      work(p, p) = c*c*app - 2.0_dp*s*c*apq + s*s*aqq
      work(q, q) = s*s*app + 2.0_dp*s*c*apq + c*c*aqq
      work(p, q) = 0.0_dp
      work(q, p) = 0.0_dp
      do i = 1, n
        vip = vectors(i, p)
        viq = vectors(i, q)
        vectors(i, p) = c*vip - s*viq
        vectors(i, q) = s*vip + c*viq
      end do
    end do
    if (sweep > nsweep) status = 2
    do i = 1, n
      values(i) = work(i, i)
    end do
    call sort_eigenpairs(values, vectors)
  end subroutine symmetric_eigen

  subroutine sort_eigenpairs(values, vectors)
    real(dp), intent(inout) :: values(:)
    real(dp), intent(inout) :: vectors(:, :)

    real(dp) :: v
    real(dp), allocatable :: col(:)
    integer :: i, j, n

    n = size(values)
    allocate(col(size(vectors, 1)))
    do i = 1, n - 1
      do j = i + 1, n
        if (values(j) < values(i)) then
          v = values(i)
          values(i) = values(j)
          values(j) = v
          col = vectors(:, i)
          vectors(:, i) = vectors(:, j)
          vectors(:, j) = col
        end if
      end do
    end do
  end subroutine sort_eigenpairs

end module cluster_linalg
