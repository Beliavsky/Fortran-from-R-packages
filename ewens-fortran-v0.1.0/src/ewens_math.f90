! SPDX-License-Identifier: MIT
module ewens_math
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_negative_inf, ieee_is_finite
  use ewens_kinds, only : dp, i8
  use copula, only : stirling_first
  implicit none
  private
  public :: log_rising_factorial, log_unsigned_stirling1
  public :: number_of_classes, class_size_spectrum

contains

  pure real(dp) function log_rising_factorial(theta, n) result(value)
    real(dp), intent(in) :: theta
    integer, intent(in) :: n
    integer :: j

    value = 0.0_dp
    do j = 0, n - 1
      value = value + log(theta + real(j, dp))
    end do
  end function log_rising_factorial

  real(dp) function log_unsigned_stirling1(n, k) result(value)
    integer, intent(in) :: n, k
    real(dp), allocatable :: old(:), new(:)
    real(dp) :: neginf
    integer(i8) :: exact
    integer :: i, j

    neginf = ieee_value(0.0_dp, ieee_negative_inf)
    if (n < 0 .or. k < 0 .or. k > n) then
      value = neginf
      return
    end if
    if (n == 0) then
      if (k == 0) then
        value = 0.0_dp
      else
        value = neginf
      end if
      return
    end if
    if (k == 0) then
      value = neginf
      return
    end if

    ! The attached copula translation returns exact int64 Stirling numbers.
    ! They are safe through n=20; use them directly in that range.
    if (n <= 20) then
      exact = stirling_first(n, k)
      if (exact > 0_i8) then
        value = log(real(exact, dp))
      else
        value = neginf
      end if
      return
    end if

    allocate(old(0:k), new(0:k))
    old = neginf
    old(0) = 0.0_dp

    do i = 1, n
      new = neginf
      do j = 1, min(i, k)
        if (.not. ieee_is_finite(old(j - 1))) then
          new(j) = old(j) + log(real(i - 1, dp))
        else if (.not. ieee_is_finite(old(j)) .or. i == 1) then
          new(j) = old(j - 1)
        else
          new(j) = logaddexp(old(j - 1), old(j) + log(real(i - 1, dp)))
        end if
      end do
      old = new
    end do
    value = old(k)
  end function log_unsigned_stirling1

  pure real(dp) function logaddexp(a, b) result(value)
    real(dp), intent(in) :: a, b
    real(dp) :: m

    m = max(a, b)
    if (.not. (m > -huge(1.0_dp))) then
      value = m
    else
      value = m + log(exp(a - m) + exp(b - m))
    end if
  end function logaddexp

  integer function number_of_classes(labels) result(k)
    integer, intent(in) :: labels(:)
    integer, allocatable :: seen(:)
    integer :: i, j
    logical :: found

    k = 0
    if (size(labels) == 0) return
    allocate(seen(size(labels)))
    seen = 0
    do i = 1, size(labels)
      found = .false.
      do j = 1, k
        if (labels(i) == seen(j)) then
          found = .true.
          exit
        end if
      end do
      if (.not. found) then
        k = k + 1
        seen(k) = labels(i)
      end if
    end do
  end function number_of_classes

  subroutine class_size_spectrum(labels, mvec, nclass)
    integer, intent(in) :: labels(:)
    integer, allocatable, intent(out) :: mvec(:)
    integer, intent(out), optional :: nclass
    integer, allocatable :: unique_labels(:), counts(:)
    integer :: i, j, k, n, max_count
    logical :: found

    n = size(labels)
    if (n == 0) then
      allocate(mvec(0))
      if (present(nclass)) nclass = 0
      return
    end if

    allocate(unique_labels(n), counts(n))
    unique_labels = 0
    counts = 0
    k = 0
    do i = 1, n
      found = .false.
      do j = 1, k
        if (labels(i) == unique_labels(j)) then
          counts(j) = counts(j) + 1
          found = .true.
          exit
        end if
      end do
      if (.not. found) then
        k = k + 1
        unique_labels(k) = labels(i)
        counts(k) = 1
      end if
    end do

    max_count = maxval(counts(1:k))
    allocate(mvec(max_count))
    mvec = 0
    do i = 1, k
      mvec(counts(i)) = mvec(counts(i)) + 1
    end do
    if (present(nclass)) nclass = k
  end subroutine class_size_spectrum

end module ewens_math
