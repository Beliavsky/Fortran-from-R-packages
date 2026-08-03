! SPDX-License-Identifier: BSD-2-Clause
module fastcluster_distances
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_value, ieee_quiet_nan
  use fastcluster_kinds, only: dp
  use fastcluster_types, only: fc_success, fc_invalid_argument, fc_nan_distance, fc_allocation_failure
  implicit none
  private

  public :: pairwise_distances
  public :: condensed_to_matrix
  public :: matrix_to_condensed

contains

  subroutine pairwise_distances(x, metric, distances, p, status, message, squared_euclidean)
    real(dp), intent(in) :: x(:, :)
    character(len=*), intent(in) :: metric
    real(dp), allocatable, intent(out) :: distances(:, :)
    real(dp), intent(in), optional :: p
    integer, intent(out), optional :: status
    character(len=:), allocatable, intent(out), optional :: message
    logical, intent(in), optional :: squared_euclidean

    character(len=:), allocatable :: metric_name
    real(dp) :: exponent, value
    integer :: i, j, n, stat
    logical :: squared

    call set_status(status, message, fc_success, 'ok')
    n = size(x, 1)
    if (n < 2 .or. size(x, 2) < 1) then
      call set_status(status, message, fc_invalid_argument, &
        'x must contain at least two observations and one variable')
      allocate(distances(0, 0))
      return
    end if

    metric_name = lower_ascii(trim(adjustl(metric)))
    squared = .false.
    if (present(squared_euclidean)) squared = squared_euclidean

    exponent = 2.0_dp
    if (present(p)) exponent = p
    if (metric_name == 'minkowski' .and. exponent <= 0.0_dp) then
      call set_status(status, message, fc_invalid_argument, &
        'the Minkowski exponent p must be positive')
      allocate(distances(0, 0))
      return
    end if
    if (squared .and. metric_name /= 'euclidean') then
      call set_status(status, message, fc_invalid_argument, &
        'squared_euclidean is valid only for the Euclidean metric')
      allocate(distances(0, 0))
      return
    end if

    allocate(distances(n, n), stat=stat)
    if (stat /= 0) then
      call set_status(status, message, fc_allocation_failure, &
        'could not allocate the distance matrix')
      return
    end if
    distances = 0.0_dp

    do j = 1, n - 1
      do i = j + 1, n
        select case (metric_name)
        case ('euclidean')
          value = euclidean_distance(x(i, :), x(j, :), squared)
        case ('maximum')
          value = maximum_distance(x(i, :), x(j, :))
        case ('manhattan')
          value = manhattan_distance(x(i, :), x(j, :))
        case ('canberra')
          value = canberra_distance(x(i, :), x(j, :))
        case ('binary')
          value = binary_distance(x(i, :), x(j, :))
        case ('minkowski')
          value = minkowski_distance(x(i, :), x(j, :), exponent)
        case default
          call set_status(status, message, fc_invalid_argument, &
            'unknown metric: '//metric_name)
          deallocate(distances)
          allocate(distances(0, 0))
          return
        end select
        if (ieee_is_nan(value)) then
          call set_status(status, message, fc_nan_distance, &
            'a pair of observations has no usable coordinates')
          deallocate(distances)
          allocate(distances(0, 0))
          return
        end if
        distances(i, j) = value
        distances(j, i) = value
      end do
    end do
  end subroutine pairwise_distances

  subroutine condensed_to_matrix(d, n, distances, status, message)
    real(dp), intent(in) :: d(:)
    integer, intent(in) :: n
    real(dp), allocatable, intent(out) :: distances(:, :)
    integer, intent(out), optional :: status
    character(len=:), allocatable, intent(out), optional :: message

    integer :: i, j, k, stat

    call set_status(status, message, fc_success, 'ok')
    if (n < 2 .or. size(d) /= n * (n - 1) / 2) then
      call set_status(status, message, fc_invalid_argument, &
        'condensed distances must have n*(n-1)/2 elements')
      allocate(distances(0, 0))
      return
    end if

    allocate(distances(n, n), stat=stat)
    if (stat /= 0) then
      call set_status(status, message, fc_allocation_failure, &
        'could not allocate the distance matrix')
      return
    end if
    distances = 0.0_dp
    k = 0
    do j = 1, n - 1
      do i = j + 1, n
        k = k + 1
        distances(i, j) = d(k)
        distances(j, i) = d(k)
      end do
    end do
  end subroutine condensed_to_matrix

  subroutine matrix_to_condensed(distances, d, status, message)
    real(dp), intent(in) :: distances(:, :)
    real(dp), allocatable, intent(out) :: d(:)
    integer, intent(out), optional :: status
    character(len=:), allocatable, intent(out), optional :: message

    integer :: i, j, k, n

    call set_status(status, message, fc_success, 'ok')
    n = size(distances, 1)
    if (n < 2 .or. size(distances, 2) /= n) then
      call set_status(status, message, fc_invalid_argument, &
        'distances must be a square matrix of order at least two')
      allocate(d(0))
      return
    end if
    allocate(d(n * (n - 1) / 2))
    k = 0
    do j = 1, n - 1
      do i = j + 1, n
        k = k + 1
        d(k) = distances(i, j)
      end do
    end do
  end subroutine matrix_to_condensed

  real(dp) function euclidean_distance(a, b, squared) result(distance)
    real(dp), intent(in) :: a(:), b(:)
    logical, intent(in) :: squared

    real(dp) :: dev, total
    integer :: count, i, n

    n = size(a)
    total = 0.0_dp
    count = 0
    do i = 1, n
      if (.not. ieee_is_nan(a(i)) .and. .not. ieee_is_nan(b(i))) then
        dev = a(i) - b(i)
        if (.not. ieee_is_nan(dev)) then
          total = total + dev * dev
          count = count + 1
        end if
      end if
    end do
    if (count == 0) then
      distance = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    if (count /= n) total = total * real(n, dp) / real(count, dp)
    if (squared) then
      distance = total
    else
      distance = sqrt(total)
    end if
  end function euclidean_distance

  real(dp) function maximum_distance(a, b) result(distance)
    real(dp), intent(in) :: a(:), b(:)

    real(dp) :: dev
    integer :: count, i

    distance = -huge(1.0_dp)
    count = 0
    do i = 1, size(a)
      if (.not. ieee_is_nan(a(i)) .and. .not. ieee_is_nan(b(i))) then
        dev = abs(a(i) - b(i))
        if (.not. ieee_is_nan(dev)) then
          distance = max(distance, dev)
          count = count + 1
        end if
      end if
    end do
    if (count == 0) distance = ieee_value(0.0_dp, ieee_quiet_nan)
  end function maximum_distance

  real(dp) function manhattan_distance(a, b) result(distance)
    real(dp), intent(in) :: a(:), b(:)

    real(dp) :: dev
    integer :: count, i, n

    n = size(a)
    distance = 0.0_dp
    count = 0
    do i = 1, n
      if (.not. ieee_is_nan(a(i)) .and. .not. ieee_is_nan(b(i))) then
        dev = abs(a(i) - b(i))
        if (.not. ieee_is_nan(dev)) then
          distance = distance + dev
          count = count + 1
        end if
      end if
    end do
    if (count == 0) then
      distance = ieee_value(0.0_dp, ieee_quiet_nan)
    else if (count /= n) then
      distance = distance * real(n, dp) / real(count, dp)
    end if
  end function manhattan_distance

  real(dp) function canberra_distance(a, b) result(distance)
    real(dp), intent(in) :: a(:), b(:)

    real(dp) :: denominator, dev, difference
    integer :: count, i, n

    n = size(a)
    distance = 0.0_dp
    count = 0
    do i = 1, n
      if (.not. ieee_is_nan(a(i)) .and. .not. ieee_is_nan(b(i))) then
        denominator = abs(a(i)) + abs(b(i))
        difference = abs(a(i) - b(i))
        if (denominator > tiny(1.0_dp) .or. difference > tiny(1.0_dp)) then
          dev = difference / denominator
          if (ieee_is_nan(dev)) then
            if (.not. ieee_is_finite(difference) .and. .not. ieee_is_finite(denominator)) dev = 1.0_dp
          end if
          if (.not. ieee_is_nan(dev)) then
            distance = distance + dev
            count = count + 1
          end if
        end if
      end if
    end do
    if (count == 0) then
      distance = ieee_value(0.0_dp, ieee_quiet_nan)
    else if (count /= n) then
      distance = distance * real(n, dp) / real(count, dp)
    end if
  end function canberra_distance

  real(dp) function binary_distance(a, b) result(distance)
    real(dp), intent(in) :: a(:), b(:)

    integer :: count, different, i, total
    logical :: a_nonzero, b_nonzero

    total = 0
    count = 0
    different = 0
    do i = 1, size(a)
      if (.not. ieee_is_nan(a(i)) .and. .not. ieee_is_nan(b(i))) then
        if (ieee_is_finite(a(i)) .and. ieee_is_finite(b(i))) then
          total = total + 1
          a_nonzero = a(i) > 0.0_dp .or. a(i) < 0.0_dp
          b_nonzero = b(i) > 0.0_dp .or. b(i) < 0.0_dp
          if (a_nonzero .or. b_nonzero) then
            count = count + 1
            if (a_nonzero .neqv. b_nonzero) different = different + 1
          end if
        end if
      end if
    end do
    if (total == 0) then
      distance = ieee_value(0.0_dp, ieee_quiet_nan)
    else if (count == 0) then
      distance = 0.0_dp
    else
      distance = real(different, dp) / real(count, dp)
    end if
  end function binary_distance

  real(dp) function minkowski_distance(a, b, p) result(distance)
    real(dp), intent(in) :: a(:), b(:), p

    real(dp) :: dev
    integer :: count, i, n

    n = size(a)
    distance = 0.0_dp
    count = 0
    do i = 1, n
      if (.not. ieee_is_nan(a(i)) .and. .not. ieee_is_nan(b(i))) then
        dev = a(i) - b(i)
        if (.not. ieee_is_nan(dev)) then
          distance = distance + abs(dev) ** p
          count = count + 1
        end if
      end if
    end do
    if (count == 0) then
      distance = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    if (count /= n) distance = distance * real(n, dp) / real(count, dp)
    distance = distance ** (1.0_dp / p)
  end function minkowski_distance

  pure function lower_ascii(text) result(lower)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: lower

    integer :: code, i

    do i = 1, len(text)
      code = iachar(text(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) then
        lower(i:i) = achar(code + iachar('a') - iachar('A'))
      else
        lower(i:i) = text(i:i)
      end if
    end do
  end function lower_ascii

  subroutine set_status(status, message, code, text)
    integer, intent(out), optional :: status
    character(len=:), allocatable, intent(out), optional :: message
    integer, intent(in) :: code
    character(len=*), intent(in) :: text

    if (present(status)) status = code
    if (present(message)) message = text
  end subroutine set_status

end module fastcluster_distances
