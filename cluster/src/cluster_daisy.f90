! SPDX-License-Identifier: GPL-2.0-or-later
module cluster_daisy
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
  use fastcluster_kinds, only: dp
  use fastcluster_distances, only: pairwise_distances
  use cluster_types, only: cluster_success, cluster_invalid_argument
  implicit none
  private

  integer, parameter, public :: variable_numeric = 0
  integer, parameter, public :: variable_binary_symmetric = 1
  integer, parameter, public :: variable_binary_asymmetric = 2
  integer, parameter, public :: variable_nominal = 3
  integer, parameter, public :: variable_ordinal = 4

  public :: daisy
  public :: daisy_mixed

contains

  subroutine daisy(x, distances, metric, weights, status, message)
    real(dp), intent(in) :: x(:, :)
    real(dp), allocatable, intent(out) :: distances(:, :)
    character(len=*), intent(in), optional :: metric
    real(dp), intent(in), optional :: weights(:)
    integer, intent(out), optional :: status
    character(len=:), allocatable, intent(out), optional :: message

    character(len=:), allocatable :: met
    real(dp), allocatable :: work(:, :)
    integer :: i, j, p, n, st
    real(dp) :: num, den, scale

    met = 'euclidean'
    if (present(metric)) met = lower_ascii(trim(adjustl(metric)))
    if (size(x, 1) < 2 .or. size(x, 2) < 1) then
      call set_status(status, message, cluster_invalid_argument, &
        'x must contain at least two observations and one variable')
      allocate(distances(0, 0))
      return
    end if
    if (present(weights)) then
      if (size(weights) /= size(x, 2) .or. any(weights < 0.0_dp) .or. sum(weights) <= 0.0_dp) then
        call set_status(status, message, cluster_invalid_argument, &
          'weights must be nonnegative and match the variables')
        allocate(distances(0, 0))
        return
      end if
    end if
    if (met /= 'gower' .and. .not. present(weights)) then
      call pairwise_distances(x, met, distances, status=st, message=message)
      if (present(status)) status = st
      return
    end if

    n = size(x, 1)
    p = size(x, 2)
    allocate(work(n, p), distances(n, n))
    work = x
    if (met == 'gower') then
      do j = 1, p
        call finite_range(x(:, j), scale)
        if (scale > 0.0_dp) then
          do i = 1, n
            if (.not. ieee_is_nan(work(i, j))) work(i, j) = work(i, j) / scale
          end do
        else
          work(:, j) = 0.0_dp
        end if
      end do
    else if (met /= 'euclidean' .and. met /= 'manhattan') then
      call set_status(status, message, cluster_invalid_argument, &
        'weighted daisy supports euclidean, manhattan, or gower metrics')
      deallocate(distances)
      allocate(distances(0, 0))
      return
    end if

    distances = 0.0_dp
    do i = 1, n - 1
      do j = i + 1, n
        num = 0.0_dp
        den = 0.0_dp
        call weighted_pair(work(i, :), work(j, :), met, weights, num, den)
        if (den <= 0.0_dp) then
          distances(i, j) = ieee_value(0.0_dp, ieee_quiet_nan)
        else if (met == 'euclidean') then
          distances(i, j) = sqrt(num * total_weight(weights, p) / den)
        else if (met == 'manhattan') then
          distances(i, j) = num * total_weight(weights, p) / den
        else
          distances(i, j) = num / den
        end if
        distances(j, i) = distances(i, j)
      end do
    end do
    call set_status(status, message, cluster_success, 'ok')
  end subroutine daisy

  subroutine daisy_mixed(x, variable_types, distances, weights, status, message)
    real(dp), intent(in) :: x(:, :)
    integer, intent(in) :: variable_types(:)
    real(dp), allocatable, intent(out) :: distances(:, :)
    real(dp), intent(in), optional :: weights(:)
    integer, intent(out), optional :: status
    character(len=:), allocatable, intent(out), optional :: message

    real(dp), allocatable :: ranges(:), w(:)
    real(dp) :: contribution, denominator, numerator
    integer :: i, j, l, n, p
    logical :: usable

    n = size(x, 1)
    p = size(x, 2)
    if (n < 2 .or. p < 1 .or. size(variable_types) /= p) then
      call set_status(status, message, cluster_invalid_argument, &
        'variable_types must contain one entry per variable')
      allocate(distances(0, 0))
      return
    end if
    if (any(variable_types < variable_numeric) .or. any(variable_types > variable_ordinal)) then
      call set_status(status, message, cluster_invalid_argument, 'unknown variable type')
      allocate(distances(0, 0))
      return
    end if
    allocate(w(p), ranges(p), distances(n, n))
    w = 1.0_dp
    if (present(weights)) then
      if (size(weights) /= p .or. any(weights < 0.0_dp) .or. sum(weights) <= 0.0_dp) then
        call set_status(status, message, cluster_invalid_argument, 'invalid variable weights')
        deallocate(distances)
        allocate(distances(0, 0))
        return
      end if
      w = weights
    end if
    ranges = 1.0_dp
    do l = 1, p
      if (variable_types(l) == variable_numeric .or. variable_types(l) == variable_ordinal) then
        call finite_range(x(:, l), ranges(l))
        if (ranges(l) <= 0.0_dp) ranges(l) = 1.0_dp
      end if
    end do
    distances = 0.0_dp
    do i = 1, n - 1
      do j = i + 1, n
        numerator = 0.0_dp
        denominator = 0.0_dp
        do l = 1, p
          usable = .not. ieee_is_nan(x(i, l)) .and. .not. ieee_is_nan(x(j, l))
          if (.not. usable) cycle
          contribution = 0.0_dp
          select case (variable_types(l))
          case (variable_numeric, variable_ordinal)
            contribution = min(1.0_dp, abs(x(i, l) - x(j, l)) / ranges(l))
          case (variable_binary_symmetric, variable_nominal)
            if (nint(x(i, l)) == nint(x(j, l))) then
              contribution = 0.0_dp
            else
              contribution = 1.0_dp
            end if
          case (variable_binary_asymmetric)
            if (abs(x(i, l)) <= epsilon(1.0_dp) .and. &
                abs(x(j, l)) <= epsilon(1.0_dp)) cycle
            if (nint(x(i, l)) == nint(x(j, l))) then
              contribution = 0.0_dp
            else
              contribution = 1.0_dp
            end if
          end select
          numerator = numerator + w(l) * contribution
          denominator = denominator + w(l)
        end do
        if (denominator > 0.0_dp) then
          distances(i, j) = numerator / denominator
        else
          distances(i, j) = ieee_value(0.0_dp, ieee_quiet_nan)
        end if
        distances(j, i) = distances(i, j)
      end do
    end do
    call set_status(status, message, cluster_success, 'ok')
  end subroutine daisy_mixed

  subroutine weighted_pair(a, b, metric, weights, numerator, denominator)
    real(dp), intent(in) :: a(:), b(:)
    character(len=*), intent(in) :: metric
    real(dp), intent(in), optional :: weights(:)
    real(dp), intent(out) :: numerator, denominator

    real(dp) :: wt, dev
    integer :: l

    numerator = 0.0_dp
    denominator = 0.0_dp
    do l = 1, size(a)
      if (ieee_is_nan(a(l)) .or. ieee_is_nan(b(l))) cycle
      wt = 1.0_dp
      if (present(weights)) wt = weights(l)
      dev = abs(a(l) - b(l))
      if (metric == 'euclidean') dev = dev * dev
      numerator = numerator + wt * dev
      denominator = denominator + wt
    end do
  end subroutine weighted_pair

  real(dp) function total_weight(weights, p) result(value)
    real(dp), intent(in), optional :: weights(:)
    integer, intent(in) :: p
    if (present(weights)) then
      value = sum(weights)
    else
      value = real(p, dp)
    end if
  end function total_weight

  subroutine finite_range(x, range_value)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: range_value

    real(dp) :: xmin, xmax
    integer :: i, count

    count = 0
    xmin = huge(1.0_dp)
    xmax = -huge(1.0_dp)
    do i = 1, size(x)
      if (.not. ieee_is_nan(x(i))) then
        xmin = min(xmin, x(i))
        xmax = max(xmax, x(i))
        count = count + 1
      end if
    end do
    if (count == 0) then
      range_value = 0.0_dp
    else
      range_value = xmax - xmin
    end if
  end subroutine finite_range

  subroutine set_status(status, message, code, text)
    integer, intent(out), optional :: status
    character(len=:), allocatable, intent(out), optional :: message
    integer, intent(in) :: code
    character(len=*), intent(in) :: text
    if (present(status)) status = code
    if (present(message)) message = text
  end subroutine set_status

  pure function lower_ascii(text) result(lower)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: lower
    integer :: i, code
    lower = text
    do i = 1, len(text)
      code = iachar(text(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code + 32)
    end do
  end function lower_ascii

end module cluster_daisy
