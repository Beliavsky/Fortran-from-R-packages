! SPDX-License-Identifier: GPL-3.0-or-later
module corpcor_weighted
  use corpcor_kinds, only : dp
  use corpcor_types, only : moments_result, scale_result, corpcor_success, &
    corpcor_invalid_argument, corpcor_dimension_error
  implicit none
  private
  public :: normalized_weights, weighted_variance, weighted_moments, weighted_scale
  public :: median_value

contains

  function normalized_weights(n, w, status) result(wn)
    integer, intent(in) :: n
    real(dp), intent(in), optional :: w(:)
    integer, intent(out), optional :: status
    real(dp), allocatable :: wn(:)
    real(dp) :: sw
    integer :: istat

    allocate(wn(max(0, n)))
    istat = corpcor_success
    if (n < 1) then
      istat = corpcor_dimension_error
    else if (present(w)) then
      if (size(w) /= n .or. any(w < 0.0_dp)) then
        wn = 0.0_dp
        istat = corpcor_invalid_argument
      else
        sw = sum(w)
        if (sw <= 0.0_dp) then
          wn = 0.0_dp
          istat = corpcor_invalid_argument
        else
          wn = w / sw
        end if
      end if
    else
      wn = 1.0_dp / real(n, dp)
    end if
    if (present(status)) status = istat
  end function normalized_weights

  function weighted_variance(x, w, status) result(v)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: w(:)
    integer, intent(out), optional :: status
    real(dp) :: v
    real(dp), allocatable :: wn(:)
    real(dp) :: mu, denom
    integer :: istat

    wn = normalized_weights(size(x), w, istat)
    v = 0.0_dp
    if (istat == corpcor_success) then
      denom = 1.0_dp - sum(wn * wn)
      if (denom <= epsilon(1.0_dp)) then
        istat = corpcor_invalid_argument
      else
        mu = sum(wn * x)
        v = sum(wn * (x - mu) ** 2) / denom
        if (v < epsilon(1.0_dp)) v = 0.0_dp
      end if
    end if
    if (present(status)) status = istat
  end function weighted_variance

  function weighted_moments(x, w) result(res)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(in), optional :: w(:)
    type(moments_result) :: res
    real(dp), allocatable :: wn(:)
    real(dp) :: denom
    integer :: n, p, j, istat

    n = size(x, 1)
    p = size(x, 2)
    allocate(res%mean(p), res%variance(p))
    res%mean = 0.0_dp
    res%variance = 0.0_dp
    wn = normalized_weights(n, w, istat)
    if (istat /= corpcor_success) then
      res%status = istat
      return
    end if
    denom = 1.0_dp - sum(wn * wn)
    if (denom <= epsilon(1.0_dp)) then
      res%status = corpcor_invalid_argument
      return
    end if
    do j = 1, p
      res%mean(j) = sum(wn * x(:, j))
      res%variance(j) = sum(wn * (x(:, j) - res%mean(j)) ** 2) / denom
      if (res%variance(j) < epsilon(1.0_dp)) res%variance(j) = 0.0_dp
    end do
    res%status = corpcor_success
  end function weighted_moments

  function weighted_scale(x, w, center, scale) result(res)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(in), optional :: w(:)
    logical, intent(in), optional :: center, scale
    type(scale_result) :: res
    type(moments_result) :: mom
    logical :: do_center, do_scale
    integer :: n, p, j

    n = size(x, 1)
    p = size(x, 2)
    allocate(res%x(n, p), res%center(p), res%scale(p), res%zero_scale(p))
    res%x = x
    res%center = 0.0_dp
    res%scale = 1.0_dp
    res%zero_scale = .false.
    do_center = .true.
    do_scale = .true.
    if (present(center)) do_center = center
    if (present(scale)) do_scale = scale

    mom = weighted_moments(x, w)
    res%status = mom%status
    if (res%status /= corpcor_success) return
    res%center = mom%mean
    res%scale = sqrt(max(0.0_dp, mom%variance))
    res%zero_scale = res%scale <= epsilon(1.0_dp)

    if (do_center) then
      do j = 1, p
        res%x(:, j) = res%x(:, j) - res%center(j)
      end do
    end if
    if (do_scale) then
      do j = 1, p
        if (res%scale(j) > 0.0_dp) then
          res%x(:, j) = res%x(:, j) / res%scale(j)
        else
          res%x(:, j) = 0.0_dp
        end if
      end do
    end if
  end function weighted_scale

  function median_value(x) result(med)
    real(dp), intent(in) :: x(:)
    real(dp) :: med
    real(dp), allocatable :: y(:)
    real(dp) :: tmp
    integer :: i, j, n

    n = size(x)
    if (n == 0) then
      med = 0.0_dp
      return
    end if
    y = x
    do i = 2, n
      tmp = y(i)
      j = i - 1
      do while (j >= 1)
        if (y(j) <= tmp) exit
        y(j + 1) = y(j)
        j = j - 1
      end do
      y(j + 1) = tmp
    end do
    if (mod(n, 2) == 1) then
      med = y((n + 1) / 2)
    else
      med = 0.5_dp * (y(n / 2) + y(n / 2 + 1))
    end if
  end function median_value

end module corpcor_weighted
