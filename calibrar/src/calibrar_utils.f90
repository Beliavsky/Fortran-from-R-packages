! SPDX-License-Identifier: GPL-2.0-only
module calibrar_utils
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
  use calibrar_kinds, only : dp
  implicit none
  private
  public :: check_bounds, repair_initial, weighted_sum, mean_value, sort_indices, clamp_vector

contains

  subroutine check_bounds(lower, upper, ok)
    real(dp), intent(in) :: lower(:), upper(:)
    logical, intent(out) :: ok
    ok = size(lower) == size(upper)
    if (ok) ok = all(lower < upper)
  end subroutine check_bounds

  subroutine repair_initial(x, lower, upper)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(in) :: lower(:), upper(:)
    integer :: i
    do i = 1, size(x)
      if (ieee_is_nan(x(i))) then
        if (is_finite(lower(i)) .and. is_finite(upper(i))) then
          x(i) = 0.5_dp*(lower(i) + upper(i))
        else if (is_finite(lower(i))) then
          x(i) = lower(i) + 0.1_dp*abs(lower(i))
        else if (is_finite(upper(i))) then
          x(i) = upper(i) - 0.1_dp*abs(upper(i))
        else
          x(i) = 0.0_dp
        end if
      end if
      x(i) = max(lower(i), min(upper(i), x(i)))
    end do
  end subroutine repair_initial

  pure function is_finite(x) result(ok)
    real(dp), intent(in) :: x
    logical :: ok
    ok = ieee_is_finite(x)
  end function is_finite

  pure function weighted_sum(x, w) result(v)
    real(dp), intent(in) :: x(:), w(:)
    real(dp) :: v
    v = sum(x*w)
  end function weighted_sum

  pure function mean_value(x) result(v)
    real(dp), intent(in) :: x(:)
    real(dp) :: v
    if (size(x) == 0) then
      v = 0.0_dp
    else
      v = sum(x)/real(size(x), dp)
    end if
  end function mean_value

  subroutine sort_indices(x, idx)
    real(dp), intent(in) :: x(:)
    integer, intent(out) :: idx(:)
    integer :: i, j, k, tmp
    if (size(idx) /= size(x)) error stop "sort_indices: size mismatch"
    idx = [(i, i=1,size(x))]
    do i = 1, size(x)-1
      k = i
      do j = i+1, size(x)
        if (x(idx(j)) < x(idx(k))) k = j
      end do
      if (k /= i) then
        tmp = idx(i); idx(i) = idx(k); idx(k) = tmp
      end if
    end do
  end subroutine sort_indices

  pure subroutine clamp_vector(x, lower, upper)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(in) :: lower(:), upper(:)
    x = max(lower, min(upper, x))
  end subroutine clamp_vector
end module calibrar_utils
