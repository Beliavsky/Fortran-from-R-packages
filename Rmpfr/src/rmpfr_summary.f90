module rmpfr_summary
  use rmpfr_kinds, only: i64
  use rmpfr_types
  implicit none
  private

  public :: mpfr_sum, mpfr_product, mpfr_minimum, mpfr_maximum, mpfr_range
  public :: mpfr_cumsum, mpfr_cumprod

contains

  function mpfr_sum(x, remove_nan) result(r)
    type(mpfr_real), intent(in) :: x(:) !! Arbitrary-precision values to add.
    logical, intent(in), optional :: remove_nan !! If true, ignore NaNs; otherwise the first NaN makes the result NaN.
    type(mpfr_real) :: r
    logical :: skip_nan
    integer :: i, p

    skip_nan = .false.
    if (present(remove_nan)) skip_nan = remove_nan
    p = max_precision_or_default(x)
    r = mpfr_zero(1, p)
    do i = 1, size(x)
      if (mpfr_is_nan(x(i))) then
        if (skip_nan) cycle
        r = mpfr_nan(p)
        return
      end if
      r = r + mpfr_copy(x(i), p)
    end do
  end function mpfr_sum

  function mpfr_product(x, remove_nan) result(r)
    type(mpfr_real), intent(in) :: x(:) !! Arbitrary-precision values to multiply.
    logical, intent(in), optional :: remove_nan !! If true, ignore NaNs; otherwise the first NaN makes the result NaN.
    type(mpfr_real) :: r
    logical :: skip_nan
    integer :: i, p

    skip_nan = .false.
    if (present(remove_nan)) skip_nan = remove_nan
    p = max_precision_or_default(x)
    r = mpfr_from_integer(1_i64, p)
    do i = 1, size(x)
      if (mpfr_is_nan(x(i))) then
        if (skip_nan) cycle
        r = mpfr_nan(p)
        return
      end if
      r = r * mpfr_copy(x(i), p)
    end do
  end function mpfr_product

  function mpfr_minimum(x, remove_nan) result(r)
    type(mpfr_real), intent(in) :: x(:) !! Arbitrary-precision values whose minimum is required.
    logical, intent(in), optional :: remove_nan !! If true, ignore NaNs; otherwise the first NaN makes the result NaN.
    type(mpfr_real) :: r
    logical :: found, skip_nan
    integer :: i, p

    skip_nan = .false.
    if (present(remove_nan)) skip_nan = remove_nan
    p = max_precision_or_default(x)
    found = .false.
    do i = 1, size(x)
      if (mpfr_is_nan(x(i))) then
        if (skip_nan) cycle
        r = mpfr_nan(p)
        return
      end if
      if (.not. found) then
        r = mpfr_copy(x(i), p)
        found = .true.
      else if (x(i) < r) then
        r = mpfr_copy(x(i), p)
      end if
    end do
    if (.not. found) r = mpfr_inf(1, p)
  end function mpfr_minimum

  function mpfr_maximum(x, remove_nan) result(r)
    type(mpfr_real), intent(in) :: x(:) !! Arbitrary-precision values whose maximum is required.
    logical, intent(in), optional :: remove_nan !! If true, ignore NaNs; otherwise the first NaN makes the result NaN.
    type(mpfr_real) :: r
    logical :: found, skip_nan
    integer :: i, p

    skip_nan = .false.
    if (present(remove_nan)) skip_nan = remove_nan
    p = max_precision_or_default(x)
    found = .false.
    do i = 1, size(x)
      if (mpfr_is_nan(x(i))) then
        if (skip_nan) cycle
        r = mpfr_nan(p)
        return
      end if
      if (.not. found) then
        r = mpfr_copy(x(i), p)
        found = .true.
      else if (x(i) > r) then
        r = mpfr_copy(x(i), p)
      end if
    end do
    if (.not. found) r = mpfr_inf(-1, p)
  end function mpfr_maximum

  subroutine mpfr_range(x, minimum, maximum, remove_nan)
    type(mpfr_real), intent(in) :: x(:) !! Arbitrary-precision values whose range endpoints are required.
    type(mpfr_real), intent(out) :: minimum !! Returned minimum, or +Inf when all values are removed.
    type(mpfr_real), intent(out) :: maximum !! Returned maximum, or -Inf when all values are removed.
    logical, intent(in), optional :: remove_nan !! If true, ignore NaNs; otherwise NaN propagates to both endpoints.

    minimum = mpfr_minimum(x, remove_nan)
    maximum = mpfr_maximum(x, remove_nan)
  end subroutine mpfr_range

  subroutine mpfr_cumsum(x, r)
    type(mpfr_real), intent(in) :: x(:) !! Arbitrary-precision values for cumulative summation.
    type(mpfr_real), allocatable, intent(out) :: r(:) !! Allocated cumulative sums with the same length as x.
    type(mpfr_real) :: running
    integer :: i, p

    p = max_precision_or_default(x)
    allocate(r(size(x)))
    running = mpfr_zero(1, p)
    do i = 1, size(x)
      running = running + mpfr_copy(x(i), p)
      r(i) = running
    end do
  end subroutine mpfr_cumsum

  subroutine mpfr_cumprod(x, r)
    type(mpfr_real), intent(in) :: x(:) !! Arbitrary-precision values for cumulative multiplication.
    type(mpfr_real), allocatable, intent(out) :: r(:) !! Allocated cumulative products with the same length as x.
    type(mpfr_real) :: running
    integer :: i, p

    p = max_precision_or_default(x)
    allocate(r(size(x)))
    running = mpfr_from_integer(1_i64, p)
    do i = 1, size(x)
      running = running * mpfr_copy(x(i), p)
      r(i) = running
    end do
  end subroutine mpfr_cumprod

  integer function max_precision_or_default(x) result(p)
    type(mpfr_real), intent(in) :: x(:) !! Arbitrary-precision array whose maximum element precision is requested.
    integer :: i

    if (size(x) == 0) then
      p = mpfr_get_default_precision()
      return
    end if
    p = 2
    do i = 1, size(x)
      p = max(p, mpfr_precision(x(i)))
    end do
  end function max_precision_or_default

end module rmpfr_summary
