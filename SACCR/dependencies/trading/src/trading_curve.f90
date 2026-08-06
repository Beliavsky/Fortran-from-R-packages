module trading_curve
  use trading_kinds, only : dp
  implicit none
  private

  type, public :: curve_t
    real(dp), allocatable :: tenors(:)
    real(dp), allocatable :: rates(:)
  contains
    procedure :: interpolate => curve_interpolate
  end type curve_t

  public :: natural_cubic_spline
  public :: linear_interpolation

contains

  subroutine curve_interpolate(self, time_points, values, methodology)
    class(curve_t), intent(in) :: self
    real(dp), intent(in) :: time_points(:)
    real(dp), intent(out) :: values(:)
    character(len=*), intent(in), optional :: methodology
    character(len=32) :: method

    if (.not. allocated(self%tenors) .or. .not. allocated(self%rates)) then
      error stop "curve_interpolate: curve is not populated"
    end if
    if (size(self%tenors) /= size(self%rates)) then
      error stop "curve_interpolate: tenors and rates have different sizes"
    end if
    if (size(values) /= size(time_points)) then
      error stop "curve_interpolate: result has the wrong size"
    end if

    method = "cubic_splines"
    if (present(methodology)) method = trim(methodology)

    select case (trim(method))
    case ("cubic_splines", "natural")
      call natural_cubic_spline(self%tenors, self%rates, time_points, values)
    case ("linear")
      call linear_interpolation(self%tenors, self%rates, time_points, values)
    case default
      error stop "curve_interpolate: methodology must be cubic_splines or linear"
    end select
  end subroutine curve_interpolate

  subroutine linear_interpolation(x, y, query, values)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: y(:)
    real(dp), intent(in) :: query(:)
    real(dp), intent(out) :: values(:)
    integer :: i

    if (size(x) < 2 .or. size(x) /= size(y)) then
      error stop "linear_interpolation: invalid knots"
    end if
    if (size(values) /= size(query)) then
      error stop "linear_interpolation: result has the wrong size"
    end if

    do i = 1, size(query)
      values(i) = interpolate_linear_point(x, y, query(i))
    end do
  end subroutine linear_interpolation

  subroutine natural_cubic_spline(x, y, query, values)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: y(:)
    real(dp), intent(in) :: query(:)
    real(dp), intent(out) :: values(:)
    real(dp), allocatable :: second_derivative(:)
    real(dp), allocatable :: u(:)
    real(dp) :: a
    real(dp) :: b
    real(dp) :: h
    real(dp) :: p
    real(dp) :: sigma
    integer :: i
    integer :: k
    integer :: n

    n = size(x)
    if (n < 2 .or. size(y) /= n) then
      error stop "natural_cubic_spline: invalid knots"
    end if
    if (size(values) /= size(query)) then
      error stop "natural_cubic_spline: result has the wrong size"
    end if
    if (any(x(2:) <= x(:n - 1))) then
      error stop "natural_cubic_spline: tenors must be strictly increasing"
    end if

    allocate(second_derivative(n), u(n))
    second_derivative = 0.0_dp
    u = 0.0_dp

    do i = 2, n - 1
      sigma = (x(i) - x(i - 1)) / (x(i + 1) - x(i - 1))
      p = sigma * second_derivative(i - 1) + 2.0_dp
      second_derivative(i) = (sigma - 1.0_dp) / p
      u(i) = (6.0_dp * ((y(i + 1) - y(i)) / (x(i + 1) - x(i)) - &
        (y(i) - y(i - 1)) / (x(i) - x(i - 1))) / &
        (x(i + 1) - x(i - 1)) - sigma * u(i - 1)) / p
    end do

    do k = n - 1, 1, -1
      second_derivative(k) = second_derivative(k) * &
        second_derivative(k + 1) + u(k)
    end do

    do i = 1, size(query)
      if (query(i) <= x(1)) then
        values(i) = interpolate_linear_point(x(:2), y(:2), query(i))
      else if (query(i) >= x(n)) then
        values(i) = interpolate_linear_point(x(n - 1:n), y(n - 1:n), query(i))
      else
        k = locate_interval(x, query(i))
        h = x(k + 1) - x(k)
        a = (x(k + 1) - query(i)) / h
        b = (query(i) - x(k)) / h
        values(i) = a * y(k) + b * y(k + 1) + &
          ((a**3 - a) * second_derivative(k) + &
          (b**3 - b) * second_derivative(k + 1)) * h**2 / 6.0_dp
      end if
    end do
  end subroutine natural_cubic_spline

  pure real(dp) function interpolate_linear_point(x, y, query) result(value)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: y(:)
    real(dp), intent(in) :: query
    integer :: k
    integer :: n

    n = size(x)
    if (query <= x(1)) then
      k = 1
    else if (query >= x(n)) then
      k = n - 1
    else
      k = locate_interval(x, query)
    end if

    value = y(k) + (query - x(k)) * (y(k + 1) - y(k)) / &
      (x(k + 1) - x(k))
  end function interpolate_linear_point

  pure integer function locate_interval(x, query) result(index_value)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: query
    integer :: left
    integer :: middle
    integer :: right

    left = 1
    right = size(x) - 1
    do while (left <= right)
      middle = (left + right) / 2
      if (query < x(middle)) then
        right = middle - 1
      else if (query >= x(middle + 1)) then
        left = middle + 1
      else
        index_value = middle
        return
      end if
    end do
    index_value = max(1, min(size(x) - 1, left))
  end function locate_interval

end module trading_curve
