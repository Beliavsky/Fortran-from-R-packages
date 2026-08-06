module polynom_orthogonal
  use polynom_kinds, only : dp
  use polynom_status, only : poly_status_t, poly_ok, poly_invalid_argument, set_status
  use polynom_core, only : polynomial_t, polylist_t, polynomial, operator(-), &
    operator(*), operator(/)
  implicit none
  private
  public :: orthogonal_polynomials

contains

  function orthogonal_polynomials(x, degree, normalized, status) result(list)
    real(dp), intent(in) :: x(:)
    integer, intent(in), optional :: degree
    logical, intent(in), optional :: normalized
    type(poly_status_t), intent(out), optional :: status
    type(polylist_t) :: list
    type(polynomial_t), allocatable :: monic(:)
    real(dp), allocatable :: values(:,:), norms(:), alpha(:)
    type(polynomial_t) :: xpoly
    logical :: norm_flag
    integer :: d, j, unique_count
    real(dp) :: beta

    unique_count = count_unique(x)
    d = max(0, unique_count - 1)
    if (present(degree)) d = degree
    norm_flag = .true.
    if (present(normalized)) norm_flag = normalized
    if (size(x) == 0 .or. d < 0 .or. d >= unique_count) then
      allocate(list%item(0))
      call set_status(status, poly_invalid_argument, &
        'degree must be between zero and the number of distinct x values minus one')
      return
    end if

    allocate(monic(0:d), values(size(x), 0:d), norms(0:d), alpha(0:max(0, d - 1)))
    xpoly = polynomial([0.0_dp, 1.0_dp])
    monic(0) = polynomial(1.0_dp)
    values(:, 0) = 1.0_dp
    norms(0) = real(size(x), dp)
    if (d >= 1) then
      alpha(0) = sum(x * values(:, 0)**2) / norms(0)
      monic(1) = xpoly - alpha(0)
      values(:, 1) = monic(1)%evaluate(x)
      norms(1) = sum(values(:, 1)**2)
    end if
    do j = 1, d - 1
      alpha(j) = sum(x * values(:, j)**2) / norms(j)
      beta = norms(j) / norms(j - 1)
      monic(j + 1) = (xpoly - alpha(j)) * monic(j) - beta * monic(j - 1)
      values(:, j + 1) = monic(j + 1)%evaluate(x)
      norms(j + 1) = sum(values(:, j + 1)**2)
    end do

    allocate(list%item(d + 1))
    do j = 0, d
      if (norm_flag) then
        list%item(j + 1) = monic(j) / sqrt(norms(j))
      else
        list%item(j + 1) = monic(j)
      end if
    end do
    call set_status(status, poly_ok, '')
  end function orthogonal_polynomials

  integer function count_unique(x)
    real(dp), intent(in) :: x(:)
    logical, allocatable :: seen(:)
    integer :: i, j
    allocate(seen(size(x)), source=.false.)
    count_unique = 0
    do i = 1, size(x)
      if (seen(i)) cycle
      count_unique = count_unique + 1
      do j = i, size(x)
        if (.not. (x(j) < x(i) .or. x(j) > x(i))) seen(j) = .true.
      end do
    end do
  end function count_unique

end module polynom_orthogonal
