! SPDX-License-Identifier: GPL-3.0-or-later
module corpcor_matrix_tools
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use corpcor_kinds, only : dp
  use corpcor_types, only : covariance_decomposition, precision_decomposition, &
    corpcor_success, corpcor_invalid_argument, corpcor_dimension_error
  use corpcor_linalg, only : covariance_to_correlation
  implicit none
  private
  public :: symmetric_matrix_to_vector, symmetric_matrix_indices
  public :: vector_to_symmetric_matrix, rebuild_covariance, decompose_covariance
  public :: rebuild_precision, decompose_precision

contains

  function symmetric_matrix_to_vector(a, include_diagonal, status) result(v)
    real(dp), intent(in) :: a(:, :)
    logical, intent(in), optional :: include_diagonal
    integer, intent(out), optional :: status
    real(dp), allocatable :: v(:)
    logical :: diag
    integer :: n, i, j, k, nv, istat

    n = size(a, 1)
    diag = .false.
    if (present(include_diagonal)) diag = include_diagonal
    istat = corpcor_success
    if (size(a, 2) /= n) then
      allocate(v(0))
      istat = corpcor_dimension_error
    else
      if (diag) then
        nv = n * (n + 1) / 2
      else
        nv = n * (n - 1) / 2
      end if
      allocate(v(nv))
      k = 0
      do j = 1, n
        if (diag) then
          do i = j, n
            k = k + 1
            v(k) = a(i, j)
          end do
        else
          do i = j + 1, n
            k = k + 1
            v(k) = a(i, j)
          end do
        end if
      end do
    end if
    if (present(status)) status = istat
  end function symmetric_matrix_to_vector

  function symmetric_matrix_indices(n, include_diagonal, status) result(index)
    integer, intent(in) :: n
    logical, intent(in), optional :: include_diagonal
    integer, intent(out), optional :: status
    integer, allocatable :: index(:, :)
    logical :: diag
    integer :: i, j, k, nv, istat

    diag = .false.
    if (present(include_diagonal)) diag = include_diagonal
    istat = corpcor_success
    if (n < 1) then
      allocate(index(0, 2))
      istat = corpcor_invalid_argument
    else
      if (diag) then
        nv = n * (n + 1) / 2
      else
        nv = n * (n - 1) / 2
      end if
      allocate(index(nv, 2))
      k = 0
      do i = 1, n
        if (diag) then
          do j = i, n
            k = k + 1
            index(k, :) = [i, j]
          end do
        else
          do j = i + 1, n
            k = k + 1
            index(k, :) = [i, j]
          end do
        end if
      end do
    end if
    if (present(status)) status = istat
  end function symmetric_matrix_indices

  function vector_to_symmetric_matrix(v, include_diagonal, order, status) result(a)
    real(dp), intent(in) :: v(:)
    logical, intent(in), optional :: include_diagonal
    integer, intent(in), optional :: order(:)
    integer, intent(out), optional :: status
    real(dp), allocatable :: a(:, :)
    real(dp), allocatable :: vin(:)
    real(dp) :: nr, nan_value
    logical :: diag
    integer :: n, i, j, k, istat

    diag = .false.
    if (present(include_diagonal)) diag = include_diagonal
    istat = corpcor_success
    if (diag) then
      nr = 0.5_dp * (-1.0_dp + sqrt(1.0_dp + 8.0_dp * real(size(v), dp)))
    else
      nr = 0.5_dp * (1.0_dp + sqrt(1.0_dp + 8.0_dp * real(size(v), dp)))
    end if
    n = nint(nr)
    if (n < 1 .or. abs(nr - real(n, dp)) > 100.0_dp * epsilon(1.0_dp)) then
      allocate(a(0, 0))
      istat = corpcor_invalid_argument
      if (present(status)) status = istat
      return
    end if
    allocate(vin(size(v)))
    if (present(order)) then
      if (size(order) /= size(v) .or. any(order < 1) .or. any(order > size(v))) then
        allocate(a(0, 0))
        istat = corpcor_invalid_argument
        if (present(status)) status = istat
        return
      end if
      vin = 0.0_dp
      do i = 1, size(v)
        vin(order(i)) = v(i)
      end do
    else
      vin = v
    end if

    allocate(a(n, n))
    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    a = nan_value
    k = 0
    do j = 1, n
      if (diag) then
        do i = j, n
          k = k + 1
          a(i, j) = vin(k)
          a(j, i) = vin(k)
        end do
      else
        do i = j + 1, n
          k = k + 1
          a(i, j) = vin(k)
          a(j, i) = vin(k)
        end do
      end if
    end do
    if (present(status)) status = istat
  end function vector_to_symmetric_matrix

  function rebuild_covariance(correlation, variance, status) result(covariance)
    real(dp), intent(in) :: correlation(:, :)
    real(dp), intent(in) :: variance(:)
    integer, intent(out), optional :: status
    real(dp), allocatable :: covariance(:, :)
    real(dp), allocatable :: sd(:)
    integer :: n, i, j, istat

    n = size(correlation, 1)
    allocate(covariance(n, size(correlation, 2)))
    covariance = 0.0_dp
    istat = corpcor_success
    if (size(correlation, 2) /= n .or. size(variance) /= n) then
      istat = corpcor_dimension_error
    else if (any(variance < 0.0_dp)) then
      istat = corpcor_invalid_argument
    else
      sd = sqrt(variance)
      do j = 1, n
        do i = 1, n
          covariance(i, j) = correlation(i, j) * sd(i) * sd(j)
        end do
      end do
    end if
    if (present(status)) status = istat
  end function rebuild_covariance

  function decompose_covariance(covariance) result(res)
    real(dp), intent(in) :: covariance(:, :)
    type(covariance_decomposition) :: res
    integer :: n, i

    n = size(covariance, 1)
    allocate(res%variance(n))
    if (size(covariance, 2) /= n) then
      allocate(res%correlation(0, 0))
      res%variance = 0.0_dp
      res%status = corpcor_dimension_error
      return
    end if
    do i = 1, n
      res%variance(i) = covariance(i, i)
    end do
    res%correlation = covariance_to_correlation(covariance, res%status)
  end function decompose_covariance

  function rebuild_precision(partial_correlation, partial_variance, status) result(precision)
    real(dp), intent(in) :: partial_correlation(:, :)
    real(dp), intent(in) :: partial_variance(:)
    integer, intent(out), optional :: status
    real(dp), allocatable :: precision(:, :)
    real(dp), allocatable :: inverse_sd(:)
    integer :: n, i, j, istat

    n = size(partial_correlation, 1)
    allocate(precision(n, size(partial_correlation, 2)))
    precision = 0.0_dp
    istat = corpcor_success
    if (size(partial_correlation, 2) /= n .or. size(partial_variance) /= n) then
      istat = corpcor_dimension_error
    else if (any(partial_variance <= 0.0_dp)) then
      istat = corpcor_invalid_argument
    else
      inverse_sd = sqrt(1.0_dp / partial_variance)
      do j = 1, n
        do i = 1, n
          precision(i, j) = -partial_correlation(i, j) * inverse_sd(i) * inverse_sd(j)
        end do
        precision(j, j) = -precision(j, j)
      end do
    end if
    if (present(status)) status = istat
  end function rebuild_precision

  function decompose_precision(precision) result(res)
    real(dp), intent(in) :: precision(:, :)
    type(precision_decomposition) :: res
    real(dp), allocatable :: work(:, :)
    integer :: n, i

    n = size(precision, 1)
    allocate(res%partial_variance(n))
    if (size(precision, 2) /= n .or. any([(precision(i, i) <= 0.0_dp, i=1,n)])) then
      allocate(res%partial_correlation(0, 0))
      res%partial_variance = 0.0_dp
      res%status = corpcor_invalid_argument
      return
    end if
    do i = 1, n
      res%partial_variance(i) = 1.0_dp / precision(i, i)
    end do
    work = -precision
    do i = 1, n
      work(i, i) = -work(i, i)
    end do
    res%partial_correlation = covariance_to_correlation(work, res%status)
  end function decompose_precision

end module corpcor_matrix_tools
