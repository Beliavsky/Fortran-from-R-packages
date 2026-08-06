! SPDX-License-Identifier: GPL-3.0-only
module mass_multivariate
  use rrcov_kinds, only : dp
  use rrcov_types, only : rrcov_success
  use rrcov_linalg, only : symmetric_eigen
  use mass_types, only : correspondence_result, mca_result, mass_success, &
    mass_invalid_argument, mass_dimension_error, mass_no_convergence
  implicit none
  private
  public :: correspondence_analysis, mca_fit, mca_predict_rows, mca_predict_factors
contains

  subroutine correspondence_analysis(table, dimensions, result)
    real(dp), intent(in) :: table(:, :)
    integer, intent(in), optional :: dimensions
    type(correspondence_result), intent(out) :: result
    real(dp), allocatable :: row_mass(:), column_mass(:), standardized(:, :)
    real(dp), allocatable :: values(:), vectors(:, :), left(:, :)
    real(dp) :: total, singular
    integer :: nr, nc, nf, i, j, k, st

    nr = size(table, 1)
    nc = size(table, 2)
    nf = 1
    if (present(dimensions)) nf = dimensions
    if (nr < 2 .or. nc < 2 .or. any(table < 0.0_dp)) then
      result%status = mass_invalid_argument
      return
    end if
    total = sum(table)
    if (total <= 0.0_dp) then
      result%status = mass_invalid_argument
      return
    end if
    allocate(row_mass(nr), column_mass(nc), standardized(nr, nc))
    row_mass = sum(table, dim=2) / total
    column_mass = sum(table, dim=1) / total
    if (any(row_mass <= 0.0_dp) .or. any(column_mass <= 0.0_dp)) then
      result%status = mass_invalid_argument
      return
    end if
    do i = 1, nr
      do j = 1, nc
        standardized(i, j) = (table(i, j) / total - &
          row_mass(i) * column_mass(j)) / sqrt(row_mass(i) * column_mass(j))
      end do
    end do
    call thin_svd(standardized, values, left, vectors, st)
    nf = min(nf, min(nr - 1, nc - 1))
    allocate(result%correlations(nf), result%row_scores(nr, nf), &
      result%column_scores(nc, nf))
    do k = 1, nf
      singular = values(k)
      result%correlations(k) = singular
      result%row_scores(:, k) = left(:, k) / sqrt(row_mass)
      result%column_scores(:, k) = vectors(:, k) / sqrt(column_mass)
    end do
    result%status = merge(mass_success, mass_no_convergence, st == rrcov_success)
  end subroutine correspondence_analysis

  subroutine mca_fit(codes, levels, dimensions, result)
    integer, intent(in) :: codes(:, :), levels(:)
    integer, intent(in), optional :: dimensions
    type(mca_result), intent(out) :: result
    real(dp), allocatable :: indicator(:, :), counts(:), scaled(:, :)
    real(dp), allocatable :: singular(:), left(:, :), right(:, :)
    integer :: n, factors, columns, nf, i, j, offset, code, k, st

    n = size(codes, 1)
    factors = size(codes, 2)
    nf = 2
    if (present(dimensions)) nf = dimensions
    if (factors /= size(levels) .or. n < 1 .or. any(levels < 2)) then
      result%status = mass_dimension_error
      return
    end if
    columns = sum(levels)
    allocate(indicator(n, columns))
    indicator = 0.0_dp
    offset = 0
    do j = 1, factors
      do i = 1, n
        code = codes(i, j)
        if (code < 1 .or. code > levels(j)) then
          result%status = mass_invalid_argument
          return
        end if
        indicator(i, offset + code) = 1.0_dp
      end do
      offset = offset + levels(j)
    end do
    counts = sum(indicator, dim=1)
    allocate(scaled(n, columns))
    do j = 1, columns
      scaled(:, j) = indicator(:, j) / sqrt(real(factors, dp) * counts(j))
    end do
    call thin_svd(scaled, singular, left, right, st)
    nf = min(nf, min(n - 1, columns - 1))
    allocate(result%singular_values(nf), result%row_scores(n, nf), &
      result%factor_scores(columns, nf))
    do k = 1, nf
      result%singular_values(k) = singular(k + 1)
      result%row_scores(:, k) = matmul(scaled, right(:, k + 1)) / real(factors, dp)
      result%factor_scores(:, k) = right(:, k + 1) / &
        sqrt(real(factors, dp) * counts)
    end do
    result%status = merge(mass_success, mass_no_convergence, st == rrcov_success)
  end subroutine mca_fit

  subroutine mca_predict_rows(codes, levels, model, scores, status)
    integer, intent(in) :: codes(:, :), levels(:)
    type(mca_result), intent(in) :: model
    real(dp), allocatable, intent(out) :: scores(:, :)
    integer, intent(out) :: status
    real(dp), allocatable :: indicator(:, :)
    integer :: n, factors, columns, i, j, offset, code
    n = size(codes, 1)
    factors = size(codes, 2)
    columns = sum(levels)
    if (size(levels) /= factors .or. size(model%factor_scores, 1) /= columns) then
      allocate(scores(0, 0))
      status = mass_dimension_error
      return
    end if
    allocate(indicator(n, columns))
    indicator = 0.0_dp
    offset = 0
    do j = 1, factors
      do i = 1, n
        code = codes(i, j)
        if (code < 1 .or. code > levels(j)) then
          allocate(scores(0, 0))
          status = mass_invalid_argument
          return
        end if
        indicator(i, offset + code) = 1.0_dp
      end do
      offset = offset + levels(j)
    end do
    scores = matmul(indicator, model%factor_scores) / real(factors, dp)
    status = mass_success
  end subroutine mca_predict_rows

  subroutine mca_predict_factors(indicator, row_scores, scores, status)
    real(dp), intent(in) :: indicator(:, :), row_scores(:, :)
    real(dp), allocatable, intent(out) :: scores(:, :)
    integer, intent(out) :: status
    real(dp), allocatable :: counts(:)
    integer :: j
    if (size(indicator, 1) /= size(row_scores, 1)) then
      allocate(scores(0, 0))
      status = mass_dimension_error
      return
    end if
    counts = sum(indicator, dim=1)
    if (any(counts <= 0.0_dp)) then
      allocate(scores(0, 0))
      status = mass_invalid_argument
      return
    end if
    allocate(scores(size(indicator, 2), size(row_scores, 2)))
    scores = matmul(transpose(indicator), row_scores)
    do j = 1, size(indicator, 2)
      scores(j, :) = scores(j, :) / counts(j)
    end do
    status = mass_success
  end subroutine mca_predict_factors

  subroutine thin_svd(a, singular, left, right, status)
    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable, intent(out) :: singular(:), left(:, :), right(:, :)
    integer, intent(out) :: status
    real(dp), allocatable :: gram(:, :), values(:), vectors(:, :)
    integer :: nr, nc, k, i, st
    nr = size(a, 1)
    nc = size(a, 2)
    gram = matmul(transpose(a), a)
    call symmetric_eigen(gram, values, vectors, st)
    k = min(nr, nc)
    allocate(singular(k), right(nc, k), left(nr, k))
    singular = sqrt(max(values(1:k), 0.0_dp))
    right = vectors(:, 1:k)
    left = 0.0_dp
    do i = 1, k
      if (singular(i) > sqrt(epsilon(1.0_dp))) then
        left(:, i) = matmul(a, right(:, i)) / singular(i)
      end if
    end do
    status = st
  end subroutine thin_svd

end module mass_multivariate
