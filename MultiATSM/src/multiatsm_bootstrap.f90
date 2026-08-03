! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module multiatsm_bootstrap
  use multiatsm_kinds, only : dp
  use multiatsm_random, only : set_random_seed, random_integer
  use multiatsm_types, only : var_model, bootstrap_result
  use multiatsm_var, only : fit_var
  implicit none
  private

  integer, parameter, public :: BOOTSTRAP_IID = 1
  integer, parameter, public :: BOOTSTRAP_WILD = 2
  integer, parameter, public :: BOOTSTRAP_BLOCK = 3

  public :: resample_residuals, simulate_var, bootstrap_var
  public :: percentile_bounds

contains

  subroutine resample_residuals(factor_residuals, yield_residuals, method, factor_draws, &
      yield_draws, info, block_length, seed)
    real(dp), intent(in) :: factor_residuals(:, :), yield_residuals(:, :)
    integer, intent(in) :: method
    real(dp), allocatable, intent(out) :: factor_draws(:, :), yield_draws(:, :)
    integer, intent(out) :: info
    integer, intent(in), optional :: block_length, seed
    integer, allocatable :: index(:)
    integer :: tf, ty, nmax, i, start, j, block, sign_value
    real(dp) :: u

    tf = size(factor_residuals, 2)
    ty = size(yield_residuals, 2)
    if (tf < 1 .or. ty < 1) then
      allocate(factor_draws(0, 0), yield_draws(0, 0))
      info = -1
      return
    end if
    if (present(seed)) call set_random_seed(seed)
    allocate(factor_draws(size(factor_residuals, 1), tf))
    allocate(yield_draws(size(yield_residuals, 1), ty))
    nmax = max(tf, ty)
    allocate(index(nmax))

    select case (method)
    case (BOOTSTRAP_IID)
      do i = 1, nmax
        index(i) = random_integer(min(tf, ty))
      end do
      do i = 1, tf
        factor_draws(:, i) = factor_residuals(:, 1 + modulo(index(i) - 1, tf))
      end do
      do i = 1, ty
        yield_draws(:, i) = yield_residuals(:, 1 + modulo(index(i) - 1, ty))
      end do
    case (BOOTSTRAP_WILD)
      do i = 1, nmax
        call random_number(u)
        sign_value = merge(1, -1, u <= 0.5_dp)
        if (i <= tf) factor_draws(:, i) = real(sign_value, dp) * factor_residuals(:, i)
        if (i <= ty) yield_draws(:, i) = real(sign_value, dp) * yield_residuals(:, i)
      end do
    case (BOOTSTRAP_BLOCK)
      block = 1
      if (present(block_length)) block = max(1, block_length)
      i = 1
      do while (i <= nmax)
        start = random_integer(max(1, min(tf, ty) - block + 1))
        do j = 0, block - 1
          if (i + j > nmax) exit
          index(i + j) = start + j
        end do
        i = i + block
      end do
      do i = 1, tf
        factor_draws(:, i) = factor_residuals(:, 1 + modulo(index(i) - 1, tf))
      end do
      do i = 1, ty
        yield_draws(:, i) = yield_residuals(:, 1 + modulo(index(i) - 1, ty))
      end do
    case default
      factor_draws = 0.0_dp
      yield_draws = 0.0_dp
      info = -2
      return
    end select
    info = 0
  end subroutine resample_residuals

  subroutine simulate_var(intercept, phi, residuals, initial, series, info)
    real(dp), intent(in) :: intercept(:), phi(:, :), residuals(:, :), initial(:)
    real(dp), allocatable, intent(out) :: series(:, :)
    integer, intent(out) :: info
    integer :: k, t, j

    k = size(intercept)
    t = size(residuals, 2) + 1
    if (size(phi, 1) /= k .or. size(phi, 2) /= k .or. &
        size(residuals, 1) /= k .or. size(initial) /= k) then
      allocate(series(0, 0))
      info = -1
      return
    end if
    allocate(series(k, t))
    series(:, 1) = initial
    do j = 2, t
      series(:, j) = intercept + matmul(phi, series(:, j - 1)) + residuals(:, j - 1)
    end do
    info = 0
  end subroutine simulate_var

  subroutine bootstrap_var(data, n_draws, method, result, info, block_length, seed)
    real(dp), intent(in) :: data(:, :)
    integer, intent(in) :: n_draws, method
    type(bootstrap_result), intent(out) :: result
    integer, intent(out) :: info
    integer, intent(in), optional :: block_length, seed
    type(var_model) :: base, fitted
    real(dp), allocatable :: residual_draw(:, :), dummy_y(:, :), dummy_draw(:, :), simulated(:, :)
    real(dp), allocatable :: matrix_draws(:, :, :)
    integer :: b, local_seed, k, t

    if (n_draws < 1 .or. size(data, 2) < 3) then
      info = -1
      return
    end if
    call fit_var(data, base, info)
    if (info /= 0) return
    k = size(data, 1)
    t = size(data, 2)
    allocate(result%phi_draws(k, k, n_draws), result%sigma_draws(k, k, n_draws))
    allocate(dummy_y(1, t - 1))
    dummy_y = 0.0_dp
    local_seed = 12345
    if (present(seed)) local_seed = seed

    do b = 1, n_draws
      if (present(block_length)) then
        call resample_residuals(base%residuals, dummy_y, method, residual_draw, dummy_draw, info, &
          block_length, local_seed + 7919 * b)
      else
        call resample_residuals(base%residuals, dummy_y, method, residual_draw, dummy_draw, info, &
          seed=local_seed + 7919 * b)
      end if
      if (info /= 0) return
      call simulate_var(base%intercept, base%phi, residual_draw, data(:, 1), simulated, info)
      if (info /= 0) return
      call fit_var(simulated, fitted, info)
      if (info /= 0) return
      result%phi_draws(:, :, b) = fitted%phi
      result%sigma_draws(:, :, b) = fitted%sigma
    end do

    allocate(matrix_draws(k, k, n_draws))
    matrix_draws = result%phi_draws
    call percentile_bounds(matrix_draws, 0.025_dp, 0.975_dp, result%lower, result%median, result%upper, info)
  end subroutine bootstrap_var

  subroutine percentile_bounds(draws, lower_probability, upper_probability, lower, median, upper, info)
    real(dp), intent(in) :: draws(:, :, :)
    real(dp), intent(in) :: lower_probability, upper_probability
    real(dp), allocatable, intent(out) :: lower(:, :), median(:, :), upper(:, :)
    integer, intent(out) :: info
    real(dp), allocatable :: values(:)
    integer :: i, j, n, il, im, iu

    n = size(draws, 3)
    if (n < 1 .or. lower_probability < 0.0_dp .or. upper_probability > 1.0_dp .or. &
        lower_probability > upper_probability) then
      allocate(lower(0, 0), median(0, 0), upper(0, 0))
      info = -1
      return
    end if
    allocate(lower(size(draws, 1), size(draws, 2)))
    allocate(median(size(draws, 1), size(draws, 2)))
    allocate(upper(size(draws, 1), size(draws, 2)), values(n))
    il = max(1, min(n, 1 + int(lower_probability * real(n - 1, dp))))
    im = max(1, min(n, 1 + int(0.5_dp * real(n - 1, dp))))
    iu = max(1, min(n, 1 + int(upper_probability * real(n - 1, dp))))
    do j = 1, size(draws, 2)
      do i = 1, size(draws, 1)
        values = draws(i, j, :)
        call sort_values(values)
        lower(i, j) = values(il)
        median(i, j) = values(im)
        upper(i, j) = values(iu)
      end do
    end do
    info = 0
  contains
    subroutine sort_values(x)
      real(dp), intent(inout) :: x(:)
      real(dp) :: key
      integer :: i, j

      do i = 2, size(x)
        key = x(i)
        j = i - 1
        do while (j >= 1)
          if (x(j) <= key) exit
          x(j + 1) = x(j)
          j = j - 1
        end do
        x(j + 1) = key
      end do
    end subroutine sort_values
  end subroutine percentile_bounds

end module multiatsm_bootstrap
