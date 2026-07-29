! SPDX-License-Identifier: GPL-2.0-or-later
module apt_cointegration
  use apt_kinds, only : dp
  use apt_regression, only : regression_result, hypothesis_result, fit_ols, &
    linear_f_test, ljung_box_test
  implicit none
  private

  integer, parameter, public :: apt_tar = 1
  integer, parameter, public :: apt_mtar = 2

  type, public :: ci_tar_fit_result
    integer :: model = apt_tar
    integer :: lag = 0
    integer :: n_total = 0
    integer :: n_cointegration = 0
    integer :: start_index = 0
    integer :: status = 0
    real(dp) :: threshold = 0.0_dp
    real(dp) :: sse = 0.0_dp
    real(dp) :: aic = 0.0_dp
    real(dp) :: bic = 0.0_dp
    real(dp), allocatable :: long_run_residuals(:)
    real(dp), allocatable :: lagged_residuals(:)
    real(dp), allocatable :: lagged_differences(:)
    real(dp), allocatable :: response(:)
    real(dp), allocatable :: design(:,:)
    type(regression_result) :: long_run
    type(regression_result) :: threshold_regression
    type(hypothesis_result) :: no_cointegration_test
    type(hypothesis_result) :: symmetry_test
  end type ci_tar_fit_result

  type, public :: ci_tar_lag_result
    integer :: model = apt_tar
    integer :: max_lag = 0
    integer :: best_lag_aic = 0
    integer :: best_lag_bic = 0
    integer :: status = 0
    real(dp) :: threshold = 0.0_dp
    real(dp) :: best_aic = 0.0_dp
    real(dp) :: best_bic = 0.0_dp
    integer, allocatable :: lag(:)
    integer, allocatable :: total_observations(:)
    integer, allocatable :: cointegration_observations(:)
    real(dp), allocatable :: sse(:)
    real(dp), allocatable :: aic(:)
    real(dp), allocatable :: bic(:)
    real(dp), allocatable :: ljung_box_4_p(:)
    real(dp), allocatable :: ljung_box_8_p(:)
    real(dp), allocatable :: ljung_box_12_p(:)
  end type ci_tar_lag_result

  type, public :: ci_tar_threshold_result
    integer :: model = apt_tar
    integer :: lag = 0
    integer :: n_total = 0
    integer :: n_cointegration = 0
    integer :: lower_index = 0
    integer :: upper_index = 0
    integer :: status = 0
    real(dp) :: trim_fraction = 0.15_dp
    real(dp) :: threshold = 0.0_dp
    real(dp) :: minimum_sse = 0.0_dp
    real(dp), allocatable :: path_threshold(:)
    real(dp), allocatable :: path_sse(:)
    real(dp), allocatable :: path_aic(:)
    real(dp), allocatable :: path_bic(:)
  end type ci_tar_threshold_result

  public :: ci_tar_fit, ci_tar_lag, ci_tar_threshold
  public :: ciTarFit, ciTarLag, ciTarThd

  interface ciTarFit
    module procedure ci_tar_fit
  end interface
  interface ciTarLag
    module procedure ci_tar_lag
  end interface
  interface ciTarThd
    module procedure ci_tar_threshold
  end interface

contains

  subroutine ci_tar_fit(y, x, result, model, lag, threshold, start_index)
    real(dp), intent(in) :: y(:), x(:)
    type(ci_tar_fit_result), intent(out) :: result
    integer, intent(in), optional :: model, lag, start_index
    real(dp), intent(in), optional :: threshold
    integer :: n, mdl, l, t0, t, j, row, p
    real(dp) :: th, lz, indv
    real(dp), allocatable :: xlr(:,:), rmat(:,:), rhs(:)

    n = size(y)
    mdl = apt_tar
    if (present(model)) mdl = model
    l = 1
    if (present(lag)) l = lag
    th = 0.0_dp
    if (present(threshold)) th = threshold

    result%model = mdl
    result%lag = l
    result%threshold = th
    result%n_total = n

    if (size(x) /= n .or. n < 4 .or. l < 0 .or. &
        (mdl /= apt_tar .and. mdl /= apt_mtar)) then
      result%status = 1
      return
    end if

    allocate(xlr(n,2))
    xlr(:,1) = 1.0_dp
    xlr(:,2) = x
    call fit_ols(y, xlr, result%long_run)
    if (result%long_run%status /= 0) then
      result%status = 2
      return
    end if
    allocate(result%long_run_residuals(n), result%lagged_residuals(n-1), &
      result%lagged_differences(max(0,n-2)))
    result%long_run_residuals = result%long_run%residuals
    result%lagged_residuals = result%long_run_residuals(1:n-1)
    if (n > 2) result%lagged_differences = &
      result%long_run_residuals(2:n-1) - result%long_run_residuals(1:n-2)

    t0 = l + 2
    if (mdl == apt_mtar) t0 = max(t0, 3)
    if (present(start_index)) t0 = max(t0, start_index)
    if (t0 > n) then
      result%status = 3
      return
    end if
    result%start_index = t0
    result%n_cointegration = n - t0 + 1
    p = l + 2
    allocate(result%response(result%n_cointegration), &
      result%design(result%n_cointegration,p))

    do t = t0, n
      row = t - t0 + 1
      result%response(row) = result%long_run_residuals(t) - &
        result%long_run_residuals(t-1)
      lz = result%long_run_residuals(t-1)
      if (mdl == apt_tar) then
        if (lz >= th) then
          indv = 1.0_dp
        else
          indv = 0.0_dp
        end if
      else
        if (result%long_run_residuals(t-1) - result%long_run_residuals(t-2) >= th) then
          indv = 1.0_dp
        else
          indv = 0.0_dp
        end if
      end if
      result%design(row,1) = lz * indv
      result%design(row,2) = lz * (1.0_dp - indv)
      do j = 1, l
        result%design(row,2+j) = result%long_run_residuals(t-j) - &
          result%long_run_residuals(t-j-1)
      end do
    end do
    call fit_ols(result%response, result%design, result%threshold_regression)
    if (result%threshold_regression%status /= 0) then
      result%status = 4
      return
    end if

    allocate(rmat(2,p), rhs(2))
    rmat = 0.0_dp
    rhs = 0.0_dp
    rmat(1,1) = 1.0_dp
    rmat(2,2) = 1.0_dp
    call linear_f_test(result%threshold_regression, rmat, rhs, &
      result%no_cointegration_test)
    deallocate(rmat, rhs)
    allocate(rmat(1,p), rhs(1))
    rmat = 0.0_dp
    rhs = 0.0_dp
    rmat(1,1) = 1.0_dp
    rmat(1,2) = -1.0_dp
    call linear_f_test(result%threshold_regression, rmat, rhs, result%symmetry_test)

    result%sse = result%threshold_regression%sse
    result%aic = result%threshold_regression%aic
    result%bic = result%threshold_regression%bic
    result%status = 0
  end subroutine ci_tar_fit

  subroutine ci_tar_lag(y, x, result, model, max_lag, threshold, adjust)
    real(dp), intent(in) :: y(:), x(:)
    type(ci_tar_lag_result), intent(out) :: result
    integer, intent(in), optional :: model, max_lag
    real(dp), intent(in), optional :: threshold
    logical, intent(in), optional :: adjust
    integer :: mdl, ml, k, common_start, best_a, best_b
    real(dp) :: th, q
    logical :: adj
    type(ci_tar_fit_result) :: fit

    mdl = apt_tar
    if (present(model)) mdl = model
    ml = 4
    if (present(max_lag)) ml = max_lag
    th = 0.0_dp
    if (present(threshold)) th = threshold
    adj = .true.
    if (present(adjust)) adj = adjust
    result%model = mdl
    result%max_lag = ml
    result%threshold = th
    if (ml < 0) then
      result%status = 1
      return
    end if
    allocate(result%lag(ml+1), result%total_observations(ml+1), &
      result%cointegration_observations(ml+1), result%sse(ml+1), &
      result%aic(ml+1), result%bic(ml+1), result%ljung_box_4_p(ml+1), &
      result%ljung_box_8_p(ml+1), result%ljung_box_12_p(ml+1))
    common_start = ml + 2
    if (mdl == apt_mtar) common_start = max(common_start, 3)
    do k = 0, ml
      if (adj) then
        call ci_tar_fit(y, x, fit, mdl, k, th, common_start)
      else
        call ci_tar_fit(y, x, fit, mdl, k, th)
      end if
      if (fit%status /= 0) then
        result%status = 2
        return
      end if
      result%lag(k+1) = k
      result%total_observations(k+1) = fit%n_total
      result%cointegration_observations(k+1) = fit%n_cointegration
      result%sse(k+1) = fit%sse
      result%aic(k+1) = fit%aic
      result%bic(k+1) = fit%bic
      call ljung_box_test(fit%threshold_regression%residuals, 4, q, &
        result%ljung_box_4_p(k+1))
      call ljung_box_test(fit%threshold_regression%residuals, 8, q, &
        result%ljung_box_8_p(k+1))
      call ljung_box_test(fit%threshold_regression%residuals, 12, q, &
        result%ljung_box_12_p(k+1))
    end do
    best_a = minloc(result%aic, dim=1)
    best_b = minloc(result%bic, dim=1)
    result%best_lag_aic = result%lag(best_a)
    result%best_lag_bic = result%lag(best_b)
    result%best_aic = result%aic(best_a)
    result%best_bic = result%bic(best_b)
    result%status = 0
  end subroutine ci_tar_lag

  subroutine ci_tar_threshold(y, x, result, model, lag, trim_fraction)
    real(dp), intent(in) :: y(:), x(:)
    type(ci_tar_threshold_result), intent(out) :: result
    integer, intent(in), optional :: model, lag
    real(dp), intent(in), optional :: trim_fraction
    integer :: mdl, l, a, b, i, np, idx
    real(dp) :: trimv
    real(dp), allocatable :: candidates(:)
    type(ci_tar_fit_result) :: base, fit

    mdl = apt_tar
    if (present(model)) mdl = model
    l = 1
    if (present(lag)) l = lag
    trimv = 0.15_dp
    if (present(trim_fraction)) trimv = trim_fraction
    result%model = mdl
    result%lag = l
    result%trim_fraction = trimv
    if (trimv < 0.0_dp .or. trimv >= 0.5_dp) then
      result%status = 1
      return
    end if
    call ci_tar_fit(y, x, base, mdl, l, 0.0_dp)
    if (base%status /= 0) then
      result%status = 2
      return
    end if
    result%n_total = base%n_total
    result%n_cointegration = base%n_cointegration
    if (mdl == apt_tar) then
      candidates = base%lagged_residuals
    else
      candidates = base%lagged_differences
    end if
    call sort_real(candidates)
    a = ceiling(real(base%n_cointegration,dp) * trimv)
    b = base%n_cointegration - floor(real(base%n_cointegration,dp) * trimv)
    a = max(1, min(a, size(candidates)))
    b = max(a, min(b, size(candidates)))
    result%lower_index = a
    result%upper_index = b
    np = b - a + 1
    allocate(result%path_threshold(np), result%path_sse(np), &
      result%path_aic(np), result%path_bic(np))
    do i = a, b
      idx = i - a + 1
      result%path_threshold(idx) = candidates(i)
      call ci_tar_fit(y, x, fit, mdl, l, candidates(i))
      if (fit%status /= 0) then
        result%status = 3
        return
      end if
      result%path_sse(idx) = fit%sse
      result%path_aic(idx) = fit%aic
      result%path_bic(idx) = fit%bic
    end do
    idx = minloc(result%path_sse, dim=1)
    result%threshold = result%path_threshold(idx)
    result%minimum_sse = result%path_sse(idx)
    result%status = 0
  end subroutine ci_tar_threshold

  recursive subroutine quicksort_real(a, left, right)
    real(dp), intent(inout) :: a(:)
    integer, intent(in) :: left, right
    integer :: i, j
    real(dp) :: pivot, tmp
    if (left >= right) return
    i = left
    j = right
    pivot = a((left+right)/2)
    do
      do while (a(i) < pivot)
        i = i + 1
      end do
      do while (a(j) > pivot)
        j = j - 1
      end do
      if (i <= j) then
        tmp = a(i)
        a(i) = a(j)
        a(j) = tmp
        i = i + 1
        j = j - 1
      end if
      if (i > j) exit
    end do
    if (left < j) call quicksort_real(a, left, j)
    if (i < right) call quicksort_real(a, i, right)
  end subroutine quicksort_real

  subroutine sort_real(a)
    real(dp), intent(inout) :: a(:)
    if (size(a) > 1) call quicksort_real(a, 1, size(a))
  end subroutine sort_real

end module apt_cointegration
