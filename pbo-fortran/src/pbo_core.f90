! SPDX-License-Identifier: MIT
! Copyright (c) 2014 Matthew R. Barry
module pbo_core
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
  use pbo_kinds, only : dp
  use pbo_types, only : pbo_result
  use pbo_combinations, only : generate_combinations
  use pbo_stats, only : first_argmax, average_rank_of, fit_line, round_significant
  implicit none
  private
  public :: performance_function, compute_pbo

  abstract interface
    subroutine performance_function(data, values)
      import dp
      real(dp), intent(in) :: data(:,:)
      real(dp), intent(out) :: values(:)
    end subroutine performance_function
  end interface
contains
  subroutine compute_pbo(data, n_subsets, metric, result, threshold, inf_sub)
    real(dp), intent(in) :: data(:,:)
    integer, intent(in) :: n_subsets
    procedure(performance_function) :: metric
    type(pbo_result), intent(out) :: result
    real(dp), intent(in), optional :: threshold, inf_sub
    integer :: t, n, subset_size, half_rows, c, b, p, row, i
    integer, allocatable :: is_indices(:), oos_indices(:)
    logical, allocatable :: in_sample(:)
    real(dp), allocatable :: is_data(:,:), oos_data(:,:), r(:), r_bar(:)
    real(dp) :: omega, legacy_a, legacy_b, legacy_r2, legacy_ar2
    real(dp) :: deg_a, deg_b, deg_r2, deg_ar2
    logical :: ok, fit_ok
    character(len=:), allocatable :: msg

    result%message = ''
    if (present(threshold)) result%threshold = threshold
    if (present(inf_sub)) result%inf_sub = inf_sub
    t = size(data,1)
    n = size(data,2)
    result%n_observations = t
    result%n_strategies = n
    result%n_subsets = n_subsets

    if (t < 2 .or. n < 1) then
      result%message = 'data must have at least two rows and one strategy'
      return
    end if
    if (n_subsets < 2 .or. mod(n_subsets,2) /= 0) then
      result%message = 'n_subsets must be a positive even integer of at least two'
      return
    end if
    if (mod(t,n_subsets) /= 0) then
      result%message = 'n_subsets must evenly divide the number of observations'
      return
    end if
    if (result%inf_sub <= 0.0_dp) then
      result%message = 'inf_sub must be positive'
      return
    end if

    call generate_combinations(n_subsets, n_subsets / 2, result%combos, ok, msg)
    if (.not. ok) then
      result%message = msg
      return
    end if
    result%n_cases = size(result%combos,2)
    subset_size = t / n_subsets
    half_rows = t / 2
    allocate(result%performance_is(n,result%n_cases))
    allocate(result%performance_oos(n,result%n_cases))
    allocate(result%selected_is(result%n_cases))
    allocate(result%selected_oos(result%n_cases))
    allocate(result%oos_rank(result%n_cases))
    allocate(result%omega_bar(result%n_cases))
    allocate(result%lambda(result%n_cases))
    allocate(result%selected_pairs(result%n_cases,2))
    allocate(is_indices(half_rows), oos_indices(half_rows), in_sample(t))
    allocate(is_data(half_rows,n), oos_data(half_rows,n), r(n), r_bar(n))

    do c = 1, result%n_cases
      in_sample = .false.
      p = 0
      do i = 1, n_subsets / 2
        b = result%combos(i,c)
        do row = (b - 1) * subset_size + 1, b * subset_size
          p = p + 1
          is_indices(p) = row
          in_sample(row) = .true.
        end do
      end do
      p = 0
      do row = 1, t
        if (.not. in_sample(row)) then
          p = p + 1
          oos_indices(p) = row
        end if
      end do
      is_data = data(is_indices,:)
      oos_data = data(oos_indices,:)
      call metric(is_data, r)
      call metric(oos_data, r_bar)
      if (any(ieee_is_nan(r)) .or. any(ieee_is_nan(r_bar))) then
        result%message = 'performance callback returned NaN'
        return
      end if
      result%performance_is(:,c) = r
      result%performance_oos(:,c) = r_bar
      result%selected_is(c) = first_argmax(r)
      result%selected_oos(c) = first_argmax(r_bar)
      result%oos_rank(c) = average_rank_of(r_bar, result%selected_is(c))
      omega = result%oos_rank(c) / real(n, dp)
      result%omega_bar(c) = omega
      if (omega >= 1.0_dp) then
        result%lambda(c) = result%inf_sub
      else
        result%lambda(c) = log(omega / (1.0_dp - omega))
      end if
      result%selected_pairs(c,1) = r(result%selected_is(c))
      result%selected_pairs(c,2) = r_bar(result%selected_is(c))
    end do

    result%phi = real(count(result%lambda <= 0.0_dp), dp) / real(result%n_cases, dp)
    result%below_threshold = real(count(result%selected_pairs(:,2) < result%threshold), dp) / &
      real(result%n_cases, dp)
    result%below_threshold = round_significant(result%below_threshold, 3)

    call fit_line(result%selected_pairs(:,2), result%selected_pairs(:,1), legacy_a, legacy_b, &
      legacy_r2, legacy_ar2, fit_ok)
    if (fit_ok) then
      result%slope = round_significant(legacy_a, 5)
      result%intercept = round_significant(legacy_b, 5)
      result%adjusted_r2 = round_significant(legacy_ar2, 2)
    end if

    call fit_line(result%selected_pairs(:,1), result%selected_pairs(:,2), deg_a, deg_b, deg_r2, &
      deg_ar2, fit_ok)
    if (fit_ok) then
      result%degradation_intercept = deg_a
      result%degradation_slope = deg_b
      result%degradation_r2 = deg_r2
    end if

    if (.not. all(ieee_is_finite(result%selected_pairs))) then
      result%message = 'selected performance pairs must be finite for regression summaries'
      return
    end if
    result%success = .true.
    result%message = 'ok'
  end subroutine compute_pbo
end module pbo_core
