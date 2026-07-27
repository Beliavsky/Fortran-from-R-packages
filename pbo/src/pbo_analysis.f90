! SPDX-License-Identifier: MIT
! Copyright (c) 2014 Matthew R. Barry
module pbo_analysis
  use pbo_kinds, only : dp
  use pbo_types, only : pbo_result, dominance_result, selection_result
  use pbo_stats, only : empirical_cdf
  implicit none
  private
  public :: dominance_curve, selection_frequencies
contains
  subroutine dominance_curve(result, increment, curve)
    type(pbo_result), intent(in) :: result
    real(dp), intent(in) :: increment
    type(dominance_result), intent(out) :: curve
    real(dp), allocatable :: all_oos(:)
    real(dp) :: lo, hi
    integer :: n_points, i, c, j, p

    curve%message = ''
    if (.not. result%success) then
      curve%message = 'pbo result is not successful'
      return
    end if
    if (increment <= 0.0_dp) then
      curve%message = 'increment must be positive'
      return
    end if
    lo = minval(result%selected_pairs(:,2))
    hi = maxval(result%selected_pairs(:,2))
    n_points = int(floor((hi - lo) / increment + 32.0_dp * epsilon(1.0_dp))) + 1
    n_points = max(n_points, 1)
    allocate(curve%performance(n_points), curve%cdf_selected(n_points))
    allocate(curve%cdf_all(n_points), curve%sd2_difference(n_points))
    allocate(curve%integrated_difference(n_points))
    allocate(all_oos(result%n_strategies * result%n_cases))
    p = 0
    do c = 1, result%n_cases
      do j = 1, result%n_strategies
        p = p + 1
        all_oos(p) = result%performance_oos(j,c)
      end do
    end do
    do i = 1, n_points
      curve%performance(i) = lo + real(i - 1, dp) * increment
      curve%cdf_selected(i) = empirical_cdf(result%selected_pairs(:,2), curve%performance(i))
      curve%cdf_all(i) = empirical_cdf(all_oos, curve%performance(i))
      curve%sd2_difference(i) = curve%cdf_all(i) - curve%cdf_selected(i)
    end do
    curve%integrated_difference(1) = 0.0_dp
    do i = 2, n_points
      curve%integrated_difference(i) = curve%integrated_difference(i - 1) + &
        0.5_dp * increment * (curve%sd2_difference(i - 1) + curve%sd2_difference(i))
    end do
    curve%success = .true.
    curve%message = 'ok'
  end subroutine dominance_curve

  subroutine selection_frequencies(result, frequencies)
    type(pbo_result), intent(in) :: result
    type(selection_result), intent(out) :: frequencies
    integer :: i, j, key_strategy, key_frequency

    allocate(frequencies%strategy(result%n_strategies))
    allocate(frequencies%frequency(result%n_strategies))
    do i = 1, result%n_strategies
      frequencies%strategy(i) = i
      frequencies%frequency(i) = count(result%selected_is == i)
    end do
    do i = 2, result%n_strategies
      key_strategy = frequencies%strategy(i)
      key_frequency = frequencies%frequency(i)
      j = i - 1
      do while (j >= 1)
        if (frequencies%frequency(j) > key_frequency) exit
        if (frequencies%frequency(j) == key_frequency .and. &
            frequencies%strategy(j) < key_strategy) exit
        frequencies%strategy(j + 1) = frequencies%strategy(j)
        frequencies%frequency(j + 1) = frequencies%frequency(j)
        j = j - 1
      end do
      frequencies%strategy(j + 1) = key_strategy
      frequencies%frequency(j + 1) = key_frequency
    end do
  end subroutine selection_frequencies
end module pbo_analysis
