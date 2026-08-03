! SPDX-License-Identifier: MIT
! PortfolioTesteR modern Fortran translation
module ptr_filters
  use ptr_kinds, only : dp
  use ptr_utils, only : is_finite, percentile
  implicit none
  private
  public :: filter_top_n, filter_threshold, filter_between, filter_by_percentile
  public :: combine_filters, apply_regime, selection_counts

contains

  subroutine filter_top_n(signal, n_select, selected, ascending)
    real(dp), intent(in) :: signal(:,:)
    integer, intent(in) :: n_select
    real(dp), allocatable, intent(out) :: selected(:,:)
    logical, intent(in), optional :: ascending
    logical :: asc
    integer :: t, j, k, best, nleft, nkeep
    logical, allocatable :: used(:)
    real(dp) :: best_value
    asc = .false.; if (present(ascending)) asc = ascending
    allocate(selected(size(signal,1),size(signal,2))); selected = 0.0_dp
    allocate(used(size(signal,2)))
    do t = 1, size(signal,1)
      used = .false.; nleft = count(is_finite(signal(t,:)))
      nkeep = min(max(n_select,0), nleft)
      do k = 1, nkeep
        best = 0
        if (asc) then
          best_value = huge(1.0_dp)
        else
          best_value = -huge(1.0_dp)
        end if
        do j = 1, size(signal,2)
          if (used(j) .or. .not. is_finite(signal(t,j))) cycle
          if (best == 0) then
            best = j; best_value = signal(t,j)
          else if (asc .and. signal(t,j) < best_value) then
            best = j; best_value = signal(t,j)
          else if (.not. asc .and. signal(t,j) > best_value) then
            best = j; best_value = signal(t,j)
          end if
        end do
        if (best > 0) then
          used(best) = .true.; selected(t,best) = 1.0_dp
        end if
      end do
    end do
  end subroutine filter_top_n

  subroutine filter_threshold(signal, value, selected, above)
    real(dp), intent(in) :: signal(:,:), value
    real(dp), allocatable, intent(out) :: selected(:,:)
    logical, intent(in), optional :: above
    logical :: use_above
    use_above = .true.; if (present(above)) use_above = above
    allocate(selected(size(signal,1),size(signal,2))); selected = 0.0_dp
    if (use_above) then
      where (is_finite(signal) .and. signal > value) selected = 1.0_dp
    else
      where (is_finite(signal) .and. signal < value) selected = 1.0_dp
    end if
  end subroutine filter_threshold

  subroutine filter_between(signal, lower, upper, selected, inside)
    real(dp), intent(in) :: signal(:,:), lower, upper
    real(dp), allocatable, intent(out) :: selected(:,:)
    logical, intent(in), optional :: inside
    logical :: use_inside
    use_inside = .true.; if (present(inside)) use_inside = inside
    allocate(selected(size(signal,1),size(signal,2))); selected = 0.0_dp
    if (use_inside) then
      where (is_finite(signal) .and. signal >= lower .and. signal <= upper) selected = 1.0_dp
    else
      where (is_finite(signal) .and. (signal < lower .or. signal > upper)) selected = 1.0_dp
    end if
  end subroutine filter_between

  subroutine filter_by_percentile(signal, prob, selected, top)
    real(dp), intent(in) :: signal(:,:), prob
    real(dp), allocatable, intent(out) :: selected(:,:)
    logical, intent(in), optional :: top
    logical :: use_top
    real(dp) :: cut
    integer :: t
    use_top = .true.; if (present(top)) use_top = top
    allocate(selected(size(signal,1),size(signal,2))); selected = 0.0_dp
    do t = 1, size(signal,1)
      if (use_top) then
        cut = percentile(signal(t,:), 1.0_dp-prob)
        where (is_finite(signal(t,:)) .and. signal(t,:) >= cut) selected(t,:) = 1.0_dp
      else
        cut = percentile(signal(t,:), prob)
        where (is_finite(signal(t,:)) .and. signal(t,:) <= cut) selected(t,:) = 1.0_dp
      end if
    end do
  end subroutine filter_by_percentile

  subroutine combine_filters(filters, selected, use_and)
    real(dp), intent(in) :: filters(:,:,:)
    real(dp), allocatable, intent(out) :: selected(:,:)
    logical, intent(in), optional :: use_and
    logical :: all_required
    integer :: t,j,k
    all_required = .true.; if (present(use_and)) all_required = use_and
    allocate(selected(size(filters,1),size(filters,2)))
    if (all_required) then
      selected = 1.0_dp
      do k = 1, size(filters,3)
        where (filters(:,:,k) <= 0.0_dp) selected = 0.0_dp
      end do
    else
      selected = 0.0_dp
      do k = 1, size(filters,3)
        where (filters(:,:,k) > 0.0_dp) selected = 1.0_dp
      end do
    end if
    do t=1,size(selected,1)
      do j=1,size(selected,2)
        if (.not. is_finite(selected(t,j))) selected(t,j)=0.0_dp
      end do
    end do
  end subroutine combine_filters

  subroutine apply_regime(selection, regime, out, partial_weight)
    real(dp), intent(in) :: selection(:,:), regime(:)
    real(dp), allocatable, intent(out) :: out(:,:)
    real(dp), intent(in), optional :: partial_weight
    real(dp) :: partial
    integer :: t
    partial = 0.0_dp; if (present(partial_weight)) partial = partial_weight
    allocate(out(size(selection,1),size(selection,2))); out = 0.0_dp
    if (size(regime) /= size(selection,1)) return
    do t=1,size(selection,1)
      if (.not. is_finite(regime(t))) cycle
      if (regime(t) > 0.0_dp) then
        out(t,:) = selection(t,:)
      else
        out(t,:) = partial * selection(t,:)
      end if
    end do
  end subroutine apply_regime

  subroutine selection_counts(selection, counts)
    real(dp), intent(in) :: selection(:,:)
    integer, allocatable, intent(out) :: counts(:)
    integer :: t
    allocate(counts(size(selection,1)))
    do t=1,size(selection,1)
      counts(t)=count(is_finite(selection(t,:)) .and. selection(t,:)>0.0_dp)
    end do
  end subroutine selection_counts

end module ptr_filters
