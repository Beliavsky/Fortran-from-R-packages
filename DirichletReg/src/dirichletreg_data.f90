! SPDX-License-Identifier: GPL-2.0-or-later
module dirichletreg_data
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
  use dirichletreg_kinds, only : dp
  implicit none
  private
  public :: prepare_composition, prepare_beta, collapse_subcomposition

contains

  subroutine prepare_composition(y_in, y, base, normalized, transformed, stat, trafo, force_transform, norm_tol)
    real(dp), intent(in) :: y_in(:,:)
    real(dp), intent(out) :: y(:,:)
    integer, intent(in), optional :: base
    logical, intent(out), optional :: normalized, transformed
    integer, intent(out), optional :: stat
    real(dp), intent(in), optional :: trafo, norm_tol
    logical, intent(in), optional :: force_transform

    integer :: i, nvalid, ibase
    real(dp) :: threshold, tol, s, nanv
    logical :: normed, trans, do_transform, valid_row

    if (present(stat)) stat = 0
    normed = .false.
    trans = .false.
    threshold = sqrt(epsilon(1.0_dp))
    if (present(trafo)) threshold = trafo
    tol = sqrt(epsilon(1.0_dp))
    if (present(norm_tol)) tol = norm_tol
    ibase = 1
    if (present(base)) ibase = base

    if (any(shape(y) /= shape(y_in)) .or. size(y_in,2) < 2 .or. ibase < 1 .or. ibase > size(y_in,2) .or. &
        threshold < 0.0_dp .or. threshold >= 0.5_dp .or. tol < 0.0_dp) then
      if (present(stat)) stat = 1
      y = 0.0_dp
      if (present(normalized)) normalized = .false.
      if (present(transformed)) transformed = .false.
      return
    end if

    y = y_in
    nanv = ieee_value(0.0_dp, ieee_quiet_nan)
    nvalid = 0

    do i = 1, size(y,1)
      valid_row = all(ieee_is_finite(y(i,:)))
      if (.not. valid_row) then
        y(i,:) = nanv
        cycle
      end if
      if (any(y(i,:) < 0.0_dp)) then
        if (present(stat)) stat = 2
        return
      end if
      nvalid = nvalid + 1
      s = sum(y(i,:))
      if (s <= 0.0_dp) then
        if (present(stat)) stat = 3
        return
      end if
      if (abs(s-1.0_dp) > tol) then
        y(i,:) = y(i,:)/s
        normed = .true.
      end if
    end do

    if (nvalid == 0) then
      if (present(stat)) stat = 4
      return
    end if

    do_transform = .false.
    if (present(force_transform)) do_transform = force_transform
    if (.not. do_transform) then
      do i = 1, size(y,1)
        if (.not. all(ieee_is_finite(y(i,:)))) cycle
        if (any(y(i,:) < threshold) .or. any(y(i,:) > 1.0_dp-threshold)) then
          do_transform = .true.
          exit
        end if
      end do
    end if

    if (do_transform) then
      do i = 1, size(y,1)
        if (.not. all(ieee_is_finite(y(i,:)))) cycle
        y(i,:) = (y(i,:)*(real(nvalid,dp)-1.0_dp) + 1.0_dp/real(size(y,2),dp))/real(nvalid,dp)
      end do
      trans = .true.
    end if

    do i = 1, size(y,1)
      if (.not. all(ieee_is_finite(y(i,:)))) cycle
      if (any(y(i,:) <= 0.0_dp) .or. any(y(i,:) >= 1.0_dp)) then
        if (present(stat)) stat = 5
        return
      end if
    end do

    if (present(normalized)) normalized = normed
    if (present(transformed)) transformed = trans
  end subroutine prepare_composition


  subroutine prepare_beta(v, y, stat, trafo, force_transform)
    real(dp), intent(in) :: v(:)
    real(dp), intent(out) :: y(:,:)
    integer, intent(out), optional :: stat
    real(dp), intent(in), optional :: trafo
    logical, intent(in), optional :: force_transform
    real(dp), allocatable :: raw(:,:)
    integer :: ierr

    if (size(y,1) /= size(v) .or. size(y,2) /= 2 .or. any(v < 0.0_dp) .or. any(v > 1.0_dp)) then
      if (present(stat)) stat = 1
      y = 0.0_dp
      return
    end if
    allocate(raw(size(v),2))
    raw(:,1) = 1.0_dp-v
    raw(:,2) = v
    if (present(trafo) .and. present(force_transform)) then
      call prepare_composition(raw, y, stat=ierr, trafo=trafo, force_transform=force_transform)
    else if (present(trafo)) then
      call prepare_composition(raw, y, stat=ierr, trafo=trafo)
    else if (present(force_transform)) then
      call prepare_composition(raw, y, stat=ierr, force_transform=force_transform)
    else
      call prepare_composition(raw, y, stat=ierr)
    end if
    if (present(stat)) stat = ierr
  end subroutine prepare_beta


  subroutine collapse_subcomposition(y, selected, ysub, stat)
    real(dp), intent(in) :: y(:,:)
    integer, intent(in) :: selected(:)
    real(dp), intent(out) :: ysub(:,:)
    integer, intent(out), optional :: stat
    logical, allocatable :: keep(:)
    integer :: d, j, k

    if (present(stat)) stat = 0
    d = size(y,2)
    if (size(selected) < 1 .or. size(selected) > d-2 .or. size(ysub,1) /= size(y,1) .or. &
        size(ysub,2) /= size(selected)+1 .or. any(selected < 1) .or. any(selected > d)) then
      if (present(stat)) stat = 1
      ysub = 0.0_dp
      return
    end if
    allocate(keep(d)); keep = .false.
    do j = 1, size(selected)
      if (keep(selected(j))) then
        if (present(stat)) stat = 2
        ysub = 0.0_dp
        return
      end if
      keep(selected(j)) = .true.
    end do
    ysub(:,1) = 0.0_dp
    do k = 1, d
      if (.not. keep(k)) ysub(:,1) = ysub(:,1) + y(:,k)
    end do
    do j = 1, size(selected)
      ysub(:,j+1) = y(:,selected(j))
    end do
  end subroutine collapse_subcomposition

end module dirichletreg_data
