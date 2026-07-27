! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2005 Antonio Fabio Di Narzo
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of a translation of tseriesChaos and is distributed
! under the GNU General Public License version 2 only.
module chaos_neighbors
  use chaos_kinds, only : dp
  use chaos_utils, only : scale_unit_interval, squared_embedding_distance, sort_pairs, &
    linear_regression, quiet_nan
  use, intrinsic :: iso_fortran_env, only : int64
  implicit none
  private
  integer, parameter :: search_direct = 1
  integer, parameter :: search_box = 2
  integer, parameter :: max_box_bins = 100

  type :: box_index_type
    integer :: nbin = 0
    integer :: npoints = 0
    integer :: offset = 0
    real(dp) :: width = 0.0_dp
    integer, allocatable :: cell_start(:)
    integer, allocatable :: members(:)
  end type box_index_type

  public :: false_nearest_fraction, false_nearest_curve
  public :: find_k_nearests, follow_neighbor_points, lyapunov_stretching
  public :: lyapunov_linear_fit
contains
  subroutine false_nearest_fraction(series, m, d, theiler, rt, eps, fraction, total, status, &
      search_method, distance_evaluations, method_used)
    real(dp), intent(in) :: series(:), rt, eps
    integer, intent(in) :: m, d, theiler
    real(dp), intent(out) :: fraction
    integer, intent(out) :: total, status
    character(len=*), intent(in), optional :: search_method
    integer(int64), intent(out), optional :: distance_evaluations
    character(len=*), intent(out), optional :: method_used
    type(box_index_type) :: box_index
    real(dp), allocatable :: x(:), distances(:)
    integer, allocatable :: ids(:)
    real(dp) :: span, eps_scaled, dist2, ratio
    integer :: scale_status, blength, i, j, nfalse, nfound, method, method_status
    integer(int64) :: evaluations

    if (present(distance_evaluations)) distance_evaluations = 0_int64
    if (present(method_used)) method_used = ""
    if (m < 1 .or. d < 1 .or. theiler < 0 .or. rt <= 0.0_dp .or. eps <= 0.0_dp) then
      fraction = quiet_nan()
      total = 0
      status = 1
      return
    end if
    call scale_unit_interval(series, x, span, scale_status)
    if (scale_status /= 0) then
      fraction = quiet_nan()
      total = 0
      status = 2
      return
    end if
    eps_scaled = eps / span
    blength = size(series) - m * d - theiler
    if (blength <= 0) then
      fraction = quiet_nan()
      total = 0
      status = 3
      return
    end if
    call resolve_search_method(search_method, blength, m, eps_scaled, method, method_status)
    if (method_status /= 0) then
      fraction = quiet_nan()
      total = 0
      status = 5
      return
    end if
    call report_method(method, method_used)
    if (method == search_box) then
      call build_box_index(x, m, d, blength, eps_scaled, box_index, method_status)
      if (method_status /= 0) then
        fraction = quiet_nan()
        total = 0
        status = 6
        return
      end if
    end if

    allocate(ids(blength), distances(blength))
    nfalse = 0
    total = 0
    evaluations = 0_int64
    do i = 1, blength
      call collect_radius_neighbors(x, m, d, theiler, eps_scaled, blength, i, method, &
        box_index, ids, distances, nfound, evaluations)
      do j = 1, nfound
        dist2 = distances(j)**2
        if (dist2 <= tiny(1.0_dp)) cycle
        ratio = (dist2 + (x(i + m * d) - x(ids(j) + m * d))**2) / dist2
        if (ratio > rt) nfalse = nfalse + 1
        total = total + 1
      end do
    end do
    if (present(distance_evaluations)) distance_evaluations = evaluations
    if (total == 0) then
      fraction = quiet_nan()
      status = 4
    else
      fraction = real(nfalse, dp) / real(total, dp)
      status = 0
    end if
  end subroutine false_nearest_fraction

  subroutine false_nearest_curve(series, max_m, d, theiler, rt, eps, fractions, totals, status, &
      search_method)
    real(dp), intent(in) :: series(:), rt, eps
    integer, intent(in) :: max_m, d, theiler
    real(dp), allocatable, intent(out) :: fractions(:)
    integer, allocatable, intent(out) :: totals(:)
    integer, intent(out) :: status
    character(len=*), intent(in), optional :: search_method
    integer :: m, local_status

    if (max_m < 1) then
      allocate(fractions(0), totals(0))
      status = 1
      return
    end if
    allocate(fractions(max_m), totals(max_m))
    status = 0
    do m = 1, max_m
      call false_nearest_fraction(series, m, d, theiler, rt, eps, fractions(m), totals(m), &
        local_status, search_method=search_method)
      if (local_status /= 0 .and. local_status /= 4) status = local_status
    end do
  end subroutine false_nearest_curve

  subroutine find_k_nearests(series, m, d, theiler, eps, nref, k, steps, nearest, distances, &
      status, search_method, distance_evaluations, method_used)
    real(dp), intent(in) :: series(:), eps
    integer, intent(in) :: m, d, theiler, nref, k, steps
    integer, allocatable, intent(out) :: nearest(:, :)
    real(dp), allocatable, intent(out) :: distances(:, :)
    integer, intent(out) :: status
    character(len=*), intent(in), optional :: search_method
    integer(int64), intent(out), optional :: distance_evaluations
    character(len=*), intent(out), optional :: method_used
    type(box_index_type) :: box_index
    real(dp), allocatable :: x(:), candidate_dist(:)
    integer, allocatable :: candidate_id(:)
    real(dp) :: span, eps_scaled
    integer :: scale_status, blength, i, nfound, nkeep, method, method_status
    integer(int64) :: evaluations

    if (present(distance_evaluations)) distance_evaluations = 0_int64
    if (present(method_used)) method_used = ""
    if (m < 1 .or. d < 1 .or. theiler < 0 .or. eps <= 0.0_dp .or. &
        nref < 1 .or. k < 1 .or. steps < 0) then
      allocate(nearest(0, 0), distances(0, 0))
      status = 1
      return
    end if
    call scale_unit_interval(series, x, span, scale_status)
    if (scale_status /= 0) then
      allocate(nearest(0, 0), distances(0, 0))
      status = 2
      return
    end if
    eps_scaled = eps / span
    blength = size(series) - (m - 1) * d - steps
    if (blength <= 0 .or. nref > blength) then
      allocate(nearest(0, 0), distances(0, 0))
      status = 3
      return
    end if
    call resolve_search_method(search_method, blength, m, eps_scaled, method, method_status)
    if (method_status /= 0) then
      allocate(nearest(0, 0), distances(0, 0))
      status = 4
      return
    end if
    call report_method(method, method_used)
    if (method == search_box) then
      call build_box_index(x, m, d, blength, eps_scaled, box_index, method_status)
      if (method_status /= 0) then
        allocate(nearest(0, 0), distances(0, 0))
        status = 5
        return
      end if
    end if

    allocate(nearest(nref, k), distances(nref, k))
    allocate(candidate_dist(blength), candidate_id(blength))
    nearest = -1
    distances = huge(1.0_dp)
    evaluations = 0_int64
    do i = 1, nref
      call collect_radius_neighbors(x, m, d, theiler, eps_scaled, blength, i, method, &
        box_index, candidate_id, candidate_dist, nfound, evaluations)
      if (nfound > 1) call sort_pairs(candidate_dist, candidate_id, nfound)
      nkeep = min(k, nfound)
      if (nkeep > 0) then
        nearest(i, 1:nkeep) = candidate_id(1:nkeep)
        distances(i, 1:nkeep) = candidate_dist(1:nkeep) * span
      end if
    end do
    if (present(distance_evaluations)) distance_evaluations = evaluations
    status = 0
  end subroutine find_k_nearests

  subroutine follow_neighbor_points(series, m, d, refs, nearest, steps, stretching, status)
    real(dp), intent(in) :: series(:)
    integer, intent(in) :: m, d, refs(:), nearest(:, :), steps
    real(dp), allocatable, intent(out) :: stretching(:)
    integer, intent(out) :: status
    integer :: time, i, j, a, b, k
    real(dp) :: mean_distance

    if (m < 1 .or. d < 1 .or. steps < 1 .or. size(refs) < 1 .or. &
        size(nearest, 1) < maxval(refs)) then
      allocate(stretching(0))
      status = 1
      return
    end if
    k = size(nearest, 2)
    if (k < 1 .or. any(nearest(refs, :) < 1)) then
      allocate(stretching(0))
      status = 2
      return
    end if
    if (max(maxval(refs), maxval(nearest(refs, :))) + steps - 1 + &
        (m - 1) * d > size(series)) then
      allocate(stretching(0))
      status = 3
      return
    end if
    allocate(stretching(steps))
    stretching = 0.0_dp
    do time = 0, steps - 1
      do i = 1, size(refs)
        mean_distance = 0.0_dp
        a = refs(i) + time
        do j = 1, k
          b = nearest(refs(i), j) + time
          mean_distance = mean_distance + &
            sqrt(squared_embedding_distance(series, a, b, m, d))
        end do
        mean_distance = mean_distance / real(k, dp)
        if (mean_distance <= 0.0_dp) then
          stretching(time + 1) = quiet_nan()
          status = 4
          return
        end if
        stretching(time + 1) = stretching(time + 1) + log(mean_distance)
      end do
      stretching(time + 1) = stretching(time + 1) / real(size(refs), dp)
    end do
    status = 0
  end subroutine follow_neighbor_points

  subroutine lyapunov_stretching(series, m, d, theiler, k, max_ref, steps, eps, stretching, &
      refs_used, status, search_method, distance_evaluations, method_used)
    real(dp), intent(in) :: series(:), eps
    integer, intent(in) :: m, d, theiler, k, max_ref, steps
    real(dp), allocatable, intent(out) :: stretching(:)
    integer, intent(out) :: refs_used, status
    character(len=*), intent(in), optional :: search_method
    integer(int64), intent(out), optional :: distance_evaluations
    character(len=*), intent(out), optional :: method_used
    integer, allocatable :: nearest(:, :), refs(:)
    real(dp), allocatable :: distances(:, :)
    integer :: nref, i, local_status
    integer(int64) :: evaluations
    character(len=16) :: selected_method

    if (present(distance_evaluations)) distance_evaluations = 0_int64
    if (present(method_used)) method_used = ""
    nref = min(max_ref, size(series) - (m - 1) * d - steps)
    if (max_ref < 0) nref = size(series) - (m - 1) * d - steps
    if (nref < 1) then
      allocate(stretching(0))
      refs_used = 0
      status = 1
      return
    end if
    call find_k_nearests(series, m, d, theiler, eps, nref, k, steps, nearest, distances, &
      local_status, search_method=search_method, distance_evaluations=evaluations, &
      method_used=selected_method)
    if (present(distance_evaluations)) distance_evaluations = evaluations
    if (present(method_used)) method_used = trim(selected_method)
    if (local_status /= 0) then
      allocate(stretching(0))
      refs_used = 0
      status = local_status
      return
    end if
    refs_used = count([(all(nearest(i, :) > 0), i = 1, nref)])
    if (refs_used < 1) then
      allocate(stretching(0))
      status = 2
      return
    end if
    allocate(refs(refs_used))
    refs_used = 0
    do i = 1, nref
      if (all(nearest(i, :) > 0)) then
        refs_used = refs_used + 1
        refs(refs_used) = i
      end if
    end do
    call follow_neighbor_points(series, m, d, refs, nearest, steps, stretching, status)
  end subroutine lyapunov_stretching

  subroutine lyapunov_linear_fit(stretching, first_step, last_step, intercept, exponent, &
      status, dt)
    real(dp), intent(in) :: stretching(:)
    integer, intent(in) :: first_step, last_step
    real(dp), intent(out) :: intercept, exponent
    integer, intent(out) :: status
    real(dp), intent(in), optional :: dt
    real(dp), allocatable :: x(:)
    real(dp) :: delta
    integer :: i, n

    if (first_step < 1 .or. last_step > size(stretching) .or. last_step <= first_step) then
      intercept = quiet_nan()
      exponent = quiet_nan()
      status = 1
      return
    end if
    delta = 1.0_dp
    if (present(dt)) delta = dt
    if (delta <= 0.0_dp) then
      intercept = quiet_nan()
      exponent = quiet_nan()
      status = 2
      return
    end if
    n = last_step - first_step + 1
    allocate(x(n))
    do i = 1, n
      x(i) = real(first_step + i - 2, dp) * delta
    end do
    call linear_regression(x, stretching(first_step:last_step), intercept, exponent, status)
  end subroutine lyapunov_linear_fit

  subroutine resolve_search_method(requested_method, npoints, m, eps, method, status)
    character(len=*), intent(in), optional :: requested_method
    integer, intent(in) :: npoints, m
    real(dp), intent(in) :: eps
    integer, intent(out) :: method, status
    character(len=16) :: requested

    requested = "auto"
    if (present(requested_method)) requested = lowercase(trim(adjustl(requested_method)))
    select case (trim(requested))
    case ("direct")
      method = search_direct
      status = 0
    case ("box")
      method = search_box
      status = 0
    case ("auto")
      if (npoints >= 256 .and. m <= 8 .and. eps < 0.5_dp) then
        method = search_box
      else
        method = search_direct
      end if
      status = 0
    case default
      method = 0
      status = 1
    end select
  end subroutine resolve_search_method

  subroutine report_method(method, method_used)
    integer, intent(in) :: method
    character(len=*), intent(out), optional :: method_used

    if (.not. present(method_used)) return
    if (method == search_box) then
      method_used = "box"
    else
      method_used = "direct"
    end if
  end subroutine report_method

  subroutine build_box_index(x, m, d, npoints, eps, index, status)
    real(dp), intent(in) :: x(:), eps
    integer, intent(in) :: m, d, npoints
    type(box_index_type), intent(out) :: index
    integer, intent(out) :: status
    integer, allocatable :: counts(:), cursor(:)
    integer :: i, cell, ncells, bx, by

    if (m < 1 .or. d < 1 .or. npoints < 1 .or. eps <= 0.0_dp) then
      status = 1
      return
    end if
    index%width = max(eps, 1.0_dp / real(max_box_bins, dp))
    index%nbin = max(1, min(max_box_bins, ceiling(1.0_dp / index%width)))
    index%npoints = npoints
    index%offset = (m - 1) * d
    ncells = index%nbin * index%nbin
    allocate(counts(ncells), cursor(ncells), index%cell_start(ncells + 1))
    allocate(index%members(npoints))
    counts = 0
    do i = 1, npoints
      bx = coordinate_bin(x(i), index%width, index%nbin)
      by = coordinate_bin(x(i + index%offset), index%width, index%nbin)
      cell = flatten_cell(bx, by, index%nbin)
      counts(cell) = counts(cell) + 1
    end do
    index%cell_start(1) = 1
    do cell = 1, ncells
      index%cell_start(cell + 1) = index%cell_start(cell) + counts(cell)
    end do
    cursor = index%cell_start(1:ncells)
    do i = 1, npoints
      bx = coordinate_bin(x(i), index%width, index%nbin)
      by = coordinate_bin(x(i + index%offset), index%width, index%nbin)
      cell = flatten_cell(bx, by, index%nbin)
      index%members(cursor(cell)) = i
      cursor(cell) = cursor(cell) + 1
    end do
    status = 0
  end subroutine build_box_index

  subroutine collect_radius_neighbors(x, m, d, theiler, eps, npoints, reference, method, &
      box_index, ids, distances, nfound, evaluations)
    real(dp), intent(in) :: x(:), eps
    integer, intent(in) :: m, d, theiler, npoints, reference, method
    type(box_index_type), intent(in) :: box_index
    integer, intent(out) :: ids(:)
    real(dp), intent(out) :: distances(:)
    integer, intent(out) :: nfound
    integer(int64), intent(inout) :: evaluations

    if (method == search_box) then
      call collect_radius_neighbors_box(x, m, d, theiler, eps, reference, box_index, &
        ids, distances, nfound, evaluations)
    else
      call collect_radius_neighbors_direct(x, m, d, theiler, eps, npoints, reference, &
        ids, distances, nfound, evaluations)
    end if
  end subroutine collect_radius_neighbors

  subroutine collect_radius_neighbors_direct(x, m, d, theiler, eps, npoints, reference, &
      ids, distances, nfound, evaluations)
    real(dp), intent(in) :: x(:), eps
    integer, intent(in) :: m, d, theiler, npoints, reference
    integer, intent(out) :: ids(:)
    real(dp), intent(out) :: distances(:)
    integer, intent(out) :: nfound
    integer(int64), intent(inout) :: evaluations
    real(dp) :: dist2, eps2
    integer :: candidate

    eps2 = eps**2
    nfound = 0
    do candidate = 1, npoints
      if (abs(candidate - reference) <= theiler) cycle
      evaluations = evaluations + 1_int64
      dist2 = squared_embedding_distance(x, reference, candidate, m, d)
      if (dist2 < eps2) then
        nfound = nfound + 1
        ids(nfound) = candidate
        distances(nfound) = sqrt(dist2)
      end if
    end do
  end subroutine collect_radius_neighbors_direct

  subroutine collect_radius_neighbors_box(x, m, d, theiler, eps, reference, index, ids, &
      distances, nfound, evaluations)
    real(dp), intent(in) :: x(:), eps
    integer, intent(in) :: m, d, theiler, reference
    type(box_index_type), intent(in) :: index
    integer, intent(out) :: ids(:)
    real(dp), intent(out) :: distances(:)
    integer, intent(out) :: nfound
    integer(int64), intent(inout) :: evaluations
    real(dp) :: dist2, eps2
    integer :: ref_bx, ref_by, bx, by, cell, pos, candidate

    eps2 = eps**2
    ref_bx = coordinate_bin(x(reference), index%width, index%nbin)
    ref_by = coordinate_bin(x(reference + index%offset), index%width, index%nbin)
    nfound = 0
    do bx = max(1, ref_bx - 1), min(index%nbin, ref_bx + 1)
      do by = max(1, ref_by - 1), min(index%nbin, ref_by + 1)
        cell = flatten_cell(bx, by, index%nbin)
        do pos = index%cell_start(cell), index%cell_start(cell + 1) - 1
          candidate = index%members(pos)
          if (abs(candidate - reference) <= theiler) cycle
          evaluations = evaluations + 1_int64
          dist2 = squared_embedding_distance(x, reference, candidate, m, d)
          if (dist2 < eps2) then
            nfound = nfound + 1
            ids(nfound) = candidate
            distances(nfound) = sqrt(dist2)
          end if
        end do
      end do
    end do
  end subroutine collect_radius_neighbors_box

  pure integer function coordinate_bin(value, width, nbin) result(bin)
    real(dp), intent(in) :: value, width
    integer, intent(in) :: nbin

    bin = int(value / width) + 1
    bin = max(1, min(nbin, bin))
  end function coordinate_bin

  pure integer function flatten_cell(bx, by, nbin) result(cell)
    integer, intent(in) :: bx, by, nbin

    cell = (bx - 1) * nbin + by
  end function flatten_cell

  pure function lowercase(text) result(lower)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: lower
    integer :: i, code

    lower = text
    do i = 1, len(text)
      code = iachar(text(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) then
        lower(i:i) = achar(code + iachar('a') - iachar('A'))
      end if
    end do
  end function lowercase
end module chaos_neighbors
