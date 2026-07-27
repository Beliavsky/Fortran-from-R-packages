! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2005 Antonio Fabio Di Narzo
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of a translation of fNonlinear/tseriesChaos computational routines and is distributed
! under the GNU General Public License version 2 or later.
module chaos_metrics
  use chaos_kinds, only : dp
  use chaos_utils, only : embedding_length, scale_unit_interval, squared_embedding_distance, quiet_nan
  use chaos_embedding, only : delay_embed
  use, intrinsic :: iso_fortran_env, only : int64
  implicit none
  private
  public :: correlation_integral, correlation_dimension_curve
  public :: average_mutual_information, recurrence_distance_matrix
  public :: space_time_separation
contains
  subroutine correlation_integral(series, m, d, theiler, eps, value, status)
    real(dp), intent(in) :: series(:), eps
    integer, intent(in) :: m, d, theiler
    real(dp), intent(out) :: value
    integer, intent(out) :: status
    integer :: nembed, i, j
    integer(int64) :: count, denom_count
    real(dp) :: eps2

    nembed = embedding_length(size(series), m, d)
    if (nembed <= 0 .or. theiler < 0 .or. theiler >= nembed .or. eps <= 0.0_dp) then
      value = quiet_nan()
      status = 1
      return
    end if
    eps2 = eps * eps
    count = 0_int64
    do i = 1, nembed - theiler
      do j = i + theiler, nembed
        if (squared_embedding_distance(series, i, j, m, d) < eps2) count = count + 1_int64
      end do
    end do
    denom_count = int(nembed - theiler + 1, int64) * int(nembed - theiler, int64) / 2_int64
    if (denom_count <= 0_int64) then
      value = quiet_nan()
      status = 2
      return
    end if
    value = real(count, dp) / real(denom_count, dp)
    status = 0
  end subroutine correlation_integral

  subroutine correlation_dimension_curve(series, max_m, d, theiler, eps_min, neps, eps, c2, status)
    real(dp), intent(in) :: series(:), eps_min
    integer, intent(in) :: max_m, d, theiler, neps
    real(dp), allocatable, intent(out) :: eps(:), c2(:, :)
    integer, intent(out) :: status
    integer :: nembed, i, j, k, w
    integer(int64) :: denom_count
    real(dp) :: eps_max, log_step, dist2

    nembed = embedding_length(size(series), max_m, d)
    if (max_m < 1 .or. d < 1 .or. theiler < 0 .or. theiler >= nembed .or. eps_min <= 0.0_dp .or. neps < 2) then
      allocate(eps(0), c2(0, 0))
      status = 1
      return
    end if
    eps_max = (maxval(series) - minval(series)) * sqrt(real(max_m, dp))
    if (eps_max <= eps_min) then
      allocate(eps(0), c2(0, 0))
      status = 2
      return
    end if
    allocate(eps(neps), c2(neps, max_m))
    log_step = log(eps_max / eps_min) / real(neps - 1, dp)
    do k = 1, neps
      eps(k) = eps_min * exp(real(k - 1, dp) * log_step)
    end do
    c2 = 0.0_dp
    do i = 1, nembed - theiler
      do j = i + theiler, nembed
        dist2 = 0.0_dp
        do w = 1, max_m
          dist2 = dist2 + (series(i + (w - 1) * d) - series(j + (w - 1) * d))**2
          do k = 1, neps
            if (dist2 < eps(k)**2) c2(k, w) = c2(k, w) + 1.0_dp
          end do
        end do
      end do
    end do
    denom_count = int(nembed - theiler + 1, int64) * int(nembed - theiler, int64) / 2_int64
    c2 = c2 / real(denom_count, dp)
    status = 0
  end subroutine correlation_dimension_curve

  subroutine average_mutual_information(series, partitions, lag_max, ami, status)
    real(dp), intent(in) :: series(:)
    integer, intent(in) :: partitions, lag_max
    real(dp), allocatable, intent(out) :: ami(:)
    integer, intent(out) :: status
    real(dp), allocatable :: x(:), joint(:, :), marginal(:)
    real(dp) :: span, total
    integer :: scale_status, lag, i, bx, by

    if (partitions < 1 .or. lag_max < 0 .or. size(series) <= lag_max) then
      allocate(ami(0))
      status = 1
      return
    end if
    call scale_unit_interval(series, x, span, scale_status)
    if (scale_status /= 0) then
      allocate(ami(0))
      status = 2
      return
    end if
    allocate(ami(0:lag_max), joint(partitions, partitions), marginal(partitions))
    do lag = 0, lag_max
      joint = 0.0_dp
      do i = 1, size(x) - lag
        bx = min(int(x(i) * real(partitions, dp)) + 1, partitions)
        by = min(int(x(i + lag) * real(partitions, dp)) + 1, partitions)
        joint(bx, by) = joint(bx, by) + 1.0_dp
      end do
      total = sum(joint)
      joint = joint / total
      marginal = sum(joint, dim=2)
      ami(lag) = entropy_term(joint) - 2.0_dp * entropy_term(marginal)
    end do
    status = 0
  contains
    pure real(dp) function entropy_term(p) result(value)
      real(dp), intent(in) :: p(..)
      integer :: ii, jj
      value = 0.0_dp
      select rank (p)
      rank (1)
        do ii = 1, size(p)
          if (p(ii) > 0.0_dp) value = value + p(ii) * log(p(ii))
        end do
      rank (2)
        do jj = 1, size(p, 2)
          do ii = 1, size(p, 1)
            if (p(ii, jj) > 0.0_dp) value = value + p(ii, jj) * log(p(ii, jj))
          end do
        end do
      rank default
        value = quiet_nan()
      end select
    end function entropy_term
  end subroutine average_mutual_information

  subroutine recurrence_distance_matrix(series, m, d, distances, status, normalize)
    real(dp), intent(in) :: series(:)
    integer, intent(in) :: m, d
    real(dp), allocatable, intent(out) :: distances(:, :)
    integer, intent(out) :: status
    logical, intent(in), optional :: normalize
    real(dp), allocatable :: embedded(:, :)
    real(dp) :: max_distance
    logical :: do_normalize
    integer :: i, j, emb_status

    call delay_embed(series, m, d, embedded, emb_status)
    if (emb_status /= 0) then
      allocate(distances(0, 0))
      status = 1
      return
    end if
    allocate(distances(size(embedded, 1), size(embedded, 1)))
    do_normalize = .true.
    if (present(normalize)) do_normalize = normalize
    do i = 1, size(embedded, 1)
      distances(i, i) = 0.0_dp
      do j = i + 1, size(embedded, 1)
        distances(i, j) = sqrt(sum((embedded(i, :) - embedded(j, :))**2))
        distances(j, i) = distances(i, j)
      end do
    end do
    if (do_normalize) then
      max_distance = maxval(distances)
      if (max_distance > 0.0_dp) distances = distances / max_distance
    end if
    status = 0
  end subroutine recurrence_distance_matrix

  subroutine space_time_separation(series, m, d, idt, steps, isolines, status)
    real(dp), intent(in) :: series(:)
    integer, intent(in) :: m, d, idt, steps
    real(dp), allocatable, intent(out) :: isolines(:, :)
    integer, intent(out) :: status
    integer, parameter :: neps = 1000, nfrac = 10
    integer :: nembed, istep, j, ibin, frac, cumulative, need, npair
    integer, allocatable :: hist(:)
    real(dp) :: eps_max2, dist2

    nembed = embedding_length(size(series), m, d)
    if (nembed <= 0 .or. idt < 1 .or. steps < 1 .or. (steps - 1) * idt >= nembed) then
      allocate(isolines(0, 0))
      status = 1
      return
    end if
    eps_max2 = ((maxval(series) - minval(series)) * sqrt(real(m, dp)))**2
    if (eps_max2 <= 0.0_dp) then
      allocate(isolines(0, 0))
      status = 2
      return
    end if
    allocate(isolines(nfrac, steps), hist(neps))
    do istep = 0, steps - 1
      hist = 0
      npair = nembed - istep * idt
      do j = 1, npair
        dist2 = squared_embedding_distance(series, j, j + istep * idt, m, d)
        ibin = min(int(dist2 * real(neps, dp) / eps_max2) + 1, neps)
        hist(ibin) = hist(ibin) + 1
      end do
      do frac = 1, nfrac
        need = ceiling(real(npair * frac, dp) / real(nfrac, dp))
        cumulative = 0
        ibin = 0
        do while (ibin < neps .and. cumulative < need)
          ibin = ibin + 1
          cumulative = cumulative + hist(ibin)
        end do
        isolines(frac, istep + 1) = sqrt(real(ibin, dp) * eps_max2 / real(neps, dp))
      end do
    end do
    status = 0
  end subroutine space_time_separation
end module chaos_metrics
