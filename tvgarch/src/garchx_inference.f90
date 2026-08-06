! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of garchx.
! Copyright (C) 2026 translation contributors.
! Original garchx package copyright (C) Genaro Sucarrat.
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 2 of the License, or
! (at your option) any later version.
module garchx_inference
   use garchx_kinds, only : dp
   use garchx_math, only : normal_cdf, student_t_quantile, rmnorm, empirical_quantile
   use garchx_linalg, only : invert_matrix, matrix_rank
   implicit none
   private
   public :: confidence_intervals, boundary_t_tests, boundary_wald_test
contains
   subroutine confidence_intervals(par, vcov, nobs, level, intervals, status)
      real(dp), intent(in) :: par(:), vcov(:, :), level
      integer, intent(in) :: nobs
      real(dp), allocatable, intent(out) :: intervals(:, :)
      integer, intent(out) :: status
      integer :: k, i, df
      real(dp) :: alpha, critical, se
      k = size(par)
      if (size(vcov, 1) /= k .or. size(vcov, 2) /= k .or. &
          level <= 0.0_dp .or. level >= 1.0_dp .or. nobs <= k) then
         status = 1
         allocate(intervals(0, 0))
         return
      end if
      df = nobs-k
      alpha = 0.5_dp*(1.0_dp-level)
      critical = student_t_quantile(1.0_dp-alpha, real(df, dp))
      allocate(intervals(k, 2))
      do i = 1, k
         se = sqrt(max(0.0_dp, vcov(i, i)))
         intervals(i, 1) = par(i)-critical*se
         intervals(i, 2) = par(i)+critical*se
      end do
      status = 0
   end subroutine confidence_intervals

   subroutine boundary_t_tests(par, vcov, indices, table, status)
      real(dp), intent(in) :: par(:), vcov(:, :)
      integer, intent(in) :: indices(:)
      real(dp), allocatable, intent(out) :: table(:, :)
      integer, intent(out) :: status
      integer :: i, idx
      real(dp) :: se, statistic
      if (size(vcov, 1) /= size(par) .or. size(vcov, 2) /= size(par) .or. &
          any(indices < 1) .or. any(indices > size(par))) then
         status = 1
         allocate(table(0, 0))
         return
      end if
      allocate(table(size(indices), 4))
      do i = 1, size(indices)
         idx = indices(i)
         se = sqrt(max(0.0_dp, vcov(idx, idx)))
         if (se <= tiny(1.0_dp)) then
            statistic = huge(1.0_dp)
         else
            statistic = par(idx)/se
         end if
         table(i, 1) = par(idx)
         table(i, 2) = se
         table(i, 3) = statistic
         table(i, 4) = 1.0_dp-normal_cdf(statistic)
      end do
      status = 0
   end subroutine boundary_t_tests

   subroutine boundary_wald_test(par, vcov, nobs, r_matrix, restrictions, levels, &
                                 n_sim, statistic, critical_values, status)
      real(dp), intent(in) :: par(:), vcov(:, :), r_matrix(:, :), restrictions(:), levels(:)
      integer, intent(in) :: nobs, n_sim
      real(dp), intent(out) :: statistic
      real(dp), allocatable, intent(out) :: critical_values(:)
      integer, intent(out) :: status
      integer :: q, k, i, j, inv_status, draw_status, rank_value
      real(dp), allocatable :: sigma(:, :), middle(:, :), middle_inv(:, :), diff(:)
      real(dp), allocatable :: draws(:, :), transformed(:), simulated_stats(:)

      k = size(par)
      q = size(r_matrix, 1)
      rank_value = matrix_rank(r_matrix)
      if (size(vcov, 1) /= k .or. size(vcov, 2) /= k .or. &
          size(r_matrix, 2) /= k .or. size(restrictions) /= q .or. &
          any(levels <= 0.0_dp) .or. any(levels >= 1.0_dp) .or. &
          nobs < 1 .or. n_sim < 1 .or. rank_value /= q) then
         status = 1
         statistic = 0.0_dp
         allocate(critical_values(0))
         return
      end if
      allocate(diff(q), middle(q, q))
      diff = matmul(r_matrix, par)-restrictions
      middle = matmul(r_matrix, matmul(vcov, transpose(r_matrix)))
      call invert_matrix(middle, middle_inv, inv_status)
      if (inv_status /= 0) then
         status = 2
         statistic = 0.0_dp
         allocate(critical_values(0))
         return
      end if
      statistic = dot_product(diff, matmul(middle_inv, diff))
      allocate(sigma(k, k))
      sigma = real(nobs, dp)*vcov
      middle = matmul(r_matrix, matmul(sigma, transpose(r_matrix)))
      call invert_matrix(middle, middle_inv, inv_status)
      if (inv_status /= 0) then
         status = 3
         allocate(critical_values(0))
         return
      end if
      call rmnorm(n_sim, vcov=sigma, draws=draws, status=draw_status)
      if (draw_status /= 0) then
         status = 4
         allocate(critical_values(0))
         return
      end if
      allocate(transformed(q), simulated_stats(n_sim), critical_values(size(levels)))
      do i = 1, n_sim
         transformed = matmul(r_matrix, draws(i, :))
         simulated_stats(i) = dot_product(transformed, matmul(middle_inv, transformed))
      end do
      do j = 1, size(levels)
         critical_values(j) = empirical_quantile(simulated_stats, 1.0_dp-levels(j))
      end do
      status = 0
   end subroutine boundary_wald_test
end module garchx_inference
