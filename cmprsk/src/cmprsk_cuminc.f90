! Copyright (C) 2000 Robert Gray
! Modern Fortran translation maintained for Fortran-from-R-packages.
! SPDX-License-Identifier: GPL-2.0-or-later
module cmprsk_cuminc
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use r_kinds, only : dp
   use r_distributions, only : r_pchisq
   use r_linalg, only : solve_system
   use cmprsk_status, only : cmprsk_success, cmprsk_invalid_argument, cmprsk_singular_matrix
   implicit none
   private

   type, public :: cuminc_curve
      real(dp), allocatable :: time(:)
      real(dp), allocatable :: estimate(:)
      real(dp), allocatable :: variance(:)
   end type cuminc_curve

   type, public :: gray_test_result
      real(dp), allocatable :: score(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp) :: statistic = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      integer :: degrees_freedom = 0
   end type gray_test_result

   public :: cumulative_incidence
   public :: gray_test
   public :: timepoint_indices
   public :: curve_timepoints

contains

   pure subroutine cumulative_incidence(time, failed_any, failed_cause, curve, status)
      real(dp), intent(in) :: time(:) !! Follow-up times sorted in nondecreasing order; values must be nonnegative.
      integer, intent(in) :: failed_any(:) !! One for failure from any cause and zero for censoring.
      integer, intent(in) :: failed_cause(:) !! One for failure from the cause of interest and zero otherwise.
      type(cuminc_curve), intent(out) :: curve !! Step-function times, cumulative-incidence estimates, and variances.
      integer, intent(out) :: status !! `cmprsk_success` or `cmprsk_invalid_argument`.

      integer :: i
      integer :: l
      integer :: ll
      integer :: lcnt
      integer :: n
      integer :: nd
      integer :: nd1
      integer :: nd2
      integer :: nfc
      integer :: rs
      real(dp) :: fk
      real(dp) :: fkn
      real(dp) :: t2
      real(dp) :: t3
      real(dp) :: t4
      real(dp) :: t5
      real(dp) :: t6
      real(dp) :: ty
      real(dp) :: v1
      real(dp) :: v2
      real(dp) :: v3

      n = size(time)
      if (size(failed_any) /= n .or. size(failed_cause) /= n .or. n < 1) then
         status = cmprsk_invalid_argument
         return
      end if
      if (any(ieee_is_nan(time)) .or. any(time < 0.0_dp) .or. any(failed_any < 0) .or. any(failed_any > 1) .or. &
          any(failed_cause < 0) .or. any(failed_cause > 1) .or. &
          any(failed_cause > failed_any)) then
         status = cmprsk_invalid_argument
         return
      end if
      do i = 2, n
         if (time(i) < time(i - 1)) then
            status = cmprsk_invalid_argument
            return
         end if
      end do

      nfc = 0
      i = 1
      do while (i <= n)
         ty = time(i)
         nd1 = 0
         do while (i <= n)
            if (time(i) /= ty) exit
            nd1 = nd1 + failed_cause(i)
            i = i + 1
         end do
         if (nd1 > 0) nfc = nfc + 1
      end do

      allocate(curve%time(2*nfc + 2))
      allocate(curve%estimate(2*nfc + 2))
      allocate(curve%variance(2*nfc + 2))
      curve%time = 0.0_dp
      curve%estimate = 0.0_dp
      curve%variance = 0.0_dp

      fk = 1.0_dp
      v1 = 0.0_dp
      v2 = 0.0_dp
      v3 = 0.0_dp
      curve%time(1) = 0.0_dp
      curve%estimate(1) = 0.0_dp
      curve%variance(1) = 0.0_dp
      lcnt = 1
      l = 1
      ll = 1
      rs = n
      ty = time(1)

      do
         l = l + 1
         do while (l <= n)
            if (time(l) /= ty) exit
            l = l + 1
         end do
         l = l - 1

         nd1 = sum(failed_cause(ll:l))
         nd2 = sum(failed_any(ll:l) - failed_cause(ll:l))
         nd = nd1 + nd2
         if (nd > 0) then
            fkn = fk * real(rs - nd, dp) / real(rs, dp)
            if (nd1 > 0) then
               lcnt = lcnt + 2
               curve%estimate(lcnt - 1) = curve%estimate(lcnt - 2)
               curve%estimate(lcnt) = curve%estimate(lcnt - 1) + &
                    fk*real(nd1, dp)/real(rs, dp)
            end if

            if (nd2 > 0 .and. fkn > 0.0_dp) then
               t5 = 1.0_dp
               if (nd2 > 1) t5 = 1.0_dp - real(nd2 - 1, dp)/real(rs - 1, dp)
               t6 = fk*fk*t5*real(nd2, dp)/real(rs*rs, dp)
               t3 = 1.0_dp/fkn
               t4 = curve%estimate(lcnt)/fkn
               v1 = v1 + t4*t4*t6
               v2 = v2 + t3*t4*t6
               v3 = v3 + t3*t3*t6
            end if

            if (nd1 > 0) then
               t5 = 1.0_dp
               if (nd1 > 1) t5 = 1.0_dp - real(nd1 - 1, dp)/real(rs - 1, dp)
               t6 = fk*fk*t5*real(nd1, dp)/real(rs*rs, dp)
               t3 = 0.0_dp
               if (fkn > 0.0_dp) t3 = 1.0_dp/fkn
               t4 = 1.0_dp + t3*curve%estimate(lcnt)
               v1 = v1 + t4*t4*t6
               v2 = v2 + t3*t4*t6
               v3 = v3 + t3*t3*t6
               t2 = curve%estimate(lcnt)
               curve%time(lcnt - 1) = time(l)
               curve%time(lcnt) = time(l)
               curve%variance(lcnt - 1) = curve%variance(lcnt - 2)
               curve%variance(lcnt) = v1 + t2*t2*v3 - 2.0_dp*t2*v2
            end if
            fk = fkn
         end if

         rs = n - l
         l = l + 1
         if (l > n) exit
         ll = l
         ty = time(l)
      end do

      lcnt = lcnt + 1
      curve%time(lcnt) = time(n)
      curve%estimate(lcnt) = curve%estimate(lcnt - 1)
      curve%variance(lcnt) = curve%variance(lcnt - 1)
      status = cmprsk_success
   end subroutine cumulative_incidence

   pure subroutine gray_test(time, status_code, group, strata, rho, result, status)
      real(dp), intent(in) :: time(:) !! Follow-up times sorted in nondecreasing order.
      integer, intent(in) :: status_code(:) !! Zero for censoring, one for the target cause, and two for competing failures.
      integer, intent(in) :: group(:) !! Group codes in `1:ng`.
      integer, intent(in) :: strata(:) !! Stratum codes in `1:nst`.
      real(dp), intent(in) :: rho !! Gray-test weight exponent in `(1 - F(t-))**rho`.
      type(gray_test_result), intent(out) :: result !! Scores, covariance, chi-square statistic, p-value, and degrees of freedom.
      integer, intent(out) :: status !! Success, invalid-argument, or singular-covariance status.

      integer :: i
      integer :: info
      integer :: k
      integer :: n
      integer :: ng
      integer :: ng1
      integer :: nst
      integer :: ns
      integer, allocatable :: group_s(:)
      integer, allocatable :: status_s(:)
      real(dp), allocatable :: covariance_s(:, :)
      real(dp), allocatable :: rhs(:)
      real(dp), allocatable :: score_s(:)
      real(dp), allocatable :: solution(:)
      real(dp), allocatable :: time_s(:)

      n = size(time)
      if (size(status_code) /= n .or. size(group) /= n .or. size(strata) /= n .or. n < 1) then
         status = cmprsk_invalid_argument
         return
      end if
      if (any(ieee_is_nan(time)) .or. ieee_is_nan(rho) .or. any(status_code < 0) .or. &
          any(status_code > 2) .or. any(group < 1) .or. any(strata < 1)) then
         status = cmprsk_invalid_argument
         return
      end if
      do i = 2, n
         if (time(i) < time(i - 1)) then
            status = cmprsk_invalid_argument
            return
         end if
      end do

      ng = maxval(group)
      nst = maxval(strata)
      if (ng < 2) then
         status = cmprsk_invalid_argument
         return
      end if
      ng1 = ng - 1
      allocate(result%score(ng1), result%covariance(ng1, ng1))
      result%score = 0.0_dp
      result%covariance = 0.0_dp

      do k = 1, nst
         ns = count(strata == k)
         if (ns == 0) cycle
         allocate(time_s(ns), status_s(ns), group_s(ns), score_s(ng1), covariance_s(ng1, ng1))
         ns = 0
         do i = 1, n
            if (strata(i) /= k) cycle
            ns = ns + 1
            time_s(ns) = time(i)
            status_s(ns) = status_code(i)
            group_s(ns) = group(i)
         end do
         call gray_test_one_stratum(time_s, status_s, group_s, rho, score_s, covariance_s)
         result%score = result%score + score_s
         result%covariance = result%covariance + covariance_s
         deallocate(time_s, status_s, group_s, score_s, covariance_s)
      end do

      allocate(rhs(ng1), solution(ng1))
      rhs = result%score
      call solve_system(result%covariance, rhs, solution, info)
      if (info /= 0) then
         result%statistic = -1.0_dp
         result%p_value = 1.0_dp
         result%degrees_freedom = ng1
         status = cmprsk_singular_matrix
         return
      end if
      result%statistic = dot_product(result%score, solution)
      result%degrees_freedom = ng1
      result%p_value = r_pchisq(result%statistic, real(ng1, dp), lower_tail=.false.)
      status = cmprsk_success
   end subroutine gray_test

   pure subroutine gray_test_one_stratum(time, status_code, group, rho, score, covariance)
      real(dp), intent(in) :: time(:) !! Sorted follow-up times for one stratum.
      integer, intent(in) :: status_code(:) !! Zero for censoring, one for target failure, and two for competing failure.
      integer, intent(in) :: group(:) !! Group codes in `1:ng` for the stratum.
      real(dp), intent(in) :: rho !! Gray-test weight exponent.
      real(dp), intent(out) :: score(:) !! Scores for groups `1:(ng-1)` relative to the last group.
      real(dp), intent(out) :: covariance(:, :) !! Estimated covariance matrix of `score`.

      integer :: d(0:2, maxval(group))
      integer :: i
      integer :: j
      integer :: k
      integer :: ll
      integer :: lu
      integer :: n
      integer :: nd1
      integer :: nd2
      integer :: ng
      integer :: ng1
      integer :: rs(maxval(group))
      real(dp) :: a(maxval(group), maxval(group))
      real(dp) :: c(maxval(group), maxval(group))
      real(dp) :: f
      real(dp) :: f1(maxval(group))
      real(dp) :: f1m(maxval(group))
      real(dp) :: fb
      real(dp) :: fm
      real(dp) :: skm(maxval(group))
      real(dp) :: skmm(maxval(group))
      real(dp) :: t1
      real(dp) :: t2
      real(dp) :: t3
      real(dp) :: t4
      real(dp) :: t5
      real(dp) :: t6
      real(dp) :: td
      real(dp) :: tq
      real(dp) :: tr
      real(dp) :: v2(size(score), maxval(group))
      real(dp) :: v3(maxval(group))

      n = size(time)
      ng = maxval(group)
      ng1 = ng - 1
      score = 0.0_dp
      covariance = 0.0_dp
      rs = 0
      do i = 1, n
         rs(group(i)) = rs(group(i)) + 1
      end do
      f1m = 0.0_dp
      f1 = 0.0_dp
      skmm = 1.0_dp
      skm = 1.0_dp
      v3 = 0.0_dp
      v2 = 0.0_dp
      c = 0.0_dp
      fm = 0.0_dp
      f = 0.0_dp
      ll = 1

      do
         lu = ll
         do while (lu < n)
            if (time(lu + 1) /= time(ll)) exit
            lu = lu + 1
         end do
         d = 0
         do i = ll, lu
            d(status_code(i), group(i)) = d(status_code(i), group(i)) + 1
         end do
         nd1 = sum(d(1, :))
         nd2 = sum(d(2, :))

         if (nd1 /= 0 .or. nd2 /= 0) then
            tr = 0.0_dp
            tq = 0.0_dp
            do i = 1, ng
               if (rs(i) <= 0) cycle
               td = real(d(1, i) + d(2, i), dp)
               skm(i) = skmm(i)*(real(rs(i), dp) - td)/real(rs(i), dp)
               f1(i) = f1m(i) + skmm(i)*real(d(1, i), dp)/real(rs(i), dp)
               tr = tr + real(rs(i), dp)/skmm(i)
               tq = tq + real(rs(i), dp)*(1.0_dp - f1m(i))/skmm(i)
            end do
            f = fm + real(nd1, dp)/tr
            fb = (1.0_dp - fm)**rho

            a = 0.0_dp
            do i = 1, ng
               if (rs(i) <= 0) cycle
               t1 = real(rs(i), dp)/skmm(i)
               a(i, i) = fb*t1*(1.0_dp - t1/tr)
               if (a(i, i) /= 0.0_dp) then
                  c(i, i) = c(i, i) + a(i, i)*real(nd1, dp)/(tr*(1.0_dp - fm))
               end if
               do j = i + 1, ng
                  if (rs(j) <= 0) cycle
                  a(i, j) = -fb*t1*real(rs(j), dp)/(skmm(j)*tr)
                  if (a(i, j) /= 0.0_dp) then
                     c(i, j) = c(i, j) + a(i, j)*real(nd1, dp)/(tr*(1.0_dp - fm))
                  end if
               end do
            end do
            do i = 2, ng
               do j = 1, i - 1
                  a(i, j) = a(j, i)
                  c(i, j) = c(j, i)
               end do
            end do

            do i = 1, ng1
               if (rs(i) <= 0) cycle
               score(i) = score(i) + fb*(real(d(1, i), dp) - &
                    real(nd1*rs(i), dp)*(1.0_dp - f1m(i))/(skmm(i)*tq))
            end do

            if (nd1 > 0) then
               do k = 1, ng
                  if (rs(k) <= 0) cycle
                  t4 = 1.0_dp
                  if (skm(k) > 0.0_dp) t4 = 1.0_dp - (1.0_dp - f)/skm(k)
                  t5 = 1.0_dp
                  if (nd1 > 1) then
                     t5 = 1.0_dp - real(nd1 - 1, dp)/(tr*skmm(k) - 1.0_dp)
                  end if
                  t3 = t5*skmm(k)*real(nd1, dp)/(tr*real(rs(k), dp))
                  v3(k) = v3(k) + t4*t4*t3
                  do i = 1, ng1
                     t1 = a(i, k) - t4*c(i, k)
                     v2(i, k) = v2(i, k) + t1*t4*t3
                     do j = 1, i
                        t2 = a(j, k) - t4*c(j, k)
                        covariance(i, j) = covariance(i, j) + t1*t2*t3
                     end do
                  end do
               end do
            end if

            if (nd2 > 0) then
               do k = 1, ng
                  if (skm(k) <= 0.0_dp .or. d(2, k) <= 0) cycle
                  t4 = (1.0_dp - f)/skm(k)
                  t5 = 1.0_dp
                  if (d(2, k) > 1) then
                     t5 = 1.0_dp - real(d(2, k) - 1, dp)/real(rs(k) - 1, dp)
                  end if
                  t6 = real(rs(k), dp)
                  t3 = t5*(skmm(k)**2)*real(d(2, k), dp)/(t6**2)
                  v3(k) = v3(k) + t4*t4*t3
                  do i = 1, ng1
                     t1 = t4*c(i, k)
                     v2(i, k) = v2(i, k) - t1*t4*t3
                     do j = 1, i
                        t2 = t4*c(j, k)
                        covariance(i, j) = covariance(i, j) + t1*t2*t3
                     end do
                  end do
               end do
            end if
         end if

         if (lu >= n) exit
         do i = ll, lu
            rs(group(i)) = rs(group(i)) - 1
         end do
         fm = f
         f1m = f1
         skmm = skm
         ll = lu + 1
      end do

      do i = 1, ng1
         do j = 1, i
            do k = 1, ng
               covariance(i, j) = covariance(i, j) + c(i, k)*c(j, k)*v3(k)
               covariance(i, j) = covariance(i, j) + c(i, k)*v2(j, k)
               covariance(i, j) = covariance(i, j) + c(j, k)*v2(i, k)
            end do
            covariance(j, i) = covariance(i, j)
         end do
      end do
   end subroutine gray_test_one_stratum

   pure subroutine timepoint_indices(curve_time, requested_time, index)
      real(dp), intent(in) :: curve_time(:) !! Step-function time coordinates in ascending order.
      real(dp), intent(in) :: requested_time(:) !! Requested times, sorted in ascending order.
      integer, intent(out) :: index(:) !! One-based curve indices, or zero when a requested time is beyond support.

      integer :: i
      integer :: k
      integer :: l
      integer :: n
      integer :: ntp

      n = size(curve_time)
      ntp = size(requested_time)
      index = 0
      if (n == 0 .or. ntp == 0) return
      l = ntp
      do i = ntp, 1, -1
         if (curve_time(n) >= requested_time(i)) exit
         index(l) = 0
         l = l - 1
      end do
      if (l <= 0) return
      if (curve_time(n) == requested_time(l)) then
         index(l) = n
         l = l - 1
      end if
      k = n - 1
      do while (l > 0)
         do i = k, 1, -1
            if (curve_time(k) <= requested_time(l)) then
               index(l) = k + 1
               l = l - 1
               exit
            else
               k = k - 1
            end if
         end do
         if (k < 1 .and. l > 0) then
            index(1:l) = 0
            exit
         end if
      end do
   end subroutine timepoint_indices

   pure subroutine curve_timepoints(curve, requested_time, estimate, variance, present)
      type(cuminc_curve), intent(in) :: curve !! Cumulative-incidence curve returned by `cumulative_incidence`.
      real(dp), intent(in) :: requested_time(:) !! Requested times in ascending order.
      real(dp), intent(out) :: estimate(:) !! Estimates at requested times; undefined where `present` is false.
      real(dp), intent(out) :: variance(:) !! Variances at requested times; undefined where `present` is false.
      logical, intent(out) :: present(:) !! True where the curve extends to the corresponding requested time.

      integer :: i
      integer :: index(size(requested_time))

      call timepoint_indices(curve%time, requested_time, index)
      estimate = 0.0_dp
      variance = 0.0_dp
      present = index > 0
      do i = 1, size(requested_time)
         if (.not. present(i)) cycle
         estimate(i) = curve%estimate(index(i))
         variance(i) = curve%variance(index(i))
      end do
   end subroutine curve_timepoints

end module cmprsk_cuminc
