! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Derived from the R package strucchange 1.6-0. See NOTICE.md and UPSTREAM.md.
module strucchange_fstats
   use r_kinds, only : dp
   use strucchange_regression, only : ols_fit
   implicit none
   private
   public :: fstats_result
   public :: compute_fstats
   public :: ave_f_statistic
   public :: exp_f_statistic
   public :: sup_f_statistic

   type :: fstats_result
      integer :: nreg = 0
      integer :: nobs = 0
      integer :: breakpoint = 0
      integer :: info = 0
      real(dp) :: lambda = 0.0_dp
      real(dp) :: rss = 0.0_dp
      integer, allocatable :: points(:)
      real(dp), allocatable :: statistics(:)
   end type fstats_result
contains
   subroutine compute_fstats(x, y, result, from_index, to_index)
      real(dp), intent(in) :: x(:, :), y(:)
      type(fstats_result), intent(out) :: result
      integer, intent(in), optional :: from_index, to_index
      real(dp), allocatable :: beta(:), e(:), beta1(:), beta2(:)
      real(dp), allocatable :: e1(:), e2(:)
      real(dp) :: full_rss, segment_rss, sigma2, rss1, rss2
      integer :: first, last, i, k, n, np, rank, rank1, rank2, info

      n = size(x, 1)
      k = size(x, 2)
      result%nobs = n
      result%nreg = k
      if (size(y) /= n .or. n <= 2 * k) then
         result%info = -1
         return
      end if
      first = floor(0.15_dp * real(n, dp))
      if (present(from_index)) first = from_index
      first = max(first, k + 1)
      last = n - first
      if (present(to_index)) last = to_index
      last = min(last, n - k - 1)
      if (first > last) then
         result%info = -2
         return
      end if

      call ols_fit(x, y, beta, e, full_rss, rank, info)
      if (info /= 0 .or. rank < k) then
         result%info = -3
         return
      end if
      np = last - first + 1
      allocate(result%points(np), result%statistics(np))
      do i = 1, np
         result%points(i) = first + i - 1
         call ols_fit(x(1:result%points(i), :), y(1:result%points(i)), &
            beta1, e1, rss1, rank1, info)
         if (info /= 0 .or. rank1 < k) then
            result%info = -4
            return
         end if
         call ols_fit(x(result%points(i) + 1:n, :), &
            y(result%points(i) + 1:n), beta2, e2, rss2, rank2, info)
         if (info /= 0 .or. rank2 < k) then
            result%info = -5
            return
         end if
         segment_rss = rss1 + rss2
         sigma2 = segment_rss / real(n - 2 * k, dp)
         if (sigma2 <= 0.0_dp) then
            result%statistics(i) = huge(1.0_dp)
         else
            result%statistics(i) = (full_rss - segment_rss) / sigma2
         end if
      end do
      i = maxloc(result%statistics, dim = 1)
      result%breakpoint = result%points(i)
      result%rss = full_rss / (1.0_dp + maxval(result%statistics) / &
         real(n - 2 * k, dp))
      result%lambda = real((n - first) * last, dp) / real(first * (n - last), dp)
      result%info = 0
   end subroutine compute_fstats

   pure real(dp) function sup_f_statistic(result) result(value)
      type(fstats_result), intent(in) :: result
      if (allocated(result%statistics)) then
         value = maxval(result%statistics)
      else
         value = 0.0_dp
      end if
   end function sup_f_statistic

   pure real(dp) function ave_f_statistic(result) result(value)
      type(fstats_result), intent(in) :: result
      if (allocated(result%statistics) .and. size(result%statistics) > 0) then
         value = sum(result%statistics) / real(size(result%statistics), dp)
      else
         value = 0.0_dp
      end if
   end function ave_f_statistic

   pure real(dp) function exp_f_statistic(result) result(value)
      type(fstats_result), intent(in) :: result
      real(dp) :: offset
      if (.not. allocated(result%statistics) .or. &
          size(result%statistics) == 0) then
         value = 0.0_dp
         return
      end if
      offset = 0.5_dp * maxval(result%statistics)
      value = offset + log(sum(exp(0.5_dp * result%statistics - offset)) / &
         real(size(result%statistics), dp))
   end function exp_f_statistic
end module strucchange_fstats
