! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Derived from the R package strucchange 1.6-0. See NOTICE.md and UPSTREAM.md.
module strucchange_breakpoints
   use r_kinds, only : dp
   use strucchange_recresid, only : recursive_residuals
   use strucchange_regression, only : ols_fit
   implicit none
   private
   public :: breakpoint_result
   public :: breakpoint_path_result
   public :: breakpoint_confint_result
   public :: compute_breakpoints
   public :: compute_breakpoint_path
   public :: segmented_fit
   public :: breakpoint_confidence_intervals
   public :: best_break_count

   type :: breakpoint_result
      integer :: nobs = 0
      integer :: nreg = 0
      integer :: n_breaks = 0
      integer :: info = 0
      real(dp) :: rss = 0.0_dp
      real(dp) :: bic = 0.0_dp
      integer, allocatable :: breakpoints(:)
   end type breakpoint_result


   type :: breakpoint_confint_result
      integer :: n_breaks = 0
      integer :: info = 0
      integer, allocatable :: intervals(:, :)
      logical, allocatable :: valid(:)
   end type breakpoint_confint_result

   type :: breakpoint_path_result
      integer :: nobs = 0
      integer :: nreg = 0
      integer :: max_breaks = 0
      integer :: info = 0
      real(dp), allocatable :: rss(:)
      real(dp), allocatable :: bic(:)
      integer, allocatable :: breakpoints(:, :)
      integer, allocatable :: break_count(:)
   end type breakpoint_path_result
contains
   subroutine build_rss_triangle(x, y, h, triangle, info)
      real(dp), intent(in) :: x(:, :), y(:)
      integer, intent(in) :: h
      real(dp), allocatable, intent(out) :: triangle(:, :)
      integer, intent(out) :: info
      real(dp), allocatable :: rr(:)
      real(dp) :: cumulative
      integer :: i, j, k, n, ierr, rindex

      n = size(x, 1)
      k = size(x, 2)
      allocate(triangle(n, n))
      triangle = huge(1.0_dp)
      if (size(y) /= n .or. h <= k .or. 2 * h > n) then
         info = -1
         return
      end if
      do i = 1, n - h + 1
         if (n - i + 1 < k + 1) exit
         call recursive_residuals(x(i:n, :), y(i:n), rr, ierr)
         if (ierr /= 0) then
            info = -2
            return
         end if
         cumulative = 0.0_dp
         do rindex = 1, size(rr)
            j = i + k + rindex - 1
            cumulative = cumulative + rr(rindex) ** 2
            if (j - i + 1 >= h) triangle(i, j) = cumulative
         end do
      end do
      info = 0
   end subroutine build_rss_triangle

   subroutine compute_breakpoints(x, y, h, n_breaks, result)
      real(dp), intent(in) :: x(:, :), y(:)
      integer, intent(in) :: h, n_breaks
      type(breakpoint_result), intent(out) :: result
      real(dp), allocatable :: triangle(:, :), cost(:, :)
      integer, allocatable :: previous(:, :)
      real(dp) :: candidate
      integer :: j, n, k, p, s, segments, ierr

      n = size(x, 1)
      k = size(x, 2)
      result%nobs = n
      result%nreg = k
      result%n_breaks = n_breaks
      if (n_breaks < 0 .or. (n_breaks + 1) * h > n) then
         result%info = -1
         return
      end if
      call build_rss_triangle(x, y, h, triangle, ierr)
      if (ierr /= 0) then
         result%info = ierr
         return
      end if
      segments = n_breaks + 1
      allocate(cost(segments, n), previous(segments, n))
      cost = huge(1.0_dp)
      previous = 0
      do j = h, n
         cost(1, j) = triangle(1, j)
      end do
      do s = 2, segments
         do j = s * h, n
            do p = (s - 1) * h, j - h
               if (cost(s - 1, p) >= huge(1.0_dp) / 4.0_dp) cycle
               if (triangle(p + 1, j) >= huge(1.0_dp) / 4.0_dp) cycle
               candidate = cost(s - 1, p) + triangle(p + 1, j)
               if (candidate < cost(s, j)) then
                  cost(s, j) = candidate
                  previous(s, j) = p
               end if
            end do
         end do
      end do
      if (cost(segments, n) >= huge(1.0_dp) / 4.0_dp) then
         result%info = -3
         return
      end if
      result%rss = cost(segments, n)
      result%bic = breakpoint_bic(n, k, n_breaks, result%rss)
      allocate(result%breakpoints(n_breaks))
      j = n
      do s = segments, 2, -1
         p = previous(s, j)
         result%breakpoints(s - 1) = p
         j = p
      end do
      result%info = 0
   end subroutine compute_breakpoints

   subroutine compute_breakpoint_path(x, y, h, max_breaks, result)
      real(dp), intent(in) :: x(:, :), y(:)
      integer, intent(in) :: h, max_breaks
      type(breakpoint_path_result), intent(out) :: result
      type(breakpoint_result) :: one
      real(dp), allocatable :: beta(:), residuals(:)
      real(dp) :: rss0
      integer :: b, info, rank, n, k

      n = size(x, 1)
      k = size(x, 2)
      result%nobs = n
      result%nreg = k
      result%max_breaks = max_breaks
      if (max_breaks < 0) then
         result%info = -1
         return
      end if
      allocate(result%rss(0:max_breaks), result%bic(0:max_breaks))
      allocate(result%break_count(0:max_breaks))
      allocate(result%breakpoints(max(1, max_breaks), 0:max_breaks))
      result%breakpoints = 0
      result%break_count = 0

      call ols_fit(x, y, beta, residuals, rss0, rank, info)
      if (info /= 0 .or. rank < k) then
         result%info = -2
         return
      end if
      result%rss(0) = rss0
      result%bic(0) = breakpoint_bic(n, k, 0, rss0)
      do b = 1, max_breaks
         call compute_breakpoints(x, y, h, b, one)
         if (one%info /= 0) then
            result%rss(b) = huge(1.0_dp)
            result%bic(b) = huge(1.0_dp)
         else
            result%rss(b) = one%rss
            result%bic(b) = one%bic
            result%break_count(b) = b
            result%breakpoints(1:b, b) = one%breakpoints
         end if
      end do
      result%info = 0
   end subroutine compute_breakpoint_path

   subroutine breakpoint_confidence_intervals(x, y, breakpoints, level, result, &
      heterogeneous_regressors, heterogeneous_errors)
      use strucchange_monitoring, only : pargmax_v
      real(dp), intent(in) :: x(:, :), y(:)
      integer, intent(in) :: breakpoints(:)
      real(dp), intent(in) :: level
      type(breakpoint_confint_result), intent(out) :: result
      logical, intent(in), optional :: heterogeneous_regressors, heterogeneous_errors
      real(dp), allocatable :: coefficients(:, :), fitted(:), residuals(:), segment_rss(:)
      real(dp), allocatable :: q_full(:, :), q1(:, :), q2(:, :)
      real(dp), allocatable :: delta(:)
      real(dp) :: a2, lower_root, upper_root, p0, phi1, phi2
      real(dp) :: qprod1, qprod2, sigma1, sigma2, sigma_full, xi
      integer, allocatable :: boundaries(:)
      integer :: first1, first2, i, k, last1, last2, n, ierr
      logical :: het_err, het_reg

      n = size(x, 1)
      k = size(x, 2)
      result%n_breaks = size(breakpoints)
      if (size(y) /= n .or. level <= 0.0_dp .or. level >= 1.0_dp .or. &
         size(breakpoints) < 1) then
         result%info = -1
         allocate(result%intervals(0, 3), result%valid(0))
         return
      end if
      het_reg = .true.
      het_err = .true.
      if (present(heterogeneous_regressors)) het_reg = heterogeneous_regressors
      if (present(heterogeneous_errors)) het_err = heterogeneous_errors

      call segmented_fit(x, y, breakpoints, coefficients, fitted, residuals, &
         segment_rss, ierr)
      if (ierr /= 0) then
         result%info = -2
         allocate(result%intervals(0, 3), result%valid(0))
         return
      end if
      allocate(result%intervals(size(breakpoints), 3), result%valid(size(breakpoints)))
      allocate(boundaries(0:size(breakpoints) + 1))
      boundaries(0) = 0
      boundaries(1:size(breakpoints)) = breakpoints
      boundaries(size(breakpoints) + 1) = n
      result%intervals = 0
      result%valid = .false.
      allocate(q_full(k, k), q1(k, k), q2(k, k), delta(k))
      q_full = matmul(transpose(x), x) / real(n, dp)
      sigma_full = sum(segment_rss) / real(n, dp)
      a2 = (1.0_dp - level) / 2.0_dp

      do i = 1, size(breakpoints)
         first1 = boundaries(i - 1) + 1
         last1 = boundaries(i)
         first2 = boundaries(i) + 1
         last2 = boundaries(i + 1)

         if (het_reg) then
            q1 = matmul(transpose(x(first1:last1, :)), x(first1:last1, :)) / &
               real(last1 - first1 + 1, dp)
            q2 = matmul(transpose(x(first2:last2, :)), x(first2:last2, :)) / &
               real(last2 - first2 + 1, dp)
         else
            q1 = q_full
            q2 = q_full
         end if
         if (het_err) then
            sigma1 = segment_rss(i) / real(last1 - first1 + 1, dp)
            sigma2 = segment_rss(i + 1) / real(last2 - first2 + 1, dp)
         else
            sigma1 = sigma_full
            sigma2 = sigma_full
         end if

         delta = coefficients(:, i + 1) - coefficients(:, i)
         qprod1 = dot_product(delta, matmul(q1, delta))
         qprod2 = dot_product(delta, matmul(q2, delta))
         if (qprod1 <= 0.0_dp .or. qprod2 <= 0.0_dp .or. &
            sigma1 <= 0.0_dp .or. sigma2 <= 0.0_dp) cycle
         if (het_reg) then
            xi = qprod2 / qprod1
         else
            xi = 1.0_dp
         end if
         phi1 = sqrt(sigma1)
         phi2 = sqrt(sigma2)
         p0 = pargmax_v(0.0_dp, xi, phi1, phi2)
         if (p0 < a2 .or. p0 > 1.0_dp - a2) cycle

         call pargmax_quantile(a2, xi, phi1, phi2, lower_root, ierr)
         if (ierr /= 0) cycle
         call pargmax_quantile(1.0_dp - a2, xi, phi1, phi2, upper_root, ierr)
         if (ierr /= 0) cycle
         lower_root = lower_root * phi1 * phi1 / qprod1
         upper_root = upper_root * phi1 * phi1 / qprod1
         result%intervals(i, 1) = breakpoints(i) - ceiling(upper_root)
         result%intervals(i, 2) = breakpoints(i)
         result%intervals(i, 3) = breakpoints(i) - floor(lower_root)
         result%valid(i) = .true.
      end do
      result%info = 0
   end subroutine breakpoint_confidence_intervals

   subroutine pargmax_quantile(probability, xi, phi1, phi2, root, info)
      use strucchange_monitoring, only : pargmax_v
      real(dp), intent(in) :: probability, xi, phi1, phi2
      real(dp), intent(out) :: root
      integer, intent(out) :: info
      real(dp) :: low, high, midpoint, p0
      integer :: iteration

      p0 = pargmax_v(0.0_dp, xi, phi1, phi2)
      if (probability < p0) then
         high = 0.0_dp
         low = -1000.0_dp
         do while (pargmax_v(low, xi, phi1, phi2) > probability)
            if (abs(low) > 1.0e12_dp) then
               info = -1
               return
            end if
            low = 2.0_dp * low
         end do
      else
         low = 0.0_dp
         high = 1000.0_dp
         do while (pargmax_v(high, xi, phi1, phi2) < probability)
            if (high > 1.0e12_dp) then
               info = -1
               return
            end if
            high = 2.0_dp * high
         end do
      end if
      do iteration = 1, 180
         midpoint = 0.5_dp * (low + high)
         if (pargmax_v(midpoint, xi, phi1, phi2) < probability) then
            low = midpoint
         else
            high = midpoint
         end if
         if (high - low <= 16.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(midpoint))) exit
      end do
      root = 0.5_dp * (low + high)
      info = 0
   end subroutine pargmax_quantile

   subroutine segmented_fit(x, y, breakpoints, coefficients, fitted, residuals, &
      segment_rss, info)
      real(dp), intent(in) :: x(:, :), y(:)
      integer, intent(in) :: breakpoints(:)
      real(dp), allocatable, intent(out) :: coefficients(:, :)
      real(dp), allocatable, intent(out) :: fitted(:), residuals(:), segment_rss(:)
      integer, intent(out) :: info
      real(dp), allocatable :: beta(:), one_residual(:)
      real(dp) :: rss
      integer :: first, i, k, last, n, rank, ierr, segments

      n = size(x, 1)
      k = size(x, 2)
      segments = size(breakpoints) + 1
      if (size(y) /= n .or. n < 1 .or. k < 1) then
         info = -1
         allocate(coefficients(0, 0), fitted(0), residuals(0), segment_rss(0))
         return
      end if
      if (size(breakpoints) > 0) then
         if (breakpoints(1) < 1 .or. breakpoints(size(breakpoints)) >= n) then
            info = -2
            allocate(coefficients(0, 0), fitted(0), residuals(0), segment_rss(0))
            return
         end if
         do i = 2, size(breakpoints)
            if (breakpoints(i) <= breakpoints(i - 1)) then
               info = -2
               allocate(coefficients(0, 0), fitted(0), residuals(0), segment_rss(0))
               return
            end if
         end do
      end if

      allocate(coefficients(k, segments), fitted(n), residuals(n), segment_rss(segments))
      first = 1
      do i = 1, segments
         if (i <= size(breakpoints)) then
            last = breakpoints(i)
         else
            last = n
         end if
         if (last - first + 1 < k) then
            info = -3
            return
         end if
         call ols_fit(x(first:last, :), y(first:last), beta, one_residual, rss, rank, ierr)
         if (ierr /= 0 .or. rank < k) then
            info = -4
            return
         end if
         coefficients(:, i) = beta
         residuals(first:last) = one_residual
         fitted(first:last) = y(first:last) - one_residual
         segment_rss(i) = rss
         first = last + 1
      end do
      info = 0
   end subroutine segmented_fit

   integer function best_break_count(path) result(n_breaks)
      type(breakpoint_path_result), intent(in) :: path
      integer :: b
      real(dp) :: best_bic

      n_breaks = -1
      if (.not. allocated(path%bic)) return
      if (size(path%bic) < 1) return
      best_bic = huge(1.0_dp)
      do b = lbound(path%bic, 1), ubound(path%bic, 1)
         if (path%bic(b) < best_bic) then
            best_bic = path%bic(b)
            n_breaks = b
         end if
      end do
   end function best_break_count

   pure real(dp) function breakpoint_bic(n, k, n_breaks, rss) result(value)
      integer, intent(in) :: n, k, n_breaks
      real(dp), intent(in) :: rss
      real(dp), parameter :: two_pi = 6.2831853071795864769_dp
      integer :: degrees

      if (rss <= 0.0_dp) then
         value = -huge(1.0_dp)
         return
      end if
      degrees = (k + 1) * (n_breaks + 1)
      value = real(n, dp) * (log(rss) + 1.0_dp - log(real(n, dp)) + &
         log(two_pi)) + log(real(n, dp)) * real(degrees, dp)
   end function breakpoint_bic
end module strucchange_breakpoints
