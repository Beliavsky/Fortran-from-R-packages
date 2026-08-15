! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module vgam_gaitd
   use vgam_kinds, only : dp
   use vgam_distributions, only : dpois_v, dnbinom_v
   implicit none
   private

   type, public :: gaitd_distribution_t
      integer :: min_support = 0
      real(dp), allocatable :: pmf(:)
      real(dp), allocatable :: cdf(:)
      real(dp) :: mean = 0.0_dp
      real(dp) :: variance = 0.0_dp
      integer :: status = 0
   contains
      procedure :: probability => gaitd_probability
      procedure :: cumulative => gaitd_cumulative
      procedure :: quantile => gaitd_quantile
      procedure :: random => gaitd_random
   end type gaitd_distribution_t

   public :: gaitd_transform_pmf, gaitd_poisson, gaitd_negative_binomial

contains

   subroutine gaitd_transform_pmf(base_pmf, min_support, dist, truncate, &
                                  altered_points, altered_probabilities, &
                                  inflated_points, inflation_probabilities, &
                                  deflated_points, deflation_factors)
      real(dp), intent(in) :: base_pmf(:)
      integer, intent(in) :: min_support
      type(gaitd_distribution_t), intent(out) :: dist
      integer, intent(in), optional :: truncate(:), altered_points(:)
      integer, intent(in), optional :: inflated_points(:), deflated_points(:)
      real(dp), intent(in), optional :: altered_probabilities(:)
      real(dp), intent(in), optional :: inflation_probabilities(:)
      real(dp), intent(in), optional :: deflation_factors(:)
      real(dp), allocatable :: work(:)
      real(dp) :: reserved, remain, total, x, ex2
      integer :: i, idx

      dist%min_support = min_support
      if (size(base_pmf) <= 0 .or. any(base_pmf < 0.0_dp)) then
         dist%status = 1
         return
      end if
      total = sum(base_pmf)
      if (total <= 0.0_dp) then
         dist%status = 2
         return
      end if
      allocate(work(size(base_pmf)))
      work = base_pmf/total

      if (present(truncate)) then
         do i = 1, size(truncate)
            idx = point_index(truncate(i), min_support, size(work))
            if (idx == 0) then
               dist%status = 3
               return
            end if
            work(idx) = 0.0_dp
         end do
      end if

      if (present(deflated_points) .neqv. present(deflation_factors)) then
         dist%status = 4
         return
      end if
      if (present(deflated_points)) then
         if (size(deflated_points) /= size(deflation_factors) .or. &
             any(deflation_factors < 0.0_dp) .or. any(deflation_factors > 1.0_dp)) then
            dist%status = 5
            return
         end if
         do i = 1, size(deflated_points)
            idx = point_index(deflated_points(i), min_support, size(work))
            if (idx == 0) then
               dist%status = 6
               return
            end if
            work(idx) = work(idx)*deflation_factors(i)
         end do
      end if

      reserved = 0.0_dp
      if (present(altered_points) .neqv. present(altered_probabilities)) then
         dist%status = 7
         return
      end if
      if (present(altered_points)) then
         if (size(altered_points) /= size(altered_probabilities) .or. &
             any(altered_probabilities < 0.0_dp) .or. any(altered_probabilities > 1.0_dp)) then
            dist%status = 8
            return
         end if
         reserved = reserved + sum(altered_probabilities)
         do i = 1, size(altered_points)
            idx = point_index(altered_points(i), min_support, size(work))
            if (idx == 0) then
               dist%status = 9
               return
            end if
            work(idx) = 0.0_dp
         end do
      end if

      if (present(inflated_points) .neqv. present(inflation_probabilities)) then
         dist%status = 10
         return
      end if
      if (present(inflated_points)) then
         if (size(inflated_points) /= size(inflation_probabilities) .or. &
             any(inflation_probabilities < 0.0_dp) .or. any(inflation_probabilities > 1.0_dp)) then
            dist%status = 11
            return
         end if
         reserved = reserved + sum(inflation_probabilities)
         do i = 1, size(inflated_points)
            idx = point_index(inflated_points(i), min_support, size(work))
            if (idx == 0) then
               dist%status = 12
               return
            end if
         end do
      end if

      if (reserved >= 1.0_dp - 100.0_dp*epsilon(1.0_dp)) then
         dist%status = 13
         return
      end if
      remain = sum(work)
      if (remain <= 0.0_dp) then
         dist%status = 14
         return
      end if
      work = work*(1.0_dp - reserved)/remain

      if (present(altered_points)) then
         do i = 1, size(altered_points)
            idx = point_index(altered_points(i), min_support, size(work))
            work(idx) = altered_probabilities(i)
         end do
      end if
      if (present(inflated_points)) then
         do i = 1, size(inflated_points)
            idx = point_index(inflated_points(i), min_support, size(work))
            work(idx) = work(idx) + inflation_probabilities(i)
         end do
      end if

      work = work/sum(work)
      dist%pmf = work
      allocate(dist%cdf(size(work)))
      dist%cdf(1) = work(1)
      do i = 2, size(work)
         dist%cdf(i) = min(1.0_dp, dist%cdf(i - 1) + work(i))
      end do
      dist%cdf(size(work)) = 1.0_dp
      dist%mean = 0.0_dp
      ex2 = 0.0_dp
      do i = 1, size(work)
         x = real(min_support + i - 1, dp)
         dist%mean = dist%mean + x*work(i)
         ex2 = ex2 + x*x*work(i)
      end do
      dist%variance = max(0.0_dp, ex2 - dist%mean*dist%mean)
      dist%status = 0
   end subroutine gaitd_transform_pmf

   subroutine gaitd_poisson(lambda, max_support, dist, truncate, altered_points, &
                            altered_probabilities, inflated_points, inflation_probabilities, &
                            deflated_points, deflation_factors)
      real(dp), intent(in) :: lambda
      integer, intent(in) :: max_support
      type(gaitd_distribution_t), intent(out) :: dist
      integer, intent(in), optional :: truncate(:), altered_points(:)
      integer, intent(in), optional :: inflated_points(:), deflated_points(:)
      real(dp), intent(in), optional :: altered_probabilities(:)
      real(dp), intent(in), optional :: inflation_probabilities(:), deflation_factors(:)
      real(dp), allocatable :: base(:)
      integer :: k

      if (lambda < 0.0_dp .or. max_support < 0) then
         dist%status = 20
         return
      end if
      allocate(base(max_support + 1))
      do k = 0, max_support
         base(k + 1) = dpois_v(k, lambda)
      end do
      call gaitd_transform_pmf(base, 0, dist, truncate, altered_points, &
                               altered_probabilities, inflated_points, inflation_probabilities, &
                               deflated_points, deflation_factors)
   end subroutine gaitd_poisson

   subroutine gaitd_negative_binomial(mu, size, max_support, dist, truncate, &
                                      altered_points, altered_probabilities, inflated_points, &
                                      inflation_probabilities, deflated_points, deflation_factors)
      real(dp), intent(in) :: mu, size
      integer, intent(in) :: max_support
      type(gaitd_distribution_t), intent(out) :: dist
      integer, intent(in), optional :: truncate(:), altered_points(:)
      integer, intent(in), optional :: inflated_points(:), deflated_points(:)
      real(dp), intent(in), optional :: altered_probabilities(:)
      real(dp), intent(in), optional :: inflation_probabilities(:), deflation_factors(:)
      real(dp), allocatable :: base(:)
      real(dp) :: prob
      integer :: k

      if (mu < 0.0_dp .or. size <= 0.0_dp .or. max_support < 0) then
         dist%status = 21
         return
      end if
      prob = size/(size + mu)
      allocate(base(max_support + 1))
      do k = 0, max_support
         base(k + 1) = dnbinom_v(k, size, prob)
      end do
      call gaitd_transform_pmf(base, 0, dist, truncate, altered_points, &
                               altered_probabilities, inflated_points, inflation_probabilities, &
                               deflated_points, deflation_factors)
   end subroutine gaitd_negative_binomial

   real(dp) function gaitd_probability(self, x) result(p)
      class(gaitd_distribution_t), intent(in) :: self
      integer, intent(in) :: x
      integer :: idx
      p = 0.0_dp
      if (.not. allocated(self%pmf)) return
      idx = point_index(x, self%min_support, size(self%pmf))
      if (idx > 0) p = self%pmf(idx)
   end function gaitd_probability

   real(dp) function gaitd_cumulative(self, x) result(p)
      class(gaitd_distribution_t), intent(in) :: self
      integer, intent(in) :: x
      integer :: idx
      if (.not. allocated(self%cdf)) then
         p = 0.0_dp
         return
      end if
      if (x < self%min_support) then
         p = 0.0_dp
      else if (x >= self%min_support + size(self%cdf) - 1) then
         p = 1.0_dp
      else
         idx = x - self%min_support + 1
         p = self%cdf(idx)
      end if
   end function gaitd_cumulative

   integer function gaitd_quantile(self, probability) result(q)
      class(gaitd_distribution_t), intent(in) :: self
      real(dp), intent(in) :: probability
      integer :: lo, hi, mid
      if (.not. allocated(self%cdf) .or. probability < 0.0_dp .or. probability > 1.0_dp) then
         q = huge(1)
         return
      end if
      if (probability <= 0.0_dp) then
         q = self%min_support
         return
      end if
      lo = 1
      hi = size(self%cdf)
      do while (lo < hi)
         mid = lo + (hi - lo)/2
         if (self%cdf(mid) >= probability) then
            hi = mid
         else
            lo = mid + 1
         end if
      end do
      q = self%min_support + lo - 1
   end function gaitd_quantile

   integer function gaitd_random(self) result(x)
      class(gaitd_distribution_t), intent(in) :: self
      real(dp) :: u
      call random_number(u)
      x = self%quantile(u)
   end function gaitd_random

   pure integer function point_index(point, min_support, n) result(idx)
      integer, intent(in) :: point, min_support, n
      idx = point - min_support + 1
      if (idx < 1 .or. idx > n) idx = 0
   end function point_index

end module vgam_gaitd
