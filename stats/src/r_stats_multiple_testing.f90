! SPDX-License-Identifier: MIT
module r_stats_multiple_testing
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
   use r_kinds, only: dp
   use r_ordering, only: r_order
   use r_optional, only: optval
   implicit none
   private
   public :: p_adjust

contains

   pure function p_adjust(p, method, n) result(adjusted)
      !! Adjusts p-values for multiple testing, preserving input order and NaN positions.
      real(dp), intent(in) :: p(:) !! P-values in `[0,1]`, optionally containing NaNs.
      character(len=*), intent(in), optional :: method !! Holm (default), hochberg, hommel, bonferroni, BH, BY, fdr, or none.
      integer, intent(in), optional :: n !! Total hypothesis count; defaults to the nonmissing input count.
      real(dp), allocatable :: adjusted(:) !! Adjusted probabilities with the same shape as `p`.
      real(dp), allocatable :: observed(:), corrected(:)
      integer, allocatable :: order(:), positions(:)
      character(len=:), allocatable :: selected
      real(dp) :: bound, harmonic, candidate
      integer :: i, k, m, total

      selected = 'holm'
      if (present(method)) selected = trim(method)
      select case (selected)
      case ('holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none')
      case default
         error stop 'p_adjust: unsupported method'
      end select
      m = count(.not. ieee_is_nan(p))
      total = optval(n, m)
      if (total < m) error stop 'p_adjust: n is smaller than the number of observed p-values'
      allocate (observed(m), positions(m), corrected(m))
      k = 0
      do i = 1, size(p)
         if (ieee_is_nan(p(i))) cycle
         if (.not. ieee_is_finite(p(i))) error stop 'p_adjust: infinite p-value'
         if (p(i) < 0.0_dp .or. p(i) > 1.0_dp) error stop 'p_adjust: p-value outside [0,1]'
         k = k + 1
         observed(k) = p(i)
         positions(k) = i
      end do
      adjusted = p
      if (m == 0 .or. total <= 1 .or. selected == 'none') return
      call r_order(observed, order)
      harmonic = 1.0_dp
      if (selected == 'BY') then
         harmonic = 0.0_dp
         do i = 1, total
            harmonic = harmonic + 1.0_dp/real(i, dp)
         end do
      end if
      select case (selected)
      case ('hommel')
         call hommel_adjust(observed, order, total, corrected)
      case ('bonferroni')
         corrected = min(1.0_dp, real(total, dp)*observed)
      case ('holm')
         bound = 0.0_dp
         do i = 1, m
            k = order(i)
            bound = max(bound, real(total - i + 1, dp)*observed(k))
            corrected(k) = min(1.0_dp, bound)
         end do
      case default
         bound = 1.0_dp
         do i = m, 1, -1
            k = order(i)
            if (selected == 'hochberg') then
               candidate = real(total - i + 1, dp)*observed(k)
            else
               candidate = harmonic*real(total, dp)/real(i, dp)*observed(k)
            end if
            bound = min(bound, candidate)
            corrected(k) = bound
         end do
      end select
      adjusted(positions) = corrected
   end function p_adjust

   pure subroutine hommel_adjust(p, order, total, corrected)
      !! Applies closed Simes testing by considering the worst subset of each size containing a hypothesis.
      real(dp), intent(in) :: p(:) !! Observed probabilities in their original order.
      integer, intent(in) :: order(:) !! Indices sorting `p` in ascending order.
      integer, intent(in) :: total !! Total number of hypotheses, including unobserved probabilities of one.
      real(dp), intent(out) :: corrected(:) !! Hommel probabilities in the original order.
      real(dp), allocatable :: sorted(:)
      real(dp) :: tail_bound, subset_bound
      integer :: subset_size, boundary, rank, m

      m = size(p)
      allocate (sorted(total), source=1.0_dp)
      sorted(:m) = p(order)
      corrected = p
      ! For fixed size, the least significant subset containing a given hypothesis
      ! consists of that hypothesis and the largest remaining p-values.
      do subset_size = 2, total
         boundary = total - subset_size + 1
         tail_bound = 1.0_dp
         do rank = boundary + 1, total
            tail_bound = min(tail_bound, real(subset_size, dp)*sorted(rank)/ &
                             real(rank - boundary + 1, dp))
         end do
         subset_bound = min(tail_bound, real(subset_size, dp)*sorted(boundary))
         do rank = 1, m
            if (rank < boundary) then
               corrected(order(rank)) = max(corrected(order(rank)), &
                  min(tail_bound, real(subset_size, dp)*sorted(rank)))
            else
               corrected(order(rank)) = max(corrected(order(rank)), subset_bound)
            end if
         end do
      end do
   end subroutine hommel_adjust

end module r_stats_multiple_testing
