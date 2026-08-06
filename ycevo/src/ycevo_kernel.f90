module ycevo_kernel
   use ycevo_kinds, only : dp
   implicit none
   private

   real(dp), parameter :: pi = acos(-1.0_dp)

   public :: epanechnikov, epanechnikov_quantile, epanechnikov_random
   public :: calc_epaker_weights, get_weights, range_nonzero, span_to_bandwidth

contains

   elemental real(dp) function epanechnikov(x) result(weight)
      real(dp), intent(in) :: x

      if (abs(x) <= 1.0_dp) then
         weight = 0.75_dp * (1.0_dp - x*x)
      else
         weight = 0.0_dp
      end if
   end function epanechnikov

   elemental real(dp) function epanechnikov_quantile(p, mu, radius) result(x)
      real(dp), intent(in) :: p, mu, radius
      real(dp) :: pp

      pp = min(max(p, 0.0_dp), 1.0_dp)
      x = 2.0_dp * sin(asin(2.0_dp*pp - 1.0_dp) / 3.0_dp) * radius + mu
   end function epanechnikov_quantile

   subroutine epanechnikov_random(x, mu, radius)
      real(dp), intent(out) :: x(:)
      real(dp), intent(in) :: mu, radius
      real(dp), allocatable :: u(:)

      allocate(u(size(x)))
      call random_number(u)
      x = epanechnikov_quantile(u, mu, radius)
   end subroutine epanechnikov_random

   subroutine calc_epaker_weights(gamma, grid, bandwidth, weights)
      real(dp), intent(in) :: gamma(:), grid(:), bandwidth(:)
      real(dp), allocatable, intent(out) :: weights(:, :)
      integer :: i, j

      if (size(bandwidth) /= 1 .and. size(bandwidth) /= size(grid)) then
         allocate(weights(0, 0))
         return
      end if

      allocate(weights(size(gamma), size(grid)))
      do j = 1, size(grid)
         do i = 1, size(gamma)
            if (size(bandwidth) == 1) then
               weights(i, j) = epanechnikov((grid(j) - gamma(i)) / bandwidth(1))
            else
               weights(i, j) = epanechnikov((grid(j) - gamma(i)) / bandwidth(j))
            end if
         end do
      end do
   end subroutine calc_epaker_weights

   subroutine get_weights(grid, bandwidth, n, weights, units)
      real(dp), intent(in) :: grid(:), bandwidth(:)
      integer, intent(in) :: n
      real(dp), allocatable, intent(out) :: weights(:, :)
      real(dp), intent(in), optional :: units
      real(dp), allocatable :: gamma(:)
      real(dp) :: scale
      integer :: i

      scale = real(n, dp)
      if (present(units)) scale = units
      allocate(gamma(n))
      do i = 1, n
         gamma(i) = real(i, dp) / scale
      end do
      call calc_epaker_weights(gamma, grid, bandwidth, weights)
   end subroutine get_weights

   subroutine range_nonzero(weights, first, last, threshold)
      real(dp), intent(in) :: weights(:, :)
      integer, allocatable, intent(out) :: first(:), last(:)
      real(dp), intent(in), optional :: threshold
      real(dp) :: cutoff
      integer :: i, j

      cutoff = 0.0_dp
      if (present(threshold)) cutoff = threshold
      allocate(first(size(weights, 2)), last(size(weights, 2)))
      first = 0
      last = 0
      do j = 1, size(weights, 2)
         do i = 1, size(weights, 1)
            if (weights(i, j) > cutoff) then
               if (first(j) == 0) first(j) = i
               last(j) = i
            end if
         end do
      end do
   end subroutine range_nonzero


   real(dp) function span_to_bandwidth(span, n, units) result(bandwidth)
      real(dp), intent(in) :: span
      integer, intent(in) :: n
      real(dp), intent(in), optional :: units
      real(dp) :: scale, gamma, centre, low, high, mid
      integer :: i, iteration, covered

      scale = real(n, dp)
      if (present(units)) scale = units
      if (n <= 0 .or. span <= 0.0_dp) then
         bandwidth = 0.0_dp
         return
      end if
      centre = (real(n/2,dp) + 0.5_dp)/scale
      low = 1.0e-8_dp
      high = max(1.0_dp, 2.0_dp*real(n,dp)/scale)
      do iteration = 1, 80
         mid = 0.5_dp*(low + high)
         covered = 0
         do i = 1, n
            gamma = real(i,dp)/scale
            if (epanechnikov((centre - gamma)/mid) > 0.0_dp) covered = covered + 1
         end do
         if (0.5_dp*real(covered,dp) < span) then
            low = mid
         else
            high = mid
         end if
      end do
      bandwidth = high
   end function span_to_bandwidth

end module ycevo_kernel
