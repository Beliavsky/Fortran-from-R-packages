! SPDX-License-Identifier: GPL-2.0-or-later
module binom_profile_module
   use binom_types
   use binom_numerics
   use r_distributions, only: r_qnorm
   use r_transforms, only: r_expm1
   use r_linalg, only: solve_system
   use r_ordering, only: r_order
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private
   public :: binom_profile
contains
   pure elemental function binom_profile(x, n, conf_level, maxsteps, del, bayes) result(ci)
      real(dp), intent(in) :: x !! Success count between zero and n.
      real(dp), intent(in) :: n !! Positive trial count.
      real(dp), intent(in), optional :: conf_level !! Confidence probability strictly between zero and one; default 0.95.
      integer, intent(in), optional :: maxsteps !! Maximum profile steps per direction, 1 through 4096; default 50.
      real(dp), intent(in), optional :: del !! Positive logit-profile step multiplier; default zmax/5, as in R.
      logical, intent(in), optional :: bayes !! Apply Jeffreys boundary adjustment; default true.
      type(binom_interval) :: ci
      real(dp) :: level, xx, nn, adjustment, alpha, zmax, step_size, eta, error, mu, z, q
      real(dp), allocatable :: probability(:), deviate(:), grid_p(:), grid_z(:)
      logical :: adjust, left_edge, right_edge
      integer :: steps, count, sign_step, step, info
      level = 0.95_dp
      if (present(conf_level)) level = conf_level
      steps = 50
      if (present(maxsteps)) steps = maxsteps
      adjust = .true.
      if (present(bayes)) adjust = bayes
      ci%mean = binom_nan()
      ci%lower = binom_nan()
      ci%upper = binom_nan()
      ci%status = binom_domain
      if (.not. valid_binomial(x, n, level) .or. steps < 1) return
      if (steps > 4096) then
         ci%status = binom_resource
         return
      end if
      ci%mean = x / n
      alpha = 1.0_dp - level
      if ((x == 0.0_dp .or. x == n) .and. .not. adjust) then
         ci%lower = 0.0_dp
         ci%upper = 1.0_dp
         if (x == 0.0_dp) ci%upper = -r_expm1(log(0.5_dp * alpha) / n)
         if (x == n) ci%lower = exp(log(0.5_dp * alpha) / n)
         ci%status = binom_ok
         return
      end if
      xx = x
      nn = n
      left_edge = x == 0.0_dp .and. adjust
      right_edge = x == n .and. adjust
      adjustment = 1.0_dp
      if (left_edge .or. right_edge) then
         xx = x + 0.5_dp
         nn = n + 1.0_dp
         adjustment = 0.5_dp
      end if
      if (nn > 1.0_dp) adjustment = min(2.0_dp / log(nn), adjustment)
      zmax = r_qnorm(0.5_dp * alpha / nn, lower_tail=.false.)
      step_size = zmax / 5.0_dp
      if (present(del)) step_size = del
      if (.not. ieee_is_finite(step_size) .or. step_size <= 0.0_dp) return
      mu = xx / nn
      eta = logit(mu)
      error = sqrt(1.0_dp / (nn * mu * (1.0_dp - mu))) * adjustment
      allocate(probability(min(2 * steps, 512)), deviate(min(2 * steps, 512)))
      count = 0
      do sign_step = -1, 1, 2
         if (left_edge .and. sign_step == -1) cycle
         if (right_edge .and. sign_step == 1) cycle
         z = 0.0_dp
         do step = 1, steps
            if (abs(z) >= zmax) exit
            if (count == size(probability)) then
               ci%status = binom_resource
               return
            end if
            ! Match stats::binomial()$linkinv endpoint clipping at |eta| > 30.
            q = eta + real(sign_step * step, dp) * step_size * error
            if (q < -30.0_dp) then
               mu = epsilon(1.0_dp)
            else if (q > 30.0_dp) then
               mu = 1.0_dp - epsilon(1.0_dp)
            else
               mu = logistic(q)
            end if
            z = real(sign_step, dp) * sqrt(binomial_deviance(mu, xx, nn))
            count = count + 1
            probability(count) = mu
            deviate(count) = z
         end do
      end do
      call sort_pairs(probability(:count), deviate(:count))
      call collapse_ties(probability, deviate, count)
      if (count < 2) then
         ci%status = binom_convergence
         return
      end if
      allocate(grid_p(3 * count), grid_z(3 * count))
      do step = 1, 3 * count
         grid_p(step) = probability(1) + (probability(count) - probability(1)) * &
            real(step - 1, dp) / real(3 * count - 1, dp)
      end do
      call profile_fmm(probability(:count), deviate(:count), grid_p, grid_z, info)
      if (info /= 0) then
         ci%status = binom_convergence
         return
      end if
      call sort_pairs(grid_z, grid_p)
      q = r_qnorm(0.5_dp * alpha, lower_tail=.false.)
      ci%lower = 0.0_dp
      ci%upper = 1.0_dp
      if (.not. left_edge) ci%lower = interpolate_profile(grid_z, grid_p, -q)
      if (.not. right_edge) ci%upper = interpolate_profile(grid_z, grid_p, q)
      ci%status = binom_ok
      if (.not. ieee_is_finite(ci%lower) .or. .not. ieee_is_finite(ci%upper)) ci%status = binom_convergence
   end function binom_profile

   pure subroutine profile_fmm(x, y, grid, value, status)
      real(dp), intent(in) :: x(:) !! Strictly increasing sampled probabilities, at least two entries.
      real(dp), intent(in) :: y(:) !! Signed square-root deviance at each probability.
      real(dp), intent(in) :: grid(:) !! Ordered interpolation points inside the sampled range.
      real(dp), intent(out) :: value(:) !! FMM cubic interpolated deviances, same size as grid.
      integer, intent(out) :: status !! Zero on success, otherwise the linear-system status.
      real(dp), allocatable :: matrix(:, :), rhs(:), second(:), h(:), slope(:)
      real(dp) :: cubic_left, cubic_right, dx, fraction
      integer :: m, i, j
      m = size(x)
      h = x(2:) - x(:m - 1)
      slope = (y(2:) - y(:m - 1)) / h
      status = 0
      if (m == 2) then
         value = y(1) + (grid - x(1)) * slope(1)
         return
      end if
      if (m == 3) then
         fraction = (slope(2) - slope(1)) / (x(3) - x(1))
         value = y(1) + (grid - x(1)) * (slope(1) + fraction * (grid - x(2)))
         return
      end if
      allocate(matrix(m, m), rhs(m), second(m))
      matrix = 0.0_dp
      rhs = 0.0_dp
      ! FMM endpoint third derivatives are those of the first/last four-point cubic.
      cubic_left = ((slope(3) - slope(2)) / (x(4) - x(2)) - &
         (slope(2) - slope(1)) / (x(3) - x(1))) / (x(4) - x(1))
      cubic_right = ((slope(m - 1) - slope(m - 2)) / (x(m) - x(m - 2)) - &
         (slope(m - 2) - slope(m - 3)) / (x(m - 1) - x(m - 3))) / (x(m) - x(m - 3))
      matrix(1, 1) = -1.0_dp
      matrix(1, 2) = 1.0_dp
      rhs(1) = 6.0_dp * h(1) * cubic_left
      matrix(m, m - 1) = -1.0_dp
      matrix(m, m) = 1.0_dp
      rhs(m) = 6.0_dp * h(m - 1) * cubic_right
      do i = 2, m - 1
         matrix(i, i - 1) = h(i - 1)
         matrix(i, i) = 2.0_dp * (h(i - 1) + h(i))
         matrix(i, i + 1) = h(i)
         rhs(i) = 6.0_dp * (slope(i) - slope(i - 1))
      end do
      call solve_system(matrix, rhs, second, status)
      if (status /= 0) return
      j = 1
      do i = 1, size(grid)
         do while (j < m - 1)
            if (grid(i) <= x(j + 1)) exit
            j = j + 1
         end do
         dx = grid(i) - x(j)
         value(i) = y(j) + dx * (slope(j) - h(j) * (2.0_dp * second(j) + second(j + 1)) / 6.0_dp + &
            dx * (0.5_dp * second(j) + dx * (second(j + 1) - second(j)) / (6.0_dp * h(j))))
      end do
   end subroutine profile_fmm

   pure subroutine sort_pairs(x, y)
      real(dp), intent(inout) :: x(:) !! Keys sorted into ascending order in place.
      real(dp), intent(inout) :: y(:) !! Same-length paired values reordered with the keys.
      integer, allocatable :: order(:)
      call r_order(x, order)
      x = x(order)
      y = y(order)
   end subroutine sort_pairs

   pure subroutine collapse_ties(x, y, count)
      real(dp), intent(inout) :: x(:) !! Sorted probability keys, collapsed in the first count entries.
      real(dp), intent(inout) :: y(:) !! Deviances averaged when repeated keys occur.
      integer, intent(inout) :: count !! Number of used entries, replaced with the distinct-key count.
      integer :: i, j, used
      real(dp) :: total
      i = 1
      used = 0
      do while (i <= count)
         j = i + 1
         total = y(i)
         do while (j <= count)
            if (x(j) /= x(i)) exit
            total = total + y(j)
            j = j + 1
         end do
         used = used + 1
         x(used) = x(i)
         y(used) = total / real(j - i, dp)
         i = j
      end do
      count = used
   end subroutine collapse_ties

   pure function interpolate_profile(x, y, point) result(value)
      real(dp), intent(in) :: x(:) !! Sorted spline deviances.
      real(dp), intent(in) :: y(:) !! Paired probabilities.
      real(dp), intent(in) :: point !! Target standard-normal deviate.
      real(dp) :: value
      integer :: j
      value = binom_nan()
      if (point < x(1) .or. point > x(size(x))) return
      do j = 1, size(x) - 1
         if (x(j + 1) < point .or. x(j + 1) == x(j)) cycle
         value = y(j) + (point - x(j)) * (y(j + 1) - y(j)) / (x(j + 1) - x(j))
         return
      end do
   end function interpolate_profile
end module binom_profile_module
