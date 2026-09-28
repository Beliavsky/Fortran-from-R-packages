! SPDX-License-Identifier: GPL-3.0-only
! Derived from R package ftsa 6.7 by Rob Hyndman and Han Lin Shang.
module ftsa_simulation
   use r_kinds, only : dp
   implicit none
   private

   public :: ftsa_sim_farma, ftsa_stop_time_sim_data

contains

   subroutine ftsa_sim_farma(n, d, no_grid, seed_number, farma, simulated, status)
      integer, intent(in) :: n !! Number of retained functional observations; must be positive.
      real(dp), intent(in) :: d !! Fractional integration parameter used in the innovation convolution; must be positive.
      integer, intent(in) :: no_grid !! Number of equally spaced functional grid points; must exceed one.
      integer, intent(in) :: seed_number !! User seed offset corresponding to upstream `123 + seed_number` initialization.
      logical, intent(in) :: farma !! Apply the moving-average operator before the autoregressive operator when true.
      real(dp), allocatable, intent(out) :: simulated(:, :) !! Simulated grid-by-time FAR or FARMA observations.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid dimensions or fractional parameter.
      real(dp), allocatable :: beta(:), brown(:, :), eta(:, :), eta_ma(:, :), series(:, :)
      integer :: i, k, t

      if (n < 1 .or. no_grid < 2 .or. d <= 0.0_dp) then
         allocate(simulated(0, 0))
         status = 1
         return
      end if
      call set_rng_seed(123 + seed_number)
      allocate(brown(no_grid, 2 * n + 100), beta(n), eta(no_grid, n + 100))
      do t = 1, size(brown, 2)
         call brownian_path(brown(:, t))
      end do
      do i = 1, n
         beta(i) = exp(log_gamma(real(i, dp) + d) - log_gamma(real(i + 1, dp)) - log_gamma(d))
      end do
      do t = n + 1, 2 * n + 100
         eta(:, t - n) = brown(:, t)
         do k = 1, n
            eta(:, t - n) = eta(:, t - n) + beta(k) * brown(:, t - k)
         end do
      end do
      if (farma) then
         allocate(eta_ma(no_grid, size(eta, 2) - 1))
         do t = 2, size(eta, 2)
            eta_ma(:, t - 1) = apply_ma_operator(eta(:, t - 1)) + eta(:, t)
         end do
         allocate(series(no_grid, size(eta_ma, 2)))
         call apply_ar_recursion(eta_ma, series)
      else
         allocate(series(no_grid, size(eta, 2)))
         call apply_ar_recursion(eta, series)
      end if
      allocate(simulated(no_grid, n))
      simulated = series(:, size(series, 2) - n + 1:size(series, 2))
      status = 0
   end subroutine ftsa_sim_farma

   subroutine ftsa_stop_time_sim_data(sample_size, omega, seed_number, grid, data, status)
      integer, intent(in) :: sample_size !! Number of raw simulated curves before the final stationarizing transformation.
      real(dp), intent(in) :: omega !! Brownian-noise scale; nonnegative values are expected.
      integer, intent(in) :: seed_number !! User seed offset corresponding to upstream `123 + seed_number` initialization.
      real(dp), allocatable, intent(out) :: grid(:) !! Fixed 101-point grid from zero through one.
      real(dp), allocatable, intent(out) :: data(:, :) !! Absolute relative first differences, with `sample_size-1` columns.
      integer, intent(out) :: status !! Zero on success; nonzero when the requested sample size or noise scale is invalid.
      real(dp), allocatable :: path(:), raw(:, :), mean_curve(:)
      integer :: half, i, t

      if (sample_size < 3 .or. omega < 0.0_dp) then
         allocate(grid(0), data(0, 0))
         status = 1
         return
      end if
      call set_rng_seed(123 + seed_number)
      allocate(grid(101), mean_curve(101), path(101), raw(101, sample_size), data(101, sample_size - 1))
      do i = 1, 101
         grid(i) = real(i - 1, dp) / 100.0_dp
         mean_curve(i) = 10.0_dp * grid(i) * (1.0_dp - grid(i))
      end do
      call brownian_path(path)
      raw(:, 1) = mean_curve + omega * path
      half = ceiling(real(sample_size, dp) / 2.0_dp)
      do t = 2, half
         call brownian_path(path)
         raw(:, t) = 0.2_dp * raw(:, t - 1) + omega * path
      end do
      do t = half + 1, sample_size
         call brownian_path(path)
         raw(:, t) = 0.9_dp * raw(:, t - 1) + omega * path
      end do
      do t = 1, sample_size - 1
         data(:, t) = abs(raw(:, t + 1) - raw(:, t)) / (abs(raw(:, t)) + 0.1_dp)
      end do
      status = 0
   end subroutine ftsa_stop_time_sim_data

   subroutine apply_ar_recursion(innovations, series)
      real(dp), intent(in) :: innovations(:, :) !! Grid-by-time functional innovations.
      real(dp), intent(out) :: series(:, :) !! FAR(1) recursion using the upstream exponential integral kernel.
      integer :: t

      series(:, 1) = innovations(:, 1)
      do t = 2, size(innovations, 2)
         series(:, t) = apply_ar_operator(series(:, t - 1)) + innovations(:, t)
      end do
   end subroutine apply_ar_recursion

   pure function apply_ar_operator(x) result(y)
      real(dp), intent(in) :: x(:) !! Functional vector acted on by the FAR kernel.
      real(dp) :: y(size(x))
      integer :: i, j, n
      real(dp) :: xi, xj

      n = size(x)
      y = 0.0_dp
      do i = 1, n
         xi = real(i, dp) / real(n, dp)
         do j = 1, n
            xj = real(j, dp) / real(n, dp)
            y(i) = y(i) + 0.34_dp * exp(0.5_dp * (xi * xi + xj * xj)) * x(j)
         end do
      end do
      y = y / real(n, dp)
   end function apply_ar_operator

   pure function apply_ma_operator(x) result(y)
      real(dp), intent(in) :: x(:) !! Functional vector acted on by the upstream Brownian-covariance MA kernel.
      real(dp) :: y(size(x))
      integer :: i, j, n

      n = size(x)
      y = 0.0_dp
      do i = 1, n
         do j = 1, n
            y(i) = y(i) + 1.5_dp * real(min(i, j), dp) / real(n, dp) * x(j)
         end do
      end do
      y = y / real(n, dp)
   end function apply_ma_operator

   subroutine brownian_path(path)
      real(dp), intent(out) :: path(:) !! Brownian-motion path on `[0,1]`, including the zero starting point.
      real(dp) :: dt, z
      integer :: i

      path(1) = 0.0_dp
      if (size(path) < 2) return
      dt = 1.0_dp / real(size(path) - 1, dp)
      do i = 2, size(path)
         call standard_normal(z)
         path(i) = path(i - 1) + sqrt(dt) * z
      end do
   end subroutine brownian_path

   subroutine standard_normal(z)
      real(dp), intent(out) :: z !! One standard-normal variate generated by the Box-Muller transform.
      real(dp) :: u1, u2
      real(dp), parameter :: two_pi = 6.2831853071795864769252867665590058_dp

      call random_number(u1)
      call random_number(u2)
      u1 = max(u1, tiny(1.0_dp))
      z = sqrt(-2.0_dp * log(u1)) * cos(two_pi * u2)
   end subroutine standard_normal

   subroutine set_rng_seed(seed)
      integer, intent(in) :: seed !! Scalar seed expanded deterministically to the processor RNG seed vector.
      integer, allocatable :: values(:)
      integer :: i, n

      call random_seed(size=n)
      allocate(values(n))
      do i = 1, n
         values(i) = modulo(seed + 104729 * i, huge(1) - 1)
         if (values(i) <= 0) values(i) = i
      end do
      call random_seed(put=values)
   end subroutine set_rng_seed

end module ftsa_simulation
