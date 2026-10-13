! Numerical algorithms placed in the public domain by Burns Statistics Ltd.
! Modern Fortran translation: see LICENSE and NOTICE.md.
module burstfin_portfolio
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
   use r_kinds, only: dp
   use r_ordering, only: r_sort_values_in_place
   implicit none
   private
   public :: dp, alpha_proxy, alpha_proxy_grid, alpha_surface, slide_weight
   public :: var_add_benchmark, var_relative_benchmark

   type :: alpha_surface
      real(dp), allocatable :: x(:), y(:), z(:, :)
      integer :: x_argument = 0, y_argument = 0
   end type alpha_surface

   interface var_add_benchmark
      module procedure add_benchmark_matrix, add_benchmark_array
   end interface
   interface var_relative_benchmark
      module procedure relative_benchmark_matrix, relative_benchmark_array
   end interface
contains
   pure elemental real(dp) function alpha_proxy(weight, vol_man, vol_bench, vol_other, cor_man, cor_bench) result(alpha)
      real(dp), intent(in), optional :: weight !! Manager portfolio weight strictly between zero and one; default 0.2.
      real(dp), intent(in), optional :: vol_man !! Positive manager volatility in decimal units; default 0.2.
      real(dp), intent(in), optional :: vol_bench !! Positive benchmark volatility in decimal units; default 0.2.
      real(dp), intent(in), optional :: vol_other !! Positive volatility of the remaining portfolio; default 0.2.
      real(dp), intent(in), optional :: cor_man !! Manager correlation with the rest, in [-1,1]; default 0.2.
      real(dp), intent(in), optional :: cor_bench !! Benchmark correlation with the rest, in [-1,1]; default 0.2.
      real(dp) :: w, vm, vb, vo, cm, cb
      w = 0.2_dp
      vm = 0.2_dp
      vb = 0.2_dp
      vo = 0.2_dp
      cm = 0.2_dp
      cb = 0.2_dp
      if (present(weight)) w = weight
      if (present(vol_man)) vm = vol_man
      if (present(vol_bench)) vb = vol_bench
      if (present(vol_other)) vo = vol_other
      if (present(cor_man)) cm = cor_man
      if (present(cor_bench)) cb = cor_bench
      alpha = ieee_value(0.0_dp, ieee_quiet_nan)
      if (.not. all(ieee_is_finite([w, vm, vb, vo, cm, cb]))) return
      if (w <= 0.0_dp .or. w >= 1.0_dp .or. min(vm, vb, vo) <= 0.0_dp) return
      if (max(abs(cm), abs(cb)) > 1.0_dp) return
      alpha = 10000.0_dp * (-w * (vm ** 2 - vb ** 2) - &
         2.0_dp * (1.0_dp - w) * vo * (vm * cm - vb * cb))
   end function alpha_proxy

   pure subroutine alpha_proxy_grid(weight, vol_man, vol_bench, vol_other, cor_man, cor_bench, surface, status, transpose_output)
      real(dp), intent(in) :: weight(:) !! One weight or the weight grid, with entries in (0,1).
      real(dp), intent(in) :: vol_man(:) !! One positive manager volatility or a volatility grid.
      real(dp), intent(in) :: vol_bench(:) !! One positive benchmark volatility or a volatility grid.
      real(dp), intent(in) :: vol_other(:) !! One positive remaining-portfolio volatility or its grid.
      real(dp), intent(in) :: cor_man(:) !! One manager correlation or a grid with entries in [-1,1].
      real(dp), intent(in) :: cor_bench(:) !! One benchmark correlation or a grid with entries in [-1,1].
      type(alpha_surface), intent(out) :: surface !! Axis values and proxy matrix, in the upstream argument order.
      integer, intent(out) :: status !! Zero on success; one for invalid values, empty input, or more than two grids.
      logical, intent(in), optional :: transpose_output !! Exchange the two variable axes, as in R; default false.
      real(dp), allocatable :: args(:, :), tmp(:), ztmp(:, :)
      integer :: lengths(6), axes(2), i, j, a, nx, ny
      real(dp) :: v(6)
      status = 1
      lengths = [size(weight), size(vol_man), size(vol_bench), size(vol_other), size(cor_man), size(cor_bench)]
      if (minval(lengths) == 0 .or. count(lengths > 1) > 2) return
      allocate(args(maxval(lengths), 6))
      args = 0.0_dp
      args(:lengths(1), 1) = weight
      args(:lengths(2), 2) = vol_man
      args(:lengths(3), 3) = vol_bench
      args(:lengths(4), 4) = vol_other
      args(:lengths(5), 5) = cor_man
      args(:lengths(6), 6) = cor_bench
      if (.not. all(ieee_is_finite(args))) return
      if (any(weight <= 0.0_dp) .or. any(weight >= 1.0_dp)) return
      if (minval(vol_man) <= 0.0_dp .or. minval(vol_bench) <= 0.0_dp .or. minval(vol_other) <= 0.0_dp) return
      if (any(abs(cor_man) > 1.0_dp) .or. any(abs(cor_bench) > 1.0_dp)) return
      axes = 0
      j = 0
      do a = 1, 6
         if (lengths(a) == 1) cycle
         j = j + 1
         axes(j) = a
         if (count(lengths > 1) == 2) call r_sort_values_in_place(args(:lengths(a), a))
      end do
      nx = 1
      ny = 1
      if (axes(1) > 0) nx = lengths(axes(1))
      if (axes(2) > 0) ny = lengths(axes(2))
      allocate(surface%x(nx), surface%y(ny), surface%z(nx, ny))
      surface%x = 0.0_dp
      surface%y = 0.0_dp
      if (axes(1) > 0) surface%x = args(:nx, axes(1))
      if (axes(2) > 0) surface%y = args(:ny, axes(2))
      do j = 1, ny
         do i = 1, nx
            v = args(1, :)
            if (axes(1) > 0) v(axes(1)) = args(i, axes(1))
            if (axes(2) > 0) v(axes(2)) = args(j, axes(2))
            surface%z(i, j) = alpha_proxy(v(1), v(2), v(3), v(4), v(5), v(6))
         end do
      end do
      if (present(transpose_output)) then
         if (transpose_output .and. axes(2) > 0) then
            tmp = surface%x
            surface%x = surface%y
            surface%y = tmp
            ztmp = surface%z
            surface%z = transpose(ztmp)
            axes = axes(2:1:-1)
         end if
      end if
      surface%x_argument = axes(1)
      surface%y_argument = axes(2)
      status = 0
   end subroutine alpha_proxy_grid

   pure subroutine slide_weight(n, weights, status, fractions, observations, locations)
      integer, intent(in) :: n !! Positive number of chronologically ordered observations.
      real(dp), allocatable, intent(out) :: weights(:) !! Nondecreasing weights of length n between zero and one.
      integer, intent(out) :: status !! Zero on success; one for invalid n, coordinates, or an all-zero specification.
      real(dp), intent(in), optional :: fractions(2) !! Ramp endpoints as fractions of n; default [0,1].
      real(dp), intent(in), optional :: observations(2) !! Endpoint ages, converted to locations by n minus age.
      real(dp), intent(in), optional :: locations(2) !! Absolute endpoints; mutually exclusive with observations.
      real(dp) :: loc(2), lo, hi, v
      integer :: i
      status = 1
      if (n < 1) return
      if (present(locations) .and. present(observations)) return
      loc = [0.0_dp, real(n, dp)]
      if (present(fractions)) loc = fractions * real(n, dp)
      if (present(observations)) loc = real(n, dp) - observations
      if (present(locations)) loc = locations
      if (.not. all(ieee_is_finite(loc))) return
      do i = 1, 2
         ! R rounds exact halves to the nearest even integer.
         v = aint(loc(i))
         if (v > loc(i)) v = v - 1.0_dp
         if (loc(i) - v > 0.5_dp) v = v + 1.0_dp
         if (loc(i) - v == 0.5_dp .and. modulo(v, 2.0_dp) /= 0.0_dp) v = v + 1.0_dp
         loc(i) = v
      end do
      lo = minval(loc)
      hi = maxval(loc)
      if (lo >= real(n, dp)) return
      allocate(weights(n))
      do i = 1, n
         weights(i) = 1.0_dp
         if (real(i, dp) >= lo .and. real(i, dp) <= hi) then
            weights(i) = (real(i, dp) - lo) / (hi - lo + 1.0_dp)
         end if
         if (lo > 1.0_dp .and. real(i, dp) <= lo) weights(i) = 0.0_dp
      end do
      status = 0
   end subroutine slide_weight

   pure subroutine add_benchmark_matrix(variance, weights, result, status, sum_to_one, asset_indices)
      real(dp), intent(in) :: variance(:, :) !! Finite square covariance matrix in caller-defined asset order.
      real(dp), intent(in) :: weights(:) !! Finite benchmark holdings, optionally mapped by asset_indices.
      real(dp), allocatable, intent(out) :: result(:, :) !! Covariance with a benchmark appended as the final asset.
      integer, intent(out) :: status !! Zero on success; one for invalid shape, holdings, mapping, or values.
      logical, intent(in), optional :: sum_to_one !! Normalize absolute holdings to one when needed; default true.
      integer, intent(in), optional :: asset_indices(:) !! Unique one-based asset positions corresponding to holdings.
      real(dp), allocatable :: w(:), bcov(:)
      real(dp) :: total, largest
      logical :: normalize
      integer :: p, j, k
      status = 1
      p = size(variance, 1)
      if (p < 1 .or. size(variance, 2) /= p) return
      if (.not. all(ieee_is_finite(variance)) .or. .not. all(ieee_is_finite(weights))) return
      if (size(weights) == 0) return
      allocate(w(p))
      w = 0.0_dp
      if (present(asset_indices)) then
         if (size(asset_indices) /= size(weights)) return
         if (any(asset_indices < 1) .or. any(asset_indices > p)) return
         do j = 1, size(weights)
            do k = 1, j - 1
               if (asset_indices(j) == asset_indices(k)) return
            end do
            w(asset_indices(j)) = weights(j)
         end do
      else
         if (size(weights) /= p) return
         w = weights
      end if
      normalize = .true.
      if (present(sum_to_one)) normalize = sum_to_one
      if (normalize) then
         largest = maxval(abs(w))
         if (largest == 0.0_dp) return
         total = sum(abs(w / largest))
         if (largest > huge(1.0_dp) / total) then
            w = (w / largest) / total
         else if (abs(largest * total - 1.0_dp) > 1.0e-10_dp) then
            w = (w / largest) / total
         end if
      end if
      bcov = matmul(variance, w)
      allocate(result(p + 1, p + 1))
      result(:p, :p) = variance
      result(:p, p + 1) = bcov
      result(p + 1, :p) = bcov
      result(p + 1, p + 1) = dot_product(w, bcov)
      if (.not. all(ieee_is_finite(result))) return
      status = 0
   end subroutine add_benchmark_matrix

   pure subroutine add_benchmark_array(variance, weights, result, status, sum_to_one, asset_indices)
      real(dp), intent(in) :: variance(:, :, :) !! Covariance cube with assets on its first two axes.
      real(dp), intent(in) :: weights(:) !! Benchmark holdings applied identically to all slices.
      real(dp), allocatable, intent(out) :: result(:, :, :) !! Cube with one benchmark row and column appended per slice.
      integer, intent(out) :: status !! Zero on success; one for invalid dimensions, mappings, or values.
      logical, intent(in), optional :: sum_to_one !! Normalize absolute holdings to one; default true.
      integer, intent(in), optional :: asset_indices(:) !! One-based asset positions for potentially sparse holdings.
      real(dp), allocatable :: slice(:, :)
      integer :: j, p
      p = size(variance, 1)
      status = 1
      if (size(variance, 3) < 1) return
      allocate(result(p + 1, p + 1, size(variance, 3)))
      do j = 1, size(variance, 3)
         call add_benchmark_matrix(variance(:, :, j), weights, slice, status, sum_to_one, asset_indices)
         if (status /= 0) return
         result(:, :, j) = slice
      end do
   end subroutine add_benchmark_array

   pure subroutine relative_benchmark_matrix(variance, benchmark, result, status)
      real(dp), intent(in) :: variance(:, :) !! Square covariance matrix including the benchmark asset.
      integer, intent(in) :: benchmark !! One-based index of the benchmark row and column to remove.
      real(dp), allocatable, intent(out) :: result(:, :) !! Covariance of each remaining asset minus the benchmark.
      integer, intent(out) :: status !! Zero on success; one for invalid dimensions or benchmark index.
      integer :: p, i, j, ii, jj
      p = size(variance, 1)
      status = 1
      if (p < 1 .or. size(variance, 2) /= p) return
      if (benchmark < 1 .or. benchmark > p) return
      allocate(result(p - 1, p - 1))
      jj = 0
      do j = 1, p
         if (j == benchmark) cycle
         jj = jj + 1
         ii = 0
         do i = 1, p
            if (i == benchmark) cycle
            ii = ii + 1
            ! Upstream uses the benchmark column on both sides, also for nonsymmetric input.
            result(ii, jj) = variance(j, i) + variance(benchmark, benchmark) - &
               variance(i, benchmark) - variance(j, benchmark)
         end do
      end do
      status = 0
   end subroutine relative_benchmark_matrix

   pure subroutine relative_benchmark_array(variance, benchmark, result, status)
      real(dp), intent(in) :: variance(:, :, :) !! Covariance cube containing the benchmark in each slice.
      integer, intent(in) :: benchmark !! One-based benchmark index shared by every slice.
      real(dp), allocatable, intent(out) :: result(:, :, :) !! Benchmark-relative covariance cube, with benchmark removed.
      integer, intent(out) :: status !! Zero on success; one for invalid shape or benchmark index.
      real(dp), allocatable :: slice(:, :)
      integer :: p, j
      p = size(variance, 1)
      status = 1
      if (p < 1 .or. size(variance, 3) < 1) return
      allocate(result(p - 1, p - 1, size(variance, 3)))
      do j = 1, size(variance, 3)
         call relative_benchmark_matrix(variance(:, :, j), benchmark, slice, status)
         if (status /= 0) return
         result(:, :, j) = slice
      end do
   end subroutine relative_benchmark_array
end module burstfin_portfolio
