! SPDX-License-Identifier: GPL-3.0-only
module kza_filters
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_quiet_nan, ieee_value
   use kza_kinds, only : dp
   use kza_utils, only : adaptive_factor, normalizer, quiet_nan, r_round_even
   implicit none
   private

   public :: kz
   public :: kza
   public :: kzsv

   interface kz
      module procedure kz_1d
      module procedure kz_2d
      module procedure kz_3d
   end interface kz

   interface kza
      module procedure kza_1d
      module procedure kza_2d
      module procedure kza_3d
   end interface kza

contains

   pure function kz_1d(x, m, k) result(filtered)
      real(dp), intent(in) :: x(:) !! One-dimensional raw data; non-finite values are omitted from each moving average.
      integer, intent(in) :: m !! Full R-level KZ window width; the internal radius is floor(m/2).
      integer, intent(in), optional :: k !! Number of moving-average iterations; defaults to 3.
      real(dp), allocatable :: filtered(:)
      real(dp), allocatable :: work(:)
      integer :: center
      integer :: iter
      integer :: niter
      integer :: q

      q = max(m / 2, 0)
      niter = 3
      if (present(k)) niter = max(k, 0)
      allocate(work(size(x)), filtered(size(x)))
      work = x
      filtered = x
      do iter = 1, niter
         do center = 1, size(x)
            filtered(center) = finite_mean_1d(work, max(1, center - q), min(size(x), center + q))
         end do
         work = filtered
      end do
   end function kz_1d

   pure function kz_2d(x, m, k, m2) result(filtered)
      real(dp), intent(in) :: x(:, :) !! Two-dimensional raw field; non-finite cells are omitted from each box average.
      integer, intent(in) :: m !! Full window width along the first array dimension; also used for dimension two when m2 is absent.
      integer, intent(in), optional :: k !! Number of moving-average iterations; defaults to 3.
      integer, intent(in), optional :: m2 !! Full window width along the second array dimension; defaults to m.
      real(dp), allocatable :: filtered(:, :)
      real(dp), allocatable :: work(:, :)
      integer :: i
      integer :: iter
      integer :: j
      integer :: niter
      integer :: q1
      integer :: q2

      q1 = max(m / 2, 0)
      q2 = q1
      if (present(m2)) q2 = max(m2 / 2, 0)
      niter = 3
      if (present(k)) niter = max(k, 0)
      allocate(work(size(x, 1), size(x, 2)), filtered(size(x, 1), size(x, 2)))
      work = x
      filtered = x
      do iter = 1, niter
         do j = 1, size(x, 2)
            do i = 1, size(x, 1)
               filtered(i, j) = finite_mean_2d(work, i, j, q1, q1, q2, q2)
            end do
         end do
         work = filtered
      end do
   end function kz_2d

   pure function kz_3d(x, m, k, m2, m3) result(filtered)
      real(dp), intent(in) :: x(:, :, :) !! Three-dimensional raw field; non-finite cells are omitted from each box average.
      integer, intent(in) :: m !! Full window width along the first array dimension; also used for other dimensions when omitted.
      integer, intent(in), optional :: k !! Number of moving-average iterations; defaults to 3.
      integer, intent(in), optional :: m2 !! Full window width along the second array dimension; defaults to m.
      integer, intent(in), optional :: m3 !! Full window width along the third array dimension; defaults to m.
      real(dp), allocatable :: filtered(:, :, :)
      real(dp), allocatable :: work(:, :, :)
      integer :: i
      integer :: iter
      integer :: j
      integer :: l
      integer :: niter
      integer :: q1
      integer :: q2
      integer :: q3

      q1 = max(m / 2, 0)
      q2 = q1
      q3 = q1
      if (present(m2)) q2 = max(m2 / 2, 0)
      if (present(m3)) q3 = max(m3 / 2, 0)
      niter = 3
      if (present(k)) niter = max(k, 0)
      allocate(work(size(x, 1), size(x, 2), size(x, 3)))
      allocate(filtered(size(x, 1), size(x, 2), size(x, 3)))
      work = x
      filtered = x
      do iter = 1, niter
         do l = 1, size(x, 3)
            do j = 1, size(x, 2)
               do i = 1, size(x, 1)
                  filtered(i, j, l) = finite_mean_3d(work, i, j, l, q1, q1, q2, q2, q3, q3)
               end do
            end do
         end do
         work = filtered
      end do
   end function kz_3d

   function kza_1d(x, m, y, k, min_size, tol, impute_tails, normalize) result(filtered)
      real(dp), intent(in) :: x(:) !! One-dimensional raw signal to filter adaptively.
      integer, intent(in) :: m !! Full R-level KZA window width; the internal per-side radius is floor(m/2).
      real(dp), intent(in), optional :: y(:) !! Optional KZ baseline; when absent, kz(x,m,k) is computed internally.
      integer, intent(in), optional :: k !! Number of adaptive-filter iterations and default-baseline KZ iterations; defaults to 3.
      integer, intent(in), optional :: min_size !! Minimum per-side adaptive radius; defaults to R round(0.05*m).
      real(dp), intent(in), optional :: tol !! Absolute derivative threshold treated as zero; defaults to 1.0e-5.
      logical, intent(in), optional :: impute_tails !! When false, reproduce R's one-dimensional tail NA marking; defaults to false.
      character(len=*), intent(in), optional :: normalize !! Shrink normalizer: "max" or "quantile" (99th percentile); defaults to "max".
      real(dp), allocatable :: filtered(:)
      real(dp), allocatable :: baseline(:)
      integer :: i
      integer :: min_window
      integer :: niter
      integer :: q
      real(dp) :: eps
      real(dp) :: nprob
      logical :: keep_tails

      niter = 3
      if (present(k)) niter = k
      min_window = r_round_even(0.05_dp * real(m, dp))
      if (present(min_size)) min_window = min_size
      eps = 1.0e-5_dp
      if (present(tol)) eps = tol
      nprob = normalization_probability(normalize)
      q = max(m / 2, 0)

      if (present(y)) then
         baseline = y
      else
         baseline = kz_1d(x, m, niter)
      end if
      filtered = kza_core_1d(x, baseline, q, niter, min_window, eps, nprob)

      keep_tails = .false.
      if (present(impute_tails)) keep_tails = impute_tails
      if (.not. keep_tails .and. size(filtered) > 0) then
         do i = 1, min(max(m, 0), size(filtered))
            filtered(i) = quiet_nan()
         end do
         do i = max(1, size(filtered) - max(m, 0)), size(filtered)
            filtered(i) = quiet_nan()
         end do
      end if
   end function kza_1d

   function kza_2d(x, m, y, k, min_size, tol, impute_tails, symmetrize, normalize, m2) result(filtered)
      real(dp), intent(in) :: x(:, :) !! Two-dimensional raw field to filter adaptively.
      integer, intent(in) :: m !! Full window width along the first array dimension; also used for dimension two when m2 is absent.
      real(dp), intent(in), optional :: y(:, :) !! Optional KZ baseline; when absent, a matching KZ baseline is computed internally.
      integer, intent(in), optional :: k !! Number of adaptive and default-baseline iterations; defaults to 3.
      integer, intent(in), optional :: min_size !! Minimum per-side adaptive radius in every direction; defaults to R round(0.05*m).
      real(dp), intent(in), optional :: tol !! Absolute derivative threshold treated as zero; defaults to 1.0e-5.
      logical, intent(in), optional :: impute_tails !! Accepted for R-interface parity; tail NA marking is intentionally one-dimensional only.
      logical, intent(in), optional :: symmetrize !! Average the adaptive filter over four 90-degree rotations; defaults to false.
      character(len=*), intent(in), optional :: normalize !! Shrink normalizer: "max" or "quantile" (99th percentile); defaults to "max".
      integer, intent(in), optional :: m2 !! Full window width along the second array dimension; defaults to m.
      real(dp), allocatable :: filtered(:, :)
      real(dp), allocatable :: acc(:, :)
      real(dp), allocatable :: back(:, :)
      real(dp), allocatable :: baseline(:, :)
      real(dp), allocatable :: xr(:, :)
      real(dp), allocatable :: yr(:, :)
      real(dp), allocatable :: zr(:, :)
      integer :: min_window
      integer :: niter
      integer :: q1
      integer :: q2
      integer :: rotation
      real(dp) :: eps
      real(dp) :: nprob
      logical :: do_symmetrize

      niter = 3
      if (present(k)) niter = k
      min_window = r_round_even(0.05_dp * real(m, dp))
      if (present(min_size)) min_window = min_size
      eps = 1.0e-5_dp
      if (present(tol)) eps = tol
      nprob = normalization_probability(normalize)
      q1 = max(m / 2, 0)
      q2 = q1
      if (present(m2)) q2 = max(m2 / 2, 0)

      if (present(y)) then
         baseline = y
      else if (present(m2)) then
         baseline = kz_2d(x, m, niter, m2)
      else
         baseline = kz_2d(x, m, niter)
      end if

      do_symmetrize = .false.
      if (present(symmetrize)) do_symmetrize = symmetrize
      if (.not. do_symmetrize) then
         filtered = kza_core_2d(x, baseline, q1, q2, niter, min_window, eps, nprob)
         return
      end if

      allocate(acc(size(x, 1), size(x, 2)))
      acc = 0.0_dp
      do rotation = 0, 3
         xr = rotate_n_2d(x, rotation)
         yr = rotate_n_2d(baseline, rotation)
         zr = kza_core_2d(xr, yr, q1, q2, niter, min_window, eps, nprob)
         back = rotate_n_2d(zr, modulo(4 - rotation, 4))
         acc = acc + back
      end do
      filtered = acc / 4.0_dp
   end function kza_2d

   function kza_3d(x, m, y, k, min_size, tol, impute_tails, normalize, m2, m3) result(filtered)
      real(dp), intent(in) :: x(:, :, :) !! Three-dimensional raw field to filter adaptively.
      integer, intent(in) :: m !! Full window width along the first array dimension; also used for other dimensions when omitted.
      real(dp), intent(in), optional :: y(:, :, :) !! Optional KZ baseline; when absent, a matching KZ baseline is computed internally.
      integer, intent(in), optional :: k !! Number of adaptive and default-baseline iterations; defaults to 3.
      integer, intent(in), optional :: min_size !! Minimum per-side adaptive radius in every direction; defaults to R round(0.05*m).
      real(dp), intent(in), optional :: tol !! Absolute derivative threshold treated as zero; defaults to 1.0e-5.
      logical, intent(in), optional :: impute_tails !! Accepted for R-interface parity; tail NA marking is intentionally one-dimensional only.
      character(len=*), intent(in), optional :: normalize !! Shrink normalizer: "max" or "quantile" (99th percentile); defaults to "max".
      integer, intent(in), optional :: m2 !! Full window width along the second array dimension; defaults to m.
      integer, intent(in), optional :: m3 !! Full window width along the third array dimension; defaults to m.
      real(dp), allocatable :: filtered(:, :, :)
      real(dp), allocatable :: baseline(:, :, :)
      integer :: min_window
      integer :: niter
      integer :: q1
      integer :: q2
      integer :: q3
      real(dp) :: eps
      real(dp) :: nprob

      niter = 3
      if (present(k)) niter = k
      min_window = r_round_even(0.05_dp * real(m, dp))
      if (present(min_size)) min_window = min_size
      eps = 1.0e-5_dp
      if (present(tol)) eps = tol
      nprob = normalization_probability(normalize)
      q1 = max(m / 2, 0)
      q2 = q1
      q3 = q1
      if (present(m2)) q2 = max(m2 / 2, 0)
      if (present(m3)) q3 = max(m3 / 2, 0)

      if (present(y)) then
         baseline = y
      else if (present(m2) .and. present(m3)) then
         baseline = kz_3d(x, m, niter, m2, m3)
      else if (present(m2)) then
         baseline = kz_3d(x, m, niter, m2)
      else if (present(m3)) then
         baseline = kz_3d(x, m, niter, m3=m3)
      else
         baseline = kz_3d(x, m, niter)
      end if
      filtered = kza_core_3d(x, baseline, q1, q2, q3, niter, min_window, eps, nprob)
   end function kza_3d

   function kzsv(kza_data, kz_data, window, min_size, tol) result(sigma)
      real(dp), intent(in) :: kza_data(:) !! Adaptive-filtered one-dimensional series whose local sample variance is estimated.
      real(dp), intent(in) :: kz_data(:) !! KZ baseline used to determine the same adaptive head/tail windows as upstream kzsv.
      integer, intent(in) :: window !! Full R-level KZA window width; the internal per-side radius is floor(window/2).
      integer, intent(in) :: min_size !! Minimum per-side adaptive radius carried from the associated kza fit.
      real(dp), intent(in), optional :: tol !! Absolute derivative threshold treated as zero; defaults to 1.0e-5.
      real(dp), allocatable :: sigma(:)
      real(dp), allocatable :: d(:)
      real(dp), allocatable :: dprime(:)
      real(dp) :: avg
      real(dp) :: eps
      real(dp) :: scale
      real(dp) :: ss
      integer :: i
      integer :: q
      integer :: qh
      integer :: qt
      integer :: size_window
      integer :: t

      if (size(kza_data) /= size(kz_data)) error stop "kzsv: KZA and KZ data sizes differ"
      q = max(window / 2, 0)
      eps = 1.0e-5_dp
      if (present(tol)) eps = tol
      call differenced_1d(kz_data, q, d, dprime)
      scale = normalizer(d, 1.0_dp)
      allocate(sigma(size(kza_data)))

      do t = 1, size(kza_data)
         call choose_window(d(t), dprime(t), q, scale, eps, min_size, qt, qh)
         qt = min(qt, t - 1)
         qh = min(qh, size(kza_data) - t)
         size_window = qt + qh + 1
         if (size_window <= 1) then
            sigma(t) = quiet_nan()
         else
            avg = sum(kza_data(t - qt:t + qh)) / real(size_window, dp)
            ss = 0.0_dp
            do i = t - qt, t + qh
               ss = ss + (kza_data(i) - avg) ** 2
            end do
            sigma(t) = ss / real(size_window - 1, dp)
         end if
      end do
   end function kzsv

   pure function kza_core_1d(x, baseline, q, niter, min_window, eps, nprob) result(filtered)
      real(dp), intent(in) :: x(:) !! Raw one-dimensional data iteratively averaged inside adaptive windows.
      real(dp), intent(in) :: baseline(:) !! Fixed KZ baseline from which adaptive differences and derivatives are computed.
      integer, intent(in) :: q !! Maximum per-side adaptive radius.
      integer, intent(in) :: niter !! Number of repeated adaptive averaging iterations.
      integer, intent(in) :: min_window !! Minimum allowed per-side adaptive radius.
      real(dp), intent(in) :: eps !! Absolute derivative threshold treated as a flat difference metric.
      real(dp), intent(in) :: nprob !! Difference normalizer probability: 1 for maximum or 0.99 for robust quantile mode.
      real(dp), allocatable :: filtered(:)
      real(dp), allocatable :: d(:)
      real(dp), allocatable :: dprime(:)
      real(dp), allocatable :: work(:)
      real(dp) :: scale
      integer :: iter
      integer :: qh
      integer :: qt
      integer :: t

      call differenced_1d(baseline, q, d, dprime)
      scale = normalizer(d, nprob)
      allocate(work(size(x)), filtered(size(x)))
      work = x
      filtered = x
      do iter = 1, max(niter, 0)
         do t = 1, size(x)
            call choose_window(d(t), dprime(t), q, scale, eps, min_window, qt, qh)
            qt = min(qt, t - 1)
            qh = min(qh, size(x) - t)
            filtered(t) = finite_mean_1d(work, t - qt, t + qh)
         end do
         work = filtered
      end do
   end function kza_core_1d

   pure function kza_core_2d(x, baseline, q1, q2, niter, min_window, eps, nprob) result(filtered)
      real(dp), intent(in) :: x(:, :) !! Raw two-dimensional field iteratively averaged inside adaptive boxes.
      real(dp), intent(in) :: baseline(:, :) !! Fixed KZ baseline used to compute directional difference metrics.
      integer, intent(in) :: q1 !! Maximum per-side radius along the first array dimension.
      integer, intent(in) :: q2 !! Maximum per-side radius along the second array dimension.
      integer, intent(in) :: niter !! Number of repeated adaptive averaging iterations.
      integer, intent(in) :: min_window !! Minimum allowed per-side adaptive radius in each direction.
      real(dp), intent(in) :: eps !! Absolute derivative threshold treated as a flat difference metric.
      real(dp), intent(in) :: nprob !! Difference normalizer probability: 1 for maximum or 0.99 for robust quantile mode.
      real(dp), allocatable :: filtered(:, :)
      real(dp), allocatable :: d1(:, :)
      real(dp), allocatable :: d2(:, :)
      real(dp), allocatable :: dp1(:, :)
      real(dp), allocatable :: dp2(:, :)
      real(dp), allocatable :: work(:, :)
      real(dp) :: scale1
      real(dp) :: scale2
      integer :: h1
      integer :: h2
      integer :: i
      integer :: iter
      integer :: j
      integer :: t1
      integer :: t2

      call differences_2d(baseline, q1, q2, d1, d2, dp1, dp2)
      scale1 = normalizer(reshape(d1, [size(d1)]), nprob)
      scale2 = normalizer(reshape(d2, [size(d2)]), nprob)
      allocate(work(size(x, 1), size(x, 2)), filtered(size(x, 1), size(x, 2)))
      work = x
      filtered = x
      do iter = 1, max(niter, 0)
         do j = 1, size(x, 2)
            do i = 1, size(x, 1)
               call choose_window(d1(i, j), dp1(i, j), q1, scale1, eps, min_window, t1, h1)
               call choose_window(d2(i, j), dp2(i, j), q2, scale2, eps, min_window, t2, h2)
               filtered(i, j) = finite_mean_2d(work, i, j, t1, h1, t2, h2)
            end do
         end do
         work = filtered
      end do
   end function kza_core_2d

   pure function kza_core_3d(x, baseline, q1, q2, q3, niter, min_window, eps, nprob) result(filtered)
      real(dp), intent(in) :: x(:, :, :) !! Raw three-dimensional field iteratively averaged inside adaptive boxes.
      real(dp), intent(in) :: baseline(:, :, :) !! Fixed KZ baseline used to compute directional difference metrics.
      integer, intent(in) :: q1 !! Maximum per-side radius along the first array dimension.
      integer, intent(in) :: q2 !! Maximum per-side radius along the second array dimension.
      integer, intent(in) :: q3 !! Maximum per-side radius along the third array dimension.
      integer, intent(in) :: niter !! Number of repeated adaptive averaging iterations.
      integer, intent(in) :: min_window !! Minimum allowed per-side adaptive radius in each direction.
      real(dp), intent(in) :: eps !! Absolute derivative threshold treated as a flat difference metric.
      real(dp), intent(in) :: nprob !! Difference normalizer probability: 1 for maximum or 0.99 for robust quantile mode.
      real(dp), allocatable :: filtered(:, :, :)
      real(dp), allocatable :: d1(:, :, :)
      real(dp), allocatable :: d2(:, :, :)
      real(dp), allocatable :: d3(:, :, :)
      real(dp), allocatable :: dp1(:, :, :)
      real(dp), allocatable :: dp2(:, :, :)
      real(dp), allocatable :: dp3(:, :, :)
      real(dp), allocatable :: work(:, :, :)
      real(dp) :: scale1
      real(dp) :: scale2
      real(dp) :: scale3
      integer :: h1
      integer :: h2
      integer :: h3
      integer :: i
      integer :: iter
      integer :: j
      integer :: l
      integer :: t1
      integer :: t2
      integer :: t3

      call differences_3d(baseline, q1, q2, q3, d1, d2, d3, dp1, dp2, dp3)
      scale1 = normalizer(reshape(d1, [size(d1)]), nprob)
      scale2 = normalizer(reshape(d2, [size(d2)]), nprob)
      scale3 = normalizer(reshape(d3, [size(d3)]), nprob)
      allocate(work(size(x, 1), size(x, 2), size(x, 3)))
      allocate(filtered(size(x, 1), size(x, 2), size(x, 3)))
      work = x
      filtered = x
      do iter = 1, max(niter, 0)
         do l = 1, size(x, 3)
            do j = 1, size(x, 2)
               do i = 1, size(x, 1)
                  call choose_window(d1(i, j, l), dp1(i, j, l), q1, scale1, eps, min_window, t1, h1)
                  call choose_window(d2(i, j, l), dp2(i, j, l), q2, scale2, eps, min_window, t2, h2)
                  call choose_window(d3(i, j, l), dp3(i, j, l), q3, scale3, eps, min_window, t3, h3)
                  filtered(i, j, l) = finite_mean_3d(work, i, j, l, t1, h1, t2, h2, t3, h3)
               end do
            end do
         end do
         work = filtered
      end do
   end function kza_core_3d

   pure subroutine differenced_1d(y, q, d, dprime)
      real(dp), intent(in) :: y(:) !! Baseline signal used to compute the adaptive absolute-difference metric.
      integer, intent(in) :: q !! Per-side difference radius; edges are clamped to the available data.
      real(dp), allocatable, intent(out) :: d(:) !! Absolute difference metric |y(t+q)-y(t-q)| with clamped endpoints.
      real(dp), allocatable, intent(out) :: dprime(:) !! Forward first difference of d, with the final value copied from its predecessor.
      integer :: left
      integer :: n
      integer :: right
      integer :: t

      n = size(y)
      allocate(d(n), dprime(n))
      if (n == 0) return
      do t = 1, n
         left = max(1, t - q)
         right = min(n, t + q)
         d(t) = abs(y(right) - y(left))
      end do
      if (n == 1) then
         dprime(1) = 0.0_dp
      else
         do t = 1, n - 1
            dprime(t) = d(t + 1) - d(t)
         end do
         dprime(n) = dprime(n - 1)
      end if
   end subroutine differenced_1d

   pure subroutine differences_2d(y, q1, q2, d1, d2, dp1, dp2)
      real(dp), intent(in) :: y(:, :) !! Two-dimensional KZ baseline used for directional adaptive metrics.
      integer, intent(in) :: q1 !! Per-side radius for differences along the first array dimension.
      integer, intent(in) :: q2 !! Per-side radius for differences along the second array dimension.
      real(dp), allocatable, intent(out) :: d1(:, :) !! Absolute difference metric along the first array dimension.
      real(dp), allocatable, intent(out) :: d2(:, :) !! Absolute difference metric along the second array dimension.
      real(dp), allocatable, intent(out) :: dp1(:, :) !! Forward derivative of d1 along the first array dimension.
      real(dp), allocatable, intent(out) :: dp2(:, :) !! Forward derivative of d2 along the second array dimension.
      integer :: i
      integer :: j
      integer :: lo
      integer :: hi
      integer :: n1
      integer :: n2

      n1 = size(y, 1)
      n2 = size(y, 2)
      allocate(d1(n1, n2), d2(n1, n2), dp1(n1, n2), dp2(n1, n2))
      do j = 1, n2
         do i = 1, n1
            lo = max(1, i - q1)
            hi = min(n1, i + q1)
            d1(i, j) = abs(y(hi, j) - y(lo, j))
            lo = max(1, j - q2)
            hi = min(n2, j + q2)
            d2(i, j) = abs(y(i, hi) - y(i, lo))
         end do
      end do

      if (n1 <= 1) then
         dp1 = 0.0_dp
      else
         do j = 1, n2
            do i = 1, n1 - 1
               dp1(i, j) = d1(i + 1, j) - d1(i, j)
            end do
            dp1(n1, j) = dp1(n1 - 1, j)
         end do
      end if
      if (n2 <= 1) then
         dp2 = 0.0_dp
      else
         do i = 1, n1
            do j = 1, n2 - 1
               dp2(i, j) = d2(i, j + 1) - d2(i, j)
            end do
            dp2(i, n2) = dp2(i, n2 - 1)
         end do
      end if
   end subroutine differences_2d

   pure subroutine differences_3d(y, q1, q2, q3, d1, d2, d3, dp1, dp2, dp3)
      real(dp), intent(in) :: y(:, :, :) !! Three-dimensional KZ baseline used for directional adaptive metrics.
      integer, intent(in) :: q1 !! Per-side radius for differences along the first array dimension.
      integer, intent(in) :: q2 !! Per-side radius for differences along the second array dimension.
      integer, intent(in) :: q3 !! Per-side radius for differences along the third array dimension.
      real(dp), allocatable, intent(out) :: d1(:, :, :) !! Absolute difference metric along the first array dimension.
      real(dp), allocatable, intent(out) :: d2(:, :, :) !! Absolute difference metric along the second array dimension.
      real(dp), allocatable, intent(out) :: d3(:, :, :) !! Absolute difference metric along the third array dimension.
      real(dp), allocatable, intent(out) :: dp1(:, :, :) !! Forward derivative of d1 along the first array dimension.
      real(dp), allocatable, intent(out) :: dp2(:, :, :) !! Forward derivative of d2 along the second array dimension.
      real(dp), allocatable, intent(out) :: dp3(:, :, :) !! Forward derivative of d3 along the third array dimension.
      integer :: i
      integer :: j
      integer :: l
      integer :: lo
      integer :: hi
      integer :: n1
      integer :: n2
      integer :: n3

      n1 = size(y, 1)
      n2 = size(y, 2)
      n3 = size(y, 3)
      allocate(d1(n1, n2, n3), d2(n1, n2, n3), d3(n1, n2, n3))
      allocate(dp1(n1, n2, n3), dp2(n1, n2, n3), dp3(n1, n2, n3))
      do l = 1, n3
         do j = 1, n2
            do i = 1, n1
               lo = max(1, i - q1)
               hi = min(n1, i + q1)
               d1(i, j, l) = abs(y(hi, j, l) - y(lo, j, l))
               lo = max(1, j - q2)
               hi = min(n2, j + q2)
               d2(i, j, l) = abs(y(i, hi, l) - y(i, lo, l))
               lo = max(1, l - q3)
               hi = min(n3, l + q3)
               d3(i, j, l) = abs(y(i, j, hi) - y(i, j, lo))
            end do
         end do
      end do

      if (n1 <= 1) then
         dp1 = 0.0_dp
      else
         do l = 1, n3
            do j = 1, n2
               do i = 1, n1 - 1
                  dp1(i, j, l) = d1(i + 1, j, l) - d1(i, j, l)
               end do
               dp1(n1, j, l) = dp1(n1 - 1, j, l)
            end do
         end do
      end if
      if (n2 <= 1) then
         dp2 = 0.0_dp
      else
         do l = 1, n3
            do i = 1, n1
               do j = 1, n2 - 1
                  dp2(i, j, l) = d2(i, j + 1, l) - d2(i, j, l)
               end do
               dp2(i, n2, l) = dp2(i, n2 - 1, l)
            end do
         end do
      end if
      if (n3 <= 1) then
         dp3 = 0.0_dp
      else
         do j = 1, n2
            do i = 1, n1
               do l = 1, n3 - 1
                  dp3(i, j, l) = d3(i, j, l + 1) - d3(i, j, l)
               end do
               dp3(i, j, n3) = dp3(i, j, n3 - 1)
            end do
         end do
      end if
   end subroutine differences_3d

   pure subroutine choose_window(d, dprime, q, scale, eps, min_window, tail, head)
      real(dp), intent(in) :: d !! Local absolute-difference metric at the current cell.
      real(dp), intent(in) :: dprime !! Forward derivative of the local difference metric.
      integer, intent(in) :: q !! Maximum per-side radius for the current direction.
      real(dp), intent(in) :: scale !! Global maximum or robust quantile used in the adaptive shrink factor.
      real(dp), intent(in) :: eps !! Derivative magnitude below which the metric is treated as flat.
      integer, intent(in) :: min_window !! Minimum per-side adaptive radius.
      integer, intent(out) :: tail !! Selected radius toward decreasing array indices.
      integer, intent(out) :: head !! Selected radius toward increasing array indices.
      integer :: shrink

      shrink = floor(real(q, dp) * adaptive_factor(d, scale))
      if (abs(dprime) < eps) then
         head = shrink
         tail = shrink
      else if (dprime < 0.0_dp) then
         head = q
         tail = shrink
      else
         head = shrink
         tail = q
      end if
      tail = max(tail, min_window)
      head = max(head, min_window)
   end subroutine choose_window

   pure real(dp) function finite_mean_1d(x, first, last) result(avg)
      real(dp), intent(in) :: x(:) !! Input signal whose finite values are averaged over the requested inclusive interval.
      integer, intent(in) :: first !! First inclusive array index of the averaging interval.
      integer, intent(in) :: last !! Last inclusive array index of the averaging interval.
      integer :: i
      integer :: nfinite
      real(dp) :: total

      total = 0.0_dp
      nfinite = 0
      do i = max(first, 1), min(last, size(x))
         if (ieee_is_finite(x(i))) then
            total = total + x(i)
            nfinite = nfinite + 1
         end if
      end do
      if (nfinite == 0) then
         avg = quiet_nan()
      else
         avg = total / real(nfinite, dp)
      end if
   end function finite_mean_1d

   pure real(dp) function finite_mean_2d(x, center1, center2, tail1, head1, tail2, head2) result(avg)
      real(dp), intent(in) :: x(:, :) !! Input field whose finite cells are averaged over the adaptive rectangle.
      integer, intent(in) :: center1 !! Center index along the first array dimension.
      integer, intent(in) :: center2 !! Center index along the second array dimension.
      integer, intent(in) :: tail1 !! Radius toward lower indices along the first dimension.
      integer, intent(in) :: head1 !! Radius toward higher indices along the first dimension.
      integer, intent(in) :: tail2 !! Radius toward lower indices along the second dimension.
      integer, intent(in) :: head2 !! Radius toward higher indices along the second dimension.
      integer :: i
      integer :: j
      integer :: nfinite
      real(dp) :: total

      total = 0.0_dp
      nfinite = 0
      do j = max(1, center2 - tail2), min(size(x, 2), center2 + head2)
         do i = max(1, center1 - tail1), min(size(x, 1), center1 + head1)
            if (ieee_is_finite(x(i, j))) then
               total = total + x(i, j)
               nfinite = nfinite + 1
            end if
         end do
      end do
      if (nfinite == 0) then
         avg = quiet_nan()
      else
         avg = total / real(nfinite, dp)
      end if
   end function finite_mean_2d

   pure real(dp) function finite_mean_3d(x, center1, center2, center3, tail1, head1, tail2, head2, tail3, head3) result(avg)
      real(dp), intent(in) :: x(:, :, :) !! Input field whose finite cells are averaged over the adaptive rectangular prism.
      integer, intent(in) :: center1 !! Center index along the first array dimension.
      integer, intent(in) :: center2 !! Center index along the second array dimension.
      integer, intent(in) :: center3 !! Center index along the third array dimension.
      integer, intent(in) :: tail1 !! Radius toward lower indices along the first dimension.
      integer, intent(in) :: head1 !! Radius toward higher indices along the first dimension.
      integer, intent(in) :: tail2 !! Radius toward lower indices along the second dimension.
      integer, intent(in) :: head2 !! Radius toward higher indices along the second dimension.
      integer, intent(in) :: tail3 !! Radius toward lower indices along the third dimension.
      integer, intent(in) :: head3 !! Radius toward higher indices along the third dimension.
      integer :: i
      integer :: j
      integer :: l
      integer :: nfinite
      real(dp) :: total

      total = 0.0_dp
      nfinite = 0
      do l = max(1, center3 - tail3), min(size(x, 3), center3 + head3)
         do j = max(1, center2 - tail2), min(size(x, 2), center2 + head2)
            do i = max(1, center1 - tail1), min(size(x, 1), center1 + head1)
               if (ieee_is_finite(x(i, j, l))) then
                  total = total + x(i, j, l)
                  nfinite = nfinite + 1
               end if
            end do
         end do
      end do
      if (nfinite == 0) then
         avg = quiet_nan()
      else
         avg = total / real(nfinite, dp)
      end if
   end function finite_mean_3d

   pure function rotate_n_2d(a, rotations) result(b)
      real(dp), intent(in) :: a(:, :) !! Matrix rotated by a multiple of 90 degrees using the upstream .rot90 orientation.
      integer, intent(in) :: rotations !! Number of 90-degree counterclockwise rotations modulo four.
      real(dp), allocatable :: b(:, :)
      integer :: i
      integer :: j
      integer :: r

      r = modulo(rotations, 4)
      select case (r)
      case (0)
         allocate(b(size(a, 1), size(a, 2)))
         b = a
      case (1)
         allocate(b(size(a, 2), size(a, 1)))
         do j = 1, size(b, 2)
            do i = 1, size(b, 1)
               b(i, j) = a(j, size(a, 2) - i + 1)
            end do
         end do
      case (2)
         allocate(b(size(a, 1), size(a, 2)))
         do j = 1, size(b, 2)
            do i = 1, size(b, 1)
               b(i, j) = a(size(a, 1) - i + 1, size(a, 2) - j + 1)
            end do
         end do
      case (3)
         allocate(b(size(a, 2), size(a, 1)))
         do j = 1, size(b, 2)
            do i = 1, size(b, 1)
               b(i, j) = a(size(a, 1) - j + 1, i)
            end do
         end do
      end select
   end function rotate_n_2d

   pure real(dp) function normalization_probability(normalize) result(prob)
      character(len=*), intent(in), optional :: normalize !! R-style normalizer name; values beginning with q select the 0.99 quantile.

      prob = 1.0_dp
      if (present(normalize)) then
         if (len_trim(normalize) > 0) then
            if (normalize(1:1) == "q" .or. normalize(1:1) == "Q") prob = 0.99_dp
         end if
      end if
   end function normalization_probability

end module kza_filters
