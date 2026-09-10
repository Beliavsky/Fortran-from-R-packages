! SPDX-License-Identifier: GPL-2.0-or-later
! Statistical and scalar/vector utilities translated from wavethresh 4.7.3.
module wavethresh_stats
   use wavethresh_types, only : dp, wd_t
   implicit none
   private

   public :: l2norm, linfnorm, ssq, logabs, mad_value, madmad
   public :: shannon_entropy, sure, newsure, dof, robust_covariance
   public :: firstdot, levarr, c2to4, lt_to_name
   public :: shrink_soft, shrink_hard, kolsmi_chi2

contains

   pure function ssq(first, second) result(value)
      real(dp), intent(in) :: first(:) !! First vector in the squared-error comparison.
      real(dp), intent(in) :: second(:) !! Second vector; must have the same length as first.
      real(dp) :: value
      if (size(first) /= size(second)) then
         value = huge(1.0_dp)
      else
         value = sum((first - second)**2)
      end if
   end function ssq

   pure function l2norm(first, second) result(value)
      real(dp), intent(in) :: first(:) !! First vector in the Euclidean-distance calculation.
      real(dp), intent(in) :: second(:) !! Second vector; must have the same length as first.
      real(dp) :: value
      value = sqrt(ssq(first, second))
   end function l2norm

   pure function linfnorm(first, second) result(value)
      real(dp), intent(in) :: first(:) !! First vector in the maximum-distance calculation.
      real(dp), intent(in) :: second(:) !! Second vector; must have the same length as first.
      real(dp) :: value
      if (size(first) /= size(second) .or. size(first) == 0) then
         value = huge(1.0_dp)
      else
         value = maxval(abs(first - second))
      end if
   end function linfnorm

   elemental function logabs(x) result(value)
      real(dp), intent(in) :: x !! Real value for which wavethresh computes log2(x**2).
      real(dp) :: value
      value = log(x * x) / log(2.0_dp)
   end function logabs

   pure function median_value(x) result(value)
      real(dp), intent(in) :: x(:) !! Sample whose ordinary median is requested.
      real(dp) :: value
      real(dp), allocatable :: work(:)
      real(dp) :: key
      integer :: i
      integer :: j
      integer :: n
      n = size(x)
      if (n == 0) then
         value = 0.0_dp
         return
      end if
      work = x
      do i = 2, n
         key = work(i)
         j = i - 1
         do while (j >= 1)
            if (work(j) <= key) exit
            work(j + 1) = work(j)
            j = j - 1
         end do
         work(j + 1) = key
      end do
      if (mod(n, 2) == 1) then
         value = work((n + 1) / 2)
      else
         value = 0.5_dp * (work(n / 2) + work(n / 2 + 1))
      end if
   end function median_value

   pure function mad_value(x) result(value)
      real(dp), intent(in) :: x(:) !! Sample for R-compatible median absolute deviation with default scaling.
      real(dp) :: value
      real(dp) :: center
      if (size(x) == 0) then
         value = 0.0_dp
         return
      end if
      center = median_value(x)
      value = 1.482602218505602_dp * median_value(abs(x - center))
   end function mad_value

   pure function madmad(x) result(value)
      real(dp), intent(in) :: x(:) !! Sample whose squared median absolute deviation is requested.
      real(dp) :: value
      value = mad_value(x)**2
   end function madmad

   pure function shannon_entropy(values, zilchtol) result(entropy)
      real(dp), intent(in) :: values(:) !! Coefficient vector entering wavethresh Shannon.entropy.
      real(dp), intent(in), optional :: zilchtol !! Squared-energy tolerance below which entropy is zero.
      real(dp) :: entropy
      real(dp) :: tolerance
      real(dp), allocatable :: square(:)
      integer :: i
      tolerance = 1.0e-300_dp
      if (present(zilchtol)) tolerance = zilchtol
      square = values**2
      if (sum(square) < tolerance) then
         entropy = 0.0_dp
         return
      end if
      entropy = 0.0_dp
      do i = 1, size(square)
         if (square(i) > 0.0_dp) entropy = entropy - square(i) * log(square(i))
      end do
   end function shannon_entropy

   pure function sure(x) result(threshold)
      real(dp), intent(in) :: x(:) !! Noise-standardized coefficients for Donoho-Johnstone SURE threshold selection.
      real(dp) :: threshold
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: cumulative(:)
      real(dp) :: best
      real(dp) :: risk
      real(dp) :: key
      integer :: d
      integer :: i
      integer :: j
      d = size(x)
      if (d == 0) then
         threshold = 0.0_dp
         return
      end if
      y = abs(x)
      do i = 2, d
         key = y(i)
         j = i - 1
         do while (j >= 1)
            if (y(j) <= key) exit
            y(j + 1) = y(j)
            j = j - 1
         end do
         y(j + 1) = key
      end do
      allocate(cumulative(0:d))
      cumulative(0) = 0.0_dp
      do i = 1, d
         cumulative(i) = cumulative(i - 1) + y(i)**2
      end do
      best = huge(1.0_dp)
      threshold = y(1)
      do i = 1, d
         risk = real(d - 2 * i, dp) + cumulative(i - 1) + real(d - i + 1, dp) * y(i)**2
         if (risk < best) then
            best = risk
            threshold = y(i)
         end if
      end do
   end function sure

   pure function newsure(sigma, x) result(threshold)
      real(dp), intent(in) :: sigma(:) !! Per-coefficient standard deviations corresponding to x.
      real(dp), intent(in) :: x(:) !! Coefficients for heteroscedastic SURE threshold selection.
      real(dp) :: threshold
      integer, allocatable :: order(:)
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: s(:)
      real(dp) :: best
      real(dp) :: risk
      integer :: d
      integer :: i
      integer :: j
      integer :: key
      d = size(x)
      if (d == 0 .or. size(sigma) /= d) then
         threshold = 0.0_dp
         return
      end if
      allocate(order(d))
      do i = 1, d
         order(i) = i
      end do
      do i = 2, d
         key = order(i)
         j = i - 1
         do while (j >= 1)
            if (abs(x(order(j))) <= abs(x(key))) exit
            order(j + 1) = order(j)
            j = j - 1
         end do
         order(j + 1) = key
      end do
      allocate(y(d), s(d))
      do i = 1, d
         y(i) = abs(x(order(i)))
         s(i) = sigma(order(i))
      end do
      best = huge(1.0_dp)
      threshold = minval(abs(x))
      do i = 1, d
         risk = real(d, dp)
         if (i > 1) then
            risk = risk - 2.0_dp * sum(s(1:i - 1)**2) + sum(y(1:i - 1)**2)
         end if
         risk = risk + real(d - i + 1, dp) * y(i)**2
         if (risk < best) then
            best = risk
            threshold = y(i)
         end if
      end do
   end function newsure

   pure function dof(object) result(count)
      type(wd_t), intent(in) :: object !! One-dimensional wavelet object whose nonzero coefficient count is requested.
      integer :: count
      integer :: level
      count = 1
      if (.not. object%ok) then
         count = 0
         return
      end if
      do level = 0, object%nlevels - 1
         count = count + count_nonzero(object%detail(level)%values)
      end do
   end function dof

   pure function count_nonzero(values) result(n_nonzero)
      real(dp), intent(in) :: values(:) !! Coefficients to count using exact zero semantics of upstream dof().
      integer :: n_nonzero
      n_nonzero = count(abs(values) > 0.0_dp)
   end function count_nonzero

   function robust_covariance(x) result(sigma)
      real(dp), intent(in) :: x(:,:) !! Matrix with variables in rows and observations in columns, matching rcov().
      real(dp), allocatable :: sigma(:,:)
      real(dp), allocatable :: a(:)
      real(dp) :: b1
      real(dp) :: b2
      real(dp) :: b3
      integer :: i
      integer :: j
      integer :: m
      m = size(x, 1)
      allocate(sigma(m, m), a(m))
      sigma = 0.0_dp
      do i = 1, m
         if (mad_value(x(i, :)) <= tiny(1.0_dp)) then
            a(i) = 0.0_dp
            sigma(i, i) = 0.0_dp
         else
            a(i) = 1.0_dp / mad_value(x(i, :))
            sigma(i, i) = 1.0_dp / a(i)**2
         end if
      end do
      do i = 2, m
         do j = 1, i - 1
            if (abs(a(i)) <= tiny(1.0_dp)) cycle
            if (abs(a(j)) <= tiny(1.0_dp)) cycle
            b1 = mad_value(a(i) * x(i, :) + a(j) * x(j, :))**2
            b2 = mad_value(a(i) * x(i, :) - a(j) * x(j, :))**2
            b3 = mad_value(a(j) * x(j, :) - a(i) * x(i, :))**2
            if (b1 + b2 > tiny(1.0_dp)) sigma(i, j) = &
               (b1 - b2) / ((b1 + b2) * a(i) * a(j))
            if (b1 + b3 > tiny(1.0_dp)) sigma(j, i) = &
               (b1 - b3) / ((b1 + b3) * a(i) * a(j))
         end do
      end do
   end function robust_covariance

   function firstdot(strings) result(position)
      character(len=*), intent(in) :: strings(:) !! Strings in which to locate the first period character.
      integer, allocatable :: position(:)
      integer :: i
      allocate(position(size(strings)))
      do i = 1, size(strings)
         position(i) = index(strings(i), ".")
      end do
   end function firstdot

   recursive function levarr(values, levels) result(reordered)
      real(dp), intent(in) :: values(:) !! Input vector to recursively split into odd and even positions.
      integer, intent(in) :: levels !! Number of recursive odd/even rearrangement stages.
      real(dp), allocatable :: reordered(:)
      real(dp), allocatable :: odd(:)
      real(dp), allocatable :: even(:)
      real(dp), allocatable :: rodd(:)
      real(dp), allocatable :: reven(:)
      if (levels <= 0 .or. size(values) <= 1) then
         reordered = values
         return
      end if
      odd = values(1:size(values):2)
      even = values(2:size(values):2)
      rodd = levarr(odd, levels - 1)
      reven = levarr(even, levels - 1)
      reordered = [rodd, reven]
   end function levarr

   pure function c2to4(index_value) result(answer)
      integer, intent(in) :: index_value !! Nonnegative integer interpreted through its base-two digits as base-four digits.
      integer :: answer
      integer :: bit
      integer :: place
      integer :: work
      answer = 0
      place = 1
      work = max(index_value, 0)
      do while (work > 0)
         bit = mod(work, 2)
         answer = answer + bit * place
         place = place * 4
         work = work / 2
      end do
   end function c2to4

   pure function lt_to_name(level, coefficient_type) result(name)
      integer, intent(in) :: level !! Two-dimensional wavelet level included in the wavethresh coefficient name.
      character(len=*), intent(in) :: coefficient_type !! Orientation code CD, DC, DD, or CC.
      character(len=32) :: name
      character(len=1) :: digit
      select case (trim(coefficient_type))
      case ("CD")
         digit = "1"
      case ("DC")
         digit = "2"
      case ("DD")
         digit = "3"
      case ("CC")
         digit = "4"
      case default
         digit = "?"
      end select
      write(name, '("w", i0, "L", a)') level, digit
   end function lt_to_name

   elemental function shrink_soft(coefficient, threshold) result(value)
      real(dp), intent(in) :: coefficient !! Coefficient to soft-threshold.
      real(dp), intent(in) :: threshold !! Nonnegative shrinkage threshold.
      real(dp) :: value
      value = sign(max(abs(coefficient) - max(threshold, 0.0_dp), 0.0_dp), coefficient)
   end function shrink_soft

   elemental function shrink_hard(coefficient, threshold) result(value)
      real(dp), intent(in) :: coefficient !! Coefficient to hard-threshold.
      real(dp), intent(in) :: threshold !! Nonnegative hard threshold.
      real(dp) :: value
      if (abs(coefficient) <= max(threshold, 0.0_dp)) then
         value = 0.0_dp
      else
         value = coefficient
      end if
   end function shrink_hard

   pure function kolsmi_chi2(data) result(value)
      real(dp), intent(in) :: data(:) !! Sequence entering wavethresh TOkolsmi.chi2 cumulative-deviation statistic.
      real(dp) :: value
      real(dp) :: total
      real(dp) :: cumulative
      integer :: i
      integer :: n
      n = size(data)
      if (n == 0) then
         value = 0.0_dp
         return
      end if
      total = sum(data)
      cumulative = 0.0_dp
      value = 0.0_dp
      do i = 1, n
         cumulative = cumulative + data(i)
         value = max(value, abs(cumulative - real(i, dp) * total / real(n, dp)))
      end do
      value = value / sqrt(2.0_dp * real(n, dp))
   end function kolsmi_chi2

end module wavethresh_stats
