! SPDX-License-Identifier: GPL-3.0-only
! Derived from sharpeRratio 1.4.3 by Damien Challet.
module sharpe_rratio_records
   use, intrinsic :: iso_fortran_env, only : int64
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use ghyp_kinds, only : dp
   implicit none
   private

   type, public :: r0_result
      real(dp) :: mean = 0.0_dp
      real(dp) :: q1 = 0.0_dp
      real(dp) :: q2 = 0.0_dp
      integer :: n = 0
      integer :: num_permutations = 0
      logical :: ok = .false.
      character(len=160) :: message = ''
   end type r0_result

   type :: rng_state
      integer(int64) :: state = 1_int64
   end type rng_state

   public :: num_records_up, num_records_down, compute_r0bar
   public :: computeR0bar

   interface computeR0bar
      module procedure compute_r0bar
   end interface computeR0bar

contains

   pure function num_records_up(x) result(number)
      real(dp), intent(in) :: x(:)
      integer :: number
      real(dp) :: xmax
      integer :: i

      if (size(x) == 0) then
         number = 0
         return
      end if
      number = 1
      xmax = x(1)
      do i = 2, size(x)
         if (x(i) > xmax) then
            number = number+1
            xmax = x(i)
         end if
      end do
   end function num_records_up

   pure function num_records_down(x) result(number)
      real(dp), intent(in) :: x(:)
      integer :: number
      real(dp) :: xmin
      integer :: i

      if (size(x) == 0) then
         number = 0
         return
      end if
      number = 1
      xmin = x(1)
      do i = 2, size(x)
         if (x(i) < xmin) then
            number = number+1
            xmin = x(i)
         end if
      end do
   end function num_records_down

   subroutine initialize_rng(rng, seed)
      type(rng_state), intent(out) :: rng
      integer(int64), intent(in), optional :: seed
      integer(int64) :: raw
      integer :: count, rate, count_max

      if (present(seed)) then
         raw = seed
      else
         call system_clock(count,count_rate=rate,count_max=count_max)
         raw = int(count,int64)+104729_int64*int(max(1,rate),int64)
      end if
      raw = modulo(raw,2147483646_int64)+1_int64
      rng%state = raw
   end subroutine initialize_rng

   function uniform_random(rng) result(value)
      type(rng_state), intent(inout) :: rng
      real(dp) :: value
      integer(int64) :: hi, lo, test

      hi = rng%state/127773_int64
      lo = modulo(rng%state,127773_int64)
      test = 16807_int64*lo-2836_int64*hi
      if (test > 0_int64) then
         rng%state = test
      else
         rng%state = test+2147483647_int64
      end if
      value = real(rng%state,dp)/2147483647.0_dp
   end function uniform_random

   subroutine shuffle_values(x, rng)
      real(dp), intent(inout) :: x(:)
      type(rng_state), intent(inout) :: rng
      real(dp) :: temporary
      integer :: i, j

      do i = size(x), 2, -1
         j = 1+int(uniform_random(rng)*real(i,dp))
         j = min(i,max(1,j))
         temporary = x(i)
         x(i) = x(j)
         x(j) = temporary
      end do
   end subroutine shuffle_values

   pure subroutine sort_real(x)
      real(dp), intent(inout) :: x(:)
      real(dp) :: key
      integer :: i, j

      do i = 2, size(x)
         key = x(i)
         j = i-1
         do while (j >= 1)
            if (x(j) <= key) exit
            x(j+1) = x(j)
            j = j-1
         end do
         x(j+1) = key
      end do
   end subroutine sort_real

   pure function source_quantile(sorted, probability) result(value)
      real(dp), intent(in) :: sorted(:), probability
      real(dp) :: value
      integer :: index

      index = floor(probability*real(size(sorted),dp))+1
      index = min(size(sorted),max(1,index))
      value = sorted(index)
   end function source_quantile

   pure function type7_quantile(sorted, probability) result(value)
      real(dp), intent(in) :: sorted(:), probability
      real(dp) :: value, h, fraction
      integer :: lo, hi

      if (size(sorted) == 1) then
         value = sorted(1)
         return
      end if
      h = 1.0_dp+real(size(sorted)-1,dp)*min(1.0_dp,max(0.0_dp,probability))
      lo = floor(h)
      hi = ceiling(h)
      fraction = h-real(lo,dp)
      value = (1.0_dp-fraction)*sorted(lo)+fraction*sorted(hi)
   end function type7_quantile

   function compute_r0bar(x, num_perm, q1, q2, seed, source_compatible) result(result)
      real(dp), intent(in) :: x(:)
      integer, intent(in), optional :: num_perm
      real(dp), intent(in), optional :: q1, q2
      integer(int64), intent(in), optional :: seed
      logical, intent(in), optional :: source_compatible
      type(r0_result) :: result
      real(dp), allocatable :: work(:), cumulative(:), r0_values(:)
      real(dp) :: lower_probability, upper_probability, r0
      integer :: permutations, permutation, i
      logical :: source_mode
      type(rng_state) :: rng

      if (size(x) < 1) then
         result%message = 'at least one observation is required'
         return
      end if
      if (.not. all(ieee_is_finite(x))) then
         result%message = 'compute_r0bar requires finite observations'
         return
      end if
      permutations = 100
      if (present(num_perm)) permutations = num_perm
      if (permutations < 1) then
         result%message = 'num_perm must be positive'
         return
      end if
      lower_probability = 0.025_dp
      upper_probability = 0.975_dp
      if (present(q1)) lower_probability = q1
      if (present(q2)) upper_probability = q2
      if (lower_probability < 0.0_dp .or. lower_probability > 1.0_dp .or. &
          upper_probability < 0.0_dp .or. upper_probability > 1.0_dp) then
         result%message = 'quantile probabilities must be in [0,1]'
         return
      end if
      source_mode = .true.
      if (present(source_compatible)) source_mode = source_compatible

      allocate(work(size(x)),cumulative(size(x)),r0_values(permutations))
      work = x
      call initialize_rng(rng,seed)
      do permutation = 1, permutations
         cumulative(1) = work(1)
         do i = 2, size(work)
            cumulative(i) = cumulative(i-1)+work(i)
         end do
         r0 = real(num_records_up(cumulative)-num_records_down(cumulative),dp)
         r0_values(permutation) = r0
         if (permutation < permutations) call shuffle_values(work,rng)
      end do
      result%mean = sum(r0_values)/real(permutations,dp)
      call sort_real(r0_values)
      if (source_mode) then
         result%q1 = source_quantile(r0_values,lower_probability)
         result%q2 = source_quantile(r0_values,upper_probability)
      else
         result%q1 = type7_quantile(r0_values,lower_probability)
         result%q2 = type7_quantile(r0_values,upper_probability)
      end if
      if (result%q1 > result%q2) then
         r0 = result%q1
         result%q1 = result%q2
         result%q2 = r0
      end if
      result%n = size(x)
      result%num_permutations = permutations
      result%ok = .true.
   end function compute_r0bar

end module sharpe_rratio_records
