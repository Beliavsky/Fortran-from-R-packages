! SPDX-License-Identifier: GPL-2.0-or-later
module rf_rng
   use r_kinds, only : dp, i64
   implicit none
   private

   integer(i64), parameter :: pm_m = 2147483647_i64
   integer(i64), parameter :: pm_a = 16807_i64
   integer(i64), parameter :: pm_q = 127773_i64
   integer(i64), parameter :: pm_r = 2836_i64

   type, public :: rf_rng_state
      integer(i64) :: state = 1976_i64
   contains
      procedure :: seed => rf_seed
      procedure :: uniform => rf_uniform
      procedure :: randint => rf_randint
   end type rf_rng_state

   public :: sample_indices, shuffle_int, shuffle_real

contains

   subroutine rf_seed(self, seed)
      class(rf_rng_state), intent(inout) :: self
      integer(i64), intent(in) :: seed

      self%state = modulo(abs(seed), pm_m - 1_i64) + 1_i64
   end subroutine rf_seed

   function rf_uniform(self) result(u)
      class(rf_rng_state), intent(inout) :: self
      real(dp) :: u
      integer(i64) :: hi, lo, test

      hi = self%state / pm_q
      lo = modulo(self%state, pm_q)
      test = pm_a * lo - pm_r * hi
      if (test > 0_i64) then
         self%state = test
      else
         self%state = test + pm_m
      end if
      u = real(self%state, dp) / real(pm_m, dp)
   end function rf_uniform

   function rf_randint(self, n) result(k)
      class(rf_rng_state), intent(inout) :: self
      integer, intent(in) :: n
      integer :: k

      if (n <= 1) then
         k = 1
      else
         k = min(n, 1 + int(self%uniform() * real(n, dp)))
      end if
   end function rf_randint

   subroutine sample_indices(rng, population_size, sample_size, replace, indices, weights, status)
      type(rf_rng_state), intent(inout) :: rng
      integer, intent(in) :: population_size, sample_size
      logical, intent(in) :: replace
      integer, intent(out) :: indices(sample_size)
      real(dp), intent(in), optional :: weights(population_size)
      integer, intent(out), optional :: status
      real(dp), allocatable :: w(:), cumulative(:)
      integer, allocatable :: pool(:)
      integer :: i, j, last
      real(dp) :: s, u

      if (present(status)) status = 0
      if (population_size <= 0 .or. sample_size < 0) then
         if (present(status)) status = 1
         if (sample_size > 0) indices = 1
         return
      end if
      if (.not. replace .and. sample_size > population_size) then
         if (present(status)) status = 2
         indices = 1
         return
      end if

      if (.not. present(weights)) then
         if (replace) then
            do i = 1, sample_size
               indices(i) = rng%randint(population_size)
            end do
         else
            allocate(pool(population_size))
            pool = [(i, i = 1, population_size)]
            last = population_size
            do i = 1, sample_size
               j = rng%randint(last)
               indices(i) = pool(j)
               pool(j) = pool(last)
               last = last - 1
            end do
         end if
         return
      end if

      allocate(w(population_size), cumulative(population_size))
      w = max(weights, 0.0_dp)
      s = sum(w)
      if (s <= 0.0_dp) then
         if (present(status)) status = 3
         indices = 1
         return
      end if
      w = w / s

      if (replace) then
         call make_cumulative(w, cumulative)
         do i = 1, sample_size
            u = rng%uniform()
            indices(i) = locate_cumulative(cumulative, u)
         end do
      else
         do i = 1, sample_size
            s = sum(w)
            if (s <= 0.0_dp) then
               if (present(status)) status = 4
               indices(i:) = 1
               return
            end if
            w = w / s
            call make_cumulative(w, cumulative)
            u = rng%uniform()
            j = locate_cumulative(cumulative, u)
            indices(i) = j
            w(j) = 0.0_dp
         end do
      end if
   end subroutine sample_indices

   subroutine make_cumulative(w, cumulative)
      real(dp), intent(in) :: w(:)
      real(dp), intent(out) :: cumulative(size(w))
      integer :: i

      cumulative(1) = w(1)
      do i = 2, size(w)
         cumulative(i) = cumulative(i - 1) + w(i)
      end do
      cumulative(size(w)) = 1.0_dp
   end subroutine make_cumulative

   integer function locate_cumulative(cumulative, u) result(pos)
      real(dp), intent(in) :: cumulative(:), u
      integer :: lo, hi, mid

      lo = 1
      hi = size(cumulative)
      do while (lo < hi)
         mid = (lo + hi) / 2
         if (u <= cumulative(mid)) then
            hi = mid
         else
            lo = mid + 1
         end if
      end do
      pos = lo
   end function locate_cumulative

   subroutine shuffle_int(rng, x)
      type(rf_rng_state), intent(inout) :: rng
      integer, intent(inout) :: x(:)
      integer :: i, j, tmp

      do i = size(x), 2, -1
         j = rng%randint(i)
         tmp = x(i)
         x(i) = x(j)
         x(j) = tmp
      end do
   end subroutine shuffle_int

   subroutine shuffle_real(rng, x)
      type(rf_rng_state), intent(inout) :: rng
      real(dp), intent(inout) :: x(:)
      integer :: i, j
      real(dp) :: tmp

      do i = size(x), 2, -1
         j = rng%randint(i)
         tmp = x(i)
         x(i) = x(j)
         x(j) = tmp
      end do
   end subroutine shuffle_real

end module rf_rng
