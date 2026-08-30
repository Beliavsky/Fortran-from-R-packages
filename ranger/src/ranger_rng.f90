! SPDX-License-Identifier: GPL-3.0-only
! Modern Fortran translation of ranger 0.18.0 computational code.
! Upstream C++ core: Copyright (c) 2014-2018 Marvin N. Wright; MIT licensed.
! The upstream R package is GPL-3. See NOTICE.md for full provenance.
module ranger_rng
   use r_kinds, only : dp, i64
   implicit none
   private

   integer, parameter :: mt_n = 312
   integer, parameter :: mt_m = 156
   integer(i64), parameter :: mt_matrix_a = int(z'B5026F5AA96619E9', i64)
   integer(i64), parameter :: mt_upper_mask = int(z'FFFFFFFF80000000', i64)
   integer(i64), parameter :: mt_lower_mask = int(z'000000007FFFFFFF', i64)
   integer(i64), parameter :: mask32 = int(z'00000000FFFFFFFF', i64)
   integer(i64), parameter :: mask16 = int(z'000000000000FFFF', i64)
   integer(i64), parameter :: mt_init_mult = int(z'5851F42D4C957F2D', i64)
   real(dp), parameter :: two32 = 4294967296.0_dp
   real(dp), parameter :: two64 = 18446744073709551616.0_dp

   type, public :: ranger_rng_state
      integer(i64) :: mt(mt_n) = 0_i64
      integer :: index = mt_n + 1
   contains
      procedure :: seed => ranger_seed
      procedure :: uniform => ranger_uniform
      procedure :: randint => ranger_randint
      procedure :: randint64 => ranger_randint_i64
      procedure, private :: next_bits => ranger_next_bits
   end type ranger_rng_state

   public :: sample_indices, shuffle_int, shuffle_real

contains

   subroutine ranger_seed(self, seed)
      class(ranger_rng_state), intent(inout) :: self
      integer(i64), intent(in) :: seed
      integer :: i
      integer(i64) :: value

      value = seed
      self%mt(1) = value
      do i = 2, mt_n
         value = ieor(self%mt(i - 1), shiftr(self%mt(i - 1), 62))
         self%mt(i) = mul_add_mod64(value, mt_init_mult, int(i - 1, i64))
      end do
      self%index = mt_n
   end subroutine ranger_seed

   function ranger_next_bits(self) result(value)
      class(ranger_rng_state), intent(inout) :: self
      integer(i64) :: value
      integer(i64) :: x, xa
      integer :: i

      if (self%index >= mt_n) then
         do i = 1, mt_n
            x = ior(iand(self%mt(i), mt_upper_mask), iand(self%mt(mod(i, mt_n) + 1), mt_lower_mask))
            xa = shiftr(x, 1)
            if (btest(x, 0)) xa = ieor(xa, mt_matrix_a)
            self%mt(i) = ieor(self%mt(mod(i - 1 + mt_m, mt_n) + 1), xa)
         end do
         self%index = 0
      end if

      self%index = self%index + 1
      value = self%mt(self%index)
      value = ieor(value, iand(shiftr(value, 29), int(z'5555555555555555', i64)))
      value = ieor(value, iand(shiftl(value, 17), int(z'71D67FFFEDA60000', i64)))
      value = ieor(value, iand(shiftl(value, 37), int(z'FFF7EEE000000000', i64)))
      value = ieor(value, shiftr(value, 43))
   end function ranger_next_bits

   function ranger_uniform(self) result(u)
      class(ranger_rng_state), intent(inout) :: self
      real(dp) :: u
      integer(i64) :: bits, hi, lo

      bits = self%next_bits()
      hi = shiftr(bits, 32)
      lo = iand(bits, mask32)
      u = (real(hi, dp) * two32 + real(lo, dp)) / two64
      if (u >= 1.0_dp) u = nearest(1.0_dp, -1.0_dp)
   end function ranger_uniform

   function ranger_randint(self, n) result(k)
      class(ranger_rng_state), intent(inout) :: self
      integer, intent(in) :: n
      integer :: k

      if (n <= 1) then
         k = 1
      else
         k = int(self%randint64(int(n, i64)))
      end if
   end function ranger_randint

   function ranger_randint_i64(self, n) result(k)
      class(ranger_rng_state), intent(inout) :: self
      integer(i64), intent(in) :: n
      integer(i64) :: k
      integer(i64) :: bits, high, low, threshold

      if (n <= 1_i64) then
         k = 1_i64
         return
      end if

      ! libstdc++ uniform_int_distribution uses Lemire's multiply-high
      ! downscaling for a full-width mt19937_64 engine. Reproduce that
      ! algorithm without requiring a nonstandard 128-bit integer kind.
      threshold = pow2_64_mod(n)
      do
         bits = self%next_bits()
         call multiply_u64(bits, n, high, low)
         if (low < 0_i64 .or. low >= threshold) exit
      end do
      k = high + 1_i64
   end function ranger_randint_i64

   pure subroutine multiply_u64(a, b, high, low)
      integer(i64), intent(in) :: a, b
      integer(i64), intent(out) :: high, low
      integer(i64) :: a0, a1, b0, b1
      integer(i64) :: h0, l0, h1, l1, h2, l2, h3, l3
      integer(i64) :: middle, carry, high_middle, high_hi

      a0 = iand(a, mask32)
      a1 = shiftr(a, 32)
      b0 = iand(b, mask32)
      b1 = shiftr(b, 32)

      call mul32(a0, b0, h0, l0)
      call mul32(a0, b1, h1, l1)
      call mul32(a1, b0, h2, l2)
      call mul32(a1, b1, h3, l3)

      middle = h0 + l1 + l2
      carry = shiftr(middle, 32)
      low = ior(shiftl(iand(middle, mask32), 32), l0)

      high_middle = l3 + h1 + h2 + carry
      high_hi = iand(h3 + shiftr(high_middle, 32), mask32)
      high = ior(shiftl(high_hi, 32), iand(high_middle, mask32))
   end subroutine multiply_u64

   pure function pow2_64_mod(n) result(value)
      integer(i64), intent(in) :: n
      integer(i64) :: value
      integer :: i

      if (n <= 1_i64) then
         value = 0_i64
         return
      end if
      value = 1_i64
      do i = 1, 64
         if (value >= n - value) then
            value = value - (n - value)
         else
            value = 2_i64 * value
         end if
      end do
   end function pow2_64_mod

   pure function mul_add_mod64(a, b, addend) result(value)
      integer(i64), intent(in) :: a, b, addend
      integer(i64) :: value
      integer(i64) :: a_hi, a_lo, b_hi, b_lo
      integer(i64) :: p_hi, p_lo, cross1_lo, cross2_lo
      integer(i64) :: dummy_hi, hi, lo, carry, add_lo, add_hi

      a_lo = iand(a, mask32)
      a_hi = shiftr(a, 32)
      b_lo = iand(b, mask32)
      b_hi = shiftr(b, 32)

      call mul32(a_lo, b_lo, p_hi, p_lo)
      call mul32(a_hi, b_lo, dummy_hi, cross1_lo)
      call mul32(a_lo, b_hi, dummy_hi, cross2_lo)
      hi = iand(p_hi + cross1_lo + cross2_lo, mask32)
      lo = p_lo

      add_lo = iand(addend, mask32)
      add_hi = shiftr(addend, 32)
      lo = lo + add_lo
      carry = shiftr(lo, 32)
      lo = iand(lo, mask32)
      hi = iand(hi + add_hi + carry, mask32)
      value = ior(shiftl(hi, 32), lo)
   end function mul_add_mod64

   pure subroutine mul32(a, b, hi, lo)
      integer(i64), intent(in) :: a, b
      integer(i64), intent(out) :: hi, lo
      integer(i64) :: a0, a1, b0, b1, p0, p1, p2, work, carry

      a0 = iand(a, mask16)
      a1 = shiftr(a, 16)
      b0 = iand(b, mask16)
      b1 = shiftr(b, 16)
      p0 = a0 * b0
      p1 = a0 * b1 + a1 * b0
      p2 = a1 * b1
      work = p0 + shiftl(iand(p1, mask16), 16)
      carry = shiftr(work, 32)
      lo = iand(work, mask32)
      hi = iand(p2 + shiftr(p1, 16) + carry, mask32)
   end subroutine mul32

   subroutine sample_indices(rng, population_size, sample_size, replace, indices, weights, status)
      type(ranger_rng_state), intent(inout) :: rng
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
         else if (sample_size < population_size / 10) then
            allocate(pool(population_size))
            pool = 0
            do i = 1, sample_size
               do
                  j = rng%randint(population_size)
                  if (pool(j) == 0) exit
               end do
               pool(j) = 1
               indices(i) = j
            end do
         else
            allocate(pool(population_size))
            pool = [(i, i = 1, population_size)]
            do i = 1, sample_size
               j = i + int(rng%uniform() * real(population_size - i + 1, dp))
               if (j > population_size) j = population_size
               last = pool(i)
               pool(i) = pool(j)
               pool(j) = last
               indices(i) = pool(i)
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
         allocate(pool(population_size))
         pool = 0
         call make_cumulative(w, cumulative)
         do i = 1, sample_size
            do
               u = rng%uniform()
               j = locate_cumulative(cumulative, u)
               if (pool(j) == 0) exit
            end do
            pool(j) = 1
            indices(i) = j
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
      type(ranger_rng_state), intent(inout) :: rng
      integer, intent(inout) :: x(:)
      integer :: i, j1, j2, tmp
      integer(i64) :: combined, swap_range

      if (size(x) <= 1) return

      i = 2
      if (mod(size(x), 2) == 0) then
         j1 = int(rng%randint64(2_i64))
         tmp = x(i)
         x(i) = x(j1)
         x(j1) = tmp
         i = i + 1
      end if

      do while (i <= size(x))
         swap_range = int(i, i64)
         combined = rng%randint64(swap_range * (swap_range + 1_i64)) - 1_i64
         j1 = int(combined / (swap_range + 1_i64)) + 1
         j2 = int(modulo(combined, swap_range + 1_i64)) + 1

         tmp = x(i)
         x(i) = x(j1)
         x(j1) = tmp
         tmp = x(i + 1)
         x(i + 1) = x(j2)
         x(j2) = tmp
         i = i + 2
      end do
   end subroutine shuffle_int

   subroutine shuffle_real(rng, x)
      type(ranger_rng_state), intent(inout) :: rng
      real(dp), intent(inout) :: x(:)
      integer :: i, j1, j2
      integer(i64) :: combined, swap_range
      real(dp) :: tmp

      if (size(x) <= 1) return

      i = 2
      if (mod(size(x), 2) == 0) then
         j1 = int(rng%randint64(2_i64))
         tmp = x(i)
         x(i) = x(j1)
         x(j1) = tmp
         i = i + 1
      end if

      do while (i <= size(x))
         swap_range = int(i, i64)
         combined = rng%randint64(swap_range * (swap_range + 1_i64)) - 1_i64
         j1 = int(combined / (swap_range + 1_i64)) + 1
         j2 = int(modulo(combined, swap_range + 1_i64)) + 1

         tmp = x(i)
         x(i) = x(j1)
         x(j1) = tmp
         tmp = x(i + 1)
         x(i + 1) = x(j2)
         x(j2) = tmp
         i = i + 2
      end do
   end subroutine shuffle_real

end module ranger_rng
