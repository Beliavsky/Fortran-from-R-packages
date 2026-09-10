module gmp_number_theory
   use, intrinsic :: iso_fortran_env, only : int64
   use gmp_bigz, only : bigz, bigz_from_int64, bigz_to_int64, bigz_fits_int64, bigz_compare, bigz_equal, &
      bigz_is_zero, bigz_is_one, bigz_is_even, bigz_add, bigz_sub, bigz_neg, bigz_abs, bigz_mul, &
      bigz_mul_small, bigz_quotient, bigz_modulo, bigz_div_small, bigz_gcd, bigz_powm
   use gmp_bigq, only : bigq, bigq_make, bigq_from_int64, bigq_add, bigq_sub, bigq_mul, bigq_pow
   implicit none
   private

   public :: factorial_z, choose_z
   public :: fibnum_z, fibnum2_z, lucnum_z, lucnum2_z
   public :: stirling1_z, stirling1_all_z, stirling2_z, stirling2_all_z
   public :: eulerian_z, eulerian_all_z, bernoulli_q, dbinom_q
   public :: isprime_z, nextprime_z, factorize_z, urand_bigz

contains

   pure elemental function factorial_z(n) result(value)
      integer(int64), intent(in) :: n !! Nonnegative factorial argument; negative values return zero as an error sentinel.
      type(bigz) :: value
      integer(int64) :: k

      if (n < 0_int64) then
         value = bigz_from_int64(0_int64)
         return
      end if
      value = bigz_from_int64(1_int64)
      do k = 2_int64, n
         value = bigz_mul_small(value, k)
      end do
   end function factorial_z

   pure elemental function choose_z(n, k) result(value)
      type(bigz), intent(in) :: n !! Generalized integer upper argument; negative values follow GMP binomial semantics.
      integer(int64), intent(in) :: k !! Nonnegative lower argument.
      type(bigz) :: value, term
      integer(int64) :: j

      if (k < 0_int64) then
         value = bigz_from_int64(0_int64)
         return
      end if
      if (n%sign >= 0) then
         if (bigz_compare(n, bigz_from_int64(k)) < 0) then
            value = bigz_from_int64(0_int64)
            return
         end if
      end if
      value = bigz_from_int64(1_int64)
      do j = 1_int64, k
         term = bigz_sub(n, bigz_from_int64(j - 1_int64))
         value = bigz_mul(value, term)
         value = bigz_quotient(value, bigz_from_int64(j))
      end do
   end function choose_z

   pure elemental function fibnum_z(n) result(value)
      integer(int64), intent(in) :: n !! Nonnegative Fibonacci index.
      type(bigz) :: value, next

      if (n < 0_int64) then
         value = bigz_from_int64(0_int64)
         return
      end if
      call fib_pair(n, value, next)
   end function fibnum_z

   pure subroutine fibnum2_z(n, previous, value)
      integer(int64), intent(in) :: n !! Nonnegative Fibonacci index.
      type(bigz), intent(out) :: previous !! Fibonacci number F(n-1), with F(-1)=1 when n=0.
      type(bigz), intent(out) :: value !! Fibonacci number F(n).
      type(bigz) :: next

      if (n < 0_int64) then
         previous = bigz_from_int64(0_int64)
         value = bigz_from_int64(0_int64)
         return
      end if
      call fib_pair(n, value, next)
      if (n == 0_int64) then
         previous = bigz_from_int64(1_int64)
      else
         previous = bigz_sub(next, value)
      end if
   end subroutine fibnum2_z

   pure elemental function lucnum_z(n) result(value)
      integer(int64), intent(in) :: n !! Nonnegative Lucas index.
      type(bigz) :: value, fn, fn1

      if (n < 0_int64) then
         value = bigz_from_int64(0_int64)
         return
      end if
      call fib_pair(n, fn, fn1)
      value = bigz_sub(bigz_mul_small(fn1, 2_int64), fn)
   end function lucnum_z

   pure subroutine lucnum2_z(n, previous, value)
      integer(int64), intent(in) :: n !! Nonnegative Lucas index.
      type(bigz), intent(out) :: previous !! Lucas number L(n-1), with L(-1)=-1 when n=0.
      type(bigz), intent(out) :: value !! Lucas number L(n).

      if (n < 0_int64) then
         previous = bigz_from_int64(0_int64)
         value = bigz_from_int64(0_int64)
         return
      end if
      value = lucnum_z(n)
      if (n == 0_int64) then
         previous = bigz_from_int64(-1_int64)
      else
         previous = lucnum_z(n - 1_int64)
      end if
   end subroutine lucnum2_z

   pure elemental function stirling1_z(n, k) result(value)
      integer, intent(in) :: n !! Nonnegative row index for signed Stirling numbers of the first kind.
      integer, intent(in) :: k !! Column index in 0 through n.
      type(bigz) :: value
      type(bigz), allocatable :: row(:), next_row(:)
      integer :: i, j

      if (n < 0 .or. k < 0 .or. k > n) then
         value = bigz_from_int64(0_int64)
         return
      end if
      allocate(row(0:n))
      row = bigz_from_int64(0_int64)
      row(0) = bigz_from_int64(1_int64)
      do i = 1, n
         allocate(next_row(0:n))
         next_row = bigz_from_int64(0_int64)
         do j = 1, i
            next_row(j) = bigz_sub(row(j - 1), bigz_mul_small(row(j), int(i - 1, int64)))
         end do
         call move_alloc(next_row, row)
      end do
      value = row(k)
   end function stirling1_z

   pure function stirling1_all_z(n) result(values)
      integer, intent(in) :: n !! Nonnegative row index; result contains k=1 through n, matching Stirling1.all.
      type(bigz), allocatable :: values(:)
      integer :: k

      if (n <= 0) then
         allocate(values(0))
         return
      end if
      allocate(values(n))
      do k = 1, n
         values(k) = stirling1_z(n, k)
      end do
   end function stirling1_all_z

   pure elemental function stirling2_z(n, k) result(value)
      integer, intent(in) :: n !! Nonnegative row index for Stirling numbers of the second kind.
      integer, intent(in) :: k !! Column index in 0 through n.
      type(bigz) :: value
      type(bigz), allocatable :: row(:), next_row(:)
      integer :: i, j

      if (n < 0 .or. k < 0 .or. k > n) then
         value = bigz_from_int64(0_int64)
         return
      end if
      allocate(row(0:n))
      row = bigz_from_int64(0_int64)
      row(0) = bigz_from_int64(1_int64)
      do i = 1, n
         allocate(next_row(0:n))
         next_row = bigz_from_int64(0_int64)
         do j = 1, i
            next_row(j) = bigz_add(row(j - 1), bigz_mul_small(row(j), int(j, int64)))
         end do
         call move_alloc(next_row, row)
      end do
      value = row(k)
   end function stirling2_z

   pure function stirling2_all_z(n) result(values)
      integer, intent(in) :: n !! Nonnegative row index; result contains k=1 through n, matching Stirling2.all.
      type(bigz), allocatable :: values(:)
      integer :: k

      if (n <= 0) then
         allocate(values(0))
         return
      end if
      allocate(values(n))
      do k = 1, n
         values(k) = stirling2_z(n, k)
      end do
   end function stirling2_all_z

   pure elemental function eulerian_z(n, k) result(value)
      integer, intent(in) :: n !! Nonnegative Eulerian row index.
      integer, intent(in) :: k !! Ascent count; valid values are 0 through n.
      type(bigz) :: value
      type(bigz), allocatable :: row(:), next_row(:)
      integer :: i, j

      if (n < 0 .or. k < 0 .or. k > n) then
         value = bigz_from_int64(0_int64)
         return
      end if
      if (n == 0) then
         value = merge(bigz_from_int64(1_int64), bigz_from_int64(0_int64), k == 0)
         return
      end if
      allocate(row(0:n))
      row = bigz_from_int64(0_int64)
      row(0) = bigz_from_int64(1_int64)
      do i = 1, n
         allocate(next_row(0:n))
         next_row = bigz_from_int64(0_int64)
         do j = 0, i - 1
            next_row(j) = bigz_mul_small(row(j), int(j + 1, int64))
            if (j > 0) then
               next_row(j) = bigz_add(next_row(j), bigz_mul_small(row(j - 1), int(i - j, int64)))
            end if
         end do
         call move_alloc(next_row, row)
      end do
      value = row(k)
   end function eulerian_z

   pure function eulerian_all_z(n) result(values)
      integer, intent(in) :: n !! Nonnegative Eulerian row index; result contains k=0 through n-1, or one for n=0.
      type(bigz), allocatable :: values(:)
      integer :: k

      if (n < 0) then
         allocate(values(0))
      else if (n == 0) then
         allocate(values(1))
         values(1) = bigz_from_int64(1_int64)
      else
         allocate(values(n))
         do k = 0, n - 1
            values(k + 1) = eulerian_z(n, k)
         end do
      end if
   end function eulerian_all_z

   pure elemental function bernoulli_q(n) result(value)
      integer, intent(in) :: n !! Nonnegative Bernoulli index using the package convention B1 = +1/2.
      type(bigq) :: value
      type(bigq), allocatable :: a(:)
      integer :: m, j

      if (n < 0) then
         value = bigq_from_int64(0_int64)
         return
      end if
      allocate(a(0:n))
      do m = 0, n
         a(m) = bigq_make(bigz_from_int64(1_int64), bigz_from_int64(int(m + 1, int64)))
         do j = m, 1, -1
            a(j - 1) = bigq_mul(bigq_from_int64(int(j, int64)), bigq_sub(a(j - 1), a(j)))
         end do
      end do
      value = a(0)
   end function bernoulli_q

   pure elemental function dbinom_q(x, size, prob) result(value)
      integer(int64), intent(in) :: x !! Number of successes; values outside 0 through size return zero.
      type(bigz), intent(in) :: size !! Nonnegative binomial trial count.
      type(bigq), intent(in) :: prob !! Exact success probability.
      type(bigq) :: value, one_minus
      type(bigz) :: sx
      integer(int64) :: n_machine

      if (x < 0_int64 .or. size%sign < 0) then
         value = bigq_from_int64(0_int64)
         return
      end if
      sx = bigz_from_int64(x)
      if (bigz_compare(sx, size) > 0) then
         value = bigq_from_int64(0_int64)
         return
      end if
      if (.not. bigz_fits_int64(size)) then
         value = bigq_from_int64(0_int64)
         return
      end if
      n_machine = bigz_to_int64(size)
      one_minus = bigq_sub(bigq_from_int64(1_int64), prob)
      value = bigq_mul(bigq_make(choose_z(size, x), bigz_from_int64(1_int64)), &
         bigq_mul(bigq_pow(prob, x), bigq_pow(one_minus, n_machine - x)))
   end function dbinom_q

   pure elemental integer function isprime_z(n, reps) result(status)
      type(bigz), intent(in) :: n !! Integer tested for compositeness or primality.
      integer, intent(in), optional :: reps !! Requested Miller-Rabin repeat count for arbitrary-size values; defaults to 40.
      integer(int64), parameter :: deterministic_bases(7) = [2_int64, 325_int64, 9375_int64, 28178_int64, &
         450775_int64, 9780504_int64, 1795265022_int64]
      integer(int64), parameter :: probable_bases(16) = [2_int64, 3_int64, 5_int64, 7_int64, 11_int64, 13_int64, &
         17_int64, 19_int64, 23_int64, 29_int64, 31_int64, 37_int64, 41_int64, 43_int64, 47_int64, 53_int64]
      type(bigz) :: nm1, d, q, base_z
      integer(int64) :: rem
      integer :: s, i, nreps, nbases
      logical :: witness_passed

      if (n%sign <= 0 .or. bigz_compare(n, bigz_from_int64(2_int64)) < 0) then
         status = 0
         return
      end if
      if (bigz_equal(n, bigz_from_int64(2_int64)) .or. bigz_equal(n, bigz_from_int64(3_int64))) then
         status = 2
         return
      end if
      if (bigz_is_even(n)) then
         status = 0
         return
      end if
      do i = 1, 16
         if (bigz_equal(n, bigz_from_int64(probable_bases(i)))) then
            status = 2
            return
         end if
         if (bigz_is_zero(bigz_modulo(n, bigz_from_int64(probable_bases(i))))) then
            status = 0
            return
         end if
      end do

      nm1 = bigz_sub(n, bigz_from_int64(1_int64))
      d = nm1
      s = 0
      do while (bigz_is_even(d))
         call bigz_div_small(d, 2_int64, q, rem)
         d = q
         s = s + 1
      end do

      if (bigz_fits_int64(n)) then
         nbases = size(deterministic_bases)
         do i = 1, nbases
            base_z = bigz_modulo(bigz_from_int64(deterministic_bases(i)), n)
            if (bigz_is_zero(base_z)) cycle
            if (.not. miller_rabin_pass(base_z, d, s, n, nm1)) then
               status = 0
               return
            end if
         end do
         status = 2
      else
         nreps = 40
         if (present(reps)) nreps = max(1, reps)
         nbases = min(nreps, size(probable_bases))
         witness_passed = .true.
         do i = 1, nbases
            base_z = bigz_from_int64(probable_bases(i))
            if (.not. miller_rabin_pass(base_z, d, s, n, nm1)) then
               witness_passed = .false.
               exit
            end if
         end do
         status = merge(1, 0, witness_passed)
      end if
   end function isprime_z

   pure elemental function nextprime_z(n) result(p)
      type(bigz), intent(in) :: n !! Integer lower bound; the returned probable prime is strictly greater than this value.
      type(bigz) :: p

      if (bigz_compare(n, bigz_from_int64(2_int64)) < 0) then
         p = bigz_from_int64(2_int64)
         return
      end if
      p = bigz_add(n, bigz_from_int64(1_int64))
      if (bigz_is_even(p)) p = bigz_add(p, bigz_from_int64(1_int64))
      do while (isprime_z(p) == 0)
         p = bigz_add(p, bigz_from_int64(2_int64))
      end do
   end function nextprime_z

   pure subroutine factorize_z(n, factors, info)
      type(bigz), intent(in) :: n !! Positive integer to factor into probable-prime factors.
      type(bigz), allocatable, intent(out) :: factors(:) !! Factors in nondecreasing order, repeated by multiplicity.
      integer, intent(out), optional :: info !! Zero on success; one for zero input; two if Pollard rho did not split a composite.
      type(bigz), allocatable :: work(:)
      integer :: stat

      allocate(work(0))
      stat = 0
      if (bigz_is_zero(n)) then
         factors = work
         stat = 1
      else
         call factor_recursive(bigz_abs(n), work, stat)
         call sort_bigz(work)
         factors = work
      end if
      if (present(info)) info = stat
   end subroutine factorize_z

   pure subroutine urand_bigz(nb, size_bits, seed, values)
      integer, intent(in) :: nb !! Number of arbitrary-precision random integers to generate.
      integer, intent(in) :: size_bits !! Bit width; each result lies in 0 through 2^size_bits-1.
      type(bigz), intent(in) :: seed !! Deterministic seed; only its residue modulo 2^31-1 is used.
      type(bigz), allocatable, intent(out) :: values(:) !! Generated integers, reproducible for a fixed seed.
      integer(int64), parameter :: modulus = 2147483647_int64
      integer(int64), parameter :: multiplier = 16807_int64
      type(bigz) :: q
      integer(int64) :: state, rem, bit
      integer :: i, j

      allocate(values(max(0, nb)))
      call bigz_div_small(bigz_abs(seed), modulus, q, rem)
      state = rem
      if (state == 0_int64) state = 1_int64
      do i = 1, size(values)
         values(i) = bigz_from_int64(0_int64)
         do j = 1, max(0, size_bits)
            state = mod(multiplier * state, modulus)
            bit = mod(state, 2_int64)
            values(i) = bigz_mul_small(values(i), 2_int64)
            if (bit == 1_int64) values(i) = bigz_add(values(i), bigz_from_int64(1_int64))
         end do
      end do
   end subroutine urand_bigz

   pure recursive subroutine fib_pair(n, fn, fn1)
      integer(int64), intent(in) :: n !! Nonnegative Fibonacci index for fast doubling.
      type(bigz), intent(out) :: fn !! Fibonacci number F(n).
      type(bigz), intent(out) :: fn1 !! Fibonacci number F(n+1).
      type(bigz) :: a, b, c, d

      if (n == 0_int64) then
         fn = bigz_from_int64(0_int64)
         fn1 = bigz_from_int64(1_int64)
         return
      end if
      call fib_pair(n / 2_int64, a, b)
      c = bigz_mul(a, bigz_sub(bigz_mul_small(b, 2_int64), a))
      d = bigz_add(bigz_mul(a, a), bigz_mul(b, b))
      if (mod(n, 2_int64) == 0_int64) then
         fn = c
         fn1 = d
      else
         fn = d
         fn1 = bigz_add(c, d)
      end if
   end subroutine fib_pair

   pure elemental logical function miller_rabin_pass(base, d, s, n, nm1) result(passed)
      type(bigz), intent(in) :: base !! Miller-Rabin witness base.
      type(bigz), intent(in) :: d !! Odd factor of n-1 such that n-1=d*2^s.
      integer, intent(in) :: s !! Exponent of two in n-1.
      type(bigz), intent(in) :: n !! Odd candidate greater than three.
      type(bigz), intent(in) :: nm1 !! Cached n-1 value.
      type(bigz) :: x
      integer :: r

      x = bigz_powm(base, d, n)
      if (bigz_is_one(x) .or. bigz_equal(x, nm1)) then
         passed = .true.
         return
      end if
      do r = 1, s - 1
         x = bigz_modulo(bigz_mul(x, x), n)
         if (bigz_equal(x, nm1)) then
            passed = .true.
            return
         end if
         if (bigz_is_one(x)) then
            passed = .false.
            return
         end if
      end do
      passed = .false.
   end function miller_rabin_pass

   pure recursive subroutine factor_recursive(n, factors, info)
      type(bigz), intent(in) :: n !! Positive composite-or-prime value to factor recursively.
      type(bigz), allocatable, intent(inout) :: factors(:) !! Accumulated prime factors.
      integer, intent(inout) :: info !! Running status, set to two when Pollard rho cannot split a composite.
      type(bigz) :: d, q
      integer :: rho_info

      if (info /= 0 .or. bigz_is_one(n)) return
      if (isprime_z(n) > 0) then
         call append_bigz(factors, n)
         return
      end if
      call pollard_rho(n, d, rho_info)
      if (rho_info /= 0 .or. bigz_is_one(d) .or. bigz_equal(d, n)) then
         info = 2
         return
      end if
      q = bigz_quotient(n, d)
      call factor_recursive(d, factors, info)
      call factor_recursive(q, factors, info)
   end subroutine factor_recursive

   pure subroutine pollard_rho(n, factor, info)
      type(bigz), intent(in) :: n !! Odd composite integer to split.
      type(bigz), intent(out) :: factor !! Nontrivial factor on success.
      integer, intent(out) :: info !! Zero on success and one after exhausting deterministic restart attempts.
      type(bigz) :: x, y, c, d, diff
      integer :: restart, iter

      if (bigz_is_even(n)) then
         factor = bigz_from_int64(2_int64)
         info = 0
         return
      end if
      do restart = 1, 24
         x = bigz_from_int64(2_int64 + int(restart, int64))
         y = x
         c = bigz_from_int64(int(restart, int64))
         d = bigz_from_int64(1_int64)
         do iter = 1, 20000
            x = rho_step(x, c, n)
            y = rho_step(rho_step(y, c, n), c, n)
            diff = bigz_abs(bigz_sub(x, y))
            d = bigz_gcd(diff, n)
            if (.not. bigz_is_one(d)) exit
         end do
         if (.not. bigz_equal(d, n) .and. .not. bigz_is_one(d)) then
            factor = d
            info = 0
            return
         end if
      end do
      factor = n
      info = 1
   end subroutine pollard_rho

   pure elemental function rho_step(x, c, n) result(y)
      type(bigz), intent(in) :: x !! Current Pollard-rho iterate.
      type(bigz), intent(in) :: c !! Polynomial constant in x^2+c.
      type(bigz), intent(in) :: n !! Composite modulus.
      type(bigz) :: y
      y = bigz_modulo(bigz_add(bigz_mul(x, x), c), n)
   end function rho_step

   pure subroutine append_bigz(values, value)
      type(bigz), allocatable, intent(inout) :: values(:) !! Dynamically growing vector of arbitrary-precision integers.
      type(bigz), intent(in) :: value !! New value appended at the end of the vector.
      type(bigz), allocatable :: tmp(:)
      integer :: n

      n = size(values)
      allocate(tmp(n + 1))
      if (n > 0) tmp(1:n) = values
      tmp(n + 1) = value
      call move_alloc(tmp, values)
   end subroutine append_bigz

   pure subroutine sort_bigz(values)
      type(bigz), intent(inout) :: values(:) !! Factor vector sorted in place by numerical value.
      type(bigz) :: key
      integer :: i, j

      do i = 2, size(values)
         key = values(i)
         j = i - 1
         do while (j >= 1)
            if (bigz_compare(values(j), key) <= 0) exit
            values(j + 1) = values(j)
            j = j - 1
         end do
         values(j + 1) = key
      end do
   end subroutine sort_bigz

end module gmp_number_theory
