program test_gmp
   use, intrinsic :: iso_fortran_env, only : int64
   use gmp_api
   implicit none

   type(bigz) :: a, b, c, qz, rz, g, s, t, seed
   type(bigz), allocatable :: factors(:), random1(:), random2(:), zv(:), zd(:), zu(:)
   type(bigq) :: qa, qb, qc
   type(bigz), allocatable :: am(:,:), bm(:,:)
   type(bigq), allocatable :: solution(:,:), inverse(:,:)
   logical, allocatable :: dup(:)
   integer :: info

   if (bigz_to_string(bigz_from_int64(int(z'8000000000000000', int64))) /= &
      '-9223372036854775808') error stop 'minimum int64 conversion'
   if (bigz_to_string(bigz_from_string('FF', 16), 10) /= '255') error stop 'radix conversion'
   if (bigq_to_string(bigq_from_string('1.25')) /= '5/4') error stop 'decimal rational parsing'

   a = bigz_from_string('123456789012345678901234567890')
   b = bigz_from_string('98765432109876543210')
   c = bigz_mul(a, b)
   if (bigz_to_string(c) /= '12193263113702179522496570642237463801111263526900') error stop 'bigz multiplication'

   qz = bigz_quotient(c, a)
   rz = bigz_modulo(c, a)
   if (.not. bigz_equal(qz, b)) error stop 'bigz exact quotient'
   if (.not. bigz_is_zero(rz)) error stop 'bigz exact remainder'

   qz = bigz_quotient(bigz_from_int64(-7_int64), bigz_from_int64(3_int64))
   rz = bigz_modulo(bigz_from_int64(-7_int64), bigz_from_int64(3_int64))
   if (bigz_to_string(qz) /= '-3') error stop 'floor division quotient'
   if (bigz_to_string(rz) /= '2') error stop 'floor division remainder'

   g = bigz_gcd(bigz_from_int64(84_int64), bigz_from_int64(30_int64))
   if (bigz_to_string(g) /= '6') error stop 'gcd'
   call bigz_gcdex(bigz_from_int64(84_int64), bigz_from_int64(30_int64), g, s, t)
   if (.not. bigz_equal(bigz_add(bigz_mul(bigz_from_int64(84_int64), s), &
      bigz_mul(bigz_from_int64(30_int64), t)), g)) error stop 'gcdex identity'

   if (bigz_to_string(bigz_pow(bigz_from_int64(2_int64), 100_int64)) /= &
      '1267650600228229401496703205376') error stop 'integer power'
   if (bigz_to_string(bigz_powm(bigz_from_int64(2_int64), bigz_from_int64(100_int64), &
      bigz_from_int64(101_int64))) /= '1') error stop 'modular power'
   if (bigz_to_string(bigz_powm(bigz_from_int64(3_int64), bigz_from_int64(-1_int64), &
      bigz_from_int64(11_int64))) /= '4') error stop 'negative modular power'

   qa = bigq_from_string('1/3')
   qb = bigq_from_string('1/6')
   qc = bigq_add(qa, qb)
   if (bigq_to_string(qc) /= '1/2') error stop 'rational addition'
   if (bigq_to_string(bigq_from_string('2/4')) /= '1/2') error stop 'rational normalization'
   if (bigz_to_string(bigq_floor(bigq_from_string('-5/2'))) /= '-3') error stop 'rational floor'
   if (bigz_to_string(bigq_trunc(bigq_from_string('-5/2'))) /= '-2') error stop 'rational truncation'
   if (bigz_to_string(bigq_round0(bigq_from_string('5/2'))) /= '2') error stop 'ties-to-even down'
   if (bigz_to_string(bigq_round0(bigq_from_string('7/2'))) /= '4') error stop 'ties-to-even up'

   if (bigz_to_string(factorial_z(20_int64)) /= '2432902008176640000') error stop 'factorial'
   if (bigz_to_string(choose_z(bigz_from_int64(100_int64), 50_int64)) /= &
      '100891344545564193334812497256') error stop 'binomial coefficient'
   if (bigz_to_string(choose_z(bigz_from_int64(-5_int64), 3_int64)) /= '-35') error stop 'negative choose'
   if (bigz_to_string(fibnum_z(100_int64)) /= '354224848179261915075') error stop 'Fibonacci'
   if (bigz_to_string(lucnum_z(10_int64)) /= '123') error stop 'Lucas'
   if (bigz_to_string(stirling1_z(5, 2)) /= '-50') error stop 'Stirling first kind'
   if (bigz_to_string(stirling2_z(5, 2)) /= '15') error stop 'Stirling second kind'
   if (bigz_to_string(eulerian_z(5, 2)) /= '66') error stop 'Eulerian number'
   if (bigq_to_string(bernoulli_q(10)) /= '5/66') error stop 'Bernoulli number'
   if (bigq_to_string(dbinom_q(2_int64, bigz_from_int64(5_int64), bigq_from_string('1/3'))) /= &
      '80/243') error stop 'exact dbinom'

   if (isprime_z(bigz_from_int64(71_int64)) /= 2) error stop 'definite prime'
   if (isprime_z(bigz_from_int64(210_int64)) /= 0) error stop 'composite'
   if (bigz_to_string(nextprime_z(bigz_from_int64(14_int64))) /= '17') error stop 'next prime'
   call factorize_z(bigz_from_int64(34455342_int64), factors, info)
   if (info /= 0) error stop 'factorization status'
   if (size(factors) /= 4) error stop 'factorization size'
   if (bigz_to_string(factors(1)) /= '2') error stop 'factor 1'
   if (bigz_to_string(factors(2)) /= '3') error stop 'factor 2'
   if (bigz_to_string(factors(3)) /= '101') error stop 'factor 3'
   if (bigz_to_string(factors(4)) /= '56857') error stop 'factor 4'

   seed = bigz_from_int64(12345_int64)
   call urand_bigz(4, 40, seed, random1)
   call urand_bigz(4, 40, seed, random2)
   if (size(random1) /= size(random2)) error stop 'random size'
   if (.not. all([(bigz_equal(random1(info), random2(info)), info = 1, size(random1))])) &
      error stop 'random reproducibility'

   allocate(am(2, 2), bm(2, 1))
   am(1, 1) = bigz_from_int64(2_int64)
   am(1, 2) = bigz_from_int64(1_int64)
   am(2, 1) = bigz_from_int64(5_int64)
   am(2, 2) = bigz_from_int64(3_int64)
   bm(1, 1) = bigz_from_int64(1_int64)
   bm(2, 1) = bigz_from_int64(2_int64)
   call bigz_solve(am, bm, solution, info)
   if (info /= 0) error stop 'matrix solve status'
   if (bigq_to_string(solution(1, 1)) /= '1') error stop 'matrix solve x1'
   if (bigq_to_string(solution(2, 1)) /= '-1') error stop 'matrix solve x2'
   call bigz_inverse(am, inverse, info)
   if (info /= 0) error stop 'matrix inverse status'
   if (bigq_to_string(inverse(1, 1)) /= '3') error stop 'inverse 11'
   if (bigq_to_string(inverse(1, 2)) /= '-1') error stop 'inverse 12'
   if (bigq_to_string(inverse(2, 1)) /= '-5') error stop 'inverse 21'
   if (bigq_to_string(inverse(2, 2)) /= '2') error stop 'inverse 22'

   allocate(zv(5))
   zv(1) = bigz_from_int64(1_int64)
   zv(2) = bigz_from_int64(3_int64)
   zv(3) = bigz_from_int64(3_int64)
   zv(4) = bigz_from_int64(8_int64)
   zv(5) = bigz_from_int64(13_int64)
   if (bigz_to_string(bigz_sum(zv)) /= '28') error stop 'vector sum'
   dup = bigz_duplicated(zv)
   if (.not. dup(3)) error stop 'duplicated'
   zu = bigz_unique(zv)
   if (size(zu) /= 4) error stop 'unique size'
   zd = bigz_diff(zv, 1, 1)
   if (size(zd) /= 4) error stop 'diff size'
   if (bigz_to_string(zd(1)) /= '2') error stop 'diff value 1'
   if (bigz_to_string(zd(2)) /= '0') error stop 'diff value 2'
   if (bigz_to_string(zd(3)) /= '5') error stop 'diff value 3'
   if (bigz_to_string(zd(4)) /= '5') error stop 'diff value 4'

   print '(a)', 'All gmp deterministic tests passed.'
end program test_gmp
