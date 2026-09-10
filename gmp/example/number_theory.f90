program number_theory
   use, intrinsic :: iso_fortran_env, only : int64
   use gmp_api
   implicit none

   type(bigz), allocatable :: factors(:)
   type(bigz) :: n
   integer :: i, info

   n = bigz_from_string('1234567891011')
   call factorize_z(n, factors, info)
   if (info /= 0) error stop 'factorization failed'

   print '(a)', 'n = '//bigz_to_string(n)
   print '(a,i0)', 'isprime(n) status = ', isprime_z(n)
   print '(a)', 'prime factors:'
   do i = 1, size(factors)
      print '(2x,a)', bigz_to_string(factors(i))
   end do
   print '(a)', 'F(200) = '//bigz_to_string(fibnum_z(200_int64))
   print '(a)', 'B(20)  = '//bigq_to_string(bernoulli_q(20))
end program number_theory
