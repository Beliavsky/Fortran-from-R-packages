program basic_arithmetic
   use, intrinsic :: iso_fortran_env, only : int64
   use gmp_api
   implicit none

   type(bigz) :: a, b
   type(bigq) :: q

   a = bigz_from_string('123456789012345678901234567890')
   b = factorial_z(50_int64)
   q = bigq_make(a, b)

   print '(a)', 'a              = '//bigz_to_string(a)
   print '(a)', '50!            = '//bigz_to_string(b)
   print '(a)', 'a / 50! exact  = '//bigq_to_string(q)
   print '(a)', 'choose(100,50) = '//bigz_to_string(choose_z(bigz_from_int64(100_int64), 50_int64))
end program basic_arithmetic
