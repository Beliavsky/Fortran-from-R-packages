program test_quasi
   use, intrinsic :: iso_fortran_env, only : int64,real64
   use randtoolbox_quasi, only : halton, torus, sobol, sobol_directions
   use randtoolbox_bits, only : int2bit,bit2int,bit2unitreal
   implicit none
   real(real64),allocatable::x(:,:)
   integer(int64),allocatable::v(:,:)
   integer,allocatable::bits(:)
   real(real64),parameter::tol=2.0e-15_real64
   x=halton(3,2,start=1_int64)
   if(maxval(abs(x-reshape([0.5_real64,0.25_real64,0.75_real64, &
      1.0_real64/3.0_real64,2.0_real64/3.0_real64,1.0_real64/9.0_real64],[3,2])))>tol) error stop 'Halton mismatch'
   x=sobol(4,2,start=1_int64)
   if(maxval(abs(x-reshape([0.5_real64,0.75_real64,0.25_real64,0.375_real64, &
      0.5_real64,0.25_real64,0.75_real64,0.375_real64],[4,2])))>tol) error stop 'Sobol mismatch'
   x=torus(2,1,start=1_int64)
   if(abs(x(1,1)-(sqrt(2.0_real64)-1.0_real64))>tol) error stop 'Torus mismatch'
   v=sobol_directions(1111,20)
   if(size(v,1)/=20 .or. size(v,2)/=1111) error stop 'Sobol direction shape mismatch'
   bits=int2bit(13_int64,8)
   if(bit2int(bits)/=13_int64) error stop 'bit roundtrip mismatch'
   if(abs(bit2unitreal([1,0,0,0])-0.5_real64)>tol) error stop 'bit2unitreal mismatch'
   print '(a)', 'quasi-random tests: PASS'
end program test_quasi
