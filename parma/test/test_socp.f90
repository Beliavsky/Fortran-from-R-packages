! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2012-2014 Alexios Galanos and Bernhard Pfaff.
! Modern Fortran translation Copyright (C) 2026 OpenAI.
! This derivative work is distributed under GPL-3.0-or-later.
program test_socp
   use parma
   implicit none
   type(socp_result) :: result
   real(dp) :: f(1),a(1,1),b(1),c(1,1),d(1),x0(1)
   integer :: ncones(1)

   f = [-1.0_dp]
   a(1,1) = 1.0_dp
   b = 0.0_dp
   c = 0.0_dp
   d = 1.0_dp
   ncones = 1
   x0 = 0.0_dp
   call socp_solve(f,a,b,c,d,ncones,result,x0=x0,max_iter=300,tol=1.0e-9_dp)
   if (result%status /= 0) then
      write(*,'(a,i0,1x,a)') 'SOCP status ',result%status,trim(result%message)
      error stop 1
   end if
   if (abs(result%x(1)-1.0_dp) > 2.0e-4_dp) then
      write(*,'(a,es24.14)') 'SOCP solution failed: ',result%x(1)
      error stop 1
   end if
   print '(a)', 'test_socp: PASS'
end program test_socp
