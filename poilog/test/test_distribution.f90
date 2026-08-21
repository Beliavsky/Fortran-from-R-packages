! SPDX-License-Identifier: GPL-3.0-only
! Derived from the GPL-3 R package poilog by Vidar Grotan and Steinar Engen.
program test_distribution
   use poilog, only : dp,dpoilog,dbipoilog
   implicit none
   integer :: n1,n2
   real(dp) :: s,p1,p2,p12


   if(abs(dpoilog(0,0.0_dp,1.0_dp)-0.38175646475548336_dp)>2.0e-9_dp) then
      print *, 'FAIL reference dpoilog n=0'
      error stop 1
   end if
   if(abs(dpoilog(10,0.0_dp,1.0_dp)-0.003846867437111658_dp)>2.0e-8_dp) then
      print *, 'FAIL reference dpoilog n=10'
      error stop 1
   end if
   if(abs(dbipoilog(2,3,0.0_dp,0.0_dp,1.0_dp,1.0_dp,0.5_dp)- &
          0.013350525377362848_dp)>3.0e-8_dp) then
      print *, 'FAIL reference dbipoilog'
      error stop 1
   end if

   s=0.0_dp
   do n1=0,500
      s=s+dpoilog(n1,0.0_dp,1.0_dp)
   end do
   if(abs(s-1.0_dp)>2.0e-6_dp) then
      print *, 'FAIL univariate normalization',s
      error stop 1
   end if

   do n1=0,8
      do n2=0,8
         p1=dpoilog(n1,0.2_dp,0.8_dp)
         p2=dpoilog(n2,-0.1_dp,1.1_dp)
         p12=dbipoilog(n1,n2,0.2_dp,-0.1_dp,0.8_dp,1.1_dp,0.0_dp)
         if(abs(p12-p1*p2)>2.0e-6_dp*max(1.0_dp,p1*p2)) then
            print *, 'FAIL rho=0 factorization',n1,n2,p12,p1*p2
            error stop 1
         end if
      end do
   end do
   print *, 'test_distribution: PASS'
end program test_distribution
