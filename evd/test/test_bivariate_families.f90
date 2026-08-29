program test_bivariate_families
   use evd
   implicit none
   real(dp) :: m1(3),m2(3),q1,q2,a,h,p,d
   integer :: fails
   fails=0
   m1=[0.0_dp,1.0_dp,0.1_dp]
   m2=[0.2_dp,1.2_dp,-0.05_dp]
   q1=0.8_dp
   q2=0.9_dp

   p=pbvlog(q1,q2,0.65_dp,m1,m2)
   d=dbvlog(q1,q2,0.65_dp,m1,m2)
   a=abvlog(0.4_dp,0.65_dp)
   h=hbvlog(0.4_dp,0.65_dp)
   call valid('log',p,d,a,h,fails)
   p=pbvalog(q1,q2,0.65_dp,[0.7_dp,0.8_dp],m1,m2)
   d=dbvalog(q1,q2,0.65_dp,[0.7_dp,0.8_dp],m1,m2)
   a=abvalog(0.4_dp,0.65_dp,0.7_dp,0.8_dp)
   h=hbvalog(0.4_dp,0.65_dp,0.7_dp,0.8_dp)
   call valid('alog',p,d,a,h,fails)
   p=pbvhr(q1,q2,1.1_dp,m1,m2)
   d=dbvhr(q1,q2,1.1_dp,m1,m2)
   a=abvhr(0.4_dp,1.1_dp)
   h=hbvhr(0.4_dp,1.1_dp)
   call valid('hr',p,d,a,h,fails)
   p=pbvneglog(q1,q2,1.2_dp,m1,m2)
   d=dbvneglog(q1,q2,1.2_dp,m1,m2)
   a=abvneglog(0.4_dp,1.2_dp)
   h=hbvneglog(0.4_dp,1.2_dp)
   call valid('neglog',p,d,a,h,fails)
   p=pbvaneglog(q1,q2,1.2_dp,[0.7_dp,0.8_dp],m1,m2)
   d=dbvaneglog(q1,q2,1.2_dp,[0.7_dp,0.8_dp],m1,m2)
   a=abvaneglog(0.4_dp,1.2_dp,0.7_dp,0.8_dp)
   h=hbvaneglog(0.4_dp,1.2_dp,0.7_dp,0.8_dp)
   call valid('aneglog',p,d,a,h,fails)
   p=pbvbilog(q1,q2,0.4_dp,0.6_dp,m1,m2)
   d=dbvbilog(q1,q2,0.4_dp,0.6_dp,m1,m2)
   a=abvbilog(0.4_dp,0.4_dp,0.6_dp)
   h=hbvbilog(0.4_dp,0.4_dp,0.6_dp)
   call valid('bilog',p,d,a,h,fails)
   p=pbvnegbilog(q1,q2,0.7_dp,1.1_dp,m1,m2)
   d=dbvnegbilog(q1,q2,0.7_dp,1.1_dp,m1,m2)
   a=abvnegbilog(0.4_dp,0.7_dp,1.1_dp)
   h=hbvnegbilog(0.4_dp,0.7_dp,1.1_dp)
   call valid('negbilog',p,d,a,h,fails)
   p=pbvct(q1,q2,1.5_dp,2.0_dp,m1,m2)
   d=dbvct(q1,q2,1.5_dp,2.0_dp,m1,m2)
   a=abvct(0.4_dp,1.5_dp,2.0_dp)
   h=hbvct(0.4_dp,1.5_dp,2.0_dp)
   call valid('ct',p,d,a,h,fails)
   p=pbvamix(q1,q2,0.2_dp,0.1_dp,m1,m2)
   d=dbvamix(q1,q2,0.2_dp,0.1_dp,m1,m2)
   a=abvamix(0.4_dp,0.2_dp,0.1_dp)
   h=hbvamix(0.4_dp,0.2_dp,0.1_dp)
   call valid('amix',p,d,a,h,fails)

   if(fails>0) then
   write(*,'(a,i0)') 'test_bivariate_families: FAIL ',fails
   error stop 1
   end if
   write(*,'(a)') 'test_bivariate_families: PASS'
contains
   subroutine valid(name,p,d,a,h,fails)
      character(len=*),intent(in)::name
      real(dp),intent(in)::p,d,a,h
      integer,intent(inout)::fails
      if(.not.(p>0.0_dp.and.p<1.0_dp.and.d>0.0_dp.and.a>=0.5_dp.and.a<=1.0_dp.and.h>=0.0_dp)) then
         write(*,'(a,4(1x,es12.4))') 'FAIL '//name,p,d,a,h
         fails=fails+1
      end if
   end subroutine
end program
