program test_general_reference
   use r_compat, only: dp
   use stabledist
   implicit none
   real(dp), parameter :: xs(5)=[-2.0_dp,-0.5_dp,0.0_dp,0.7_dp,3.0_dp]
   real(dp), parameter :: d15(5)=[ &
      0.07522409466433215_dp,0.26727525241221917_dp,0.28534475886171384_dp, &
      0.23477983286148052_dp,0.04065007198675813_dp]
   real(dp), parameter :: p15(5)=[ &
      0.07357990368713119_dp,0.32953891457654283_dp,0.46981551097802410_dp, &
      0.65589313415423840_dp,0.92675473333275080_dp]
   real(dp), parameter :: d07(5)=[ &
      0.07519071206839882_dp,0.18537431089644340_dp,0.26645679057558735_dp, &
      0.32454484263032307_dp,0.006012353186245266_dp]
   real(dp), parameter :: p07(5)=[ &
      0.33635516436113455_dp,0.51452880562210670_dp,0.62626590922235800_dp, &
      0.84993525313346100_dp,0.97249281969983460_dp]
   real(dp), parameter :: d10(5)=[ &
      0.035278035136051315_dp,0.29398603437056536_dp,0.28565927321974630_dp, &
      0.19576738871731153_dp,0.04845721359236821_dp]
   real(dp), parameter :: p10(5)=[ &
      0.06046679249142643_dp,0.27462438788927860_dp,0.42367655151451520_dp, &
      0.59302914930048170_dp,0.82828731611827460_dp]
   integer :: i

   do i=1,5
      call chk(dstable(xs(i),1.5_dp,0.4_dp),d15(i),3e-8_dp,'d 1.5/.4')
      call chk(pstable(xs(i),1.5_dp,0.4_dp),p15(i),3e-8_dp,'p 1.5/.4')
      call chk(dstable(xs(i),0.7_dp,-0.8_dp),d07(i),3e-8_dp,'d .7/-.8')
      call chk(pstable(xs(i),0.7_dp,-0.8_dp),p07(i),3e-8_dp,'p .7/-.8')
      call chk(dstable(xs(i),1.0_dp,0.6_dp),d10(i),3e-8_dp,'d 1/.6')
      call chk(pstable(xs(i),1.0_dp,0.6_dp),p10(i),3e-8_dp,'p 1/.6')
   end do

   ! Independent Nolan S1 reference points (scale/location non-unit).
   call chk(dstable(0.7_dp,1.5_dp,0.4_dp,gamma=2.0_dp,delta=0.3_dp,pm=1), &
      0.12285268974890633_dp,3e-8_dp,'S1 density alpha 1.5')
   call chk(pstable(0.7_dp,1.5_dp,0.4_dp,gamma=2.0_dp,delta=0.3_dp,pm=1), &
      0.63186262900968060_dp,3e-8_dp,'S1 cdf alpha 1.5')
   call chk(dstable(0.7_dp,1.0_dp,0.6_dp,gamma=2.3_dp,delta=-0.2_dp,pm=1), &
      0.12078890945656855_dp,3e-8_dp,'S1 density alpha 1')
   call chk(pstable(0.7_dp,1.0_dp,0.6_dp,gamma=2.3_dp,delta=-0.2_dp,pm=1), &
      0.44429334179558955_dp,3e-8_dp,'S1 cdf alpha 1')

   print *, 'test_general_reference: PASS'
contains
   subroutine chk(a,b,tol,msg)
      real(dp),intent(in)::a,b,tol
      character(len=*),intent(in)::msg
      if(abs(a-b)>tol*max(1.0_dp,abs(b)))then
         print *,trim(msg),a,b,abs(a-b)
         error stop 1
      end if
   end subroutine chk
end program test_general_reference
