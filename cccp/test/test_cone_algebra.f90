! SPDX-License-Identifier: GPL-3.0-or-later
program test_cone_algebra
   use cccp
   implicit none
   real(dp) :: s3(3), z3(3), p3(3), e3(3), jp3(3)
   real(dp) :: ss(4), zz(4), pp(4), ee(4), jpp(4)
   integer :: info, fails
   fails=0

   s3=[2.0_dp,0.3_dp,-0.2_dp]
   z3=[1.5_dp,0.1_dp,0.2_dp]
   call cone_identity(cone_socc,3,e3)
   call check(maxval(abs(e3-[1.0_dp,0.0_dp,0.0_dp]))<1e-14_dp,'SOC identity',fails)
   call cone_jordan_product(cone_socc,3,s3,z3,jp3)
   call cone_jordan_inverse(cone_socc,3,jp3,z3,p3,info)
   call check(info==0.and.maxval(abs(p3-s3))<1e-11_dp,'SOC inverse',fails)
   call check(cone_max_step(cone_socc,3,s3)<0.0_dp,'SOC interior step',fails)

   ss=reshape(reshape([2.0_dp,0.4_dp,0.4_dp,1.0_dp],[2,2]),[4])
   zz=reshape(reshape([1.5_dp,0.0_dp,0.0_dp,0.8_dp],[2,2]),[4])
   call cone_identity(cone_psdc,2,ee)
   call check(maxval(abs(ee-[1.0_dp,0.0_dp,0.0_dp,1.0_dp]))<1e-14_dp,'PSD identity',fails)
   call cone_jordan_product(cone_psdc,2,ss,zz,jpp)
   call cone_jordan_inverse(cone_psdc,2,jpp,zz,pp,info)
   call check(info==0.and.maxval(abs(pp-ss))<1e-11_dp,'PSD inverse',fails)
   call check(cone_inner_product(cone_psdc,ss,ss,2)>0.0_dp,'PSD inner product',fails)

   if(fails>0)error stop 1
   print '(a)','test_cone_algebra: PASS'
contains
   subroutine check(ok,name,fails)
      logical,intent(in)::ok
      character(len=*),intent(in)::name
      integer,intent(inout)::fails
      if(.not.ok)then;print '(a,a)','FAIL: ',name;fails=fails+1;end if
   end subroutine check
end program test_cone_algebra
