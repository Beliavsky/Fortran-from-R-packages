! SPDX-License-Identifier: GPL-3.0-only
program test_smoothing
   use smoots
   implicit none
   real(dp),parameter::tol=2.0e-10_dp
   real(dp)::y(21)
   real(dp),allocatable::g(:),w(:,:),k(:)
   integer::i,status
   y=[1.0_dp,2.0_dp,4.0_dp,7.0_dp,11.0_dp,16.0_dp,22.0_dp,29.0_dp,37.0_dp,46.0_dp,56.0_dp, &
      67.0_dp,79.0_dp,92.0_dp,106.0_dp,121.0_dp,137.0_dp,154.0_dp,172.0_dp,191.0_dp,211.0_dp]
   call gsmooth(y,0,1,1,0.2_dp,1,g,w,status)
   call assert_close(g(1),-2.421369863013699_dp,tol,'gsmooth left')
   call assert_close(g(11),58.4_dp,tol,'gsmooth center')
   call assert_close(g(21),207.57863013698636_dp,tol,'gsmooth right')
   call assert_close(w(1,1),0.43442092154420936_dp,tol,'weight')
   call gsmooth(y,1,2,1,0.2_dp,1,g,w,status)
   call assert_close(g(1),10.5_dp,1.0e-8_dp,'derivative left')
   call assert_close(g(11),220.5_dp,1.0e-8_dp,'derivative center')
   call gsmooth(y,0,3,2,0.2_dp,1,g,w,status)
   do i=1,21
      call assert_close(g(i),y(i),1.0e-8_dp,'cubic polynomial reproduction')
   end do
   call knsmooth(y,1,0.2_dp,0,k,status)
   call assert_close(k(1),3.498245614035088_dp,tol,'kernel left')
   call assert_close(k(11),58.073619631901835_dp,tol,'kernel center')
   call assert_close(k(21),183.28771929824563_dp,tol,'kernel right')
   print '(a)','test_smoothing: PASS'
contains
   subroutine assert_close(a,b,t,msg)
      real(dp),intent(in)::a,b,t
      character(len=*),intent(in)::msg
      if(abs(a-b)>t)then
         print *,trim(msg),a,b,abs(a-b);error stop 1
      end if
   end subroutine assert_close
end program test_smoothing
