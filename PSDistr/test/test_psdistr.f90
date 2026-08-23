program test_psdistr
   use psdistr, only : dp, dtppn,ptppn,qtppn, dpc,ppc,qpc, ddsn,pdsn,qdsn, den,pen,qen, dspc,pspc,qspc, deck,peck,qeck
   implicit none
   integer :: fails
   fails=0
   call chk(dtppn(2.0_dp,1.0_dp,1.0_dp,1.0_dp,2.0_dp),0.4839414490_dp,5e-8_dp,'dtppn')
   call chk(ptppn(2.0_dp,1.0_dp,1.0_dp,1.0_dp,2.0_dp),0.8413447461_dp,5e-8_dp,'ptppn')
   call chk(qtppn(0.5_dp,1.0_dp,1.0_dp,1.0_dp,2.0_dp),1.0_dp,5e-8_dp,'qtppn')
   call chk(dpc(0.0_dp,1.0_dp,2.0_dp,2.0_dp),0.1933340584_dp,5e-8_dp,'dpc')
   call chk(ppc(0.0_dp,1.0_dp,2.0_dp,2.0_dp),0.4012936743_dp,5e-8_dp,'ppc')
   call chk(qpc(0.5_dp,1.0_dp,2.0_dp,2.0_dp),1.0_dp,5e-8_dp,'qpc')
   call chk(ddsn(-0.5_dp,2.0_dp,2.0_dp,2.0_dp,0.0_dp),1.054_dp,6e-4_dp,'ddsn')
   call chk(pdsn(-0.5_dp,2.0_dp,2.0_dp,2.0_dp,0.0_dp),0.7734_dp,6e-5_dp,'pdsn')
   call chk(qdsn(0.5_dp,2.0_dp,2.0_dp,2.0_dp,0.0_dp),-0.6823_dp,6e-5_dp,'qdsn')
   call chk(den(1.0_dp,1.0_dp,2.0_dp,2.0_dp,2.0_dp,1.0_dp),0.2666_dp,6e-5_dp,'den')
   call chk(pen(1.0_dp,1.0_dp,2.0_dp,2.0_dp,2.0_dp,1.0_dp),0.7279_dp,6e-5_dp,'pen')
   call chk(qen(0.5_dp,1.0_dp,2.0_dp,2.0_dp,2.0_dp,1.0_dp),0.2910_dp,6e-5_dp,'qen')
   call chk(dspc(0.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,0.0_dp),0.2420_dp,6e-5_dp,'dspc')
   call chk(pspc(0.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,0.0_dp),0.8413_dp,6e-5_dp,'pspc')
   call chk(qspc(0.5_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,0.0_dp),-0.6823_dp,6e-5_dp,'qspc')
   call chk(deck(1.0_dp,2.0_dp,3.0_dp),0.2307_dp,6e-5_dp,'deck')
   call chk(peck(1.0_dp,2.0_dp,3.0_dp),0.9294_dp,6e-5_dp,'peck')
   call chk(qeck(0.1_dp,1.0_dp,1.0_dp),-0.6084_dp,2e-4_dp,'qeck')
   call inversion_tests()
   if(fails/=0)then;write(*,'(a,i0)')'test_psdistr: FAIL ',fails;error stop 1;end if
   write(*,'(a)')'test_psdistr: PASS'
contains
   subroutine chk(x,ref,tol,name)
      real(dp),intent(in)::x,ref,tol;character(*),intent(in)::name
      if(abs(x-ref)>tol)then;fails=fails+1;write(*,'(a,2es18.8)')trim(name)//' mismatch ',x,ref;end if
   end subroutine
   subroutine inversion_tests()
      real(dp),parameter::pv(5)=[0.01_dp,0.1_dp,0.5_dp,0.9_dp,0.99_dp]
      integer::i
      do i=1,size(pv)
         call chk(ptppn(qtppn(pv(i),0.3_dp,1.2_dp,0.7_dp,2.3_dp),0.3_dp,1.2_dp,0.7_dp,2.3_dp),pv(i),2e-10_dp,'tppn inv')
         call chk(ppc(qpc(pv(i),-0.2_dp,1.4_dp,1.7_dp),-0.2_dp,1.4_dp,1.7_dp),pv(i),2e-10_dp,'pc inv')
         call chk(pdsn(qdsn(pv(i),1.3_dp,0.8_dp,-0.4_dp,0.5_dp),1.3_dp,0.8_dp,-0.4_dp,0.5_dp),pv(i),2e-10_dp,'dsn inv')
         call chk(pen(qen(pv(i),1.0_dp,2.0_dp,2.0_dp,1.5_dp,0.4_dp),1.0_dp,2.0_dp,2.0_dp,1.5_dp,0.4_dp),pv(i),2e-10_dp,'en inv')
         call chk(pspc(qspc(pv(i),0.9_dp,1.1_dp,0.2_dp,1.5_dp,-0.4_dp),0.9_dp,1.1_dp,0.2_dp,1.5_dp,-0.4_dp),pv(i),2e-10_dp,'spc inv')
         call chk(peck(qeck(pv(i),2.0_dp,0.7_dp),2.0_dp,0.7_dp),pv(i),2e-10_dp,'eck inv')
      end do
   end subroutine
end program test_psdistr
