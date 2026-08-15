program test_v02
   use gamlss_dist
   implicit none
   integer, parameter :: n=400
   real(dp) :: y(n),x(n,1),p,err,s
   real(dp) :: start2(2)
   integer :: i,k
   type(gamlss_fit_result_t) :: fit

   call assert_close(dGIG(1.3_dp,2.0_dp,0.7_dp,0.4_dp),0.4318142010370882_dp,2.0e-10_dp,'dGIG')
   call assert_close(pGIG(1.3_dp,2.0_dp,0.7_dp,0.4_dp),0.3576180638843589_dp,2.0e-8_dp,'pGIG')
   call assert_close(dSHASHo(0.7_dp,0.2_dp,1.3_dp,-0.4_dp,1.2_dp), &
      0.3009104106443373_dp,2.0e-12_dp,'dSHASHo')
   call assert_close(pSHASHo(0.7_dp,0.2_dp,1.3_dp,-0.4_dp,1.2_dp), &
      0.8307922812691645_dp,2.0e-12_dp,'pSHASHo')
   call assert_close(dSHASH(0.7_dp,0.2_dp,1.3_dp,0.8_dp,1.2_dp), &
      0.3253728432508289_dp,2.0e-12_dp,'dSHASH')
   call assert_close(pSHASH(0.7_dp,0.2_dp,1.3_dp,0.8_dp,1.2_dp), &
      0.6607936039254347_dp,2.0e-12_dp,'pSHASH')
   call assert_close(dSIMPLEX(0.4_dp,0.55_dp,0.7_dp),1.016861650374601_dp,2.0e-12_dp,'dSIMPLEX')
   call assert_close(pSIMPLEX(0.4_dp,0.55_dp,0.7_dp),0.03723123337508907_dp,3.0e-8_dp,'pSIMPLEX')
   call assert_close(dSEP(0.3_dp,-0.1_dp,1.1_dp,0.6_dp,1.7_dp), &
      0.4083872521785789_dp,2.0e-12_dp,'dSEP')
   call assert_close(pSEP(0.3_dp,-0.1_dp,1.1_dp,0.6_dp,1.7_dp), &
      0.4696533975770092_dp,4.0e-8_dp,'pSEP')
   call assert_close(dSEP1(0.3_dp,-0.1_dp,1.1_dp,0.6_dp,1.7_dp), &
      0.3946506172192555_dp,2.0e-12_dp,'dSEP1')
   call assert_close(pSEP1(0.3_dp,-0.1_dp,1.1_dp,0.6_dp,1.7_dp), &
      0.4805087667992761_dp,4.0e-8_dp,'pSEP1')
   call assert_close(dST1(0.3_dp,-0.1_dp,1.1_dp,0.6_dp,4.5_dp), &
      0.3686652156813058_dp,2.0e-12_dp,'dST1')
   call assert_close(pST1(0.3_dp,-0.1_dp,1.1_dp,0.6_dp,4.5_dp), &
      0.4652277976678861_dp,5.0e-8_dp,'pST1')
   call assert_close(dST2(0.3_dp,-0.1_dp,1.1_dp,0.6_dp,4.5_dp), &
      0.3737818492347098_dp,2.0e-12_dp,'dST2')
   call assert_close(pST2(0.3_dp,-0.1_dp,1.1_dp,0.6_dp,4.5_dp), &
      0.4736714565527808_dp,5.0e-8_dp,'pST2')
   call assert_close(dST3(0.3_dp,-0.1_dp,1.1_dp,1.4_dp,4.5_dp), &
      0.3116534850376537_dp,2.0e-12_dp,'dST3')
   call assert_close(pST3(0.3_dp,-0.1_dp,1.1_dp,1.4_dp,4.5_dp), &
      0.4659513549320484_dp,2.0e-12_dp,'pST3')
   call assert_close(dST4(0.3_dp,-0.1_dp,1.1_dp,1.4_dp,4.5_dp), &
      0.2990948831474671_dp,2.0e-12_dp,'dST4')
   call assert_close(pST4(0.3_dp,-0.1_dp,1.1_dp,1.4_dp,4.5_dp), &
      0.6543435663871513_dp,2.0e-12_dp,'pST4')
   call assert_close(dSEP3(0.3_dp,-0.1_dp,1.1_dp,1.4_dp,1.7_dp), &
      0.3047413893671900_dp,2.0e-12_dp,'dSEP3')
   call assert_close(pSEP3(0.3_dp,-0.1_dp,1.1_dp,1.4_dp,1.7_dp), &
      0.4636906230124702_dp,2.0e-12_dp,'pSEP3')
   call assert_close(dSEP4(0.3_dp,-0.1_dp,1.1_dp,1.4_dp,1.7_dp), &
      0.4213681013346156_dp,2.0e-12_dp,'dSEP4')
   call assert_close(pSEP4(0.3_dp,-0.1_dp,1.1_dp,1.4_dp,1.7_dp), &
      0.6942558259048140_dp,2.0e-12_dp,'pSEP4')
   call assert_close(dST5(0.3_dp,-0.1_dp,1.1_dp,0.6_dp,1.7_dp), &
      0.2629346487267686_dp,2.0e-12_dp,'dST5')
   call assert_close(pST5(0.3_dp,-0.1_dp,1.1_dp,0.6_dp,1.7_dp), &
      0.4237991765476417_dp,2.0e-12_dp,'pST5')
   call assert_close(dNET(0.7_dp,0.2_dp,1.3_dp,1.5_dp,2.2_dp), &
      0.2675174013723437_dp,2.0e-12_dp,'dNET')
   call assert_close(pNET(0.7_dp,0.2_dp,1.3_dp,1.5_dp,2.2_dp), &
      0.6405536258188026_dp,2.0e-12_dp,'pNET')

   call assert_close(dGPO(4.0_dp,2.3_dp,0.4_dp),0.06695282038356563_dp,2.0e-12_dp,'dGPO')
   call assert_close(dDPO(4.0_dp,2.3_dp,0.7_dp),0.1142632106929987_dp,2.0e-12_dp,'dDPO')
   call assert_close(dDEL(4.0_dp,2.3_dp,0.6_dp,0.35_dp),0.09790325815865808_dp,2.0e-12_dp,'dDEL')
   call assert_close(dSI(4.0_dp,2.3_dp,0.6_dp,-0.4_dp),0.08529596521792872_dp,2.0e-10_dp,'dSI')
   call assert_close(dSICHEL(4.0_dp,2.3_dp,0.6_dp,-0.4_dp),0.08190988895640468_dp,2.0e-10_dp,'dSICHEL')
   call assert_close(dYULE(4.0_dp,2.2_dp),0.02625118575558054_dp,2.0e-12_dp,'dYULE')
   call assert_close(dWARING(4.0_dp,2.2_dp,0.7_dp),0.04392544753986090_dp,2.0e-12_dp,'dWARING')
   call assert_close(dZIPF(4.0_dp,1.5_dp),0.02329504050902428_dp,2.0e-12_dp,'dZIPF')
   if (dGPO(2.5_dp,2.3_dp,0.4_dp)/=0.0_dp) error stop 'dGPO noninteger support'
   if (dDPO(2.5_dp,2.3_dp,0.7_dp)/=0.0_dp) error stop 'dDPO noninteger support'

   s=0.0_dp
   do k=0,120
      s=s+dGPO(real(k,dp),2.3_dp,0.4_dp)
   end do
   call assert_close(s,1.0_dp,2.0e-9_dp,'GPO normalization')
   s=0.0_dp
   do k=0,120
      s=s+dDEL(real(k,dp),2.3_dp,0.6_dp,0.35_dp)
   end do
   call assert_close(s,1.0_dp,2.0e-10_dp,'DEL normalization')
   s=0.0_dp
   do k=0,120
      s=s+dSICHEL(real(k,dp),2.3_dp,0.6_dp,-0.4_dp)
   end do
   call assert_close(s,1.0_dp,2.0e-8_dp,'SICHEL normalization')

   p=0.73_dp
   call assert_close(pSHASHo(qSHASHo(p,0.2_dp,1.3_dp,-0.4_dp,1.2_dp), &
      0.2_dp,1.3_dp,-0.4_dp,1.2_dp),p,2.0e-10_dp,'SHASHo round trip')
   call assert_close(pST5(qST5(p,-0.1_dp,1.1_dp,0.6_dp,1.7_dp), &
      -0.1_dp,1.1_dp,0.6_dp,1.7_dp),p,2.0e-10_dp,'ST5 round trip')
   call assert_close(pST3(qST3(p,-0.1_dp,1.1_dp,1.4_dp,4.5_dp), &
      -0.1_dp,1.1_dp,1.4_dp,4.5_dp),p,3.0e-10_dp,'ST3 round trip')
   call assert_close(pST4(qST4(p,-0.1_dp,1.1_dp,1.4_dp,4.5_dp), &
      -0.1_dp,1.1_dp,1.4_dp,4.5_dp),p,3.0e-10_dp,'ST4 round trip')
   call assert_close(pSEP3(qSEP3(p,-0.1_dp,1.1_dp,1.4_dp,1.7_dp), &
      -0.1_dp,1.1_dp,1.4_dp,1.7_dp),p,3.0e-10_dp,'SEP3 round trip')
   call assert_close(pSEP4(qSEP4(p,-0.1_dp,1.1_dp,1.4_dp,1.7_dp), &
      -0.1_dp,1.1_dp,1.4_dp,1.7_dp),p,3.0e-10_dp,'SEP4 round trip')

   x=1.0_dp
   do i=1,n
      p=(real(i,dp)-0.5_dp)/real(n,dp)
      y(i)=real(qGPO(p,2.3_dp,0.4_dp),dp)
   end do
   start2=[log(2.0_dp),log(0.5_dp)]
   call fit_gamlss(y,x,GAMLSS_GPO,fit,x_sigma=x,start=start2,max_iter=120,tol=1.0e-8_dp)
   if (.not.fit%converged) error stop 'GPO fit did not converge'
   err=abs(fit%fitted_mu(1)-2.3_dp)
   if (err>0.12_dp) error stop 'GPO fitted mu'
   err=abs(fit%fitted_sigma(1)-0.4_dp)
   if (err>0.12_dp) error stop 'GPO fitted sigma'

   print '(a)', 'test_v02: PASS'
contains
   subroutine assert_close(a,b,tol,msg)
      real(dp), intent(in) :: a,b,tol
      character(*), intent(in) :: msg
      if (abs(a-b)>tol) then
         print '(a,2es24.15)', trim(msg)//' failed: ',a,b
         error stop 1
      end if
   end subroutine assert_close
end program test_v02
