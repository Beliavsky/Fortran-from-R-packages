program test_v03
   use gamlss_dist
   implicit none
   integer, parameter :: n=600
   real(dp) :: y(n),x1(n,1),p,s,alpha
   integer :: bd(n),i,k
   type(gamlss_fit_result_t) :: fit

   call assert_close(dSN1(0.7_dp,0.1_dp,1.2_dp,0.8_dp),0.3845854493807751_dp,2.0e-14_dp,'dSN1')
   call assert_close(pSN2(0.4_dp,-0.2_dp,0.9_dp,1.7_dp),0.4837060949973497_dp,2.0e-14_dp,'pSN2')
   call assert_close(dGT(1.1_dp,0.2_dp,1.3_dp,4.5_dp,1.7_dp),0.2369355429705677_dp,3.0e-14_dp,'dGT')
   call assert_close(pGT(1.1_dp,0.2_dp,1.3_dp,4.5_dp,1.7_dp),0.8086840502459027_dp,3.0e-14_dp,'pGT')
   call assert_close(dexGAUS(4.2_dp,3.5_dp,0.8_dp,1.1_dp),0.3501771713522542_dp,2.0e-14_dp,'dexGAUS')
   call assert_close(pexGAUS(4.2_dp,3.5_dp,0.8_dp,1.1_dp),0.4240181586600088_dp,2.0e-14_dp,'pexGAUS')
   call assert_close(dPARETO(2.4_dp,1.7_dp),1.7_dp*2.4_dp**(-2.7_dp),2.0e-15_dp,'dPARETO')
   call assert_close(pPARETO1(1.5_dp,2.2_dp),1.0_dp-2.5_dp**(-2.2_dp),2.0e-15_dp,'pPARETO1')
   call assert_close(pPARETO2(qPARETO2(0.73_dp,1.4_dp,0.6_dp),1.4_dp,0.6_dp),0.73_dp,2.0e-14_dp,'PARETO2 rt')
   call assert_close(pPARETO2o(qPARETO2o(0.73_dp,1.4_dp,2.1_dp),1.4_dp,2.1_dp),0.73_dp,2.0e-14_dp,'PARETO2o rt')
   call assert_close(dST3C(0.3_dp,-0.1_dp,1.1_dp,1.4_dp,4.5_dp), &
      dST3(0.3_dp,-0.1_dp,1.1_dp,1.4_dp,4.5_dp),2.0e-15_dp,'ST3C alias')
   call assert_close(pSST(qSST(0.67_dp,0.4_dp,1.2_dp,0.8_dp,7.0_dp),0.4_dp,1.2_dp,0.8_dp,7.0_dp), &
      0.67_dp,4.0e-12_dp,'SST rt')
   call assert_close(pSN1(qSN1(0.63_dp,0.2_dp,1.1_dp,-1.3_dp),0.2_dp,1.1_dp,-1.3_dp), &
      0.63_dp,3.0e-8_dp,'SN1 rt')

   call assert_close(dGAF(1.7_dp,2.1_dp,0.45_dp,1.6_dp), &
      dgamma_v(1.7_dp,1.0_dp/(0.45_dp*2.1_dp**(-0.2_dp))**2, &
      2.1_dp*(0.45_dp*2.1_dp**(-0.2_dp))**2),2.0e-14_dp,'GAF map')
   call assert_close(dNBF(3.0_dp,2.4_dp,0.35_dp,1.7_dp), &
      dnbinom_v(3,1.0_dp/(0.35_dp*2.4_dp**(-0.3_dp)), &
      1.0_dp/(1.0_dp+2.4_dp*0.35_dp*2.4_dp**(-0.3_dp))),2.0e-14_dp,'NBF map')
   call assert_close(dZINBF(0.0_dp,2.4_dp,0.35_dp,1.7_dp,0.18_dp), &
      0.18_dp+0.82_dp*dNBF(0.0_dp,2.4_dp,0.35_dp,1.7_dp),2.0e-14_dp,'ZINBF zero')

   call assert_close(dDBI(3.0_dp,0.35_dp,0.6_dp,8.0_dp),0.3633253187926869_dp,3.0e-14_dp,'dDBI')
   s=0.0_dp
   do k=0,8
      s=s+dDBI(real(k,dp),0.35_dp,0.6_dp,8.0_dp)
   end do
   call assert_close(s,1.0_dp,3.0e-14_dp,'DBI norm')
   alpha=1.0_dp/(sqrt(2.3_dp**2+0.7_dp**2)-2.3_dp)
   call assert_close(dPIG2(4.0_dp,2.3_dp,0.7_dp),dPIG(4.0_dp,2.3_dp,alpha),2.0e-14_dp,'PIG2 map')

   call assert_close(dZIPIG(0.0_dp,2.3_dp,0.6_dp,0.2_dp), &
      0.2_dp+0.8_dp*dPIG(0.0_dp,2.3_dp,0.6_dp),2.0e-14_dp,'ZIPIG zero')
   call assert_close(dZAPIG(0.0_dp,2.3_dp,0.6_dp,0.2_dp),0.2_dp,2.0e-15_dp,'ZAPIG zero')
   call assert_close(dZISICHEL(0.0_dp,2.3_dp,0.6_dp,-0.4_dp,0.15_dp), &
      0.15_dp+0.85_dp*dSICHEL(0.0_dp,2.3_dp,0.6_dp,-0.4_dp),2.0e-13_dp,'ZISICHEL zero')
   call assert_close(dZASICHEL(0.0_dp,2.3_dp,0.6_dp,-0.4_dp,0.15_dp),0.15_dp,2.0e-14_dp,'ZASICHEL zero')
   call assert_close(dZIBB(0.0_dp,0.4_dp,0.3_dp,0.12_dp,10), &
      0.12_dp+0.88_dp*dBB(0.0_dp,0.4_dp,0.3_dp,10),2.0e-14_dp,'ZIBB zero')
   call assert_close(dZABB(0.0_dp,0.4_dp,0.3_dp,0.12_dp,10),0.12_dp,2.0e-15_dp,'ZABB zero')
   call assert_close(dZIBNB(0.0_dp,2.0_dp,0.5_dp,1.2_dp,0.1_dp), &
      0.1_dp+0.9_dp*dBNB(0.0_dp,2.0_dp,0.5_dp,1.2_dp),2.0e-14_dp,'ZIBNB zero')
   call assert_close(dZABNB(0.0_dp,2.0_dp,0.5_dp,1.2_dp,0.1_dp),0.1_dp,2.0e-15_dp,'ZABNB zero')
   call assert_close(dZAZIPF(0.0_dp,1.5_dp,0.17_dp),0.17_dp,2.0e-15_dp,'ZAZIPF zero')

   s=0.0_dp
   do k=0,150
      s=s+dZAPIG(real(k,dp),2.3_dp,0.6_dp,0.2_dp)
   end do
   call assert_close(s,1.0_dp,3.0e-9_dp,'ZAPIG norm')
   s=0.0_dp
   do k=0,150
      s=s+dZASICHEL(real(k,dp),2.3_dp,0.6_dp,-0.4_dp,0.15_dp)
   end do
   call assert_close(s,1.0_dp,3.0e-8_dp,'ZASICHEL norm')

   x1=1.0_dp
   do i=1,n
      p=(real(i,dp)-0.5_dp)/real(n,dp)
      y(i)=qPARETO1(p,2.4_dp)
   end do
   call fit_gamlss(y,x1,GAMLSS_PARETO1,fit,max_iter=100,tol=1.0e-9_dp)
   if(.not.fit%converged)error stop 'PARETO1 fit'
   if(abs(fit%fitted_mu(1)-2.4_dp)>0.04_dp)error stop 'PARETO1 recovery'

   do i=1,n
      p=(real(i,dp)-0.5_dp)/real(n,dp)
      y(i)=real(qDBI(p,0.42_dp,0.72_dp,10.0_dp),dp)
      bd(i)=10
   end do
   call fit_dbi(y,bd,x1,x1,fit,max_iter=150,tol=1.0e-9_dp)
   if(.not.fit%converged)error stop 'DBI fit'
   if(abs(fit%fitted_mu(1)-0.42_dp)>0.025_dp)error stop 'DBI mu recovery'
   if(abs(fit%fitted_sigma(1)-0.72_dp)>0.08_dp)error stop 'DBI sigma recovery'

   do i=1,n
      p=(real(i,dp)-0.5_dp)/real(n,dp)
      y(i)=real(qZIBB(p,0.42_dp,0.28_dp,0.14_dp,10),dp)
   end do
   call fit_zibb(y,bd,x1,x1,x1,fit,max_iter=180,tol=1.0e-9_dp)
   if(.not.fit%converged)error stop 'ZIBB fit'
   if(abs(fit%fitted_mu(1)-0.42_dp)>0.03_dp)error stop 'ZIBB mu recovery'
   if(abs(fit%fitted_sigma(1)-0.28_dp)>0.07_dp)error stop 'ZIBB sigma recovery'
   if(abs(fit%fitted_nu(1)-0.14_dp)>0.04_dp)error stop 'ZIBB nu recovery'

   print '(a)', 'test_v03: PASS'
contains
   subroutine assert_close(a,b,tol,msg)
      real(dp),intent(in)::a,b,tol
      character(*),intent(in)::msg
      if(abs(a-b)>tol)then
         print '(a,2es24.15)',trim(msg)//' failed: ',a,b
         error stop 1
      end if
   end subroutine assert_close
end program test_v03
