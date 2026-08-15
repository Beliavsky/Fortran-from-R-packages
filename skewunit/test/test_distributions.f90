program test_distributions
   use skewunit
   implicit none
   integer :: failures

   failures = 0
   call check_close('dasin',dasin(0.2_dp),0.7957747154594766_dp,2.0e-13_dp)
   call check_close('pasin',pasin(0.2_dp),0.2951672353008665_dp,2.0e-12_dp)
   call check_close('dtriang',dtriang(0.2_dp),0.8_dp,1.0e-14_dp)
   call check_close('ptriang',ptriang(0.2_dp),0.08_dp,1.0e-14_dp)
   call check_close('duquad',duquad(0.2_dp),1.08_dp,1.0e-14_dp)
   call check_close('puquad',puquad(0.2_dp),0.392_dp,1.0e-14_dp)
   call check_close('djsb',djsb(0.2_dp,1.3_dp),0.6389532477240347_dp,2.0e-13_dp)
   call check_close('pjsb',pjsb(0.2_dp,1.3_dp),0.03575833521737554_dp,2.0e-13_dp)
   call check_close('dsbeta',dsbeta(0.2_dp,1.3_dp),1.0242770466978477_dp,2.0e-12_dp)
   call check_close('psbeta',psbeta(0.2_dp,1.3_dp),0.1624710378665525_dp,2.0e-12_dp)

   call check_close('skew density asin/asin', &
      dskewunit(0.2_dp,0.5_dp,1.0_dp,1.0_dp,family_asin,family_asin), &
      0.6414156136996497_dp,3.0e-12_dp)
   call check_close('skew cdf asin/asin', &
      pskewunit(0.2_dp,0.5_dp,1.0_dp,1.0_dp,family_asin,family_asin), &
      0.21116644729194248_dp,2.0e-9_dp)
   call check_close('skew density triang/JSB', &
      dskewunit(0.2_dp,-0.4_dp,1.4_dp,1.0_dp,family_triang,family_jsb), &
      1.2055103701167003_dp,3.0e-12_dp)
   call check_close('skew cdf triang/JSB', &
      pskewunit(0.8_dp,-0.4_dp,1.4_dp,1.0_dp,family_triang,family_jsb), &
      0.9680118905515998_dp,2.0e-9_dp)
   call check_close('skew density sbeta/JSB', &
      dskewunit(0.35_dp,0.7_dp,1.3_dp,0.8_dp,family_sbeta,family_jsb), &
      0.8344650190510489_dp,5.0e-12_dp)
   call check_close('skew cdf sbeta/JSB', &
      pskewunit(0.35_dp,0.7_dp,1.3_dp,0.8_dp,family_sbeta,family_jsb), &
      0.15286416555030066_dp,3.0e-9_dp)

   call check_close('lambda zero reduces to baseline', &
      dskewunit(0.37_dp,0.0_dp,1.0_dp,1.0_dp,family_triang,family_asin), &
      dtriang(0.37_dp),5.0e-14_dp)
   call check_close('upper log JSB',pjsb(0.2_dp,1.3_dp,.false.,.true.), &
      log(1.0_dp-0.03575833521737554_dp),3.0e-13_dp)

   if (failures == 0) then
      print '(a)', 'test_distributions: PASS'
   else
      print '(a,i0)', 'test_distributions: FAIL ', failures
      error stop 1
   end if

contains

   subroutine check_close(name, got, expected, tol)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: got, expected, tol
      if (abs(got-expected) > tol*max(1.0_dp,abs(expected))) then
         failures = failures+1
         print '(a,2(1x,es24.16))', trim(name)//' got/expected:',got,expected
      end if
   end subroutine check_close

end program test_distributions
