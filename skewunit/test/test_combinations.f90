program test_combinations
   use skewunit
   implicit none
   real(dp), parameter :: dref(5,5) = reshape([ &
      0.5935490385505621_dp,0.6567908447883526_dp,0.4696387071150669_dp, &
      0.4500236212497260_dp,0.5391494401817583_dp, &
      0.1825767894769627_dp,0.20203008683520002_dp,0.14446174080000002_dp, &
      0.13842810386352336_dp,0.16584337172387112_dp, &
      1.3324144399699445_dp,1.47438130432_dp,1.0542572799999999_dp, &
      1.0102248210947462_dp,1.2102967956179942_dp, &
      1.5766680459479514_dp,1.744659784771503_dp,1.2475200776280218_dp, &
      1.4033372563406499_dp,1.5132771206582156_dp, &
      1.0323375942215753_dp,1.1423316972617847_dp,0.8168249993975837_dp, &
      0.91884769962544_dp,0.9908318153245723_dp ],[5,5])
   real(dp), parameter :: pref(5,5) = reshape([ &
      0.29323916616505064_dp,0.3705211977568351_dp,0.13675662543389933_dp, &
      0.10244892403740752_dp,0.20605292414226042_dp, &
      0.34205211340748765_dp,0.4372286814519039_dp,0.15114007313920003_dp, &
      0.10823031393678921_dp,0.23566707837161133_dp, &
      0.21958415601571898_dp,0.26319543949568003_dp,0.1356317584_dp, &
      0.11769188996432801_dp,0.17762181350458583_dp, &
      0.20231173972753505_dp,0.238520284213793_dp,0.13303289182465203_dp, &
      0.162021601874411_dp,0.1874113023658084_dp, &
      0.26778505135686936_dp,0.32780549768737516_dp,0.15086086914499475_dp, &
      0.19424956952057507_dp,0.24017432935682578_dp ],[5,5])
   real(dp) :: got, xv(3), pv(3)
   integer :: f1, f2

   do f1 = 1, 5
      do f2 = 1, 5
         got = dskewunit(0.37_dp,0.6_dp,1.3_dp,0.8_dp,f1,f2)
         call check_close('density combination',got,dref(f2,f1),3.0e-11_dp)
         got = pskewunit(0.37_dp,0.6_dp,1.3_dp,0.8_dp,f1,f2)
         call check_close('cdf combination',got,pref(f2,f1),4.0e-9_dp)
      end do
   end do

   xv = [0.1_dp,0.4_dp,0.9_dp]
   call pskewunit_vec(xv,pv,lambda=-0.3_dp,delta=1.2_dp, &
      family1=family_triang,family2=family_sbeta)
   do f1 = 1, 3
      call check_close('vector cdf',pv(f1),pskewunit(xv(f1),-0.3_dp,1.2_dp, &
         family1=family_triang,family2=family_sbeta),1.0e-13_dp)
   end do

   print '(a)', 'test_combinations: PASS'

contains

   subroutine check_close(name, got, expected, tol)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: got, expected, tol
      if (abs(got-expected) > tol*max(1.0_dp,abs(expected))) then
         print '(a,2(1x,es24.16))', trim(name)//' got/expected:',got,expected
         error stop 1
      end if
   end subroutine check_close

end program test_combinations
