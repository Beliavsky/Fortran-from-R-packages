! SPDX-License-Identifier: GPL-3.0-only
program test_moments
   use mcrp_module
   implicit none

   real(dp) :: r(6, 3), w(3)
   real(dp), allocatable :: a2(:, :), a3(:, :), a4(:, :)
   real(dp), allocatable :: v(:)
   integer :: failures

   failures = 0
   r = transpose(reshape([ &
      0.012_dp, -0.004_dp,  0.009_dp, &
     -0.006_dp,  0.011_dp,  0.003_dp, &
      0.018_dp,  0.002_dp, -0.005_dp, &
     -0.009_dp, -0.007_dp,  0.014_dp, &
      0.004_dp,  0.015_dp, -0.002_dp, &
      0.021_dp, -0.003_dp,  0.006_dp], [3, 6]))
   w = [0.40_dp, 0.35_dp, 0.25_dp]

   a2 = m2(r)
   call check(size(a2, 1) == 3 .and. size(a2, 2) == 3, 'M2 shape')
   call close(a2(1, 1), 1.5506666666666668e-4_dp, 1.0e-16_dp, 'M2(1,1)')
   call close(a2(1, 2), -2.2266666666666668e-5_dp, 1.0e-16_dp, 'M2(1,2)')
   call close(a2(3, 3), 4.9366666666666669e-5_dp, 1.0e-16_dp, 'M2(3,3)')

   a3 = m3(r)
   call check(size(a3, 1) == 3 .and. size(a3, 2) == 9, 'M3 shape')
   call close(a3(1, 1), -2.2407407407407452e-7_dp, 1.0e-18_dp, 'M3(1,1)')
   call close(a3(2, 5), 2.4407407407407414e-7_dp, 1.0e-18_dp, 'M3(2,5)')
   call close(a3(3, 9), 1.0592592592592556e-8_dp, 1.0e-18_dp, 'M3(3,9)')

   a4 = m4(r)
   call check(size(a4, 1) == 3 .and. size(a4, 2) == 27, 'M4 shape')
   call close(a4(1, 1), 2.4258407407407414e-8_dp, 1.0e-19_dp, 'M4(1,1)')
   call close(a4(1, 14), -7.4285185185185230e-10_dp, 1.0e-20_dp, 'M4(1,14)')
   call close(a4(2, 6), -1.3169444444444443e-9_dp, 1.0e-20_dp, 'M4(2,6)')
   call close(a4(3, 27), 3.0692476851851853e-9_dp, 1.0e-20_dp, 'M4(3,27)')

   call close(pm2(r, w), 1.6655750000000007e-5_dp, 1.0e-16_dp, 'pm2')
   v = dm2(r, w)
   call close(v(1), 9.1000000000000016e-5_dp, 1.0e-16_dp, 'dm2(1)')
   call close(sum(cm2(r, w, .true.)), 1.0_dp, 1.0e-12_dp, 'cm2 percentage sum')
   call close(sum(cm2(r, w, .false.)), pm2(r, w), 1.0e-16_dp, &
      'cm2 absolute sum')
   call close(port_risk(r, w), pm2(r, w), 1.0e-18_dp, 'PortRisk alias')
   v = port_risk_deriv(r, w)
   call close(v(2), 1.5239999999999996e-5_dp, 1.0e-18_dp, 'PortRiskDeriv alias')
   call close(sum(port_risk_contrib(r, w, .true.)), 1.0_dp, 1.0e-12_dp, &
      'PortRiskContrib percentage sum')

   call close(pm3(r, w), -4.4785312500000002e-8_dp, 1.0e-18_dp, 'pm3')
   v = dm3(r, w)
   call close(v(2), -2.4371666666666664e-7_dp, 1.0e-18_dp, 'dm3(2)')
   call close(sum(cm3(r, w, .true.)), 1.0_dp, 1.0e-12_dp, 'cm3 percentage sum')
   call close(sum(cm3(r, w, .false.)), pm3(r, w), 1.0e-18_dp, &
      'cm3 absolute sum')
   call close(port_skew(r, w), -0.65885419981693472_dp, 1.0e-13_dp, 'PortSkew')
   v = port_skew_deriv(r, w)
   call close(v(1), -0.33939553989581689_dp, 1.0e-12_dp, 'PortSkewDeriv(1)')
   call close(sum(port_skew_contrib(r, w, .true.)), 1.0_dp, 1.0e-11_dp, &
      'PortSkewContrib percentage sum')
   call close(sum(port_skew_contrib(r, w, .false.)), port_skew(r, w), &
      1.0e-12_dp, 'PortSkewContrib absolute sum')

   call close(pm4(r, w), 4.8629542018229172e-10_dp, 1.0e-20_dp, 'pm4')
   v = dm4(r, w)
   call close(v(1), 4.6447822916666678e-9_dp, 1.0e-19_dp, 'dm4(1)')
   call close(sum(cm4(r, w, .true.)), 1.0_dp, 1.0e-12_dp, 'cm4 percentage sum')
   call close(sum(cm4(r, w, .false.)), pm4(r, w), 1.0e-19_dp, &
      'cm4 absolute sum')
   call close(port_kurt(r, w), 1.7529591370625068_dp, 1.0e-13_dp, 'PortKurt')
   v = port_kurt_deriv(r, w)
   call close(v(3), -2.3581219458255740_dp, 1.0e-11_dp, 'PortKurtDeriv(3)')
   call close(sum(port_kurt_contrib(r, w, .true.)), 1.0_dp, 1.0e-11_dp, &
      'PortKurtContrib percentage sum')
   call close(sum(port_kurt_contrib(r, w, .false.)), port_kurt(r, w), &
      1.0e-11_dp, 'PortKurtContrib absolute sum')

   if (failures /= 0) then
      write(*, '(a, i0)') 'test_moments failures: ', failures
      error stop 1
   end if
   write(*, '(a)') 'test_moments: all tests passed'

contains

   subroutine close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: label
      if (abs(actual - expected) > tolerance) then
         failures = failures + 1
         write(*, '(a, 2es24.15)') trim(label)//': ', actual, expected
      end if
   end subroutine close

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         failures = failures + 1
         write(*, '(a)') trim(label)//': failed'
      end if
   end subroutine check

end program test_moments
