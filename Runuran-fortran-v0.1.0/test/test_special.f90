! SPDX-License-Identifier: GPL-2.0-or-later
program test_special
  use runuran
  implicit none
  type(continuous_distribution)::d
  real(dp)::f0,f1
  real(dp)::v
  d=udig(1.0_dp,2.0_dp)
  v=d%pdf(1.0_dp)
  if(abs(v-sqrt(1.0_dp/pi))>2e-13_dp)then;print *,'IG pdf',v;error stop 1;end if
  d=udhyperbolic(2.0_dp,0.0_dp,1.0_dp,0.0_dp)
  f0=d%pdf(0.0_dp)
  f1=d%pdf(1.0_dp)
  if(.not.(f0>0.0_dp .and. f1<f0))then
    print *,'hyperbolic failed'
    error stop 1
  end if
  d=udgig(1.0_dp,2.0_dp,2.0_dp)
  if(abs(d%pdf(1.0_dp)-0.4838037750126474_dp)>3e-10_dp)then
    print *,'GIG failed',d%pdf(1.0_dp);error stop 1
  end if
  d=udghyp(1.2_dp,2.1_dp,0.4_dp,1.3_dp,-0.2_dp)
  if(abs(d%pdf(0.6_dp)-0.3532239542303128_dp)>2e-8_dp)then
    print *,'GHYP failed',d%pdf(0.6_dp);error stop 1
  end if
  d=udvg(1.3_dp,2.2_dp,0.5_dp,0.1_dp)
  if(abs(d%pdf(0.9_dp)-0.3089434959699380_dp)>2e-8_dp)then
    print *,'VG failed',d%pdf(0.9_dp);error stop 1
  end if
  d=udmeixner(2.0_dp,0.2_dp,1.3_dp,0.0_dp)
  if(abs(d%pdf(0.7_dp)-0.2511528890595145_dp)>2e-10_dp)then
    print *,'Meixner failed',d%pdf(0.7_dp);error stop 1
  end if
  d=udplanck(1.0_dp)
  if(abs(d%pdf(1.0_dp)-6.0_dp/(pi*pi*(exp(1.0_dp)-1.0_dp)))>2e-10_dp)then
    print *,'Planck failed',d%pdf(1.0_dp);error stop 1
  end if
  print *,'test_special: PASS'
end program
