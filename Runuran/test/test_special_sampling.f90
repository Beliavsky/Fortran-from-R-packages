! SPDX-License-Identifier: GPL-2.0-or-later
program test_special_sampling
  use runuran
  implicit none
  type(rng_state)::rng
  type(continuous_distribution)::d
  real(dp)::x(8000),m,expected
  integer::fails
  fails=0
  call rng_seed(rng,987654321_i8)

  d=udgig(1.0_dp,2.0_dp,2.0_dp)
  call d%sample_n(rng,x)
  m=sum(x)/real(size(x),dp)
  expected=1.8143077587637895_dp
  if(abs(m-expected)>0.07_dp)then
    print *,'GIG sample mean',m,expected
    fails=fails+1
  end if

  d=udvg(1.3_dp,2.2_dp,0.5_dp,0.1_dp)
  call d%sample_n(rng,x)
  m=sum(x)/real(size(x),dp)
  expected=0.1_dp+2.0_dp*0.5_dp*1.3_dp/(2.2_dp**2-0.5_dp**2)
  if(abs(m-expected)>0.06_dp)then
    print *,'VG sample mean',m,expected
    fails=fails+1
  end if

  d=udghyp(1.2_dp,2.1_dp,0.4_dp,1.3_dp,-0.2_dp)
  call d%sample_n(rng,x)
  m=sum(x)/real(size(x),dp)
  ! E[W] for GIG(lambda, gamma^2, delta^2) times beta, plus mu.
  ! Fixed reference from high-precision Bessel-K evaluation.
  expected=0.2276448909327812_dp
  if(abs(m-expected)>0.07_dp)then
    print *,'GHYP sample mean',m,expected
    fails=fails+1
  end if

  if(fails==0)then
    print *,'test_special_sampling: PASS'
  else
    error stop 1
  end if
end program test_special_sampling
