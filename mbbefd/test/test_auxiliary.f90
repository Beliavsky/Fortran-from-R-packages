program test_auxiliary
  use mbbefd, only : dp, doiunif, poiunif, qoiunif, moibeta, tloibeta, &
    make_eecf, eecf_t, etl, theil_empirical, swiss_re, trans_01, invt_01
  implicit none
  real(dp) :: x(5),bg(2),z
  type(eecf_t) :: ec
  x=[0.1_dp,0.2_dp,0.5_dp,1.0_dp,1.0_dp]
  if(abs(etl(x)-0.4_dp)>1.0e-14_dp) error stop 'etl'
  ec=make_eecf(x)
  z=(sum(min(x,0.5_dp))/5.0_dp)/(sum(x)/5.0_dp)
  if(abs(ec%evaluate(0.5_dp)-z)>1.0e-13_dp) error stop 'eecf'
  if(abs(doiunif(1.0_dp,0.3_dp)-0.3_dp)>1.0e-14_dp) error stop 'oi atom'
  if(abs(poiunif(0.5_dp,0.3_dp)-0.35_dp)>1.0e-14_dp) error stop 'oi cdf'
  if(abs(qoiunif(0.8_dp,0.3_dp)-1.0_dp)>1.0e-14_dp) error stop 'oi quantile'
  if(abs(tloibeta(2.0_dp,3.0_dp,0.2_dp)-0.2_dp)>1.0e-14_dp) error stop 'oi tl'
  if(abs(moibeta(1.0_dp,2.0_dp,3.0_dp,0.2_dp)-0.52_dp)>5.0e-12_dp) error stop 'oi moment'
  bg=swiss_re(2.0_dp);if(any(bg<=0.0_dp)) error stop 'swiss'
  z=0.37_dp;if(abs(invt_01(trans_01(z))-z)>1.0e-13_dp) error stop 'transform'
  if(theil_empirical(x)<0.0_dp) error stop 'theil'
  print '(a)', 'test_auxiliary: PASS'
end program test_auxiliary
