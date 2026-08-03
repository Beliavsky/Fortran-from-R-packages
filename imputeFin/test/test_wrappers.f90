! SPDX-License-Identifier: GPL-3.0-only
program test_wrappers
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_nan
  use imputefin
  implicit none
  integer,parameter::n=120
  real(dp)::ohlc(n,4),vol(n),nanv,c,spread
  real(dp),allocatable::oi(:,:),vi(:),rolling(:)
  logical,allocatable::mask(:)
  integer::i,st
  nanv=ieee_value(0.0_dp,ieee_quiet_nan)
  do i=1,n
    c=100.0_dp*exp(0.002_dp*real(i,dp)+0.01_dp*sin(0.15_dp*real(i,dp)))
    spread=0.01_dp+0.002_dp*cos(0.1_dp*real(i,dp))
    ohlc(i,4)=c;ohlc(i,1)=c*exp(0.002_dp*sin(real(i,dp)))
    ohlc(i,2)=c*exp(spread);ohlc(i,3)=c*exp(-spread);vol(i)=1.0e6_dp*exp(0.1_dp*sin(0.05_dp*real(i,dp)))
  end do
  ohlc(30:34,:)=nanv;ohlc(70,4)=nanv;vol(40:45)=nanv
  call impute_ohlc(ohlc,oi,rolling_window=200,seed=42_8,status=st)
  call check(st==impute_ok,'OHLC status')
  call check(.not.any(ieee_is_nan(oi)),'OHLC complete')
  call check(all(oi(:,2)>=max(oi(:,1),oi(:,4))),'OHLC high constraint')
  call check(all(oi(:,3)<=min(oi(:,1),oi(:,4))),'OHLC low constraint')
  call impute_vol(vol,vi,rolling_window=200,seed=43_8,status=st)
  call check(st==impute_ok.and..not.any(ieee_is_nan(vi)),'volume imputation')
  mask=is_inner_na([nanv,nanv,1.0_dp,nanv,2.0_dp,nanv])
  call check(.not.mask(1).and..not.mask(2).and.mask(4).and..not.mask(6),'inner NA mask')
  call impute_rolling_ar1_gaussian(vol,rolling,rolling_window=200,seed=44_8,status=st)
  call check(st==impute_ok,'rolling API')
  print '(a)', 'test_wrappers: PASS'
contains
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(*),intent(in)::msg
    if(.not.ok)then;write(*,'(a)')'FAIL: '//trim(msg);error stop 1;end if
  end subroutine check
end program test_wrappers
