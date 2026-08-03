! SPDX-License-Identifier: GPL-3.0-only
program demo_imputefin
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use imputefin
  implicit none
  real(dp)::prices(10),vol(10),ohlc(10,4)
  real(dp),allocatable::prices_imp(:),vol_imp(:),ohlc_imp(:,:)
  integer::i,status
  do i=1,10
    prices(i)=100.0_dp*exp(0.01_dp*real(i,dp))
    vol(i)=1.0e6_dp*(1.0_dp+0.02_dp*real(i,dp))
    ohlc(i,:)=[prices(i)*0.998_dp,prices(i)*1.01_dp,prices(i)*0.99_dp,prices(i)]
  end do
  prices(5)=ieee_value(0.0_dp,ieee_quiet_nan)
  vol(6)=ieee_value(0.0_dp,ieee_quiet_nan)
  ohlc(4,:)=ieee_value(0.0_dp,ieee_quiet_nan)
  call impute_rolling_ar1_gaussian(prices,prices_imp,rolling_window=20,status=status)
  call impute_vol(vol,vol_imp,rolling_window=20,status=status)
  call impute_ohlc(ohlc,ohlc_imp,rolling_window=20,status=status)
  write(*,'(a,*(f10.3,1x))')'prices: ',prices_imp
  write(*,'(a,*(f10.0,1x))')'volume: ',vol_imp
  write(*,'(a,4f12.4)')'imputed OHLC row 4: ',ohlc_imp(4,:)
end program demo_imputefin
