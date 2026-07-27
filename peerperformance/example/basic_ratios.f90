! SPDX-License-Identifier: GPL-2.0-or-later
program basic_ratios
  use peerperformance, only: dp, sharpe, modified_sharpe
  implicit none
  real(dp) :: returns(12,3), sr(3), msr(3)
  integer :: i
  do i=1,12
    returns(i,1)=0.010_dp+0.020_dp*sin(0.7_dp*real(i,dp))
    returns(i,2)=0.006_dp+0.018_dp*cos(0.5_dp*real(i,dp))
    returns(i,3)=-0.001_dp+0.022_dp*sin(0.4_dp*real(i,dp)+0.3_dp)
  end do
  call sharpe(returns,sr)
  call modified_sharpe(returns,0.95_dp,msr,na_negative=.false.)
  print '(a,3(1x,f10.6))','Sharpe:',sr
  print '(a,3(1x,f10.6))','Modified Sharpe:',msr
end program basic_ratios
