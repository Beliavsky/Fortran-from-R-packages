! SPDX-License-Identifier: GPL-2.0-or-later
program peer_screening
  use peerperformance, only: dp, peer_control, screening_result, sharpe_screening
  implicit none
  real(dp) :: returns(60,4), t
  type(peer_control) :: control
  type(screening_result) :: result
  integer :: i
  do i=1,60
    t=real(i,dp)
    returns(i,1)=0.006_dp+0.018_dp*sin(0.31_dp*t)
    returns(i,2)=0.003_dp+0.020_dp*cos(0.27_dp*t)
    returns(i,3)=0.001_dp+0.021_dp*sin(0.21_dp*t+0.5_dp)
    returns(i,4)=-0.001_dp+0.019_dp*cos(0.25_dp*t+0.2_dp)
  end do
  control%has_lambda=.true.; control%lambda=0.5_dp; control%min_obs=20
  call sharpe_screening(returns,control,result)
  if (result%status/=0) error stop trim(result%message)
  print '(a)', 'fund  Sharpe      pi+       pi0       pi-'
  do i=1,4
    print '(i4,4(1x,f9.5))',i,result%estimate(1,i),result%pipos(1,i), &
         result%pizero(1,i),result%pineg(1,i)
  end do
end program peer_screening
