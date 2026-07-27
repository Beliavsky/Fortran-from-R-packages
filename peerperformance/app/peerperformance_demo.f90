! SPDX-License-Identifier: GPL-2.0-or-later
program peerperformance_demo
  use peerperformance, only: dp, peer_control, screening_result, alpha_screening
  implicit none
  integer, parameter :: n=72, p=5
  real(dp) :: returns(n,p), factors(n,2), t
  type(peer_control) :: control
  type(screening_result) :: result
  integer :: i, j

  do i=1,n
    t=real(i,dp)
    factors(i,1)=0.012_dp*sin(0.13_dp*t)
    factors(i,2)=0.009_dp*cos(0.09_dp*t)
    do j=1,p
      returns(i,j)=0.0015_dp*real(4-j,dp)+0.15_dp*real(j,dp)*factors(i,1)+ &
           0.08_dp*real(p-j,dp)*factors(i,2)+0.018_dp*sin((0.19_dp+0.03_dp*real(j,dp))*t+0.2_dp*real(j,dp))
    end do
  end do
  control%has_lambda=.true.; control%lambda=0.5_dp
  control%min_obs=24; control%screen_beta=.true.
  call alpha_screening(returns,control,result,factors=factors)
  if (result%status/=0) error stop trim(result%message)
  print '(a)', 'fund   alpha        pi+          pi0          pi-'
  do j=1,p
    print '(i4,4(1x,f11.6))',j,result%estimate(1,j),result%pipos(1,j), &
         result%pizero(1,j),result%pineg(1,j)
  end do
end program peerperformance_demo
