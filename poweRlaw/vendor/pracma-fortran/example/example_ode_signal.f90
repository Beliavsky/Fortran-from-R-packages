! SPDX-License-Identifier: GPL-3.0-or-later
module example_ode_functions
   use pracma_kinds, only : dp
   implicit none
contains
   subroutine oscillator(t,y,dydt)
      real(dp),intent(in)::t,y(:)
      real(dp),intent(out)::dydt(:)
      dydt(1)=y(2)
      dydt(2)=-y(1)-0.1_dp*y(2)+0.0_dp*t
   end subroutine oscillator
end module example_ode_functions

program example_ode_signal
   use pracma
   use example_ode_functions
   implicit none
   type(ode_result) :: sol
   type(peak_result) :: peaks
   real(dp), allocatable :: grid(:), values(:), smooth(:)

   sol=ode45(oscillator,[0.0_dp,20.0_dp],[1.0_dp,0.0_dp])
   grid=linspace(0.0_dp,4.0_dp,21)
   values=sin(grid)+0.05_dp*cos(13.0_dp*grid)
   smooth=savgol(values,7,3)
   peaks=findpeaks(smooth)

   print '(a,i0)','accepted ODE steps: ',sol%accepted_steps
   print '(a,2f12.6)','final oscillator state: ',sol%y(:,size(sol%t))
   print '(a,i0)','smoothed peaks: ',size(peaks%indices)
end program example_ode_signal
