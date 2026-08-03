! SPDX-License-Identifier: Artistic-2.0
program forecast_example
   use ldhmm
   implicit none
   type(ldhmm_model) :: model
   real(dp) :: param(2,3), gamma_matrix(2,2), delta(2), x(8)
   real(dp), allocatable :: forecast(:, :)
   integer :: i

   param(1,:) = [0.002936740_dp,0.01977561_dp,1.141693_dp]
   param(2,:) = [-0.001707031_dp,0.03718047_dp,1.324177_dp]
   gamma_matrix(1,:) = [0.98083875_dp,0.01916125_dp]
   gamma_matrix(2,:) = [0.04931245_dp,0.95068755_dp]
   delta = [0.7201662_dp,0.2798338_dp]
   x = [0.01_dp,-0.02_dp,0.005_dp,0.03_dp,-0.01_dp,0.0_dp,0.04_dp,-0.03_dp]

   model = ldhmm_create(2,param,gamma_matrix,delta)
   forecast = ldhmm_forecast_state(model,x,10)
   do i = 1, size(forecast,2)
      write(*,'(i3,2(1x,f11.7))') i, forecast(:,i)
   end do
end program forecast_example
