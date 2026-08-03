! SPDX-License-Identifier: Artistic-2.0
program demo_ldhmm
   use ldhmm
   implicit none
   type(ldhmm_model) :: model, decoded
   real(dp) :: param(2,3), gamma_matrix(2,2), delta(2), x(10)
   real(dp), allocatable :: state_forecast(:, :), statistics(:, :)
   integer :: i, status

   param(1,:) = [0.002_dp,0.015_dp,1.0_dp]
   param(2,:) = [-0.003_dp,0.040_dp,1.4_dp]
   gamma_matrix(1,:) = [0.97_dp,0.03_dp]
   gamma_matrix(2,:) = [0.08_dp,0.92_dp]
   delta = [0.70_dp,0.30_dp]
   x = [0.010_dp,-0.020_dp,0.005_dp,0.030_dp,-0.010_dp, &
      0.008_dp,0.012_dp,-0.035_dp,0.006_dp,0.018_dp]

   model = ldhmm_create(2,param,gamma_matrix,delta,status=status)
   decoded = ldhmm_decode(model,x,status=status)
   state_forecast = ldhmm_forecast_state(model,x,5,status)
   statistics = ldhmm_ld_stats(model,annualize=.true.)

   write(*,'(a,es14.6)') 'minus log-likelihood: ', ldhmm_mllk(model,x)
   write(*,'(a)') 'local states:'
   write(*,'(*(i0,1x))') decoded%states_local
   write(*,'(a)') 'five-step state forecasts:'
   do i = 1, size(state_forecast,2)
      write(*,'(i3,2(1x,f10.6))') i, state_forecast(:,i)
   end do
   write(*,'(a)') 'annualized state statistics: R, V, kurtosis, allocation'
   do i = 1, model%m
      write(*,'(i3,4(1x,es13.5))') i, statistics(i,:)
   end do
end program demo_ldhmm
