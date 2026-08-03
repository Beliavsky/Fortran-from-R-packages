! SPDX-License-Identifier: Artistic-2.0
program fit_example
   use ldhmm
   implicit none
   type(ldhmm_model) :: model, fitted
   type(ldhmm_fit_control) :: control
   real(dp) :: param(1,2), gamma_matrix(1,1), delta(1), x(21)
   integer :: status

   param(1,:) = [0.05_dp,0.08_dp]
   gamma_matrix = 1.0_dp
   delta = 1.0_dp
   x = [-0.030_dp,-0.020_dp,-0.015_dp,-0.010_dp,-0.008_dp,-0.005_dp, &
      -0.003_dp,0.0_dp,0.002_dp,0.004_dp,0.006_dp,0.008_dp,0.010_dp, &
       0.012_dp,0.014_dp,0.016_dp,0.018_dp,0.020_dp,0.023_dp,0.027_dp,0.032_dp]
   model = ldhmm_create(1,param,gamma_matrix,delta)

   control%optimizer = 'bfgs'
   control%max_iterations = 100
   fitted = ldhmm_fit(model,x,control,status)
   write(*,'(a,i0)') 'status: ', status
   write(*,'(a,i0)') 'iterations: ', fitted%iterations
   write(*,'(a,es14.6)') 'mllk: ', fitted%mllk
   write(*,'(a,2(1x,es14.6))') 'mu, sigma: ', fitted%param(1,:)
end program fit_example
