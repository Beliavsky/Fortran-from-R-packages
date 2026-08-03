! SPDX-License-Identifier: Artistic-2.0
program test_fit
   use ldhmm
   implicit none
   type(ldhmm_model) :: model, fitted
   type(ldhmm_fit_control) :: control
   real(dp) :: param(1,2), gamma_matrix(1,1), delta(1), x(21), initial_value
   integer :: status

   param(1,:) = [0.05_dp,0.08_dp]
   gamma_matrix = 1.0_dp
   delta = 1.0_dp
   x = [-0.030_dp,-0.020_dp,-0.015_dp,-0.010_dp,-0.008_dp,-0.005_dp, &
      -0.003_dp,0.0_dp,0.002_dp,0.004_dp,0.006_dp,0.008_dp,0.010_dp, &
       0.012_dp,0.014_dp,0.016_dp,0.018_dp,0.020_dp,0.023_dp,0.027_dp,0.032_dp]
   model = ldhmm_create(1,param,gamma_matrix,delta,status=status)
   initial_value = ldhmm_mllk(model,x)

   control%optimizer = 'bfgs'
   control%max_iterations = 100
   control%tolerance = 1.0e-8_dp
   fitted = ldhmm_fit(model,x,control,status)
   call assert_true(status == LDHMM_SUCCESS,'BFGS status')
   call assert_true(fitted%mllk < initial_value-1.0_dp,'BFGS improves likelihood')
   call assert_true(abs(fitted%param(1,1)) < 0.02_dp,'BFGS location')
   call assert_true(fitted%param(1,2) > 0.0_dp,'BFGS scale')

   control%optimizer = 'nelder-mead'
   control%max_iterations = 300
   fitted = ldhmm_fit(model,x,control,status)
   call assert_true(status == LDHMM_SUCCESS,'Nelder-Mead status')
   call assert_true(fitted%mllk < initial_value-1.0_dp,'Nelder-Mead improves likelihood')
   print '(a)', 'test_fit: PASS'

contains

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         write(*,'(a)') 'FAIL: '//message
         error stop 1
      end if
   end subroutine assert_true

end program test_fit
