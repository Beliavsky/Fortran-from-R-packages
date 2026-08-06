program demo_tvgarchkf
   use fgarch_kinds, only : dp
   use tvgarchkf
   implicit none
   type(tvgarch_spec) :: truth, initial
   type(tvgarch_simulation_result) :: simulation
   type(tvgarch_fit_result) :: fit
   integer :: n

   n = 500
   truth = make_tvgarch_spec( &
      make_tv_function([0.04_dp,0.03_dp],tv_polynomial), &
      make_tv_function([0.08_dp,0.03_dp],tv_polynomial), &
      make_tv_function([0.78_dp,-0.08_dp],tv_polynomial))
   simulation = tvgarch_simulate(n,truth,seed=20260804,corrected_constraints=.true.)
   if (simulation%status /= 0) error stop trim(simulation%message)

   initial = make_tvgarch_spec( &
      make_tv_function([0.05_dp,0.01_dp],tv_polynomial), &
      make_tv_function([0.10_dp,0.01_dp],tv_polynomial), &
      make_tv_function([0.72_dp,-0.02_dp],tv_polynomial))
   fit = tvgarch_kalman_fit(simulation%returns,initial,corrected_constraints=.true.,max_iterations=1800)
   if (fit%filter%status /= 0) error stop trim(fit%message)

   write(*,'(a,f14.5)') 'Kalman criterion: ',fit%criterion
   write(*,'(a,*(f10.5,1x))') 'Estimated coefficients: ',fit%rounded_parameters
   write(*,'(a,3f12.6)') 'Final omega/alpha/beta: ',fit%filter%omega(n),fit%filter%alpha(n),fit%filter%beta(n)
   write(*,'(a,f12.6)') 'Final conditional sigma: ',fit%filter%sigma(n)
end program demo_tvgarchkf
