program demo_tsmarch
  use ghyp_kinds, only : dp, i8
  use tsgarch, only : garch_spec, fit_options
  use tsmarch
  implicit none
  type(dcc_spec) :: dspec
  type(dcc_parameters) :: dpar
  type(dcc_simulation) :: sim
  type(dcc_fit) :: fit
  type(dcc_forecast) :: pred
  type(garch_spec) :: marginal
  type(fit_options) :: options
  real(dp) :: qbar(3,3)
  real(dp), allocatable :: data(:, :)
  integer :: i

  qbar = reshape([1.0_dp,0.45_dp,0.20_dp, &
                  0.45_dp,1.0_dp,0.30_dp, &
                  0.20_dp,0.30_dp,1.0_dp],[3,3])
  dspec%distribution='mvn'
  dspec%alpha_order=1
  dspec%gamma_order=1
  dspec%beta_order=1
  allocate(dpar%alpha(1),dpar%gamma(1),dpar%beta(1))
  dpar%alpha=0.04_dp
  dpar%gamma=0.02_dp
  dpar%beta=0.90_dp
  dpar%shape=8.0_dp

  sim=simulate_dcc_innovations(dspec,dpar,qbar,250,paths=1,burn=120,seed=20260804_i8)
  allocate(data(250,3))
  do i=1,250
    data(i,:)=[0.012_dp,0.016_dp,0.010_dp]*sim%innovations(i,:,1)
  end do

  marginal%model='ewma'
  marginal%distribution='norm'
  marginal%constant=.false.
  dspec%constant_correlation=.true.
  options%max_iterations=400
  options%compute_inference=.false.
  fit=estimate_dcc(data,marginal,dspec,options)
  pred=forecast_dcc(data,fit,5,paths=500,seed=17_i8)

  write(*,'(a,f12.4)') 'multivariate log likelihood: ',fit%log_likelihood
  write(*,'(a,3f10.5)') 'last fitted correlations: ', &
    fit%filtered%correlation(2,1,size(data,1)), &
    fit%filtered%correlation(3,1,size(data,1)), &
    fit%filtered%correlation(3,2,size(data,1))
  write(*,'(a,3f10.5)') 'one-step sigma forecast: ',pred%sigma(1,:)
end program demo_tsmarch
