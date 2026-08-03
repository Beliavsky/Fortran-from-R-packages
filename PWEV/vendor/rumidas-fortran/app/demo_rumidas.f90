program demo_rumidas
  use rumidas
  implicit none
  integer,parameter::n=90,k=3
  real(dp)::ret(n),z(n),g(n),mv(k+1,n),truth(5),start(5),tau
  type(garch_midas_spec)::spec
  type(rumidas_fit_result)::fit
  type(rumidas_fit_control)::control
  integer::i,status

  truth=[0.08_dp,0.84_dp,log(1.0e-4_dp),0.0_dp,2.5_dp]
  tau=exp(truth(3))
  do i=1,n
    z(i)=sin(1.31_dp*real(i,dp))+0.35_dp*cos(0.57_dp*real(i,dp))
    mv(:,i)=[0.1_dp,0.2_dp,0.3_dp,0.4_dp]
  end do
  z=z/sqrt(sum(z*z)/real(n,dp))
  g(1)=1.0_dp;ret(1)=sqrt(tau)*z(1)
  do i=2,n
    g(i)=1.0_dp-truth(1)-truth(2)+truth(1)*ret(i-1)**2/tau+truth(2)*g(i-1)
    ret(i)=sqrt(tau*g(i))*z(i)
  end do

  spec=garch_midas_spec(RUMIDAS_GM,RUMIDAS_NORMAL,RUMIDAS_BETA_LAG,k,0,.false.)
  start=[0.05_dp,0.75_dp,log(sum(ret*ret)/real(n,dp)),0.0_dp,2.0_dp]
  control%random_starts=1;control%max_iterations=250
  call fit_garch_midas(spec,ret,mv,fit,status,start=start,control=control)
  print '(a)', 'rumidas modern Fortran demonstration'
  print '(a,l1)', 'converged: ',fit%converged
  print '(a,f12.4)', 'log likelihood: ',fit%loglik
  print '(a,5f12.6)', 'coefficients: ',fit%coefficients
  print '(a,es14.5)', 'last fitted variance: ',fit%conditional(n)
end program demo_rumidas
