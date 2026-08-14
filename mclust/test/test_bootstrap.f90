program test_bootstrap
  use mclust
  implicit none
  integer,parameter::n=80
  real(dp)::x(n,2)
  type(mclust_fit)::fit
  type(parameter_bootstrap_result)::pb
  type(bootstrap_lrt_result)::lrt
  integer::i,st
  do i=1,n/2
    x(i,:)=[-2.5_dp+0.25_dp*sin(0.43_dp*i),-2.0_dp+0.2_dp*cos(0.37_dp*i)]
  end do
  do i=n/2+1,n
    x(i,:)=[2.5_dp+0.25_dp*sin(0.47_dp*i),2.0_dp+0.2_dp*cos(0.41_dp*i)]
  end do
  call fit_model(x,2,'EII',fit); if(fit%status<0)error stop 'fit'
  call mclust_parameter_bootstrap(fit,x,4,pb,parametric=.true.,status=st)
  if(st/=0 .or. pb%nboot/=4 .or. size(pb%pro,1)/=4)error stop 'parameter bootstrap'
  call bootstrap_lrt(x,'EII',1,4,lrt,status=st)
  if(st/=0 .or. lrt%p_value<0.0_dp .or. lrt%p_value>1.0_dp)error stop 'lrt bootstrap'
  print *, 'test_bootstrap PASS ',lrt%p_value
end program test_bootstrap
