program basic_mclust
  use mclust
  implicit none
  real(dp)::x(100,2)
  type(mclust_fit)::fit
  integer::i,status
  do i=1,50
    x(i,:)=[-2.0_dp+0.3_dp*sin(0.6_dp*i),-1.5_dp+0.2_dp*cos(0.4_dp*i)]
  end do
  do i=51,100
    x(i,:)=[2.5_dp+0.3_dp*sin(0.5_dp*i),2.0_dp+0.2_dp*cos(0.7_dp*i)]
  end do
  call mclust_select(x,fit,g_values=[1,2,3],model_names=[character(len=3)::'EII','EEE','VVV'],status=status)
  if(status/=0) error stop 'mclust_select failed'
  print '(a,i0,a,a,a,f12.4)', 'selected G=',fit%g,' model=',trim(fit%model_name),' BIC=',fit%bic
end program basic_mclust
