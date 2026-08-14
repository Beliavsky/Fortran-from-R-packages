program test_models
  use mclust
  implicit none
  integer,parameter :: n=120,d=2
  real(dp)::x(n,d)
  character(len=3),parameter :: mods(14)=[character(len=3) :: &
    'EII','VII','EEI','VEI','EVI','VVI','EEE', &
    'EVE','VEE','VVE','EEV','VEV','EVV','VVV']
  type(mclust_fit)::fit,best
  real(dp),allocatable::z0(:,:)
  integer::i,j,st
  integer,parameter::gs(3)=[1,2,3]
  character(len=3),parameter::selmods(4)=[character(len=3)::'EII','VII','EEE','VVV']
  do i=1,n/2
    x(i,1)=-4.0_dp+0.35_dp*sin(real(i,dp)*0.71_dp)
    x(i,2)=-2.5_dp+0.30_dp*cos(real(i,dp)*0.43_dp)
  end do
  do i=n/2+1,n
    x(i,1)=4.0_dp+0.40_dp*sin(real(i,dp)*0.67_dp)
    x(i,2)=3.0_dp+0.25_dp*cos(real(i,dp)*0.39_dp)
  end do
  call hc_responsibilities(x,2,z0,'VVV',st)
  if(st/=0 .or. size(z0,1)/=n .or. size(z0,2)/=2) error stop 'hc initialization failed'
  if(maxval(abs(sum(z0,dim=2)-1.0_dp))>1e-12_dp) error stop 'hc z invalid'
  do j=1,size(mods)
    call fit_model(x,2,mods(j),fit,z_init=z0)
    if(fit%status<0) then
      print *, 'failed model ',mods(j),fit%status
      error stop 'model fit failure'
    end if
    if(abs(sum(fit%pro)-1.0_dp)>1e-10_dp) error stop 'bad proportions'
    if(.not.(fit%loglik>-huge(1.0_dp))) error stop 'bad loglik'
    if(any(fit%sigma(1,1,:)<=0.0_dp) .or. any(fit%sigma(2,2,:)<=0.0_dp)) error stop 'bad covariance'
  end do
  call mclust_select(x,best,g_values=gs,model_names=selmods,status=st)
  if(st/=0) error stop 'selection failed'
  if(best%g/=2) then
    print *, 'selected g=',best%g,' model=',best%model_name,' bic=',best%bic
    error stop 'expected two clusters'
  end if
  if(n_var_params('VVV',2,2)/=6) error stop 'parameter count VVV'
  if(n_mclust_params('VVV',2,2)/=11) error stop 'total parameter count VVV'
  print *, 'test_models PASS ', best%model_name, best%bic
end program test_models
