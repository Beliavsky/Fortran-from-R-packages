program diabetes_reference
  use mclust
  implicit none
  real(dp) :: x(145,3)
  type(mclust_fit) :: fit
  integer :: u,i,st
  integer,parameter :: gs(9)=[1,2,3,4,5,6,7,8,9]
  open(newunit=u,file='validation/diabetes_x.dat',status='old',action='read')
  do i=1,145
    read(u,*) x(i,:)
  end do
  close(u)
  call mclust_select(x,fit,g_values=gs,status=st)
  if(st/=0) error stop 'diabetes: model selection failed'
  if(fit%g/=3 .or. fit%model_name/='VVV') error stop 'diabetes: selected model mismatch'
  if(any([(count(fit%classification==i),i=1,fit%g)]/=[82,35,28])) &
    error stop 'diabetes: classification counts mismatch'
  if(maxval(abs(fit%pro-[0.5553630_dp,0.2479432_dp,0.1966939_dp]))>2.0e-5_dp) &
    error stop 'diabetes: mixing proportions mismatch'
  if(abs(fit%loglik-(-2295.118_dp))>0.05_dp) error stop 'diabetes: loglik mismatch'
  if(abs(fit%bic-(-4734.561_dp))>0.10_dp) error stop 'diabetes: BIC mismatch'
  if(abs(fit%icl-(-4749.16_dp))>0.20_dp) error stop 'diabetes: ICL mismatch'
  print *, 'diabetes_reference PASS',fit%model_name,fit%g,fit%loglik,fit%bic,fit%icl
end program diabetes_reference
