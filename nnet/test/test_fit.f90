program test_fit
use nnet, only: dp,nnet_model_t,nnet_fit,nnet_predict
implicit none
type(nnet_model_t)::m
real(dp)::x(8,1),y(8,1),init(2),dec(1)
real(dp),allocatable::p(:,:)
integer::i
do i=1,8
x(i,1)=real(i-4,dp)/2.0_dp
end do
y(:,1)=1.25_dp+2.5_dp*x(:,1)
init=0.0_dp
dec=0.0_dp
call nnet_fit(m,x,y,0,initial_wts=init,linout=.true.,skip=.true.,rang=0.0_dp,decay=dec,maxit=200,reltol=1e-10_dp,hessian=.true.)
p=nnet_predict(m,x)
if(maxval(abs(p-y))>1e-6_dp) then
print *,m%wts,maxval(abs(p-y))
error stop 'linear fit'
end if
if(maxval(abs(m%hessian-transpose(m%hessian)))>1e-12_dp) error stop 'fit hessian symmetry'
print *,'test_fit passed',m%wts
end program
