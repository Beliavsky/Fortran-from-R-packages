program basic_nnet
use nnet, only: dp, nnet_model_t, nnet_fit, nnet_predict
implicit none
type(nnet_model_t) :: fit
real(dp) :: x(8,1), y(8,1), decay(1), init(2)
real(dp), allocatable :: pred(:,:)
integer :: i
do i=1,8
   x(i,1)=real(i-4,dp)/2.0_dp
end do
y(:,1)=1.25_dp+2.5_dp*x(:,1)
decay=0.0_dp
init=0.0_dp
call nnet_fit(fit,x,y,hidden_size=0,initial_wts=init,linout=.true.,skip=.true., &
   rang=0.0_dp,decay=decay,maxit=200,reltol=1.0e-10_dp)
pred=nnet_predict(fit,x)
print '(a,f12.6)','objective: ',fit%value
print '(a,2f12.6)','weights:   ',fit%wts
print '(a,es12.4)','max error: ',maxval(abs(pred-y))
end program
