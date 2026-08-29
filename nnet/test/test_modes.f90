program test_modes
use nnet, only: dp,nnet_model_t,build_network,nnet_objective,nnet_gradient,nnet_hessian_exact,nnet_predict_raw
implicit none
type(nnet_model_t)::m
real(dp)::x(2,2),yb(2,1),ys(2,3),yc(2,3),cw(2),wb(3),ws(9),f
real(dp),allocatable::p(:,:),g(:),h(:,:),gp(:),gm(:),wp(:),wm(:)
real(dp)::eps,err
integer::j
x=reshape([1.0_dp,0.5_dp,-1.0_dp,2.0_dp],[2,2])
cw=1.0_dp
! binary logistic / entropy
call build_network(m,2,0,1,entropy=.true.,skip=.true.)
wb=[0.1_dp,0.4_dp,-0.2_dp]
yb(:,1)=[1.0_dp,0.0_dp]
p=nnet_predict_raw(m,x,wb)
if(any(abs(p(:,1)-1.0_dp/(1.0_dp+exp(-[0.7_dp,-0.1_dp])))>1e-14_dp)) error stop 'entropy prediction'
f=nnet_objective(m,x,yb,cw,wb)
g=nnet_gradient(m,x,yb,cw,wb)
h=nnet_hessian_exact(m,x,yb,cw,wb)
if(.not.(f>0.0_dp).or.any(g/=g)) error stop 'entropy objective'
eps=1e-6_dp
err=0.0_dp
allocate(wp(3),wm(3))
do j=1,3
   wp=wb
   wm=wb
   wp(j)=wp(j)+eps
   wm(j)=wm(j)-eps
   gp=nnet_gradient(m,x,yb,cw,wp)
   gm=nnet_gradient(m,x,yb,cw,wm)
   err=max(err,maxval(abs(h(:,j)-(gp-gm)/(2*eps))))
end do
if(err>2e-6_dp) then
print *,'entropy hess err',err
error stop 'entropy hessian'
end if
deallocate(wp,wm)
! ordinary softmax
call build_network(m,2,0,3,softmax=.true.,skip=.true.)
ws=[0.0_dp,0.2_dp,-0.1_dp, 0.1_dp,-0.3_dp,0.4_dp, -0.2_dp,0.5_dp,0.2_dp]
ys=0.0_dp
ys(1,2)=1.0_dp
ys(2,3)=1.0_dp
p=nnet_predict_raw(m,x,ws)
if(maxval(abs(sum(p,dim=2)-1.0_dp))>2e-15_dp) error stop 'softmax sum'
h=nnet_hessian_exact(m,x,ys,cw,ws)
allocate(wp(9),wm(9))
err=0.0_dp
do j=1,9
   wp=ws
   wm=ws
   wp(j)=wp(j)+eps
   wm(j)=wm(j)-eps
   gp=nnet_gradient(m,x,ys,cw,wp)
   gm=nnet_gradient(m,x,ys,cw,wm)
   err=max(err,maxval(abs(h(:,j)-(gp-gm)/(2*eps))))
end do
if(err>3e-6_dp) then
print *,'softmax hess err',err
error stop 'softmax hessian'
end if
deallocate(wp,wm)
! censored softmax: first row allowed classes 1/2, second row 2/3
yc=0.0_dp
yc(1,1:2)=1.0_dp
yc(2,2:3)=1.0_dp
call build_network(m,2,0,3,softmax=.true.,censored=.true.,skip=.true.)
f=nnet_objective(m,x,yc,cw,ws)
if(abs(f+log(sum(p(1,1:2)))+log(sum(p(2,2:3))))>2e-14_dp) error stop 'censored likelihood'
h=nnet_hessian_exact(m,x,yc,cw,ws)
allocate(wp(9),wm(9))
err=0.0_dp
do j=1,9
   wp=ws
   wm=ws
   wp(j)=wp(j)+eps
   wm(j)=wm(j)-eps
   gp=nnet_gradient(m,x,yc,cw,wp)
   gm=nnet_gradient(m,x,yc,cw,wm)
   err=max(err,maxval(abs(h(:,j)-(gp-gm)/(2*eps))))
end do
if(err>4e-6_dp) then
print *,'censored hess err',err
error stop 'censored hessian'
end if
print *,'test_modes passed'
end program
