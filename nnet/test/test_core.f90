program test_core
use nnet, only: dp,nnet_model_t,build_network,nnet_objective,nnet_gradient,nnet_hessian_exact,nnet_predict_raw
implicit none
type(nnet_model_t)::m
real(dp)::x(3,2),y(3,1),cw(3),w(9),f,eps,maxg,maxh
real(dp),allocatable::g(:),gp(:),gm(:),h(:,:),p(:,:),wp(:),wm(:)
integer::j
x=reshape([0.2_dp,-0.1_dp,0.7_dp, 0.4_dp,0.8_dp,-0.3_dp],[3,2])
y(:,1)=[0.3_dp,0.8_dp,0.1_dp]
cw=[1.0_dp,2.0_dp,0.5_dp]
w=[0.1_dp,0.2_dp,-0.3_dp,-0.2_dp,0.4_dp,0.1_dp,0.05_dp,0.7_dp,-0.5_dp]
call build_network(m,2,2,1,linout=.false.,entropy=.false.,skip=.false.)
m%decay=0.01_dp
f=nnet_objective(m,x,y,cw,w)
g=nnet_gradient(m,x,y,cw,w)
h=nnet_hessian_exact(m,x,y,cw,w)
p=nnet_predict_raw(m,x,w)
if(any(p<=0.0_dp).or.any(p>=1.0_dp)) error stop 'core prediction range'
eps=1.0e-6_dp
maxg=0.0_dp
maxh=0.0_dp
allocate(wp(9),wm(9))
do j=1,9
 wp=w
 wm=w
 wp(j)=wp(j)+eps
 wm(j)=wm(j)-eps
 maxg=max(maxg,abs(g(j)-(nnet_objective(m,x,y,cw,wp)-nnet_objective(m,x,y,cw,wm))/(2*eps)))
 gp=nnet_gradient(m,x,y,cw,wp)
 gm=nnet_gradient(m,x,y,cw,wm)
 maxh=max(maxh,maxval(abs(h(:,j)-(gp-gm)/(2*eps))))
end do
if(maxg>2e-7_dp) then
print *, 'gradient err',maxg
error stop 'gradient'
end if
if(maxh>3e-5_dp) then
print *, 'hessian err',maxh
error stop 'hessian'
end if
if(abs(f-0.0_dp)<tiny(1.0_dp)) error stop 'objective unexpectedly zero'
print '(a,es12.4,a,es12.4)','test_core passed; max grad err=',maxg,' hess err=',maxh
end program
