! SPDX-License-Identifier: GPL-2.0-or-later
module gmm_linear
use r_compat, only: dp, optim_result_t, optim_bfgs, chol, chol2inv, pchisq
use gmm_linalg, only: invert_matrix, least_squares
use gmm_covariance, only: moment_covariance, hac_covariance, bw_andrews
implicit none
private
public :: linear_gmm_result_t, linear_gmm_fit, tsls_fit, linear_moments, linear_gradient
public :: LINEAR_TWO_STEP, LINEAR_ITERATIVE, LINEAR_CUE
integer,parameter :: LINEAR_TWO_STEP=1,LINEAR_ITERATIVE=2,LINEAR_CUE=3

type :: linear_gmm_result_t
   real(dp),allocatable :: coefficients(:),initial(:),fitted(:),residuals(:),moments(:,:),gradient(:,:)
   real(dp),allocatable :: omega(:,:),weight(:,:),vcov(:,:)
   real(dp)::objective=0.0_dp,j_stat=0.0_dp,j_pvalue=1.0_dp,bandwidth=0.0_dp
   integer::n=0,k=0,q=0,df=0,convergence=0,iterations=0,outer_iterations=0
end type

real(dp),allocatable,save :: cue_y(:),cue_x(:,:),cue_z(:,:)
integer,save :: cue_cov=1
real(dp),save :: cue_bw=0.0_dp
character(len=32),save :: cue_kernel='Quadratic Spectral'
logical,save :: cue_centered=.true.

contains

pure function linear_moments(theta,y,x,z) result(gt)
real(dp),intent(in)::theta(:),y(:),x(:,:),z(:,:)
real(dp)::gt(size(y),size(z,2)),e(size(y))
e=y-matmul(x,theta)
gt=z*spread(e,2,size(z,2))
end function linear_moments

pure function linear_gradient(x,z) result(g)
real(dp),intent(in)::x(:,:),z(:,:)
real(dp)::g(size(z,2),size(x,2))
g=-matmul(transpose(z),x)/real(size(x,1),dp)
end function linear_gradient

subroutine linear_gmm_fit(y,x,z,res,method,covariance,kernel,bandwidth,centered,identity_weight,tol,maxouter)
real(dp),intent(in)::y(:),x(:,:),z(:,:)
type(linear_gmm_result_t),intent(out)::res
integer,intent(in),optional::method,covariance,maxouter
character(len=*),intent(in),optional::kernel
real(dp),intent(in),optional::bandwidth,tol
logical,intent(in),optional::centered,identity_weight
integer::meth,covcode,mo,o,info
real(dp)::bw,eps,delta
logical::ctr,ident
character(len=32)::kern
real(dp),allocatable::theta(:),newtheta(:),w(:,:),s(:,:),gt(:,:),winv(:,:),g(:,:),a(:,:),ainv(:,:),mid(:,:)
type(optim_result_t)::op
meth=LINEAR_TWO_STEP
if(present(method)) meth=method
covcode=1
if(present(covariance)) covcode=covariance
mo=100
if(present(maxouter)) mo=maxouter
eps=1.0e-7_dp
if(present(tol)) eps=tol
ctr=.true.
if(present(centered)) ctr=centered
ident=.false.
if(present(identity_weight)) ident=identity_weight
kern='Quadratic Spectral'
if(present(kernel)) kern=kernel
res%n=size(y)
res%k=size(x,2)
res%q=size(z,2)
res%df=res%q-res%k
allocate(w(res%q,res%q))
w=0
call set_identity(w)
call weighted_iv(y,x,z,w,theta,info)
res%initial=theta
bw=0.0_dp
if(covcode==3) then
   gt=linear_moments(theta,y,x,z)
   if(present(bandwidth)) then
   bw=bandwidth
   else
   bw=bw_andrews(gt,kern,'AR(1)')
   end if
end if
if(ident .or. res%q==res%k) then
   res%outer_iterations=1
else if(meth==LINEAR_TWO_STEP) then
   gt=linear_moments(theta,y,x,z)
   s=linear_cov(gt,y-matmul(x,theta),z,covcode,bw,kern,ctr)
   winv=spd_inverse(s)
   call weighted_iv(y,x,z,winv,newtheta,info)
   theta=newtheta
   res%outer_iterations=2
else if(meth==LINEAR_ITERATIVE) then
   do o=1,mo
      gt=linear_moments(theta,y,x,z)
      s=linear_cov(gt,y-matmul(x,theta),z,covcode,bw,kern,ctr)
      winv=spd_inverse(s)
      call weighted_iv(y,x,z,winv,newtheta,info)
      delta=maxval(abs(newtheta-theta))
      theta=newtheta
      if(delta<=eps) exit
   end do
   res%outer_iterations=o
else
   cue_y=y
   cue_x=x
   cue_z=z
   cue_cov=covcode
   cue_bw=bw
   cue_kernel=kern
   cue_centered=ctr
   op=optim_bfgs(linear_cue_obj,theta,maxit=500,reltol=1.0e-11_dp)
   theta=op%par
   res%convergence=op%convergence
   res%iterations=op%counts(1)
   res%outer_iterations=1
   if(allocated(cue_y)) deallocate(cue_y,cue_x,cue_z)
end if
res%coefficients=theta
res%fitted=matmul(x,theta)
res%residuals=y-res%fitted
res%moments=linear_moments(theta,y,x,z)
res%gradient=linear_gradient(x,z)
res%bandwidth=bw
res%omega=linear_cov(res%moments,res%residuals,z,covcode,bw,kern,ctr)
if(ident) then
   allocate(res%weight(res%q,res%q))
   res%weight=0
   call set_identity(res%weight)
else
   res%weight=spd_inverse(res%omega)
end if
block
   real(dp),allocatable::gb(:)
   allocate(gb(res%q))
   gb=sum(res%moments,dim=1)/real(res%n,dp)
   res%objective=dot_product(gb,matmul(res%weight,gb))
   res%j_stat=real(res%n,dp)*res%objective
end block
if(res%df>0) res%j_pvalue=1.0_dp-pchisq(res%j_stat,real(res%df,dp))
g=res%gradient
a=matmul(transpose(g),matmul(res%weight,g))
call invert_matrix(a,ainv,info)
if(info==0) then
   if(.not.ident) then
      res%vcov=ainv/real(res%n,dp)
   else
      mid=matmul(transpose(g),matmul(res%weight,matmul(res%omega,matmul(res%weight,g))))
      res%vcov=matmul(ainv,matmul(mid,ainv))/real(res%n,dp)
   end if
else
   allocate(res%vcov(res%k,res%k))
   res%vcov=huge(1.0_dp)
end if
end subroutine linear_gmm_fit

subroutine tsls_fit(y,x,z,res)
real(dp),intent(in)::y(:),x(:,:),z(:,:)
type(linear_gmm_result_t),intent(out)::res
real(dp),allocatable::zz(:,:),w(:,:),beta(:)
integer::info
zz=matmul(transpose(z),z)/real(size(z,1),dp)
w=spd_inverse(zz)
call weighted_iv(y,x,z,w,beta,info)
res%n=size(y)
res%k=size(x,2)
res%q=size(z,2)
res%df=res%q-res%k
res%coefficients=beta
res%fitted=matmul(x,beta)
res%residuals=y-res%fitted
res%moments=linear_moments(beta,y,x,z)
res%gradient=linear_gradient(x,z)
res%omega=zz*sum((res%residuals-sum(res%residuals)/real(res%n,dp))**2)/real(res%n,dp)
res%weight=w
block
 real(dp),allocatable::a(:,:),ainv(:,:)
 real(dp)::sig2
 integer::ii
 sig2=dot_product(res%residuals,res%residuals)/real(max(1,res%n-res%k),dp)
 a=matmul(transpose(x),matmul(z,matmul(w,matmul(transpose(z),x))))
 call invert_matrix(a,ainv,ii)
 if(ii==0) then
 res%vcov=sig2*ainv
 else
 allocate(res%vcov(res%k,res%k))
 res%vcov=huge(1.0_dp)
 end if
end block
end subroutine tsls_fit

subroutine weighted_iv(y,x,z,w,theta,info)
real(dp),intent(in)::y(:),x(:,:),z(:,:),w(:,:)
real(dp),allocatable,intent(out)::theta(:)
integer,intent(out)::info
real(dp),allocatable::zx(:,:),zy(:),a(:,:),b(:),ainv(:,:)
zx=matmul(transpose(z),x)/real(size(y),dp)
zy=matmul(transpose(z),y)/real(size(y),dp)
a=matmul(transpose(zx),matmul(w,zx))
b=matmul(transpose(zx),matmul(w,zy))
call invert_matrix(a,ainv,info)
if(info==0) then
theta=matmul(ainv,b)
else
allocate(theta(size(x,2)))
theta=0.0_dp
end if
end subroutine weighted_iv

function linear_cov(gt,e,z,covcode,bw,kernel,ctr) result(s)
real(dp),intent(in)::gt(:,:),e(:),z(:,:),bw
integer,intent(in)::covcode
character(len=*),intent(in)::kernel
logical,intent(in)::ctr
real(dp),allocatable::s(:,:)
real(dp)::ec(size(e)),sig
select case(covcode)
case(2) ! iid: upstream linear one-equation specialization
   ec=e-sum(e)/real(size(e),dp)
   sig=dot_product(ec,ec)/real(size(e),dp)
   s=sig*matmul(transpose(z),z)/real(size(e),dp)
case(3)
   s=hac_covariance(gt,bw,kernel,ctr)
case default
   s=moment_covariance(gt,ctr)
end select
end function linear_cov

pure function linear_cue_obj(theta) result(v)
real(dp),intent(in)::theta(:)
real(dp)::v
real(dp),allocatable::gt(:,:),s(:,:),w(:,:)
real(dp)::gb(size(cue_z,2)),e(size(cue_y))
e=cue_y-matmul(cue_x,theta)
gt=cue_z*spread(e,2,size(cue_z,2))
gb=sum(gt,dim=1)/real(size(gt,1),dp)
select case(cue_cov)
case(2)
   s=sum((e-sum(e)/real(size(e),dp))**2)/real(size(e),dp)*matmul(transpose(cue_z),cue_z)/real(size(e),dp)
case(3)
   s=hac_covariance(gt,cue_bw,cue_kernel,cue_centered)
case default
   s=moment_covariance(gt,cue_centered)
end select
w=spd_inverse(s)
v=dot_product(gb,matmul(w,gb))
end function linear_cue_obj

pure function spd_inverse(a) result(ai)
real(dp),intent(in)::a(:,:)
real(dp),allocatable::ai(:,:),r(:,:),aa(:,:)
real(dp)::ridge
integer::i
allocate(aa(size(a,1),size(a,2)))
aa=0.5_dp*(a+transpose(a))
ridge=1.0e-12_dp*max(1.0_dp,maxval(abs(aa)))
do i=1,size(aa,1)
aa(i,i)=aa(i,i)+ridge
end do
r=chol(aa)
ai=chol2inv(r)
end function spd_inverse

pure subroutine set_identity(a)
real(dp),intent(inout)::a(:,:)
integer::i
a=0.0_dp
do i=1,min(size(a,1),size(a,2))
a(i,i)=1.0_dp
end do
end subroutine set_identity

end module gmm_linear
