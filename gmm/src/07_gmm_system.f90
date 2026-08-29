! SPDX-License-Identifier: GPL-2.0-or-later
module gmm_system
use r_compat, only: dp, chol, chol2inv, pchisq
use gmm_linalg, only: invert_matrix
use gmm_covariance, only: moment_covariance, hac_covariance, bw_andrews
implicit none
private
public :: system_gmm_result_t, system_gmm_fit, sur_fit, three_sls_fit, five_fit, random_effect_fit
public :: SYS_MDS, SYS_HAC, SYS_CONDHOM
integer,parameter :: SYS_MDS=1,SYS_HAC=2,SYS_CONDHOM=3

type :: system_gmm_result_t
   real(dp),allocatable :: coefficients(:),initial(:),fitted(:,:),residuals(:,:),moments(:,:),gradient(:,:)
   real(dp),allocatable :: omega(:,:),weight(:,:),vcov(:,:),sigma(:,:)
   real(dp)::objective=0.0_dp,j_stat=0.0_dp,j_pvalue=1.0_dp,bandwidth=0.0_dp
   integer::n=0,neq=0,k=0,q=0,npar=0,df=0,convergence=0
   logical::common_coef=.false.
end type

contains

subroutine system_gmm_fit(y,x,z,res,covariance,common_coef,kernel,bandwidth,centered,identity_weight)
real(dp),intent(in)::y(:,:),x(:,:,:),z(:,:,:)
type(system_gmm_result_t),intent(out)::res
integer,intent(in),optional::covariance
logical,intent(in),optional::common_coef,centered,identity_weight
character(len=*),intent(in),optional::kernel
real(dp),intent(in),optional::bandwidth
integer::covcode,info
logical::common,ctr,ident
character(len=32)::kern
real(dp)::bw
real(dp),allocatable::d(:,:),b(:),theta(:),theta2(:),w(:,:),omega(:,:),ainv(:,:),a(:,:)
real(dp),allocatable::gt(:,:),g(:,:),mid(:,:)
covcode=SYS_CONDHOM
if(present(covariance)) covcode=covariance
common=.false.
if(present(common_coef)) common=common_coef
ctr=.true.
if(present(centered)) ctr=centered
ident=.false.
if(present(identity_weight)) ident=identity_weight
kern='Quadratic Spectral'
if(present(kernel)) kern=kernel
res%n=size(y,1)
res%neq=size(y,2)
res%k=size(x,2)
res%q=size(z,2)
res%common_coef=common
if(common) then
res%npar=res%k
else
res%npar=res%k*res%neq
end if
res%df=res%q*res%neq-res%npar
call system_design(y,x,z,common,d,b)
allocate(w(size(d,1),size(d,1)))
w=0
call set_identity(w)
call weighted_system(d,b,w,theta,info)
res%initial=theta
call system_outputs(theta,y,x,z,common,res%fitted,res%residuals,gt,g)
bw=0.0_dp
if(covcode==SYS_HAC) then
if(present(bandwidth)) then
bw=bandwidth
else
bw=bw_andrews(gt,kern,'AR(1)')
end if
end if
if(.not.ident .and. size(d,1)>size(d,2)) then
   omega=system_cov(gt,res%residuals,z,covcode,bw,kern,ctr)
   w=spd_inverse(omega)
   call weighted_system(d,b,w,theta2,info)
   theta=theta2
   call system_outputs(theta,y,x,z,common,res%fitted,res%residuals,gt,g)
end if
res%coefficients=theta
res%moments=gt
res%gradient=g
res%bandwidth=bw
res%sigma=matmul(transpose(res%residuals-spread(sum(res%residuals,dim=1)/real(res%n,dp),1,res%n)), &
                 res%residuals-spread(sum(res%residuals,dim=1)/real(res%n,dp),1,res%n))/real(res%n,dp)
res%omega=system_cov(gt,res%residuals,z,covcode,bw,kern,ctr)
if(ident) then
allocate(res%weight(size(res%omega,1),size(res%omega,2)))
res%weight=0
call set_identity(res%weight)
else
res%weight=spd_inverse(res%omega)
end if
block
 real(dp),allocatable::gb(:)
 allocate(gb(size(gt,2)))
 gb=sum(gt,dim=1)/real(res%n,dp)
 res%objective=dot_product(gb,matmul(res%weight,gb))
 res%j_stat=real(res%n,dp)*res%objective
end block
if(res%df>0)res%j_pvalue=1.0_dp-pchisq(res%j_stat,real(res%df,dp))
a=matmul(transpose(g),matmul(res%weight,g))
call invert_matrix(a,ainv,info)
if(info==0)then
 if(.not.ident)then
 res%vcov=ainv/real(res%n,dp)
 else
   mid=matmul(transpose(g),matmul(res%weight,matmul(res%omega,matmul(res%weight,g))))
   res%vcov=matmul(ainv,matmul(mid,ainv))/real(res%n,dp)
 end if
else
allocate(res%vcov(res%npar,res%npar))
res%vcov=huge(1.0_dp)
end if
end subroutine system_gmm_fit

subroutine sur_fit(y,x,res,common_coef)
real(dp),intent(in)::y(:,:),x(:,:,:)
type(system_gmm_result_t),intent(out)::res
logical,intent(in),optional::common_coef
real(dp),allocatable::z(:,:,:)
logical::cc
cc=.false.
if(present(common_coef))cc=common_coef
z=x
call system_gmm_fit(y,x,z,res,covariance=SYS_CONDHOM,common_coef=cc)
end subroutine sur_fit

subroutine three_sls_fit(y,x,z_common,res,common_coef)
real(dp),intent(in)::y(:,:),x(:,:,:),z_common(:,:)
type(system_gmm_result_t),intent(out)::res
logical,intent(in),optional::common_coef
real(dp),allocatable::z(:,:,:)
integer::j
allocate(z(size(z_common,1),size(z_common,2),size(y,2)))
do j=1,size(y,2)
z(:,:,j)=z_common
end do
if (present(common_coef)) then
   call system_gmm_fit(y,x,z,res,SYS_CONDHOM,common_coef)
else
   call system_gmm_fit(y,x,z,res,SYS_CONDHOM)
end if
end subroutine three_sls_fit

subroutine five_fit(y,x,z,res,common_coef)
real(dp),intent(in)::y(:,:),x(:,:,:),z(:,:,:)
type(system_gmm_result_t),intent(out)::res
logical,intent(in),optional::common_coef
if (present(common_coef)) then
   call system_gmm_fit(y,x,z,res,SYS_CONDHOM,common_coef)
else
   call system_gmm_fit(y,x,z,res,SYS_CONDHOM)
end if
end subroutine five_fit

subroutine random_effect_fit(y,x,res)
real(dp),intent(in)::y(:,:),x(:,:,:)
type(system_gmm_result_t),intent(out)::res
call sur_fit(y,x,res,.true.)
end subroutine random_effect_fit

subroutine system_design(y,x,z,common,d,b)
real(dp),intent(in)::y(:,:),x(:,:,:),z(:,:,:)
logical,intent(in)::common
real(dp),allocatable,intent(out)::d(:,:),b(:)
integer::m,k,q,j,r1,r2,c1,c2
m=size(y,2)
k=size(x,2)
q=size(z,2)
if(common)then
allocate(d(m*q,k))
else
allocate(d(m*q,m*k))
end if
allocate(b(m*q))
d=0
b=0
do j=1,m
 r1=(j-1)*q+1
 r2=j*q
 b(r1:r2)=matmul(transpose(z(:,:,j)),y(:,j))/real(size(y,1),dp)
 if(common)then
   d(r1:r2,:)=matmul(transpose(z(:,:,j)),x(:,:,j))/real(size(y,1),dp)
 else
   c1=(j-1)*k+1
   c2=j*k
   d(r1:r2,c1:c2)=matmul(transpose(z(:,:,j)),x(:,:,j))/real(size(y,1),dp)
 end if
end do
end subroutine system_design

subroutine weighted_system(d,b,w,theta,info)
real(dp),intent(in)::d(:,:),b(:),w(:,:)
real(dp),allocatable,intent(out)::theta(:)
integer,intent(out)::info
real(dp),allocatable::a(:,:),rhs(:),ai(:,:)
a=matmul(transpose(d),matmul(w,d))
rhs=matmul(transpose(d),matmul(w,b))
call invert_matrix(a,ai,info)
if(info==0)then
theta=matmul(ai,rhs)
else
allocate(theta(size(d,2)))
theta=0
end if
end subroutine weighted_system

subroutine system_outputs(theta,y,x,z,common,fitted,resid,gt,g)
real(dp),intent(in)::theta(:),y(:,:),x(:,:,:),z(:,:,:)
logical,intent(in)::common
real(dp),allocatable,intent(out)::fitted(:,:),resid(:,:),gt(:,:),g(:,:)
integer::m,n,k,q,j,r1,r2,c1,c2
real(dp),allocatable::tj(:)
n=size(y,1)
m=size(y,2)
k=size(x,2)
q=size(z,2)
allocate(fitted(n,m),resid(n,m),gt(n,m*q))
gt=0
if(common)then
allocate(g(m*q,k))
g=0
else
allocate(g(m*q,m*k))
g=0
end if
do j=1,m
 if(common)then
 tj=theta
 else
 c1=(j-1)*k+1
 c2=j*k
 tj=theta(c1:c2)
 end if
 fitted(:,j)=matmul(x(:,:,j),tj)
 resid(:,j)=y(:,j)-fitted(:,j)
 r1=(j-1)*q+1
 r2=j*q
 gt(:,r1:r2)=z(:,:,j)*spread(resid(:,j),2,q)
 if(common)then
 g(r1:r2,:)=-matmul(transpose(z(:,:,j)),x(:,:,j))/real(n,dp)
 else
 g(r1:r2,c1:c2)=-matmul(transpose(z(:,:,j)),x(:,:,j))/real(n,dp)
 end if
end do
end subroutine system_outputs

function system_cov(gt,resid,z,covcode,bw,kernel,ctr) result(s)
real(dp),intent(in)::gt(:,:),resid(:,:),z(:,:,:),bw
integer,intent(in)::covcode
character(len=*),intent(in)::kernel
logical,intent(in)::ctr
real(dp),allocatable::s(:,:)
real(dp),allocatable::ec(:,:),sig(:,:)
integer::m,q,i,j,r1,r2,c1,c2,n
n=size(gt,1)
m=size(resid,2)
q=size(z,2)
select case(covcode)
case(SYS_HAC);s=hac_covariance(gt,bw,kernel,ctr)
case(SYS_CONDHOM)
 ec=resid-spread(sum(resid,dim=1)/real(n,dp),1,n)
 sig=matmul(transpose(ec),ec)/real(n,dp)
 allocate(s(m*q,m*q))
 s=0
 do i=1,m
 do j=1,m
   r1=(i-1)*q+1
   r2=i*q
   c1=(j-1)*q+1
   c2=j*q
   s(r1:r2,c1:c2)=sig(i,j)*matmul(transpose(z(:,:,i)),z(:,:,j))/real(n,dp)
 end do
 end do
case default;s=moment_covariance(gt,ctr)
end select
end function system_cov

pure function spd_inverse(a) result(ai)
use r_compat, only: chol,chol2inv
real(dp),intent(in)::a(:,:)
real(dp),allocatable::ai(:,:),aa(:,:),r(:,:)
integer::i
allocate(aa(size(a,1),size(a,2)))
aa=0.5_dp*(a+transpose(a))
do i=1,size(aa,1)
aa(i,i)=aa(i,i)+1.0e-12_dp*max(1.0_dp,maxval(abs(aa)))
end do
r=chol(aa)
ai=chol2inv(r)
end function spd_inverse

pure subroutine set_identity(a)
real(dp),intent(inout)::a(:,:)
integer::i
a=0
do i=1,min(size(a,1),size(a,2))
a(i,i)=1
end do
end subroutine
end module gmm_system
