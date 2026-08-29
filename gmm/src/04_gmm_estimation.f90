! SPDX-License-Identifier: GPL-2.0-or-later
module gmm_estimation
use r_compat, only: dp, optim_result_t, optim_bfgs, chol, chol2inv, pchisq
use gmm_linalg, only: colmeans_mat, invert_matrix
use gmm_covariance, only: moment_covariance, hac_covariance, bw_andrews
implicit none
private
public :: gmm_fit_result_t, gmm_fit, gmm_evaluate, gmm_fit_fixed_weight, gmm_moment_function, gmm_gradient_function
public :: GMM_TWO_STEP, GMM_ITERATIVE, GMM_CUE, COV_MDS, COV_IID, COV_HAC, COV_IDENT

integer,parameter :: GMM_TWO_STEP=1,GMM_ITERATIVE=2,GMM_CUE=3
integer,parameter :: COV_MDS=1,COV_IID=2,COV_HAC=3,COV_IDENT=4

type :: gmm_fit_result_t
   real(dp), allocatable :: theta(:), theta_initial(:), moments(:,:), gradient(:,:), omega(:,:), weight(:,:), vcov(:,:)
   real(dp) :: objective=0.0_dp, j_stat=0.0_dp, j_pvalue=1.0_dp, bandwidth=0.0_dp
   integer :: n=0,q=0,k=0,df=0,convergence=0,iterations=0,outer_iterations=0
end type gmm_fit_result_t

abstract interface
   pure function gmm_moment_function(theta,data) result(gt)
      import :: dp
      real(dp),intent(in)::theta(:),data(:,:)
      real(dp),allocatable::gt(:,:)
   end function gmm_moment_function
   pure function gmm_gradient_function(theta,data) result(g)
      import :: dp
      real(dp),intent(in)::theta(:),data(:,:)
      real(dp),allocatable::g(:,:)
   end function gmm_gradient_function
end interface

real(dp),allocatable,save :: opt_data(:,:), opt_weight(:,:)
procedure(gmm_moment_function),pointer,save :: opt_moment=>null()
procedure(gmm_gradient_function),pointer,save :: opt_gradient=>null()
integer,save :: opt_cov=COV_MDS
logical,save :: opt_cue=.false.,opt_have_grad=.false.
real(dp),save :: opt_bw=0.0_dp
character(len=32),save :: opt_kernel='Quadratic Spectral'
logical,save :: opt_centered=.true.

contains

subroutine gmm_fit(moment,data,theta0,res,method,covariance,kernel,bandwidth,centered,tol,maxit,maxouter,gradient)
procedure(gmm_moment_function) :: moment
real(dp),intent(in)::data(:,:),theta0(:)
type(gmm_fit_result_t),intent(out)::res
integer,intent(in),optional::method,covariance,maxit,maxouter
character(len=*),intent(in),optional::kernel
real(dp),intent(in),optional::bandwidth,tol
logical,intent(in),optional::centered
procedure(gmm_gradient_function),optional :: gradient
integer::meth,covcode,mi,mo,o
real(dp)::eps,delta,bw
real(dp),allocatable::theta(:),gt(:,:),omega(:,:),w(:,:)
type(optim_result_t)::op
logical::ctr
character(len=32)::kern
meth=GMM_TWO_STEP
if(present(method)) meth=method
covcode=COV_MDS
if(present(covariance)) covcode=covariance
mi=500
if(present(maxit)) mi=maxit
mo=100
if(present(maxouter)) mo=maxouter
eps=1.0e-7_dp
if(present(tol)) eps=tol
ctr=.true.
if(present(centered)) ctr=centered
kern='Quadratic Spectral'
if(present(kernel)) kern=kernel

theta=theta0
gt=moment(theta,data)
res%n=size(gt,1)
res%q=size(gt,2)
res%k=size(theta)
res%df=res%q-res%k
bw=0.0_dp
if(covcode==COV_HAC) then
   if(present(bandwidth)) then
      bw=bandwidth
   else
      bw=bw_andrews(gt,kern,'AR(1)')
   end if
end if
call setup_optimizer(moment,data,covcode,kern,bw,ctr,.false.,gradient)
allocate(opt_weight(res%q,res%q))
opt_weight=0.0_dp
call set_identity(opt_weight)
if(meth==GMM_CUE .and. covcode/=COV_IDENT) then
   ! Follow the package's CUE structure: obtain a stable first-stage start, then
   ! minimize with the covariance updated at each candidate theta.
   op=run_bfgs(theta,mi,.false.)
   theta=op%par
   if(covcode==COV_HAC .and. .not.present(bandwidth)) then
      gt=moment(theta,data)
      bw=bw_andrews(gt,kern,'AR(1)')
      opt_bw=bw
   end if
   opt_cue=.true.
   op=run_bfgs(theta,mi,.false.)
   theta=op%par
   res%outer_iterations=1
else
   op=run_bfgs(theta,mi,present(gradient))
   theta=op%par
   if(meth==GMM_TWO_STEP .and. covcode/=COV_IDENT) then
      gt=moment(theta,data)
      omega=make_cov(gt,covcode,bw,kern,ctr)
      w=spd_inverse(omega)
      opt_weight=w
      op=run_bfgs(theta,mi,present(gradient))
      theta=op%par
      res%outer_iterations=2
   else if(meth==GMM_ITERATIVE .and. covcode/=COV_IDENT) then
      do o=1,mo
         gt=moment(theta,data)
         omega=make_cov(gt,covcode,bw,kern,ctr)
         w=spd_inverse(omega)
         opt_weight=w
         res%theta_initial=theta
         op=run_bfgs(theta,mi,present(gradient))
         delta=maxval(abs(op%par-theta))
         theta=op%par
         if(delta<=eps) exit
      end do
      res%outer_iterations=o
   else
      res%outer_iterations=1
   end if
end if
res%theta=theta
res%convergence=op%convergence
res%iterations=op%counts(1)
res%bandwidth=bw
call finalize_result(moment,data,res,covcode,kern,bw,ctr,gradient,meth)
call clear_optimizer()
end subroutine gmm_fit

subroutine gmm_fit_fixed_weight(moment,data,theta0,weight,res,gradient,centered,true_fixed,maxit)
procedure(gmm_moment_function)::moment
real(dp),intent(in)::data(:,:),theta0(:),weight(:,:)
type(gmm_fit_result_t),intent(out)::res
procedure(gmm_gradient_function),optional::gradient
logical,intent(in),optional::centered,true_fixed
integer,intent(in),optional::maxit
logical::ctr,tf
integer::mi,info
type(optim_result_t)::op
real(dp),allocatable::gt(:,:),g(:,:),s(:,:),a(:,:),ai(:,:),t1(:,:),mid(:,:)
real(dp),allocatable::gb(:)
ctr=.true.
if(present(centered))ctr=centered
tf=.false.
if(present(true_fixed))tf=true_fixed
mi=500
if(present(maxit))mi=maxit
call setup_optimizer(moment,data,COV_MDS,'Quadratic Spectral',0.0_dp,ctr,.false.,gradient)
allocate(opt_weight(size(weight,1),size(weight,2)))
opt_weight=weight
op=run_bfgs(theta0,mi,present(gradient))
res%theta=op%par
res%convergence=op%convergence
res%iterations=op%counts(1)
gt=moment(res%theta,data)
res%moments=gt
res%n=size(gt,1)
res%q=size(gt,2)
res%k=size(res%theta)
res%df=res%q-res%k
allocate(gb(res%q))
gb=sum(gt,dim=1)/real(res%n,dp)
res%objective=dot_product(gb,matmul(weight,gb))
res%j_stat=real(res%n,dp)*res%objective
if(res%df>0)res%j_pvalue=1.0_dp-pchisq(res%j_stat,real(res%df,dp))
res%weight=weight
s=moment_covariance(gt,ctr)
res%omega=s
if(present(gradient))then
g=gradient(res%theta,data)
else
g=numeric_gradient(moment,res%theta,data)
end if
res%gradient=g
a=matmul(transpose(g),matmul(weight,g))
call invert_matrix(a,ai,info)
if(info==0)then
 if(tf)then
  res%vcov=ai/real(res%n,dp)
 else
  t1=matmul(ai,matmul(transpose(g),weight))
  mid=matmul(t1,matmul(s,transpose(t1)))
  res%vcov=mid/real(res%n,dp)
 end if
else
 allocate(res%vcov(res%k,res%k))
 res%vcov=huge(1.0_dp)
end if
call clear_optimizer()
end subroutine gmm_fit_fixed_weight

subroutine gmm_evaluate(moment,data,theta,res,covariance,kernel,bandwidth,centered,gradient)
procedure(gmm_moment_function)::moment
real(dp),intent(in)::data(:,:),theta(:)
type(gmm_fit_result_t),intent(out)::res
integer,intent(in),optional::covariance
character(len=*),intent(in),optional::kernel
real(dp),intent(in),optional::bandwidth
logical,intent(in),optional::centered
procedure(gmm_gradient_function),optional::gradient
integer::covcode
real(dp)::bw
logical::ctr
character(len=32)::kern
real(dp),allocatable::gt(:,:)
covcode=COV_MDS
if(present(covariance)) covcode=covariance
ctr=.true.
if(present(centered)) ctr=centered
kern='Quadratic Spectral'
if(present(kernel)) kern=kernel
gt=moment(theta,data)
bw=0.0_dp
if(covcode==COV_HAC) then
if(present(bandwidth)) then
bw=bandwidth
else
bw=bw_andrews(gt,kern,'AR(1)')
end if
end if
res%theta=theta
res%n=size(gt,1)
res%q=size(gt,2)
res%k=size(theta)
res%df=res%q-res%k
res%bandwidth=bw
call finalize_result(moment,data,res,covcode,kern,bw,ctr,gradient,GMM_TWO_STEP)
end subroutine gmm_evaluate

subroutine setup_optimizer(moment,data,covcode,kernel,bw,ctr,cue,gradient)
procedure(gmm_moment_function)::moment
real(dp),intent(in)::data(:,:),bw
integer,intent(in)::covcode
character(len=*),intent(in)::kernel
logical,intent(in)::ctr,cue
procedure(gmm_gradient_function),optional::gradient
if(allocated(opt_data)) deallocate(opt_data)
opt_data=data
opt_moment=>moment
opt_cov=covcode
opt_kernel=kernel
opt_bw=bw
opt_centered=ctr
opt_cue=cue
opt_have_grad=.false.
nullify(opt_gradient)
if(present(gradient)) then
opt_gradient=>gradient
opt_have_grad=.true.
end if
end subroutine setup_optimizer

subroutine clear_optimizer()
if(allocated(opt_data)) deallocate(opt_data)
if(allocated(opt_weight)) deallocate(opt_weight)
nullify(opt_moment)
nullify(opt_gradient)
opt_have_grad=.false.
opt_cue=.false.
end subroutine clear_optimizer

function run_bfgs(theta,maxit,use_grad) result(op)
real(dp),intent(in)::theta(:)
integer,intent(in)::maxit
logical,intent(in)::use_grad
type(optim_result_t)::op
if(use_grad .and. opt_have_grad .and. .not.opt_cue) then
   op=optim_bfgs(gmm_obj,theta,maxit=maxit,reltol=1.0e-10_dp,gr=gmm_obj_grad)
else
   op=optim_bfgs(gmm_obj,theta,maxit=maxit,reltol=1.0e-10_dp)
end if
end function run_bfgs

pure function gmm_obj(theta) result(v)
real(dp),intent(in)::theta(:)
real(dp)::v
real(dp),allocatable::gt(:,:),om(:,:),w(:,:)
real(dp)::gb(size(opt_weight,1))
gt=opt_moment(theta,opt_data)
gb=sum(gt,dim=1)/real(size(gt,1),dp)
if(opt_cue) then
   om=make_cov_pure(gt,opt_cov,opt_bw,opt_kernel,opt_centered)
   w=spd_inverse(om)
   v=dot_product(gb,matmul(w,gb))
else
   v=dot_product(gb,matmul(opt_weight,gb))
end if
if(.not.(v<huge(1.0_dp))) v=huge(1.0_dp)/100.0_dp
end function gmm_obj

pure function gmm_obj_grad(theta) result(gv)
real(dp),intent(in)::theta(:)
real(dp),allocatable::gv(:)
real(dp),allocatable::gt(:,:),g(:,:)
real(dp)::gb(size(opt_weight,1))
gt=opt_moment(theta,opt_data)
gb=sum(gt,dim=1)/real(size(gt,1),dp)
g=opt_gradient(theta,opt_data)
allocate(gv(size(theta)))
gv=2.0_dp*matmul(transpose(g),matmul(opt_weight,gb))
end function gmm_obj_grad

subroutine finalize_result(moment,data,res,covcode,kernel,bw,ctr,gradient,meth)
procedure(gmm_moment_function)::moment
real(dp),intent(in)::data(:,:),bw
integer,intent(in)::covcode,meth
character(len=*),intent(in)::kernel
logical,intent(in)::ctr
procedure(gmm_gradient_function),optional::gradient
type(gmm_fit_result_t),intent(inout)::res
real(dp),allocatable::s(:,:),w(:,:),g(:,:),a(:,:),ainv(:,:),mid(:,:),gt(:,:)
real(dp)::gb(res%q)
integer::info
gt=moment(res%theta,data)
res%moments=gt
gb=sum(gt,dim=1)/real(res%n,dp)
s=make_cov(gt,covcode,bw,kernel,ctr)
res%omega=s
if(covcode==COV_IDENT) then
   allocate(w(res%q,res%q))
   w=0
   call set_identity(w)
else
   w=spd_inverse(s)
end if
res%weight=w
res%objective=dot_product(gb,matmul(w,gb))
res%j_stat=real(res%n,dp)*res%objective
if(res%df>0) then
res%j_pvalue=1.0_dp-pchisq(res%j_stat,real(res%df,dp))
else
res%j_pvalue=1.0_dp
end if
if(present(gradient)) then
g=gradient(res%theta,data)
else
g=numeric_gradient(moment,res%theta,data)
end if
res%gradient=g
a=matmul(transpose(g),matmul(w,g))
call invert_matrix(a,ainv,info)
if(info/=0) then
   allocate(res%vcov(res%k,res%k))
   res%vcov=huge(1.0_dp)
else
   if(covcode/=COV_IDENT .and. meth/=GMM_CUE) then
      ! Optimal GMM simplification, matching FinRes when W=S^{-1}.
      res%vcov=ainv/real(res%n,dp)
   else
      mid=matmul(transpose(g),matmul(w,matmul(s,matmul(w,g))))
      res%vcov=matmul(ainv,matmul(mid,ainv))/real(res%n,dp)
   end if
end if
end subroutine finalize_result

function numeric_gradient(moment,theta,data) result(g)
procedure(gmm_moment_function)::moment
real(dp),intent(in)::theta(:),data(:,:)
real(dp),allocatable::g(:,:)
real(dp),allocatable::gp(:,:),gm(:,:)
real(dp)::tp(size(theta)),tm(size(theta)),h
integer::j,q
q=size(moment(theta,data),2)
allocate(g(q,size(theta)))
do j=1,size(theta)
   h=epsilon(1.0_dp)**(1.0_dp/3.0_dp)*max(1.0_dp,abs(theta(j)))
   tp=theta
   tm=theta
   tp(j)=tp(j)+h
   tm(j)=tm(j)-h
   gp=moment(tp,data)
   gm=moment(tm,data)
   g(:,j)=(sum(gp,dim=1)-sum(gm,dim=1))/(2.0_dp*h*real(size(data,1),dp))
end do
end function numeric_gradient

function make_cov(gt,covcode,bw,kernel,ctr) result(s)
real(dp),intent(in)::gt(:,:),bw
integer,intent(in)::covcode
character(len=*),intent(in)::kernel
logical,intent(in)::ctr
real(dp),allocatable::s(:,:)
allocate(s(size(gt,2),size(gt,2)))
select case(covcode)
case(COV_HAC); s=hac_covariance(gt,bw,kernel,ctr)
case(COV_IDENT)
s=0
call set_identity(s)
case default; s=moment_covariance(gt,ctr)
end select
end function make_cov

pure function make_cov_pure(gt,covcode,bw,kernel,ctr) result(s)
real(dp),intent(in)::gt(:,:),bw
integer,intent(in)::covcode
character(len=*),intent(in)::kernel
logical,intent(in)::ctr
real(dp),allocatable::s(:,:)
allocate(s(size(gt,2),size(gt,2)))
select case(covcode)
case(COV_HAC); s=hac_covariance(gt,bw,kernel,ctr)
case(COV_IDENT)
s=0
call set_identity(s)
case default; s=moment_covariance(gt,ctr)
end select
end function make_cov_pure

pure function spd_inverse(a) result(ai)
real(dp),intent(in)::a(:,:)
real(dp),allocatable::ai(:,:),r(:,:),aa(:,:)
real(dp)::ridge
integer::i
allocate(aa(size(a,1),size(a,2)))
aa=0.5_dp*(a+transpose(a))
ridge=max(1.0e-12_dp,epsilon(1.0_dp)*max(1.0_dp,maxval(abs(aa))))
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

end module gmm_estimation
