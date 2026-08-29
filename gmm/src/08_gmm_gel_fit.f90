! SPDX-License-Identifier: GPL-2.0-or-later
module gmm_gel_fit
use r_compat, only: dp, optim_result_t, optim_bfgs, chol, chol2inv, pchisq
use gmm_estimation, only: gmm_moment_function, gmm_gradient_function
use gmm_gel, only: GEL_EL,GEL_ET,GEL_CUE,GEL_ETEL,GEL_HD,GEL_ETHD,GEL_RCUE,gel_rho,gel_implied_prob
implicit none
private
public :: gel_fit_result_t, gel_fit, gel_evaluate

type :: gel_fit_result_t
   real(dp),allocatable :: coefficients(:),lambda(:),prob(:),moments(:,:),gradient(:,:)
   real(dp),allocatable :: vcov_par(:,:),vcov_lambda(:,:),khat(:,:)
   real(dp)::objective=0.0_dp,lr_stat=0.0_dp,lm_stat=0.0_dp,j_stat=0.0_dp
   real(dp)::lr_pvalue=1.0_dp,lm_pvalue=1.0_dp,j_pvalue=1.0_dp,k_scale=1.0_dp,bandwidth=1.0_dp
   integer::n=0,q=0,k=0,df=0,conv_par=0,conv_lambda=0,iterations=0
end type

real(dp),allocatable,save :: gel_data(:,:)
procedure(gmm_moment_function),pointer,save :: gel_moment=>null()
integer,save :: gel_type=GEL_EL
real(dp),save :: gel_kscale=1.0_dp
integer,save :: gel_inner_maxit=100
real(dp),save :: gel_inner_tol=1.0e-9_dp

contains

subroutine gel_fit(moment,data,theta0,type,res,gradient,tol,maxit,maxiterlam,k_scale,bandwidth)
procedure(gmm_moment_function)::moment
real(dp),intent(in)::data(:,:),theta0(:)
integer,intent(in)::type
type(gel_fit_result_t),intent(out)::res
procedure(gmm_gradient_function),optional::gradient
real(dp),intent(in),optional::tol,k_scale,bandwidth
integer,intent(in),optional::maxit,maxiterlam
integer::mi
real(dp)::rtol
type(optim_result_t)::op
gel_data=data
gel_moment=>moment
gel_type=type
gel_kscale=1.0_dp
if(present(k_scale))gel_kscale=k_scale
gel_inner_maxit=100
if(present(maxiterlam))gel_inner_maxit=maxiterlam
gel_inner_tol=1.0e-9_dp
if(present(tol))gel_inner_tol=tol
mi=500
if(present(maxit))mi=maxit
rtol=1.0e-10_dp
if(present(tol))rtol=max(tol,1.0e-12_dp)
op=optim_bfgs(gel_theta_objective,theta0,maxit=mi,reltol=rtol)
res%coefficients=op%par
res%conv_par=op%convergence
res%iterations=op%counts(1)
res%k_scale=gel_kscale
res%bandwidth=1.0_dp
if(present(bandwidth))res%bandwidth=bandwidth
call finish_gel(moment,data,type,res,gradient)
if(allocated(gel_data))deallocate(gel_data)
nullify(gel_moment)
end subroutine gel_fit

subroutine gel_evaluate(moment,data,theta,type,res,gradient,k_scale,bandwidth)
procedure(gmm_moment_function)::moment
real(dp),intent(in)::data(:,:),theta(:)
integer,intent(in)::type
type(gel_fit_result_t),intent(out)::res
procedure(gmm_gradient_function),optional::gradient
real(dp),intent(in),optional::k_scale,bandwidth
res%coefficients=theta
res%k_scale=1.0_dp
if(present(k_scale))res%k_scale=k_scale
res%bandwidth=1.0_dp
if(present(bandwidth))res%bandwidth=bandwidth
call finish_gel(moment,data,type,res,gradient)
end subroutine gel_evaluate

subroutine finish_gel(moment,data,type,res,gradient)
procedure(gmm_moment_function)::moment
real(dp),intent(in)::data(:,:)
integer,intent(in)::type
type(gel_fit_result_t),intent(inout)::res
procedure(gmm_gradient_function),optional::gradient
real(dp),allocatable::g(:,:),lambda(:),gt(:,:)
integer::conv
real(dp)::obj
gt=moment(res%coefficients,data)
res%moments=gt
res%n=size(gt,1)
res%q=size(gt,2)
res%k=size(res%coefficients)
res%df=res%q-res%k
call gel_eval_inner(gt,type,res%k_scale,lambda,obj,conv)
res%lambda=lambda
res%objective=obj
res%conv_lambda=conv
res%prob=gel_implied_prob(gt,lambda,type,res%k_scale)
if(present(gradient))then
g=gradient(res%coefficients,data)
else
g=numeric_gradient(moment,res%coefficients,data)
end if
res%gradient=g
call gel_covariance(gt,g,res%prob,res%k_scale,res%bandwidth,res%vcov_par,res%vcov_lambda,res%khat)
call gel_spec_stats(res,type)
end subroutine finish_gel

pure function gel_theta_objective(theta) result(obj)
real(dp),intent(in)::theta(:)
real(dp)::obj
real(dp),allocatable::gt(:,:),lambda(:)
integer::conv
gt=gel_moment(theta,gel_data)
call gel_eval_inner(gt,gel_type,gel_kscale,lambda,obj,conv)
if(conv/=0 .or. .not.(obj<huge(1.0_dp)))obj=huge(1.0_dp)/100.0_dp
end function gel_theta_objective

pure subroutine gel_eval_inner(gt,type,k,lambda,obj,conv)
real(dp),intent(in)::gt(:,:),k
integer,intent(in)::type
real(dp),allocatable,intent(out)::lambda(:)
real(dp),intent(out)::obj
integer,intent(out)::conv
real(dp),allocatable::r(:)
integer::st
st=type
if(type==GEL_ETEL.or.type==GEL_ETHD)st=GEL_ET
if(type==GEL_CUE)then
   call cue_lambda_pure(gt,lambda,conv)
else if(type==GEL_RCUE)then
   call rcue_lambda_pure(gt,k,lambda,obj,conv)
   return
else
   call newton_lambda_pure(gt,st,k,lambda,conv)
end if
if(conv/=0)then
obj=huge(1.0_dp)/100.0_dp
return
end if
r=gel_rho(gt,lambda,0,type,k)
select case(type)
case(GEL_EL,GEL_CUE);obj=sum(r)/real(size(r),dp)
case(GEL_ET);obj=sum(r)/real(size(r),dp)+1.0_dp
case(GEL_HD);obj=sum(r)/real(size(r),dp)+2.0_dp
case(GEL_ETEL);obj=sum(r)/real(size(r),dp)
case(GEL_ETHD);obj=sum(r)
end select
end subroutine gel_eval_inner

pure subroutine newton_lambda_pure(gt,type,k,lambda,conv)
real(dp),intent(in)::gt(:,:),k
integer,intent(in)::type
real(dp),allocatable,intent(out)::lambda(:)
integer,intent(out)::conv
real(dp),allocatable::r1(:),r2(:),f(:),jmat(:,:),invneg(:,:),step(:),trial(:)
real(dp)::oldnorm
integer::i,half
allocate(lambda(size(gt,2)))
lambda=0
conv=1
do i=1,gel_inner_maxit
   r1=gel_rho(gt,lambda,1,type,k)
   r2=gel_rho(gt,lambda,2,type,k)
   f=-sum(gt*spread(r1,2,size(gt,2)),dim=1)/real(size(gt,1),dp)
   oldnorm=sum(abs(f))
   if(oldnorm<gel_inner_tol)then
   conv=0
   return
   end if
   jmat=matmul(transpose(gt*spread(r2,2,size(gt,2))),gt)
   invneg=spd_inverse(-jmat)
   step=-matmul(invneg,f)
   if(sum(abs(step))<gel_inner_tol)then
   conv=0
   return
   end if
   trial=lambda+step
   if(type==GEL_EL)then
      do half=1,60
         if(all(1.0_dp-k*matmul(gt,trial)>0.0_dp))exit
         step=step/2.0_dp
         trial=lambda+step
      end do
      if(.not.all(1.0_dp-k*matmul(gt,trial)>0.0_dp))return
   end if
   lambda=trial
end do
end subroutine newton_lambda_pure

pure subroutine cue_lambda_pure(gt,lambda,conv)
real(dp),intent(in)::gt(:,:)
real(dp),allocatable,intent(out)::lambda(:)
integer,intent(out)::conv
real(dp),allocatable::a(:,:),rhs(:),ai(:,:)
a=matmul(transpose(gt),gt)
rhs=-sum(gt,dim=1)
ai=spd_inverse(a)
lambda=matmul(ai,rhs)
conv=0
end subroutine cue_lambda_pure

pure subroutine rcue_lambda_pure(gt,k,lambda,obj,conv)
real(dp),intent(in)::gt(:,:),k
real(dp),allocatable,intent(out)::lambda(:)
real(dp),intent(out)::obj
integer,intent(out)::conv
logical::keep(size(gt,1)),newkeep(size(gt,1))
real(dp),allocatable::sub(:,:),lam(:),r(:)
real(dp)::raw(size(gt,1))
integer::i,n1,n0,c
keep=.true.
conv=1
obj=0.0_dp
allocate(lambda(size(gt,2)))
lambda=0
DO i=1,gel_inner_maxit
   n1=count(keep)
   if(n1<size(gt,2)+1)then
   conv=2
   return
   end if
   sub=pack_rows(gt,keep)
   call cue_lambda_pure(sub,lam,c)
   if(c/=0)then
   conv=3
   return
   end if
   lambda=lam
   raw=1.0_dp+k*matmul(gt,lambda)
   newkeep=raw>0.0_dp
   if(all(newkeep.eqv.keep))then
      r=gel_rho(sub,lambda,0,GEL_CUE,k)
      n0=size(gt,1)-n1
      obj=sum(r)/real(size(gt,1),dp)+real(n0,dp)/(2.0_dp*real(size(gt,1),dp))
      conv=0
      return
   end if
   keep=newkeep
END DO
end subroutine rcue_lambda_pure

pure function pack_rows(a,mask) result(b)
real(dp),intent(in)::a(:,:)
logical,intent(in)::mask(:)
real(dp)::b(count(mask),size(a,2))
integer::i,j
j=0
do i=1,size(a,1)
if(mask(i))then
j=j+1
b(j,:)=a(i,:)
end if
end do
end function pack_rows

subroutine gel_covariance(gt,g,pt,k1,bw,vcpar,vclam,khat)
real(dp),intent(in)::gt(:,:),g(:,:),pt(:),k1,bw
real(dp),allocatable,intent(out)::vcpar(:,:),vclam(:,:),khat(:,:)
real(dp),allocatable::gw(:,:),kinv(:,:),a(:,:),ainv(:,:),proj(:,:)
integer::n,q,k
n=size(gt,1)
q=size(gt,2)
k=size(g,2)
gw=gt*spread(sqrt(max(pt,0.0_dp)*bw),2,q)
khat=matmul(transpose(gw),gw)
kinv=spd_inverse(khat)
a=matmul(transpose(g/k1),matmul(kinv,g/k1))
ainv=spd_inverse(a)
vcpar=ainv/real(n,dp)
! Algebraic equivalent of the projection used in .vcovGel.
proj=matmul(g/k1,matmul(ainv,matmul(transpose(g/k1),kinv)))
! Covariance of lambda: K^{-1}(I-G(G'K^{-1}G)^{-1}G'K^{-1}) / n, scaled as upstream.
vclam=matmul(kinv,identity_minus(proj))*bw*bw/real(n,dp)
vclam=0.5_dp*(vclam+transpose(vclam))
end subroutine gel_covariance

pure function identity_minus(a) result(b)
real(dp),intent(in)::a(:,:)
real(dp)::b(size(a,1),size(a,2))
integer::i
b=-a
do i=1,min(size(a,1),size(a,2))
b(i,i)=b(i,i)+1.0_dp
end do
end function identity_minus

subroutine gel_spec_stats(res,type)
type(gel_fit_result_t),intent(inout)::res
integer,intent(in)::type
real(dp)::gb(res%q)
gb=sum(res%moments,dim=1)/real(res%n,dp)
res%lr_stat=2.0_dp*res%objective*real(res%n,dp)/(res%bandwidth*res%k_scale**2)
if(type==GEL_ETHD)res%lr_stat=2.0_dp*res%lr_stat
res%lm_stat=real(res%n,dp)*dot_product(res%lambda,matmul(res%khat,res%lambda))/(res%bandwidth**2)
res%j_stat=real(res%n,dp)*dot_product(gb,matmul(spd_inverse(res%khat),gb))/(res%k_scale**2)
if(res%df>0)then
 res%lr_pvalue=1-pchisq(res%lr_stat,real(res%df,dp))
 res%lm_pvalue=1-pchisq(res%lm_stat,real(res%df,dp))
 res%j_pvalue=1-pchisq(res%j_stat,real(res%df,dp))
end if
end subroutine gel_spec_stats

function numeric_gradient(moment,theta,data) result(g)
procedure(gmm_moment_function)::moment
real(dp),intent(in)::theta(:),data(:,:)
real(dp),allocatable::g(:,:),gp(:,:),gm(:,:)
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
 g(:,j)=(sum(gp,dim=1)-sum(gm,dim=1))/(2*h*real(size(data,1),dp))
end do
end function numeric_gradient

pure function spd_inverse(a) result(ai)
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

end module gmm_gel_fit
