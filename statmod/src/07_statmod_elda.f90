module statmod_elda
use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf
use r_compat, only: dp, qnorm, pchisq
use statmod_linalg, only: weighted_least_squares, symmetric_inverse
implicit none
private
public :: elda_one_group_result_t, elda_result_t, elda_one_group, elda_fit
real(dp),parameter::pi=acos(-1.0_dp)

type :: elda_one_group_result_t
   real(dp),allocatable::p(:)
   real(dp)::lambda=0,alpha=0,ci_alpha(3)=0,ci_frequency(3)=0,deviance=0
   real(dp)::z_score_logdose=0,z_score_dose=0,dloglik_logdose=0,fisher_logdose=0
   real(dp)::dloglik_dose=0,fisher_dose=0
   integer::iter=0
end type

type :: elda_result_t
   real(dp),allocatable::ci(:,:) ! groups x lower/estimate/upper
   real(dp)::difference_chisq=0,difference_p=1
   integer::difference_df=0
   real(dp)::slope_estimate=1,slope_se=0,slope_wald_z=0,slope_wald_p=1
   real(dp)::slope_lr_z=0,slope_lr_p=1
   real(dp)::score_logdose_z=0,score_logdose_p=1,score_dose_z=0,score_dose_p=1
end type
contains

function limdil_allpos(tested,dose,confidence,observed) result(lambda)
real(dp),intent(in)::tested(:),dose(:),confidence
logical,intent(in)::observed
real(dp)::lambda,alpha,dosem,testedsum,beta,f,deriv,step
integer::iter
alpha=1-confidence
dosem=minval(dose)
testedsum=sum(tested,mask=dose==dosem)
beta=log(-log(1-alpha**(1/testedsum)))-log(dosem)
lambda=exp(beta)
if(observed)lambda=-expm1_local(lambda)
do iter=1,1000
   if(observed)then
      f=sum(tested*log(1-(1-lambda)**dose))-log(alpha)
      deriv=sum(tested*(-dose)*(1-lambda)**(dose-1)/(1-(1-lambda)**dose))
   else
      f=sum(tested*log(1-exp(-lambda*dose)))-log(alpha)
      deriv=sum(tested*dose*exp(-dose*lambda)/(1-exp(-dose*lambda)))
   end if
   step=f/deriv
   lambda=lambda-step
   if(-step<1e-6_dp)exit
end do
end function

function elda_one_group(response,dose,tested,observed,confidence,tol,maxit) result(out)
real(dp),intent(in)::response(:),dose(:),tested(:)
logical,intent(in),optional::observed
real(dp),intent(in),optional::confidence,tol
integer,intent(in),optional::maxit
type(elda_one_group_result_t)::out
real(dp)::conf,tolerance,sizev,nall,pmean,lambda,ph,step,d1,d2,alpha,z,se,fi
real(dp),allocatable::p(:),onemp(:),v(:),x(:),mueta(:),infoa(:),mubeta(:)
integer::iter,mx
logical::obs
obs=.false.
if(present(observed))obs=observed
conf=0.95_dp
if(present(confidence))conf=confidence
tolerance=1e-8_dp
if(present(tol))tolerance=tol
mx=100
if(present(maxit))mx=maxit
sizev=1-conf
allocate(p(size(response)))
if(all(response<1e-14_dp))then
   nall=sum(dose*tested)
   if(obs)then
   ph=1-sizev**(1/nall)
   else
   ph=-log(sizev)/nall
   end if
   out%ci_frequency=[ieee_value(ph,ieee_positive_inf),ieee_value(ph,ieee_positive_inf),1/ph]
   allocate(out%p(size(response)))
   out%p=0
   return
end if
if(all(response/tested>1-1e-14_dp))then
   ph=limdil_allpos(tested,dose,conf,obs)
   out%ci_frequency=[1/ph,1.0_dp,1.0_dp]
   allocate(out%p(size(response)))
   out%p=1
   return
end if
pmean=sum(response)/sum(tested)
lambda=-log(1-pmean)/maxval(dose)
do iter=1,mx
   p=-expm1_local(-lambda*dose)
   onemp=exp(-lambda*dose)
   d1=sum(tested*dose*(response/tested-p)/p)/size(p)
   d2=-sum(tested*(response/tested)*dose*dose*onemp/(p*p))/size(p)
   step=d1/d2
   lambda=lambda-step
   if(abs(step)<tolerance)exit
end do
out%iter=iter
out%lambda=lambda
alpha=log(lambda)
out%alpha=alpha
p=-expm1_local(-lambda*dose)
onemp=exp(-lambda*dose)
out%p=p
fi=sum(tested*dose*dose*onemp/p)*lambda**2
se=1/sqrt(fi)
z=qnorm((1-conf)/2,lower_tail=.false.)
out%ci_alpha=[alpha-z*se,alpha,alpha+z*se]
if(obs)then
out%ci_frequency=-1/expm1_local(-exp(out%ci_alpha))
else
out%ci_frequency=exp(-out%ci_alpha)
end if
out%deviance=binomial_deviance(response,tested,p)
v=p*onemp/tested
x=log(dose)
mueta=lambda*dose*onemp
infoa=mueta*mueta/v
ph=sum(x*infoa)/sum(infoa)
mubeta=(x-ph)*mueta
out%dloglik_logdose=sum(mubeta*(response/tested-p)/v)
out%fisher_logdose=sum(mubeta*mubeta/v)
out%z_score_logdose=out%dloglik_logdose/sqrt(out%fisher_logdose)
x=dose
ph=sum(x*infoa)/sum(infoa)
mubeta=(x-ph)*mueta
out%dloglik_dose=sum(mubeta*(response/tested-p)/v)
out%fisher_dose=sum(mubeta*mubeta/v)
out%z_score_dose=out%dloglik_dose/sqrt(out%fisher_dose)
end function

pure function binomial_deviance(y,n,p) result(dev)
real(dp),intent(in)::y(:),n(:),p(:)
real(dp)::dev,yy,nn,term
integer::i
dev=0
do i=1,size(y)
   yy=y(i)
   nn=n(i)
   term=0
   if(yy>0)term=term+yy*log((yy/nn)/p(i))
   if(yy<nn)term=term+(nn-yy)*log(((nn-yy)/nn)/(1-p(i)))
   dev=dev+2*term
end do
end function

function elda_fit(response,dose,tested,group,observed,confidence,test_unit_slope) result(out)
real(dp),intent(in)::response(:),dose(:),tested(:)
integer,intent(in)::group(:)
logical,intent(in),optional::observed,test_unit_slope
real(dp),intent(in),optional::confidence
type(elda_result_t)::out
type(elda_one_group_result_t)::fg,feq
integer,allocatable::levels(:),idx(:)
integer::ng,i,j,n
real(dp)::dev0,dl1,fi1,dl2,fi2,conf
logical::obs,ts
obs=.false.
if(present(observed))obs=observed
ts=.false.
if(present(test_unit_slope))ts=test_unit_slope
conf=0.95_dp
if(present(confidence))conf=confidence
call unique_int(group,levels)
ng=size(levels)
allocate(out%ci(ng,3))
dev0=0
dl1=0
fi1=0
dl2=0
fi2=0
n=size(group)
do i=1,ng
   idx=pack([(j,j=1,n)],group==levels(i))
   fg=elda_one_group(response(idx),dose(idx),tested(idx),obs,conf)
   out%ci(i,:)=max(fg%ci_frequency,1.0_dp)
   dev0=dev0+fg%deviance
   dl1=dl1+fg%dloglik_logdose
   fi1=fi1+fg%fisher_logdose
   dl2=dl2+fg%dloglik_dose
   fi2=fi2+fg%fisher_dose
end do
if(ng>1)then
   feq=elda_one_group(response,dose,tested,obs,conf)
   out%difference_chisq=max(feq%deviance-dev0,0.0_dp)
   out%difference_df=ng-1
   out%difference_p=1-pchisq(out%difference_chisq,real(ng-1,dp))
end if
if(fi1>1e-15_dp)then
   out%score_logdose_z=dl1/sqrt(fi1)
   out%score_logdose_p=erfc(abs(out%score_logdose_z)/sqrt(2.0_dp))
   out%score_dose_z=dl2/sqrt(fi2)
   out%score_dose_p=erfc(abs(out%score_dose_z)/sqrt(2.0_dp))
end if
if(ts)call elda_slope_fit(response,dose,tested,group,levels,dev0,out)
end function

subroutine elda_slope_fit(response,dose,tested,group,levels,dev0,out)
real(dp),intent(in)::response(:),dose(:),tested(:),dev0
integer,intent(in)::group(:),levels(:)
type(elda_result_t),intent(inout)::out
real(dp),allocatable::x(:,:),beta(:),oldbeta(:),mu(:),eta(:),dmu(:),w(:),z(:),fit(:),res(:),xtwx(:,:),ainv(:,:)
integer::n,ng,i,j,iter,rank,info
real(dp)::dev,se
n=size(response)
ng=size(levels)
allocate(x(n,ng+1))
x=0
do i=1,n
   do j=1,ng
   if(group(i)==levels(j))x(i,j)=1
   end do
   x(i,ng+1)=log(dose(i))
end do
allocate(beta(ng+1))
beta=0
beta(ng+1)=1
! initialize intercepts from group one-group alphas
 do j=1,ng
   block
      integer,allocatable::ix(:)
      type(elda_one_group_result_t)::fg
      ix=pack([(i,i=1,n)],group==levels(j))
      fg=elda_one_group(response(ix),dose(ix),tested(ix))
      beta(j)=fg%alpha
   end block
 end do
do iter=1,50
   oldbeta=beta
   eta=matmul(x,beta)
   mu=1-exp(-exp(eta))
   dmu=exp(eta-exp(eta))
   w=tested*dmu*dmu/max(mu*(1-mu),1e-14_dp)
   z=eta+(response/tested-mu)/max(dmu,1e-14_dp)
   call weighted_least_squares(x,z,w,beta,fit,res,rank,info)
   if(info/=0)exit
   if(maxval(abs(beta-oldbeta))<1e-10_dp)exit
end do
eta=matmul(x,beta)
mu=1-exp(-exp(eta))
dev=binomial_deviance(response,tested,mu)
out%slope_estimate=beta(ng+1)
out%slope_lr_z=sqrt(max(dev0-dev,0.0_dp))*sign(1.0_dp,beta(ng+1)-1)
out%slope_lr_p=1-pchisq(max(dev0-dev,0.0_dp),1.0_dp)
dmu=exp(eta-exp(eta))
w=tested*dmu*dmu/max(mu*(1-mu),1e-14_dp)
allocate(xtwx(ng+1,ng+1))
xtwx=0
do i=1,n
   do j=1,ng+1
   xtwx(j,:)=xtwx(j,:)+w(i)*x(i,j)*x(i,:)
   end do
end do
call symmetric_inverse(xtwx,ainv,info)
if(info==0)then
   se=sqrt(ainv(ng+1,ng+1))
   out%slope_se=se
   out%slope_wald_z=(out%slope_estimate-1)/se
   out%slope_wald_p=erfc(abs(out%slope_wald_z)/sqrt(2.0_dp))
end if
end subroutine

subroutine unique_int(x,u)
integer,intent(in)::x(:)
integer,allocatable,intent(out)::u(:)
integer::i,n
allocate(u(size(x)))
n=0
do i=1,size(x)
if(n==0.or..not.any(u(1:n)==x(i)))then
n=n+1
u(n)=x(i)
end if
end do
u=u(:n)
end subroutine


pure elemental function expm1_local(x) result(v)
real(dp),intent(in)::x
real(dp)::v
if(abs(x)<1.0e-5_dp)then
   v=x*(1.0_dp+x*(0.5_dp+x*(1.0_dp/6.0_dp+x*(1.0_dp/24.0_dp+x/120.0_dp))))
else
   v=exp(x)-1.0_dp
end if
end function expm1_local

end module statmod_elda
