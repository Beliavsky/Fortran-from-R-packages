! GPL-2.0-or-later. Likelihood/profile fitting corresponding to spatialProcess/mKrigMLE*.
module fields_spatial_process
use fields_kinds, only: dp
use fields_covariance, only: stationary_covariance
use fields_kriging, only: krig_fit, krig_fit_covariance
implicit none
private
public :: spatial_process_fit, spatial_profile_grid_result
public :: spatial_profile_grid, fit_spatial_process

type :: spatial_process_fit
   character(len=32) :: model='exponential'
   real(dp) :: a_range=1.0_dp
   real(dp) :: smoothness=0.5_dp
   real(dp) :: power=1.0_dp
   real(dp) :: objective=huge(1.0_dp)
   type(krig_fit) :: fit
end type

type :: spatial_profile_grid_result
   real(dp), allocatable :: a_range(:),lambda(:),objective(:,:)
   integer :: best_range=0,best_lambda=0
end type

contains

function spatial_profile_grid(x,y,a_ranges,lambdas,model,smoothness,power,t,weights,criterion) result(out)
real(dp),intent(in)::x(:,:),y(:),a_ranges(:),lambdas(:)
character(len=*),intent(in),optional::model,criterion
real(dp),intent(in),optional::smoothness,power,t(:,:),weights(:)
type(spatial_profile_grid_result)::out
real(dp),allocatable::k(:,:)
type(krig_fit)::fit
real(dp)::nu,pp,best
character(len=32)::mod,crit
integer::i,j
nu=0.5_dp;if(present(smoothness))nu=smoothness
pp=1.0_dp;if(present(power))pp=power
mod='exponential';if(present(model))mod=lower(trim(model))
crit='reml';if(present(criterion))crit=lower(trim(criterion))
if(any(a_ranges<=0.0_dp).or.any(lambdas<0.0_dp))error stop 'spatial_profile_grid: invalid grid'
out%a_range=a_ranges;out%lambda=lambdas;allocate(out%objective(size(a_ranges),size(lambdas)))
best=huge(1.0_dp)
do i=1,size(a_ranges)
   k=stationary_covariance(x,x,model=mod,a_range=a_ranges(i),smoothness=nu,p=pp)
   do j=1,size(lambdas)
      fit=krig_fit_covariance(y,k,lambdas(j),t,weights)
      if(fit%info/=0)then
         out%objective(i,j)=huge(1.0_dp)
      else
         select case(trim(crit))
         case('gcv');out%objective(i,j)=fit%gcv
         case('mle','ml');out%objective(i,j)=fit%profile_mle
         case default;out%objective(i,j)=fit%profile_reml
         end select
      end if
      if(out%objective(i,j)<best)then;best=out%objective(i,j);out%best_range=i;out%best_lambda=j;end if
   end do
end do
end function spatial_profile_grid

function fit_spatial_process(x,y,model,a_range_start,lambda_start,smoothness_start,power,t,weights,criterion, &
                            range_bounds,lambda_bounds,smoothness_bounds,optimize_smoothness,tol,maxiter) result(out)
real(dp),intent(in)::x(:,:),y(:)
character(len=*),intent(in),optional::model,criterion
real(dp),intent(in),optional::a_range_start,lambda_start,smoothness_start,power,t(:,:),weights(:)
real(dp),intent(in),optional::range_bounds(2),lambda_bounds(2),smoothness_bounds(2),tol
logical,intent(in),optional::optimize_smoothness
integer,intent(in),optional::maxiter
type(spatial_process_fit)::out
real(dp)::lr,ll,lnu,rb(2),lb(2),sb(2),eps,old,obj
real(dp)::pp
integer::it,nit
logical::opts
character(len=32)::mod,crit
mod='exponential';if(present(model))mod=lower(trim(model))
crit='reml';if(present(criterion))crit=lower(trim(criterion))
pp=1.0_dp;if(present(power))pp=power
lr=log(1.0_dp);if(present(a_range_start))lr=log(a_range_start)
ll=log(0.1_dp);if(present(lambda_start))ll=log(max(lambda_start,1.0e-14_dp))
lnu=log(0.5_dp);if(present(smoothness_start))lnu=log(smoothness_start)
rb=[1.0e-4_dp,1.0e4_dp];if(present(range_bounds))rb=range_bounds
lb=[1.0e-10_dp,1.0e10_dp];if(present(lambda_bounds))lb=lambda_bounds
sb=[0.05_dp,10.0_dp];if(present(smoothness_bounds))sb=smoothness_bounds
if(any(rb<=0.0_dp).or.any(lb<=0.0_dp).or.any(sb<=0.0_dp))error stop 'fit_spatial_process: bounds must be positive'
eps=1.0e-5_dp;if(present(tol))eps=tol;nit=30;if(present(maxiter))nit=maxiter
opts=.false.;if(present(optimize_smoothness))opts=optimize_smoothness
obj=spatial_objective(lr,ll,lnu,x,y,mod,pp,t,weights,crit)
do it=1,nit
   old=obj
   lr=golden_one_parameter(1,lr,ll,lnu,log(rb(1)),log(rb(2)),x,y,mod,pp,t,weights,crit,eps)
   ll=golden_one_parameter(2,lr,ll,lnu,log(lb(1)),log(lb(2)),x,y,mod,pp,t,weights,crit,eps)
   if(opts.and.trim(mod)=='matern') lnu=golden_one_parameter(3,lr,ll,lnu,log(sb(1)),log(sb(2)),x,y,mod,pp,t,weights,crit,eps)
   obj=spatial_objective(lr,ll,lnu,x,y,mod,pp,t,weights,crit)
   if(abs(old-obj)<=eps*(1.0_dp+abs(obj)))exit
end do
out%model=mod;out%a_range=exp(lr);out%smoothness=exp(lnu);out%power=pp;out%objective=obj
out%fit=krig_fit_covariance(y,stationary_covariance(x,x,model=mod,a_range=out%a_range,smoothness=out%smoothness,p=pp), &
                            exp(ll),t,weights)
end function fit_spatial_process

real(dp) function golden_one_parameter(which,lr,ll,lnu,a,b,x,y,model,p,t,weights,criterion,tol) result(best)
integer,intent(in)::which
real(dp),intent(in)::lr,ll,lnu,a,b,x(:,:),y(:),p,tol
real(dp),intent(in),optional::t(:,:),weights(:)
character(len=*),intent(in)::model,criterion
real(dp)::left,right,c,d,fc,fd,gr,r1,l1,n1
integer::it
left=a;right=b;gr=(sqrt(5.0_dp)-1.0_dp)/2.0_dp
c=right-gr*(right-left);d=left+gr*(right-left)
r1=lr;l1=ll;n1=lnu;call set_one(which,c,r1,l1,n1);fc=spatial_objective(r1,l1,n1,x,y,model,p,t,weights,criterion)
r1=lr;l1=ll;n1=lnu;call set_one(which,d,r1,l1,n1);fd=spatial_objective(r1,l1,n1,x,y,model,p,t,weights,criterion)
do it=1,100
   if(abs(right-left)<=tol*(1.0_dp+abs(left)+abs(right)))exit
   if(fc<=fd)then
      right=d;d=c;fd=fc;c=right-gr*(right-left);r1=lr;l1=ll;n1=lnu;call set_one(which,c,r1,l1,n1)
      fc=spatial_objective(r1,l1,n1,x,y,model,p,t,weights,criterion)
   else
      left=c;c=d;fc=fd;d=left+gr*(right-left);r1=lr;l1=ll;n1=lnu;call set_one(which,d,r1,l1,n1)
      fd=spatial_objective(r1,l1,n1,x,y,model,p,t,weights,criterion)
   end if
end do
if(fc<=fd)then;best=c;else;best=d;end if
end function golden_one_parameter

subroutine set_one(which,v,lr,ll,lnu)
integer,intent(in)::which
real(dp),intent(in)::v
real(dp),intent(inout)::lr,ll,lnu
select case(which);case(1);lr=v;case(2);ll=v;case(3);lnu=v;end select
end subroutine set_one

real(dp) function spatial_objective(lr,ll,lnu,x,y,model,p,t,weights,criterion) result(v)
real(dp),intent(in)::lr,ll,lnu,x(:,:),y(:),p
real(dp),intent(in),optional::t(:,:),weights(:)
character(len=*),intent(in)::model,criterion
real(dp),allocatable::k(:,:)
type(krig_fit)::fit
k=stationary_covariance(x,x,model=model,a_range=exp(lr),smoothness=exp(lnu),p=p)
fit=krig_fit_covariance(y,k,exp(ll),t,weights)
if(fit%info/=0)then;v=huge(1.0_dp);return;end if
select case(trim(criterion));case('gcv');v=fit%gcv;case('mle','ml');v=fit%profile_mle;case default;v=fit%profile_reml;end select
end function spatial_objective

pure function lower(s) result(t)
character(len=*),intent(in)::s
character(len=len(s))::t
integer::i,c
t=s;do i=1,len(s);c=iachar(t(i:i));if(c>=65.and.c<=90)t(i:i)=achar(c+32);end do
end function lower

end module fields_spatial_process
