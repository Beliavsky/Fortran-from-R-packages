! Additional computational kernels derived from statmod R/tweedie.R,
! R/glmscoretest.R and R/growthcurve.R.
! Upstream license: GPL-2 | GPL-3. See LICENSE and NOTICE.md.
module statmod_glm_misc
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use r_compat, only: dp
use statmod_linalg, only: least_squares
use statmod_utils, only: compare_two_growth_curves
implicit none
private
public :: tweedie_variance, tweedie_linkfun, tweedie_linkinv, tweedie_mu_eta
public :: tweedie_deviance_residual, glm_scoretest
public :: growth_comparison_t, compare_growth_curves, p_adjust_holm

type, public :: growth_comparison_t
   integer, allocatable :: group1(:), group2(:)
   real(dp), allocatable :: statistic(:), p_value(:), adjusted_p_value(:)
end type growth_comparison_t

contains

pure elemental function tweedie_variance(mu,var_power) result(v)
real(dp),intent(in)::mu,var_power
real(dp)::v
v=mu**var_power
end function tweedie_variance

pure elemental function tweedie_linkfun(mu,link_power) result(eta)
real(dp),intent(in)::mu,link_power
real(dp)::eta
if(link_power==0.0_dp)then
   eta=log(mu)
else
   eta=mu**link_power
end if
end function tweedie_linkfun

pure elemental function tweedie_linkinv(eta,link_power) result(mu)
real(dp),intent(in)::eta,link_power
real(dp)::mu
if(link_power==0.0_dp)then
   mu=max(exp(eta),epsilon(1.0_dp))
else
   mu=eta**(1.0_dp/link_power)
end if
end function tweedie_linkinv

pure elemental function tweedie_mu_eta(eta,link_power) result(v)
real(dp),intent(in)::eta,link_power
real(dp)::v
if(link_power==0.0_dp)then
   v=max(exp(eta),epsilon(1.0_dp))
else
   v=(1.0_dp/link_power)*eta**(1.0_dp/link_power-1.0_dp)
end if
end function tweedie_mu_eta

pure elemental function tweedie_deviance_residual(y,mu,wt,var_power) result(dev)
real(dp),intent(in)::y,mu,wt,var_power
real(dp)::dev,y1,theta,kappa
! This is the dev.resids component of statmod::tweedie().
y1=y+0.1_dp*merge(1.0_dp,0.0_dp,y==0.0_dp)
if(var_power==1.0_dp)then
   theta=log(y1/mu)
else
   theta=(y1**(1.0_dp-var_power)-mu**(1.0_dp-var_power))/(1.0_dp-var_power)
end if
if(var_power==2.0_dp)then
   kappa=log(y1/mu)
else
   kappa=(y**(2.0_dp-var_power)-mu**(2.0_dp-var_power))/(2.0_dp-var_power)
end if
dev=2.0_dp*wt*(y*theta-kappa)
end function tweedie_deviance_residual

function glm_scoretest(residuals,weights,x_existing,x_new,dispersion) result(score)
! Numerical kernel of glm.scoretest after extracting a fitted GLM's residuals,
! weights and existing design matrix. Columns of x_new are tested separately.
real(dp),intent(in)::residuals(:),weights(:),x_existing(:,:),x_new(:,:)
real(dp),intent(in),optional::dispersion
real(dp),allocatable::score(:)
real(dp),allocatable::xe(:,:),xn(:,:),r(:),w(:),ws(:),beta(:),fit(:),rr(:),xres(:,:)
real(dp)::disp,den
integer::i,j,nkeep,rank,info,df
logical,allocatable::keep(:)
keep=weights>0.0_dp
nkeep=count(keep)
allocate(score(size(x_new,2)))
score=0.0_dp
if(nkeep==0)return
r=pack(residuals,keep)
w=pack(weights,keep)
ws=sqrt(w)
allocate(xe(nkeep,size(x_existing,2)),xn(nkeep,size(x_new,2)))
do j=1,size(x_existing,2)
xe(:,j)=pack(x_existing(:,j),keep)*ws
end do
do j=1,size(x_new,2)
xn(:,j)=pack(x_new(:,j),keep)*ws
end do
if(present(dispersion))then
   disp=dispersion
else
   df=nkeep-size(x_existing,2)
   if(df<=0)return
   disp=sum(w*r*r)/real(df,dp)
end if
allocate(xres(nkeep,size(x_new,2)))
do j=1,size(x_new,2)
   if(size(x_existing,2)>0)then
      call least_squares(xe,xn(:,j),beta,fit,rr,rank,info)
      if(info/=0)return
      xres(:,j)=rr
   else
      xres(:,j)=xn(:,j)
   end if
   den=sqrt(sum(xres(:,j)*xres(:,j))*disp)
   if(den>0.0_dp)score(j)=sum(xres(:,j)*(ws*r))/den
end do
end function glm_scoretest

function compare_growth_curves(group,y,nsim,n0) result(out)
! Pairwise computational kernel of compareGrowthCurves for integer group labels.
integer,intent(in)::group(:)
real(dp),intent(in)::y(:,:)
integer,intent(in),optional::nsim
real(dp),intent(in),optional::n0
type(growth_comparison_t)::out
integer,allocatable::lev(:),sel(:)
integer::i,j,k,nlev,np,ns
real(dp)::base,stat,pv
call unique_int(group,lev)
nlev=size(lev)
np=nlev*(nlev-1)/2
allocate(out%group1(np),out%group2(np),out%statistic(np),out%p_value(np),out%adjusted_p_value(np))
ns=100
if(present(nsim))ns=nsim
base=0.5_dp
if(present(n0))base=n0
k=0
do i=1,nlev-1
   do j=i+1,nlev
      k=k+1
      sel=pack([(np,np=1,size(group))],(group==lev(i)).or.(group==lev(j)))
      call compare_two_growth_curves(group(sel),y(sel,:),ns,stat,pv,base)
      out%group1(k)=lev(i)
      out%group2(k)=lev(j)
      out%statistic(k)=stat
      out%p_value(k)=pv
   end do
end do
out%adjusted_p_value=p_adjust_holm(out%p_value)
end function compare_growth_curves

function p_adjust_holm(p) result(adj)
real(dp),intent(in)::p(:)
real(dp)::adj(size(p)),ps(size(p)),running
integer::ord(size(p)),i,j,n,tmpi
n=size(p)
ps=p
ord=[(i,i=1,n)]
! Stable insertion sort of indices by p value.
do i=2,n
   tmpi=ord(i)
   j=i-1
   do while(j>=1)
      if(p(ord(j))<=p(tmpi))exit
      ord(j+1)=ord(j)
      j=j-1
   end do
   ord(j+1)=tmpi
end do
running=0.0_dp
adj=1.0_dp
do i=1,n
   running=max(running,real(n-i+1,dp)*p(ord(i)))
   adj(ord(i))=min(1.0_dp,running)
end do
end function p_adjust_holm

subroutine unique_int(x,u)
integer,intent(in)::x(:)
integer,allocatable,intent(out)::u(:)
integer,allocatable::tmp(:)
integer::i,n
allocate(tmp(size(x)))
n=0
do i=1,size(x)
   if(n==0)then
      n=1
      tmp(n)=x(i)
   else if(.not.any(tmp(1:n)==x(i)))then
      n=n+1
      tmp(n)=x(i)
   end if
end do
allocate(u(n))
u=tmp(:n)
end subroutine unique_int

end module statmod_glm_misc
