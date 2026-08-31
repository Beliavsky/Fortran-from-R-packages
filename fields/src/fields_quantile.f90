! GPL-2.0-or-later. Quantile smoothing workflows translated from QSreg/QTps.
module fields_quantile
use fields_kinds, only: dp
use fields_spline1d, only: spline_fit, smoothing_spline, smoothing_spline_gcv, spline_predict
use fields_tps, only: tps_fit_type, tps_fit, tps_fit_gcv, tps_predict
use fields_stats, only: quantile_type7
implicit none
private

public :: quantile_spline_fit, quantile_tps_fit
public :: quantile_smoothing_spline, quantile_thin_plate_spline
public :: quantile_spline_predict, quantile_tps_predict
public :: qsreg_psi, qsreg_sigma

type :: quantile_spline_fit
   type(spline_fit) :: fit
   logical :: converged=.false.
   integer :: iterations=0
   real(dp) :: psi_scale=0.0_dp
   real(dp) :: alpha=0.5_dp
   real(dp) :: cv_pseudo=0.0_dp
   real(dp), allocatable :: y_raw(:), pseudo_y(:), convergence(:)
end type quantile_spline_fit

type :: quantile_tps_fit
   type(tps_fit_type) :: fit
   logical :: converged=.false.
   integer :: iterations=0
   real(dp) :: psi_scale=0.0_dp
   real(dp) :: alpha=0.5_dp
   real(dp) :: cv_pseudo=0.0_dp
   real(dp), allocatable :: y_raw(:), pseudo_y(:), convergence(:)
end type quantile_tps_fit

contains

pure elemental real(dp) function qsreg_psi(r,alpha,c) result(v)
real(dp),intent(in)::r
real(dp),intent(in),optional::alpha,c
real(dp)::a,cc,t
a=0.5_dp;if(present(alpha))a=alpha
cc=1.0_dp;if(present(c))cc=c
if(cc<=0.0_dp .or. a<=0.0_dp .or. a>=1.0_dp) then
   v=0.0_dp; return
end if
if(r<0.0_dp) then
   t=2.0_dp*(1.0_dp-a)*r/cc
else
   t=2.0_dp*a*r/cc
end if
v=min(2.0_dp*a,max(-2.0_dp*(1.0_dp-a),t))
end function qsreg_psi

pure elemental real(dp) function qsreg_sigma(r,alpha,c) result(v)
real(dp),intent(in)::r
real(dp),intent(in),optional::alpha,c
real(dp)::a,cc
a=0.5_dp;if(present(alpha))a=alpha
cc=1.0_dp;if(present(c))cc=c
if(cc<=0.0_dp .or. a<=0.0_dp .or. a>=1.0_dp) then
   v=huge(1.0_dp); return
end if
if(r<0.0_dp) then
   v=(1.0_dp-a)*r*r/cc
else
   v=a*r*r/cc
end if
if(r>cc) v=2.0_dp*a*r-a*cc
if(r< -cc) v=-2.0_dp*(1.0_dp-a)*r-(1.0_dp-a)*cc
end function qsreg_sigma

function quantile_smoothing_spline(x,y,lambda,f_start,psi_scale,cost,alpha,maxiter,tolerance) result(out)
real(dp),intent(in)::x(:),y(:)
real(dp),intent(in),optional::lambda,f_start(:),psi_scale,cost,alpha,tolerance
integer,intent(in),optional::maxiter
type(quantile_spline_fit)::out
real(dp),allocatable::f(:),fnew(:),yp(:),fcv(:)
real(dp)::ps,c,a,tol,denom,change
integer::n,it,nit
logical::fixed_lambda
n=size(y);if(size(x)/=n .or. n<3)error stop 'quantile_smoothing_spline: invalid x/y'
a=0.5_dp;if(present(alpha))a=alpha
c=1.0_dp;if(present(cost))c=cost
tol=1.0e-3_dp;if(present(tolerance))tol=tolerance
nit=100;if(present(maxiter))nit=maxiter
fixed_lambda=present(lambda)
allocate(f(n),fnew(n),yp(n),out%convergence(nit))
if(present(f_start))then
   if(size(f_start)/=n)error stop 'quantile_smoothing_spline: f_start mismatch'
   f=f_start
else
   f=quantile_type7(y,0.5_dp)
end if
ps=mad_scale(y-f);if(present(psi_scale))ps=psi_scale
if(ps<=0.0_dp)ps=max(1.0e-12_dp,sqrt(sum((y-sum(y)/real(n,dp))**2)/real(max(1,n-1),dp))*0.05_dp)
out%convergence=huge(1.0_dp)
do it=1,nit
   yp=f+c*ps*qsreg_psi((y-f)/ps,c=c,alpha=a)
   if(fixed_lambda)then
      out%fit=smoothing_spline(x,yp,lambda,find_diag=.true.)
   else
      out%fit=smoothing_spline_gcv(x,yp)
   end if
   fnew=out%fit%fitted
   denom=max(1.0e-14_dp,sum(abs(f))/real(n,dp))
   change=sum(abs(fnew-f))/real(n,dp)/denom
   out%convergence(it)=change
   f=fnew
   if(change<=tol)then;out%converged=.true.;exit;end if
end do
out%iterations=it
! Reproduce QSreg's final pseudo-data/refit and pseudo-data LOO criterion.
yp=f+c*ps*qsreg_psi((y-f)/ps,c=c,alpha=a)
if(fixed_lambda)then
   out%fit=smoothing_spline(x,yp,lambda,find_diag=.true.)
else
   out%fit=smoothing_spline_gcv(x,yp)
end if
allocate(fcv(n))
where(abs(1.0_dp-out%fit%diag_a)>1.0e-12_dp)
   fcv=out%fit%fitted/(1.0_dp-out%fit%diag_a)-out%fit%diag_a*yp/(1.0_dp-out%fit%diag_a)
elsewhere
   fcv=out%fit%fitted
end where
out%cv_pseudo=sum(qsreg_sigma(y-fcv,alpha=a,c=ps))/real(n,dp)
out%psi_scale=ps;out%alpha=a;out%y_raw=y;out%pseudo_y=yp
end function quantile_smoothing_spline

function quantile_thin_plate_spline(x,y,lambda,m,f_start,psi_scale,cost,alpha,maxiter,tolerance) result(out)
real(dp),intent(in)::x(:,:),y(:)
real(dp),intent(in),optional::lambda,f_start(:),psi_scale,cost,alpha,tolerance
integer,intent(in),optional::m,maxiter
type(quantile_tps_fit)::out
real(dp),allocatable::f(:),fnew(:),yp(:),fcv(:)
real(dp)::ps,c,a,tol,denom,change
integer::n,it,nit,mm
logical::fixed_lambda
n=size(y);if(size(x,1)/=n)error stop 'quantile_thin_plate_spline: x/y mismatch'
a=0.5_dp;if(present(alpha))a=alpha
c=1.0_dp;if(present(cost))c=cost
tol=1.0e-3_dp;if(present(tolerance))tol=tolerance
nit=100;if(present(maxiter))nit=maxiter
mm=max(2,int(ceiling(0.5_dp*real(size(x,2),dp)+0.1_dp)));if(present(m))mm=m
fixed_lambda=present(lambda)
allocate(f(n),fnew(n),yp(n),out%convergence(nit))
if(present(f_start))then
   if(size(f_start)/=n)error stop 'quantile_thin_plate_spline: f_start mismatch'
   f=f_start
else
   f=quantile_type7(y,0.5_dp)
end if
ps=mad_scale(y-f);if(present(psi_scale))ps=psi_scale
if(ps<=0.0_dp)ps=max(1.0e-12_dp,sqrt(sum((y-sum(y)/real(n,dp))**2)/real(max(1,n-1),dp))*0.05_dp)
out%convergence=huge(1.0_dp)
do it=1,nit
   yp=f+c*ps*qsreg_psi((y-f)/ps,c=c,alpha=a)
   if(fixed_lambda)then
      out%fit=tps_fit(x,yp,lambda,mm)
   else
      out%fit=tps_fit_gcv(x,yp,mm)
   end if
   if(out%fit%info/=0)exit
   fnew=out%fit%fitted
   denom=max(1.0e-14_dp,sum(abs(f))/real(n,dp))
   change=sum(abs(fnew-f))/real(n,dp)/denom
   out%convergence(it)=change
   f=fnew
   if(change<=tol)then;out%converged=.true.;exit;end if
end do
out%iterations=it
yp=f+c*ps*qsreg_psi((y-f)/ps,c=c,alpha=a)
if(fixed_lambda)then;out%fit=tps_fit(x,yp,lambda,mm);else;out%fit=tps_fit_gcv(x,yp,mm);end if
allocate(fcv(n))
where(abs(1.0_dp-out%fit%diag_a)>1.0e-12_dp)
   fcv=out%fit%fitted/(1.0_dp-out%fit%diag_a)-out%fit%diag_a*yp/(1.0_dp-out%fit%diag_a)
elsewhere
   fcv=out%fit%fitted
end where
out%cv_pseudo=sum(qsreg_sigma(y-fcv,alpha=a,c=ps))/real(n,dp)
out%psi_scale=ps;out%alpha=a;out%y_raw=y;out%pseudo_y=yp
end function quantile_thin_plate_spline

function quantile_spline_predict(out,xnew,derivative) result(v)
type(quantile_spline_fit),intent(in)::out
real(dp),intent(in)::xnew(:)
integer,intent(in),optional::derivative
real(dp),allocatable::v(:)
if(present(derivative))then
   v=spline_predict(out%fit,xnew,derivative)
else
   v=spline_predict(out%fit,xnew)
end if
end function quantile_spline_predict

function quantile_tps_predict(out,xnew) result(v)
type(quantile_tps_fit),intent(in)::out
real(dp),intent(in)::xnew(:,:)
real(dp),allocatable::v(:)
v=tps_predict(out%fit,xnew)
end function quantile_tps_predict

real(dp) function mad_scale(x) result(v)
real(dp),intent(in)::x(:)
real(dp),allocatable::a(:)
real(dp)::med
med=quantile_type7(x,0.5_dp);allocate(a(size(x)));a=abs(x-med)
v=1.482602218505602_dp*quantile_type7(a,0.5_dp)
end function mad_scale

end module fields_quantile
