! GPL-2.0-or-later. Modern wrappers around the fields cubic spline kernels.
module fields_spline1d
use fields_kinds, only: dp
use fields_native, only: css, rcss
implicit none
private

public :: spline_fit, robust_spline_fit
public :: smoothing_spline, smoothing_spline_gcv, spline_predict
public :: robust_smoothing_spline, robust_spline_predict
public :: spline_trace, spline_df_to_lambda

type :: spline_fit
   real(dp) :: lambda = 1.0_dp
   real(dp) :: trace = 0.0_dp
   real(dp) :: gcv = 0.0_dp
   integer :: ierr = 0
   real(dp), allocatable :: x(:)
   real(dp), allocatable :: y(:)
   real(dp), allocatable :: weights(:)
   real(dp), allocatable :: fitted(:)
   real(dp), allocatable :: diag_a(:)
end type spline_fit

type :: robust_spline_fit
   real(dp) :: lambda = 1.0_dp
   real(dp) :: trace = 0.0_dp
   real(dp) :: cv = 0.0_dp
   integer :: iterations = 0
   real(dp) :: convergence_test = 0.0_dp
   integer :: ierr = 0
   real(dp), allocatable :: x(:)
   real(dp), allocatable :: y(:)
   real(dp), allocatable :: robust_weights(:)
   real(dp), allocatable :: fitted(:)
   real(dp), allocatable :: diag_a(:)
end type robust_spline_fit

contains

function smoothing_spline(x,y,lambda,weights,find_diag) result(fit)
real(dp), intent(in) :: x(:),y(:),lambda
real(dp), intent(in), optional :: weights(:)
logical, intent(in), optional :: find_diag
type(spline_fit) :: fit
real(dp), allocatable :: xw(:),yw(:),sdw(:),sy(:),diag(:),xg(:),yg(:)
integer :: n,job(3),ierr
real(dp) :: trace,gcv,h
logical :: fd
n=size(x)
if(size(y)/=n .or. n<2) error stop 'smoothing_spline: x/y length mismatch or too few points'
if(lambda<0.0_dp) error stop 'smoothing_spline: lambda must be nonnegative'
allocate(xw(n),yw(n),sdw(n),sy(n),diag(n),xg(1),yg(1))
xw=x; yw=y; sdw=1.0_dp
if(present(weights)) then
   if(size(weights)/=n .or. any(weights<=0.0_dp)) error stop 'smoothing_spline: invalid weights'
   sdw=1.0_dp/sqrt(weights)
end if
fd=.true.; if(present(find_diag)) fd=find_diag
job=[merge(3,1,fd),0,0]
diag=0.0_dp
if(job(1)==3) then
   diag(1)=1.0_dp
   if(n>=2) diag(2)=0.0_dp
end if
xg=0.0_dp; yg=0.0_dp; trace=0.0_dp; gcv=0.0_dp
if(lambda==0.0_dp) then
   h=-1001.0_dp
else
   h=log(lambda)
end if
ierr=0
call css(h,n,xw,yw,sdw,sy,trace,diag,gcv,1,xg,yg,job,0,ierr)
fit%lambda=lambda; fit%trace=trace; fit%gcv=gcv; fit%ierr=ierr
fit%x=x; fit%y=y; allocate(fit%weights(n)); fit%weights=1.0_dp
if(present(weights)) fit%weights=weights
fit%fitted=sy
if(fd) fit%diag_a=diag
end function smoothing_spline

function smoothing_spline_gcv(x,y,weights,lambda_min,lambda_max,cost,offset,tol,maxiter) result(fit)
real(dp), intent(in) :: x(:),y(:)
real(dp), intent(in), optional :: weights(:),lambda_min,lambda_max,cost,offset,tol
integer, intent(in), optional :: maxiter
type(spline_fit) :: fit
real(dp) :: a,b,c,d,fc,fd,gr,eps,cost0,off0,lamlo,lamhi
integer :: it,nit
lamlo=1.0e-12_dp; if(present(lambda_min)) lamlo=lambda_min
lamhi=1.0e12_dp; if(present(lambda_max)) lamhi=lambda_max
if(lamlo<=0.0_dp .or. lamhi<=lamlo) error stop 'smoothing_spline_gcv: invalid lambda interval'
cost0=1.0_dp; if(present(cost)) cost0=cost
off0=0.0_dp; if(present(offset)) off0=offset
eps=1.0e-6_dp; if(present(tol)) eps=tol
nit=100; if(present(maxiter)) nit=maxiter
a=log(lamlo); b=log(lamhi); gr=(sqrt(5.0_dp)-1.0_dp)/2.0_dp
c=b-gr*(b-a); d=a+gr*(b-a)
fc=spline_gcv_at_loglambda(c,x,y,weights,cost0,off0)
fd=spline_gcv_at_loglambda(d,x,y,weights,cost0,off0)
do it=1,nit
   if(abs(b-a)<=eps*(1.0_dp+abs(a)+abs(b))) exit
   if(fc<=fd) then
      b=d; d=c; fd=fc; c=b-gr*(b-a)
      fc=spline_gcv_at_loglambda(c,x,y,weights,cost0,off0)
   else
      a=c; c=d; fc=fd; d=a+gr*(b-a)
      fd=spline_gcv_at_loglambda(d,x,y,weights,cost0,off0)
   end if
end do
if(fc<=fd) then
   fit=smoothing_spline(x,y,exp(c),weights,find_diag=.true.)
else
   fit=smoothing_spline(x,y,exp(d),weights,find_diag=.true.)
end if
end function smoothing_spline_gcv

real(dp) function spline_gcv_at_loglambda(h,x,y,weights,cost,offset) result(gcv)
real(dp), intent(in) :: h,x(:),y(:),cost,offset
real(dp), intent(in), optional :: weights(:)
real(dp), allocatable :: xw(:),yw(:),sdw(:),sy(:),diag(:),xg(:),yg(:)
real(dp) :: trace
integer :: n,job(3),ierr
n=size(x); allocate(xw(n),yw(n),sdw(n),sy(n),diag(n),xg(1),yg(1))
xw=x; yw=y; sdw=1.0_dp
if(present(weights)) then
   if(size(weights)/=n .or. any(weights<=0.0_dp)) error stop 'spline_gcv: invalid weights'
   sdw=1.0_dp/sqrt(weights)
end if
diag=0.0_dp; diag(1)=cost; if(n>=2) diag(2)=offset
trace=0.0_dp; gcv=0.0_dp; xg=0.0_dp; yg=0.0_dp; job=[3,0,0]
ierr=0
call css(h,n,xw,yw,sdw,sy,trace,diag,gcv,1,xg,yg,job,0,ierr)
if(ierr/=0 .or. .not.(gcv>=0.0_dp)) gcv=huge(1.0_dp)
end function spline_gcv_at_loglambda

function spline_predict(fit,xnew,derivative) result(ynew)
type(spline_fit), intent(in) :: fit
real(dp), intent(in) :: xnew(:)
integer, intent(in), optional :: derivative
real(dp), allocatable :: ynew(:)
real(dp), allocatable :: xw(:),yw(:),sdw(:),sy(:),diag(:),xg(:)
real(dp) :: trace,gcv,h
integer :: n,job(3),ierr,ider
n=size(fit%x); if(n<2) error stop 'spline_predict: invalid fit'
ider=0; if(present(derivative)) ider=derivative
if(ider<0 .or. ider>2) error stop 'spline_predict: derivative must be 0, 1, or 2'
allocate(xw(n),yw(n),sdw(n),sy(n),diag(n),xg(size(xnew)),ynew(size(xnew)))
xw=fit%x; yw=fit%y; sdw=1.0_dp/sqrt(fit%weights); diag=0.0_dp; xg=xnew
trace=0.0_dp; gcv=0.0_dp; job=[0,3,0]
if(fit%lambda==0.0_dp) then; h=-1001.0_dp; else; h=log(fit%lambda); end if
ierr=0
call css(h,n,xw,yw,sdw,sy,trace,diag,gcv,size(xnew),xg,ynew,job,ider,ierr)
if(ierr/=0) error stop 'spline_predict: css failed'
end function spline_predict

real(dp) function spline_trace(x,weights,lambda) result(tr)
real(dp), intent(in) :: x(:),lambda
real(dp), intent(in), optional :: weights(:)
real(dp), allocatable :: y(:),xw(:),sdw(:),sy(:),diag(:),xg(:),yg(:)
real(dp) :: gcv,h
integer :: n,job(3),ierr
n=size(x); allocate(y(n),xw(n),sdw(n),sy(n),diag(n),xg(1),yg(1))
y=0.0_dp; xw=x; sdw=1.0_dp
if(present(weights)) sdw=1.0_dp/sqrt(weights)
diag=0.0_dp; tr=0.0_dp; gcv=0.0_dp; xg=0.0_dp; yg=0.0_dp; job=[1,0,0]
if(lambda==0.0_dp) then; h=-1001.0_dp; else; h=log(lambda); end if
ierr=0
call css(h,n,xw,y,sdw,sy,tr,diag,gcv,1,xg,yg,job,0,ierr)
if(ierr/=0) tr=huge(1.0_dp)
end function spline_trace

real(dp) function spline_df_to_lambda(df,x,weights,tol) result(lambda)
real(dp), intent(in) :: df,x(:)
real(dp), intent(in), optional :: weights(:),tol
real(dp) :: lo,hi,mid,t,eps
integer :: it,n
n=size(x); if(df<2.0_dp .or. df>real(n,dp)) error stop 'spline_df_to_lambda: df outside feasible range'
eps=1.0e-6_dp; if(present(tol)) eps=tol
lo=-30.0_dp; hi=30.0_dp
do it=1,100
   t=spline_trace(x,weights,exp(lo)); if(t>=df) exit
   lo=lo-5.0_dp
end do
do it=1,100
   t=spline_trace(x,weights,exp(hi)); if(t<=df) exit
   hi=hi+5.0_dp
end do
do it=1,120
   mid=0.5_dp*(lo+hi); t=spline_trace(x,weights,exp(mid))
   if(abs(t-df)<=eps*max(1.0_dp,df)) exit
   if(t>df) then; lo=mid; else; hi=mid; end if
end do
lambda=exp(mid)
end function spline_df_to_lambda

function robust_smoothing_spline(x,y,lambda,c_scale,alpha,maxiter,tol,cost,offset) result(fit)
real(dp), intent(in) :: x(:),y(:),lambda
real(dp), intent(in), optional :: c_scale,alpha,tol,cost,offset
integer, intent(in), optional :: maxiter
type(robust_spline_fit) :: fit
real(dp), allocatable :: xw(:),yw(:),wt(:),sy(:),diag(:),xg(:),yg(:)
real(dp) :: h,trace,cv,din(10),dout(10)
integer :: n,job(3),ideriv,ierr
n=size(x); if(size(y)/=n .or. n<2) error stop 'robust_smoothing_spline: invalid input'
allocate(xw(n),yw(n),wt(n),sy(n),diag(n),xg(1),yg(1))
xw=x; yw=y; wt=1.0_dp; sy=0.0_dp; diag=0.0_dp; xg=0.0_dp; yg=0.0_dp
din=0.0_dp; din(1:6)=[1.0_dp,0.0_dp,50.0_dp,1.0e-5_dp,1.345_dp,1.0_dp]
if(present(cost)) din(1)=cost
if(present(offset)) din(2)=offset
if(present(maxiter)) din(3)=real(maxiter,dp)
if(present(tol)) din(4)=tol
if(present(c_scale)) din(5)=c_scale
if(present(alpha)) din(6)=alpha
if(lambda==0.0_dp) then; h=-1001.0_dp; else; h=log(lambda); end if
trace=0.0_dp; cv=0.0_dp; dout=0.0_dp; job=[3,0,0]; ideriv=0
ierr=0
call rcss(h,n,xw,yw,wt,sy,trace,diag,cv,1,xg,yg,job,ideriv,din,dout,ierr)
fit%lambda=lambda; fit%trace=trace; fit%cv=cv; fit%ierr=ierr
fit%iterations=nint(dout(1)); fit%convergence_test=dout(2)
fit%x=x; fit%y=y; fit%robust_weights=wt; fit%fitted=sy; fit%diag_a=diag
end function robust_smoothing_spline

function robust_spline_predict(fit,xnew,derivative,c_scale,alpha) result(ynew)
type(robust_spline_fit), intent(in) :: fit
real(dp), intent(in) :: xnew(:)
integer, intent(in), optional :: derivative
real(dp), intent(in), optional :: c_scale,alpha
real(dp), allocatable :: ynew(:)
real(dp), allocatable :: xw(:),yw(:),wt(:),sy(:),diag(:),xg(:)
real(dp) :: h,trace,cv,din(10),dout(10)
integer :: n,job(3),ider,ierr
n=size(fit%x); allocate(xw(n),yw(n),wt(n),sy(n),diag(n),xg(size(xnew)),ynew(size(xnew)))
xw=fit%x; yw=fit%y; wt=fit%robust_weights; sy=0.0_dp; diag=0.0_dp; xg=xnew
din=0.0_dp; din(1:6)=[1.0_dp,0.0_dp,1.0_dp,1.0e-12_dp,1.345_dp,1.0_dp]
if(present(c_scale)) din(5)=c_scale
if(present(alpha)) din(6)=alpha
if(fit%lambda==0.0_dp) then; h=-1001.0_dp; else; h=log(fit%lambda); end if
ider=0; if(present(derivative)) ider=derivative
trace=0.0_dp; cv=0.0_dp; dout=0.0_dp; job=[0,3,0]
ierr=0
call rcss(h,n,xw,yw,wt,sy,trace,diag,cv,size(xnew),xg,ynew,job,ider,din,dout,ierr)
if(ierr/=0) error stop 'robust_spline_predict: rcss failed'
end function robust_spline_predict

end module fields_spline1d
