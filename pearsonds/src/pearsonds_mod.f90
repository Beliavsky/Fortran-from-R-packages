! SPDX-License-Identifier: GPL-2.0-or-later
!
! Modern Fortran translation of the computational code in the R package
! PearsonDS 1.3.2 (2025-03-24), by Martin Becker, Stefan Kloessner, and
! contributors. See NOTICE and LICENSES/ for attribution and license details.
! Pearson-IV normalization and rejection sampling are ports of code whose upstream
! comment credits Joel Heinrich, "A Guide to the Pearson Type IV Distribution"
! (University of Pennsylvania, 2004-12-21).
module pearsonds_mod
use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_is_finite
use r_compat, only: dp, dnorm, normal_cdf, qnorm, dbeta, pbeta, qbeta, dgamma, pgamma, qgamma, &
   dt, pt, qt, df, pf, qf, rnorm_vec, rbeta, rgamma, rt_vec, rf_rng, runif1, &
   optim_result_t, optim_bfgs, integrate_result_t, integrate
implicit none
private

integer, parameter, public :: pearson_type0 = 0
integer, parameter, public :: pearson_type_i = 1
integer, parameter, public :: pearson_type_ii = 2
integer, parameter, public :: pearson_type_iii = 3
integer, parameter, public :: pearson_type_iv = 4
integer, parameter, public :: pearson_type_v = 5
integer, parameter, public :: pearson_type_vi = 6
integer, parameter, public :: pearson_type_vii = 7

type, public :: pearson_params_t
   integer :: family = -1
   integer :: npar = 0
   real(kind=dp) :: par(4) = 0.0_dp
   integer :: status = 0
   character(len=160) :: message = ""
end type pearson_params_t

type, public :: pearson_ml_result_t
   type(pearson_params_t) :: fit
   real(kind=dp) :: objective = huge(1.0_dp)
   integer :: convergence = 1
   integer :: counts(2) = 0
   character(len=160) :: message = ""
end type pearson_ml_result_t

type, public :: pearson_msc_result_t
   real(kind=dp) :: loglik(8) = 0.0_dp
   real(kind=dp) :: criteria(5,8) = 0.0_dp
   type(pearson_params_t) :: fits(8)
   type(pearson_params_t) :: best(5)
end type pearson_msc_result_t

interface dpearson
   module procedure dpearson_scalar
   module procedure dpearson_vec
end interface dpearson

interface ppearson
   module procedure ppearson_scalar
   module procedure ppearson_vec
end interface ppearson

interface qpearson
   module procedure qpearson_scalar
   module procedure qpearson_vec
end interface qpearson

interface ppearsoniv
   module procedure ppearsoniv_scalar
   module procedure ppearsoniv_vec
end interface ppearsoniv

interface qpearsoniv
   module procedure qpearsoniv_scalar
   module procedure qpearsoniv_vec
end interface qpearsoniv

public :: dp
public :: dpearson0, ppearson0, qpearson0, rpearson0
public :: dpearsoni, ppearsoni, qpearsoni, rpearsoni
public :: dpearsonii, ppearsonii, qpearsonii, rpearsonii
public :: dpearsoniii, ppearsoniii, qpearsoniii, rpearsoniii
public :: dpearsoniv, ppearsoniv, qpearsoniv, rpearsoniv
public :: dpearsonv, ppearsonv, qpearsonv, rpearsonv
public :: dpearsonvi, ppearsonvi, qpearsonvi, rpearsonvi
public :: dpearsonvii, ppearsonvii, qpearsonvii, rpearsonvii
public :: dpearson, ppearson, qpearson, rpearson
public :: pearson0moments, pearsonimoments, pearsoniimoments, pearsoniiimoments
public :: pearsonivmoments, pearsonvmoments, pearsonvimoments, pearsonviimoments
public :: pearson_moments, emp_moments, pearson_fit_m, match_moments
public :: pearson_fit_ml_type, pearson_fit_ml, pearson_msc
public :: pearson_iv_norm, log_pearson_iv_norm, hypergeom_2f1

contains

pure elemental logical function close_enough(a, b) result(ok)
real(kind=dp), intent(in) :: a, b
real(kind=dp) :: sc
sc = max(1.0_dp, abs(a), abs(b))
ok = abs(a-b) <= 100.0_dp*epsilon(1.0_dp)*sc
end function close_enough

pure elemental real(kind=dp) function nan_dp() result(x)
x = ieee_value(0.0_dp, ieee_quiet_nan)
end function nan_dp

pure elemental real(kind=dp) function pos_inf() result(x)
x = ieee_value(0.0_dp, ieee_positive_inf)
end function pos_inf

pure elemental real(kind=dp) function tail_finish(p, lower_tail, log_p) result(out)
real(kind=dp), intent(in) :: p
logical, intent(in), optional :: lower_tail, log_p
logical :: lt, lp
lt = .true.
lp = .false.
if (present(lower_tail)) lt = lower_tail
if (present(log_p)) lp = log_p
out = p
if (.not. lt) out = 1.0_dp-out
if (lp) out = log(out)
end function tail_finish

pure elemental real(kind=dp) function probability_input(p, lower_tail, log_p) result(out)
real(kind=dp), intent(in) :: p
logical, intent(in), optional :: lower_tail, log_p
logical :: lt, lp
lt = .true.
lp = .false.
if (present(lower_tail)) lt = lower_tail
if (present(log_p)) lp = log_p
out = p
if (lp) out = exp(out)
if (.not. lt) out = 1.0_dp-out
end function probability_input

pure elemental real(kind=dp) function dpearson0(x, mean, sd, log_) result(out)
real(kind=dp), intent(in) :: x, mean, sd
logical, intent(in), optional :: log_
out = dnorm(x, mean=mean, sd=sd, log_=log_)
end function dpearson0

pure elemental real(kind=dp) function ppearson0(q, mean, sd, lower_tail, log_p) result(out)
real(kind=dp), intent(in) :: q, mean, sd
logical, intent(in), optional :: lower_tail, log_p
real(kind=dp) :: p
if (sd <= 0.0_dp) then
   out = nan_dp()
   return
end if
p = normal_cdf((q-mean)/sd)
out = tail_finish(p, lower_tail, log_p)
end function ppearson0

pure elemental real(kind=dp) function qpearson0(p, mean, sd, lower_tail, log_p) result(out)
real(kind=dp), intent(in) :: p, mean, sd
logical, intent(in), optional :: lower_tail, log_p
real(kind=dp) :: pp
pp = probability_input(p, lower_tail, log_p)
out = qnorm(pp, mean=mean, sd=sd)
end function qpearson0

function rpearson0(n, mean, sd) result(out)
integer, intent(in) :: n
real(kind=dp), intent(in) :: mean, sd
real(kind=dp), allocatable :: out(:)
out = mean + sd*rnorm_vec(n)
end function rpearson0

pure elemental real(kind=dp) function dpearsoni(x, a, b, location, scale, log_) result(out)
real(kind=dp), intent(in) :: x, a, b, location, scale
logical, intent(in), optional :: log_
logical :: lg
lg = .false.
if (present(log_)) lg = log_
if (a <= 0.0_dp .or. b <= 0.0_dp .or. scale == 0.0_dp) then
   out = nan_dp()
   return
end if
out = dbeta((x-location)/scale, a, b, log_=lg)
if (lg) then
   out = out-log(abs(scale))
else
   out = out/abs(scale)
end if
end function dpearsoni

pure elemental real(kind=dp) function ppearsoni(q, a, b, location, scale, lower_tail, log_p) result(out)
real(kind=dp), intent(in) :: q, a, b, location, scale
logical, intent(in), optional :: lower_tail, log_p
real(kind=dp) :: p
if (a <= 0.0_dp .or. b <= 0.0_dp .or. scale == 0.0_dp) then
   out = nan_dp()
   return
end if
p = pbeta((q-location)/scale, a, b)
if (scale < 0.0_dp) p = 1.0_dp-p
out = tail_finish(p, lower_tail, log_p)
end function ppearsoni

pure elemental real(kind=dp) function qpearsoni(p, a, b, location, scale, lower_tail, log_p) result(out)
real(kind=dp), intent(in) :: p, a, b, location, scale
logical, intent(in), optional :: lower_tail, log_p
real(kind=dp) :: pp
pp = probability_input(p, lower_tail, log_p)
if (scale < 0.0_dp) pp = 1.0_dp-pp
out = location + scale*qbeta(pp, a, b)
end function qpearsoni

function rpearsoni(n, a, b, location, scale) result(out)
integer, intent(in) :: n
real(kind=dp), intent(in) :: a, b, location, scale
real(kind=dp), allocatable :: out(:)
out = location + scale*rbeta(n, a, b)
end function rpearsoni

pure elemental real(kind=dp) function dpearsonii(x, a, location, scale, log_) result(out)
real(kind=dp), intent(in) :: x, a, location, scale
logical, intent(in), optional :: log_
out = dpearsoni(x, a, a, location, scale, log_)
end function dpearsonii

pure elemental real(kind=dp) function ppearsonii(q, a, location, scale, lower_tail, log_p) result(out)
real(kind=dp), intent(in) :: q, a, location, scale
logical, intent(in), optional :: lower_tail, log_p
out = ppearsoni(q, a, a, location, scale, lower_tail, log_p)
end function ppearsonii

pure elemental real(kind=dp) function qpearsonii(p, a, location, scale, lower_tail, log_p) result(out)
real(kind=dp), intent(in) :: p, a, location, scale
logical, intent(in), optional :: lower_tail, log_p
out = qpearsoni(p, a, a, location, scale, lower_tail, log_p)
end function qpearsonii

function rpearsonii(n, a, location, scale) result(out)
integer, intent(in) :: n
real(kind=dp), intent(in) :: a, location, scale
real(kind=dp), allocatable :: out(:)
out = location + scale*rbeta(n, a, a)
end function rpearsonii

pure elemental real(kind=dp) function dpearsoniii(x, shape, location, scale, log_) result(out)
real(kind=dp), intent(in) :: x, shape, location, scale
logical, intent(in), optional :: log_
real(kind=dp) :: y
if (shape <= 0.0_dp .or. scale == 0.0_dp) then
   out = nan_dp()
   return
end if
y = sign(1.0_dp, scale)*(x-location)
out = dgamma(y, shape, rate=1.0_dp/abs(scale), log_=log_)
end function dpearsoniii

pure elemental real(kind=dp) function ppearsoniii(q, shape, location, scale, lower_tail, log_p) result(out)
real(kind=dp), intent(in) :: q, shape, location, scale
logical, intent(in), optional :: lower_tail, log_p
real(kind=dp) :: p
if (shape <= 0.0_dp .or. scale == 0.0_dp) then
   out = nan_dp()
   return
end if
p = pgamma(sign(1.0_dp,scale)*(q-location), shape, rate=1.0_dp/abs(scale))
if (scale < 0.0_dp) p = 1.0_dp-p
out = tail_finish(p, lower_tail, log_p)
end function ppearsoniii

pure elemental real(kind=dp) function qpearsoniii(p, shape, location, scale, lower_tail, log_p) result(out)
real(kind=dp), intent(in) :: p, shape, location, scale
logical, intent(in), optional :: lower_tail, log_p
real(kind=dp) :: pp
pp = probability_input(p, lower_tail, log_p)
if (scale < 0.0_dp) pp = 1.0_dp-pp
out = location + sign(1.0_dp,scale)*qgamma(pp, shape, rate=1.0_dp/abs(scale))
end function qpearsoniii

function rpearsoniii(n, shape, location, scale) result(out)
integer, intent(in) :: n
real(kind=dp), intent(in) :: shape, location, scale
real(kind=dp), allocatable :: out(:)
out = location + sign(1.0_dp,scale)*rgamma(n, shape, scale=abs(scale))
end function rpearsoniii

pure real(kind=dp) function log_gammar2(x_in, y) result(out)
! Port of PearsonDS src/pearsonIV.c: loggammar2().
real(kind=dp), intent(in) :: x_in, y
real(kind=dp) :: x, y2, xmin, r, s, p, f, t
x = x_in
y2 = y*y
xmin = max(2.0_dp*y2, 10.0_dp)
r = 0.0_dp
s = 1.0_dp
p = 1.0_dp
f = 0.0_dp
do while (x < xmin)
   t = y/x
   r = r + log(1.0_dp+t*t)
   x = x + 1.0_dp
end do
do while (p > s*epsilon(1.0_dp))
   p = p*(y2+f*f)
   f = f+1.0_dp
   p = p/(x*f)
   x = x+1.0_dp
   s = s+p
end do
out = -r-log(s)
end function log_gammar2

pure elemental real(kind=dp) function log_pearson_iv_norm(m, nu, scale) result(out)
real(kind=dp), intent(in) :: m, nu, scale
real(kind=dp), parameter :: half = 0.5_dp
if (m <= half .or. scale <= 0.0_dp) then
   out = nan_dp()
   return
end if
out = -half*log(acos(-1.0_dp)) + log_gammar2(m, half*nu) + &
      log_gamma(m) - log_gamma(m-half) - log(scale)
end function log_pearson_iv_norm

pure elemental real(kind=dp) function pearson_iv_norm(m, nu, scale) result(out)
real(kind=dp), intent(in) :: m, nu, scale
out = exp(log_pearson_iv_norm(m, nu, scale))
end function pearson_iv_norm

pure elemental real(kind=dp) function dpearsoniv(x, m, nu, location, scale, log_) result(out)
real(kind=dp), intent(in) :: x, m, nu, location, scale
logical, intent(in), optional :: log_
logical :: lg
real(kind=dp) :: z, ld
lg = .false.
if (present(log_)) lg = log_
if (scale <= 0.0_dp .or. m <= 0.5_dp) then
   out = nan_dp()
   return
end if
z = (x-location)/scale
ld = log_pearson_iv_norm(m,nu,scale) - m*log(1.0_dp+z*z) - nu*atan(z)
if (lg) then
   out = ld
else
   out = exp(ld)
end if
end function dpearsoniv

function ppearsoniv_scalar(q, m, nu, location, scale, lower_tail, log_p, tol) result(out)
real(kind=dp), intent(in) :: q, m, nu, location, scale
logical, intent(in), optional :: lower_tail, log_p
real(kind=dp), intent(in), optional :: tol
real(kind=dp) :: out, eps, mode, pin
real(kind=dp) :: ninf, pinf
logical :: lt, lp
integer :: nsub
type(integrate_result_t) :: ir
if (q /= q) then
   out = q
   return
end if
if (.not. ieee_is_finite(q)) then
   out = tail_finish(merge(1.0_dp,0.0_dp,q>0.0_dp),lower_tail,log_p)
   return
end if
if (scale <= 0.0_dp .or. m <= 0.5_dp) then
   out = nan_dp()
   return
end if
lt = .true.
lp = .false.
if (present(lower_tail)) lt = lower_tail
if (present(log_p)) lp = log_p
eps = 1.0e-8_dp
if (present(tol)) eps = max(tol, epsilon(1.0_dp))
nsub = 160
mode = location-scale*nu/(2.0_dp*m)
pinf = pos_inf()
ninf = -pinf
if (q > mode) then
   ir = integrate(pdf, q, pinf, rel_tol=eps, subdivisions=nsub)
   pin = 1.0_dp-ir%value
else
   ir = integrate(pdf, ninf, q, rel_tol=eps, subdivisions=nsub)
   pin = ir%value
end if
pin = max(0.0_dp, min(1.0_dp, pin))
if (.not. lt) pin = 1.0_dp-pin
if (lp) then
   out = log(pin)
else
   out = pin
end if
contains
   function pdf(z) result(v)
   real(kind=dp), intent(in) :: z
   real(kind=dp) :: v
   v = dpearsoniv(z,m,nu,location,scale)
   end function pdf
end function ppearsoniv_scalar

function ppearsoniv_vec(q, m, nu, location, scale, lower_tail, log_p, tol) result(out)
real(kind=dp), intent(in) :: q(:), m, nu, location, scale
logical, intent(in), optional :: lower_tail, log_p
real(kind=dp), intent(in), optional :: tol
real(kind=dp), allocatable :: out(:)
integer :: i
allocate(out(size(q)))
do i=1,size(q)
   out(i) = ppearsoniv_scalar(q(i),m,nu,location,scale,lower_tail,log_p,tol)
end do
end function ppearsoniv_vec

function qpearsoniv_scalar(p, m, nu, location, scale, lower_tail, log_p, tol) result(out)
real(kind=dp), intent(in) :: p, m, nu, location, scale
logical, intent(in), optional :: lower_tail, log_p
real(kind=dp), intent(in), optional :: tol
real(kind=dp) :: out, pp, eps, xold, xnew, dens, cdf, ttol
integer :: it
if (scale <= 0.0_dp .or. m <= 0.5_dp) then
   out = nan_dp()
   return
end if
pp = probability_input(p, lower_tail, log_p)
if (pp /= pp .or. pp < 0.0_dp .or. pp > 1.0_dp) then
   out = nan_dp()
   return
else if (pp == 0.0_dp) then
   out = -pos_inf()
   return
else if (pp == 1.0_dp) then
   out = pos_inf()
   return
end if
eps = 1.0e-8_dp
if (present(tol)) eps = max(tol, epsilon(1.0_dp))
xold = location-scale*nu/(2.0_dp*m)
do it=1,30
   dens = dpearsoniv(xold,m,nu,location,scale)
   if (.not. ieee_is_finite(dens) .or. dens <= tiny(1.0_dp)) exit
   ttol = max(eps*1.0e-2_dp*dens, epsilon(1.0_dp))
   cdf = ppearsoniv_scalar(xold,m,nu,location,scale,tol=ttol)
   xnew = xold-(cdf-pp)/dens
   if (.not. ieee_is_finite(xnew)) exit
   if (abs(xnew-xold) < eps) then
      out = xnew
      return
   end if
   xold = xnew
end do
! Robust bracketed fallback if the R-style Newton iteration does not settle.
call iv_quantile_bisect(pp,m,nu,location,scale,eps,out)
end function qpearsoniv_scalar

subroutine iv_quantile_bisect(p,m,nu,location,scale,tol,out)
real(kind=dp), intent(in) :: p,m,nu,location,scale,tol
real(kind=dp), intent(out) :: out
real(kind=dp) :: lo,hi,mid,pl,ph,pm,step
integer :: it
lo = location-scale*nu/(2.0_dp*m)
hi = lo
pl = ppearsoniv_scalar(lo,m,nu,location,scale,tol=max(tol,1.0e-7_dp))
ph = pl
step = scale
if (p < pl) then
   do it=1,80
      lo = lo-step
      pl = ppearsoniv_scalar(lo,m,nu,location,scale,tol=max(tol,1.0e-7_dp))
      if (pl <= p) exit
      step = 2.0_dp*step
   end do
else
   do it=1,80
      hi = hi+step
      ph = ppearsoniv_scalar(hi,m,nu,location,scale,tol=max(tol,1.0e-7_dp))
      if (ph >= p) exit
      step = 2.0_dp*step
   end do
end if
if (p < ppearsoniv_scalar(location-scale*nu/(2.0_dp*m),m,nu,location,scale, &
   tol=max(tol,1.0e-7_dp))) hi = location-scale*nu/(2.0_dp*m)
if (p >= ppearsoniv_scalar(location-scale*nu/(2.0_dp*m),m,nu,location,scale, &
   tol=max(tol,1.0e-7_dp))) lo = location-scale*nu/(2.0_dp*m)
do it=1,80
   mid = 0.5_dp*(lo+hi)
   pm = ppearsoniv_scalar(mid,m,nu,location,scale,tol=max(tol,1.0e-7_dp))
   if (pm < p) then
      lo = mid
   else
      hi = mid
   end if
   if (abs(hi-lo) <= tol*max(1.0_dp,abs(mid))) exit
end do
out = 0.5_dp*(lo+hi)
end subroutine iv_quantile_bisect

function qpearsoniv_vec(p, m, nu, location, scale, lower_tail, log_p, tol) result(out)
real(kind=dp), intent(in) :: p(:), m, nu, location, scale
logical, intent(in), optional :: lower_tail, log_p
real(kind=dp), intent(in), optional :: tol
real(kind=dp), allocatable :: out(:)
integer :: i
allocate(out(size(p)))
do i=1,size(p)
   out(i) = qpearsoniv_scalar(p(i),m,nu,location,scale,lower_tail,log_p,tol)
end do
end function qpearsoniv_vec

function rpearsoniv(n, m, nu, location, scale) result(out)
! Translation of the rejection sampler in PearsonDS src/pearsonIV.c.
integer, intent(in) :: n
real(kind=dp), intent(in) :: m, nu, location, scale
real(kind=dp), allocatable :: out(:)
real(kind=dp) :: b, mode, cosm, r, rc, x, z, logk, u
integer :: i, s
allocate(out(max(0,n)))
if (scale <= 0.0_dp .or. m < 1.0_dp) then
   out = nan_dp()
   return
end if
b = 2.0_dp*m-2.0_dp
if (abs(b) <= tiny(1.0_dp)) then
   ! The original sampler is singular exactly at m=1; use inverse CDF here.
   do i=1,size(out)
      out(i) = qpearsoniv_scalar(runif1(),m,nu,location,scale,tol=1.0e-7_dp)
   end do
   return
end if
mode = atan(-nu/b)
cosm = b/sqrt(b*b+nu*nu)
r = b*log(cosm)-nu*mode
logk = log_pearson_iv_norm(m,nu,1.0_dp)
rc = exp(-r-logk)
do i=1,size(out)
   do
      s = 0
      z = 0.0_dp
      x = 4.0_dp*runif1()
      if (x > 2.0_dp) then
         x = x-2.0_dp
         s = 1
      end if
      if (x > 1.0_dp) then
         z = log(x-1.0_dp)
         x = 1.0_dp-z
      end if
      if (s == 1) then
         x = mode+rc*x
      else
         x = mode-rc*x
      end if
      if (abs(x) >= 0.5_dp*acos(-1.0_dp)) cycle
      u = max(runif1(),tiny(1.0_dp))
      if (z+log(u) <= b*log(cos(x))-nu*x-r) exit
   end do
   out(i) = scale*tan(x)+location
end do
end function rpearsoniv

pure elemental real(kind=dp) function dpearsonv(x, shape, location, scale, log_) result(out)
real(kind=dp), intent(in) :: x, shape, location, scale
logical, intent(in), optional :: log_
logical :: lg
real(kind=dp) :: y
lg=.false.
if(present(log_)) lg=log_
if (shape <= 0.0_dp .or. scale == 0.0_dp) then
   out=nan_dp()
   return
end if
y=sign(1.0_dp,scale)*(x-location)
if (y <= 0.0_dp) then
   if (lg) then
   out=-pos_inf()
   else
   out=0.0_dp
   end if
   return
end if
if (lg) then
   out=dgamma(1.0_dp/y,shape,rate=abs(scale),log_=.true.)-2.0_dp*log(y)
else
   out=dgamma(1.0_dp/y,shape,rate=abs(scale))/(x-location)**2
end if
end function dpearsonv

pure elemental real(kind=dp) function ppearsonv(q, shape, location, scale, lower_tail, log_p) result(out)
real(kind=dp), intent(in) :: q, shape, location, scale
logical, intent(in), optional :: lower_tail, log_p
real(kind=dp) :: y,p
if (shape <= 0.0_dp .or. scale == 0.0_dp) then
   out=nan_dp()
   return
end if
y=sign(1.0_dp,scale)*(q-location)
if (y <= 0.0_dp) then
   p=merge(0.0_dp,1.0_dp,scale>0.0_dp)
else
   p=pgamma(1.0_dp/y,shape,rate=abs(scale))
   if (scale > 0.0_dp) p=1.0_dp-p
end if
out=tail_finish(p,lower_tail,log_p)
end function ppearsonv

pure elemental real(kind=dp) function qpearsonv(p, shape, location, scale, lower_tail, log_p) result(out)
real(kind=dp), intent(in) :: p, shape, location, scale
logical, intent(in), optional :: lower_tail, log_p
real(kind=dp) :: pp,g
pp=probability_input(p,lower_tail,log_p)
if (scale > 0.0_dp) pp=1.0_dp-pp
g=qgamma(pp,shape,rate=abs(scale))
out=location+sign(1.0_dp,scale)/g
end function qpearsonv

function rpearsonv(n,shape,location,scale) result(out)
integer,intent(in)::n
real(kind=dp),intent(in)::shape,location,scale
real(kind=dp),allocatable::out(:)
real(kind=dp),allocatable::g(:)
g=rgamma(n,shape,rate=abs(scale))
out=location+sign(1.0_dp,scale)/g
end function rpearsonv

pure elemental real(kind=dp) function dpearsonvi(x,a,b,location,scale,log_) result(out)
real(kind=dp),intent(in)::x,a,b,location,scale
logical,intent(in),optional::log_
logical::lg
real(kind=dp)::ns
lg=.false.
if(present(log_))lg=log_
if(a<=0.0_dp.or.b<=0.0_dp.or.scale==0.0_dp)then
out=nan_dp()
return
end if
ns=scale*a/b
out=df((x-location)/ns,2.0_dp*a,2.0_dp*b,log_=lg)
if(lg)then
out=out-log(abs(ns))
else
out=out/abs(ns)
end if
end function dpearsonvi

pure elemental real(kind=dp) function ppearsonvi(q,a,b,location,scale,lower_tail,log_p) result(out)
real(kind=dp),intent(in)::q,a,b,location,scale
logical,intent(in),optional::lower_tail,log_p
real(kind=dp)::ns,p
if(a<=0.0_dp.or.b<=0.0_dp.or.scale==0.0_dp)then
out=nan_dp()
return
end if
ns=scale*a/b
p=pf((q-location)/ns,2.0_dp*a,2.0_dp*b)
if(ns<0.0_dp)p=1.0_dp-p
out=tail_finish(p,lower_tail,log_p)
end function ppearsonvi

pure elemental real(kind=dp) function qpearsonvi(p,a,b,location,scale,lower_tail,log_p) result(out)
real(kind=dp),intent(in)::p,a,b,location,scale
logical,intent(in),optional::lower_tail,log_p
real(kind=dp)::ns,pp
ns=scale*a/b
pp=probability_input(p,lower_tail,log_p)
if(ns<0.0_dp)pp=1.0_dp-pp
out=location+ns*qf(pp,2.0_dp*a,2.0_dp*b)
end function qpearsonvi

function rpearsonvi(n,a,b,location,scale) result(out)
integer,intent(in)::n
real(kind=dp),intent(in)::a,b,location,scale
real(kind=dp),allocatable::out(:)
real(kind=dp)::ns
ns=scale*a/b
out=location+ns*rf_rng(n,2.0_dp*a,2.0_dp*b)
end function rpearsonvi

pure elemental real(kind=dp) function dpearsonvii(x,dof,location,scale,log_) result(out)
real(kind=dp),intent(in)::x,dof,location,scale
logical,intent(in),optional::log_
logical::lg
lg=.false.
if(present(log_))lg=log_
if(dof<=0.0_dp.or.scale==0.0_dp)then
out=nan_dp()
return
end if
out=dt((x-location)/scale,dof,log_=lg)
if(lg)then
out=out-log(abs(scale))
else
out=out/abs(scale)
end if
end function dpearsonvii

pure elemental real(kind=dp) function ppearsonvii(q,dof,location,scale,lower_tail,log_p) result(out)
real(kind=dp),intent(in)::q,dof,location,scale
logical,intent(in),optional::lower_tail,log_p
real(kind=dp)::p
if(dof<=0.0_dp.or.scale==0.0_dp)then
out=nan_dp()
return
end if
p=pt((q-location)/scale,dof)
if(scale<0.0_dp)p=1.0_dp-p
out=tail_finish(p,lower_tail,log_p)
end function ppearsonvii

pure elemental real(kind=dp) function qpearsonvii(p,dof,location,scale,lower_tail,log_p) result(out)
real(kind=dp),intent(in)::p,dof,location,scale
logical,intent(in),optional::lower_tail,log_p
real(kind=dp)::pp
pp=probability_input(p,lower_tail,log_p)
if(scale<0.0_dp)pp=1.0_dp-pp
out=location+scale*qt(pp,dof)
end function qpearsonvii

function rpearsonvii(n,dof,location,scale) result(out)
integer,intent(in)::n
real(kind=dp),intent(in)::dof,location,scale
real(kind=dp),allocatable::out(:)
out=location+scale*rt_vec(n,dof)
end function rpearsonvii

function dpearson_scalar(x,params,log_) result(out)
real(kind=dp),intent(in)::x
type(pearson_params_t),intent(in)::params
logical,intent(in),optional::log_
real(kind=dp)::out
real(kind=dp),allocatable::tmp(:)
tmp=dpearson_vec([x],params,log_)
out=tmp(1)
end function dpearson_scalar

function ppearson_scalar(q,params,lower_tail,log_p,tol) result(out)
real(kind=dp),intent(in)::q
type(pearson_params_t),intent(in)::params
logical,intent(in),optional::lower_tail,log_p
real(kind=dp),intent(in),optional::tol
real(kind=dp)::out
real(kind=dp),allocatable::tmp(:)
tmp=ppearson_vec([q],params,lower_tail,log_p,tol)
out=tmp(1)
end function ppearson_scalar

function qpearson_scalar(p,params,lower_tail,log_p,tol) result(out)
real(kind=dp),intent(in)::p
type(pearson_params_t),intent(in)::params
logical,intent(in),optional::lower_tail,log_p
real(kind=dp),intent(in),optional::tol
real(kind=dp)::out
real(kind=dp),allocatable::tmp(:)
tmp=qpearson_vec([p],params,lower_tail,log_p,tol)
out=tmp(1)
end function qpearson_scalar

function dpearson_vec(x, params, log_) result(out)
real(kind=dp),intent(in)::x(:)
type(pearson_params_t),intent(in)::params
logical,intent(in),optional::log_
real(kind=dp),allocatable::out(:)
select case(params%family)
case(0); out=dpearson0(x,params%par(1),params%par(2),log_)
case(1); out=dpearsoni(x,params%par(1),params%par(2),params%par(3),params%par(4),log_)
case(2); out=dpearsonii(x,params%par(1),params%par(2),params%par(3),log_)
case(3); out=dpearsoniii(x,params%par(1),params%par(2),params%par(3),log_)
case(4); out=dpearsoniv(x,params%par(1),params%par(2),params%par(3),params%par(4),log_)
case(5); out=dpearsonv(x,params%par(1),params%par(2),params%par(3),log_)
case(6); out=dpearsonvi(x,params%par(1),params%par(2),params%par(3),params%par(4),log_)
case(7); out=dpearsonvii(x,params%par(1),params%par(2),params%par(3),log_)
case default
allocate(out(size(x)))
out=nan_dp()
end select
end function dpearson_vec

function ppearson_vec(q,params,lower_tail,log_p,tol) result(out)
real(kind=dp),intent(in)::q(:)
type(pearson_params_t),intent(in)::params
logical,intent(in),optional::lower_tail,log_p
real(kind=dp),intent(in),optional::tol
real(kind=dp),allocatable::out(:)
select case(params%family)
case(0); out=ppearson0(q,params%par(1),params%par(2),lower_tail,log_p)
case(1); out=ppearsoni(q,params%par(1),params%par(2),params%par(3),params%par(4),lower_tail,log_p)
case(2); out=ppearsonii(q,params%par(1),params%par(2),params%par(3),lower_tail,log_p)
case(3); out=ppearsoniii(q,params%par(1),params%par(2),params%par(3),lower_tail,log_p)
case(4); out=ppearsoniv(q,params%par(1),params%par(2),params%par(3),params%par(4),lower_tail,log_p,tol)
case(5); out=ppearsonv(q,params%par(1),params%par(2),params%par(3),lower_tail,log_p)
case(6); out=ppearsonvi(q,params%par(1),params%par(2),params%par(3),params%par(4),lower_tail,log_p)
case(7); out=ppearsonvii(q,params%par(1),params%par(2),params%par(3),lower_tail,log_p)
case default
allocate(out(size(q)))
out=nan_dp()
end select
end function ppearson_vec

function qpearson_vec(p,params,lower_tail,log_p,tol) result(out)
real(kind=dp),intent(in)::p(:)
type(pearson_params_t),intent(in)::params
logical,intent(in),optional::lower_tail,log_p
real(kind=dp),intent(in),optional::tol
real(kind=dp),allocatable::out(:)
select case(params%family)
case(0); out=qpearson0(p,params%par(1),params%par(2),lower_tail,log_p)
case(1); out=qpearsoni(p,params%par(1),params%par(2),params%par(3),params%par(4),lower_tail,log_p)
case(2); out=qpearsonii(p,params%par(1),params%par(2),params%par(3),lower_tail,log_p)
case(3); out=qpearsoniii(p,params%par(1),params%par(2),params%par(3),lower_tail,log_p)
case(4); out=qpearsoniv(p,params%par(1),params%par(2),params%par(3),params%par(4),lower_tail,log_p,tol)
case(5); out=qpearsonv(p,params%par(1),params%par(2),params%par(3),lower_tail,log_p)
case(6); out=qpearsonvi(p,params%par(1),params%par(2),params%par(3),params%par(4),lower_tail,log_p)
case(7); out=qpearsonvii(p,params%par(1),params%par(2),params%par(3),lower_tail,log_p)
case default
allocate(out(size(p)))
out=nan_dp()
end select
end function qpearson_vec

function rpearson(n,params) result(out)
integer,intent(in)::n
type(pearson_params_t),intent(in)::params
real(kind=dp),allocatable::out(:)
select case(params%family)
case(0); out=rpearson0(n,params%par(1),params%par(2))
case(1); out=rpearsoni(n,params%par(1),params%par(2),params%par(3),params%par(4))
case(2); out=rpearsonii(n,params%par(1),params%par(2),params%par(3))
case(3); out=rpearsoniii(n,params%par(1),params%par(2),params%par(3))
case(4); out=rpearsoniv(n,params%par(1),params%par(2),params%par(3),params%par(4))
case(5); out=rpearsonv(n,params%par(1),params%par(2),params%par(3))
case(6); out=rpearsonvi(n,params%par(1),params%par(2),params%par(3),params%par(4))
case(7); out=rpearsonvii(n,params%par(1),params%par(2),params%par(3))
case default
allocate(out(max(0,n)))
out=nan_dp()
end select
end function rpearson

pure function pearson0moments(mean,sd) result(z)
real(kind=dp),intent(in)::mean,sd
real(kind=dp)::z(4)
z=[mean,sd*sd,0.0_dp,3.0_dp]
end function pearson0moments

pure function pearsonimoments(a,b,location,scale) result(z)
real(kind=dp),intent(in)::a,b,location,scale
real(kind=dp)::z(4)
z(1)=scale*a/(a+b)+location
z(2)=scale**2*a*b/((a+b)**2*(a+b+1.0_dp))
z(3)=sign(1.0_dp,scale)*2.0_dp*(b-a)*sqrt(a+b+1.0_dp)/((a+b+2.0_dp)*sqrt(a*b))
z(4)=3.0_dp+6.0_dp*(a**3-a*a*(2.0_dp*b-1.0_dp)+b*b*(b+1.0_dp)-2.0_dp*a*b*(b+2.0_dp))/ &
     (a*b*(a+b+2.0_dp)*(a+b+3.0_dp))
end function pearsonimoments

pure function pearsoniimoments(a,location,scale) result(z)
real(kind=dp),intent(in)::a,location,scale
real(kind=dp)::z(4)
z(1)=scale/2.0_dp+location
z(2)=scale**2*a*a/((2.0_dp*a)**2*(2.0_dp*a+1.0_dp))
z(3)=0.0_dp
z(4)=3.0_dp+6.0_dp*(-2.0_dp*(a+1.0_dp))/((2.0_dp*a+2.0_dp)*(2.0_dp*a+3.0_dp))
end function pearsoniimoments

pure function pearsoniiimoments(shape,location,scale) result(z)
real(kind=dp),intent(in)::shape,location,scale
real(kind=dp)::z(4)
z=[scale*shape+location,shape*scale**2,sign(1.0_dp,scale)*2.0_dp/sqrt(shape),3.0_dp+6.0_dp/shape]
end function pearsoniiimoments

pure function pearsonivmoments(m,nu,location,scale) result(z)
real(kind=dp),intent(in)::m,nu,location,scale
real(kind=dp)::z(4),r
r=2.0_dp*(m-1.0_dp)
if(m>1.0_dp)then
z(1)=location-scale*nu/r
else
z(1)=nan_dp()
end if
if(m>1.5_dp)then
z(2)=scale**2*(r*r+nu*nu)/(r*r*(r-1.0_dp))
else
z(2)=pos_inf()
end if
if(m>2.0_dp)then
z(3)=-4.0_dp*nu/(r-2.0_dp)*sqrt((r-1.0_dp)/(r*r+nu*nu))
else
z(3)=nan_dp()
end if
if(m>2.5_dp)then
 z(4)=3.0_dp*(r-1.0_dp)*((r+6.0_dp)*(r*r+nu*nu)-8.0_dp*r*r)/((r-2.0_dp)*(r-3.0_dp)*(r*r+nu*nu))
else
 z(4)=nan_dp()
end if
end function pearsonivmoments

pure function pearsonvmoments(shape,location,scale) result(z)
real(kind=dp),intent(in)::shape,location,scale
real(kind=dp)::z(4)
if(shape>1.0_dp)then
z(1)=location+scale/(shape-1.0_dp)
else
z(1)=nan_dp()
end if
if(shape>2.0_dp)then
z(2)=scale**2/((shape-1.0_dp)**2*(shape-2.0_dp))
else
z(2)=nan_dp()
end if
if(shape>3.0_dp)then
z(3)=4.0_dp*sign(1.0_dp,scale)*sqrt(shape-2.0_dp)/(shape-3.0_dp)
else
z(3)=nan_dp()
end if
if(shape>4.0_dp)then
z(4)=3.0_dp+(30.0_dp*shape-66.0_dp)/((shape-3.0_dp)*(shape-4.0_dp))
else
z(4)=nan_dp()
end if
end function pearsonvmoments

pure function pearsonvimoments(a,b,location,scale) result(z)
real(kind=dp),intent(in)::a,b,location,scale
real(kind=dp)::z(4),mmm,vvn,vvv,ssn,sss,kkn,kkk
mmm=scale*a/(b-1.0_dp)
vvn=scale**2*(a+1.0_dp)*a/((b-1.0_dp)*(b-2.0_dp))
vvv=vvn-mmm**2
ssn=scale**3*(a+2.0_dp)*(a+1.0_dp)*a/((b-1.0_dp)*(b-2.0_dp)*(b-3.0_dp))
sss=ssn-3.0_dp*mmm*vvn+2.0_dp*mmm**3
kkn=scale**4*(a+3.0_dp)*(a+2.0_dp)*(a+1.0_dp)*a/((b-1.0_dp)*(b-2.0_dp)*(b-3.0_dp)*(b-4.0_dp))
kkk=kkn-4.0_dp*mmm*ssn+6.0_dp*mmm*mmm*vvn-3.0_dp*mmm**4
if(b>1.0_dp)then
z(1)=mmm+location
else
z(1)=nan_dp()
end if
if(b>2.0_dp)then
z(2)=scale**2*a*(a+b-1.0_dp)/((b-2.0_dp)*(b-1.0_dp)**2)
else
z(2)=nan_dp()
end if
if(b>3.0_dp)then
z(3)=sss/(vvv**1.5_dp)
else
z(3)=nan_dp()
end if
if(b>4.0_dp)then
z(4)=kkk/(vvv**2)
else
z(4)=nan_dp()
end if
end function pearsonvimoments

pure function pearsonviimoments(dof,location,scale) result(z)
real(kind=dp),intent(in)::dof,location,scale
real(kind=dp)::z(4)
if(dof>1.0_dp)then
z(1)=location
else
z(1)=nan_dp()
end if
if(dof>2.0_dp)then
z(2)=scale**2*dof/(dof-2.0_dp)
else
z(2)=pos_inf()
end if
if(dof>3.0_dp)then
z(3)=0.0_dp
else
z(3)=nan_dp()
end if
if(dof>4.0_dp)then
z(4)=3.0_dp+6.0_dp/(dof-4.0_dp)
else
z(4)=nan_dp()
end if
end function pearsonviimoments

pure function pearson_moments(params) result(z)
type(pearson_params_t),intent(in)::params
real(kind=dp)::z(4)
select case(params%family)
case(0);z=pearson0moments(params%par(1),params%par(2))
case(1);z=pearsonimoments(params%par(1),params%par(2),params%par(3),params%par(4))
case(2);z=pearsoniimoments(params%par(1),params%par(2),params%par(3))
case(3);z=pearsoniiimoments(params%par(1),params%par(2),params%par(3))
case(4);z=pearsonivmoments(params%par(1),params%par(2),params%par(3),params%par(4))
case(5);z=pearsonvmoments(params%par(1),params%par(2),params%par(3))
case(6);z=pearsonvimoments(params%par(1),params%par(2),params%par(3),params%par(4))
case(7);z=pearsonviimoments(params%par(1),params%par(2),params%par(3))
case default;z=nan_dp()
end select
end function pearson_moments

pure function emp_moments(x) result(z)
real(kind=dp),intent(in)::x(:)
real(kind=dp)::z(4),mu,v
integer::n
n=size(x)
if(n==0.or.any(x/=x))then
z=nan_dp()
return
end if
mu=sum(x)/real(n,dp)
v=sum((x-mu)**2)/real(n,dp)
z(1)=mu
z(2)=v
if(v>0.0_dp)then
 z(3)=sum((x-mu)**3)/(real(n,dp)*v**1.5_dp)
 z(4)=sum((x-mu)**4)/(real(n,dp)*v**2)
else
 z(3)=0.0_dp
 z(4)=0.0_dp
end if
end function emp_moments

pure function pearson_fit_m(mean,variance,skewness,kurtosis) result(res)
real(kind=dp),intent(in)::mean,variance,skewness,kurtosis
type(pearson_params_t)::res
real(kind=dp)::mmm,vvv,sss,kkk,c0,c1,c2,kap,a1,a2,a,m1,m2,sca,loc,r,nu,disc,c1x
mmm=mean
vvv=variance
sss=skewness
kkk=kurtosis
res%family=-1
res%npar=0
res%par=0.0_dp
res%status=0
res%message=""
if(.not.ieee_is_finite(vvv).or.vvv<=0.0_dp.or.any([mmm,sss,kkk]/=[mmm,sss,kkk]))then
 res%status=1
 res%message="invalid moments"
 return
end if
if(close_enough(sss*sss,kkk-1.0_dp))then
 res%status=2
 res%message="boundary two-point distribution is not in the Pearson continuous system"
 return
end if
if(sss*sss>kkk-1.0_dp)then
 res%status=3
 res%message="no probability distribution has the requested moments"
 return
end if
c0=(4.0_dp*kkk-3.0_dp*sss*sss)/(10.0_dp*kkk-12.0_dp*sss*sss-18.0_dp)*vvv
c1=sss*(kkk+3.0_dp)/(10.0_dp*kkk-12.0_dp*sss*sss-18.0_dp)*sqrt(vvv)
c2=(2.0_dp*kkk-3.0_dp*sss*sss-6.0_dp)/(10.0_dp*kkk-12.0_dp*sss*sss-18.0_dp)
if(close_enough(sss,0.0_dp))then
 if(close_enough(kkk,3.0_dp))then
  res%family=0
  res%npar=2
  res%par(1:2)=[mmm,sqrt(vvv)]
  return
 else if(kkk<3.0_dp)then
  a1=sqrt(vvv)/2.0_dp*(-sqrt(-16.0_dp*kkk*(2.0_dp*kkk-6.0_dp))/(2.0_dp*kkk-6.0_dp))
  m1=-1.0_dp/(2.0_dp*c2)
  sca=2.0_dp*a1
  loc=mmm-sca/2.0_dp
  res%family=2
  res%npar=3
  res%par(1:3)=[1.0_dp+m1,loc,sca]
  return
 else
  r=6.0_dp*(kkk-1.0_dp)/(2.0_dp*kkk-6.0_dp)
  a=sqrt(vvv*(r-1.0_dp))
  res%family=7
  res%npar=3
  res%par(1:3)=[1.0_dp+r,mmm,a/sqrt(1.0_dp+r)]
  return
 end if
else if(.not.close_enough(2.0_dp*kkk-3.0_dp*sss*sss-6.0_dp,0.0_dp))then
 kap=0.25_dp*sss*sss*(kkk+3.0_dp)**2/((4.0_dp*kkk-3.0_dp*sss*sss)*(2.0_dp*kkk-3.0_dp*sss*sss-6.0_dp))
 disc=sss*sss*(kkk+3.0_dp)**2-4.0_dp*(4.0_dp*kkk-3.0_dp*sss*sss)*(2.0_dp*kkk-3.0_dp*sss*sss-6.0_dp)
 if(kap<0.0_dp)then
  a1=sqrt(vvv)/2.0_dp*((-sss*(kkk+3.0_dp)-sqrt(disc))/(2.0_dp*kkk-3.0_dp*sss*sss-6.0_dp))
  a2=sqrt(vvv)/2.0_dp*((-sss*(kkk+3.0_dp)+sqrt(disc))/(2.0_dp*kkk-3.0_dp*sss*sss-6.0_dp))
  if(a1>0.0_dp)then
  a=a1
  a1=a2
  a2=a
  end if
  m1=-(sss*(kkk+3.0_dp)+a1*(10.0_dp*kkk-12.0_dp*sss*sss-18.0_dp)/sqrt(vvv))/sqrt(disc)
  m2=-(-sss*(kkk+3.0_dp)-a2*(10.0_dp*kkk-12.0_dp*sss*sss-18.0_dp)/sqrt(vvv))/sqrt(disc)
  sca=a2-a1
  loc=mmm-sca*(m1+1.0_dp)/(m1+m2+2.0_dp)
  res%family=1
  res%npar=4
  res%par=[1.0_dp+m1,1.0_dp+m2,loc,sca]
  return
 else if(close_enough(kap,1.0_dp))then
  c1x=c1/(2.0_dp*c2)
  sca=-(c1-c1x)/c2
  res%family=5
  res%npar=3
  res%par(1:3)=[1.0_dp/c2-1.0_dp,mmm-c1x,sca]
  return
 else if(kap>1.0_dp)then
  a1=sqrt(vvv)/2.0_dp*((-sss*(kkk+3.0_dp)-sqrt(disc))/(2.0_dp*kkk-3.0_dp*sss*sss-6.0_dp))
  a2=sqrt(vvv)/2.0_dp*((-sss*(kkk+3.0_dp)+sqrt(disc))/(2.0_dp*kkk-3.0_dp*sss*sss-6.0_dp))
  if(a1>0.0_dp)then
  a=a1
  a1=a2
  a2=a
  end if
  a=c1
  m1=(a+a1)/(c2*(a2-a1))
  m2=-(a+a2)/(c2*(a2-a1))
  sca=a2-a1
  loc=mmm+sca*(m2+1.0_dp)/(m2+m1+2.0_dp)
  res%family=6
  res%npar=4
  res%par=[1.0_dp+m2,-m2-m1-1.0_dp,loc,sca]
  return
 else if(kap>0.0_dp.and.kap<1.0_dp)then
  r=6.0_dp*(kkk-sss*sss-1.0_dp)/(2.0_dp*kkk-3.0_dp*sss*sss-6.0_dp)
  nu=-r*(r-2.0_dp)*sss/sqrt(16.0_dp*(r-1.0_dp)-sss*sss*(r-2.0_dp)**2)
  sca=sqrt(vvv*(16.0_dp*(r-1.0_dp)-sss*sss*(r-2.0_dp)**2))/4.0_dp
  loc=mmm-((r-2.0_dp)*sss*sqrt(vvv))/4.0_dp
  res%family=4
  res%npar=4
  res%par=[1.0_dp+r/2.0_dp,nu,loc,sca]
  return
 end if
else
 a=c0/(c1*c1)-1.0_dp
 loc=mmm-c0/c1
 res%family=3
 res%npar=3
 res%par(1:3)=[a+1.0_dp,loc,c1]
 return
end if
res%status=4
res%message="moments did not map to a Pearson family"
end function pearson_fit_m

pure function adjusted_start(mom,family) result(res)
real(kind=dp),intent(in)::mom(4)
integer,intent(in)::family
type(pearson_params_t)::res
real(kind=dp)::mmm,vvv,sss,kkk,slim1,slim2
mmm=mom(1)
vvv=mom(2)
sss=mom(3)
kkk=mom(4)
select case(family)
case(0)
 res=pearson_fit_m(mmm,vvv,0.0_dp,3.0_dp)
case(1)
 if(close_enough(kkk,1.0_dp).or.kkk<1.0_dp)kkk=2.0_dp
 slim1=max(0.0_dp,-2.0_dp+2.0_dp*kkk/3.0_dp)
 slim2=kkk-1.0_dp
 if(close_enough(sss*sss,slim1).or.close_enough(sss*sss,slim2).or.sss*sss>slim2.or.sss*sss<slim1) &
  sss=sign(1.0_dp,sss)*0.5_dp*(sqrt(slim1)+sqrt(slim2))
 res=pearson_fit_m(mmm,vvv,sss,kkk)
case(2)
 if(close_enough(kkk,1.0_dp).or.close_enough(kkk,3.0_dp).or.kkk<1.0_dp.or.kkk>3.0_dp)kkk=2.0_dp
 res=pearson_fit_m(mmm,vvv,0.0_dp,kkk)
case(3)
 if(close_enough(kkk,3.0_dp).or.kkk<3.0_dp)kkk=4.0_dp
 res=pearson_fit_m(mmm,vvv,sign(1.0_dp,sss)*sqrt(-2.0_dp+2.0_dp*kkk/3.0_dp),kkk)
case(4)
 if(close_enough(kkk,3.0_dp).or.kkk<3.0_dp)kkk=4.0_dp
 slim2=(kkk**2+78.0_dp*kkk-63.0_dp-sqrt(kkk+147.0_dp)*sqrt(kkk+3.0_dp)**3)/72.0_dp
 if(close_enough(sss,0.0_dp).or.close_enough(sss*sss,slim2).or.sss*sss>slim2) &
  sss=sign(1.0_dp,sss)*0.5_dp*sqrt(max(0.0_dp,slim2))
 if(close_enough(sss,0.0_dp))sss=0.25_dp*sqrt(max(0.0_dp,slim2))
 res=pearson_fit_m(mmm,vvv,sss,kkk)
case(5)
 if(close_enough(kkk,3.0_dp).or.kkk<3.0_dp)kkk=4.0_dp
 slim1=(kkk**2+78.0_dp*kkk-63.0_dp-sqrt(kkk+147.0_dp)*sqrt(kkk+3.0_dp)**3)/72.0_dp
 res=pearson_fit_m(mmm,vvv,sign(1.0_dp,sss)*sqrt(max(0.0_dp,slim1)),kkk)
case(6)
 if(close_enough(kkk,3.0_dp).or.kkk<3.0_dp)kkk=4.0_dp
 slim1=(kkk**2+78.0_dp*kkk-63.0_dp-sqrt(kkk+147.0_dp)*sqrt(kkk+3.0_dp)**3)/72.0_dp
 slim2=-2.0_dp+2.0_dp*kkk/3.0_dp
 if(close_enough(sss*sss,slim1).or.close_enough(sss*sss,slim2).or.sss*sss>slim2.or.sss*sss<slim1) &
  sss=sign(1.0_dp,sss)*0.5_dp*(sqrt(max(0.0_dp,slim1))+sqrt(max(0.0_dp,slim2)))
 res=pearson_fit_m(mmm,vvv,sss,kkk)
case(7)
 if(close_enough(kkk,3.0_dp).or.kkk<3.0_dp)kkk=4.0_dp
 res=pearson_fit_m(mmm,vvv,0.0_dp,kkk)
case default
 res%status=1
 res%message="invalid Pearson family"
end select
end function adjusted_start

function pearson_fit_ml_type(x,family,maxit,reltol) result(res)
! Maximum-likelihood fit. PearsonDS uses nlminb() with direct parameters.
! Here the same log-likelihoods are optimized after smooth transformations that
! enforce positivity and sample-support constraints, avoiding invalid finite-
! difference steps at hard support boundaries.
real(kind=dp),intent(in)::x(:)
integer,intent(in)::family
integer,intent(in),optional::maxit
real(kind=dp),intent(in),optional::reltol
type(pearson_ml_result_t)::res
type(pearson_params_t)::start
real(kind=dp)::mom(4),rt,xmin,xmax,sgn,low0,high0,a0,b0,scale0
real(kind=dp),allocatable::theta(:)
real(kind=dp)::decoded(4)
type(optim_result_t)::opt
integer::mi,n
logical::ok
n=size(x)
mi=300
if(present(maxit))mi=maxit
rt=1.0e-8_dp
if(present(reltol))rt=reltol
if(n<2.or.any(x/=x))then
res%fit%status=1
res%message="invalid data"
return
end if
mom=emp_moments(x)
xmin=minval(x)
xmax=maxval(x)
if(family==0)then
 res%fit%family=0
 res%fit%npar=2
 res%fit%par(1)=mom(1)
 res%fit%par(2)=sqrt(mom(2))
 res%objective=-sum(dpearson0(x,res%fit%par(1),res%fit%par(2),log_=.true.))
 res%convergence=0
 res%counts=[1,0]
 return
end if
start=adjusted_start(mom,family)
if(start%status/=0.or.start%family/=family)then
 res%fit=start
 res%fit%status=2
 res%message="could not construct family-specific moment start"
 return
end if
sgn=1.0_dp
select case(family)
case(1)
 a0=start%par(1)
 b0=start%par(2)
 if(start%par(4)>0.0_dp)then
  low0=start%par(3)
  high0=start%par(3)+start%par(4)
 else
  low0=start%par(3)+start%par(4)
  high0=start%par(3)
  scale0=a0
  a0=b0
  b0=scale0
 end if
 low0=min(low0,xmin-0.1_dp)
 high0=max(high0,xmax+0.1_dp)
 theta=[log(max(a0,1.0e-8_dp)),log(max(b0,1.0e-8_dp)), &
        log(max(xmin-low0,1.0e-8_dp)),log(max(high0-xmax,1.0e-8_dp))]
case(2)
 low0=start%par(2)
 high0=start%par(2)+start%par(3)
 if(start%par(3)<0.0_dp)then
  low0=start%par(2)+start%par(3)
  high0=start%par(2)
 end if
 low0=min(low0,xmin-0.01_dp)
 high0=max(high0,xmax+0.01_dp)
 theta=[log(max(start%par(1),1.0e-8_dp)),log(max(xmin-low0,1.0e-8_dp)), &
        log(max(high0-xmax,1.0e-8_dp))]
case(3)
 sgn=sign(1.0_dp,start%par(3))
 scale0=max(abs(start%par(3)),1.0e-8_dp)
 if(sgn>0.0_dp)then
  low0=min(start%par(2),xmin-0.01_dp)
  theta=[log(max(start%par(1),1.0e-8_dp)), &
       log(max(xmin-low0,1.0e-8_dp)),log(scale0)]
 else
  high0=max(start%par(2),xmax+0.01_dp)
  theta=[log(max(start%par(1),1.0e-8_dp)), &
       log(max(high0-xmax,1.0e-8_dp)),log(scale0)]
 end if
case(4)
 theta=[log(max(start%par(1)-0.5_dp,1.0e-8_dp)),start%par(2),start%par(3), &
        log(max(start%par(4),1.0e-8_dp))]
case(5)
 sgn=sign(1.0_dp,start%par(3))
 scale0=max(abs(start%par(3)),1.0e-8_dp)
 if(sgn>0.0_dp)then
  low0=min(start%par(2),xmin-0.01_dp)
  theta=[log(max(start%par(1),1.0e-8_dp)), &
       log(max(xmin-low0,1.0e-8_dp)),log(scale0)]
 else
  high0=max(start%par(2),xmax+0.01_dp)
  theta=[log(max(start%par(1),1.0e-8_dp)), &
       log(max(high0-xmax,1.0e-8_dp)),log(scale0)]
 end if
case(6)
 sgn=sign(1.0_dp,start%par(4))
 scale0=max(abs(start%par(4)),1.0e-8_dp)
 if(sgn>0.0_dp)then
  low0=min(start%par(3),xmin-0.01_dp)
  theta=[log(max(start%par(1),1.0e-8_dp)), &
       log(max(start%par(2),1.0e-8_dp)),log(max(xmin-low0,1.0e-8_dp)),log(scale0)]
 else
  high0=max(start%par(3),xmax+0.01_dp)
  theta=[log(max(start%par(1),1.0e-8_dp)), &
       log(max(start%par(2),1.0e-8_dp)),log(max(high0-xmax,1.0e-8_dp)),log(scale0)]
 end if
case(7)
 theta=[log(max(start%par(1),1.0e-8_dp)),start%par(2),log(max(abs(start%par(3)),1.0e-8_dp))]
end select
opt=optim_bfgs(objective,theta,maxit=mi,reltol=rt)
call decode_theta(opt%par,decoded,ok)
res%fit%family=family
res%fit%npar=parameter_count(family)
res%fit%par=0.0_dp
if(ok)res%fit%par(1:res%fit%npar)=decoded(1:res%fit%npar)
res%objective=opt%value
res%convergence=opt%convergence
res%counts=opt%counts
if(allocated(opt%message))res%message=opt%message
if(.not.ok.or..not.ieee_is_finite(res%objective).or.res%objective>=1.0e90_dp)then
 res%convergence=max(10,res%convergence)
 res%fit%status=3
 if(len_trim(res%message)==0)res%message="non-finite or invalid likelihood fit"
end if
contains
 pure integer function parameter_count(fam) result(k)
 integer,intent(in)::fam
 integer,parameter::kk(0:7)=[2,4,3,3,4,3,4,3]
 k=kk(fam)
 end function parameter_count

 pure subroutine decode_theta(t,pout,valid)
 real(kind=dp),intent(in)::t(:)
 real(kind=dp),intent(out)::pout(4)
 logical,intent(out)::valid
 real(kind=dp)::e1,e2,e3,e4,lo,hi
 pout=0.0_dp
 valid=.false.
 if(any(.not.ieee_is_finite(t)))return
 select case(family)
 case(1)
  if(any(abs(t([1,2,3,4]))>80.0_dp))return
  e1=exp(t(1))
  e2=exp(t(2))
  e3=exp(t(3))
  e4=exp(t(4))
  lo=xmin-e3
  hi=xmax+e4
  pout=[e1,e2,lo,hi-lo]
 case(2)
  if(any(abs(t([1,2,3]))>80.0_dp))return
  e1=exp(t(1))
  e2=exp(t(2))
  e3=exp(t(3))
  lo=xmin-e2
  hi=xmax+e3
  pout(1:3)=[e1,lo,hi-lo]
 case(3)
  if(any(abs(t([1,2,3]))>80.0_dp))return
  e1=exp(t(1))
  e2=exp(t(2))
  e3=exp(t(3))
  if(sgn>0.0_dp)then
  lo=xmin-e2
  else
  lo=xmax+e2
  end if
  pout(1:3)=[e1,lo,sgn*e3]
 case(4)
  if(abs(t(1))>80.0_dp.or.abs(t(4))>80.0_dp)return
  pout=[0.5_dp+exp(t(1)),t(2),t(3),exp(t(4))]
 case(5)
  if(any(abs(t([1,2,3]))>80.0_dp))return
  e1=exp(t(1))
  e2=exp(t(2))
  e3=exp(t(3))
  if(sgn>0.0_dp)then
  lo=xmin-e2
  else
  lo=xmax+e2
  end if
  pout(1:3)=[e1,lo,sgn*e3]
 case(6)
  if(any(abs(t([1,2,3,4]))>80.0_dp))return
  e1=exp(t(1))
  e2=exp(t(2))
  e3=exp(t(3))
  e4=exp(t(4))
  if(sgn>0.0_dp)then
  lo=xmin-e3
  else
  lo=xmax+e3
  end if
  pout=[e1,e2,lo,sgn*e4]
 case(7)
  if(abs(t(1))>80.0_dp.or.abs(t(3))>80.0_dp)return
  pout(1:3)=[exp(t(1)),t(2),exp(t(3))]
 end select
 valid=all(ieee_is_finite(pout(1:parameter_count(family))))
 end subroutine decode_theta

 pure function objective(t) result(v)
 real(kind=dp),intent(in)::t(:)
 real(kind=dp)::v,pv(4)
 real(kind=dp),allocatable::ld(:)
 logical::valid
 call decode_theta(t,pv,valid)
 if(.not.valid)then
  v=1.0e100_dp
  return
 end if
 select case(family)
 case(1);ld=dpearsoni(x,pv(1),pv(2),pv(3),pv(4),log_=.true.)
 case(2);ld=dpearsonii(x,pv(1),pv(2),pv(3),log_=.true.)
 case(3);ld=dpearsoniii(x,pv(1),pv(2),pv(3),log_=.true.)
 case(4);ld=dpearsoniv(x,pv(1),pv(2),pv(3),pv(4),log_=.true.)
 case(5);ld=dpearsonv(x,pv(1),pv(2),pv(3),log_=.true.)
 case(6);ld=dpearsonvi(x,pv(1),pv(2),pv(3),pv(4),log_=.true.)
 case(7);ld=dpearsonvii(x,pv(1),pv(2),pv(3),log_=.true.)
 end select
 if(any(.not.ieee_is_finite(ld)))then
 v=1.0e100_dp
 else
 v=-sum(ld)
 end if
 end function objective
end function pearson_fit_ml_type

function pearson_fit_ml(x,maxit,reltol) result(best)
real(kind=dp),intent(in)::x(:)
integer,intent(in),optional::maxit
real(kind=dp),intent(in),optional::reltol
type(pearson_ml_result_t)::best,tmp
integer::fam
best%objective=huge(1.0_dp)
do fam=0,7
 tmp=pearson_fit_ml_type(x,fam,maxit,reltol)
 if(tmp%objective<best%objective)best=tmp
end do
end function pearson_fit_ml

function pearson_msc(x,maxit,reltol) result(out)
real(kind=dp),intent(in)::x(:)
integer,intent(in),optional::maxit
real(kind=dp),intent(in),optional::reltol
type(pearson_msc_result_t)::out
type(pearson_ml_result_t)::ml
integer,parameter::k(8)=[2,4,3,3,4,3,4,3]
integer::j,b
real(kind=dp)::n
n=real(size(x),dp)
do j=1,8
 ml=pearson_fit_ml_type(x,j-1,maxit,reltol)
 out%fits(j)=ml%fit
 out%loglik(j)=-ml%objective
 out%criteria(1,j)=-2.0_dp*out%loglik(j)
 out%criteria(2,j)=2.0_dp*k(j)-2.0_dp*out%loglik(j)
 out%criteria(3,j)=2.0_dp*k(j)*n/(n-k(j)-1.0_dp)-2.0_dp*out%loglik(j)
 out%criteria(4,j)=k(j)*log(n)-2.0_dp*out%loglik(j)
 out%criteria(5,j)=2.0_dp*k(j)*log(log(n))-2.0_dp*out%loglik(j)
end do
do j=1,5
 b=minloc(out%criteria(j,:),dim=1)
 out%best(j)=out%fits(b)
end do
end function pearson_msc

function match_moments(mean,variance,skewness,kurtosis,family,skewness_sign,return_distribution) result(res)
real(kind=dp),intent(in)::mean,variance,skewness,kurtosis
integer,intent(in)::family
integer,intent(in),optional::skewness_sign
logical,intent(in),optional::return_distribution
type(pearson_params_t)::res
type(pearson_params_t)::fit
real(kind=dp)::sss,kkk,sf
logical::rd
sss=skewness
kkk=kurtosis
sf=1.0_dp
rd=.false.
if(present(skewness_sign))sf=merge(1.0_dp,-1.0_dp,skewness_sign>=0)
if(present(return_distribution))rd=return_distribution
if(kkk/=kkk)then
 select case(family)
 case(0);kkk=3.0_dp
 case(3);kkk=3.0_dp+1.5_dp*sss*sss
 case(5)
  if(sss*sss<32.0_dp)then
   kkk=-3.0_dp*(16.0_dp+13.0_dp*sss*sss+2.0_dp* &
      sqrt(64.0_dp+48.0_dp*sss*sss+12.0_dp*sss**4+sss**6))/(sss*sss-32.0_dp)
  else
   res%status=1
   res%message="moments and distribution type do not fit"
   return
  end if
 case default
 res%status=1
 res%message="kurtosis must be provided for this type"
 return
 end select
end if
if(sss/=sss)then
 select case(family)
 case(0,2,7);sss=0.0_dp
 case(3);sss=sf*sqrt((kkk-3.0_dp)/1.5_dp)
 case(5)
  ! PearsonDS 1.3.2 references sss while sss is NA in this branch. Preserve
  ! that unsupported source behavior as an explicit status rather than inventing a formula.
  res%status=2
  res%message="PearsonDS type-V missing-skewness branch is undefined in the source"
  return
 case default
 res%status=1
 res%message="skewness must be provided for this type"
 return
 end select
end if
fit=pearson_fit_m(mean,variance,sss,kkk)
if(fit%status/=0.or.fit%family/=family)then
 res=fit
 res%status=3
 res%message="moments and distribution type do not fit"
 return
end if
if(rd)then
 res=fit
else
 res%family=family
 res%npar=4
 res%par=[mean,variance,sss,kkk]
end if
end function match_moments

pure function hypergeom_2f1(a,b,c,z,tol,maxit) result(res)
! Double-complex Gauss hypergeometric series. This is a compact translation of
! the package's internal F21 functionality for |z| comfortably below 1. Public
! Pearson-IV CDF evaluation uses the package's no-GSL numerical-integration path.
complex(kind=dp),intent(in)::a,b,c,z
real(kind=dp),intent(in),optional::tol
integer,intent(in),optional::maxit
complex(kind=dp)::res,term
real(kind=dp)::eps
integer::n,nmax
eps=1.0e-12_dp
if(present(tol))eps=max(tol,epsilon(1.0_dp))
nmax=100000
if(present(maxit))nmax=maxit
res=(1.0_dp,0.0_dp)
term=(1.0_dp,0.0_dp)
do n=1,nmax
 term=term*(a+real(n-1,dp))*(b+real(n-1,dp))*z/((c+real(n-1,dp))*real(n,dp))
 res=res+term
 if(abs(term)<=eps*max(1.0_dp,abs(res)))exit
end do
end function hypergeom_2f1

end module pearsonds_mod
