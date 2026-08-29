module truncnorm_core
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan, &
   ieee_positive_inf, ieee_negative_inf
use r_compat, only: dp, dnorm, normal_cdf, qnorm, runif1, rnorm1
implicit none
private

real(dp), parameter :: t1 = 0.15_dp
real(dp), parameter :: t2 = 2.18_dp
real(dp), parameter :: t3 = 0.725_dp
real(dp), parameter :: t4 = 0.45_dp

public :: dtruncnorm_scalar, ptruncnorm_scalar, qtruncnorm_scalar
public :: etruncnorm_scalar, vtruncnorm_scalar, rtruncnorm_scalar
public :: dtruncnorm_recycle, ptruncnorm_recycle, qtruncnorm_recycle
public :: dtruncnorm_vec, ptruncnorm_vec, qtruncnorm_vec
public :: etruncnorm_recycle, vtruncnorm_recycle, rtruncnorm_recycle
public :: rtruncnorm_n

contains

pure function log1p_stable(x) result(y)
real(dp), intent(in) :: x
real(dp) :: y
if (abs(x) < 1.0e-8_dp) then
   y = x * (1.0_dp + x * (-0.5_dp + x * (1.0_dp / 3.0_dp + x * (-0.25_dp + 0.2_dp * x))))
else
   y = log(1.0_dp + x)
end if
end function log1p_stable

pure function nan_dp() result(x)
real(dp) :: x
x = ieee_value(0.0_dp, ieee_quiet_nan)
end function nan_dp

pure function posinf_dp() result(x)
real(dp) :: x
x = ieee_value(0.0_dp, ieee_positive_inf)
end function posinf_dp

pure function neginf_dp() result(x)
real(dp) :: x
x = ieee_value(0.0_dp, ieee_negative_inf)
end function neginf_dp

pure function logspace_sub(logx, logy) result(z)
! log(exp(logx)-exp(logy)), requiring logx >= logy.
real(dp), intent(in) :: logx, logy
real(dp) :: z
if (logy > logx) then
   z = nan_dp()
else if (logy == neginf_dp()) then
   z = logx
else if (logx == logy) then
   z = neginf_dp()
else
   z = logx + log1p_stable(-exp(logy - logx))
end if
end function logspace_sub

pure function log_normal_cdf(z) result(lp)
! Stable log Phi(z). This supplies the log-tail form missing from r_compat.
real(dp), intent(in) :: z
real(dp) :: lp, zz, invz2, corr
real(dp), parameter :: half_log_2pi = 0.918938533204672741780329736406_dp
if (z /= z) then
   lp = nan_dp()
else if (z == neginf_dp()) then
   lp = neginf_dp()
else if (z == posinf_dp()) then
   lp = 0.0_dp
else if (z < -37.0_dp) then
   zz = -z
   invz2 = 1.0_dp / (zz * zz)
   corr = 1.0_dp - invz2 + 3.0_dp * invz2**2 - 15.0_dp * invz2**3 + 105.0_dp * invz2**4
   lp = -0.5_dp * z * z - log(zz) - half_log_2pi + log(corr)
else if (z < 8.0_dp) then
   lp = log(0.5_dp * erfc(-z / sqrt(2.0_dp)))
else
   lp = log1p_stable(-0.5_dp * erfc(z / sqrt(2.0_dp)))
end if
end function log_normal_cdf

pure function log_normal_sf(z) result(lp)
real(dp), intent(in) :: z
real(dp) :: lp
lp = log_normal_cdf(-z)
end function log_normal_sf

pure function standard_interval_logprob(alpha, beta) result(lp)
! log(Phi(beta)-Phi(alpha)), stable in both tails.
real(dp), intent(in) :: alpha, beta
real(dp) :: lp, la, lb
if (alpha /= alpha .or. beta /= beta .or. beta < alpha) then
   lp = nan_dp()
else if (alpha == beta) then
   lp = neginf_dp()
else if (alpha == neginf_dp() .and. beta == posinf_dp()) then
   lp = 0.0_dp
else if (beta <= 0.0_dp) then
   lb = log_normal_cdf(beta)
   la = log_normal_cdf(alpha)
   lp = logspace_sub(lb, la)
else if (alpha >= 0.0_dp) then
   la = log_normal_sf(alpha)
   lb = log_normal_sf(beta)
   lp = logspace_sub(la, lb)
else
   lp = log(normal_cdf(beta) - normal_cdf(alpha))
end if
end function standard_interval_logprob

pure function standard_interval_prob(alpha, beta) result(p)
real(dp), intent(in) :: alpha, beta
real(dp) :: p, lp
lp = standard_interval_logprob(alpha, beta)
if (lp == neginf_dp()) then
   p = 0.0_dp
else
   p = exp(lp)
end if
end function standard_interval_prob

pure function ptruncnorm_stable(q, a, b, mean, sd) result(p)
real(dp), intent(in) :: q, a, b, mean, sd
real(dp) :: p, alpha, beta, z, lnumer, ldenom
if (sd <= 0.0_dp .or. sd /= sd .or. a > b) then
   p = nan_dp()
   return
end if
if (q < a) then
   p = 0.0_dp
   return
else if (q > b) then
   p = 1.0_dp
   return
end if
if (a == b) then
   if (q < a) then
      p = 0.0_dp
   else
      p = 1.0_dp
   end if
   return
end if
alpha = (a - mean) / sd
beta = (b - mean) / sd
z = (q - mean) / sd
ldenom = standard_interval_logprob(alpha, beta)
lnumer = standard_interval_logprob(alpha, min(z, beta))
if (.not. ieee_is_finite(ldenom)) then
   ! Match the upstream finite-precision behavior in extremely remote,
   ! ultra-narrow intervals by falling back to a uniform approximation.
   p = (q - a) / (b - a)
else if (lnumer == neginf_dp()) then
   p = 0.0_dp
else
   p = exp(lnumer - ldenom)
   p = max(0.0_dp, min(1.0_dp, p))
end if
end function ptruncnorm_stable

pure function dtruncnorm_scalar(x, a, b, mean, sd) result(y)
real(dp), intent(in) :: x
real(dp), intent(in), optional :: a, b, mean, sd
real(dp) :: y, aa, bb, mu, sig, alpha, beta, logz, logphi, rawz

aa = neginf_dp()
bb = posinf_dp()
mu = 0.0_dp
sig = 1.0_dp
if (present(a)) aa = a
if (present(b)) bb = b
if (present(mean)) mu = mean
if (present(sd)) sig = sd
if (sig <= 0.0_dp .or. sig /= sig .or. aa > bb) then
   y = nan_dp()
   return
end if
if (x < aa .or. x > bb) then
   y = 0.0_dp
   return
end if
if (aa == bb) then
   y = posinf_dp()
   return
end if
alpha = (aa - mu) / sig
beta = (bb - mu) / sig
rawz = normal_cdf(beta) - normal_cdf(alpha)
logz = standard_interval_logprob(alpha, beta)
logphi = dnorm((x - mu) / sig, 0.0_dp, 1.0_dp, log_=.true.)
if (rawz <= 0.0_dp .or. .not. ieee_is_finite(rawz) .or. .not. ieee_is_finite(logz)) then
   ! Upstream falls back to a uniform density when finite-precision pnorm
   ! cannot resolve the truncation mass. Keep that documented behavior.
   y = 1.0_dp / (bb - aa)
else
   y = exp(logphi - log(sig) - logz)
end if
end function dtruncnorm_scalar

pure function ptruncnorm_scalar(q, a, b, mean, sd) result(p)
real(dp), intent(in) :: q
real(dp), intent(in), optional :: a, b, mean, sd
real(dp) :: p, aa, bb, mu, sig
aa = neginf_dp()
bb = posinf_dp()
mu = 0.0_dp
sig = 1.0_dp
if (present(a)) aa = a
if (present(b)) bb = b
if (present(mean)) mu = mean
if (present(sd)) sig = sd
p = ptruncnorm_stable(q, aa, bb, mu, sig)
end function ptruncnorm_scalar

function qtruncnorm_scalar(p, a, b, mean, sd) result(q)
real(dp), intent(in) :: p
real(dp), intent(in), optional :: a, b, mean, sd
real(dp) :: q, aa, bb, mu, sig, lo, hi, flo, fhi
integer :: iter

aa = neginf_dp()
bb = posinf_dp()
mu = 0.0_dp
sig = 1.0_dp
if (present(a)) aa = a
if (present(b)) bb = b
if (present(mean)) mu = mean
if (present(sd)) sig = sd
if (sig <= 0.0_dp .or. sig /= sig .or. aa > bb .or. p /= p .or. p < 0.0_dp .or. p > 1.0_dp) then
   q = nan_dp()
   return
end if
if (p == 0.0_dp) then
   q = aa
   return
else if (p == 1.0_dp) then
   q = bb
   return
else if (aa == neginf_dp() .and. bb == posinf_dp()) then
   q = qnorm(p, mean=mu, sd=sig)
   return
end if

lo = aa
hi = bb
if (lo == neginf_dp()) then
   lo = min(-1.0_dp, mu - sig)
   do
      flo = ptruncnorm_stable(lo, aa, bb, mu, sig) - p
      if (flo < 0.0_dp) exit
      lo = 2.0_dp * lo
      if (.not. ieee_is_finite(lo)) then
         lo = mu - 40.0_dp * sig
         exit
      end if
   end do
end if
if (hi == posinf_dp()) then
   hi = max(1.0_dp, mu + sig)
   do
      fhi = ptruncnorm_stable(hi, aa, bb, mu, sig) - p
      if (fhi > 0.0_dp) exit
      hi = 2.0_dp * hi
      if (.not. ieee_is_finite(hi)) then
         hi = mu + 40.0_dp * sig
         exit
      end if
   end do
end if
q = brent_quantile(lo, hi, aa, bb, mu, sig, p)
end function qtruncnorm_scalar

function brent_quantile(ax, bx, alower, bupper, mean, sd, target) result(root)
! Specialized port of upstream truncnorm_zeroin(), itself derived from R's
! Brent/zeroin implementation and NETLIB c/brent.shar.
real(dp), intent(in) :: ax, bx, alower, bupper, mean, sd, target
real(dp) :: root
real(dp) :: a, b, c, fa, fb, fc, prev_step, tol_act, pp, qq, new_step
real(dp) :: t1i, cb, t2i
integer :: iter

a = ax
b = bx
c = a
fa = ptruncnorm_stable(a, alower, bupper, mean, sd) - target
fb = ptruncnorm_stable(b, alower, bupper, mean, sd) - target
fc = fa
if (fa == 0.0_dp) then
   root = a
   return
end if
if (fb == 0.0_dp) then
   root = b
   return
end if

do iter = 1, 201
   prev_step = b - a
   if (abs(fc) < abs(fb)) then
      ! Keep the sequential assignments of the original C routine.
      a = b
      b = c
      c = a
      fa = fb
      fb = fc
      fc = fa
   end if
   tol_act = 2.0_dp * epsilon(1.0_dp) * abs(b)
   new_step = 0.5_dp * (c - b)
   if (abs(new_step) <= tol_act .or. fb == 0.0_dp) then
      root = b
      return
   end if
   if (abs(prev_step) >= tol_act .and. abs(fa) > abs(fb)) then
      cb = c - b
      if (a == c) then
         t1i = fb / fa
         pp = cb * t1i
         qq = 1.0_dp - t1i
      else
         qq = fa / fc
         t1i = fb / fc
         t2i = fb / fa
         pp = t2i * (cb * qq * (qq - t1i) - (b - a) * (t1i - 1.0_dp))
         qq = (qq - 1.0_dp) * (t1i - 1.0_dp) * (t2i - 1.0_dp)
      end if
      if (pp > 0.0_dp) then
         qq = -qq
      else
         pp = -pp
      end if
      if (pp < 0.75_dp * cb * qq - abs(tol_act * qq) / 2.0_dp .and. &
          pp < abs(prev_step * qq / 2.0_dp)) new_step = pp / qq
   end if
   if (abs(new_step) < tol_act) then
      if (new_step > 0.0_dp) then
         new_step = tol_act
      else
         new_step = -tol_act
      end if
   end if
   a = b
   fa = fb
   b = b + new_step
   fb = ptruncnorm_stable(b, alower, bupper, mean, sd) - target
   if ((fb > 0.0_dp .and. fc > 0.0_dp) .or. (fb < 0.0_dp .and. fc < 0.0_dp)) then
      c = a
      fc = fa
   end if
end do
root = b
end function brent_quantile

pure function e_lefttruncnorm(a, mean, sd) result(e)
real(dp), intent(in) :: a, mean, sd
real(dp) :: e, alpha, logphi, logtail
alpha = (a - mean) / sd
logphi = dnorm(alpha, 0.0_dp, 1.0_dp, log_=.true.)
logtail = log_normal_sf(alpha)
e = mean + sd * exp(logphi - logtail)
end function e_lefttruncnorm

pure function e_righttruncnorm(b, mean, sd) result(e)
real(dp), intent(in) :: b, mean, sd
real(dp) :: e, beta, logphi, logcdf
beta = (b - mean) / sd
logphi = dnorm(beta, 0.0_dp, 1.0_dp, log_=.true.)
logcdf = log_normal_cdf(beta)
e = mean - sd * exp(logphi - logcdf)
end function e_righttruncnorm

pure function e_twosided(a, b, mean, sd) result(e)
real(dp), intent(in) :: a, b, mean, sd
real(dp) :: e, alpha, beta, lpa, lpb, lz, logphi_a, logphi_b, ratio
if (b < mean - 6.0_dp * sd .or. a > mean + 6.0_dp * sd) then
   e = 0.5_dp * (a + b)
   return
end if
alpha = (a - mean) / sd
beta = (b - mean) / sd
logphi_a = dnorm(alpha, 0.0_dp, 1.0_dp, log_=.true.)
logphi_b = dnorm(beta, 0.0_dp, 1.0_dp, log_=.true.)
lz = standard_interval_logprob(alpha, beta)
if (logphi_a >= logphi_b) then
   ratio = exp(logspace_sub(logphi_a, logphi_b) - lz)
else
   ratio = -exp(logspace_sub(logphi_b, logphi_a) - lz)
end if
e = mean + sd * ratio
end function e_twosided

pure function etruncnorm_scalar(a, b, mean, sd) result(e)
real(dp), intent(in), optional :: a, b, mean, sd
real(dp) :: e, aa, bb, mu, sig
aa = neginf_dp()
bb = posinf_dp()
mu = 0.0_dp
sig = 1.0_dp
if (present(a)) aa = a
if (present(b)) bb = b
if (present(mean)) mu = mean
if (present(sd)) sig = sd
if (sig <= 0.0_dp .or. sig /= sig .or. aa > bb) then
   e = nan_dp()
else if (ieee_is_finite(aa) .and. ieee_is_finite(bb)) then
   e = e_twosided(aa, bb, mu, sig)
else if (aa == neginf_dp() .and. ieee_is_finite(bb)) then
   e = e_righttruncnorm(bb, mu, sig)
else if (ieee_is_finite(aa) .and. bb == posinf_dp()) then
   e = e_lefttruncnorm(aa, mu, sig)
else if (aa == neginf_dp() .and. bb == posinf_dp()) then
   e = mu
else
   e = nan_dp()
end if
end function etruncnorm_scalar

pure function v_lefttruncnorm(a, mean, sd) result(v)
real(dp), intent(in) :: a, mean, sd
real(dp) :: v, alpha, lambda
alpha = (a - mean) / sd
lambda = exp(dnorm(alpha, 0.0_dp, 1.0_dp, log_=.true.) - log_normal_sf(alpha))
v = sd * sd * (1.0_dp - lambda * (lambda - alpha))
end function v_lefttruncnorm

pure function v_righttruncnorm(b, mean, sd) result(v)
real(dp), intent(in) :: b, mean, sd
real(dp) :: v
v = v_lefttruncnorm(-b, -mean, sd)
end function v_righttruncnorm

pure function v_twosided(a, b, mean, sd) result(vout)
real(dp), intent(in) :: a, b, mean, sd
real(dp) :: vout, vv, pi1, pi2, pi3, e1, e2, e3, v1, v3, c1, c3
if (b < mean - 6.0_dp * sd .or. a > mean + 6.0_dp * sd) then
   vout = (b - a)**2 / 12.0_dp
   return
end if
vv = sd * sd
pi1 = normal_cdf((a - mean) / sd)
pi2 = normal_cdf((b - mean) / sd) - pi1
pi3 = 1.0_dp - normal_cdf((b - mean) / sd)
e1 = e_righttruncnorm(a, mean, sd)
e2 = e_twosided(a, b, mean, sd)
e3 = e_lefttruncnorm(b, mean, sd)
v1 = v_righttruncnorm(a, mean, sd)
v3 = v_lefttruncnorm(b, mean, sd)
c1 = pi1 * (v1 + (e1 - mean)**2)
c3 = pi3 * (v3 + (e3 - mean)**2)
vout = (vv - c1 - c3) / pi2 - (e2 - mean)**2
end function v_twosided

pure function vtruncnorm_scalar(a, b, mean, sd) result(v)
real(dp), intent(in), optional :: a, b, mean, sd
real(dp) :: v, aa, bb, mu, sig
aa = neginf_dp()
bb = posinf_dp()
mu = 0.0_dp
sig = 1.0_dp
if (present(a)) aa = a
if (present(b)) bb = b
if (present(mean)) mu = mean
if (present(sd)) sig = sd
if (sig <= 0.0_dp .or. sig /= sig .or. aa > bb) then
   v = nan_dp()
else if (ieee_is_finite(aa) .and. ieee_is_finite(bb)) then
   v = v_twosided(aa, bb, mu, sig)
else if (aa == neginf_dp() .and. ieee_is_finite(bb)) then
   v = v_righttruncnorm(bb, mu, sig)
else if (ieee_is_finite(aa) .and. bb == posinf_dp()) then
   v = v_lefttruncnorm(aa, mu, sig)
else if (aa == neginf_dp() .and. bb == posinf_dp()) then
   v = sig * sig
else
   v = nan_dp()
end if
end function vtruncnorm_scalar

function exponential_draw(rate) result(x)
real(dp), intent(in) :: rate
real(dp) :: x
x = -log(max(tiny(1.0_dp), 1.0_dp - runif1())) / rate
end function exponential_draw

function ers_a_inf(a) result(x)
real(dp), intent(in) :: a
real(dp) :: x, rho
do
   x = exponential_draw(a) + a
   rho = exp(-0.5_dp * (x - a)**2)
   if (runif1() <= rho) exit
end do
end function ers_a_inf

function ers_a_b(a, b) result(x)
real(dp), intent(in) :: a, b
real(dp) :: x, rho
do
   x = exponential_draw(a) + a
   rho = exp(-0.5_dp * (x - a)**2)
   if (runif1() <= rho .and. x <= b) exit
end do
end function ers_a_b

function nrs_a_b(a, b) result(x)
real(dp), intent(in) :: a, b
real(dp) :: x
do
   x = rnorm1()
   if (x >= a .and. x <= b) exit
end do
end function nrs_a_b

function nrs_a_inf(a) result(x)
real(dp), intent(in) :: a
real(dp) :: x
do
   x = rnorm1()
   if (x >= a) exit
end do
end function nrs_a_inf

function hnrs_a_b(a, b) result(x)
real(dp), intent(in) :: a, b
real(dp) :: x
do
   x = abs(rnorm1())
   if (x >= a .and. x <= b) exit
end do
end function hnrs_a_b

function urs_a_b(a, b) result(x)
real(dp), intent(in) :: a, b
real(dp) :: x, u, phi_a, ub
real(dp), parameter :: inv_sqrt_2pi = 0.398942280401432677939946059934_dp
phi_a = dnorm(a, 0.0_dp, 1.0_dp)
if (a < 0.0_dp .and. b > 0.0_dp) then
   ub = inv_sqrt_2pi
else
   ub = phi_a
end if
do
   x = a + (b - a) * runif1()
   u = runif1()
   if (u * ub <= dnorm(x, 0.0_dp, 1.0_dp)) exit
end do
end function urs_a_b

function r_lefttruncnorm(a, mean, sd) result(x)
real(dp), intent(in) :: a, mean, sd
real(dp) :: x, alpha
alpha = (a - mean) / sd
if (alpha < t4) then
   x = mean + sd * nrs_a_inf(alpha)
else
   x = mean + sd * ers_a_inf(alpha)
end if
end function r_lefttruncnorm

function r_righttruncnorm(b, mean, sd) result(x)
real(dp), intent(in) :: b, mean, sd
real(dp) :: x, beta
beta = (b - mean) / sd
x = mean - sd * r_lefttruncnorm(-beta, 0.0_dp, 1.0_dp)
end function r_righttruncnorm

function r_twosided(a, b, mean, sd) result(x)
real(dp), intent(in) :: a, b, mean, sd
real(dp) :: x, alpha, beta, phi_a, phi_b
alpha = (a - mean) / sd
beta = (b - mean) / sd
phi_a = dnorm(alpha, 0.0_dp, 1.0_dp)
phi_b = dnorm(beta, 0.0_dp, 1.0_dp)
if (beta <= alpha) then
   x = nan_dp()
else if (alpha <= 0.0_dp .and. beta >= 0.0_dp) then
   if (phi_a <= t1 .or. phi_b <= t1) then
      x = mean + sd * nrs_a_b(alpha, beta)
   else
      x = mean + sd * urs_a_b(alpha, beta)
   end if
else if (alpha > 0.0_dp) then
   if (phi_a / phi_b <= t2) then
      x = mean + sd * urs_a_b(alpha, beta)
   else if (alpha < t3) then
      x = mean + sd * hnrs_a_b(alpha, beta)
   else
      x = mean + sd * ers_a_b(alpha, beta)
   end if
else
   if (phi_b / phi_a <= t2) then
      x = mean - sd * urs_a_b(-beta, -alpha)
   else if (beta > -t3) then
      x = mean - sd * hnrs_a_b(-beta, -alpha)
   else
      x = mean - sd * ers_a_b(-beta, -alpha)
   end if
end if
end function r_twosided

function rtruncnorm_scalar(a, b, mean, sd) result(x)
real(dp), intent(in), optional :: a, b, mean, sd
real(dp) :: x, aa, bb, mu, sig
aa = neginf_dp()
bb = posinf_dp()
mu = 0.0_dp
sig = 1.0_dp
if (present(a)) aa = a
if (present(b)) bb = b
if (present(mean)) mu = mean
if (present(sd)) sig = sd
if (sig <= 0.0_dp .or. sig /= sig .or. aa > bb) then
   x = nan_dp()
else if (ieee_is_finite(aa) .and. ieee_is_finite(bb)) then
   x = r_twosided(aa, bb, mu, sig)
else if (aa == neginf_dp() .and. ieee_is_finite(bb)) then
   x = r_righttruncnorm(bb, mu, sig)
else if (ieee_is_finite(aa) .and. bb == posinf_dp()) then
   x = r_lefttruncnorm(aa, mu, sig)
else if (aa == neginf_dp() .and. bb == posinf_dp()) then
   x = mu + sig * rnorm1()
else
   x = nan_dp()
end if
end function rtruncnorm_scalar

function dtruncnorm_vec(x, a, b, mean, sd) result(y)
real(dp), intent(in) :: x(:)
real(dp), intent(in), optional :: a, b, mean, sd
real(dp), allocatable :: y(:)
integer :: i
allocate(y(size(x)))
do i = 1, size(x)
   y(i) = dtruncnorm_scalar(x(i), a, b, mean, sd)
end do
end function dtruncnorm_vec

function ptruncnorm_vec(q, a, b, mean, sd) result(y)
real(dp), intent(in) :: q(:)
real(dp), intent(in), optional :: a, b, mean, sd
real(dp), allocatable :: y(:)
integer :: i
allocate(y(size(q)))
do i = 1, size(q)
   y(i) = ptruncnorm_scalar(q(i), a, b, mean, sd)
end do
end function ptruncnorm_vec

function qtruncnorm_vec(p, a, b, mean, sd) result(y)
real(dp), intent(in) :: p(:)
real(dp), intent(in), optional :: a, b, mean, sd
real(dp), allocatable :: y(:)
integer :: i
allocate(y(size(p)))
do i = 1, size(p)
   y(i) = qtruncnorm_scalar(p(i), a, b, mean, sd)
end do
end function qtruncnorm_vec

function rtruncnorm_n(n, a, b, mean, sd) result(y)
integer, intent(in) :: n
real(dp), intent(in), optional :: a, b, mean, sd
real(dp), allocatable :: y(:)
integer :: i
allocate(y(max(0,n)))
do i = 1, size(y)
   y(i) = rtruncnorm_scalar(a, b, mean, sd)
end do
end function rtruncnorm_n

function dtruncnorm_recycle(x, a, b, mean, sd) result(y)
real(dp), intent(in) :: x(:), a(:), b(:), mean(:), sd(:)
real(dp), allocatable :: y(:)
integer :: n, i
if (min(size(x), size(a), size(b), size(mean), size(sd)) == 0) then
   allocate(y(0))
   return
end if
n = max(size(x), size(a), size(b), size(mean), size(sd))
allocate(y(n))
do i = 1, n
   y(i) = dtruncnorm_scalar(x(1+mod(i-1,size(x))), a(1+mod(i-1,size(a))), &
      b(1+mod(i-1,size(b))), mean(1+mod(i-1,size(mean))), sd(1+mod(i-1,size(sd))))
end do
end function dtruncnorm_recycle

function ptruncnorm_recycle(q, a, b, mean, sd) result(y)
real(dp), intent(in) :: q(:), a(:), b(:), mean(:), sd(:)
real(dp), allocatable :: y(:)
integer :: n, i
if (min(size(q), size(a), size(b), size(mean), size(sd)) == 0) then
   allocate(y(0))
   return
end if
n = max(size(q), size(a), size(b), size(mean), size(sd))
allocate(y(n))
do i = 1, n
   y(i) = ptruncnorm_scalar(q(1+mod(i-1,size(q))), a(1+mod(i-1,size(a))), &
      b(1+mod(i-1,size(b))), mean(1+mod(i-1,size(mean))), sd(1+mod(i-1,size(sd))))
end do
end function ptruncnorm_recycle

function qtruncnorm_recycle(p, a, b, mean, sd) result(y)
real(dp), intent(in) :: p(:), a(:), b(:), mean(:), sd(:)
real(dp), allocatable :: y(:)
integer :: n, i
if (min(size(p), size(a), size(b), size(mean), size(sd)) == 0) then
   allocate(y(0))
   return
end if
n = max(size(p), size(a), size(b), size(mean), size(sd))
allocate(y(n))
do i = 1, n
   y(i) = qtruncnorm_scalar(p(1+mod(i-1,size(p))), a(1+mod(i-1,size(a))), &
      b(1+mod(i-1,size(b))), mean(1+mod(i-1,size(mean))), sd(1+mod(i-1,size(sd))))
end do
end function qtruncnorm_recycle

function etruncnorm_recycle(a, b, mean, sd) result(y)
real(dp), intent(in) :: a(:), b(:), mean(:), sd(:)
real(dp), allocatable :: y(:)
integer :: n, i
if (min(size(a), size(b), size(mean), size(sd)) == 0) then
   allocate(y(0))
   return
end if
n = max(size(a), size(b), size(mean), size(sd))
allocate(y(n))
do i = 1, n
   y(i) = etruncnorm_scalar(a(1+mod(i-1,size(a))), b(1+mod(i-1,size(b))), &
      mean(1+mod(i-1,size(mean))), sd(1+mod(i-1,size(sd))))
end do
end function etruncnorm_recycle

function vtruncnorm_recycle(a, b, mean, sd) result(y)
real(dp), intent(in) :: a(:), b(:), mean(:), sd(:)
real(dp), allocatable :: y(:)
integer :: n, i
if (min(size(a), size(b), size(mean), size(sd)) == 0) then
   allocate(y(0))
   return
end if
n = max(size(a), size(b), size(mean), size(sd))
allocate(y(n))
do i = 1, n
   y(i) = vtruncnorm_scalar(a(1+mod(i-1,size(a))), b(1+mod(i-1,size(b))), &
      mean(1+mod(i-1,size(mean))), sd(1+mod(i-1,size(sd))))
end do
end function vtruncnorm_recycle

function rtruncnorm_recycle(n, a, b, mean, sd) result(y)
integer, intent(in) :: n
real(dp), intent(in) :: a(:), b(:), mean(:), sd(:)
real(dp), allocatable :: y(:)
integer :: nn, i
if (min(size(a), size(b), size(mean), size(sd)) == 0 .or. n <= 0) then
   allocate(y(0))
   return
end if
nn = max(n, size(a), size(b), size(mean), size(sd))
allocate(y(nn))
do i = 1, nn
   y(i) = rtruncnorm_scalar(a(1+mod(i-1,size(a))), b(1+mod(i-1,size(b))), &
      mean(1+mod(i-1,size(mean))), sd(1+mod(i-1,size(sd))))
end do
end function rtruncnorm_recycle

end module truncnorm_core
