! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from / supporting R package tweedie 3.1.0 by Peter K. Dunn.
module tweedie_likelihood_mod
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
use r_compat, only: dp, r_digamma, solve_real, dnorm, dgamma, dpois
use tweedie_distribution_mod, only: dtweedie, dtweedie_saddle, tweedie_dev
use tweedie_series_mod, only: dtweedie_series, dtweedie_logw_smallp, &
   dtweedie_jw_smallp, dtweedie_logv_bigp, dtweedie_kv_bigp
use tweedie_inversion_mod, only: dtweedie_inversion
implicit none
private
integer, parameter, public :: tweedie_method_auto = 0
integer, parameter, public :: tweedie_method_series = 1
integer, parameter, public :: tweedie_method_inversion = 2
integer, parameter, public :: tweedie_method_saddle = 3
integer, parameter, public :: tweedie_phi_mle_method = 0
integer, parameter, public :: tweedie_phi_saddle_method = 1
public :: tweedie_loglik, tweedie_aic, dtweedie_dlogfdphi, dtweedie_dldphi
public :: tweedie_phi_mle, tweedie_glm_fit, tweedie_profile_grid

type, public :: tweedie_profile_result
   real(dp), allocatable :: power(:)
   real(dp), allocatable :: phi(:)
   real(dp), allocatable :: loglik(:)
   real(dp), allocatable :: beta(:,:)
   logical, allocatable :: converged(:)
end type tweedie_profile_result

contains

function density_for_method(y, mu, phi, power, method) result(d)
real(dp), intent(in) :: y, mu, phi, power
integer, intent(in) :: method
real(dp) :: d
select case(method)
case(tweedie_method_series)
   d = dtweedie_series(y, power, mu, phi)
case(tweedie_method_inversion)
   d = dtweedie_inversion(y, mu, phi, power)
case(tweedie_method_saddle)
   d = dtweedie_saddle(y, mu, phi, power)
case default
   d = dtweedie(y, mu, phi, power)
end select
end function density_for_method

function tweedie_loglik(y, mu, phi, power, method, weights) result(loglik)
real(dp), intent(in) :: y(:), mu(:), phi, power
integer, intent(in), optional :: method
real(dp), intent(in), optional :: weights(:)
real(dp) :: loglik, d, w
integer :: i, im
im = tweedie_method_auto
if (present(method)) im = method
if (size(mu) /= size(y) .or. phi <= 0.0_dp) then
   loglik = ieee_value(0.0_dp, ieee_quiet_nan)
   return
end if
if (present(weights)) then
   if (size(weights) /= size(y)) then
      loglik = ieee_value(0.0_dp, ieee_quiet_nan)
      return
   end if
end if
loglik = 0.0_dp
do i = 1, size(y)
   w = 1.0_dp
   if (present(weights)) w = weights(i)
   if (w == 0.0_dp) cycle
   if (power == 0.0_dp) then
      d = dnorm(y(i), mean=mu(i), sd=sqrt(phi))
   else
      d = density_for_method(y(i), mu(i), phi, power, im)
   end if
   if (.not. ieee_is_finite(d) .or. d <= 0.0_dp) then
      loglik = -huge(1.0_dp)
      return
   end if
   loglik = loglik + w * log(d)
end do
end function tweedie_loglik

pure elemental function tweedie_aic(loglik, npar, k) result(aic)
real(dp), intent(in) :: loglik
integer, intent(in) :: npar
real(dp), intent(in), optional :: k
real(dp) :: aic, penalty
penalty = 2.0_dp
if (present(k)) penalty = k
aic = -2.0_dp * loglik + penalty * real(npar, dp)
end function tweedie_aic

function dtweedie_dlogfdphi(y, mu, phi, power) result(f)
real(dp), intent(in) :: y, mu, phi, power
real(dp) :: f, a, aa, bb, kv, dv, logv, jw, dw, logw
real(dp) :: delta, f1, f2
if (phi <= 0.0_dp .or. mu <= 0.0_dp) then
   f = ieee_value(0.0_dp, ieee_quiet_nan)
   return
end if
if (power == 0.0_dp) then
   f = -0.5_dp / phi + (y - mu)**2 / (2.0_dp * phi**2)
   return
end if
if (power == 1.0_dp) then
   f = mu - y - y * log(mu / phi) + y * r_digamma(1.0_dp + y / phi)
   f = f / phi**2
   return
end if
if (power == 2.0_dp) then
   if (y <= 0.0_dp) then
      f = ieee_value(0.0_dp, ieee_quiet_nan)
   else
      f = -log(y) + y / mu + r_digamma(1.0_dp / phi) - 1.0_dp + log(mu * phi)
      f = f / phi**2
   end if
   return
end if
if (y < 0.0_dp) then
   f = ieee_value(0.0_dp, ieee_quiet_nan)
   return
end if

a = (2.0_dp - power) / (1.0_dp - power)
aa = y * mu**(1.0_dp-power) / (phi**2 * (power-1.0_dp))
bb = mu**(2.0_dp-power) / (phi**2 * (2.0_dp-power))

if (power > 2.0_dp) then
   if (y <= 0.0_dp) then
      f = ieee_value(0.0_dp, ieee_quiet_nan)
      return
   end if
   kv = dtweedie_kv_bigp(y, phi, power)
   dv = kv * (a - 1.0_dp) / phi
   logv = dtweedie_logv_bigp(y, phi, power)
   if (.not. ieee_is_finite(logv) .or. y < 1.0_dp) then
      delta = max(1.0e-7_dp, 1.0e-5_dp * phi)
      f1 = dtweedie(y, mu, phi, power)
      f2 = dtweedie(y, mu, phi + delta, power)
      if (f1 <= 0.0_dp .or. f2 <= 0.0_dp) then
         f = ieee_value(0.0_dp, ieee_quiet_nan)
      else
         f = (log(f2) - log(f1)) / delta
      end if
   else
      f = aa + bb + dv / exp(logv)
   end if
else
   if (y == 0.0_dp) then
      f = mu**(2.0_dp-power) / (phi**2 * (2.0_dp-power))
   else
      jw = dtweedie_jw_smallp(y, phi, power)
      dw = jw * (a - 1.0_dp) / phi
      logw = dtweedie_logw_smallp(y, phi, power)
      f = aa + bb + dw / exp(logw)
   end if
end if
end function dtweedie_dlogfdphi

function dtweedie_dldphi(phi, mu, power, y, weights) result(score)
real(dp), intent(in) :: phi, mu(:), power, y(:)
real(dp), intent(in), optional :: weights(:)
real(dp) :: score, f, w, scale
integer :: i
if (size(mu) /= size(y) .or. phi <= 0.0_dp) then
   score = ieee_value(0.0_dp, ieee_quiet_nan)
   return
end if
if (present(weights)) then
   if (size(weights) /= size(y)) then
      score = ieee_value(0.0_dp, ieee_quiet_nan)
      return
   end if
end if
score = 0.0_dp
! Preserve the stabilization used by the R package when phi^(1/(p-2)) is in (0,1).
if (power /= 1.0_dp .and. power /= 2.0_dp .and. power /= 0.0_dp) then
   scale = phi**(1.0_dp/(power-2.0_dp))
else
   scale = 1.0_dp
end if
do i = 1, size(y)
   w = 1.0_dp
   if (present(weights)) w = weights(i)
   if (w == 0.0_dp) cycle
   if (power /= 0.0_dp .and. power /= 1.0_dp .and. power /= 2.0_dp .and. &
       scale > 0.0_dp .and. scale < 1.0_dp) then
      f = dtweedie_dlogfdphi(scale*y(i), scale*mu(i), 1.0_dp, power)
      f = f * scale**(2.0_dp-power)
   else
      f = dtweedie_dlogfdphi(y(i), mu(i), phi, power)
   end if
   score = score - 2.0_dp * w * f
end do
end function dtweedie_dldphi

function tweedie_phi_mle(y, mu, power, method, weights, phi_start, lower, upper, &
   tol, maxit) result(phi_hat)
real(dp), intent(in) :: y(:), mu(:), power
integer, intent(in), optional :: method, maxit
real(dp), intent(in), optional :: weights(:), phi_start, lower, upper, tol
real(dp) :: phi_hat, phi0, lo, hi, x1, x2, f1, f2, gr, eps
integer :: im, it, mit, expand
if (size(mu) /= size(y) .or. size(y) == 0) then
   phi_hat = ieee_value(0.0_dp, ieee_quiet_nan)
   return
end if
if (present(weights)) then
   if (size(weights) /= size(y)) then
      phi_hat = ieee_value(0.0_dp, ieee_quiet_nan)
      return
   end if
end if
im = tweedie_method_auto
if (present(method)) im = method
phi0 = sum(tweedie_dev(y, mu, power)) / real(size(y), dp)
if (present(phi_start)) phi0 = phi_start
if (.not. ieee_is_finite(phi0) .or. phi0 <= 0.0_dp) phi0 = 1.0_dp
lo = min(1.0e-3_dp, 0.5_dp*phi0)
hi = 10.0_dp*phi0
if (present(lower)) lo = lower
if (present(upper)) hi = upper
lo = max(lo, tiny(1.0_dp))
if (hi <= lo) hi = 10.0_dp*lo
! Expand a user-unspecified upper bound if the likelihood is still increasing there.
if (.not. present(upper)) then
   do expand = 1, 8
      if (objective(hi) >= objective(hi/2.0_dp)) exit
      hi = 2.0_dp*hi
   end do
end if
eps = 1.0e-8_dp
if (present(tol)) eps = tol
mit = 100
if (present(maxit)) mit = maxit
gr = (sqrt(5.0_dp)-1.0_dp)/2.0_dp
x1 = hi - gr*(hi-lo)
x2 = lo + gr*(hi-lo)
f1 = objective(x1)
f2 = objective(x2)
do it = 1, mit
   if (abs(hi-lo) <= eps*max(1.0_dp,abs(x1)+abs(x2))) exit
   if (f1 > f2) then
      lo=x1
      x1=x2
      f1=f2
      x2=lo+gr*(hi-lo)
      f2=objective(x2)
   else
      hi=x2
      x2=x1
      f2=f1
      x1=hi-gr*(hi-lo)
      f1=objective(x1)
   end if
end do
phi_hat = 0.5_dp*(lo+hi)
contains
function objective(phi) result(v)
real(dp), intent(in) :: phi
real(dp) :: v, ll
ll = tweedie_loglik(y, mu, phi, power, im, weights)
if (.not. ieee_is_finite(ll)) then
   v = huge(1.0_dp)
else
   v = -2.0_dp*ll
end if
end function objective
end function tweedie_phi_mle

subroutine tweedie_glm_fit(x, y, power, link_power, beta, mu, converged, iterations, &
   weights, offset, tol, maxit)
real(dp), intent(in) :: x(:,:), y(:), power, link_power
real(dp), allocatable, intent(out) :: beta(:), mu(:)
logical, intent(out) :: converged
integer, intent(out) :: iterations
real(dp), intent(in), optional :: weights(:), offset(:), tol
integer, intent(in), optional :: maxit
real(dp), allocatable :: eta(:), z(:), ww(:), target(:), xtwx(:,:), xtwz(:), bnew(:)
real(dp) :: wt, off, dm, eps, ybar, delta
integer :: n, k, i, j, l, mit
n=size(y)
k=size(x,2)
allocate(beta(k),mu(n),eta(n),z(n),ww(n),target(n),xtwx(k,k),xtwz(k),bnew(k))
beta=0.0_dp
mu=0.0_dp
converged=.false.
iterations=0
if (size(x,1)/=n .or. n==0 .or. k==0) return
if (present(weights)) then
   if (size(weights)/=n) return
end if
if (present(offset)) then
   if (size(offset)/=n) return
end if
ybar=max(sum(max(y,0.0_dp))/real(n,dp),1.0e-6_dp)
do i=1,n
   if(link_power==0.0_dp)then
      target(i)=log(max(y(i),0.1_dp*ybar))
   else
      target(i)=max(y(i),0.1_dp*ybar)**link_power
   end if
   if(present(offset))target(i)=target(i)-offset(i)
end do
xtwx=matmul(transpose(x),x)
xtwz=matmul(transpose(x),target)
beta=solve_real(xtwx,xtwz)
eps=1.0e-9_dp
if(present(tol))eps=tol
mit=25
if(present(maxit))mit=maxit
do iterations=1,mit
   eta=matmul(x,beta)
   if(present(offset))eta=eta+offset
   do i=1,n
      if(link_power==0.0_dp)then
         mu(i)=exp(min(eta(i),700.0_dp))
         dm=mu(i)
      else
         if(eta(i)<=0.0_dp)then
            mu(i)=max(ybar,1.0e-8_dp)
            dm=mu(i)**(1.0_dp-link_power)/link_power
         else
            mu(i)=eta(i)**(1.0_dp/link_power)
            dm=mu(i)**(1.0_dp-link_power)/link_power
         end if
      end if
      mu(i)=max(mu(i),tiny(1.0_dp))
      if(abs(dm)<=tiny(1.0_dp))dm=sign(tiny(1.0_dp),dm+tiny(1.0_dp))
      z(i)=eta(i)+(y(i)-mu(i))/dm
      off=0.0_dp
      if(present(offset))off=offset(i)
      target(i)=z(i)-off
      wt=1.0_dp
      if(present(weights))wt=weights(i)
      ww(i)=max(0.0_dp,wt)*dm*dm/(mu(i)**power)
   end do
   xtwx=0.0_dp
   xtwz=0.0_dp
   do i=1,n
      do j=1,k
         xtwz(j)=xtwz(j)+ww(i)*x(i,j)*target(i)
         do l=1,k
            xtwx(j,l)=xtwx(j,l)+ww(i)*x(i,j)*x(i,l)
         end do
      end do
   end do
   bnew=solve_real(xtwx,xtwz)
   delta=maxval(abs(bnew-beta))
   beta=bnew
   if(delta<=eps*(1.0_dp+maxval(abs(beta))))then
      converged=.true.
      exit
   end if
end do
eta=matmul(x,beta)
if(present(offset))eta=eta+offset
do i=1,n
   if(link_power==0.0_dp)then
      mu(i)=exp(min(eta(i),700.0_dp))
   else if(eta(i)>0.0_dp)then
      mu(i)=eta(i)**(1.0_dp/link_power)
   else
      mu(i)=max(ybar,1.0e-8_dp)
   end if
end do
end subroutine tweedie_glm_fit

subroutine tweedie_profile_grid(x, y, powers, link_power, result, weights, offset, &
   method, phi_method, tol, maxit)
real(dp), intent(in) :: x(:,:), y(:), powers(:), link_power
class(tweedie_profile_result), intent(out) :: result
real(dp), intent(in), optional :: weights(:), offset(:), tol
integer, intent(in), optional :: method, phi_method, maxit
real(dp), allocatable :: beta(:), mu(:)
real(dp) :: phi, phi_start
integer :: i, im, iphi, nit, mit
logical :: ok
im=tweedie_method_auto
if(present(method))im=method
iphi=tweedie_phi_mle_method
if(present(phi_method))iphi=phi_method
mit=25
if(present(maxit))mit=maxit
allocate(result%power(size(powers)),result%phi(size(powers)),result%loglik(size(powers)))
allocate(result%beta(size(x,2),size(powers)),result%converged(size(powers)))
result%power=powers
result%phi=ieee_value(0.0_dp,ieee_quiet_nan)
result%loglik=-huge(1.0_dp)
result%beta=0.0_dp
result%converged=.false.
phi_start=ieee_value(0.0_dp,ieee_quiet_nan)
do i=1,size(powers)
   call tweedie_glm_fit(x,y,powers(i),link_power,beta,mu,ok,nit,weights,offset,tol,mit)
   result%converged(i)=ok
   if(.not.ok)cycle
   if(iphi==tweedie_phi_saddle_method)then
      phi=sum(tweedie_dev(y,mu,powers(i)))/real(size(y),dp)
   else
      if(ieee_is_finite(phi_start).and.phi_start>0.0_dp)then
         phi=tweedie_phi_mle(y,mu,powers(i),im,phi_start=phi_start)
      else
         phi=tweedie_phi_mle(y,mu,powers(i),im)
      end if
   end if
   result%phi(i)=phi
   phi_start=phi
   result%loglik(i)=tweedie_loglik(y,mu,phi,powers(i),im)
   result%beta(:,i)=beta
end do
end subroutine tweedie_profile_grid

end module tweedie_likelihood_mod
