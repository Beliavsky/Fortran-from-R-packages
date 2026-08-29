! SPDX-License-Identifier: GPL-2.0-or-later
!
! Modern Fortran translation of the computational core of the R package
! CompQuadForm (GPL >= 2). See NOTICE and upstream metadata for attribution.
module compquadform_mod
use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
use r_compat, only: dp, normal_cdf, pgamma, integrate_result_t, integrate
implicit none
private

real(kind=dp), parameter :: pi = acos(-1.0_dp)

type, public :: davies_result_t
   real(kind=dp) :: trace(7) = 0.0_dp
   integer :: ifault = 0
   real(kind=dp) :: qq = 0.0_dp
end type davies_result_t

type, public :: farebrother_result_t
   real(kind=dp) :: dnsty = 0.0_dp
   integer :: ifault = 0
   real(kind=dp) :: qq = 0.0_dp
end type farebrother_result_t

type, public :: imhof_result_t
   real(kind=dp) :: qq = 0.0_dp
   real(kind=dp) :: abserr = 0.0_dp
end type imhof_result_t

public :: davies, farebrother, imhof, liu

contains

function davies(q, lambda, h, delta, sigma, lim, acc) result(out)
! Survival probability for a quadratic form using Davies's AS 155 algorithm.
real(kind=dp), intent(in) :: q
real(kind=dp), intent(in) :: lambda(:)
integer, intent(in), optional :: h(:)
real(kind=dp), intent(in), optional :: delta(:)
real(kind=dp), intent(in), optional :: sigma
integer, intent(in), optional :: lim
real(kind=dp), intent(in), optional :: acc
type(davies_result_t) :: out
integer, allocatable :: nv(:), th(:)
real(kind=dp), allocatable :: lb(:), nc(:)
integer :: r, limv, count, j, nj, nt, ntm, i
integer :: rats(4)
real(kind=dp) :: sigmav, accv, sigsq, lmax, lmin, meanv, cv
real(kind=dp) :: intl, ersm, acc1, almx, xlim, xnt, xntm
real(kind=dp) :: utx, tausq, sdv, intv, intv1, x, up, un, d1, d2
real(kind=dp) :: lj, ncj, qfval, denom, tmp1, tmp2
logical :: ndtsrt, fail, limit_hit, do_main

r = size(lambda)
out%trace = 0.0_dp
out%ifault = 0
out%qq = 2.0_dp
if (r < 1) then
   out%ifault = 3
   return
end if
if (present(h)) then
   if (size(h) /= r) then
      out%ifault = 3
      return
   end if
end if
if (present(delta)) then
   if (size(delta) /= r) then
      out%ifault = 3
      return
   end if
end if

allocate(lb(r), nc(r), nv(r), th(r))
lb = lambda
if (present(h)) then
   nv = h
else
   nv = 1
end if
if (present(delta)) then
   nc = delta
else
   nc = 0.0_dp
end if
sigmav = 0.0_dp
if (present(sigma)) sigmav = sigma
limv = 10000
if (present(lim)) limv = lim
accv = 1.0e-4_dp
if (present(acc)) accv = acc

cv = q
count = 0
intl = 0.0_dp
ersm = 0.0_dp
qfval = -1.0_dp
acc1 = accv
ndtsrt = .true.
fail = .false.
limit_hit = .false.
xlim = real(limv, kind=dp)
rats = [1, 2, 4, 8]

sigsq = sigmav * sigmav
sdv = sigsq
lmax = 0.0_dp
lmin = 0.0_dp
meanv = 0.0_dp
do j = 1, r
   nj = nv(j)
   lj = lb(j)
   ncj = nc(j)
   if (nj < 0 .or. ncj < 0.0_dp) then
      out%ifault = 3
      call finish()
      return
   end if
   sdv = sdv + lj * lj * (2.0_dp * real(nj, dp) + 4.0_dp * ncj)
   meanv = meanv + lj * (real(nj, dp) + ncj)
   if (lmax < lj) then
      lmax = lj
   else if (lmin > lj) then
      lmin = lj
   end if
end do
if (sdv == 0.0_dp) then
   if (cv > 0.0_dp) then
      qfval = 1.0_dp
   else
      qfval = 0.0_dp
   end if
   call finish()
   return
end if
if (lmin == 0.0_dp .and. lmax == 0.0_dp .and. sigmav == 0.0_dp) then
   out%ifault = 3
   call finish()
   return
end if

sdv = sqrt(sdv)
almx = max(-lmin, lmax)
utx = 16.0_dp / sdv
up = 4.5_dp / sdv
un = -up
call findu(utx, 0.5_dp * acc1)
if (limit_hit) then
   out%ifault = 4
   call finish()
   return
end if

if (cv /= 0.0_dp .and. almx > 0.07_dp * sdv) then
   denom = cfe(cv)
   if (limit_hit) then
      out%ifault = 4
      call finish()
      return
   end if
   if (denom /= 0.0_dp) then
      tausq = 0.25_dp * acc1 / denom
      if (fail) then
         fail = .false.
      else if (truncation(utx, tausq) < 0.2_dp * acc1) then
         if (limit_hit) then
            out%ifault = 4
            call finish()
            return
         end if
         sigsq = sigsq + tausq
         call findu(utx, 0.25_dp * acc1)
         if (limit_hit) then
            out%ifault = 4
            call finish()
            return
         end if
         out%trace(6) = sqrt(tausq)
      end if
   end if
end if
out%trace(5) = utx
acc1 = 0.5_dp * acc1

do
   d1 = ctff(acc1, up) - cv
   if (limit_hit) then
      out%ifault = 4
      call finish()
      return
   end if
   if (d1 < 0.0_dp) then
      qfval = 1.0_dp
      call finish()
      return
   end if
   d2 = cv - ctff(acc1, un)
   if (limit_hit) then
      out%ifault = 4
      call finish()
      return
   end if
   if (d2 < 0.0_dp) then
      qfval = 0.0_dp
      call finish()
      return
   end if
   intv = 2.0_dp * pi / max(d1, d2)
   xnt = utx / intv
   xntm = 3.0_dp / sqrt(acc1)
   do_main = .true.
   if (xnt > 1.5_dp * xntm) then
      do_main = .false.
      if (xntm > xlim) then
         out%ifault = 1
         call finish()
         return
      end if
      ntm = int(floor(xntm + 0.5_dp))
      intv1 = utx / real(ntm, dp)
      x = 2.0_dp * pi / intv1
      if (x <= abs(cv)) then
         do_main = .true.
      else
         tmp1 = cfe(cv - x)
         tmp2 = cfe(cv + x)
         if (limit_hit) then
            out%ifault = 4
            call finish()
            return
         end if
         if (fail) then
            do_main = .true.
         else
            tausq = 0.33_dp * acc1 / (1.1_dp * (tmp1 + tmp2))
            acc1 = 0.67_dp * acc1
            call davies_integrate(ntm, intv1, tausq, .false.)
            xlim = xlim - xntm
            sigsq = sigsq + tausq
            out%trace(3) = out%trace(3) + 1.0_dp
            out%trace(2) = out%trace(2) + real(ntm + 1, dp)
            call findu(utx, 0.25_dp * acc1)
            if (limit_hit) then
               out%ifault = 4
               call finish()
               return
            end if
            acc1 = 0.75_dp * acc1
         end if
      end if
   end if
   if (do_main) exit
end do

out%trace(4) = intv
if (xnt > xlim) then
   out%ifault = 1
   call finish()
   return
end if
nt = int(floor(xnt + 0.5_dp))
call davies_integrate(nt, intv, 0.0_dp, .true.)
out%trace(3) = out%trace(3) + 1.0_dp
out%trace(2) = out%trace(2) + real(nt + 1, dp)
qfval = 0.5_dp - intl
out%trace(1) = ersm
up = ersm
x = up + accv / 10.0_dp
do i = 1, 4
   if (real(rats(i), dp) * x == real(rats(i), dp) * up) out%ifault = 2
end do
call finish()

contains

subroutine finish()
out%trace(7) = real(count, dp)
out%qq = 1.0_dp - qfval
end subroutine finish

subroutine counter_step()
count = count + 1
if (count > limv) limit_hit = .true.
end subroutine counter_step

pure real(kind=dp) function square(v) result(y)
real(kind=dp), intent(in) :: v
y = v * v
end function square

pure real(kind=dp) function cube(v) result(y)
real(kind=dp), intent(in) :: v
y = v * v * v
end function cube

real(kind=dp) function exp1(v) result(y)
real(kind=dp), intent(in) :: v
if (v < -50.0_dp) then
   y = 0.0_dp
else
   y = exp(v)
end if
end function exp1

real(kind=dp) function log1(v, first) result(yout)
real(kind=dp), intent(in) :: v
logical, intent(in) :: first
real(kind=dp) :: s, s1, term, y, kreal
if (abs(v) > 0.1_dp) then
   if (first) then
      yout = log(1.0_dp + v)
   else
      yout = log(1.0_dp + v) - v
   end if
   return
end if
y = v / (2.0_dp + v)
term = 2.0_dp * cube(y)
kreal = 3.0_dp
if (first) then
   s = 2.0_dp * y
else
   s = -v * y
end if
y = square(y)
do
   s1 = s + term / kreal
   if (s1 == s) exit
   kreal = kreal + 2.0_dp
   term = term * y
   s = s1
end do
yout = s
end function log1

subroutine order_lambdas()
integer :: jj, kk
real(kind=dp) :: ljj
do jj = 1, r
   ljj = abs(lb(jj))
   kk = jj - 1
   do while (kk >= 1)
      if (ljj > abs(lb(th(kk)))) then
         th(kk + 1) = th(kk)
         kk = kk - 1
      else
         exit
      end if
   end do
   th(kk + 1) = jj
end do
ndtsrt = .false.
end subroutine order_lambdas

real(kind=dp) function errbd(u, cx) result(yout)
real(kind=dp), intent(in) :: u
real(kind=dp), intent(out) :: cx
real(kind=dp) :: uu, sum1, xloc, yloc, ljj, ncjj
integer :: jj, njj
call counter_step()
if (limit_hit) then
   cx = meanv
   yout = 1.0_dp
   return
end if
cx = u * sigsq
sum1 = u * cx
uu = 2.0_dp * u
do jj = r, 1, -1
   njj = nv(jj)
   ljj = lb(jj)
   ncjj = nc(jj)
   xloc = uu * ljj
   yloc = 1.0_dp - xloc
   cx = cx + ljj * (ncjj / yloc + real(njj, dp)) / yloc
   sum1 = sum1 + ncjj * square(xloc / yloc) + &
      real(njj, dp) * (square(xloc) / yloc + log1(-xloc, .false.))
end do
yout = exp1(-0.5_dp * sum1)
end function errbd

real(kind=dp) function ctff(accx, upn) result(c2out)
real(kind=dp), intent(in) :: accx
real(kind=dp), intent(inout) :: upn
real(kind=dp) :: u1, u2, u, rb, xconst, c1loc, c2loc, eb
u2 = upn
u1 = 0.0_dp
c1loc = meanv
rb = 2.0_dp * merge(lmax, lmin, u2 > 0.0_dp)
do
   u = u2 / (1.0_dp + u2 * rb)
   eb = errbd(u, c2loc)
   if (limit_hit) then
      c2out = c2loc
      return
   end if
   if (eb <= accx) exit
   u1 = u2
   c1loc = c2loc
   u2 = 2.0_dp * u2
end do
do
   if (c2loc == meanv) exit
   u = (c1loc - meanv) / (c2loc - meanv)
   if (u >= 0.9_dp) exit
   u = 0.5_dp * (u1 + u2)
   eb = errbd(u / (1.0_dp + u * rb), xconst)
   if (limit_hit) then
      c2out = c2loc
      return
   end if
   if (eb > accx) then
      u1 = u
      c1loc = xconst
   else
      u2 = u
      c2loc = xconst
   end if
end do
upn = u2
c2out = c2loc
end function ctff

real(kind=dp) function truncation(uin, tausq_in) result(yout)
real(kind=dp), intent(in) :: uin, tausq_in
real(kind=dp) :: u, sum1, sum2, prod1, prod2, prod3
real(kind=dp) :: ljj, ncjj, xloc, yloc, err1, err2
integer :: jj, njj, ss
call counter_step()
if (limit_hit) then
   yout = 1.0_dp
   return
end if
sum1 = 0.0_dp
prod2 = 0.0_dp
prod3 = 0.0_dp
ss = 0
sum2 = (sigsq + tausq_in) * square(uin)
prod1 = 2.0_dp * sum2
u = 2.0_dp * uin
do jj = 1, r
   ljj = lb(jj)
   ncjj = nc(jj)
   njj = nv(jj)
   xloc = square(u * ljj)
   sum1 = sum1 + ncjj * xloc / (1.0_dp + xloc)
   if (xloc > 1.0_dp) then
      prod2 = prod2 + real(njj, dp) * log(xloc)
      prod3 = prod3 + real(njj, dp) * log1(xloc, .true.)
      ss = ss + njj
   else
      prod1 = prod1 + real(njj, dp) * log1(xloc, .true.)
   end if
end do
sum1 = 0.5_dp * sum1
prod2 = prod1 + prod2
prod3 = prod1 + prod3
xloc = exp1(-sum1 - 0.25_dp * prod2) / pi
yloc = exp1(-sum1 - 0.25_dp * prod3) / pi
if (ss == 0) then
   err1 = 1.0_dp
else
   err1 = 2.0_dp * xloc / real(ss, dp)
end if
if (prod3 > 1.0_dp) then
   err2 = 2.5_dp * yloc
else
   err2 = 1.0_dp
end if
err1 = min(err1, err2)
xloc = 0.5_dp * sum2
if (xloc <= yloc) then
   err2 = 1.0_dp
else
   err2 = yloc / xloc
end if
yout = min(err1, err2)
end function truncation

subroutine findu(utx_inout, accx)
real(kind=dp), intent(inout) :: utx_inout
real(kind=dp), intent(in) :: accx
real(kind=dp) :: u, ut, tv
real(kind=dp), parameter :: divis(4) = [2.0_dp, 1.4_dp, 1.2_dp, 1.1_dp]
integer :: ii
ut = utx_inout
u = ut / 4.0_dp
tv = truncation(u, 0.0_dp)
if (limit_hit) return
if (tv > accx) then
   u = ut
   do
      tv = truncation(u, 0.0_dp)
      if (limit_hit) return
      if (tv <= accx) exit
      ut = ut * 4.0_dp
      u = ut
   end do
else
   ut = u
   u = u / 4.0_dp
   do
      tv = truncation(u, 0.0_dp)
      if (limit_hit) return
      if (tv > accx) exit
      ut = u
      u = u / 4.0_dp
   end do
end if
do ii = 1, 4
   u = ut / divis(ii)
   tv = truncation(u, 0.0_dp)
   if (limit_hit) return
   if (tv <= accx) ut = u
end do
utx_inout = ut
end subroutine findu

subroutine davies_integrate(nterm, interv, tausq_in, mainx)
integer, intent(in) :: nterm
real(kind=dp), intent(in) :: interv, tausq_in
logical, intent(in) :: mainx
real(kind=dp) :: inpi, u, sum1, sum2, sum3, xloc, yloc, zloc
integer :: k, jj, njj
inpi = interv / pi
do k = nterm, 0, -1
   u = (real(k, dp) + 0.5_dp) * interv
   sum1 = -2.0_dp * u * cv
   sum2 = abs(sum1)
   sum3 = -0.5_dp * sigsq * square(u)
   do jj = r, 1, -1
      njj = nv(jj)
      xloc = 2.0_dp * lb(jj) * u
      yloc = square(xloc)
      sum3 = sum3 - 0.25_dp * real(njj, dp) * log1(yloc, .true.)
      yloc = nc(jj) * xloc / (1.0_dp + yloc)
      zloc = real(njj, dp) * atan(xloc) + yloc
      sum1 = sum1 + zloc
      sum2 = sum2 + abs(zloc)
      sum3 = sum3 - 0.5_dp * xloc * yloc
   end do
   xloc = inpi * exp1(sum3) / u
   if (.not. mainx) xloc = xloc * (1.0_dp - exp1(-0.5_dp * tausq_in * square(u)))
   sum1 = sin(0.5_dp * sum1) * xloc
   sum2 = 0.5_dp * sum2 * xloc
   intl = intl + sum1
   ersm = ersm + sum2
end do
end subroutine davies_integrate

real(kind=dp) function cfe(xin) result(yout)
real(kind=dp), intent(in) :: xin
real(kind=dp), parameter :: log28 = 0.0866_dp
real(kind=dp) :: axl, axl1, axl2, sxl, sum1, ljj
integer :: jj, kk, t
call counter_step()
if (limit_hit) then
   yout = 1.0_dp
   return
end if
if (ndtsrt) call order_lambdas()
axl = abs(xin)
if (xin > 0.0_dp) then
   sxl = 1.0_dp
else
   sxl = -1.0_dp
end if
sum1 = 0.0_dp
do jj = r, 1, -1
   t = th(jj)
   if (lb(t) * sxl > 0.0_dp) then
      ljj = abs(lb(t))
      axl1 = axl - ljj * (real(nv(t), dp) + nc(t))
      axl2 = ljj / log28
      if (axl1 > axl2) then
         axl = axl1
      else
         if (axl > axl2) axl = axl2
         sum1 = (axl - axl1) / ljj
         do kk = jj - 1, 1, -1
            sum1 = sum1 + real(nv(th(kk)), dp) + nc(th(kk))
         end do
         exit
      end if
   end if
end do
if (sum1 > 100.0_dp) then
   fail = .true.
   yout = 1.0_dp
else
   yout = 2.0_dp**(sum1 / 4.0_dp) / (pi * square(axl))
end if
end function cfe

end function davies

function farebrother(q, lambda, h, delta, maxit, eps, mode) result(out)
! Survival probability using Farebrother's AS 204 / Ruben expansion.
real(kind=dp), intent(in) :: q
real(kind=dp), intent(in) :: lambda(:)
integer, intent(in), optional :: h(:)
real(kind=dp), intent(in), optional :: delta(:)
integer, intent(in), optional :: maxit
real(kind=dp), intent(in), optional :: eps, mode
type(farebrother_result_t) :: out
integer, allocatable :: mult(:)
real(kind=dp), allocatable :: nc(:), gamma(:), theta(:), avec(:), bvec(:)
integer :: n, maxitv, i, k, m, j
real(kind=dp) :: epsv, modev, ao, aoinv, z, bbeta, eps2, hold, hold2
real(kind=dp) :: sumv, sum1, dans, lans, pans, prbty, tol, resv
logical :: converged

n = size(lambda)
out%dnsty = 0.0_dp
out%ifault = 0
out%qq = 0.0_dp
if (n < 1) then
   out%ifault = 2
   out%qq = 3.0_dp
   return
end if
if (present(h)) then
   if (size(h) /= n) then
      out%ifault = 2
      out%qq = 3.0_dp
      return
   end if
end if
if (present(delta)) then
   if (size(delta) /= n) then
      out%ifault = 2
      out%qq = 3.0_dp
      return
   end if
end if
maxitv = 100000
if (present(maxit)) maxitv = maxit
epsv = 1.0e-10_dp
if (present(eps)) epsv = eps
modev = 1.0_dp
if (present(mode)) modev = mode
allocate(mult(n), nc(n))
if (present(h)) then
   mult = h
else
   mult = 1
end if
if (present(delta)) then
   nc = delta
else
   nc = 0.0_dp
end if

if (q <= 0.0_dp .or. maxitv < 1 .or. epsv <= 0.0_dp) then
   resv = -2.0_dp
   out%ifault = 2
   out%qq = 1.0_dp - resv
   return
end if

tol = -200.0_dp
sumv = lambda(1)
bbeta = sumv
do i = 1, n
   hold = lambda(i)
   if (hold <= 0.0_dp .or. mult(i) < 1 .or. nc(i) < 0.0_dp) then
      resv = -7.0_dp
      out%ifault = -i
      out%qq = 1.0_dp - resv
      return
   end if
   if (bbeta > hold) bbeta = hold
   if (sumv < hold) sumv = hold
end do
if (modev > 0.0_dp) then
   bbeta = modev * bbeta
else
   bbeta = 2.0_dp / (1.0_dp / bbeta + 1.0_dp / sumv)
end if

allocate(gamma(n), theta(n), avec(maxitv), bvec(maxitv))
k = 0
sumv = 1.0_dp
sum1 = 0.0_dp
do i = 1, n
   hold = bbeta / lambda(i)
   gamma(i) = 1.0_dp - hold
   sumv = sumv * hold**mult(i)
   sum1 = sum1 + nc(i)
   k = k + mult(i)
   theta(i) = 1.0_dp
end do
ao = exp(0.5_dp * (log(sumv) - sum1))
if (ao <= 0.0_dp) then
   resv = 0.0_dp
   out%dnsty = 0.0_dp
   out%ifault = 1
   out%qq = 1.0_dp
   return
end if

z = q / bbeta
if (mod(k, 2) == 0) then
   i = 2
   lans = -0.5_dp * z
   dans = exp(lans)
   pans = 1.0_dp - dans
else
   i = 1
   lans = -0.5_dp * (z + log(z)) - 0.22579135264473_dp
   dans = exp(lans)
   pans = normal_cdf(sqrt(z)) - normal_cdf(-sqrt(z))
end if
k = k - 2
do j = i, k, 2
   if (lans < tol) then
      lans = lans + log(z / real(j, dp))
      dans = exp(lans)
   else
      dans = dans * z / real(j, dp)
   end if
   pans = pans - dans
end do

prbty = pans
out%dnsty = dans
eps2 = epsv / ao
aoinv = 1.0_dp / ao
sumv = aoinv - 1.0_dp
converged = .false.
do m = 1, maxitv
   sum1 = 0.0_dp
   do i = 1, n
      hold = theta(i)
      hold2 = hold * gamma(i)
      theta(i) = hold2
      sum1 = sum1 + hold2 * real(mult(i), dp) + real(m, dp) * nc(i) * (hold - hold2)
   end do
   sum1 = 0.5_dp * sum1
   bvec(m) = sum1
   do i = m - 1, 1, -1
      sum1 = sum1 + bvec(i) * avec(m - i)
   end do
   sum1 = sum1 / real(m, dp)
   avec(m) = sum1
   k = k + 2
   if (lans < tol) then
      lans = lans + log(z / real(k, dp))
      dans = exp(lans)
   else
      dans = dans * z / real(k, dp)
   end if
   pans = pans - dans
   sumv = sumv - sum1
   out%dnsty = out%dnsty + dans * sum1
   sum1 = pans * sum1
   prbty = prbty + sum1
   if (prbty < -aoinv) then
      resv = -3.0_dp
      out%ifault = 3
      out%qq = 1.0_dp - resv
      return
   end if
   if (abs(pans * sumv) < eps2 .and. abs(sum1) < eps2) then
      converged = .true.
      exit
   end if
end do
if (.not. converged .or. m >= maxitv) out%ifault = 4
out%dnsty = ao * out%dnsty / (2.0_dp * bbeta)
prbty = ao * prbty
if (prbty < 0.0_dp .or. prbty > 1.0_dp) then
   out%ifault = out%ifault + 5
else if (out%dnsty < 0.0_dp) then
   out%ifault = out%ifault + 6
end if
resv = prbty
out%qq = 1.0_dp - resv
end function farebrother

function imhof(q, lambda, h, delta, epsabs, epsrel, limit) result(out)
! Survival probability using Imhof's inversion integral.
real(kind=dp), intent(in) :: q
real(kind=dp), intent(in) :: lambda(:)
integer, intent(in), optional :: h(:)
real(kind=dp), intent(in), optional :: delta(:)
real(kind=dp), intent(in), optional :: epsabs, epsrel
integer, intent(in), optional :: limit
type(imhof_result_t) :: out
integer, allocatable :: hv(:)
real(kind=dp), allocatable :: dv(:)
real(kind=dp) :: epsabsv, epsrelv, tol, nanv
integer :: r, limitv
integer :: j
type(integrate_result_t) :: integ

r = size(lambda)
nanv = ieee_value(0.0_dp, ieee_quiet_nan)
out%qq = nanv
out%abserr = nanv
if (r < 1) return
if (present(h)) then
   if (size(h) /= r) return
end if
if (present(delta)) then
   if (size(delta) /= r) return
end if
allocate(hv(r), dv(r))
if (present(h)) then
   hv = h
else
   hv = 1
end if
if (present(delta)) then
   dv = delta
else
   dv = 0.0_dp
end if
if (any(dv < 0.0_dp)) return

epsabsv = 1.0e-6_dp
if (present(epsabs)) epsabsv = epsabs
epsrelv = 1.0e-6_dp
if (present(epsrel)) epsrelv = epsrel
limitv = 10000
if (present(limit)) limitv = limit
if (limitv >= 99999999 .or. limitv < 2) return
if (epsabsv <= 0.0_dp .or. epsrelv <= 0.0_dp) return

tol = min(epsabsv, epsrelv)
integ = integrate(imhof_integrand, 0.0_dp, huge(1.0_dp), &
   rel_tol=tol, subdivisions=limitv)
out%qq = 0.5_dp + integ%value / pi
out%abserr = integ%abs_error

contains

function imhof_integrand(u) result(v)
real(kind=dp), intent(in) :: u
real(kind=dp) :: v, theta_v, rho_v, x, sumlim
if (abs(u) <= sqrt(epsilon(1.0_dp))) then
   sumlim = 0.0_dp
   do j = 1, r
      sumlim = sumlim + lambda(j) * (real(hv(j), dp) + dv(j))
   end do
   v = 0.5_dp * (sumlim - q)
   return
end if
theta_v = 0.0_dp
rho_v = 1.0_dp
do j = 1, r
   x = lambda(j) * u
   theta_v = theta_v + real(hv(j), dp) * atan(x) + dv(j) * x / (1.0_dp + x * x)
   rho_v = rho_v * (1.0_dp + x * x)**(0.25_dp * real(hv(j), dp)) * &
      exp(0.5_dp * dv(j) * x * x / (1.0_dp + x * x))
end do
theta_v = 0.5_dp * theta_v - 0.5_dp * q * u
v = sin(theta_v) / (u * rho_v)
if (.not. ieee_is_finite(v)) v = 0.0_dp
end function imhof_integrand

end function imhof

function liu(q, lambda, h, delta) result(qq)
! Liu-Tang-Zhang moment-matching approximation to the survival probability.
real(kind=dp), intent(in) :: q
real(kind=dp), intent(in) :: lambda(:)
integer, intent(in), optional :: h(:)
real(kind=dp), intent(in), optional :: delta(:)
real(kind=dp) :: qq
integer, allocatable :: hv(:)
real(kind=dp), allocatable :: dv(:)
integer :: r
real(kind=dp) :: c1, c2, c3, c4, s1, s2, muq, sigmaq, tstar
real(kind=dp) :: a, ncp, ell, mux, sigmax, xstar, nanv

r = size(lambda)
nanv = ieee_value(0.0_dp, ieee_quiet_nan)
qq = nanv
if (r < 1) return
if (present(h)) then
   if (size(h) /= r) return
end if
if (present(delta)) then
   if (size(delta) /= r) return
end if
allocate(hv(r), dv(r))
if (present(h)) then
   hv = h
else
   hv = 1
end if
if (present(delta)) then
   dv = delta
else
   dv = 0.0_dp
end if
if (any(dv < 0.0_dp)) return

c1 = sum(lambda * real(hv, dp)) + sum(lambda * dv)
c2 = sum(lambda**2 * real(hv, dp)) + 2.0_dp * sum(lambda**2 * dv)
c3 = sum(lambda**3 * real(hv, dp)) + 3.0_dp * sum(lambda**3 * dv)
c4 = sum(lambda**4 * real(hv, dp)) + 4.0_dp * sum(lambda**4 * dv)
if (c2 <= 0.0_dp) return
s1 = c3 / c2**1.5_dp
s2 = c4 / c2**2
muq = c1
sigmaq = sqrt(2.0_dp * c2)
tstar = (q - muq) / sigmaq
if (s1 * s1 > s2) then
   a = 1.0_dp / (s1 - sqrt(s1 * s1 - s2))
   ncp = s1 * a**3 - a**2
   ell = a**2 - 2.0_dp * ncp
else
   a = 1.0_dp / s1
   ncp = 0.0_dp
   ell = c2**3 / c3**2
end if
mux = ell + ncp
sigmax = sqrt(2.0_dp) * a
xstar = tstar * sigmax + mux
qq = 1.0_dp - precise_pchisq(xstar, ell, ncp)
end function liu

pure function precise_pchisq(x, df, ncp) result(p)
! Exact-enough noncentral chi-square CDF built from r_compat's exact pgamma.
! The supplied r_compat pchisq uses a Wilson-Hilferty approximation for central
! chi-square terms; that approximation materially degrades Liu parity.  This
! package-local helper adds the missing precise Poisson-mixture evaluation.
real(kind=dp), intent(in) :: x, df, ncp
real(kind=dp) :: p
real(kind=dp) :: lam, w, sumw, term
integer :: j, jm
if (df <= 0.0_dp .or. ncp < 0.0_dp) then
   p = ieee_value(0.0_dp, ieee_quiet_nan)
   return
end if
if (x <= 0.0_dp) then
   p = 0.0_dp
   return
end if
if (.not. ieee_is_finite(x)) then
   p = 1.0_dp
   return
end if
if (ncp == 0.0_dp) then
   p = pgamma(0.5_dp * x, 0.5_dp * df)
   return
end if
lam = 0.5_dp * ncp
jm = int(floor(lam))
w = exp(-lam + real(jm, dp) * log(lam) - log_gamma(real(jm + 1, dp)))
p = w * pgamma(0.5_dp * x, 0.5_dp * df + real(jm, dp))
sumw = w
term = w
j = jm
do while (j > 0)
   term = term * real(j, dp) / lam
   j = j - 1
   p = p + term * pgamma(0.5_dp * x, 0.5_dp * df + real(j, dp))
   sumw = sumw + term
end do
term = w
j = jm
do j = jm + 1, jm + max(1000, int(50.0_dp * sqrt(lam + 1.0_dp)))
   term = term * lam / real(j, dp)
   p = p + term * pgamma(0.5_dp * x, 0.5_dp * df + real(j, dp))
   sumw = sumw + term
   if (term <= 1.0e-15_dp .and. sumw >= 1.0_dp - 1.0e-14_dp) exit
end do
p = max(0.0_dp, min(1.0_dp, p))
end function precise_pchisq

end module compquadform_mod
