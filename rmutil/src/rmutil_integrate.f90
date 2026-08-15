! rmutil computational translation
! Copyright (C) 1998-2001 J.K. Lindsey
! Copyright (C) 2026 OpenAI (modern Fortran translation)
! SPDX-License-Identifier: GPL-2.0-or-later
module rmutil_integrate
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use rmutil_kinds, only : dp, pi
   implicit none
   private
   public :: integrate_romberg, integrate_adaptive, integrate_2d
   public :: toms614_integrate

   abstract interface
      function scalar_function(x) result(y)
         import dp
         real(dp), intent(in) :: x
         real(dp) :: y
      end function scalar_function
      function function_2d(x, y) result(z)
         import dp
         real(dp), intent(in) :: x, y
         real(dp) :: z
      end function function_2d
   end interface

contains

   recursive real(dp) function integrate_adaptive(f, a, b, tol, max_depth) result(ans)
      procedure(scalar_function) :: f
      real(dp), intent(in) :: a, b
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: max_depth
      real(dp) :: eps, fa, fb, fc, c, s
      integer :: depth
      eps = 1.0e-10_dp
      if (present(tol)) eps = tol
      depth = 24
      if (present(max_depth)) depth = max_depth
      c = 0.5_dp*(a+b)
      fa = f(a)
      fb = f(b)
      fc = f(c)
      s = (b-a)*(fa+4.0_dp*fc+fb)/6.0_dp
      ans = recurse(a, b, fa, fb, fc, s, eps, depth)
   contains
      recursive function recurse(a0,b0,fa0,fb0,fc0,s0,eps0,lev) result(v)
         real(dp), intent(in) :: a0,b0,fa0,fb0,fc0,s0,eps0
         integer, intent(in) :: lev
         real(dp) :: v, c0, d0, e0, fd, fe, sl, sr, s2
         c0 = 0.5_dp*(a0+b0)
         d0 = 0.5_dp*(a0+c0)
         e0 = 0.5_dp*(c0+b0)
         fd = f(d0)
         fe = f(e0)
         sl = (c0-a0)*(fa0+4.0_dp*fd+fc0)/6.0_dp
         sr = (b0-c0)*(fc0+4.0_dp*fe+fb0)/6.0_dp
         s2 = sl + sr
         if (lev <= 0 .or. abs(s2-s0) <= 15.0_dp*eps0) then
            v = s2 + (s2-s0)/15.0_dp
         else
            v = recurse(a0,c0,fa0,fc0,fd,sl,eps0/2.0_dp,lev-1) + &
               recurse(c0,b0,fc0,fb0,fe,sr,eps0/2.0_dp,lev-1)
         end if
      end function recurse
   end function integrate_adaptive

   real(dp) function romberg_finite(f, a, b, tol, max_order) result(ans)
      procedure(scalar_function) :: f
      real(dp), intent(in) :: a, b, tol
      integer, intent(in) :: max_order
      real(dp), allocatable :: r(:,:)
      real(dp) :: h, sumv
      integer :: i, j, k, nnew
      allocate(r(max_order,max_order))
      r = 0.0_dp
      h = b-a
      r(1,1) = 0.5_dp*h*(f(a)+f(b))
      do i = 2, max_order
         nnew = 2**(i-2)
         sumv = 0.0_dp
         do k = 1, nnew
            sumv = sumv + f(a + (real(2*k-1,dp)*h)/(2.0_dp**real(i-1,dp)))
         end do
         r(i,1) = 0.5_dp*r(i-1,1) + h*sumv/(2.0_dp**real(i-1,dp))
         do j = 2, i
            r(i,j) = r(i,j-1) + (r(i,j-1)-r(i-1,j-1))/(4.0_dp**real(j-1,dp)-1.0_dp)
         end do
         if (abs(r(i,i)-r(i-1,i-1)) <= tol*max(1.0_dp,abs(r(i,i)))) then
            ans = r(i,i)
            return
         end if
      end do
      ans = r(max_order,max_order)
   end function romberg_finite

   real(dp) function integrate_romberg(f, a, b, tol, max_order) result(ans)
      procedure(scalar_function) :: f
      real(dp), intent(in) :: a, b
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: max_order
      real(dp) :: eps
      integer :: nmax
      eps = 1.0e-8_dp
      if (present(tol)) eps = tol
      nmax = 16
      if (present(max_order)) nmax = max_order
      if (ieee_is_finite(a) .and. ieee_is_finite(b)) then
         ans = romberg_finite(f,a,b,eps,nmax)
      else if (.not.ieee_is_finite(a) .and. .not.ieee_is_finite(b)) then
         ans = romberg_finite(left_transform,-1.0_dp,-sqrt(epsilon(1.0_dp)),eps,nmax) + &
            romberg_finite(f,-1.0_dp,1.0_dp,eps,nmax) + &
            romberg_finite(right_transform,sqrt(epsilon(1.0_dp)),1.0_dp,eps,nmax)
      else if (.not.ieee_is_finite(b)) then
         if (a > 0.0_dp) then
            ans = romberg_finite(right_transform,sqrt(epsilon(1.0_dp)),1.0_dp/a,eps,nmax)
         else
            ans = romberg_finite(f,a,1.0_dp,eps,nmax) + &
               romberg_finite(right_transform,sqrt(epsilon(1.0_dp)),1.0_dp,eps,nmax)
         end if
      else
         if (b < 0.0_dp) then
            ans = romberg_finite(left_transform,1.0_dp/b,-sqrt(epsilon(1.0_dp)),eps,nmax)
         else
            ans = romberg_finite(f,-1.0_dp,b,eps,nmax) + &
               romberg_finite(left_transform,-1.0_dp,-sqrt(epsilon(1.0_dp)),eps,nmax)
         end if
      end if
   contains
      function right_transform(x) result(y)
         real(dp), intent(in) :: x
         real(dp) :: y
         y = f(1.0_dp/x)/(x*x)
      end function right_transform
      function left_transform(x) result(y)
         real(dp), intent(in) :: x
         real(dp) :: y
         y = f(1.0_dp/x)/(x*x)
      end function left_transform
   end function integrate_romberg

   real(dp) function integrate_2d(f, ax, bx, ay, by, tol) result(ans)
      procedure(function_2d) :: f
      real(dp), intent(in) :: ax, bx, ay, by
      real(dp), intent(in), optional :: tol
      real(dp), allocatable :: nodes(:), weights(:)
      real(dp) :: tx, ty, x, y, jx, jy, eps
      integer :: i, j, nquad
      eps = 1.0e-6_dp
      if (present(tol)) eps = max(tol,epsilon(1.0_dp))
      if (eps >= 1.0e-5_dp) then
         nquad = 24
      else if (eps >= 1.0e-9_dp) then
         nquad = 48
      else
         nquad = 72
      end if
      allocate(nodes(nquad),weights(nquad))
      call gauss_legendre_rule(nquad,nodes,weights)
      ans = 0.0_dp
      do i = 1, nquad
         tx = nodes(i)
         call map_interval(tx,ax,bx,x,jx)
         do j = 1, nquad
            ty = nodes(j)
            call map_interval(ty,ay,by,y,jy)
            ans = ans + weights(i)*weights(j)*f(x,y)*jx*jy
         end do
      end do
   end function integrate_2d

   subroutine gauss_legendre_rule(n,x,w)
      integer, intent(in) :: n
      real(dp), intent(out) :: x(n), w(n)
      integer :: i, j, m
      real(dp) :: z, z1, p1, p2, p3, pp
      m = (n+1)/2
      do i = 1, m
         z = cos(pi*(real(i,dp)-0.25_dp)/(real(n,dp)+0.5_dp))
         do
            p1 = 1.0_dp
            p2 = 0.0_dp
            do j = 1, n
               p3 = p2
               p2 = p1
               p1 = ((2.0_dp*real(j,dp)-1.0_dp)*z*p2 - real(j-1,dp)*p3)/real(j,dp)
            end do
            pp = real(n,dp)*(z*p1-p2)/(z*z-1.0_dp)
            z1 = z
            z = z1-p1/pp
            if (abs(z-z1) <= 4.0_dp*epsilon(1.0_dp)) exit
         end do
         x(i) = -z
         x(n+1-i) = z
         w(i) = 2.0_dp/((1.0_dp-z*z)*pp*pp)
         w(n+1-i) = w(i)
      end do
   end subroutine gauss_legendre_rule

   subroutine map_interval(t,a,b,x,jac)
      real(dp), intent(in) :: t,a,b
      real(dp), intent(out) :: x,jac
      real(dp) :: u
      if (ieee_is_finite(a) .and. ieee_is_finite(b)) then
         x = 0.5_dp*((b-a)*t+a+b)
         jac = 0.5_dp*(b-a)
      else if (ieee_is_finite(a)) then
         u = (1.0_dp+t)/(1.0_dp-t)
         x = a+u
         jac = 2.0_dp/(1.0_dp-t)**2
      else if (ieee_is_finite(b)) then
         u = (1.0_dp+t)/(1.0_dp-t)
         x = b-u
         jac = 2.0_dp/(1.0_dp-t)**2
      else
         x = tan(0.5_dp*pi*t)
         jac = 0.5_dp*pi/(cos(0.5_dp*pi*t)**2)
      end if
   end subroutine map_interval

   subroutine toms614_integrate(f, a, b, dpar, m, pnorm, eps, inf, quadr)
      ! Modern module/interface port of ACM Algorithm 614 (INTHP), as
      ! distributed in rmutil. The algorithm is due to K. Sikorski,
      ! F. Stenger, and J. Schwing, ACM TOMS 10 (1984), 152-160.
      procedure(scalar_function) :: f
      real(dp), intent(in) :: a, b, dpar, pnorm
      integer, intent(inout) :: m, inf
      real(dp), intent(inout) :: eps
      real(dp), intent(out) :: quadr
      integer :: i, i1, k, l, l1, m1, m2, n, n1
      real(dp) :: alfa, exph, exph0, c, h, s, t, u, v, w
      real(dp) :: c0, e1, h0, h1, s1, v0, v1, v2
      real(dp) :: w1, w2, w3, w4, ba, sr, sq2, cor, sumv
      real(dp) :: eps3, sum1, sum2
      logical :: inf1, inf2

      quadr = 0.0_dp
      if (m < 3) goto 270
      if (pnorm < 1.0_dp .and. pnorm /= 0.0_dp) goto 280
      if (pnorm >= 1.0_dp .and. (dpar <= 0.0_dp .or. dpar > pi/2.0_dp)) goto 280
      if (inf == 4 .and. a >= b) goto 290
      sq2 = sqrt(2.0_dp)
      i1 = inf - 2
      ba = b - a
      n1 = 0
      u = 1.0_dp
10    u = u/10.0_dp
      t = 1.0_dp + u
      if (1.0_dp /= t) goto 10
      u = u*10.0_dp
      if (eps < u) eps = u
      if (pnorm == 0.0_dp) goto 40
      if (pnorm == 1.0_dp) alfa = 1.0_dp
      if (pnorm > 1.0_dp) alfa = (pnorm-1.0_dp)/pnorm
      c = 2.0_dp*pi/(1.0_dp-1.0_dp/exp(pi*sqrt(alfa))) + 4.0_dp**alfa/alfa
      w = log(c/eps)
      w1 = w*w/(pi*pi*alfa)
      n = int(w1)
      if (w1 > real(n,dp)) n = n + 1
      if (w1 == 0.0_dp) n = 1
      n1 = 2*n + 1
      sr = sqrt(alfa*real(n,dp))
      if (n1 <= m) goto 20
      n1 = 1
      n = (m-1)/2
      sr = sqrt(alfa*real(n,dp))
      m = 2*n + 1
      eps = c/exp(pi*sr)
      goto 30
20    m = n1
      n1 = 0
30    h = 2.0_dp*dpar/sr
      sum2 = 0.0_dp
      l1 = n
      k = n
      inf1 = .false.
      inf2 = .false.
      h0 = h
      goto 50
40    h = 1.0_dp
      h0 = 1.0_dp
      eps3 = eps/3.0_dp
      sr = sqrt(eps)
      v1 = eps*10.0_dp
      v2 = v1
      m1 = m - 1
      n = m1/2
      m2 = n
      l1 = 0
      inf1 = .true.
      inf2 = .false.
50    i = 0
      if (inf == 1) sumv = f(0.0_dp)
      if (inf == 2) sumv = f(a+1.0_dp)
      if (inf == 3) sumv = f(a+log(1.0_dp+sq2))/sq2
      if (inf == 4) sumv = f((a+b)/2.0_dp)*ba/4.0_dp
60    exph = exp(h)
      exph0 = exp(h0)
      h1 = h0
      e1 = exph0
      u = 0.0_dp
      cor = 0.0_dp
70    if (i1 < 0) goto 80
      if (i1 == 0) goto 90
      goto 100
80    v = f(h1)
      h1 = h1 + h
      goto 150
90    v = e1*f(a+e1)
      e1 = e1*exph
      goto 150
100   if (inf == 4) goto 140
      w1 = sqrt(e1+1.0_dp/e1)
      w2 = sqrt(e1)
      if (e1 < 0.1_dp) goto 110
      s = log(e1+w1*w2)
      goto 130
110   w3 = e1
      w4 = e1*e1
      c0 = 1.0_dp
      s = e1
      s1 = e1
      t = 0.0_dp
120   c0 = -c0*(0.5_dp+t)*(2.0_dp*t+1.0_dp)/(2.0_dp*t+3.0_dp)/(t+1.0_dp)
      t = t + 1.0_dp
      w3 = w3*w4
      s = s + c0*w3
      if (s == s1) goto 130
      s1 = s
      goto 120
130   v = w2/w1*f(a+s)
      e1 = e1*exph
      goto 150
140   w1 = e1 + 1.0_dp
      v = e1/(w1*w1)*f((a+b*e1)/w1)*ba
      e1 = e1*exph
150   i = i + 1
      sum1 = u + v
      if (abs(u) < abs(v)) goto 160
      cor = v - (sum1-u) + cor
      goto 170
160   cor = u - (sum1-v) + cor
170   u = sum1
      if (i < l1) goto 70
      if (inf1) goto 190
      if (inf2) goto 210
      l1 = k
180   inf2 = .true.
      i = 0
      exph = 1.0_dp/exph
      h0 = -h0
      e1 = 1.0_dp/exph0
      h1 = h0
      h = -h
      goto 70
190   v0 = v1
      v1 = v2
      v2 = abs(v)
      if (v0+v1+v2 <= eps3) goto 200
      if (i < m2) goto 70
      n1 = 5
200   if (inf2) k = i
      if (.not.inf2) l = i
      v1 = 10.0_dp*eps
      v2 = v1
      m2 = m1 - l
      if (.not.inf2) goto 180
      if (n1 == 5) goto 260
      sum2 = sum1 + cor + sumv
      m2 = 2*(k+l)
      if (m2 > m1) goto 240
      inf1 = .false.
      inf2 = .false.
      l1 = l
      i = 0
      h = -h
      h0 = h/2.0_dp
      goto 60
210   if (pnorm >= 1.0_dp) goto 220
      h = -h
      sum1 = (sum1+cor)*h
      w1 = (sum1+sum2)/2.0_dp
      if (abs(sum1-sum2) <= sr) goto 230
      m2 = 2*m2
      if (m2 > m1) goto 250
      i = 0
      k = 2*k
      l = 2*l
      l1 = l
      h = h/2.0_dp
      h0 = h/2.0_dp
      sum2 = w1
      inf2 = .false.
      goto 60
220   quadr = -h*(sum1+cor+sumv)
      inf = n1
      return
230   quadr = w1
      inf = 2
      m = m2 + 1
      return
240   quadr = sum2
      inf = 3
      m = k + l + 1
      return
250   quadr = w1
      inf = 3
      m = m2/2 + 1
      return
260   quadr = u + cor + sumv
      inf = 4
      m = k + l + 1
      return
270   inf = 10
      return
280   inf = 11
      return
290   inf = 12
      return
   end subroutine toms614_integrate

end module rmutil_integrate
