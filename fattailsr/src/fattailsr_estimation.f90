! FatTailsR modern Fortran translation
! Copyright (C) 2014-2026 Patrice Kiener
! Licensed under GPL-2.0-only. See COPYING.
module fattailsr_estimation
   use fattailsr_kinds, only : dp
   use fattailsr_math, only : pi, sqrt3, logit
   use fattailsr_params, only : kiener_parameters, make_k4, kd2a, kd2w, kd2e
   use fattailsr_distributions, only : qlkiener
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
   implicit none
   private

   public :: five_probs, seven_probs, eleven_probs, quantile_type6
   public :: estimkiener5, estimkiener7, estimkiener11, fit_kiener_k4

contains

   pure subroutine sort_real(x)
      real(dp), intent(inout) :: x(:)
      integer :: i, j
      real(dp) :: key
      do i = 2, size(x)
         key = x(i)
         j = i - 1
         do while (j >= 1)
            if (x(j) <= key) exit
            x(j+1) = x(j)
            j = j - 1
         end do
         x(j+1) = key
      end do
   end subroutine sort_real

   pure function probability_grid_value(n, offset) result(p)
      integer, intent(in) :: n, offset
      real(dp) :: p
      real(dp), parameter :: listp(26) = [ &
         0.25e-8_dp, 0.50e-8_dp, 1.00e-8_dp, &
         0.25e-7_dp, 0.50e-7_dp, 1.00e-7_dp, &
         0.25e-6_dp, 0.50e-6_dp, 1.00e-6_dp, &
         0.25e-5_dp, 0.50e-5_dp, 1.00e-5_dp, &
         0.25e-4_dp, 0.50e-4_dp, 1.00e-4_dp, &
         0.25e-3_dp, 0.50e-3_dp, 1.00e-3_dp, &
         0.25e-2_dp, 0.50e-2_dp, 1.00e-2_dp, &
         0.25e-1_dp, 0.50e-1_dp, 1.00e-1_dp, 0.15_dp, 0.20_dp ]
      integer :: i, idx
      idx = 0
      do i = 1, size(listp)
         if (listp(i) <= 1.0_dp/real(n,dp)) idx = i
      end do
      idx = min(idx + offset, size(listp))
      p = listp(idx)
   end function probability_grid_value

   pure function five_probs(n, i_extreme) result(p)
      integer, intent(in) :: n
      integer, intent(in), optional :: i_extreme
      real(dp) :: p(5)
      integer :: i
      i = 4
      if (present(i_extreme)) i = i_extreme
      if (i == 0) then
         p = [probability_grid_value(n,1), 0.25_dp, 0.50_dp, 0.75_dp, &
              1.0_dp-probability_grid_value(n,1)]
      else
         p = [real(i,dp)/real(n+1,dp), 0.25_dp, 0.50_dp, 0.75_dp, &
              real(n+1-i,dp)/real(n+1,dp)]
      end if
   end function five_probs

   pure function seven_probs(n) result(p)
      integer, intent(in) :: n
      real(dp) :: p(7), p1, p2
      p1 = probability_grid_value(n,1)
      p2 = probability_grid_value(n,2)
      p = [p1, p2, 0.25_dp, 0.50_dp, 0.75_dp, 1.0_dp-p2, 1.0_dp-p1]
   end function seven_probs

   pure function eleven_probs(n) result(p)
      integer, intent(in) :: n
      real(dp) :: p(11), p1, p2, p3
      p1 = probability_grid_value(n,1)
      p2 = probability_grid_value(n,2)
      p3 = probability_grid_value(n,3)
      p = [p1,p2,p3,0.25_dp,0.35_dp,0.50_dp,0.65_dp,0.75_dp,&
           1.0_dp-p3,1.0_dp-p2,1.0_dp-p1]
   end function eleven_probs

   pure function quantile_type6(x, p) result(q)
      real(dp), intent(in) :: x(:), p
      real(dp) :: q
      real(dp), allocatable :: y(:)
      real(dp) :: h, frac
      integer :: lo
      allocate(y(size(x)))
      y = x
      call sort_real(y)
      h = real(size(y)+1,dp)*p
      if (h <= 1.0_dp) then
         q = y(1)
      else if (h >= real(size(y),dp)) then
         q = y(size(y))
      else
         lo = int(floor(h))
         frac = h - real(lo,dp)
         q = y(lo) + frac*(y(lo+1)-y(lo))
      end if
   end function quantile_type6

   pure function objective_k(k, lq, lp, target) result(v)
      real(dp), intent(in) :: k, lq, lp, target
      real(dp) :: v
      v = (target - sinh(lp/k)/sinh(lq/k))**2
   end function objective_k

   pure function golden_k(lq, lp, target, lo, hi) result(kbest)
      real(dp), intent(in) :: lq, lp, target, lo, hi
      real(dp) :: kbest
      real(dp), parameter :: gr = 0.6180339887498948482_dp
      real(dp) :: a, b, c, d, fc, fd
      integer :: iter
      a = lo
      b = hi
      c = b - gr*(b-a)
      d = a + gr*(b-a)
      fc = objective_k(c,lq,lp,target)
      fd = objective_k(d,lq,lp,target)
      do iter = 1, 120
         if (abs(b-a) < 1.0e-6_dp*max(1.0_dp,abs(a)+abs(b))) exit
         if (fc < fd) then
            b=d; d=c; fd=fc; c=b-gr*(b-a); fc=objective_k(c,lq,lp,target)
         else
            a=c; c=d; fc=fd; d=a+gr*(b-a); fd=objective_k(d,lq,lp,target)
         end if
      end do
      kbest = 0.5_dp*(a+b)
   end function golden_k

   pure function invalid_parameters() result(par)
      type(kiener_parameters) :: par
      real(dp) :: nan
      nan = ieee_value(nan, ieee_quiet_nan)
      par%m=nan; par%g=nan; par%a=nan; par%k=nan; par%w=nan; par%d=nan; par%e=nan
   end function invalid_parameters

   pure function estimkiener5(x5, p5, maxk, maxe) result(par)
      real(dp), intent(in) :: x5(5), p5(5)
      real(dp), intent(in), optional :: maxk, maxe
      type(kiener_parameters) :: par
      real(dp) :: km, em, dx(5), lp(5), d, target, k, e, g, threshold
      integer :: i
      km = 20.0_dp
      em = 0.90_dp
      if (present(maxk)) km = maxk
      if (present(maxe)) em = maxe
      if (any(x5(2:5) <= x5(1:4))) then
         par = invalid_parameters()
         return
      end if
      do i=1,5
         dx(i)=abs(x5(i)-x5(3)); lp(i)=logit(p5(i))
      end do
      d = log((x5(5)-x5(3))/(x5(3)-x5(1)))/(2.0_dp*lp(5))
      target = (x5(5)-x5(1))/(x5(4)-x5(2))*cosh(d*lp(4))/cosh(d*lp(5))
      threshold = sinh(lp(5)/km)/sinh(lp(4)/km)
      if (target <= threshold) then
         k = km
      else
         k = golden_k(lp(4),lp(5),target,0.1_dp,km)
      end if
      if (abs(d) > tiny(1.0_dp)) k = min(k,abs(1.0_dp/d)*em)
      e = kd2e(k,d)
      g = (x5(4)-x5(2))*pi/sqrt3/(2.0_dp*k*sinh(lp(4)/k)*cosh(d*lp(4)))
      par = kiener_parameters(m=x5(3),g=g,a=kd2a(k,d),k=k,w=kd2w(k,d),d=d,e=e)
   end function estimkiener5

   pure function estim_kappa6(lg, lp, dx1p, dx1g, dxg, dxp, maxk) result(k)
      real(dp), intent(in) :: lg,lp,dx1p,dx1g,dxg,dxp,maxk
      real(dp) :: k, h, lgh, rss, psi, den
      h = lp/lg
      lgh = lg*sqrt(max(0.0_dp,(-7.0_dp+3.0_dp*h*h)/30.0_dp))
      rss = 1.2_dp - 1.6_dp/(-1.0_dp+h*h)
      psi = sqrt(dx1p*dxp/(dx1g*dxg))/h
      if (psi <= 1.0_dp) then
         k = maxk
      else
         den = -1.0_dp + sqrt(max(0.0_dp,1.0_dp+rss*(-1.0_dp+psi)))
         if (den <= 0.0_dp) then
            k = maxk
         else
            k = lgh/sqrt(den)
         end if
      end if
      if (.not. ieee_is_finite(k)) k=maxk
      k=min(k,maxk)
   end function estim_kappa6

   pure function estimkiener7(x7, p7, maxk) result(par)
      real(dp), intent(in) :: x7(7), p7(7)
      real(dp), intent(in), optional :: maxk
      type(kiener_parameters) :: par
      real(dp) :: km, dx(7), lp(7), k, d, g, e
      integer :: i
      km=10.0_dp
      if (present(maxk)) km=maxk
      if (any(x7(2:7) <= x7(1:6))) then
         par=invalid_parameters(); return
      end if
      do i=1,7
         dx(i)=abs(x7(i)-x7(4)); lp(i)=logit(p7(i))
      end do
      k=(estim_kappa6(lp(5),lp(7),dx(1),dx(3),dx(5),dx(7),km)+&
         estim_kappa6(lp(5),lp(6),dx(2),dx(3),dx(5),dx(6),km))/2.0_dp
      d=log(dx(7)/dx(1))/(4.0_dp*lp(7))+log(dx(6)/dx(2))/(4.0_dp*lp(6))
      if (abs(d)>0.90_dp/k) d=sign(0.90_dp/k,d)
      g=sqrt(dx(3)*dx(5))*pi/sqrt3/(k*sinh(lp(5)/k))
      e=kd2e(k,d)
      par=kiener_parameters(m=x7(4),g=g,a=kd2a(k,d),k=k,w=kd2w(k,d),d=d,e=e)
   end function estimkiener7

   pure function estimkiener11(x11,p11,ord,maxk) result(par)
      real(dp), intent(in) :: x11(11),p11(11)
      integer, intent(in), optional :: ord
      real(dp), intent(in), optional :: maxk
      type(kiener_parameters) :: par
      integer :: o, i
      real(dp) :: km, dx(11),lp(11), kvals(6),k,d1,d2,d3,d,g65,g75,g,e
      o=7; km=10.0_dp
      if(present(ord)) o=ord
      if(present(maxk)) km=maxk
      if(any(x11(2:11)<=x11(1:10)) .or. o<1 .or. o>12) then
         par=invalid_parameters(); return
      end if
      do i=1,11
         dx(i)=abs(x11(i)-x11(6)); lp(i)=logit(p11(i))
      end do
      kvals=[&
        estim_kappa6(lp(7),lp(11),dx(1),dx(5),dx(7),dx(11),km),&
        estim_kappa6(lp(7),lp(10),dx(2),dx(5),dx(7),dx(10),km),&
        estim_kappa6(lp(7),lp(9), dx(3),dx(5),dx(7),dx(9), km),&
        estim_kappa6(lp(8),lp(11),dx(1),dx(4),dx(8),dx(11),km),&
        estim_kappa6(lp(8),lp(10),dx(2),dx(4),dx(8),dx(10),km),&
        estim_kappa6(lp(8),lp(9), dx(3),dx(4),dx(8),dx(9), km)]
      select case(o)
      case(1); k=kvals(1)
      case(2); k=kvals(2)
      case(3); k=sum(kvals(1:2))/2
      case(4); k=sum(kvals(1:3))/3
      case(5); k=kvals(4)
      case(6); k=kvals(5)
      case(7); k=sum(kvals(4:5))/2
      case(8); k=sum(kvals(4:6))/3
      case(9); k=(kvals(1)+kvals(4))/2
      case(10); k=(kvals(2)+kvals(5))/2
      case(11); k=(kvals(1)+kvals(2)+kvals(4)+kvals(5))/4
      case default; k=sum(kvals)/6
      end select
      d1=log(dx(11)/dx(1))/(2*lp(11)); d2=log(dx(10)/dx(2))/(2*lp(10)); d3=log(dx(9)/dx(3))/(2*lp(9))
      select case(o)
      case(1,5,9); d=d1
      case(2,6,10); d=d2
      case(3,7,11); d=(d1+d2)/2
      case default; d=(d1+d2+d3)/3
      end select
      if(abs(d)>0.90_dp/k) d=sign(0.90_dp/k,d)
      g75=sqrt(dx(4)*dx(8))*pi/sqrt3/(k*sinh(lp(8)/k))
      g65=sqrt(dx(5)*dx(7))*pi/sqrt3/(k*sinh(lp(7)/k))
      if(o<=4) then; g=g65
      else if(o<=8) then; g=g75
      else; g=(g65+g75)/2
      end if
      e=kd2e(k,d)
      par=kiener_parameters(m=x11(6),g=g,a=kd2a(k,d),k=k,w=kd2w(k,d),d=d,e=e)
   end function estimkiener11

   pure function ppoint(i,n,a) result(p)
      integer,intent(in)::i,n
      real(dp),intent(in)::a
      real(dp)::p
      p=(real(i,dp)-a)/(real(n,dp)+1.0_dp-2.0_dp*a)
   end function ppoint

   pure function fit_objective(y, m, g, k, e) result(sse)
      real(dp),intent(in)::y(:),m,g,k,e
      real(dp)::sse,lp
      integer::i,n
      type(kiener_parameters)::par
      par=make_k4(m,g,k,e); n=size(y); sse=0.0_dp
      do i=1,n
         lp=logit(ppoint(i,n,0.0_dp))
         sse=sse+(y(i)-qlkiener(lp,par))**2
      end do
      sse=sse/real(n,dp)
   end function fit_objective

   pure function fit_kiener_k4(x,maxk,mink,maxe) result(par)
      real(dp),intent(in)::x(:)
      real(dp),intent(in),optional::maxk,mink,maxe
      type(kiener_parameters)::par
      real(dp),allocatable::y(:)
      real(dp)::km,kn,em,m,g,k,e,best,cand,step(3),trial(3),theta(3)
      real(dp)::p5(5),x5(5)
      integer::i,j,iter,sgn
      km=20.0_dp; kn=1.53_dp; em=0.5_dp
      if(present(maxk))km=maxk; if(present(mink))kn=mink; if(present(maxe))em=maxe
      if(size(x)<15) then; par=invalid_parameters(); return; end if
      allocate(y(size(x))); y=x; call sort_real(y)
      m=quantile_type6(y,0.5_dp); p5=five_probs(size(y),4)
      do i=1,5; x5(i)=quantile_type6(y,p5(i)); end do
      par=estimkiener5(x5,p5,km,min(0.9_dp,em))
      if(.not.ieee_is_finite(par%g)) then
         g=sqrt(sum((y-sum(y)/real(size(y),dp))**2)/real(size(y)-1,dp))*pi/sqrt3/2
         k=4.0_dp; e=0.0_dp
      else
         g=par%g; k=min(max(par%k,kn),km); e=min(max(par%e,-em),em)
      end if
      theta=[g,k,e]; step=[max(0.1_dp*g,1.0e-4_dp),0.5_dp,0.05_dp]
      best=fit_objective(y,m,theta(1),theta(2),theta(3))
      do iter=1,300
         do j=1,3
            do sgn=-1,1,2
               trial=theta; trial(j)=trial(j)+real(sgn,dp)*step(j)
               trial(1)=max(trial(1),1.0e-12_dp)
               trial(2)=min(max(trial(2),kn),km)
               trial(3)=min(max(trial(3),-em),em)
               cand=fit_objective(y,m,trial(1),trial(2),trial(3))
               if(cand<best) then; theta=trial; best=cand; end if
            end do
         end do
         step=step*0.94_dp
         if(maxval(step)<1.0e-7_dp) exit
      end do
      par=make_k4(m,theta(1),theta(2),theta(3))
   end function fit_kiener_k4

end module fattailsr_estimation
