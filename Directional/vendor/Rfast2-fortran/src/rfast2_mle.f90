module rfast2_mle
   use rfast_special, only : dp, pi, digamma_r, trigamma_r
   use rfast_arrays, only : mean_r, variance_r, median_r
   use rfast_mle, only : mle_result, beta_mle, gamma_mle
   use rfast_extra_mle, only : weibull_mle
   use rfast2_arrays, only : quantile_rfast2
   implicit none
   private

   public :: cauchy0_mle, halfcauchy_mle, kumar_mle, powerlaw_mle, gammapois_mle
   public :: zigamma_mle, zil_mle, ziweibull_mle, simplex_mle, gnormal0_mle
   public :: unitweibull_mle, cbern_mle, sp_mle, truncexp_mle

contains

   function cauchy0_mle(x,tol,maxiter) result(res)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxiter
      type(mle_result) :: res
      real(dp) :: eps,s,logs,der,der2,oldll,newll
      real(dp), allocatable :: down(:),x2(:)
      integer :: n,it,mi

      n = size(x)
      eps = 1.0e-7_dp
      if (present(tol)) eps = tol
      mi = 100
      if (present(maxiter)) mi = maxiter
      s = 0.5_dp*(quantile_rfast2(x,0.75_dp)-quantile_rfast2(x,0.25_dp))
      s = max(s,sqrt(tiny(1.0_dp)))
      logs = log(s)
      allocate(x2(n),down(n))
      x2 = x*x
      oldll = -huge(1.0_dp)
      do it = 1, mi
         s = exp(logs)
         down = 1.0_dp/(x2+s*s)
         newll = real(n,dp)*logs+sum(log(down))-real(n,dp)*log(pi)
         der = real(n,dp)-2.0_dp*s*s*sum(down)
         der2 = -4.0_dp*s**4*sum(down*down)
         if (abs(der2) <= tiny(1.0_dp)) exit
         if (it > 1 .and. abs(newll-oldll) <= eps) exit
         oldll = newll
         logs = logs-der/der2
      end do
      s = exp(logs)
      down = 1.0_dp/(x2+s*s)
      allocate(res%param(1))
      res%param = [s]
      res%loglik = real(n,dp)*log(s)+sum(log(down))-real(n,dp)*log(pi)
      res%iters = it
   end function cauchy0_mle

   function halfcauchy_mle(x,tol,maxiter) result(res)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxiter
      type(mle_result) :: res

      if (any(x < 0.0_dp)) then
         res%status = 1
         return
      end if
      if (present(tol) .and. present(maxiter)) then
         res = cauchy0_mle(x,tol,maxiter)
      else if (present(tol)) then
         res = cauchy0_mle(x,tol=tol)
      else if (present(maxiter)) then
         res = cauchy0_mle(x,maxiter=maxiter)
      else
         res = cauchy0_mle(x)
      end if
      if (res%status == 0) res%loglik = res%loglik + real(size(x),dp)*log(2.0_dp)
   end function halfcauchy_mle

   function powerlaw_mle(x) result(res)
      real(dp), intent(in) :: x(:)
      type(mle_result) :: res
      real(dp) :: xmin,com,a
      integer :: n

      n = size(x)
      if (n == 0 .or. any(x <= 0.0_dp)) then
         res%status = 1
         return
      end if
      xmin = minval(x)
      com = sum(log(x))-real(n,dp)*log(xmin)
      a = 1.0_dp+real(n,dp)/com
      allocate(res%param(2))
      res%param = [a,xmin]
      res%loglik = real(n,dp)*log((a-1.0_dp)/xmin)-a*com
      res%iters = 1
   end function powerlaw_mle

   function sp_mle(x) result(res)
      real(dp), intent(in) :: x(:)
      type(mle_result) :: res
      real(dp) :: slx,b
      integer :: n

      n = size(x)
      if (n == 0 .or. any(x <= 0.0_dp) .or. any(x > 1.0_dp)) then
         res%status = 1
         return
      end if
      slx = sum(log(x))
      b = -real(n,dp)/slx
      allocate(res%param(1))
      res%param = [b]
      res%loglik = real(n,dp)*log(b)+(b-1.0_dp)*slx
      res%iters = 1
   end function sp_mle

   function kumar_mle(x,tol,maxiter) result(res)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxiter
      type(mle_result) :: res
      type(mle_result) :: ini
      real(dp), allocatable :: lx(:),xa(:),ya(:),com(:)
      real(dp) :: eps,a,b,slx,dera,derb,dera2,derb2,derab,det
      real(dp) :: old(2),new(2),scom
      integer :: n,it,mi

      if (any(x <= 0.0_dp) .or. any(x >= 1.0_dp)) then
         res%status = 1
         return
      end if
      n = size(x)
      eps = 1.0e-7_dp
      if (present(tol)) eps = tol
      mi = 50
      if (present(maxiter)) mi = maxiter
      ini = beta_mle(x)
      if (ini%status /= 0) then
         res%status = ini%status
         return
      end if
      a = ini%param(1)
      b = ini%param(2)
      old = log([a,b])
      allocate(lx(n),xa(n),ya(n),com(n))
      lx = log(x)
      slx = sum(lx)
      do it = 1, mi
         a = exp(old(1))
         b = exp(old(2))
         xa = x**a
         ya = max(tiny(1.0_dp),1.0_dp-xa)
         com = xa*lx/ya
         scom = sum(com)
         derab = -b*a*scom
         dera = real(n,dp)+a*slx+(1.0_dp-1.0_dp/b)*derab
         dera2 = a*slx-(b-1.0_dp)*a*a*sum(com*lx/ya)
         derb2 = b*sum(log(ya))
         derb = real(n,dp)+derb2
         det = dera2*derb2-derab*derab
         if (abs(det) <= tiny(1.0_dp)) then
            res%status = 2
            exit
         end if
         new(1) = old(1)-(derb2*dera-derab*derb)/det
         new(2) = old(2)-(-derab*dera+dera2*derb)/det
         if (sum(abs(new-old)) <= eps) then
            old = new
            exit
         end if
         old = new
      end do
      a = exp(old(1))
      b = exp(old(2))
      xa = x**a
      ya = max(tiny(1.0_dp),1.0_dp-xa)
      allocate(res%param(2))
      res%param = [a,b]
      res%loglik = real(n,dp)*log(a*b)+(a-1.0_dp)*slx+(b-1.0_dp)*sum(log(ya))
      res%iters = it
   end function kumar_mle

   function gammapois_mle(x,tol,maxiter) result(res)
      integer, intent(in) :: x(:)
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxiter
      type(mle_result) :: res
      real(dp) :: eps,m,m2,p,a,b,la,lb,da,db,daa,dbb,dab,det,llold,llnew
      integer :: n,it,mi

      n = size(x)
      if (n == 0 .or. any(x < 0)) then
         res%status = 1
         return
      end if
      eps = 1.0e-7_dp
      if (present(tol)) eps = tol
      mi = 100
      if (present(maxiter)) mi = maxiter
      m = real(sum(x),dp)/real(n,dp)
      m2 = sum(real(x,dp)**2)/real(n,dp)
      if (m <= 0.0_dp .or. m2-m*m <= m) then
         a = max(1.0_dp,m*100.0_dp+1.0_dp)
      else
         p = max(1.0e-8_dp,1.0_dp-m/(m2-m*m))
         a = max(1.0e-8_dp,m/p-m)
      end if
      b = max(1.0e-8_dp,a/max(m,1.0e-8_dp))
      la = log(a)
      lb = log(b)
      llold = -huge(1.0_dp)
      do it = 1, mi
         a = exp(la)
         b = exp(lb)
         call gammapois_derivatives(x,a,b,da,db,daa,dbb,dab,llnew)
         det = daa*dbb-dab*dab
         if (abs(det) <= tiny(1.0_dp)) then
            res%status = 2
            exit
         end if
         if (it > 1 .and. abs(llnew-llold) <= eps) exit
         llold = llnew
         la = la-(dbb*da-dab*db)/det
         lb = lb-(-dab*da+daa*db)/det
         la = max(-30.0_dp,min(30.0_dp,la))
         lb = max(-30.0_dp,min(30.0_dp,lb))
      end do
      a = exp(la)
      b = exp(lb)
      call gammapois_derivatives(x,a,b,da,db,daa,dbb,dab,llnew)
      allocate(res%param(2))
      res%param = [a,b]
      res%loglik = llnew
      res%iters = it
   end function gammapois_mle

   subroutine gammapois_derivatives(x,a,b,da,db,daa,dbb,dab,ll)
      integer, intent(in) :: x(:)
      real(dp), intent(in) :: a,b
      real(dp), intent(out) :: da,db,daa,dbb,dab,ll
      real(dp) :: sx,p
      integer :: n,i

      n = size(x)
      sx = real(sum(x),dp)
      p = b/(1.0_dp+b)
      ll = sx*log(p)-real(n,dp)*a*log(1.0_dp+b)
      da = -real(n,dp)*digamma_r(a)-real(n,dp)*log(1.0_dp+b)
      daa = -real(n,dp)*trigamma_r(a)
      do i = 1, n
         ll = ll+log_gamma(real(x(i),dp)+a)-log_gamma(real(x(i)+1,dp))-log_gamma(a)
         da = da+digamma_r(real(x(i),dp)+a)
         daa = daa+trigamma_r(real(x(i),dp)+a)
      end do
      da = a*da
      db = sx*(1.0_dp-p)-real(n,dp)*a*p
      dab = -real(n,dp)*a*p
      daa = da+a*a*daa
      dbb = -p*(1.0_dp-p)*(sx+real(n,dp)*a)
   end subroutine gammapois_derivatives

   function zigamma_mle(x,tol) result(res)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in), optional :: tol
      type(mle_result) :: res
      type(mle_result) :: g
      real(dp), allocatable :: xp(:)
      real(dp) :: prob,ll0
      integer :: n,n1,n0

      n = size(x)
      xp = pack(x,x>0.0_dp)
      n1 = size(xp)
      n0 = n-n1
      if (n1 == 0) then
         res%status = 1
         return
      end if
      prob = real(n1,dp)/real(n,dp)
      ll0 = 0.0_dp
      if (n0 > 0) ll0 = ll0+real(n0,dp)*log(1.0_dp-prob)
      if (n1 > 0) ll0 = ll0+real(n1,dp)*log(prob)
      if (present(tol)) then
         g = gamma_mle(xp,tol)
      else
         g = gamma_mle(xp)
      end if
      allocate(res%param(3))
      res%param = [prob,g%param]
      res%loglik = ll0+g%loglik
      res%iters = g%iters
      res%status = g%status
   end function zigamma_mle

   function zil_mle(x) result(res)
      real(dp), intent(in) :: x(:)
      type(mle_result) :: res
      real(dp), allocatable :: xp(:),y(:)
      real(dp) :: prob,m,v,ll0
      integer :: n,n1,n0

      n = size(x)
      xp = pack(x,x>0.0_dp)
      n1 = size(xp)
      n0 = n-n1
      if (n1 < 2 .or. any(xp >= 1.0_dp)) then
         res%status = 1
         return
      end if
      prob = real(n1,dp)/real(n,dp)
      allocate(y(n1))
      y = log(xp)-log(1.0_dp-xp)
      m = mean_r(y)
      v = sum((y-m)**2)/real(n1,dp)
      ll0 = 0.0_dp
      if (n0 > 0) ll0 = ll0+real(n0,dp)*log(1.0_dp-prob)
      ll0 = ll0+real(n1,dp)*log(prob)
      allocate(res%param(3))
      res%param = [prob,m,v*real(n1,dp)/real(n1-1,dp)]
      res%loglik = ll0-0.5_dp*real(n1,dp)*(log(2.0_dp*pi*v)+1.0_dp)-sum(y)
      res%iters = 1
   end function zil_mle

   function ziweibull_mle(x,tol) result(res)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in), optional :: tol
      type(mle_result) :: res
      type(mle_result) :: w
      real(dp), allocatable :: xp(:)
      real(dp) :: prob,ll0
      integer :: n,n1,n0

      n = size(x)
      xp = pack(x,x>0.0_dp)
      n1 = size(xp)
      n0 = n-n1
      if (n1 == 0) then
         res%status = 1
         return
      end if
      prob = real(n1,dp)/real(n,dp)
      ll0 = real(n1,dp)*log(prob)
      if (n0 > 0) ll0 = ll0+real(n0,dp)*log(1.0_dp-prob)
      if (present(tol)) then
         w = weibull_mle(xp,tol)
      else
         w = weibull_mle(xp)
      end if
      allocate(res%param(3))
      res%param = [prob,w%param]
      res%loglik = ll0+w%loglik
      res%iters = w%iters
      res%status = w%status
   end function ziweibull_mle

   function simplex_mle(x,tol) result(res)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in), optional :: tol
      type(mle_result) :: res
      real(dp) :: a,b,c,d,fc,fd,eps,m,s,obj
      integer :: it,n

      n = size(x)
      if (any(x <= 0.0_dp) .or. any(x >= 1.0_dp)) then
         res%status = 1
         return
      end if
      eps = 1.0e-7_dp
      if (present(tol)) eps = tol
      a = max(1.0e-10_dp,minval(x)*0.5_dp)
      b = min(1.0_dp-1.0e-10_dp,0.5_dp*(1.0_dp+maxval(x)))
      c = b-(b-a)/1.6180339887498948482_dp
      d = a+(b-a)/1.6180339887498948482_dp
      fc = simplex_obj(c,x)
      fd = simplex_obj(d,x)
      do it = 1, 300
         if (abs(b-a) <= eps) exit
         if (fc < fd) then
            b = d
            d = c
            fd = fc
            c = b-(b-a)/1.6180339887498948482_dp
            fc = simplex_obj(c,x)
         else
            a = c
            c = d
            fc = fd
            d = a+(b-a)/1.6180339887498948482_dp
            fd = simplex_obj(d,x)
         end if
      end do
      m = 0.5_dp*(a+b)
      obj = simplex_obj(m,x)
      s = sqrt(obj/real(n,dp))
      allocate(res%param(2))
      res%param = [m,s]
      res%loglik = -0.5_dp*real(n,dp)*log(2.0_dp*pi)-1.5_dp*sum(log(x*(1.0_dp-x))) &
                   -real(n,dp)*log(s)-0.5_dp*real(n,dp)
      res%iters = it
   end function simplex_mle

   real(dp) function simplex_obj(m,x) result(v)
      real(dp), intent(in) :: m,x(:)
      v = sum((x-m)**2/(x*(1.0_dp-x)))/(m*m*(1.0_dp-m)**2)
   end function simplex_obj

   function gnormal0_mle(x,tol,maxiter) result(res)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxiter
      type(mle_result) :: res
      real(dp) :: lo,hi,mid,flo,fmid,eps,b,a
      integer :: it,mi,n

      n = size(x)
      if (n == 0 .or. all(abs(x) <= tiny(1.0_dp))) then
         res%status = 1
         return
      end if
      eps = 1.0e-6_dp
      if (present(tol)) eps = tol
      mi = 200
      if (present(maxiter)) mi = maxiter
      lo = 1.0e-4_dp
      hi = 30.0_dp
      flo = gnormal_root(lo,x)
      mid = 0.5_dp*(lo+hi)
      do it = 1, mi
         mid = 0.5_dp*(lo+hi)
         fmid = gnormal_root(mid,x)
         if (abs(fmid) <= eps .or. abs(hi-lo) <= eps) exit
         if (flo*fmid <= 0.0_dp) then
            hi = mid
         else
            lo = mid
            flo = fmid
         end if
      end do
      b = mid
      a = (b/real(n,dp)*sum(abs(x)**b))**(1.0_dp/b)
      allocate(res%param(2))
      res%param = [a,b]
      res%loglik = real(n,dp)*log(b)-sum(abs(x)**b)/a**b-real(n,dp)*log(2.0_dp*a)-real(n,dp)*log_gamma(1.0_dp/b)
      res%iters = it
   end function gnormal0_mle

   real(dp) function gnormal_root(b,x) result(f)
      real(dp), intent(in) :: b,x(:)
      real(dp) :: y(size(x)),sy

      y = abs(x)**b
      sy = sum(y)
      f = 1.0_dp+digamma_r(1.0_dp/b)/b-sum(y*log(max(abs(x),tiny(1.0_dp))))/sy &
          +log(b*sy/real(size(x),dp))/b
   end function gnormal_root

   function unitweibull_mle(x,tol,maxiter) result(res)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxiter
      type(mle_result) :: res,w
      real(dp), allocatable :: lx(:)
      real(dp) :: a,b
      integer :: n

      if (any(x <= 0.0_dp) .or. any(x >= 1.0_dp)) then
         res%status = 1
         return
      end if
      lx = -log(x)
      if (present(tol) .and. present(maxiter)) then
         w = weibull_mle(lx,tol,maxiter)
      else if (present(tol)) then
         w = weibull_mle(lx,tol)
      else if (present(maxiter)) then
         w = weibull_mle(lx,maxiter=maxiter)
      else
         w = weibull_mle(lx)
      end if
      res = w
      if (res%status /= 0) return
      a = res%param(1)
      b = res%param(2)
      n = size(x)
      res%loglik = sum(lx)+real(n,dp)*log(a*b)+(b-1.0_dp)*sum(log(lx))-a*sum(lx**b)
   end function unitweibull_mle

   function cbern_mle(x,tol) result(res)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in), optional :: tol
      type(mle_result) :: res
      real(dp) :: a,b,c,d,fc,fd,eps,lam
      integer :: it

      if (any(x < 0.0_dp) .or. any(x > 1.0_dp)) then
         res%status = 1
         return
      end if
      eps = 1.0e-6_dp
      if (present(tol)) eps = tol
      a = 1.0e-10_dp
      b = 1.0_dp-1.0e-10_dp
      c = b-(b-a)/1.6180339887498948482_dp
      d = a+(b-a)/1.6180339887498948482_dp
      fc = cbern_ll(c,x)
      fd = cbern_ll(d,x)
      do it = 1, 300
         if (abs(b-a) <= eps) exit
         if (fc > fd) then
            b = d
            d = c
            fd = fc
            c = b-(b-a)/1.6180339887498948482_dp
            fc = cbern_ll(c,x)
         else
            a = c
            c = d
            fc = fd
            d = a+(b-a)/1.6180339887498948482_dp
            fd = cbern_ll(d,x)
         end if
      end do
      lam = 0.5_dp*(a+b)
      allocate(res%param(1))
      res%param = [lam]
      res%loglik = real(size(x),dp)*log(2.0_dp)+cbern_ll(lam,x)
      res%iters = it
   end function cbern_mle

   real(dp) function cbern_ll(lam,x) result(v)
      real(dp), intent(in) :: lam,x(:)
      real(dp) :: c,delta

      delta = 1.0_dp-2.0_dp*lam
      if (abs(delta) < 1.0e-8_dp) then
         c = 0.0_dp
      else
         c = real(size(x),dp)*(log(abs(atanh(delta)))-log(abs(delta)))
      end if
      v = c+sum(x)*log(lam)+sum(1.0_dp-x)*log(1.0_dp-lam)
   end function cbern_ll

   function truncexp_mle(x,bound,tol) result(res)
      real(dp), intent(in) :: x(:),bound
      real(dp), intent(in), optional :: tol
      type(mle_result) :: res
      real(dp) :: a,b,c,d,fc,fd,eps,lam
      integer :: it

      eps = 1.0e-7_dp
      if (present(tol)) eps = tol
      a = max(tiny(1.0_dp),bound*1.0e-8_dp)
      b = bound
      c = b-(b-a)/1.6180339887498948482_dp
      d = a+(b-a)/1.6180339887498948482_dp
      fc = truncexp_ll(c,x,bound)
      fd = truncexp_ll(d,x,bound)
      do it = 1, 300
         if (abs(b-a) <= eps) exit
         if (fc > fd) then
            b = d
            d = c
            fd = fc
            c = b-(b-a)/1.6180339887498948482_dp
            fc = truncexp_ll(c,x,bound)
         else
            a = c
            c = d
            fc = fd
            d = a+(b-a)/1.6180339887498948482_dp
            fd = truncexp_ll(d,x,bound)
         end if
      end do
      lam = 0.5_dp*(a+b)
      allocate(res%param(1))
      res%param = [lam]
      res%loglik = truncexp_ll(lam,x,bound)
      res%iters = it
   end function truncexp_mle

   real(dp) function truncexp_ll(lam,x,bound) result(v)
      real(dp), intent(in) :: lam,x(:),bound
      v = -real(size(x),dp)*log(lam)-sum(x)/lam-real(size(x),dp)*log(1.0_dp-exp(-bound/lam))
   end function truncexp_ll

end module rfast2_mle
