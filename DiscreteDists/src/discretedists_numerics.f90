module discretedists_numerics
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_is_finite
   use discretedists_kinds, only : dp
   implicit none
   private

   real(dp), parameter :: pi = acos(-1.0_dp)

   type, public :: sum_result_t
      real(dp) :: value = 0.0_dp
      real(dp) :: abs_error = 0.0_dp
      integer :: terms = 0
      logical :: converged = .true.
   end type sum_result_t

   abstract interface
      function scalar_series_fun(x) result(y)
         import dp
         real(dp), intent(in) :: x
         real(dp) :: y
      end function scalar_series_fun
      function objective_fun(theta, context) result(value)
         import dp
         real(dp), intent(in) :: theta(:)
         class(*), intent(in) :: context
         real(dp) :: value
      end function objective_fun
   end interface

   public :: add_series, stopping, lambert_wm1, normal_cdf, normal_quantile
   public :: logistic, logit, nelder_mead, mean_sample, variance_sample
   public :: finite_difference_1, finite_difference_2, objective_fun

contains

   pure real(dp) function mean_sample(x) result(m)
      real(dp), intent(in) :: x(:)
      if (size(x) == 0) then
         m = ieee_value(0.0_dp, ieee_quiet_nan)
      else
         m = sum(x)/real(size(x),dp)
      end if
   end function mean_sample

   pure real(dp) function variance_sample(x) result(v)
      real(dp), intent(in) :: x(:)
      real(dp) :: m
      if (size(x) < 2) then
         v = 0.0_dp
      else
         m = sum(x)/real(size(x),dp)
         v = sum((x-m)**2)/real(size(x)-1,dp)
      end if
   end function variance_sample

   pure logical function stopping(x,tol) result(done)
      real(dp), intent(in) :: x(:),tol
      done = all(abs(x) <= tol)
   end function stopping

   function add_series(f, lower, upper, abs_tol, max_terms) result(ans)
      procedure(scalar_series_fun) :: f
      real(dp), intent(in) :: lower, upper
      real(dp), intent(in), optional :: abs_tol
      integer, intent(in), optional :: max_terms
      type(sum_result_t) :: ans
      real(dp) :: tol,x,t
      integer :: k,m,n
      tol = epsilon(1.0_dp)
      if (present(abs_tol)) tol=max(abs_tol,0.0_dp)
      m=1000000
      if (present(max_terms)) m=max(1,max_terms)
      if (lower >= upper) then
         ans%value=ieee_value(0.0_dp,ieee_quiet_nan);ans%converged=.false.;return
      end if
      ans%value=0.0_dp;ans%abs_error=0.0_dp;ans%terms=0;ans%converged=.true.
      if (ieee_is_finite(lower) .and. ieee_is_finite(upper)) then
         n=max(0,int(floor(upper-lower))+1)
         do k=0,n-1
            ans%value=ans%value+f(lower+real(k,dp))
         end do
         ans%terms=n
         return
      end if
      if (ieee_is_finite(lower)) then
         do k=0,min(300,m-1)
            ans%value=ans%value+f(lower+real(k,dp));ans%terms=ans%terms+1
         end do
         x=lower+real(ans%terms,dp)
         do while(ans%terms<m)
            t=f(x);ans%value=ans%value+t;ans%terms=ans%terms+1;ans%abs_error=abs(t)
            if (abs(t)<tol) return
            x=x+1.0_dp
         end do
         ans%converged=.false.;return
      end if
      if (ieee_is_finite(upper)) then
         do k=0,min(300,m-1)
            ans%value=ans%value+f(upper-real(k,dp));ans%terms=ans%terms+1
         end do
         x=upper-real(ans%terms,dp)
         do while(ans%terms<m)
            t=f(x);ans%value=ans%value+t;ans%terms=ans%terms+1;ans%abs_error=abs(t)
            if (abs(t)<tol) return
            x=x-1.0_dp
         end do
         ans%converged=.false.;return
      end if
      ! Match the upstream strategy: seed with -100:100, then add symmetric tails.
      do k=-100,100
         ans%value=ans%value+f(real(k,dp));ans%terms=ans%terms+1
      end do
      x=101.0_dp
      do while(ans%terms+2<=m)
         t=f(x)+f(-x);ans%value=ans%value+t;ans%terms=ans%terms+2;ans%abs_error=abs(t)
         if (abs(t)<tol) return
         x=x+1.0_dp
      end do
      ans%converged=.false.
   end function add_series

   pure real(dp) function log1p_local(x) result(v)
      real(dp), intent(in) :: x
      if (abs(x) < 1.0e-8_dp) then
         v=x-x*x/2.0_dp+x**3/3.0_dp-x**4/4.0_dp+x**5/5.0_dp
      else
         v=log(1.0_dp+x)
      end if
   end function log1p_local

   pure real(dp) function logistic(x) result(p)
      real(dp), intent(in) :: x
      if (x >= 0.0_dp) then
         p=1.0_dp/(1.0_dp+exp(-x))
      else
         p=exp(x)/(1.0_dp+exp(x))
      end if
   end function logistic

   pure real(dp) function logit(p) result(x)
      real(dp), intent(in) :: p
      x=log(p)-log1p_local(-p)
   end function logit

   pure real(dp) function normal_cdf(x) result(p)
      real(dp), intent(in) :: x
      p=0.5_dp*erfc(-x/sqrt(2.0_dp))
   end function normal_cdf

   pure real(dp) function normal_quantile(p) result(x)
      real(dp), intent(in) :: p
      real(dp), parameter :: a(6) = [ &
         -3.969683028665376e+01_dp, 2.209460984245205e+02_dp, -2.759285104469687e+02_dp, &
          1.383577518672690e+02_dp,-3.066479806614716e+01_dp, 2.506628277459239e+00_dp ]
      real(dp), parameter :: b(5) = [ &
         -5.447609879822406e+01_dp, 1.615858368580409e+02_dp, -1.556989798598866e+02_dp, &
          6.680131188771972e+01_dp,-1.328068155288572e+01_dp ]
      real(dp), parameter :: c(6) = [ &
         -7.784894002430293e-03_dp,-3.223964580411365e-01_dp,-2.400758277161838e+00_dp, &
         -2.549732539343734e+00_dp, 4.374664141464968e+00_dp, 2.938163982698783e+00_dp ]
      real(dp), parameter :: d(4) = [ &
          7.784695709041462e-03_dp, 3.224671290700398e-01_dp, 2.445134137142996e+00_dp, &
          3.754408661907416e+00_dp ]
      real(dp), parameter :: plow=0.02425_dp, phigh=1.0_dp-plow
      real(dp) :: q,r,e,u
      integer :: it
      if (p <= 0.0_dp) then
         x=-ieee_value(0.0_dp,ieee_positive_inf);return
      else if (p >= 1.0_dp) then
         x=ieee_value(0.0_dp,ieee_positive_inf);return
      else if (p < plow) then
         q=sqrt(-2.0_dp*log(p))
         x=(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
           ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
      else if (p > phigh) then
         q=sqrt(-2.0_dp*log1p_local(-p))
         x=-(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
           ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
      else
         q=p-0.5_dp;r=q*q
         x=(((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q/ &
           (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
      end if
      ! Two Newton refinements substantially improve the Acklam seed.
      do it=1,2
         e=normal_cdf(x)-p
         u=e*sqrt(2.0_dp*pi)*exp(0.5_dp*x*x)
         x=x-u/(1.0_dp+0.5_dp*x*u)
      end do
   end function normal_quantile

   function lambert_wm1(x) result(w)
      real(dp), intent(in) :: x
      real(dp) :: w,em1,q,ew,f,den,dw
      integer :: it
      em1=-exp(-1.0_dp)
      if (x < em1 .or. x >= 0.0_dp) then
         w=ieee_value(0.0_dp,ieee_quiet_nan);return
      end if
      if (abs(x-em1) <= 8.0_dp*epsilon(1.0_dp)) then
         w=-1.0_dp;return
      end if
      if (x < -0.05_dp) then
         q=sqrt(max(0.0_dp,2.0_dp*(1.0_dp+exp(1.0_dp)*x)))
         w=-1.0_dp-q-q*q/3.0_dp-11.0_dp*q**3/72.0_dp
      else
         w=log(-x)-log(-log(-x))
         if (w>-1.0_dp) w=-1.1_dp
      end if
      do it=1,100
         ew=exp(w);f=w*ew-x
         den=ew*(w+1.0_dp)-(w+2.0_dp)*f/(2.0_dp*w+2.0_dp)
         if (den==0.0_dp) exit
         dw=f/den;w=w-dw
         if (abs(dw)<=8.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(w))) exit
         if (w>-1.0_dp) w=-1.0_dp-epsilon(1.0_dp)
      end do
   end function lambert_wm1

   subroutine nelder_mead(fun,context,start,best,fbest,status,iterations,max_iter,tol,step)
      procedure(objective_fun) :: fun
      class(*), intent(in) :: context
      real(dp), intent(in) :: start(:)
      real(dp), intent(out) :: best(size(start)),fbest
      integer, intent(out) :: status,iterations
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: tol,step
      integer :: n,m,it,i,j,ilo,ihi,inhi,mit
      real(dp) :: ftol,delta,fr,fe,fc,spread,diam
      real(dp), allocatable :: x(:,:),fv(:),cent(:),xr(:),xe(:),xc(:)
      n=size(start);m=n+1;mit=5000;if(present(max_iter))mit=max_iter
      ftol=1.0e-9_dp;if(present(tol))ftol=tol
      delta=0.2_dp;if(present(step))delta=step
      allocate(x(n,m),fv(m),cent(n),xr(n),xe(n),xc(n))
      x(:,1)=start
      do j=2,m
         x(:,j)=start
         x(j-1,j)=x(j-1,j)+delta*max(1.0_dp,abs(start(j-1)))
      end do
      do j=1,m;fv(j)=fun(x(:,j),context);end do
      status=1
      do it=1,mit
         ilo=minloc(fv,dim=1);ihi=maxloc(fv,dim=1)
         inhi=1
         do j=1,m
            if(j==ihi)cycle
            if(inhi==ihi .or. fv(j)>fv(inhi))inhi=j
         end do
         spread=maxval(abs(fv-fv(ilo)))
         diam=0.0_dp
         do j=1,m;diam=max(diam,maxval(abs(x(:,j)-x(:,ilo))));end do
         if(spread<=ftol*(1.0_dp+abs(fv(ilo))) .and. diam<=sqrt(ftol)*(1.0_dp+maxval(abs(x(:,ilo)))))then
            status=0;exit
         end if
         cent=0.0_dp
         do j=1,m;if(j/=ihi)cent=cent+x(:,j);end do
         cent=cent/real(n,dp)
         xr=cent+(cent-x(:,ihi));fr=fun(xr,context)
         if(fr<fv(ilo))then
            xe=cent+2.0_dp*(xr-cent);fe=fun(xe,context)
            if(fe<fr)then;x(:,ihi)=xe;fv(ihi)=fe;else;x(:,ihi)=xr;fv(ihi)=fr;end if
         else if(fr<fv(inhi))then
            x(:,ihi)=xr;fv(ihi)=fr
         else
            if(fr<fv(ihi))then
               xc=cent+0.5_dp*(xr-cent)
            else
               xc=cent+0.5_dp*(x(:,ihi)-cent)
            end if
            fc=fun(xc,context)
            if(fc<min(fr,fv(ihi)))then
               x(:,ihi)=xc;fv(ihi)=fc
            else
               do j=1,m
                  if(j==ilo)cycle
                  x(:,j)=x(:,ilo)+0.5_dp*(x(:,j)-x(:,ilo));fv(j)=fun(x(:,j),context)
               end do
            end if
         end if
      end do
      iterations=min(it,mit);ilo=minloc(fv,dim=1);best=x(:,ilo);fbest=fv(ilo)
   end subroutine nelder_mead

   function finite_difference_1(fun,x,h) result(d)
      interface
         function fun(z) result(v)
            import dp
            real(dp),intent(in)::z
            real(dp)::v
         end function fun
      end interface
      real(dp),intent(in)::x,h
      real(dp)::d
      d=(fun(x+h)-fun(x-h))/(2.0_dp*h)
   end function finite_difference_1

   function finite_difference_2(fun,x,h) result(d2)
      interface
         function fun(z) result(v)
            import dp
            real(dp),intent(in)::z
            real(dp)::v
         end function fun
      end interface
      real(dp),intent(in)::x,h
      real(dp)::d2
      d2=(fun(x+h)-2.0_dp*fun(x)+fun(x-h))/(h*h)
   end function finite_difference_2

end module discretedists_numerics
