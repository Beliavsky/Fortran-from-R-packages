module suppdists_special
   use suppdists_kinds, only : dp, pi, sqrt2, sqrt2pi
   implicit none
   private
   public :: normal_pdf, normal_cdf, normal_quantile, beta_inc, beta_pdf, beta_quantile
   public :: gamma_p, gamma_q, chisq_pdf, chisq_cdf, randn, rand_gamma, rand_chisq
   public :: adaptive_integral, clamp01

   abstract interface
      function scalar_fun(x) result(y)
         import dp
         real(dp), intent(in) :: x
         real(dp) :: y
      end function scalar_fun
   end interface
contains
   pure real(dp) function clamp01(x) result(y)
      real(dp), intent(in) :: x
      y = max(0.0_dp, min(1.0_dp, x))
   end function clamp01

   pure real(dp) function normal_pdf(x) result(y)
      real(dp), intent(in) :: x
      y = exp(-0.5_dp*x*x)/sqrt2pi
   end function normal_pdf

   pure real(dp) function normal_cdf(x) result(p)
      real(dp), intent(in) :: x
      p = 0.5_dp*erfc(-x/sqrt2)
   end function normal_cdf

   pure real(dp) function normal_quantile(p) result(x)
      real(dp), intent(in) :: p
      real(dp), parameter :: a(6) = [ &
         -3.969683028665376e1_dp, 2.209460984245205e2_dp, &
         -2.759285104469687e2_dp, 1.383577518672690e2_dp, &
         -3.066479806614716e1_dp, 2.506628277459239_dp ]
      real(dp), parameter :: b(5) = [ &
         -5.447609879822406e1_dp, 1.615858368580409e2_dp, &
         -1.556989798598866e2_dp, 6.680131188771972e1_dp, &
         -1.328068155288572e1_dp ]
      real(dp), parameter :: c(6) = [ &
         -7.784894002430293e-3_dp, -3.223964580411365e-1_dp, &
         -2.400758277161838_dp, -2.549732539343734_dp, &
          4.374664141464968_dp, 2.938163982698783_dp ]
      real(dp), parameter :: d(4) = [ &
          7.784695709041462e-3_dp, 3.224671290700398e-1_dp, &
          2.445134137142996_dp, 3.754408661907416_dp ]
      real(dp), parameter :: plow=0.02425_dp, phigh=1.0_dp-plow
      real(dp) :: q, r, e
      if (p <= 0.0_dp) then
         x = -huge(1.0_dp)
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
      else if (p < plow) then
         q = sqrt(-2.0_dp*log(p))
         x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
             ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
      else if (p <= phigh) then
         q = p-0.5_dp
         r = q*q
         x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q / &
             (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
      else
         q = sqrt(-2.0_dp*log(1.0_dp-p))
         x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
              ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
      end if
      if (abs(x) < huge(1.0_dp)/2.0_dp) then
         e = normal_cdf(x)-p
         x = x - e/max(normal_pdf(x),tiny(1.0_dp))
      end if
   end function normal_quantile

   pure real(dp) function betacf(a,b,x) result(cf)
      real(dp), intent(in) :: a,b,x
      integer, parameter :: maxit=300
      real(dp), parameter :: eps=3.0e-14_dp, fpmin=1.0e-300_dp
      real(dp) :: qab,qap,qam,c,d,h,aa,del
      integer :: m,m2
      qab=a+b; qap=a+1.0_dp; qam=a-1.0_dp
      c=1.0_dp
      d=1.0_dp-qab*x/qap
      if (abs(d)<fpmin) d=fpmin
      d=1.0_dp/d; h=d
      do m=1,maxit
         m2=2*m
         aa=real(m,dp)*(b-real(m,dp))*x/((qam+real(m2,dp))*(a+real(m2,dp)))
         d=1.0_dp+aa*d; if (abs(d)<fpmin) d=fpmin
         c=1.0_dp+aa/c; if (abs(c)<fpmin) c=fpmin
         d=1.0_dp/d; h=h*d*c
         aa=-(a+real(m,dp))*(qab+real(m,dp))*x/ &
            ((a+real(m2,dp))*(qap+real(m2,dp)))
         d=1.0_dp+aa*d; if (abs(d)<fpmin) d=fpmin
         c=1.0_dp+aa/c; if (abs(c)<fpmin) c=fpmin
         d=1.0_dp/d; del=d*c; h=h*del
         if (abs(del-1.0_dp) <= eps) exit
      end do
      cf=h
   end function betacf

   pure real(dp) function beta_inc(x,a,b) result(p)
      real(dp), intent(in) :: x,a,b
      real(dp) :: bt
      if (a<=0.0_dp .or. b<=0.0_dp) then
         p = 0.0_dp; return
      end if
      if (x<=0.0_dp) then
         p=0.0_dp; return
      else if (x>=1.0_dp) then
         p=1.0_dp; return
      end if
      bt=exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b)+a*log(x)+b*log(1.0_dp-x))
      if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
         p=bt*betacf(a,b,x)/a
      else
         p=1.0_dp-bt*betacf(b,a,1.0_dp-x)/b
      end if
      p=clamp01(p)
   end function beta_inc

   pure real(dp) function beta_pdf(x,a,b) result(y)
      real(dp), intent(in) :: x,a,b
      if (x<=0.0_dp .or. x>=1.0_dp .or. a<=0.0_dp .or. b<=0.0_dp) then
         y=0.0_dp
      else
         y=exp((a-1.0_dp)*log(x)+(b-1.0_dp)*log(1.0_dp-x) &
            -log_gamma(a)-log_gamma(b)+log_gamma(a+b))
      end if
   end function beta_pdf

   pure real(dp) function beta_quantile(p,a,b) result(x)
      real(dp), intent(in) :: p,a,b
      real(dp) :: lo,hi,mid
      integer :: i
      if (p<=0.0_dp) then; x=0.0_dp; return; end if
      if (p>=1.0_dp) then; x=1.0_dp; return; end if
      lo=0.0_dp; hi=1.0_dp
      do i=1,100
         mid=0.5_dp*(lo+hi)
         if (beta_inc(mid,a,b)<p) then
            lo=mid
         else
            hi=mid
         end if
      end do
      x=0.5_dp*(lo+hi)
   end function beta_quantile

   pure real(dp) function gamma_p(a,x) result(p)
      real(dp), intent(in) :: a,x
      integer, parameter :: itmax=300
      real(dp), parameter :: eps=3e-14_dp, fpmin=1e-300_dp
      real(dp) :: ap,del,sumv,b,c,d,h,an
      integer :: n
      if (a<=0.0_dp .or. x<0.0_dp) then; p=0.0_dp; return; end if
      if (x<=0.0_dp) then; p=0.0_dp; return; end if
      if (x<a+1.0_dp) then
         ap=a; sumv=1.0_dp/a; del=sumv
         do n=1,itmax
            ap=ap+1.0_dp; del=del*x/ap; sumv=sumv+del
            if (abs(del)<abs(sumv)*eps) exit
         end do
         p=sumv*exp(-x+a*log(x)-log_gamma(a))
      else
         b=x+1.0_dp-a; c=1.0_dp/fpmin; d=1.0_dp/b; h=d
         do n=1,itmax
            an=-real(n,dp)*(real(n,dp)-a)
            b=b+2.0_dp
            d=an*d+b; if(abs(d)<fpmin)d=fpmin
            c=b+an/c; if(abs(c)<fpmin)c=fpmin
            d=1.0_dp/d; del=d*c; h=h*del
            if(abs(del-1.0_dp)<eps)exit
         end do
         p=1.0_dp-exp(-x+a*log(x)-log_gamma(a))*h
      end if
      p=clamp01(p)
   end function gamma_p

   pure real(dp) function gamma_q(a,x) result(q)
      real(dp), intent(in) :: a,x
      q=1.0_dp-gamma_p(a,x)
   end function gamma_q

   pure real(dp) function chisq_pdf(x,df) result(y)
      real(dp), intent(in) :: x,df
      real(dp) :: a
      if(x<=0.0_dp .or. df<=0.0_dp)then; y=0.0_dp; return; end if
      a=0.5_dp*df
      y=exp((a-1.0_dp)*log(x)-0.5_dp*x-a*log(2.0_dp)-log_gamma(a))
   end function chisq_pdf

   pure real(dp) function chisq_cdf(x,df) result(p)
      real(dp), intent(in) :: x,df
      if(x<=0.0_dp)then;p=0.0_dp;else;p=gamma_p(0.5_dp*df,0.5_dp*x);end if
   end function chisq_cdf

   real(dp) function randn() result(z)
      real(dp) :: u1,u2
      call random_number(u1); call random_number(u2)
      u1=max(u1,tiny(1.0_dp))
      z=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
   end function randn

   recursive real(dp) function rand_gamma(shape) result(g)
      real(dp), intent(in) :: shape
      real(dp) :: d,c,x,v,u
      if(shape<=0.0_dp)then; g=0.0_dp; return; end if
      if(shape<1.0_dp)then
         call random_number(u)
         g=rand_gamma(shape+1.0_dp)*u**(1.0_dp/shape)
         return
      end if
      d=shape-1.0_dp/3.0_dp; c=1.0_dp/sqrt(9.0_dp*d)
      do
         do
            x=randn(); v=1.0_dp+c*x
            if(v>0.0_dp)exit
         end do
         v=v**3; call random_number(u)
         if(u<1.0_dp-0.0331_dp*x**4)exit
         if(log(u)<0.5_dp*x*x+d*(1.0_dp-v+log(v)))exit
      end do
      g=d*v
   end function rand_gamma

   real(dp) function rand_chisq(df) result(x)
      real(dp), intent(in) :: df
      x=2.0_dp*rand_gamma(0.5_dp*df)
   end function rand_chisq

   function adaptive_integral(f,a,b,tol) result(ans)
      procedure(scalar_fun) :: f
      real(dp), intent(in) :: a,b,tol
      real(dp) :: ans,fa,fb,fc,s
      real(dp) :: c
      c=0.5_dp*(a+b); fa=f(a); fb=f(b); fc=f(c)
      s=(b-a)*(fa+4.0_dp*fc+fb)/6.0_dp
      ans=asr(f,a,b,tol,s,fa,fb,fc,20)
   contains
      recursive function asr(fun,aa,bb,eps,whole,faa,fbb,fcc,depth) result(v)
         procedure(scalar_fun) :: fun
         real(dp), intent(in) :: aa,bb,eps,whole,faa,fbb,fcc
         integer, intent(in) :: depth
         real(dp) :: v,cc,ld,rd,fl,fr,left,right
         cc=0.5_dp*(aa+bb); ld=0.5_dp*(aa+cc); rd=0.5_dp*(cc+bb)
         fl=fun(ld); fr=fun(rd)
         left=(cc-aa)*(faa+4.0_dp*fl+fcc)/6.0_dp
         right=(bb-cc)*(fcc+4.0_dp*fr+fbb)/6.0_dp
         if(depth<=0 .or. abs(left+right-whole)<=15.0_dp*eps)then
            v=left+right+(left+right-whole)/15.0_dp
         else
            v=asr(fun,aa,cc,0.5_dp*eps,left,faa,fcc,fl,depth-1)+ &
              asr(fun,cc,bb,0.5_dp*eps,right,fcc,fbb,fr,depth-1)
         end if
      end function asr
   end function adaptive_integral
end module suppdists_special
