module lmomco_distributions
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use lmomco_kinds, only : dp, pi
   use lmomco_types, only : lmomco_params
   use lmomco_math, only : normal_cdf, normal_pdf, normal_quantile, gamma_p, gamma_q, gamma_quantile, &
      student_t_cdf, student_t_quantile, bessel_i_nu, clamp01
   implicit none
   private
   public :: lmomco_pdf, lmomco_cdf, lmomco_quantile, lmomco_rng
   public :: params_valid

contains

   pure function lower_string(s) result(t)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: t
      integer :: i, k
      t=s
      do i=1,len(s)
         k=iachar(t(i:i))
         if(k>=iachar('A') .and. k<=iachar('Z')) t(i:i)=achar(k+32)
      end do
   end function lower_string

   pure logical function params_valid(par) result(ok)
      type(lmomco_params), intent(in) :: par
      character(len=:), allocatable :: f
      f=trim(lower_string(par%family)); ok=.true.
      select case(f)
      case('cau','exp','gum','lap','nor','ray','revgum','sla','lmrq')
         ok=par%npar>=2 .and. par%p(2)>0.0_dp
      case('gam','kur','emu','kmu','rice')
         ok=par%npar>=2 .and. all(par%p(1:2)>0.0_dp)
         if(f=='rice') ok=par%p(1)>=0.0_dp .and. par%p(2)>0.0_dp
         if(f=='kmu') ok=par%p(1)>=0.0_dp .and. par%p(2)>0.0_dp
         if(f=='emu') ok=par%p(1)>0.0_dp .and. par%p(1)<1.0_dp .and. par%p(2)>0.0_dp
      case('gev','glo','gno','gpa','ln3','pdq3','pdq4','pe3','st3','wei','gep')
         ok=par%npar>=3 .and. par%p(2)>0.0_dp
      case('gov')
         ok=par%npar>=3 .and. par%p(2)>0.0_dp .and. par%p(3)>-1.0_dp
      case('tri')
         ok=par%npar>=3 .and. par%p(1)<=par%p(2) .and. par%p(2)<=par%p(3) .and. par%p(3)>par%p(1)
      case('smd')
         ok=par%npar>=4 .and. par%p(2)>0.0_dp .and. par%p(3)>0.0_dp .and. par%p(4)>0.0_dp
      case('aep4')
         ok=par%npar>=4 .and. par%p(2)>0.0_dp .and. par%p(3)>0.0_dp .and. par%p(4)>0.0_dp
      case('gld')
         ok=par%npar>=4 .and. abs(par%p(2))>tiny(1.0_dp)
      case('kap')
         ok=par%npar>=4 .and. par%p(2)>0.0_dp
      case('wak')
         ok=par%npar>=5 .and. par%p(2)>=0.0_dp .and. par%p(4)>=0.0_dp
      case('texp')
         ok=par%npar>=3 .and. ((par%p(3)>0.0_dp) .or. par%p(2)>0.0_dp)
      case default
         ok=.false.
      end select
   end function params_valid

   real(dp) function lmomco_quantile(prob, par) result(x)
      real(dp), intent(in) :: prob
      type(lmomco_params), intent(in) :: par
      character(len=:), allocatable :: f
      real(dp) :: p,u,a,k,h,g,y,z,b,c,d,alpha,beta,lf,kf2,bu,bb
      if(.not.params_valid(par)) then
         x=ieee_value(0.0_dp,ieee_quiet_nan); return
      end if
      p=max(0.0_dp,min(1.0_dp,prob)); f=trim(lower_string(par%family))
      u=par%p(1); a=par%p(2)
      select case(f)
      case('cau')
         if(p<=0.0_dp) then; x=-huge(1.0_dp)
         else if(p>=1.0_dp) then; x=huge(1.0_dp)
         else; x=u+a*tan(pi*(p-0.5_dp)); end if
      case('exp')
         if(p>=1.0_dp) then; x=huge(1.0_dp); else; x=u-a*log(1.0_dp-p); end if
      case('gam')
         x=a*gamma_quantile(p,u)
      case('gev')
         k=par%p(3)
         if(p<=0.0_dp) then
            if(k<0.0_dp) then; x=u+a/k; else; x=-huge(1.0_dp); end if
         else if(p>=1.0_dp) then
            if(k>0.0_dp) then; x=u+a/k; else; x=huge(1.0_dp); end if
         else
            y=-log(-log(p))
            if(abs(k)>sqrt(epsilon(1.0_dp))) y=(1.0_dp-exp(-k*y))/k
            x=u+a*y
         end if
      case('glo')
         k=par%p(3)
         if(p<=0.0_dp) then
            if(k<0.0_dp) then; x=u+a/k; else; x=-huge(1.0_dp); end if
         else if(p>=1.0_dp) then
            if(k>0.0_dp) then; x=u+a/k; else; x=huge(1.0_dp); end if
         else
            y=log(p/(1.0_dp-p)); if(abs(k)>sqrt(epsilon(1.0_dp))) y=(1.0_dp-exp(-k*y))/k
            x=u+a*y
         end if
      case('gno')
         k=par%p(3); y=normal_quantile(p)
         if(abs(k)>sqrt(epsilon(1.0_dp)) .and. abs(y)<huge(1.0_dp)/100.0_dp) y=(1.0_dp-exp(-k*y))/k
         x=u+a*y
      case('gov')
         b=par%p(3); x=u+a*((b+1.0_dp)*p**b-b*p**(b+1.0_dp))
      case('gpa')
         k=par%p(3)
         if(p<=0.0_dp) then; x=u
         else if(p>=1.0_dp) then
            if(k>0.0_dp) then; x=u+a/k; else; x=huge(1.0_dp); end if
         else
            y=-log(1.0_dp-p); if(abs(k)>sqrt(epsilon(1.0_dp))) y=(1.0_dp-exp(-k*y))/k
            x=u+a*y
         end if
      case('gum')
         if(p<=0.0_dp) then; x=-huge(1.0_dp)
         else if(p>=1.0_dp) then; x=huge(1.0_dp)
         else; x=u-a*log(-log(p)); end if
      case('revgum')
         if(p<=0.0_dp) then; x=-huge(1.0_dp)
         else if(p>=1.0_dp) then; x=huge(1.0_dp)
         else; x=u+a*log(-log(1.0_dp-p)); end if
      case('lap')
         if(p<=0.0_dp) then; x=-huge(1.0_dp)
         else if(p>=1.0_dp) then; x=huge(1.0_dp)
         else if(p<=0.5_dp) then; x=u+a*log(2.0_dp*p)
         else; x=u-a*log(2.0_dp*(1.0_dp-p)); end if
      case('ln3')
         x=exp(normal_quantile(p)*par%p(3)+a)+u
      case('nor')
         x=u+a*normal_quantile(p)
      case('ray')
         if(p>=1.0_dp) then; x=huge(1.0_dp); else; x=u+sqrt(-2.0_dp*a*a*log(1.0_dp-p)); end if
      case('tri')
         b=par%p(2); c=par%p(3); a=c-u; d=b-u; z=(c-b); y=d/a
         if(p<y) then; x=u+sqrt(p*a*d)
         else if(p>y) then; x=c-sqrt((1.0_dp-p)*a*z)
         else; x=b; end if
      case('kur')
         x=(1.0_dp-(1.0_dp-p)**(1.0_dp/a))**(1.0_dp/u)
      case('smd')
         if(p>=1.0_dp) then; x=huge(1.0_dp)
         else; x=u+a*((1.0_dp-p)**(-1.0_dp/par%p(4))-1.0_dp)**(1.0_dp/par%p(3)); end if
      case('st3')
         x=u+a*student_t_quantile(p,max(1.001_dp,par%p(3)))
      case('pe3')
         g=par%p(3)
         if(abs(g)<=sqrt(epsilon(1.0_dp))) then; x=u+a*normal_quantile(p)
         else
            alpha=4.0_dp/(g*g); beta=abs(0.5_dp*a*g)
            if(g>0.0_dp) then; x=u-alpha*beta+beta*gamma_quantile(p,alpha)
            else; x=u+alpha*beta-beta*gamma_quantile(1.0_dp-p,alpha); end if
         end if
      case('gld')
         x=u+a*(p**par%p(3)-(1.0_dp-p)**par%p(4))
      case('lmrq')
         if(p>=1.0_dp) then; x=huge(1.0_dp)
         else; x=-(u+a)*log(1.0_dp-p)-2.0_dp*a*p; end if
      case('pdq3')
         k=max(-(1.0_dp-sqrt(epsilon(1.0_dp))),min(1.0_dp-sqrt(epsilon(1.0_dp)),par%p(3)))
         z=max(epsilon(1.0_dp),min(1.0_dp-epsilon(1.0_dp),p)); lf=log(z/(1.0_dp-z))
         x=u+a*(lf+k*log((1.0_dp-k*(2.0_dp*z-1.0_dp))**2/(4.0_dp*z*(1.0_dp-z))))
      case('pdq4')
         k=min(1.0_dp-sqrt(epsilon(1.0_dp)),par%p(3)); z=max(epsilon(1.0_dp),min(1.0_dp-epsilon(1.0_dp),p))
         lf=log(z/(1.0_dp-z)); kf2=k*(2.0_dp*z-1.0_dp)
         if(abs(k)<sqrt(epsilon(1.0_dp))) then; x=u+a*lf
         else if(k>0.0_dp) then; x=u+a*(lf-2.0_dp*k*atanh(kf2))
         else; x=u+a*(lf+2.0_dp*k*atan(kf2)); end if
      case('wak')
         b=par%p(3); c=par%p(4); d=par%p(5)
         if(p<=0.0_dp) then; x=u
         else if(p>=1.0_dp) then; x=huge(1.0_dp)
         else
            z=-log(1.0_dp-p)
            if(abs(b)<sqrt(epsilon(1.0_dp))) then; y=z; else; y=(1.0_dp-exp(-b*z))/b; end if
            if(abs(d)<sqrt(epsilon(1.0_dp))) then; bb=z; else; bb=(1.0_dp-exp(d*z))/(-d); end if
            x=u+a*y+c*bb
         end if
      case('wei')
         b=par%p(2); d=par%p(3)
         if(p>=1.0_dp) then; x=huge(1.0_dp)
         else; x=u+b*(-log(1.0_dp-p))**(1.0_dp/d); end if
      case('gep')
         k=par%p(2); h=par%p(3)
         if(p>=1.0_dp) then; x=huge(1.0_dp)
         else; x=-u*log(1.0_dp+(1.0_dp/h)*log(1.0_dp-p**(1.0_dp/k)*(1.0_dp-exp(-h)))); end if
      case('kap')
         g=par%p(3); h=par%p(4)
         if(p<=0.0_dp .or. p>=1.0_dp) then
            x=quantile_by_cdf(p,par)
         else
            y=-log(p); if(abs(h)>sqrt(epsilon(1.0_dp))) y=(1.0_dp-exp(-h*y))/h
            y=-log(y); if(abs(g)>sqrt(epsilon(1.0_dp))) y=(1.0_dp-exp(-g*y))/g
            x=u+a*y
         end if
      case('aep4')
         k=par%p(3); h=par%p(4); y=lmomco_cdf(u,par)
         if(p<y) then
            z=gamma_quantile(1.0_dp-p*(1.0_dp+k*k)/(k*k),1.0_dp/h)
            x=u-a*k*z**(1.0_dp/h)
         else
            z=gamma_quantile(1.0_dp-(1.0_dp-p)*(1.0_dp+k*k),1.0_dp/h)
            x=u+a/k*z**(1.0_dp/h)
         end if
      case('texp')
         if(par%p(3)>0.0_dp) then; x=par%p(3)*p
         else if(u<=0.0_dp) then; x=-a*log(1.0_dp-p)
         else; bu=1.0_dp-exp(-u/a); x=-a*log(1.0_dp-p*bu); end if
      case('rice','emu','kmu','sla')
         x=quantile_by_cdf(p,par)
      case default
         x=ieee_value(0.0_dp,ieee_quiet_nan)
      end select
   end function lmomco_quantile

   real(dp) function lmomco_cdf(x, par) result(p)
      real(dp), intent(in) :: x
      type(lmomco_params), intent(in) :: par
      character(len=:), allocatable :: f
      real(dp) :: u,a,k,z,y,h,iga,igb
      if(.not.params_valid(par)) then; p=ieee_value(0.0_dp,ieee_quiet_nan); return; end if
      f=trim(lower_string(par%family)); u=par%p(1); a=par%p(2)
      select case(f)
      case('cau')
         p=0.5_dp+atan((x-u)/a)/pi
      case('exp')
         if(x<u) then; p=0.0_dp; else; p=1.0_dp-exp(-(x-u)/a); end if
      case('gam')
         if(x<=0.0_dp) then; p=0.0_dp; else; p=gamma_p(u,x/a); end if
      case('nor')
         p=normal_cdf((x-u)/a)
      case('ln3')
         if(x<=u) then; p=0.0_dp; else; p=normal_cdf((log(x-u)-a)/par%p(3)); end if
      case('lap')
         if(x<=u) then; p=0.5_dp*exp((x-u)/a); else; p=1.0_dp-0.5_dp*exp(-(x-u)/a); end if
      case('ray')
         if(x<=u) then; p=0.0_dp; else; p=1.0_dp-exp(-0.5_dp*((x-u)/a)**2); end if
      case('gum')
         p=exp(-exp(-(x-u)/a))
      case('revgum')
         p=1.0_dp-exp(-exp((x-u)/a))
      case('gpa')
         k=par%p(3); z=(x-u)/a
         if(z<=0.0_dp) then; p=0.0_dp
         else if(abs(k)<sqrt(epsilon(1.0_dp))) then; p=1.0_dp-exp(-z)
         else if(1.0_dp-k*z<=0.0_dp) then; p=1.0_dp
         else; p=1.0_dp-(1.0_dp-k*z)**(1.0_dp/k); end if
      case('gev')
         k=par%p(3); z=(x-u)/a
         if(abs(k)<sqrt(epsilon(1.0_dp))) then; p=exp(-exp(-z))
         else
            y=1.0_dp-k*z
            if(y<=0.0_dp) then
               if(k>0.0_dp) then; p=1.0_dp; else; p=0.0_dp; end if
            else; p=exp(-y**(1.0_dp/k)); end if
         end if
      case('glo')
         k=par%p(3); z=(x-u)/a
         if(abs(k)<sqrt(epsilon(1.0_dp))) then; p=1.0_dp/(1.0_dp+exp(-z))
         else
            y=1.0_dp-k*z
            if(y<=0.0_dp) then
               if(k>0.0_dp) then; p=1.0_dp; else; p=0.0_dp; end if
            else; p=1.0_dp/(1.0_dp+y**(1.0_dp/k)); end if
         end if
      case('gno')
         k=par%p(3); z=(x-u)/a
         if(abs(k)<sqrt(epsilon(1.0_dp))) then; p=normal_cdf(z)
         else
            y=1.0_dp-k*z
            if(y<=0.0_dp) then
               if(k>0.0_dp) then; p=1.0_dp; else; p=0.0_dp; end if
            else; p=normal_cdf(-log(y)/k); end if
         end if
      case('tri')
         if(x<=par%p(1)) then; p=0.0_dp
         else if(x>=par%p(3)) then; p=1.0_dp
         else if(x<=par%p(2)) then; p=(x-par%p(1))**2/((par%p(3)-par%p(1))*(par%p(2)-par%p(1)))
         else; p=1.0_dp-(par%p(3)-x)**2/((par%p(3)-par%p(1))*(par%p(3)-par%p(2))); end if
      case('kur')
         if(x<=0.0_dp) then; p=0.0_dp
         else if(x>=1.0_dp) then; p=1.0_dp
         else; p=1.0_dp-(1.0_dp-x**u)**a; end if
      case('smd')
         if(x<=u) then; p=0.0_dp
         else; p=1.0_dp-(1.0_dp+((x-u)/a)**par%p(3))**(-par%p(4)); end if
      case('st3')
         p=student_t_cdf((x-u)/a,max(1.001_dp,par%p(3)))
      case('sla')
         z=(x-u)/a
         if(abs(z)<1.0e-10_dp) then; p=0.5_dp
         else; p=normal_cdf(z)-(normal_pdf(0.0_dp)-normal_pdf(z))/z; end if
      case('aep4')
         k=par%p(3); h=par%p(4); iga=1.0_dp/h
         if(x<u) then; igb=((u-x)/(a*k))**h; p=k*k/(1.0_dp+k*k)*gamma_q(iga,igb)
         else; igb=(k*(x-u)/a)**h; p=1.0_dp-gamma_q(iga,igb)/(1.0_dp+k*k); end if
      case('rice','emu','kmu')
         p=cdf_special_positive(x,par)
      case default
         p=cdf_from_quantile(x,par)
      end select
      p=clamp01(p)
   end function lmomco_cdf

   real(dp) function lmomco_pdf(x, par) result(v)
      real(dp), intent(in) :: x
      type(lmomco_params), intent(in) :: par
      character(len=:), allocatable :: f
      real(dp) :: u,a,z,k,h,del
      if(.not.params_valid(par)) then; v=ieee_value(0.0_dp,ieee_quiet_nan); return; end if
      f=trim(lower_string(par%family)); u=par%p(1); a=par%p(2)
      select case(f)
      case('cau'); z=(x-u)/a; v=1.0_dp/(pi*a*(1.0_dp+z*z))
      case('exp'); if(x<u) then; v=0.0_dp; else; v=exp(-(x-u)/a)/a; end if
      case('gam'); if(x<=0.0_dp) then; v=0.0_dp; else; v=exp((u-1.0_dp)*log(x/a)-x/a-log_gamma(u))/a; end if
      case('nor'); v=normal_pdf((x-u)/a)/a
      case('ln3'); if(x<=u) then; v=0.0_dp; else; z=(log(x-u)-a)/par%p(3); v=normal_pdf(z)/(par%p(3)*(x-u)); end if
      case('lap'); v=0.5_dp*exp(-abs(x-u)/a)/a
      case('ray'); if(x<=u) then; v=0.0_dp; else; z=(x-u)/a; v=z*exp(-0.5_dp*z*z)/a; end if
      case('gum'); z=(x-u)/a; v=exp(-z-exp(-z))/a
      case('revgum'); z=(x-u)/a; v=exp(z-exp(z))/a
      case('tri')
         if(x<par%p(1).or.x>par%p(3)) then; v=0.0_dp
         else if(x<=par%p(2)) then; v=2.0_dp*(x-par%p(1))/((par%p(3)-par%p(1))*(par%p(2)-par%p(1)))
         else; v=2.0_dp*(par%p(3)-x)/((par%p(3)-par%p(1))*(par%p(3)-par%p(2))); end if
      case('kur'); if(x<=0.0_dp.or.x>=1.0_dp) then; v=0.0_dp; else; v=u*a*x**(u-1.0_dp)*(1.0_dp-x**u)**(a-1.0_dp); end if
      case('smd')
         if(x<=u) then
            v=0.0_dp
         else
            z=(x-u)/a
            v=par%p(3)*par%p(4)*z**(par%p(3)-1.0_dp) / &
               (a*(1.0_dp+z**par%p(3))**(par%p(4)+1.0_dp))
         end if
      case('sla')
         z=(x-u)/a; del=normal_pdf(0.0_dp)-normal_pdf(z)
         if(abs(z)<1.0e-8_dp) then; v=1.0_dp/(2.0_dp*sqrt(2.0_dp*pi)*a); else; v=del/(z*z*a); end if
      case('aep4')
         k=par%p(3); h=par%p(4); z=h*k/(a*(1.0_dp+k*k)*gamma(1.0_dp/h)); yabs: block
            real(dp) :: yy, sk
            yy=abs(x-u)/a
            if(x-u>=0.0_dp) then; sk=k; else; sk=1.0_dp/k; end if
            v=z*exp(-(sk*yy)**h)
         end block yabs
      case('rice','emu','kmu')
         v=pdf_special_positive(x,par)
      case default
         v=pdf_from_quantile(x,par)
      end select
      if(v<0.0_dp .and. abs(v)<1.0e-12_dp) v=0.0_dp
   end function lmomco_pdf

   subroutine lmomco_rng(par, x)
      type(lmomco_params), intent(in) :: par
      real(dp), intent(out) :: x(:)
      real(dp) :: u
      integer :: i
      do i=1,size(x)
         call random_number(u)
         u=max(epsilon(1.0_dp),min(1.0_dp-epsilon(1.0_dp),u))
         x(i)=lmomco_quantile(u,par)
      end do
   end subroutine lmomco_rng

   real(dp) function cdf_from_quantile(x,par) result(p)
      real(dp), intent(in) :: x
      type(lmomco_params), intent(in) :: par
      real(dp) :: lo,hi,mid,q
      integer :: i
      lo=epsilon(1.0_dp); hi=1.0_dp-epsilon(1.0_dp)
      if(x<=lmomco_quantile(lo,par)) then; p=0.0_dp; return; end if
      if(x>=lmomco_quantile(hi,par)) then; p=1.0_dp; return; end if
      do i=1,100
         mid=0.5_dp*(lo+hi); q=lmomco_quantile(mid,par)
         if(q<x) then; lo=mid; else; hi=mid; end if
      end do
      p=0.5_dp*(lo+hi)
   end function cdf_from_quantile

   real(dp) function pdf_from_quantile(x,par) result(v)
      real(dp), intent(in) :: x
      type(lmomco_params), intent(in) :: par
      real(dp) :: p,h,q1,q2
      p=cdf_from_quantile(x,par); h=max(1.0e-7_dp,min(1.0e-4_dp,0.05_dp*min(p,1.0_dp-p)))
      if(h<=0.0_dp) then; v=0.0_dp; return; end if
      q1=lmomco_quantile(max(epsilon(1.0_dp),p-h),par)
      q2=lmomco_quantile(min(1.0_dp-epsilon(1.0_dp),p+h),par)
      if(q2<=q1) then; v=0.0_dp; else; v=(2.0_dp*h)/(q2-q1); end if
   end function pdf_from_quantile

   real(dp) function quantile_by_cdf(p,par) result(x)
      real(dp), intent(in) :: p
      type(lmomco_params), intent(in) :: par
      real(dp) :: lo,hi,mid,scale
      integer :: i
      if(p<=0.0_dp) then; x=-huge(1.0_dp); return; end if
      if(p>=1.0_dp) then; x=huge(1.0_dp); return; end if
      scale=max(1.0_dp,maxval(abs(par%p(1:max(1,par%npar)))))
      lo=-scale; hi=scale
      if(trim(lower_string(par%family))=='rice' .or. trim(lower_string(par%family))=='emu' .or. &
         trim(lower_string(par%family))=='kmu') lo=0.0_dp
      do while(lmomco_cdf(lo,par)>p .and. abs(lo)<1.0e12_dp*scale); lo=2.0_dp*lo-scale; end do
      do while(lmomco_cdf(hi,par)<p .and. hi<1.0e12_dp*scale); hi=2.0_dp*hi+scale; end do
      do i=1,110
         mid=0.5_dp*(lo+hi)
         if(lmomco_cdf(mid,par)<p) then; lo=mid; else; hi=mid; end if
      end do
      x=0.5_dp*(lo+hi)
   end function quantile_by_cdf

   real(dp) function pdf_special_positive(x,par) result(v)
      real(dp), intent(in) :: x
      type(lmomco_params), intent(in) :: par
      character(len=:), allocatable :: f
      real(dp) :: a,k,mu,eta,m,h,hh,nu,toi,tmpb,bv
      f=trim(lower_string(par%family)); v=0.0_dp
      if(x<0.0_dp) return
      select case(f)
      case('rice')
         k=par%p(1); a=par%p(2)
         if(abs(k)<=tiny(1.0_dp)) then; v=x/(a*a)*exp(-x*x/(2.0_dp*a*a))
         else
            toi=x*k/(a*a)
            if(toi<80.0_dp) then
               v=x/(a*a)*exp(-(x*x+k*k)/(2.0_dp*a*a))*bessel_i_nu(toi,0.0_dp)
            else
               v=x/(a*a)*exp(-(x-k)**2/(2.0_dp*a*a))/sqrt(2.0_dp*pi*toi)
            end if
         end if
      case('emu')
         eta=par%p(1); m=par%p(2); h=1.0_dp/(1.0_dp-eta*eta); hh=eta/(1.0_dp-eta*eta); nu=m-0.5_dp
         if(eta<=1.0e-6_dp) then
            v=2.0_dp*(2.0_dp*m)**(2.0_dp*m)/gamma(2.0_dp*m)*x**(4.0_dp*m-1.0_dp)*exp(-2.0_dp*m*x*x)
         else
            toi=2.0_dp*m*hh*x*x
            tmpb=4.0_dp*sqrt(pi)*m**(m+0.5_dp)*h**m/(gamma(m)*hh**nu)
            v=tmpb*x**(2.0_dp*m)*exp(-2.0_dp*m*h*x*x)*bessel_i_nu(toi,nu)
         end if
      case('kmu')
         k=par%p(1); mu=par%p(2)
         if(k<=1.0e-12_dp) then
            v=2.0_dp*mu**mu/gamma(mu)*x**(2.0_dp*mu-1.0_dp)*exp(-mu*x*x)
         else
            tmpb=2.0_dp*mu*(1.0_dp+k)**((mu+1.0_dp)/2.0_dp)/(k**((mu-1.0_dp)/2.0_dp)*exp(mu*k))
            toi=2.0_dp*mu*sqrt(k*(k+1.0_dp))*x
            bv=bessel_i_nu(toi,mu-1.0_dp)
            v=tmpb*x**mu*exp(-mu*(1.0_dp+k)*x*x)*bv
         end if
      end select
      if(.not.(v>=0.0_dp)) v=0.0_dp
   end function pdf_special_positive

   real(dp) function cdf_special_positive(x,par) result(p)
      real(dp), intent(in) :: x
      type(lmomco_params), intent(in) :: par
      integer, parameter :: n=1200
      integer :: i
      real(dp) :: h,s,xx
      if(x<=0.0_dp) then; p=0.0_dp; return; end if
      h=x/real(n,dp); s=pdf_special_positive(0.0_dp,par)+pdf_special_positive(x,par)
      do i=1,n-1
         xx=h*real(i,dp)
         if(mod(i,2)==0) then; s=s+2.0_dp*pdf_special_positive(xx,par)
         else; s=s+4.0_dp*pdf_special_positive(xx,par); end if
      end do
      p=clamp01(s*h/3.0_dp)
   end function cdf_special_positive

end module lmomco_distributions
