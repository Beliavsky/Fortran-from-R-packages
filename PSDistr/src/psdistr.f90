module psdistr
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use psdistr_kinds, only : dp
   use psdistr_math, only : normal_pdf, normal_cdf, normal_quantile, signed_root, regularized_beta, beta_quantile
   implicit none
   private
   public :: dp
   public :: dtppn, ptppn, qtppn, rtppn
   public :: dpc, ppc, qpc, rpc
   public :: ddsn, pdsn, qdsn, rdsn
   public :: den, pen, qen, ren
   public :: dspc, pspc, qspc, rspc
   public :: deck, peck, qeck, reck

   interface dtppn
      module procedure dtppn_scalar, dtppn_vec
   end interface
   interface ptppn
      module procedure ptppn_scalar, ptppn_vec
   end interface
   interface qtppn
      module procedure qtppn_scalar, qtppn_vec
   end interface
   interface dpc
      module procedure dpc_scalar, dpc_vec
   end interface
   interface ppc
      module procedure ppc_scalar, ppc_vec
   end interface
   interface qpc
      module procedure qpc_scalar, qpc_vec
   end interface
   interface ddsn
      module procedure ddsn_scalar, ddsn_vec
   end interface
   interface pdsn
      module procedure pdsn_scalar, pdsn_vec
   end interface
   interface qdsn
      module procedure qdsn_scalar, qdsn_vec
   end interface
   interface den
      module procedure den_scalar, den_vec
   end interface
   interface pen
      module procedure pen_scalar, pen_vec
   end interface
   interface qen
      module procedure qen_scalar, qen_vec
   end interface
   interface dspc
      module procedure dspc_scalar, dspc_vec
   end interface
   interface pspc
      module procedure pspc_scalar, pspc_vec
   end interface
   interface qspc
      module procedure qspc_scalar, qspc_vec
   end interface
   interface deck
      module procedure deck_scalar, deck_vec
   end interface
   interface peck
      module procedure peck_scalar, peck_vec
   end interface
   interface qeck
      module procedure qeck_scalar, qeck_vec
   end interface

contains
   elemental real(dp) function nanv() result(y)
      y=ieee_value(0.0_dp,ieee_quiet_nan)
   end function nanv

   elemental real(dp) function dtppn_scalar(x,teta,s1,s2,c) result(y)
      real(dp),intent(in)::x,teta,s1,s2,c
      real(dp)::z
      if(s1<=0.or.s2<=0.or.c<1) then; y=nanv(); return; end if
      if(x<teta) then
         z=(teta-x)/s1
         if(abs(z)<=tiny(1.0_dp) .and. abs(c-1.0_dp)<=epsilon(1.0_dp)) then; y=normal_pdf(0.0_dp)/s1
         else; y=c*z**(c-1.0_dp)*normal_pdf(z**c)/s1; end if
      else
         z=(x-teta)/s2
         if(abs(z)<=tiny(1.0_dp) .and. abs(c-1.0_dp)<=epsilon(1.0_dp)) then; y=normal_pdf(0.0_dp)/s2
         else; y=c*z**(c-1.0_dp)*normal_pdf(z**c)/s2; end if
      end if
   end function
   pure function dtppn_vec(x,teta,s1,s2,c) result(y)
      real(dp),intent(in)::x(:),teta,s1,s2,c; real(dp)::y(size(x)); integer::i
      do i=1,size(x); y(i)=dtppn_scalar(x(i),teta,s1,s2,c); end do
   end function

   elemental real(dp) function ptppn_scalar(x,teta,s1,s2,c) result(y)
      real(dp),intent(in)::x,teta,s1,s2,c
      if(s1<=0.or.s2<=0.or.c<1) then; y=nanv(); return; end if
      if(x<teta) then; y=normal_cdf(-((teta-x)/s1)**c)
      else; y=normal_cdf(((x-teta)/s2)**c); end if
   end function
   pure function ptppn_vec(x,teta,s1,s2,c) result(y)
      real(dp),intent(in)::x(:),teta,s1,s2,c; real(dp)::y(size(x)); integer::i
      do i=1,size(x); y(i)=ptppn_scalar(x(i),teta,s1,s2,c); end do
   end function

   elemental real(dp) function qtppn_scalar(p,teta,s1,s2,c) result(x)
      real(dp),intent(in)::p,teta,s1,s2,c; real(dp)::q
      if(p<=0.or.p>=1.or.s1<=0.or.s2<=0.or.c<1) then; x=nanv(); return; end if
      q=normal_quantile(p)
      if(p<0.5_dp) then; x=teta-s1*signed_root(-q,c)
      else; x=teta+s2*signed_root(q,c); end if
   end function
   pure function qtppn_vec(p,teta,s1,s2,c) result(x)
      real(dp),intent(in)::p(:),teta,s1,s2,c; real(dp)::x(size(p)); integer::i
      do i=1,size(p); x(i)=qtppn_scalar(p(i),teta,s1,s2,c); end do
   end function
   subroutine rtppn(n,teta,s1,s2,c,x)
      integer,intent(in)::n; real(dp),intent(in)::teta,s1,s2,c; real(dp),intent(out)::x(n)
      real(dp)::u,q; integer::i
      if(s1<=0.or.s2<=0.or.c<1) then; x=nanv(); return; end if
      do i=1,n; call random_number(u); u=max(tiny(1.0_dp),min(1.0_dp-epsilon(1.0_dp),u)); q=normal_quantile(u)
         if(q<0) then; x(i)=teta-s1*(-q)**(1.0_dp/c); else; x(i)=teta+s2*q**(1.0_dp/c); end if
      end do
   end subroutine

   elemental real(dp) function dpc_scalar(x,teta,s2,c) result(y)
      real(dp),intent(in)::x,teta,s2,c; real(dp)::z,a
      if(s2<=0.or.c<1) then; y=nanv(); return; end if
      z=(x-teta)/s2; a=abs(z)
      if(abs(a)<=tiny(1.0_dp)) then
         if(abs(c-1.0_dp)<=epsilon(1.0_dp)) then; y=normal_pdf(0.0_dp)/s2; else; y=0.0_dp; end if
      else; y=c/s2*a**(c-1.0_dp)*normal_pdf(a**c); end if
   end function
   pure function dpc_vec(x,teta,s2,c) result(y)
      real(dp),intent(in)::x(:),teta,s2,c; real(dp)::y(size(x)); integer::i
      do i=1,size(x); y(i)=dpc_scalar(x(i),teta,s2,c); end do
   end function
   elemental real(dp) function ppc_scalar(x,teta,s2,c) result(y)
      real(dp),intent(in)::x,teta,s2,c; real(dp)::e
      if(s2<=0.or.c<1) then; y=nanv(); return; end if
      e=(x-teta)/s2; y=normal_cdf(sign(abs(e)**c,e))
   end function
   pure function ppc_vec(x,teta,s2,c) result(y)
      real(dp),intent(in)::x(:),teta,s2,c; real(dp)::y(size(x));integer::i
      do i=1,size(x);y(i)=ppc_scalar(x(i),teta,s2,c);end do
   end function
   elemental real(dp) function qpc_scalar(p,teta,s2,c) result(x)
      real(dp),intent(in)::p,teta,s2,c; real(dp)::q
      if(p<=0.or.p>=1.or.s2<=0.or.c<1) then;x=nanv();return;end if
      q=normal_quantile(p); x=teta+s2*signed_root(q,c)
   end function
   pure function qpc_vec(p,teta,s2,c) result(x)
      real(dp),intent(in)::p(:),teta,s2,c;real(dp)::x(size(p));integer::i
      do i=1,size(p);x(i)=qpc_scalar(p(i),teta,s2,c);end do
   end function
   subroutine rpc(n,teta,s2,c,x)
      integer,intent(in)::n;real(dp),intent(in)::teta,s2,c;real(dp),intent(out)::x(n);real(dp)::u;integer::i
      if(s2<=0.or.c<1)then;x=nanv();return;end if
      do i=1,n;call random_number(u);u=max(tiny(1.0_dp),min(1.0_dp-epsilon(1.0_dp),u));x(i)=qpc_scalar(u,teta,s2,c);end do
   end subroutine

   elemental real(dp) function ddsn_scalar(x,a,b,c,teta) result(y)
      real(dp),intent(in)::x,a,b,c,teta;real(dp)::z,u
      if(a<0.or.b<0.or.a+b<=0)then;y=nanv();return;end if
      z=x-teta;u=a*z**3+b*z+c;y=(3*a*z*z+b)*normal_pdf(u)
   end function
   pure function ddsn_vec(x,a,b,c,teta) result(y)
      real(dp),intent(in)::x(:),a,b,c,teta;real(dp)::y(size(x));integer::i
      do i=1,size(x);y(i)=ddsn_scalar(x(i),a,b,c,teta);end do
   end function
   elemental real(dp) function pdsn_scalar(x,a,b,c,teta) result(y)
      real(dp),intent(in)::x,a,b,c,teta;real(dp)::z
      if(a<0.or.b<0.or.a+b<=0)then;y=nanv();return;end if
      z=x-teta;y=normal_cdf(a*z**3+b*z+c)
   end function
   pure function pdsn_vec(x,a,b,c,teta) result(y)
      real(dp),intent(in)::x(:),a,b,c,teta;real(dp)::y(size(x));integer::i
      do i=1,size(x);y(i)=pdsn_scalar(x(i),a,b,c,teta);end do
   end function
   elemental real(dp) function solve_cubic_monotone(w,a,b,c,teta) result(x)
      real(dp),intent(in)::w,a,b,c,teta;real(dp)::q,t,u
      if(abs(a)<=tiny(1.0_dp))then;x=(w-c)/b+teta
      else if(abs(b)<=tiny(1.0_dp))then;x=signed_root((w-c)/a,3.0_dp)+teta
      else
         q=(c-w)/a;t=0.5_dp*(-q+sqrt(q*q+4.0_dp*(b/a)**3/27.0_dp));u=signed_root(t,3.0_dp)
         if(abs(u)<tiny(1.0_dp))then;x=teta;else;x=u-b/(3.0_dp*a*u)+teta;end if
      end if
   end function
   elemental real(dp) function qdsn_scalar(p,a,b,c,teta) result(x)
      real(dp),intent(in)::p,a,b,c,teta
      if(p<=0.or.p>=1.or.a<0.or.b<0.or.a+b<=0)then;x=nanv();return;end if
      x=solve_cubic_monotone(normal_quantile(p),a,b,c,teta)
   end function
   pure function qdsn_vec(p,a,b,c,teta) result(x)
      real(dp),intent(in)::p(:),a,b,c,teta;real(dp)::x(size(p));integer::i
      do i=1,size(p);x(i)=qdsn_scalar(p(i),a,b,c,teta);end do
   end function
   subroutine rdsn(n,a,b,c,teta,x)
      integer,intent(in)::n;real(dp),intent(in)::a,b,c,teta;real(dp),intent(out)::x(n);real(dp)::u;integer::i
      if(a<0.or.b<0.or.a+b<=0)then;x=nanv();return;end if
      do i=1,n;call random_number(u);u=max(tiny(1.0_dp),min(1.0_dp-epsilon(1.0_dp),u));x(i)=qdsn_scalar(u,a,b,c,teta);end do
   end subroutine

   elemental real(dp) function den_scalar(x,a1,b1,a2,b2,c) result(y)
      real(dp),intent(in)::x,a1,b1,a2,b2,c;real(dp)::e1,e2,u
      if(b1<=0.or.b2<=0)then;y=nanv();return;end if
      e1=exp((a1-x)/b1);e2=exp((x-a2)/b2);u=c-e1+e2;y=(e1/b1+e2/b2)*normal_pdf(u)
   end function
   pure function den_vec(x,a1,b1,a2,b2,c) result(y)
      real(dp),intent(in)::x(:),a1,b1,a2,b2,c;real(dp)::y(size(x));integer::i
      do i=1,size(x);y(i)=den_scalar(x(i),a1,b1,a2,b2,c);end do
   end function
   elemental real(dp) function pen_scalar(x,a1,b1,a2,b2,c) result(y)
      real(dp),intent(in)::x,a1,b1,a2,b2,c
      if(b1<=0.or.b2<=0)then;y=nanv();return;end if
      y=normal_cdf(c-exp((a1-x)/b1)+exp((x-a2)/b2))
   end function
   pure function pen_vec(x,a1,b1,a2,b2,c) result(y)
      real(dp),intent(in)::x(:),a1,b1,a2,b2,c;real(dp)::y(size(x));integer::i
      do i=1,size(x);y(i)=pen_scalar(x(i),a1,b1,a2,b2,c);end do
   end function
   real(dp) function qen_scalar(p,a1,b1,a2,b2,c) result(x)
      real(dp),intent(in)::p,a1,b1,a2,b2,c;real(dp)::lo,hi,mid;integer::i
      if(p<=0.or.p>=1.or.b1<=0.or.b2<=0)then;x=nanv();return;end if
      lo=min(a1,a2)-max(b1,b2);hi=max(a1,a2)+max(b1,b2)
      do while(pen_scalar(lo,a1,b1,a2,b2,c)>p);lo=lo-2.0_dp*max(b1,b2);end do
      do while(pen_scalar(hi,a1,b1,a2,b2,c)<p);hi=hi+2.0_dp*max(b1,b2);end do
      do i=1,120;mid=0.5_dp*(lo+hi);if(pen_scalar(mid,a1,b1,a2,b2,c)<p)then;lo=mid;else;hi=mid;end if;end do
      x=0.5_dp*(lo+hi)
   end function
   function qen_vec(p,a1,b1,a2,b2,c) result(x)
      real(dp),intent(in)::p(:),a1,b1,a2,b2,c;real(dp)::x(size(p));integer::i
      do i=1,size(p);x(i)=qen_scalar(p(i),a1,b1,a2,b2,c);end do
   end function
   subroutine ren(n,a1,b1,a2,b2,c,x)
      integer,intent(in)::n;real(dp),intent(in)::a1,b1,a2,b2,c;real(dp),intent(out)::x(n);real(dp)::u;integer::i
      if(b1<=0.or.b2<=0)then;x=nanv();return;end if
      do i=1,n;call random_number(u);u=max(tiny(1.0_dp),min(1.0_dp-epsilon(1.0_dp),u));x(i)=qen_scalar(u,a1,b1,a2,b2,c);end do
   end subroutine

   elemental real(dp) function dspc_scalar(x,a,b,c,d,teta) result(y)
      real(dp),intent(in)::x,a,b,c,d,teta;real(dp)::z,u
      if(a<0.or.b<0.or.a+b<=0.or.d<1)then;y=nanv();return;end if
      z=x-teta;u=a*z**3+b*z+c
      if(abs(u)<=tiny(1.0_dp).and.d>1.0_dp)then;y=0.0_dp
      else;y=d*(3*a*z*z+b)*abs(u)**(d-1.0_dp)*normal_pdf(abs(u)**d);end if
   end function
   pure function dspc_vec(x,a,b,c,d,teta) result(y)
      real(dp),intent(in)::x(:),a,b,c,d,teta;real(dp)::y(size(x));integer::i
      do i=1,size(x);y(i)=dspc_scalar(x(i),a,b,c,d,teta);end do
   end function
   elemental real(dp) function pspc_scalar(x,a,b,c,d,teta) result(y)
      real(dp),intent(in)::x,a,b,c,d,teta;real(dp)::z,u,w
      if(a<0.or.b<0.or.a+b<=0.or.d<1)then;y=nanv();return;end if
      z=x-teta;u=a*z**3+b*z+c;w=sign(abs(u)**d,u);y=normal_cdf(w)
   end function
   pure function pspc_vec(x,a,b,c,d,teta) result(y)
      real(dp),intent(in)::x(:),a,b,c,d,teta;real(dp)::y(size(x));integer::i
      do i=1,size(x);y(i)=pspc_scalar(x(i),a,b,c,d,teta);end do
   end function
   elemental real(dp) function qspc_scalar(p,a,b,c,d,teta) result(x)
      real(dp),intent(in)::p,a,b,c,d,teta;real(dp)::w
      if(p<=0.or.p>=1.or.a<0.or.b<0.or.a+b<=0.or.d<1)then;x=nanv();return;end if
      w=signed_root(normal_quantile(p),d);x=solve_cubic_monotone(w,a,b,c,teta)
   end function
   pure function qspc_vec(p,a,b,c,d,teta) result(x)
      real(dp),intent(in)::p(:),a,b,c,d,teta;real(dp)::x(size(p));integer::i
      do i=1,size(p);x(i)=qspc_scalar(p(i),a,b,c,d,teta);end do
   end function
   subroutine rspc(n,a,b,c,d,teta,x)
      integer,intent(in)::n;real(dp),intent(in)::a,b,c,d,teta;real(dp),intent(out)::x(n);real(dp)::u;integer::i
      if(a<0.or.b<0.or.a+b<=0.or.d<1)then;x=nanv();return;end if
      do i=1,n;call random_number(u);u=max(tiny(1.0_dp),min(1.0_dp-epsilon(1.0_dp),u));x(i)=qspc_scalar(u,a,b,c,d,teta);end do
   end subroutine

   elemental real(dp) function deck_scalar(x,a,p) result(y)
      real(dp),intent(in)::x,a,p;real(dp)::norm
      if(a<=0.or.p<=-1)then;y=nanv();return;end if
      if(abs(x)>a)then;y=0.0_dp;return;end if
      norm=a*sqrt(acos(-1.0_dp))*exp(log_gamma(p+1.0_dp)-log_gamma(p+1.5_dp))
      y=(max(0.0_dp,1.0_dp-(x/a)**2))**p/norm
   end function
   pure function deck_vec(x,a,p) result(y)
      real(dp),intent(in)::x(:),a,p;real(dp)::y(size(x));integer::i
      do i=1,size(x);y(i)=deck_scalar(x(i),a,p);end do
   end function
   elemental real(dp) function peck_scalar(x,a,p) result(y)
      real(dp),intent(in)::x,a,p
      if(a<=0.or.p<=-1)then;y=nanv();return;end if
      if(x<=-a)then;y=0.0_dp;else if(x>=a)then;y=1.0_dp
      else if(abs(x)<=tiny(1.0_dp))then;y=0.5_dp
      else;y=0.5_dp+0.5_dp*sign(1.0_dp,x)*regularized_beta((x/a)**2,0.5_dp,p+1.0_dp);end if
   end function
   pure function peck_vec(x,a,p) result(y)
      real(dp),intent(in)::x(:),a,p;real(dp)::y(size(x));integer::i
      do i=1,size(x);y(i)=peck_scalar(x(i),a,p);end do
   end function
   elemental real(dp) function qeck_scalar(q,a,p) result(x)
      real(dp),intent(in)::q,a,p;real(dp)::z
      if(q<=0.or.q>=1.or.a<=0.or.p<=-1)then;x=nanv();return;end if
      if(q<0.5_dp)then;z=beta_quantile(1.0_dp-2.0_dp*q,0.5_dp,p+1.0_dp);x=-a*sqrt(z)
      else;z=beta_quantile(2.0_dp*q-1.0_dp,0.5_dp,p+1.0_dp);x=a*sqrt(z);end if
   end function
   pure function qeck_vec(q,a,p) result(x)
      real(dp),intent(in)::q(:),a,p;real(dp)::x(size(q));integer::i
      do i=1,size(q);x(i)=qeck_scalar(q(i),a,p);end do
   end function
   subroutine reck(n,a,p,x)
      integer,intent(in)::n;real(dp),intent(in)::a,p;real(dp),intent(out)::x(n);real(dp)::u;integer::i
      if(a<=0.or.p<=-1)then;x=nanv();return;end if
      do i=1,n;call random_number(u);u=max(tiny(1.0_dp),min(1.0_dp-epsilon(1.0_dp),u));x(i)=qeck_scalar(u,a,p);end do
   end subroutine
end module psdistr
