module suppdists_rank
   use suppdists_kinds, only : dp, i8
   use suppdists_special, only : normal_pdf, normal_cdf, beta_inc, beta_pdf, beta_quantile
   use suppdists_stats, only : dist_stats
   use suppdists_friedman_tables, only : get_friedman_table
   implicit none
   private
   public :: dkendall, pkendall, qkendall, rkendall, skendall
   public :: dfriedman, pfriedman, qfriedman, rfriedman, sfriedman
   public :: dspearman, pspearman, qspearman, rspearman, sspearman
   public :: dkruskalwallis, pkruskalwallis, qkruskalwallis, rkruskalwallis, skruskalwallis
   public :: dnormscore, pnormscore, qnormscore, rnormscore, snormscore, norm_order
contains
   pure real(dp) function kendall_exact(n,k,density) result(p)
      integer,intent(in)::n,k
      logical,intent(in)::density
      integer(i8), allocatable :: a(:),b(:)
      integer::m,j,t,maxk
      integer(i8)::s
      if(k<0)then;p=0.0_dp;return;end if
      allocate(a(0:k),b(0:k));a=0_i8;a(0)=1_i8
      do m=2,n
         b=0_i8;maxk=min(k,m*(m-1)/2)
         do j=0,maxk
            s=0_i8
            do t=0,min(m-1,j);s=s+a(j-t);end do
            b(j)=s
         end do
         a=b
      end do
      if(density)then
         s=a(k)
      else
         s=sum(a)
      end if
      p=real(s,dp)/exp(log_gamma(real(n+1,dp)))
   end function kendall_exact

   pure real(dp) function kendall_edge_cdf(n,k) result(p)
      integer,intent(in)::n,k
      real(dp)::dn,mn,s,f,g,sd,l4,l6,x,z,ph3,ph5,ph7
      dn=real(n,dp);mn=dn*(dn-1.0_dp)/4.0_dp
      s=dn*(dn+1.0_dp)*(2.0_dp*dn+1.0_dp)/6.0_dp
      f=((dn+1.0_dp)*3.0_dp*dn-1.0_dp)/5.0_dp
      g=(((dn*dn+2.0_dp)*dn-1.0_dp)*3.0_dp*dn+1.0_dp)/7.0_dp
      sd=s-dn;l4=-1.2_dp*(s*f-dn)/(sd*sd);l6=(48.0_dp/7.0_dp)*(s*g-dn)/(sd**3)
      sd=sqrt(sd/12.0_dp);x=(real(k,dp)+0.5_dp-mn)/sd;z=normal_pdf(x)
      ph3=z*x*(3.0_dp-x*x)
      ph5=-z*x*((x*x-10.0_dp)*x*x+15.0_dp)
      ph7=z*x*(((-x*x+21.0_dp)*x*x-105.0_dp)*x*x+105.0_dp)
      p=normal_cdf(x)+((35.0_dp*l4*l4*ph7/56.0_dp+l6*ph5)/30.0_dp+l4*ph3)/24.0_dp
      p=max(0.0_dp,min(1.0_dp,p))
   end function kendall_edge_cdf

   pure real(dp) function pkendall(tau,n) result(p)
      real(dp),intent(in)::tau; integer,intent(in)::n
      integer::k,m
      if(n<2 .or. tau < -1.0_dp .or. tau>1.0_dp)then;p=0.0_dp;return;end if
      m=n*(n-1)/2;k=nint((1.0_dp+tau)*real(m,dp)/2.0_dp)
      if(k<0)then;p=0.0_dp;else if(k>=m)then;p=1.0_dp
      else if(n<=12)then;p=kendall_exact(n,k,.false.)
      else;p=kendall_edge_cdf(n,k);end if
   end function pkendall

   pure real(dp) function dkendall(tau,n) result(p)
      real(dp),intent(in)::tau;integer,intent(in)::n
      integer::k,m
      m=n*(n-1)/2;k=nint((1.0_dp+tau)*real(m,dp)/2.0_dp)
      if(k<0 .or. k>m)then;p=0.0_dp
      else if(n<=12)then;p=kendall_exact(n,k,.true.)
      else;p=max(0.0_dp,kendall_edge_cdf(n,k)-merge(kendall_edge_cdf(n,k-1),0.0_dp,k>0));end if
   end function dkendall

   pure real(dp) function qkendall(prob,n) result(tau)
      real(dp),intent(in)::prob;integer,intent(in)::n
      integer::k,m
      m=n*(n-1)/2
      do k=0,m
         tau=4.0_dp*real(k,dp)/(real(n,dp)*real(n-1,dp))-1.0_dp
         if(pkendall(tau,n)>=prob)return
      end do
      tau=1.0_dp
   end function qkendall

   real(dp) function rkendall(n) result(tau)
      integer,intent(in)::n;real(dp)::u
      call random_number(u);tau=qkendall(u,n)
   end function rkendall

   function skendall(n) result(s)
      integer,intent(in)::n;type(dist_stats)::s
      integer::k,m;real(dp)::tau,p
      s%variance=real(4*n+10,dp)/real(9*n*(n-1),dp)
      s%mean=0.0_dp;s%median=0.0_dp;s%mode=0.0_dp;s%third_central=0.0_dp
      m=n*(n-1)/2;s%fourth_central=0.0_dp
      do k=0,m
         tau=4.0_dp*real(k,dp)/(real(n,dp)*real(n-1,dp))-1.0_dp;p=dkendall(tau,n)
         s%fourth_central=s%fourth_central+p*tau**4
      end do
   end function skendall

   function friedman_exact_cdf(x,r,n,found) result(p)
      real(dp), intent(in) :: x
      integer, intent(in) :: r,n
      logical, intent(out) :: found
      real(dp) :: p
      integer, allocatable :: st(:)
      real(dp), allocatable :: qt(:)
      integer :: ss,i
      call get_friedman_table(r,n,st,qt,found)
      if(.not.found)then;p=0.0_dp;return;end if
      ss=nint(x*real(n*r*(r+1),dp)/12.0_dp)
      if(mod(r,2)==0)ss=4*ss
      p=1.0_dp
      do i=1,size(st)
         if(st(i)>ss)then
            p=1.0_dp-qt(i)
            return
         end if
      end do
   end function friedman_exact_cdf

   function friedman_exact_pmf(x,r,n,found) result(p)
      real(dp), intent(in) :: x
      integer, intent(in) :: r,n
      logical, intent(out) :: found
      real(dp) :: p
      integer, allocatable :: st(:)
      real(dp), allocatable :: qt(:)
      integer :: ss,i
      call get_friedman_table(r,n,st,qt,found)
      if(.not.found)then;p=0.0_dp;return;end if
      ss=nint(x*real(n*r*(r+1),dp)/12.0_dp)
      if(mod(r,2)==0)ss=4*ss
      p=0.0_dp
      do i=1,size(st)
         if(st(i)==ss)then
            if(i<size(st))then
               p=qt(i)-qt(i+1)
            else
               p=qt(i)
            end if
            return
         end if
      end do
   end function friedman_exact_pmf

   function friedman_exact_quantile(prob,r,n,found) result(x)
      real(dp), intent(in) :: prob
      integer, intent(in) :: r,n
      logical, intent(out) :: found
      real(dp) :: x,cdf,ss0
      integer, allocatable :: st(:)
      real(dp), allocatable :: qt(:)
      integer :: i
      call get_friedman_table(r,n,st,qt,found)
      if(.not.found)then;x=0.0_dp;return;end if
      do i=1,size(st)
         if(i<size(st))then;cdf=1.0_dp-qt(i+1);else;cdf=1.0_dp;end if
         if(cdf>=prob)then
            ss0=real(st(i),dp)
            if(mod(r,2)==0)ss0=ss0/4.0_dp
            x=12.0_dp*ss0/real(n*r*(r+1),dp)
            return
         end if
      end do
      x=real(n*(r-1),dp)
   end function friedman_exact_quantile

   real(dp) function pfriedman(x,r,n) result(p)
      real(dp),intent(in)::x
      integer,intent(in)::r,n
      real(dp)::m,s,w,a,b,pe
      logical :: found
      if(r<3 .or. n<2 .or. x<0.0_dp)then;p=0.0_dp;return;end if
      pe=friedman_exact_cdf(x,r,n,found)
      if(found)then;p=pe;return;end if
      m=real(n*n*r*(r*r-1),dp)/12.0_dp
      s=x*real(n*r*(r+1),dp)/12.0_dp
      if(s>=m)then;p=1.0_dp;return;end if
      s=max(1.0_dp,2.0_dp*real(ceiling(s)/2,dp))
      w=(s-1.0_dp)/(m+2.0_dp);a=real(r-1,dp)-2.0_dp/real(n,dp);b=a*real(n-1,dp)
      p=1.0_dp-beta_inc(1.0_dp-w,b/2.0_dp,a/2.0_dp)
   end function pfriedman

   real(dp) function dfriedman(x,r,n) result(p)
      real(dp),intent(in)::x
      integer,intent(in)::r,n
      real(dp)::step,pe
      logical :: found
      pe=friedman_exact_pmf(x,r,n,found)
      if(found)then;p=pe;return;end if
      step=12.0_dp/real(n*r*(r+1),dp)
      p=max(0.0_dp,pfriedman(x,r,n)-pfriedman(max(0.0_dp,x-step),r,n))
   end function dfriedman

   real(dp) function qfriedman(prob,r,n) result(x)
      real(dp),intent(in)::prob
      integer,intent(in)::r,n
      real(dp)::step,mx,xe
      logical :: found
      xe=friedman_exact_quantile(prob,r,n,found)
      if(found)then;x=xe;return;end if
      step=12.0_dp/real(n*r*(r+1),dp);mx=real(n*(r-1),dp);x=0.0_dp
      do while(x<mx .and. pfriedman(x,r,n)<prob);x=x+step;end do
      x=min(x,mx)
   end function qfriedman

   real(dp) function rfriedman(r,n) result(x)
      integer,intent(in)::r,n;real(dp)::u
      call random_number(u);x=qfriedman(u,r,n)
   end function rfriedman

   function sfriedman(r,n) result(s)
      integer,intent(in)::r,n;type(dist_stats)::s
      real(dp)::rr,nn
      rr=real(r,dp);nn=real(n,dp)
      s%mean=rr-1.0_dp;s%median=qfriedman(0.5_dp,r,n);s%mode=max(0.0_dp,rr-3.0_dp)
      s%variance=2.0_dp*(rr-1.0_dp)*(nn-1.0_dp)/nn
      s%third_central=s%variance*4.0_dp*(nn-2.0_dp)/nn
      s%fourth_central=3.0_dp*s%variance*s%variance
   end function sfriedman

   real(dp) function pspearman(rho,r) result(p)
      real(dp),intent(in)::rho;integer,intent(in)::r
      p=pfriedman((rho+1.0_dp)*real(r-1,dp),r,2)
   end function pspearman
   real(dp) function dspearman(rho,r) result(p)
      real(dp),intent(in)::rho;integer,intent(in)::r
      p=dfriedman((rho+1.0_dp)*real(r-1,dp),r,2)
   end function dspearman
   real(dp) function qspearman(prob,r) result(rho)
      real(dp),intent(in)::prob;integer,intent(in)::r
      rho=qfriedman(prob,r,2)/real(r-1,dp)-1.0_dp
   end function qspearman
   real(dp) function rspearman(r) result(rho)
      integer,intent(in)::r;real(dp)::u
      call random_number(u);rho=qspearman(u,r)
   end function rspearman
   function sspearman(r) result(s)
      integer,intent(in)::r;type(dist_stats)::s
      s%variance=1.0_dp/real(r-1,dp)
      s%mean=0.0_dp;s%median=0.0_dp;s%mode=0.0_dp;s%third_central=0.0_dp
      s%fourth_central=3.0_dp*s%variance*s%variance
   end function sspearman

   pure real(dp) function kw_variance(c,n,u,normal_score) result(v)
      integer,intent(in)::c,n;real(dp),intent(in)::u;logical,intent(in)::normal_score
      real(dp)::cc,nn,alpha,np,nm,nc,cm,den,e,e2,e4,g
      integer::i,m
      cc=real(c,dp);nn=real(n,dp)
      if(.not.normal_score)then
         v=2.0_dp*(cc-1.0_dp)-0.4_dp*(3.0_dp*cc*cc-6.0_dp*cc+ &
            nn*(2.0_dp*cc*cc-6.0_dp*cc+1.0_dp))/(nn*(nn+1.0_dp))-1.2_dp*u
      else
         alpha=0.375_dp;np=nn+1.0_dp;nm=nn-1.0_dp;nc=nn-cc;cm=cc-1.0_dp;den=1.0_dp-2.0_dp*alpha
         m=n/2;e2=0.0_dp;e4=0.0_dp
         do i=1,m
            e=normal_quantile_local((real(i,dp)-alpha)/(nn+den));e=e*e;e2=e2+e;e4=e4+e*e
         end do
         e2=4.0_dp*e2*e2
         g=(np*nn*nm*nm*2.0_dp*e4-3.0_dp*nm**3*e2)/(nm*(nn-2.0_dp)*(nn-3.0_dp)*e2)
         v=2.0_dp*cm*nc/np-g*(np*cc*cc+2.0_dp*cm*nc-nn*np*u)/(nn*np)
      end if
   contains
      pure real(dp) function normal_quantile_local(p) result(x)
         use suppdists_special, only: normal_quantile
         real(dp),intent(in)::p;x=normal_quantile(p)
      end function normal_quantile_local
   end function kw_variance

   pure real(dp) function pkw(x,c,n,u,normal_score) result(p)
      real(dp),intent(in)::x,u;integer,intent(in)::c,n;logical,intent(in)::normal_score
      real(dp)::v,d,a,b,nn,cc
      if(x<=0.0_dp)then;p=0.0_dp;return;end if
      nn=real(n,dp);cc=real(c,dp);v=kw_variance(c,n,u,normal_score)
      d=((nn-cc)*(cc-1.0_dp)-v)/((nn-1.0_dp)*v);a=(cc-1.0_dp)*d;b=(nn-cc)*d
      p=beta_inc(min(1.0_dp,x/(nn-1.0_dp)),a,b)
   end function pkw
   pure real(dp) function dkw(x,c,n,u,normal_score) result(p)
      real(dp),intent(in)::x,u;integer,intent(in)::c,n;logical,intent(in)::normal_score
      real(dp)::v,d,a,b,nn,cc,z
      nn=real(n,dp);cc=real(c,dp);v=kw_variance(c,n,u,normal_score)
      d=((nn-cc)*(cc-1.0_dp)-v)/((nn-1.0_dp)*v);a=(cc-1.0_dp)*d;b=(nn-cc)*d;z=x/(nn-1.0_dp)
      p=beta_pdf(z,a,b)/(nn-1.0_dp)
   end function dkw
   pure real(dp) function qkw(prob,c,n,u,normal_score) result(x)
      real(dp),intent(in)::prob,u;integer,intent(in)::c,n;logical,intent(in)::normal_score
      real(dp)::v,d,a,b,nn,cc
      nn=real(n,dp);cc=real(c,dp);v=kw_variance(c,n,u,normal_score)
      d=((nn-cc)*(cc-1.0_dp)-v)/((nn-1.0_dp)*v);a=(cc-1.0_dp)*d;b=(nn-cc)*d
      x=(nn-1.0_dp)*beta_quantile(prob,a,b)
   end function qkw
   real(dp) function rkw(c,n,u,normal_score) result(x)
      integer,intent(in)::c,n;real(dp),intent(in)::u;logical,intent(in)::normal_score;real(dp)::p
      call random_number(p);x=qkw(p,c,n,u,normal_score)
   end function rkw
   function skw(c,n,u,normal_score) result(s)
      integer,intent(in)::c,n;real(dp),intent(in)::u;logical,intent(in)::normal_score;type(dist_stats)::s
      real(dp)::v,d,a,b,scale,meanb,varb,skew,kurt
      v=kw_variance(c,n,u,normal_score);scale=real(n-1,dp);s%mean=real(c-1,dp);s%variance=v
      s%median=qkw(0.5_dp,c,n,u,normal_score)
      d=((real(n-c,dp)*real(c-1,dp)-v)/(real(n-1,dp)*v));a=real(c-1,dp)*d;b=real(n-c,dp)*d
      if(a>1 .and. b>1)then;s%mode=scale*(a-1)/(a+b-2);else;s%mode=0;end if
      meanb=a/(a+b);varb=a*b/((a+b)**2*(a+b+1))
      skew=2*(b-a)*sqrt(a+b+1)/((a+b+2)*sqrt(a*b))
      kurt=6*((a-b)**2*(a+b+1)-a*b*(a+b+2))/(a*b*(a+b+2)*(a+b+3))
      s%third_central=skew*(scale*sqrt(varb))**3
      s%fourth_central=(kurt+3.0_dp)*(scale*scale*varb)**2
   end function skw

   pure real(dp) function pkruskalwallis(x,c,n,u) result(p)
      real(dp), intent(in) :: x,u
      integer, intent(in) :: c,n
      p=pkw(x,c,n,u,.false.)
   end function pkruskalwallis

   pure real(dp) function dkruskalwallis(x,c,n,u) result(p)
      real(dp), intent(in) :: x,u
      integer, intent(in) :: c,n
      p=dkw(x,c,n,u,.false.)
   end function dkruskalwallis

   pure real(dp) function qkruskalwallis(prob,c,n,u) result(x)
      real(dp), intent(in) :: prob,u
      integer, intent(in) :: c,n
      x=qkw(prob,c,n,u,.false.)
   end function qkruskalwallis

   real(dp) function rkruskalwallis(c,n,u) result(x)
      integer, intent(in) :: c,n
      real(dp), intent(in) :: u
      x=rkw(c,n,u,.false.)
   end function rkruskalwallis

   function skruskalwallis(c,n,u) result(s)
      integer, intent(in) :: c,n
      real(dp), intent(in) :: u
      type(dist_stats) :: s
      s=skw(c,n,u,.false.)
   end function skruskalwallis

   pure real(dp) function pnormscore(x,c,n,u) result(p)
      real(dp), intent(in) :: x,u
      integer, intent(in) :: c,n
      p=pkw(x,c,n,u,.true.)
   end function pnormscore

   pure real(dp) function dnormscore(x,c,n,u) result(p)
      real(dp), intent(in) :: x,u
      integer, intent(in) :: c,n
      p=dkw(x,c,n,u,.true.)
   end function dnormscore

   pure real(dp) function qnormscore(prob,c,n,u) result(x)
      real(dp), intent(in) :: prob,u
      integer, intent(in) :: c,n
      x=qkw(prob,c,n,u,.true.)
   end function qnormscore

   real(dp) function rnormscore(c,n,u) result(x)
      integer, intent(in) :: c,n
      real(dp), intent(in) :: u
      x=rkw(c,n,u,.true.)
   end function rnormscore

   function snormscore(c,n,u) result(s)
      integer, intent(in) :: c,n
      real(dp), intent(in) :: u
      type(dist_stats) :: s
      s=skw(c,n,u,.true.)
   end function snormscore

   subroutine norm_order(n,s)
      use suppdists_special, only : normal_quantile
      integer,intent(in)::n;real(dp),allocatable,intent(out)::s(:)
      real(dp),parameter::eps(4)=[.419885_dp,.450536_dp,.456936_dp,.468488_dp]
      real(dp),parameter::dl1(4)=[.112063_dp,.12177_dp,.239299_dp,.215159_dp]
      real(dp),parameter::dl2(4)=[.080122_dp,.111348_dp,-.211867_dp,-.115049_dp]
      real(dp),parameter::gam(4)=[.474798_dp,.469051_dp,.208597_dp,.259784_dp]
      real(dp),parameter::lam(4)=[.282765_dp,.304856_dp,.407708_dp,.414093_dp]
      real(dp),parameter::bb=-.283833_dp,d=-.106136_dp
      real(dp),allocatable::half(:);real(dp)::an,ai,e1,e2;integer::m,i,j
      m=n/2;allocate(half(m));an=real(n,dp)
      do i=1,m
         j=min(i,4);ai=real(i,dp);e1=(ai-eps(j))/(an+gam(j))
         if(i<=3)then;e2=e1**lam(j);else;e2=e1**(lam(4)+bb/(ai+d));end if
         half(i)=e1+e2*(dl1(j)+e2*dl2(j))/an-correc(i,n)
         half(i)=-normal_quantile(half(i))
      end do
      allocate(s(n))
      do i=1,m;s(i)=-half(i);s(n-i+1)=half(i);end do
      if(mod(n,2)==1)s(m+1)=0.0_dp
   contains
      pure real(dp) function correc(i,n) result(c)
         integer,intent(in)::i,n
         real(dp),parameter::c1(7)=[9.5_dp,28.7_dp,1.9_dp,0.0_dp,-7.0_dp,-6.2_dp,-1.6_dp]
         real(dp),parameter::c2(7)=[-6195.0_dp,-9569.0_dp,-6728.0_dp,-17614.0_dp,-8278.0_dp,-3570.0_dp,1075.0_dp]
         real(dp),parameter::c3(7)=[93380.0_dp,175160.0_dp,410400.0_dp,2157600.0_dp,2376000.0_dp,2065000.0_dp,2065000.0_dp]
         real(dp)::an
         if(i*n==4)then;c=1.9e-5_dp;return;end if
         if(i<1 .or. i>7 .or. (i/=4 .and. n>20) .or. (i==4 .and. n>40))then;c=0.0_dp;return;end if
         an=1.0_dp/real(n*n,dp);c=(c1(i)+an*(c2(i)+an*c3(i)))*1.0e-6_dp
      end function correc
   end subroutine norm_order
end module suppdists_rank
