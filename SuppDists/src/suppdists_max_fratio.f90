module suppdists_max_fratio
   use suppdists_kinds, only : dp
   use suppdists_special, only : chisq_pdf, chisq_cdf, rand_chisq, adaptive_integral
   use suppdists_stats, only : dist_stats
   implicit none
   private
   public :: dmaxfratio, pmaxfratio, qmaxfratio, rmaxfratio, smaxfratio
contains
   real(dp) function pmaxfratio(ratio,df,k) result(p)
      real(dp),intent(in)::ratio;integer,intent(in)::df,k
      if(ratio<1.0_dp .or. df<1 .or. k<1)then;p=0.0_dp;return;end if
      if(k==1)then;p=1.0_dp;return;end if
      p=real(k,dp)*adaptive_integral(fun,0.0_dp,1.0_dp-1e-10_dp,2e-8_dp)
      p=max(0.0_dp,min(1.0_dp,p))
   contains
      function fun(t) result(y)
         real(dp),intent(in)::t;real(dp)::y,x,diff,jac
         if(t<=0.0_dp .or. t>=1.0_dp)then;y=0.0_dp;return;end if
         x=t/(1.0_dp-t);jac=1.0_dp/(1.0_dp-t)**2
         diff=max(0.0_dp,chisq_cdf(ratio*x,real(df,dp))-chisq_cdf(x,real(df,dp)))
         y=chisq_pdf(x,real(df,dp))*diff**(k-1)*jac
      end function fun
   end function pmaxfratio

   real(dp) function dmaxfratio(ratio,df,k) result(f)
      real(dp),intent(in)::ratio;integer,intent(in)::df,k
      if(ratio<=1.0_dp .or. df<1 .or. k<2)then;f=0.0_dp;return;end if
      f=real(k*(k-1),dp)*adaptive_integral(fun,0.0_dp,1.0_dp-1e-10_dp,3e-8_dp)
   contains
      function fun(t) result(y)
         real(dp),intent(in)::t;real(dp)::y,x,diff,jac
         if(t<=0.0_dp .or. t>=1.0_dp)then;y=0.0_dp;return;end if
         x=t/(1.0_dp-t);jac=1.0_dp/(1.0_dp-t)**2
         diff=max(0.0_dp,chisq_cdf(ratio*x,real(df,dp))-chisq_cdf(x,real(df,dp)))
         y=chisq_pdf(x,real(df,dp))*diff**(k-2)*chisq_pdf(ratio*x,real(df,dp))*x*jac
      end function fun
   end function dmaxfratio

   real(dp) function qmaxfratio(prob,df,k) result(x)
      real(dp),intent(in)::prob;integer,intent(in)::df,k
      real(dp)::lo,hi,mid;integer::i
      if(prob<=0.0_dp)then;x=1.0_dp;return;end if
      if(prob>=1.0_dp)then;x=huge(1.0_dp);return;end if
      lo=1.0_dp;hi=2.0_dp
      do while(pmaxfratio(hi,df,k)<prob .and. hi<1e8_dp);hi=hi*2.0_dp;end do
      do i=1,70
         mid=0.5_dp*(lo+hi)
         if(pmaxfratio(mid,df,k)<prob)then;lo=mid;else;hi=mid;end if
      end do
      x=0.5_dp*(lo+hi)
   end function qmaxfratio

   real(dp) function rmaxfratio(df,k) result(r)
      integer,intent(in)::df,k
      real(dp)::x,mn,mx;integer::i
      mn=huge(1.0_dp);mx=0.0_dp
      do i=1,k;x=rand_chisq(real(df,dp));mn=min(mn,x);mx=max(mx,x);end do
      r=mx/mn
   end function rmaxfratio

   function smaxfratio(df,k) result(s)
      integer,intent(in)::df,k;type(dist_stats)::s
      real(dp)::lo,hi,norm
      lo=qmaxfratio(1e-4_dp,df,k);hi=qmaxfratio(0.9999_dp,df,k)
      norm=adaptive_integral(pdf,lo,hi,1e-7_dp)
      s%mean=adaptive_integral(m1,lo,hi,1e-7_dp)/norm
      s%median=qmaxfratio(0.5_dp,df,k);s%mode=mode_search(lo,hi)
      s%variance=adaptive_integral(m2,lo,hi,1e-7_dp)/norm
      s%third_central=adaptive_integral(m3,lo,hi,1e-7_dp)/norm
      s%fourth_central=adaptive_integral(m4,lo,hi,1e-7_dp)/norm
   contains
      function pdf(x) result(y);real(dp),intent(in)::x;real(dp)::y;y=dmaxfratio(x,df,k);end function
      function m1(x) result(y);real(dp),intent(in)::x;real(dp)::y;y=x*dmaxfratio(x,df,k);end function
      function m2(x) result(y);real(dp),intent(in)::x;real(dp)::y;y=(x-s%mean)**2*dmaxfratio(x,df,k);end function
      function m3(x) result(y);real(dp),intent(in)::x;real(dp)::y;y=(x-s%mean)**3*dmaxfratio(x,df,k);end function
      function m4(x) result(y);real(dp),intent(in)::x;real(dp)::y;y=(x-s%mean)**4*dmaxfratio(x,df,k);end function
      function mode_search(a,b) result(xm)
         real(dp),intent(in)::a,b;real(dp)::xm,x,y,best;integer::j
         best=-1.0_dp;xm=a
         do j=0,500;x=a+(b-a)*real(j,dp)/500.0_dp;y=dmaxfratio(x,df,k);if(y>best)then;best=y;xm=x;end if;end do
      end function mode_search
   end function smaxfratio
end module suppdists_max_fratio
