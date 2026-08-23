module rfast_score_tests
   use rfast_special, only : dp, digamma_r, trigamma_r, reg_gamma_q, student_t_cdf
   use rfast_arrays, only : mean_r, variance_r
   use rfast_linalg, only : solve_linear
   use rfast_mle, only : mle_result, beta_mle, gamma_mle, geometric_mle, invgauss_mle, negbin_mle, ztp_mle
   use rfast_extra_mle, only : weibull_mle
   implicit none
   private

   type, public :: score_result
      real(dp), allocatable :: statistic(:)
      real(dp), allocatable :: pvalue(:)
   end type score_result

   public :: score_betaregs, score_expregs, score_gammaregs, score_geomregs
   public :: score_glms, score_invgaussregs, score_multinomregs, score_negbinregs
   public :: score_weibregs, score_ztpregs

contains

   pure real(dp) function chisq_upper(x, df, logged) result(p)
      real(dp), intent(in) :: x, df
      logical, intent(in) :: logged
      real(dp) :: q
      q = reg_gamma_q(0.5_dp*df, 0.5_dp*max(0.0_dp,x))
      if (logged) then
         p = log(max(tiny(1.0_dp),q))
      else
         p = q
      end if
   end function chisq_upper

   pure real(dp) function t_twosided(t, df, logged) result(p)
      real(dp), intent(in) :: t, df
      logical, intent(in) :: logged
      real(dp) :: q
      q = 2.0_dp * max(0.0_dp, 1.0_dp - student_t_cdf(abs(t),df))
      q = min(1.0_dp,q)
      if (logged) then
         p = log(max(tiny(1.0_dp),q))
      else
         p = q
      end if
   end function t_twosided

   function score_betaregs(y,x,logged) result(res)
      real(dp), intent(in) :: y(:), x(:,:)
      logical, intent(in), optional :: logged
      type(score_result) :: res
      type(mle_result) :: fit
      real(dp) :: m1,m2,u(size(x,2)),den(size(x,2)),z(size(y))
      integer :: j
      logical :: lg
      lg=.false.;if(present(logged))lg=logged
      fit=beta_mle(y)
      allocate(res%statistic(size(x,2)),res%pvalue(size(x,2)))
      if(fit%status/=0.or..not.allocated(fit%param))then
         res%statistic=huge(1.0_dp);res%pvalue=merge(-huge(1.0_dp),0.0_dp,lg);return
      end if
      m1=digamma_r(fit%param(1))-digamma_r(fit%param(2))
      m2=trigamma_r(fit%param(1))+trigamma_r(fit%param(2))
      z=log(y)-log(1.0_dp-y)-m1
      do j=1,size(x,2)
         u(j)=dot_product(x(:,j),z);den(j)=sum(x(:,j)**2)*m2
         res%statistic(j)=u(j)*u(j)/max(tiny(1.0_dp),den(j))
         res%pvalue(j)=chisq_upper(res%statistic(j),1.0_dp,lg)
      end do
   end function score_betaregs

   function score_expregs(y,x,logged) result(res)
      real(dp), intent(in) :: y(:), x(:,:)
      logical, intent(in), optional :: logged
      type(score_result) :: res
      real(dp) :: lam,u,vu
      integer :: j
      logical :: lg
      lg=.false.;if(present(logged))lg=logged
      lam=mean_r(y);allocate(res%statistic(size(x,2)),res%pvalue(size(x,2)))
      do j=1,size(x,2)
         u=dot_product(x(:,j),y)*lam-sum(x(:,j));vu=sum(x(:,j)**2)*lam**4
         res%statistic(j)=u*u/max(tiny(1.0_dp),vu);res%pvalue(j)=chisq_upper(res%statistic(j),1.0_dp,lg)
      end do
   end function score_expregs

   function score_gammaregs(y,x,logged) result(res)
      real(dp), intent(in) :: y(:), x(:,:)
      logical, intent(in), optional :: logged
      type(score_result) :: res
      type(mle_result)::fit
      real(dp)::m,u,vb
      integer::j
      logical::lg
      lg=.false.;if(present(logged))lg=logged;fit=gamma_mle(y);allocate(res%statistic(size(x,2)),res%pvalue(size(x,2)))
      if(fit%status/=0.or..not.allocated(fit%param))then;res%statistic=huge(1.0_dp);res%pvalue=0.0_dp;return;end if
      m=fit%param(1)/fit%param(2)
      do j=1,size(x,2)
         u=sum(x(:,j))-dot_product(x(:,j),y)/m;vb=sum(x(:,j)**2)/fit%param(1)
         res%statistic(j)=u*u/max(tiny(1.0_dp),vb);res%pvalue(j)=chisq_upper(res%statistic(j),1.0_dp,lg)
      end do
   end function score_gammaregs

   function score_geomregs(y,x,logged) result(res)
      integer, intent(in) :: y(:)
      real(dp), intent(in) :: x(:,:)
      logical, intent(in), optional :: logged
      type(score_result)::res
      type(mle_result)::fit
      real(dp)::p,u,vb
      integer::j
      logical::lg
      lg=.false.;if(present(logged))lg=logged;fit=geometric_mle(y,1);p=fit%param(1)
      allocate(res%statistic(size(x,2)),res%pvalue(size(x,2)))
      do j=1,size(x,2)
         u=(1.0_dp-p)*sum(x(:,j))-p*sum(real(y,dp)*x(:,j))
         vb=(1.0_dp-p)*sum(x(:,j)**2);res%statistic(j)=u*u/max(tiny(1.0_dp),vb)
         res%pvalue(j)=chisq_upper(res%statistic(j),1.0_dp,lg)
      end do
   end function score_geomregs

   function score_glms(y,x,binomial,logged) result(res)
      real(dp), intent(in) :: y(:), x(:,:)
      logical, intent(in), optional :: binomial, logged
      type(score_result)::res
      real(dp)::my,sy,mx,sx,r,stat,fac
      integer::j,n
      logical::bn,lg
      n=size(y);bn=.false.;if(present(binomial))bn=binomial;lg=.false.;if(present(logged))lg=logged
      my=mean_r(y);sy=sqrt(max(0.0_dp,variance_r(y)));allocate(res%statistic(size(x,2)),res%pvalue(size(x,2)))
      if(.not.bn)fac=sy/sqrt(max(tiny(1.0_dp),sum(y)/real(n,dp)))*sqrt(real(n-1,dp))
      do j=1,size(x,2)
         mx=sum(x(:,j))/real(n,dp);sx=sqrt(sum((x(:,j)-mx)**2)/max(1.0_dp,real(n-1,dp)))
         if(sy<=0.0_dp.or.sx<=0.0_dp)then;r=0.0_dp;else;r=sum((y-my)*(x(:,j)-mx))/real(n-1,dp)/(sy*sx);end if
         if(bn)then;stat=r*sqrt(real(n,dp));else;stat=fac*r;end if
         res%statistic(j)=stat;res%pvalue(j)=t_twosided(stat,real(n-2,dp),lg)
      end do
   end function score_glms

   function score_invgaussregs(y,x,logged) result(res)
      real(dp), intent(in) :: y(:), x(:,:)
      logical, intent(in), optional :: logged
      type(score_result)::res
      real(dp)::m,lambda,u,vu
      integer::j,n
      logical::lg
      n=size(y);lg=.false.;if(present(logged))lg=logged;m=sum(y)/real(n,dp)
      lambda=1.0_dp/(sum(1.0_dp/y)/real(n,dp)-1.0_dp/m)
      allocate(res%statistic(size(x,2)),res%pvalue(size(x,2)))
      do j=1,size(x,2)
         u=dot_product(x(:,j),m-y)*lambda;vu=m**3*sum(x(:,j)**2)
         res%statistic(j)=u*u/max(tiny(1.0_dp),vu);res%pvalue(j)=chisq_upper(res%statistic(j),1.0_dp,lg)
      end do
   end function score_invgaussregs

   function score_multinomregs(y,x,logged) result(res)
      integer, intent(in) :: y(:)
      real(dp), intent(in) :: x(:,:)
      logical, intent(in), optional :: logged
      type(score_result)::res
      integer::k,dof,j,c,n,info
      real(dp),allocatable::m(:),vp(:,:),u(:),sol(:)
      real(dp)::sx,sx2,stat
      logical::lg
      n=size(y);k=maxval(y);dof=k-1;lg=.false.;if(present(logged))lg=logged
      if(dof<=1)then
         block
            real(dp)::yy(n)
            yy=merge(1.0_dp,0.0_dp,y==k)
            res=score_glms(yy,x,.true.,lg)
         end block
         return
      end if
      allocate(m(dof),vp(dof,dof),u(dof),sol(dof),res%statistic(size(x,2)),res%pvalue(size(x,2)))
      do c=1,dof;m(c)=real(count(y==c+1),dp)/real(n,dp);end do
      vp=0.0_dp
      do c=1,dof;vp(c,c)=m(c);end do
      do c=1,dof;vp(c,:)=vp(c,:)-m(c)*m;end do
      do j=1,size(x,2)
         sx=sum(x(:,j));sx2=sum(x(:,j)**2);u=0.0_dp
         do c=1,dof;u(c)=sum(x(:,j),mask=(y==c+1))-sx*m(c);end do
         call solve_linear(vp,u,sol,info)
         if(info/=0.or.sx2<=0.0_dp)then;stat=huge(1.0_dp);else;stat=dot_product(u,sol)/sx2;end if
         res%statistic(j)=stat;res%pvalue(j)=chisq_upper(stat,real(dof,dp),lg)
      end do
   end function score_multinomregs

   function score_negbinregs(y,x,logged) result(res)
      integer, intent(in) :: y(:)
      real(dp), intent(in) :: x(:,:)
      logical, intent(in), optional :: logged
      type(score_result)::res
      type(mle_result)::fit
      real(dp)::r,p,my,u,vu
      integer::j
      logical::lg
      lg=.false.;if(present(logged))lg=logged;fit=negbin_mle(y);p=fit%param(1);r=fit%param(2);my=fit%param(3)
      allocate(res%statistic(size(x,2)),res%pvalue(size(x,2)))
      do j=1,size(x,2)
         u=p*sum(real(y,dp)*x(:,j))-(1.0_dp-p)*r*sum(x(:,j))
         vu=sum(x(:,j)**2)*(p*p*(my+my*my/r));res%statistic(j)=u*u/max(tiny(1.0_dp),vu)
         res%pvalue(j)=chisq_upper(res%statistic(j),1.0_dp,lg)
      end do
   end function score_negbinregs

   function score_weibregs(y,x,logged) result(res)
      real(dp), intent(in) :: y(:),x(:,:)
      logical, intent(in), optional :: logged
      type(score_result)::res
      type(mle_result)::fit
      real(dp)::k,lam,u,vu
      integer::j
      logical::lg
      lg=.false.;if(present(logged))lg=logged;fit=weibull_mle(y);k=fit%param(1);lam=fit%param(2)
      allocate(res%statistic(size(x,2)),res%pvalue(size(x,2)))
      do j=1,size(x,2)
         u=sum(x(:,j)*y**k)/lam**k-sum(x(:,j));vu=sum(x(:,j)**2)
         res%statistic(j)=u*u/max(tiny(1.0_dp),vu);res%pvalue(j)=chisq_upper(res%statistic(j),1.0_dp,lg)
      end do
   end function score_weibregs

   function score_ztpregs(y,x,logged) result(res)
      integer, intent(in) :: y(:)
      real(dp), intent(in) :: x(:,:)
      logical, intent(in), optional :: logged
      type(score_result)::res
      type(mle_result)::fit
      real(dp)::lam,elam,ey,u,vu
      integer::j
      logical::lg
      lg=.false.;if(present(logged))lg=logged;fit=ztp_mle(y);lam=fit%param(1);elam=exp(lam);ey=lam*elam/(elam-1.0_dp)
      allocate(res%statistic(size(x,2)),res%pvalue(size(x,2)))
      do j=1,size(x,2)
         u=sum(real(y,dp)*x(:,j))-ey*sum(x(:,j));vu=sum(x(:,j)**2)*(ey*(1.0_dp+lam-ey))
         res%statistic(j)=u*u/max(tiny(1.0_dp),vu);res%pvalue(j)=chisq_upper(res%statistic(j),1.0_dp,lg)
      end do
   end function score_ztpregs

end module rfast_score_tests
