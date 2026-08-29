module urca_regression
   use urca_kinds, only : dp
   use urca_types, only : lm_result
   use urca_linalg, only : invert_spd
   implicit none
   private
   public :: lm_fit, lm_fit_multi, add_intercept, coefficient_t
   public :: nested_f_stat, aic_lm, bic_lm

   interface
      subroutine dgels(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
         import dp
         character(len=1), intent(in) :: trans
         integer,intent(in)::m,n,nrhs,lda,ldb,lwork
         real(dp),intent(inout)::a(lda,*),b(ldb,*),work(*)
         integer,intent(out)::info
      end subroutine dgels
   end interface
contains
   function add_intercept(x) result(z)
      real(dp),intent(in)::x(:,:)
      real(dp),allocatable::z(:,:)
      allocate(z(size(x,1),size(x,2)+1))
      z(:,1)=1.0_dp
      z(:,2:)=x
   end function add_intercept

   function lm_fit(x,y) result(fit)
      real(dp),intent(in)::x(:,:),y(:)
      type(lm_result)::fit
      real(dp),allocatable::a(:,:),b(:,:),work(:),xtx(:,:),xtxi(:,:)
      real(dp)::q(1),pi
      integer::n,p,ldb,lwork,ii
      n=size(x,1)
      p=size(x,2)
      fit%nobs=n
      fit%rank=p
      fit%df_resid=n-p
      if(size(y)/=n .or. n<p .or. p<1) then
      fit%info=-1
      return
      end if
      ldb=max(n,p)
      allocate(a(n,p),b(ldb,1))
      a=x
      b=0
      b(1:n,1)=y
      call dgels('N',n,p,1,a,n,b,ldb,q,-1,fit%info)
      if(fit%info/=0)return
      lwork=max(1,int(q(1)))
      allocate(work(lwork))
      a=x
      b=0
      b(1:n,1)=y
      call dgels('N',n,p,1,a,n,b,ldb,work,lwork,fit%info)
      if(fit%info/=0)return
      allocate(fit%beta(p),fit%fitted(n),fit%residuals(n),fit%vcov(p,p))
      fit%beta=b(1:p,1)
      fit%fitted=matmul(x,fit%beta)
      fit%residuals=y-fit%fitted
      fit%rss=sum(fit%residuals**2)
      if(fit%df_resid>0)fit%sigma2=fit%rss/real(fit%df_resid,dp)
      allocate(xtx(p,p))
      xtx=matmul(transpose(x),x)
      call invert_spd(xtx,xtxi,ii)
      if(ii==0)then
      fit%vcov=fit%sigma2*xtxi
      else
      fit%vcov=0
      fit%info=ii
      end if
      pi=acos(-1.0_dp)
      if(fit%rss>0)fit%loglik=-0.5_dp*real(n,dp)*(log(2*pi)+1+log(fit%rss/real(n,dp)))
   end function lm_fit

   subroutine lm_fit_multi(x,y,beta,residuals,sigma,info)
      real(dp),intent(in)::x(:,:),y(:,:)
      real(dp),allocatable,intent(out)::beta(:,:),residuals(:,:),sigma(:,:)
      integer,intent(out)::info
      real(dp),allocatable::a(:,:),b(:,:),work(:)
      real(dp)::q(1)
      integer::n,p,m,ldb,lwork
      n=size(x,1)
      p=size(x,2)
      m=size(y,2)
      ldb=max(n,p)
      if(size(y,1)/=n .or. n<p)then
      info=-1
      allocate(beta(0,0),residuals(0,0),sigma(0,0))
      return
      end if
      allocate(a(n,p),b(ldb,m))
      a=x
      b=0
      b(1:n,:)=y
      call dgels('N',n,p,m,a,n,b,ldb,q,-1,info)
      if(info/=0)return
      lwork=max(1,int(q(1)))
      allocate(work(lwork))
      a=x
      b=0
      b(1:n,:)=y
      call dgels('N',n,p,m,a,n,b,ldb,work,lwork,info)
      if(info/=0)return
      allocate(beta(p,m),residuals(n,m),sigma(m,m))
      beta=b(1:p,:)
      residuals=y-matmul(x,beta)
      sigma=matmul(transpose(residuals),residuals)/real(n,dp)
   end subroutine lm_fit_multi

   real(dp) function coefficient_t(fit,j) result(t)
      type(lm_result),intent(in)::fit
      integer,intent(in)::j
      if(j<1 .or. j>size(fit%beta) .or. fit%vcov(j,j)<=0)then
      t=0
      else
      t=fit%beta(j)/sqrt(fit%vcov(j,j))
      end if
   end function coefficient_t

   real(dp) function nested_f_stat(rss_r,rss_u,df_r,df_u) result(f)
      real(dp),intent(in)::rss_r,rss_u
      integer,intent(in)::df_r,df_u
      integer::q
      q=df_r-df_u
      if(q<=0 .or. df_u<=0 .or. rss_u<=0)then
      f=0
      else
      f=((rss_r-rss_u)/real(q,dp))/(rss_u/real(df_u,dp))
      end if
   end function nested_f_stat

   real(dp) function aic_lm(fit,kpen) result(aic)
      type(lm_result),intent(in)::fit
      real(dp),intent(in),optional::kpen
      real(dp)::k
      k=2.0_dp
      if(present(kpen))k=kpen
      aic=-2.0_dp*fit%loglik+k*real(fit%rank+1,dp)
   end function aic_lm
   real(dp) function bic_lm(fit) result(bic)
      type(lm_result),intent(in)::fit
      bic=-2.0_dp*fit%loglik+log(real(fit%nobs,dp))*real(fit%rank+1,dp)
   end function bic_lm
end module urca_regression
