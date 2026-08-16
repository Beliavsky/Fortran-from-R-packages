module rfast_regression
   use rfast_special, only : dp, pi, log1p_r
   use rfast_arrays, only : mean_r
   use rfast_linalg, only : solve_linear, inverse_matrix
   implicit none
   private
   type, public :: regression_result
      real(dp), allocatable :: beta(:), fitted(:), residuals(:), covariance(:,:)
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: deviance = huge(1.0_dp)
      real(dp) :: dispersion = 1.0_dp
      integer :: iterations = 0
      integer :: status = 0
   end type regression_result
   public :: lmfit, glm_logistic, glm_poisson, poisson_only, logistic_only
   public :: weighted_least_squares, ridge_regression, ar1_fit

contains

   function weighted_least_squares(x,y,w) result(res)
      real(dp),intent(in)::x(:,:),y(:),w(:);type(regression_result)::res
      real(dp)::xtwx(size(x,2),size(x,2)),xtwy(size(x,2)),b(size(x,2)),inv(size(x,2),size(x,2))
      real(dp)::wx(size(x,1),size(x,2)),rss;integer::j,info,n,p
      n=size(x,1);p=size(x,2)
      do j=1,p;wx(:,j)=w*x(:,j);end do
      xtwx=matmul(transpose(x),wx);xtwy=matmul(transpose(x),w*y);call solve_linear(xtwx,xtwy,b,info)
      allocate(res%beta(p),res%fitted(n),res%residuals(n),res%covariance(p,p));res%status=info
      if(info/=0)then;res%beta=0;res%fitted=0;res%residuals=y;res%covariance=0;return;end if
      res%beta=b;res%fitted=matmul(x,b);res%residuals=y-res%fitted;rss=sum(w*res%residuals**2)
      res%dispersion=rss/max(1.0_dp,real(n-p,dp));call inverse_matrix(xtwx,inv,info);res%covariance=res%dispersion*inv
      res%iterations=1
   end function weighted_least_squares

   function lmfit(x,y,intercept) result(res)
      real(dp),intent(in)::x(:,:),y(:);logical,intent(in),optional::intercept;type(regression_result)::res
      logical::inc;real(dp),allocatable::xx(:,:),w(:);integer::n,p
      n=size(x,1);p=size(x,2);inc=.false.;if(present(intercept))inc=intercept
      if(inc)then;allocate(xx(n,p+1));xx(:,1)=1.0_dp;xx(:,2:)=x;else;allocate(xx(n,p));xx=x;end if
      allocate(w(n));w=1.0_dp;res=weighted_least_squares(xx,y,w)
      res%deviance=sum(res%residuals**2)
      if(res%deviance>0)res%loglik=-0.5_dp*n*(log(2.0_dp*pi*res%deviance/n)+1.0_dp)
   end function lmfit

   function glm_logistic(x,y,tol,maxiter) result(res)
      real(dp),intent(in)::x(:,:),y(:);real(dp),intent(in),optional::tol;integer,intent(in),optional::maxiter
      type(regression_result)::res;integer::n,p,i,mi,j,info;real(dp)::eps,eta(size(y)),mu(size(y)),w(size(y)),z(size(y))
      real(dp)::b(size(x,2)),bold(size(x,2)),xtwx(size(x,2),size(x,2)),xtwz(size(x,2)),wx(size(x,1),size(x,2))
      real(dp)::inv(size(x,2),size(x,2)),ll
      n=size(x,1);p=size(x,2);eps=1e-8_dp;if(present(tol))eps=tol;mi=100;if(present(maxiter))mi=maxiter;b=0.0_dp
      do i=1,mi
         bold = b
         eta = matmul(x,b)
         where (eta >= 0.0_dp)
            mu = 1.0_dp/(1.0_dp + exp(-eta))
         elsewhere
            mu = exp(eta)/(1.0_dp + exp(eta))
         end where
         mu=max(1e-12_dp,min(1.0_dp-1e-12_dp,mu));w=mu*(1-mu);z=eta+(y-mu)/w
         do j=1,p;wx(:,j)=w*x(:,j);end do;xtwx=matmul(transpose(x),wx);xtwz=matmul(transpose(x),w*z)
         call solve_linear(xtwx,xtwz,b,info);if(info/=0)exit;if(maxval(abs(b-bold))<eps)exit
      end do
      allocate(res%beta(p),res%fitted(n),res%residuals(n),res%covariance(p,p));res%beta=b;res%iterations=i;res%status=info
      eta = matmul(x,b)
      where (eta >= 0.0_dp)
         mu = 1.0_dp/(1.0_dp + exp(-eta))
      elsewhere
         mu = exp(eta)/(1.0_dp + exp(eta))
      end where
      mu=max(1e-15_dp,min(1.0_dp-1e-15_dp,mu));res%fitted=mu;res%residuals=y-mu
      ll=sum(y*log(mu)+(1-y)*log1p_r(-mu));res%loglik=ll;res%deviance=-2.0_dp*ll;w=mu*(1-mu)
      do j=1,p;wx(:,j)=w*x(:,j);end do;xtwx=matmul(transpose(x),wx);call inverse_matrix(xtwx,inv,info);res%covariance=inv
   end function glm_logistic

   function glm_poisson(x,y,tol,maxiter) result(res)
      real(dp),intent(in)::x(:,:),y(:);real(dp),intent(in),optional::tol;integer,intent(in),optional::maxiter
      type(regression_result)::res;integer::n,p,i,mi,j,info,k;real(dp)::eps,eta(size(y)),mu(size(y)),w(size(y)),z(size(y))
      real(dp)::b(size(x,2)),bold(size(x,2)),xtwx(size(x,2),size(x,2)),xtwz(size(x,2)),wx(size(x,1),size(x,2))
      real(dp)::inv(size(x,2),size(x,2)),ll
      n=size(x,1);p=size(x,2);eps=1e-8_dp;if(present(tol))eps=tol;mi=100;if(present(maxiter))mi=maxiter;b=0.0_dp
      b(1)=log(max(1e-8_dp,mean_r(y)))
      do i=1,mi
         bold=b;eta=max(-700.0_dp,min(700.0_dp,matmul(x,b)));mu=exp(eta);w=max(mu,1e-12_dp);z=eta+(y-mu)/w
         do j=1,p;wx(:,j)=w*x(:,j);end do;xtwx=matmul(transpose(x),wx);xtwz=matmul(transpose(x),w*z)
         call solve_linear(xtwx,xtwz,b,info);if(info/=0)exit;if(maxval(abs(b-bold))<eps)exit
      end do
      allocate(res%beta(p),res%fitted(n),res%residuals(n),res%covariance(p,p));res%beta=b;res%iterations=i;res%status=info
      eta=max(-700.0_dp,min(700.0_dp,matmul(x,b)));mu=exp(eta);res%fitted=mu;res%residuals=y-mu;ll=0.0_dp
      do k=1,n
         if(y(k)>0)then;ll=ll+y(k)*log(mu(k))-mu(k)-log_gamma(y(k)+1.0_dp);else;ll=ll-mu(k);end if
      end do
      res%loglik=ll;res%deviance=0.0_dp
      do k=1,n
         if(y(k)>0)then;res%deviance=res%deviance+2*(y(k)*log(y(k)/mu(k))-(y(k)-mu(k)))
         else;res%deviance=res%deviance+2*mu(k);end if
      end do
      w=max(mu,1e-12_dp);do j=1,p;wx(:,j)=w*x(:,j);end do;xtwx=matmul(transpose(x),wx)
      call inverse_matrix(xtwx,inv,info);res%covariance=inv
   end function glm_poisson

   function poisson_only(x,y) result(beta)
      real(dp),intent(in)::x(:,:),y(:);real(dp),allocatable::beta(:);type(regression_result)::r
      r=glm_poisson(x,y);beta=r%beta
   end function poisson_only

   function logistic_only(x,y) result(beta)
      real(dp),intent(in)::x(:,:),y(:);real(dp),allocatable::beta(:);type(regression_result)::r
      r=glm_logistic(x,y);beta=r%beta
   end function logistic_only

   function ridge_regression(x,y,lambda) result(res)
      real(dp),intent(in)::x(:,:),y(:),lambda;type(regression_result)::res
      real(dp)::a(size(x,2),size(x,2)),b(size(x,2)),coef(size(x,2));integer::i,info,n,p
      n=size(x,1);p=size(x,2);a=matmul(transpose(x),x);do i=1,p;a(i,i)=a(i,i)+lambda;end do
      b = matmul(transpose(x),y)
      call solve_linear(a,b,coef,info)
      allocate(res%beta(p), res%fitted(n), res%residuals(n), res%covariance(p,p))
      res%beta=coef;res%fitted=matmul(x,coef);res%residuals=y-res%fitted;res%status=info;res%deviance=sum(res%residuals**2)
   end function ridge_regression

   function ar1_fit(y) result(res)
      real(dp),intent(in)::y(:);type(regression_result)::res;real(dp)::x(size(y)-1,2),yy(size(y)-1)
      x(:,1)=1.0_dp;x(:,2)=y(1:size(y)-1);yy=y(2:size(y));res=lmfit(x,yy,.false.)
   end function ar1_fit

end module rfast_regression
