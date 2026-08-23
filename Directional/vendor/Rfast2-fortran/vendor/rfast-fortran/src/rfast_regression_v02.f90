module rfast_regression_v02
   use rfast_special, only : dp, pi
   use rfast_linalg, only : solve_linear, inverse_matrix
   use rfast_regression, only : regression_result, glm_logistic, glm_poisson, lmfit, weighted_least_squares
   implicit none
   private

   type, public :: multinomial_result
      real(dp), allocatable :: beta(:,:)
      real(dp) :: loglik = -huge(1.0_dp)
      integer :: iterations = 0
      integer :: status = 0
   end type multinomial_result

   type, public :: multivariate_regression_result
      real(dp), allocatable :: beta(:,:)
      real(dp) :: objective = huge(1.0_dp)
      integer :: iterations = 0
      integer :: status = 0
   end type multivariate_regression_result

   public :: gamma_regression, invgauss_regression, quasipoisson_regression
   public :: proportion_regression, multinomial_regression, spatial_median_regression

contains

   subroutine add_intercept(x,xx)
      real(dp),intent(in)::x(:,:)
      real(dp),allocatable,intent(out)::xx(:,:)
      allocate(xx(size(x,1),size(x,2)+1));xx(:,1)=1.0_dp;xx(:,2:)=x
   end subroutine add_intercept

   function gamma_regression(y,x,tol,maxiter) result(res)
      real(dp),intent(in)::y(:),x(:,:)
      real(dp),intent(in),optional::tol
      integer,intent(in),optional::maxiter
      type(regression_result)::res
      real(dp),allocatable::xx(:,:),wx(:,:)
      real(dp),allocatable::b(:),bn(:),score(:),step(:),eta(:),mu(:),con(:)
      real(dp),allocatable::info_mat(:,:),inv(:,:)
      real(dp)::eps,oldobj,obj,phi
      integer::n,p,it,mi,j,info
      if(any(y<=0.0_dp))then;res%status=1;return;end if
      call add_intercept(x,xx);n=size(y);p=size(xx,2);eps=1e-8_dp;if(present(tol))eps=tol;mi=100;if(present(maxiter))mi=maxiter
      allocate(wx(n,p),b(p),bn(p),score(p),step(p),eta(n),mu(n),con(n),info_mat(p,p),inv(p,p))
      b=0.0_dp
      b(1)=log(sum(y)/real(n,dp))
      oldobj=huge(1.0_dp)
      do it=1,mi
         eta=matmul(xx,b);mu=exp(max(-700.0_dp,min(700.0_dp,eta)));con=y/mu
         score=matmul(transpose(xx),con-1.0_dp)
         do j=1,p;wx(:,j)=con*xx(:,j);end do;info_mat=matmul(transpose(xx),wx)
         call solve_linear(info_mat,score,step,info);if(info/=0)then;res%status=info;return;end if
         bn=b+step;obj=sum(y/exp(max(-700.0_dp,min(700.0_dp,matmul(xx,bn))))+matmul(xx,bn))
         if(maxval(abs(bn-b))<=eps.or.abs(oldobj-obj)<=eps)then;b=bn;exit;end if;b=bn;oldobj=obj
      end do
      eta=matmul(xx,b);mu=exp(max(-700.0_dp,min(700.0_dp,eta)));con=y/mu
      phi=sum((con-1.0_dp)**2)/real(max(1,n-p),dp)
      allocate(res%beta(p),res%fitted(n),res%residuals(n),res%covariance(p,p));res%beta=b;res%fitted=mu;res%residuals=y-mu
      res%deviance=-2.0_dp*sum(log(con))+2.0_dp*sum(con)-2.0_dp*real(n,dp);res%dispersion=phi;res%iterations=it
      do j=1,p;wx(:,j)=con*xx(:,j);end do;info_mat=matmul(transpose(xx),wx);call inverse_matrix(info_mat,inv,info)
      if(info==0)res%covariance=phi*inv
   end function gamma_regression

   function invgauss_regression(y,x,tol,maxiter) result(res)
      real(dp),intent(in)::y(:),x(:,:)
      real(dp),intent(in),optional::tol
      integer,intent(in),optional::maxiter
      type(regression_result)::res
      real(dp),allocatable::xx(:,:),b(:),bn(:),f(:),step(:),mu(:),w(:),wx(:,:),h(:,:),inv(:,:)
      real(dp)::eps,lambda,sy
      integer::n,p,it,mi,j,info
      if(any(y<=0.0_dp))then;res%status=1;return;end if
      call add_intercept(x,xx);n=size(y);p=size(xx,2);eps=1e-8_dp;if(present(tol))eps=tol;mi=100;if(present(maxiter))mi=maxiter
      allocate(b(p),bn(p),f(p),step(p),mu(n),w(n),wx(n,p),h(p,p),inv(p,p));b=0.0_dp;b(1)=sum(log(y))/real(n,dp)
      do it=1,mi
         mu=exp(max(-700.0_dp,min(700.0_dp,matmul(xx,b))));f=matmul(transpose(xx),(mu-y)/(mu*mu))
         w=2.0_dp*y/(mu*mu)-1.0_dp/mu;do j=1,p;wx(:,j)=w*xx(:,j);end do;h=matmul(transpose(xx),wx)
         call solve_linear(h,f,step,info);if(info/=0)then;res%status=info;return;end if;bn=b-step
         if(sum(abs(bn-b))<=eps)then;b=bn;exit;end if;b=bn
      end do
      mu=exp(max(-700.0_dp,min(700.0_dp,matmul(xx,b))));sy=sum(1.0_dp/y);lambda=1.0_dp/(sy/real(n,dp)-sum(1.0_dp/mu)/real(n,dp))
      allocate(res%beta(p),res%fitted(n),res%residuals(n),res%covariance(p,p))
      res%beta=b
      res%fitted=mu
      res%residuals=y-mu
      res%iterations=it
      res%loglik=0.5_dp*real(n,dp)*log(lambda/(2.0_dp*pi))-1.5_dp*sum(log(y)) &
          -0.5_dp*lambda*sum(-y/(mu*mu)+sy/real(n,dp));res%deviance=real(n,dp)/lambda
      res%dispersion=sum((y-mu)**2/(mu**3))/real(max(1,n-p),dp)
      w=2.0_dp*y/(mu*mu)-1.0_dp/mu;do j=1,p;wx(:,j)=w*xx(:,j);end do;h=matmul(transpose(xx),wx)
      call inverse_matrix(h,inv,info);if(info==0)res%covariance=inv
   end function invgauss_regression

   function quasipoisson_regression(y,x,tol,maxiter) result(res)
      real(dp),intent(in)::y(:),x(:,:)
      real(dp),intent(in),optional::tol
      integer,intent(in),optional::maxiter
      type(regression_result)::res
      real(dp),allocatable::xx(:,:)
      integer::n,p
      call add_intercept(x,xx);res=glm_poisson(xx,y,tol,maxiter);n=size(y);p=size(xx,2)
      if(res%status==0)then
         res%dispersion=sum((y-res%fitted)**2/max(res%fitted,tiny(1.0_dp)))/real(max(1,n-p),dp)
         res%covariance=res%dispersion*res%covariance
      end if
   end function quasipoisson_regression

   function proportion_regression(y,x,quasi,tol,maxiter) result(res)
      real(dp),intent(in)::y(:),x(:,:)
      logical,intent(in),optional::quasi
      real(dp),intent(in),optional::tol
      integer,intent(in),optional::maxiter
      type(regression_result)::res
      real(dp),allocatable::xx(:,:)
      real(dp)::phi
      integer::n,p
      logical::q
      q=.true.;if(present(quasi))q=quasi;call add_intercept(x,xx);res=glm_logistic(xx,y,tol,maxiter);n=size(y);p=size(xx,2)
      if(res%status==0.and.q)then
         phi=sum((y-res%fitted)**2/max(tiny(1.0_dp),res%fitted*(1.0_dp-res%fitted)))/real(max(1,n-p),dp)
         res%dispersion=phi;res%covariance=phi*res%covariance
      end if
   end function proportion_regression

   function multinomial_regression(y,x,tol,maxiter) result(res)
      integer,intent(in)::y(:)
      real(dp),intent(in)::x(:,:)
      real(dp),intent(in),optional::tol
      integer,intent(in),optional::maxiter
      type(multinomial_result)::res
      real(dp),allocatable::xx(:,:),b(:,:),bn(:,:),prob(:,:),eta(:,:),score(:),step(:),h(:,:)
      real(dp)::eps,den,pr
      integer::n,p,k,d,it,mi,i,j,c,e,ii,jj,info
      n=size(y);k=maxval(y);if(minval(y)<1.or.k<2)then;res%status=1;return;end if
      call add_intercept(x,xx);p=size(xx,2);d=k-1;eps=1e-8_dp;if(present(tol))eps=tol;mi=100;if(present(maxiter))mi=maxiter
      allocate(b(p,d),bn(p,d),prob(n,d),eta(n,d),score(p*d),step(p*d),h(p*d,p*d));b=0.0_dp
      do c=1,d;b(1,c)=log(max(tiny(1.0_dp),real(count(y==c+1),dp))/max(1.0_dp,real(count(y==1),dp)));end do
      do it=1,mi
         eta=matmul(xx,b);prob=0.0_dp
         do i=1,n
            pr=max(0.0_dp,maxval(eta(i,:)));den=exp(-pr)+sum(exp(eta(i,:)-pr));prob(i,:)=exp(eta(i,:)-pr)/den
         end do
         score=0.0_dp;h=0.0_dp
         do c=1,d
            ii=(c-1)*p
            do j=1,p;score(ii+j)=sum(xx(:,j)*(merge(1.0_dp,0.0_dp,y==c+1)-prob(:,c)));end do
            do e=1,d
               jj=(e-1)*p
               do j=1,p
                  h(ii+j,jj+1:jj+p)=h(ii+j,jj+1:jj+p) &
                      +matmul((prob(:,c)*(merge(1.0_dp,0.0_dp,c==e)-prob(:,e))*xx(:,j)),xx)
               end do
            end do
         end do
         call solve_linear(h,score,step,info);if(info/=0)then;res%status=info;return;end if
         bn=b+reshape(step,[p,d]);if(maxval(abs(bn-b))<=eps)then;b=bn;exit;end if;b=bn
      end do
      eta=matmul(xx,b);res%loglik=0.0_dp
      do i=1,n
         pr=max(0.0_dp,maxval(eta(i,:)));den=exp(-pr)+sum(exp(eta(i,:)-pr))
         if(y(i)==1)then;res%loglik=res%loglik-pr-log(den);else;res%loglik=res%loglik+eta(i,y(i)-1)-pr-log(den);end if
      end do
      allocate(res%beta(p,d));res%beta=b;res%iterations=it
   end function multinomial_regression

   function spatial_median_regression(y,x,tol,maxiter) result(res)
      real(dp),intent(in)::y(:,:),x(:,:)
      real(dp),intent(in),optional::tol
      integer,intent(in),optional::maxiter
      type(multivariate_regression_result)::res
      real(dp),allocatable::xx(:,:),b(:,:),bn(:,:),fit(:,:),e(:,:),w(:),a(:,:),rhs(:),sol(:)
      real(dp)::eps,obj,oldobj
      integer::n,p,d,it,mi,i,j,info
      call add_intercept(x,xx);n=size(y,1);d=size(y,2);p=size(xx,2);if(size(xx,1)/=n)then;res%status=1;return;end if
      eps=1e-8_dp;if(present(tol))eps=tol;mi=200;if(present(maxiter))mi=maxiter
      allocate(b(p,d),bn(p,d),fit(n,d),e(n,d),w(n),a(p,p),rhs(p),sol(p));a=matmul(transpose(xx),xx)
      do j=1,d
         rhs=matmul(transpose(xx),y(:,j))
         call solve_linear(a,rhs,sol,info)
         if(info/=0)then;res%status=info;return;end if
         b(:,j)=sol
      end do
      oldobj=huge(1.0_dp)
      do it=1,mi
         fit=matmul(xx,b);e=y-fit;do i=1,n;w(i)=1.0_dp/max(1e-12_dp,sqrt(sum(e(i,:)**2)));end do
         a=matmul(transpose(xx),spread(w,2,p)*xx)
         do j=1,d;rhs=matmul(transpose(xx),w*y(:,j));call solve_linear(a,rhs,sol,info);if(info/=0)exit;bn(:,j)=sol;end do
         if(info/=0)then;res%status=info;return;end if
         fit=matmul(xx,bn);e=y-fit;obj=sum([(sqrt(sum(e(i,:)**2)),i=1,n)])
         if(maxval(abs(bn-b))<=eps.or.abs(oldobj-obj)<=eps)then;b=bn;oldobj=obj;exit;end if;b=bn;oldobj=obj
      end do
      allocate(res%beta(p,d));res%beta=b;res%objective=oldobj;res%iterations=it
   end function spatial_median_regression

end module rfast_regression_v02
