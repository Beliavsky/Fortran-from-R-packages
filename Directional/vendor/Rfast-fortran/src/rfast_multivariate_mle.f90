module rfast_multivariate_mle
   use rfast_special, only : dp, pi, digamma_r, trigamma_r
   use rfast_arrays, only : colmeans
   use rfast_linalg, only : covariance_matrix, determinant_matrix, logdet_spd, mahalanobis_sq, solve_linear
   implicit none
   private

   type, public :: multivariate_mle_result
      real(dp), allocatable :: location(:)
      real(dp), allocatable :: scatter(:,:)
      real(dp), allocatable :: mean_original(:)
      real(dp), allocatable :: covariance_original(:,:)
      real(dp), allocatable :: param(:)
      real(dp) :: loglik = -huge(1.0_dp)
      integer :: iterations = 0
      integer :: status = 0
   end type multivariate_mle_result

   public :: mvnorm_mle, mvlognormal_mle, mvt_mle, dirichlet_mle, inverse_dirichlet_mle

contains

   function mvnorm_mle(x) result(res)
      real(dp), intent(in) :: x(:,:)
      type(multivariate_mle_result) :: res
      integer :: n,d,info
      real(dp) :: ld
      n=size(x,1);d=size(x,2)
      allocate(res%location(d),res%scatter(d,d))
      res%location=colmeans(x)
      res%scatter=covariance_matrix(x,.true.)
      ld=logdet_spd(res%scatter,info)
      if(info/=0)then;res%status=info;return;end if
      res%loglik=-0.5_dp*real(n,dp)*(real(d,dp)*log(2.0_dp*pi)+ld+real(d,dp))
      res%iterations=1
   end function mvnorm_mle

   function mvlognormal_mle(x) result(res)
      real(dp), intent(in) :: x(:,:)
      type(multivariate_mle_result) :: res
      real(dp), allocatable :: y(:,:)
      real(dp) :: si(size(x,2))
      integer :: n,d,i,j,info
      real(dp) :: ld
      n=size(x,1);d=size(x,2)
      if(any(x<=0.0_dp))then;res%status=1;return;end if
      allocate(y(n,d),res%location(d),res%scatter(d,d),res%mean_original(d),res%covariance_original(d,d))
      y=log(x);res%location=colmeans(y);res%scatter=covariance_matrix(y,.true.)
      ld=logdet_spd(res%scatter,info);if(info/=0)then;res%status=info;return;end if
      res%loglik=-0.5_dp*real(n,dp)*(real(d,dp)*log(2.0_dp*pi)+ld+real(d,dp))-sum(y)
      do i=1,d;si(i)=res%scatter(i,i);res%mean_original(i)=exp(res%location(i)+0.5_dp*si(i));end do
      do j=1,d
         do i=1,d
            res%covariance_original(i,j)=res%mean_original(i)*res%mean_original(j)*(exp(res%scatter(i,j))-1.0_dp)
         end do
      end do
      res%iterations=1
   end function mvlognormal_mle

   function mvt_mle(x,v,tol,maxiter) result(res)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(in), optional :: v,tol
      integer, intent(in), optional :: maxiter
      type(multivariate_mle_result) :: res
      integer :: n,p,it,mi,info,i
      real(dp) :: nu,eps,con,el1,el2,sumw,ld
      real(dp), allocatable :: w(:),dis(:),z(:,:),m(:),r(:,:)
      n=size(x,1);p=size(x,2);nu=5.0_dp;if(present(v))nu=v;eps=1e-7_dp;if(present(tol))eps=tol
      mi=200;if(present(maxiter))mi=maxiter
      if(nu<=0.0_dp.or.n<=p)then;res%status=1;return;end if
      allocate(w(n),dis(n),z(n,p),m(p),r(p,p))
      m=colmeans(x);r=covariance_matrix(x,.false.);if(nu/=1.0_dp)r=abs(nu-1.0_dp)/nu*r
      con=real(n,dp)*log_gamma(0.5_dp*(nu+real(p,dp)))-real(n,dp)*log_gamma(0.5_dp*nu) &
          -0.5_dp*real(n*p,dp)*log(pi*nu)
      el1=-huge(1.0_dp);el2=-huge(1.0_dp)
      do it=1,mi
         dis=mahalanobis_sq(x,m,r);w=(nu+real(p,dp))/(nu+dis);sumw=sum(w)
         do i=1,n;z(i,:)=sqrt(w(i))*(x(i,:)-m);end do
         r=matmul(transpose(z),z)/sumw
         m=matmul(transpose(x),w)/sumw
         dis=mahalanobis_sq(x,m,r);ld=logdet_spd(r,info)
         if(info/=0)then;res%status=info;return;end if
         el2=-real(n,dp)*ld-(nu+real(p,dp))*sum(log(1.0_dp+dis/nu))
         if(it>1.and.el2-el1<=eps)exit
         el1=el2
      end do
      allocate(res%location(p),res%scatter(p,p));res%location=m;res%scatter=r
      res%loglik=0.5_dp*el2+con;res%iterations=it
   end function mvt_mle

   function dirichlet_mle(x,tol,maxiter) result(res)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxiter
      type(multivariate_mle_result) :: res
      integer :: n,p,it,mi,j,info
      real(dp) :: eps,sa,phi
      real(dp), allocatable :: a(:),an(:),g(:),h(:,:),step(:),sl(:),m(:),v(:)
      n=size(x,1);p=size(x,2);eps=1e-8_dp;if(present(tol))eps=tol;mi=200;if(present(maxiter))mi=maxiter
      if(any(x<=0.0_dp).or.any(x>=1.0_dp))then;res%status=1;return;end if
      allocate(a(p),an(p),g(p),h(p,p),step(p),sl(p),m(p),v(p));m=colmeans(x)
      do j=1,p
         v(j)=sum((x(:,j)-m(j))**2)/real(n,dp)
      end do
      phi=huge(1.0_dp)
      do j=1,p
         if(v(j)>0.0_dp.and.m(j)>0.0_dp.and.m(j)<1.0_dp)phi=min(phi,m(j)*(1.0_dp-m(j))/v(j)-1.0_dp)
      end do
      if(phi<=0.0_dp.or.phi>=huge(1.0_dp)/2.0_dp)phi=real(p,dp)
      a=max(1e-3_dp,m*phi);sl=0.0_dp
      do j=1,p;sl(j)=sum(log(x(:,j)));end do
      do it=1,mi
         sa=sum(a);g=real(n,dp)*(digamma_r(sa)-digamma_r(a))+sl
         h=real(n,dp)*trigamma_r(sa)
         do j=1,p;h(j,j)=h(j,j)-real(n,dp)*trigamma_r(a(j));end do
         call solve_linear(h,g,step,info);if(info/=0)then;res%status=info;return;end if
         an=a-step
         do while(any(an<=0.0_dp));step=0.5_dp*step;an=a-step;end do
         if(sum(abs(an-a))<=eps*max(1.0_dp,sum(abs(a))))exit
         a=an
      end do
      a=an;allocate(res%param(p));res%param=a;res%iterations=it
      res%loglik=real(n,dp)*log_gamma(sum(a))-real(n,dp)*sum(log_gamma(a))
      do j=1,p;res%loglik=res%loglik+(a(j)-1.0_dp)*sl(j);end do
   end function dirichlet_mle

   function inverse_dirichlet_mle(x,tol,maxiter) result(res)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxiter
      type(multivariate_mle_result) :: res
      integer :: n,p,d,it,mi,j,info
      real(dp) :: eps,phi,mv,vv,aD,lik1,lik2,sx2,sumlogs
      real(dp),allocatable :: a(:),an(:),g(:),h(:,:),step(:),com(:),clog(:),cm(:),cv(:),rs(:)
      n=size(x,1);p=size(x,2);d=p+1;eps=1e-8_dp;if(present(tol))eps=tol;mi=200;if(present(maxiter))mi=maxiter
      if(any(x<=0.0_dp))then;res%status=1;return;end if
      allocate(a(d),an(d),g(d),h(d,d),step(d),com(d),clog(p),cm(p),cv(p),rs(n))
      rs=sum(x,dim=2);sx2=sum(log(1.0_dp+rs));sumlogs=0.0_dp
      do j=1,p
         clog(j)=sum(log(x(:,j)));cm(j)=sum(x(:,j))/real(n,dp);cv(j)=sum((x(:,j)-cm(j))**2)/real(n,dp)
         sumlogs=sumlogs+clog(j)
      end do
      mv=sum(cm)/real(p,dp);vv=sum(cv)/real(p,dp)
      aD=0.5_dp*(mv*mv+mv)/max(vv,1e-8_dp)+1.0_dp
      a(1:p)=0.5_dp*abs(cm*(aD-1.0_dp));a(d)=0.5_dp*abs(aD);a=max(a,1e-3_dp)
      com(1:p)=clog-sx2;com(d)=-sx2;lik1=-huge(1.0_dp);lik2=-huge(1.0_dp)
      do it=1,mi
         phi=sum(a);g=real(n,dp)*(digamma_r(phi)-digamma_r(a))+com
         h=real(n,dp)*trigamma_r(phi)
         do j=1,d;h(j,j)=h(j,j)-real(n,dp)*trigamma_r(a(j));end do
         call solve_linear(h,g,step,info);if(info/=0)then;res%status=info;return;end if
         an=a-step;do while(any(an<=0.0_dp));step=0.5_dp*step;an=a-step;end do
         phi=sum(an);lik2=real(n,dp)*log_gamma(phi)-real(n,dp)*sum(log_gamma(an))-phi*sx2
         do j=1,p;lik2=lik2+an(j)*clog(j);end do
         if(it>1.and.abs(lik2-lik1)<=eps)exit
         a=an;lik1=lik2
      end do
      allocate(res%param(d));res%param=an;res%iterations=it;res%loglik=lik2-sumlogs
   end function inverse_dirichlet_mle

end module rfast_multivariate_mle
