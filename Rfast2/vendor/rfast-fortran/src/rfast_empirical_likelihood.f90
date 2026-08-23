module rfast_empirical_likelihood
   use rfast_special, only : dp, reg_gamma_q
   use rfast_linalg, only : solve_linear
   implicit none
   private

   type, public :: el_result
      real(dp) :: lambda = 0.0_dp
      real(dp) :: statistic = huge(1.0_dp)
      real(dp) :: pvalue = 0.0_dp
      real(dp), allocatable :: weights(:)
      integer :: iterations = 0
      integer :: status = 0
   end type el_result

   type, public :: el_two_sample_result
      real(dp) :: mean = 0.0_dp
      real(dp) :: lambda = 0.0_dp
      real(dp) :: statistic = huge(1.0_dp)
      real(dp) :: pvalue = 0.0_dp
      real(dp), allocatable :: weights1(:), weights2(:)
      integer :: iterations = 0
      integer :: status = 0
   end type el_two_sample_result

   type, public :: mv_el_result
      real(dp), allocatable :: lambda(:), weights1(:), weights2(:)
      real(dp) :: statistic = huge(1.0_dp)
      real(dp) :: pvalue = 0.0_dp
      integer :: df = 0
      integer :: iterations = 0
      integer :: status = 0
   end type mv_el_result

   public :: el_test1, el_test2, eel_test1, eel_test2, mv_eeltest1, mv_eeltest2

contains

   pure real(dp) function chisq_sf(x,df) result(p)
      real(dp),intent(in)::x
      integer,intent(in)::df
      p=reg_gamma_q(0.5_dp*real(df,dp),0.5_dp*max(0.0_dp,x))
   end function chisq_sf

   function el_test1(x,mu,tol) result(res)
      real(dp),intent(in)::x(:),mu
      real(dp),intent(in),optional::tol
      type(el_result)::res
      real(dp)::y(size(x)),lo,hi,mid,f,eps,delta
      integer::it,n
      n=size(x);eps=1e-9_dp;if(present(tol))eps=tol;y=x-mu
      if(maxval(y)<=0.0_dp.or.minval(y)>=0.0_dp)then;res%status=1;return;end if
      lo=-1.0_dp/maxval(y);hi=-1.0_dp/minval(y)
      delta=100.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(lo),abs(hi));lo=lo+delta;hi=hi-delta
      do it=1,200
         mid=0.5_dp*(lo+hi);f=sum(y/(1.0_dp+mid*y))
         if(abs(f)<=eps.or.abs(hi-lo)<=eps*max(1.0_dp,abs(mid)))exit
         if(f>0.0_dp)then;lo=mid;else;hi=mid;end if
      end do
      res%lambda=mid;allocate(res%weights(n));res%weights=1.0_dp/(1.0_dp+mid*y);res%weights=res%weights/sum(res%weights)
      res%statistic=2.0_dp*sum(log(1.0_dp+mid*y));res%pvalue=chisq_sf(res%statistic,1);res%iterations=it
   end function el_test1

   real(dp) function el_profile_stat(mu,x,y,tol) result(v)
      real(dp),intent(in)::mu,x(:),y(:),tol
      type(el_result)::a,b
      a=el_test1(x,mu,tol);b=el_test1(y,mu,tol)
      if(a%status/=0.or.b%status/=0)then;v=huge(1.0_dp);else;v=a%statistic+b%statistic;end if
   end function el_profile_stat

   function el_test2(x,y,tol) result(res)
      real(dp),intent(in)::x(:),y(:)
      real(dp),intent(in),optional::tol
      type(el_two_sample_result)::res
      real(dp)::eps,a,b,c,d,fc,fd,phi,mu
      integer::it
      type(el_result)::r1,r2
      eps=1e-8_dp;if(present(tol))eps=tol;a=min(sum(x)/size(x),sum(y)/size(y))-1.0_dp
      b=max(sum(x)/size(x),sum(y)/size(y))+1.0_dp;phi=(sqrt(5.0_dp)-1.0_dp)/2.0_dp
      c=b-phi*(b-a);d=a+phi*(b-a);fc=el_profile_stat(c,x,y,eps);fd=el_profile_stat(d,x,y,eps)
      do it=1,200
         if(abs(b-a)<=eps*max(1.0_dp,abs(a)+abs(b)))exit
         if(fc<fd)then;b=d;d=c;fd=fc;c=b-phi*(b-a);fc=el_profile_stat(c,x,y,eps)
         else;a=c;c=d;fc=fd;d=a+phi*(b-a);fd=el_profile_stat(d,x,y,eps);end if
      end do
      mu=0.5_dp*(a+b);r1=el_test1(x,mu,eps);r2=el_test1(y,mu,eps)
      res%mean=mu;res%statistic=r1%statistic+r2%statistic;res%pvalue=chisq_sf(res%statistic,1);res%iterations=it
      if(r1%status/=0.or.r2%status/=0)then;res%status=1;return;end if
      allocate(res%weights1(size(x)),res%weights2(size(y)));res%weights1=r1%weights;res%weights2=r2%weights
   end function el_test2

   function eel_test1(x,mu,tol,maxiter) result(res)
      real(dp),intent(in)::x(:),mu
      real(dp),intent(in),optional::tol
      integer,intent(in),optional::maxiter
      type(el_result)::res
      real(dp)::lam,lam2,eps,mx,s0,s1,s2,f,der,z(size(x)),w(size(x))
      integer::it,mi,n
      n=size(x);eps=1e-10_dp;if(present(tol))eps=tol;mi=200;if(present(maxiter))mi=maxiter;lam=0.0_dp;lam2=lam
      do it=1,mi
         z=lam*x;mx=maxval(z);w=exp(z-mx);s0=sum(w);s1=sum(x*w);s2=sum(x*x*w)
         f=s1/s0-mu;der=s2/s0-(s1/s0)**2
         if(der<=tiny(1.0_dp))then;res%status=2;return;end if
         lam2=lam-f/der;if(abs(lam2-lam)<=eps*max(1.0_dp,abs(lam)))exit;lam=lam2
      end do
      z=lam2*x;mx=maxval(z);w=exp(z-mx);w=w/sum(w);allocate(res%weights(n));res%weights=w
      res%lambda=lam2;res%statistic=-2.0_dp*sum(log(real(n,dp)*w));res%pvalue=chisq_sf(res%statistic,1);res%iterations=it
   end function eel_test1

   function eel_test2(x,y,tol,maxiter) result(res)
      real(dp),intent(in)::x(:),y(:)
      real(dp),intent(in),optional::tol
      integer,intent(in),optional::maxiter
      type(el_two_sample_result)::res
      real(dp)::lam,lam2,eps,mx,my,sx0,sx1,sx2,sy0,sy1,sy2,f,der
      real(dp)::wx(size(x)),wy(size(y)),zx(size(x)),zy(size(y))
      integer::it,mi,n1,n2
      n1=size(x);n2=size(y);eps=1e-10_dp;if(present(tol))eps=tol;mi=200;if(present(maxiter))mi=maxiter;lam=0.0_dp;lam2=lam
      do it=1,mi
         zx=lam*x;mx=maxval(zx);wx=exp(zx-mx);sx0=sum(wx);sx1=sum(x*wx);sx2=sum(x*x*wx)
         zy=-lam*y;my=maxval(zy);wy=exp(zy-my);sy0=sum(wy);sy1=sum(y*wy);sy2=sum(y*y*wy)
         f=sx1/sx0-sy1/sy0;der=sx2/sx0-(sx1/sx0)**2+sy2/sy0-(sy1/sy0)**2
         if(der<=tiny(1.0_dp))then;res%status=2;return;end if
         lam2=lam-f/der;if(abs(lam2-lam)<=eps*max(1.0_dp,abs(lam)))exit;lam=lam2
      end do
      zx=lam2*x;mx=maxval(zx);wx=exp(zx-mx);wx=wx/sum(wx);zy=-lam2*y;my=maxval(zy);wy=exp(zy-my);wy=wy/sum(wy)
      allocate(res%weights1(n1),res%weights2(n2));res%weights1=wx;res%weights2=wy;res%lambda=lam2
      res%statistic=-2.0_dp*sum(log(real(n1,dp)*wx))-2.0_dp*sum(log(real(n2,dp)*wy))
      res%pvalue=chisq_sf(res%statistic,1);res%iterations=it
   end function eel_test2

   function mv_eeltest1(x,mu,tol,maxiter) result(res)
      real(dp),intent(in)::x(:,:),mu(:)
      real(dp),intent(in),optional::tol
      integer,intent(in),optional::maxiter
      type(mv_el_result)::res
      integer::n,d,it,mi,i,j,k,info
      real(dp)::eps,mx,s0
      real(dp),allocatable::lam(:),ln(:),w(:),mean(:),f(:),der(:,:),step(:),z(:)
      n=size(x,1);d=size(x,2);eps=1e-8_dp;if(present(tol))eps=tol;mi=200;if(present(maxiter))mi=maxiter
      if(size(mu)/=d)then;res%status=1;return;end if
      allocate(lam(d),ln(d),w(n),mean(d),f(d),der(d,d),step(d),z(n));lam=0.0_dp
      do it=1,mi
         z=matmul(x,lam);mx=maxval(z);w=exp(z-mx);s0=sum(w);w=w/s0;mean=matmul(transpose(x),w);f=mean-mu;der=0.0_dp
         do i=1,n
            do j=1,d
               do k=1,d;der(j,k)=der(j,k)+w(i)*(x(i,j)-mean(j))*(x(i,k)-mean(k));end do
            end do
         end do
         call solve_linear(der,f,step,info);if(info/=0)then;res%status=info;return;end if
         ln=lam-step;if(sum(abs(ln-lam))<=eps*max(1.0_dp,sum(abs(lam))))exit;lam=ln
      end do
      z=matmul(x,ln);mx=maxval(z);w=exp(z-mx);w=w/sum(w);allocate(res%lambda(d),res%weights1(n));res%lambda=ln;res%weights1=w
      res%statistic=-2.0_dp*sum(log(real(n,dp)*w));res%df=d;res%pvalue=chisq_sf(res%statistic,d);res%iterations=it
   end function mv_eeltest1

   function mv_eeltest2(x,y,tol,maxiter) result(res)
      real(dp),intent(in)::x(:,:),y(:,:)
      real(dp),intent(in),optional::tol
      integer,intent(in),optional::maxiter
      type(mv_el_result)::res
      integer::n1,n2,d,it,mi,i,j,k,info
      real(dp)::eps,mx,my
      real(dp),allocatable::lam(:),ln(:),wx(:),wy(:),meanx(:),meany(:),f(:),der(:,:),step(:),zx(:),zy(:)
      n1=size(x,1);n2=size(y,1);d=size(x,2);eps=1e-8_dp;if(present(tol))eps=tol;mi=200;if(present(maxiter))mi=maxiter
      if(size(y,2)/=d)then;res%status=1;return;end if
      allocate(lam(d),ln(d),wx(n1),wy(n2),meanx(d),meany(d),f(d),der(d,d),step(d),zx(n1),zy(n2));lam=0.0_dp
      do it=1,mi
         zx=matmul(x,lam);mx=maxval(zx);wx=exp(zx-mx);wx=wx/sum(wx)
         zy=-matmul(y,lam);my=maxval(zy);wy=exp(zy-my);wy=wy/sum(wy)
         meanx=matmul(transpose(x),wx);meany=matmul(transpose(y),wy);f=meanx-meany;der=0.0_dp
         do i=1,n1;do j=1,d;do k=1,d;der(j,k)=der(j,k)+wx(i)*(x(i,j)-meanx(j))*(x(i,k)-meanx(k));end do;end do;end do
         do i=1,n2;do j=1,d;do k=1,d;der(j,k)=der(j,k)+wy(i)*(y(i,j)-meany(j))*(y(i,k)-meany(k));end do;end do;end do
         call solve_linear(der,f,step,info);if(info/=0)then;res%status=info;return;end if
         ln=lam-step;if(sum(abs(ln-lam))<=eps*max(1.0_dp,sum(abs(lam))))exit;lam=ln
      end do
      zx=matmul(x,ln);mx=maxval(zx);wx=exp(zx-mx);wx=wx/sum(wx);zy=-matmul(y,ln);my=maxval(zy);wy=exp(zy-my);wy=wy/sum(wy)
      allocate(res%lambda(d),res%weights1(n1),res%weights2(n2));res%lambda=ln;res%weights1=wx;res%weights2=wy
      res%statistic=-2.0_dp*sum(log(real(n1,dp)*wx))-2.0_dp*sum(log(real(n2,dp)*wy));res%df=d
      res%pvalue=chisq_sf(res%statistic,d);res%iterations=it
   end function mv_eeltest2

end module rfast_empirical_likelihood
