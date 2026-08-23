module aep_fit
   use aep_special, only: dp, golden_min, bisect_root, digamma_aep
   use aep_distribution, only: daep, paep
   implicit none
   private
   public :: aep_fit_result, fitaep, aep_reg_result, regaep

   type :: aep_fit_result
      real(dp) :: alpha = 1.0_dp
      real(dp) :: sigma = 1.0_dp
      real(dp) :: mu = 0.0_dp
      real(dp) :: epsilon = 0.0_dp
      real(dp) :: log_likelihood = 0.0_dp
      real(dp) :: aic = 0.0_dp
      real(dp) :: caic = 0.0_dp
      real(dp) :: bic = 0.0_dp
      real(dp) :: hqic = 0.0_dp
      real(dp) :: ad = 0.0_dp
      real(dp) :: cvm = 0.0_dp
      real(dp) :: ks = 0.0_dp
      real(dp) :: inv_ofim(4,4) = 0.0_dp
      integer :: iterations = 0
      logical :: converged = .false.
   end type

   type :: aep_reg_result
      real(dp), allocatable :: beta(:)
      real(dp) :: alpha = 1.0_dp
      real(dp) :: sigma = 1.0_dp
      real(dp) :: epsilon = 0.0_dp
      real(dp) :: r2 = 0.0_dp
      real(dp) :: adjusted_r2 = 0.0_dp
      real(dp) :: f_stat = 0.0_dp
      real(dp), allocatable :: inv_ofim(:,:)
      real(dp), allocatable :: residuals(:)
      integer :: iterations = 0
      logical :: converged = .false.
   end type
contains
   subroutine fitaep(x, fit, starts, max_iter, tol)
      real(dp), intent(in) :: x(:)
      type(aep_fit_result), intent(out) :: fit
      real(dp), intent(in), optional :: starts(4)
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: tol
      integer :: n, it, mit, i
      real(dp) :: alpha,sigma,mu,epsi,old(4),newp(4),conv,ttol,varx,m4
      real(dp), allocatable :: e(:), xs(:), cdf(:), pdf(:), d(:,:)
      logical :: ok
      n=size(x); if(n<4) error stop "fitaep: need at least four observations"
      mit=6000; if(present(max_iter))mit=max_iter
      ttol=1.0e-4_dp; if(present(tol))ttol=tol
      if(present(starts))then
         alpha=starts(1);sigma=starts(2);mu=starts(3);epsi=starts(4)
      else
         mu=sample_median(x); varx=sample_var(x)
         m4=sum((x-sum(x)/real(n,dp))**4)/real(n,dp)
         block
            real(dp) :: fa,fb
            fa=kurt_eq(0.05_dp);fb=kurt_eq(2.0_dp)
            if(fa*fb<0.0_dp)then
               alpha=bisect_root(kurt_eq,0.05_dp,2.0_dp,1.0e-9_dp,ok)
            else
               alpha=1.0_dp
            end if
         end block
         epsi=1.0_dp-2.0_dp*real(count(x<mu),dp)/real(n,dp)
         epsi=max(-0.999_dp,min(0.999_dp,epsi))
         sigma=sqrt(max(tiny(1.0_dp),varx*gamma(1.0_dp/alpha)/gamma(3.0_dp/alpha)))
      end if
      allocate(e(n))
      fit%converged=.false.
      do it=1,mit
         old=[alpha,sigma,mu,epsi]
         do i=1,n
            if(abs(x(i)-mu)<=1.0e-8_dp)then
               e(i)=max(1.0e-12_dp,abs(mu))
            else
               e(i)=alpha/2.0_dp*(abs(x(i)-mu)/sigma)**(alpha-2.0_dp) * &
                    abs(1.0_dp+sign(1.0_dp,x(i)-mu)*epsi)**(2.0_dp-alpha)
            end if
         end do
         mu=sum(x*e/(1.0_dp+sgnvec(x-mu)*epsi)**2)/sum(e/(1.0_dp+sgnvec(x-mu)*epsi)**2)
         sigma=sqrt(max(tiny(1.0_dp),2.0_dp/real(n,dp)*sum((x-mu)**2*e/(1.0_dp+sgnvec(x-mu)*epsi)**2)))
         epsi=golden_min(eps_obj,-0.999_dp,0.999_dp,1.0e-9_dp)
         alpha=golden_min(alpha_obj,0.01_dp,2.0_dp,1.0e-9_dp)
         newp=[alpha,sigma,mu,epsi]
         conv=sum(abs(newp-old))
         if(conv<ttol)then;fit%converged=.true.;exit;end if
      end do
      fit%alpha=alpha;fit%sigma=sigma;fit%mu=mu;fit%epsilon=epsi;fit%iterations=it
      allocate(xs(n),cdf(n),pdf(n),d(n,4));xs=x;call sort_real(xs)
      cdf=paep(xs,alpha,sigma,mu,epsi);pdf=daep(xs,alpha,sigma,mu,epsi)
      fit%log_likelihood=sum(log(max(pdf,tiny(1.0_dp))))
      fit%ks=0.0_dp;fit%cvm=1.0_dp/(12.0_dp*real(n,dp));fit%ad=-real(n,dp)
      do i=1,n
         cdf(i)=min(1.0_dp-1.0e-15_dp,max(1.0e-15_dp,cdf(i)))
         fit%ks=max(fit%ks,max(real(i,dp)/n-cdf(i),cdf(i)-real(i-1,dp)/n))
         fit%cvm=fit%cvm+(cdf(i)-(2.0_dp*real(i,dp)-1.0_dp)/(2.0_dp*n))**2
         fit%ad=fit%ad-( (2.0_dp*i-1.0_dp)*log(cdf(i)) + &
              (2.0_dp*n+1.0_dp-2.0_dp*i)*log(1.0_dp-cdf(i)) )/real(n,dp)
      end do
      fit%aic=-2.0_dp*fit%log_likelihood+8.0_dp
      if(n>5)then
         fit%caic=fit%aic+40.0_dp/real(n-5,dp)
      else
         fit%caic=huge(1.0_dp)
      end if
      fit%bic=-2.0_dp*fit%log_likelihood+4.0_dp*log(real(n,dp))
      fit%hqic=-2.0_dp*fit%log_likelihood+8.0_dp*log(log(real(n,dp)))
      call score_matrix(x,alpha,sigma,mu,epsi,d)
      call inv_xtx(d,fit%inv_ofim)
   contains
      real(dp) function kurt_eq(a) result(v)
         real(dp),intent(in)::a
         v=-m4/max(tiny(1.0_dp),((real(n-1,dp)/n)*varx)**2) + &
           gamma(5.0_dp/a)*gamma(1.0_dp/a)/gamma(3.0_dp/a)**2
      end function
      real(dp) function eps_obj(z) result(v)
         real(dp),intent(in)::z
         v=sum((x-mu)**2*e/(sigma*sigma*(1.0_dp+sgnvec(x-mu)*z)**2))
      end function
      real(dp) function alpha_obj(a) result(v)
         real(dp),intent(in)::a
         v=real(n,dp)*log_gamma(1.0_dp+1.0_dp/a)+sum((abs(x-mu)/(sigma*(1.0_dp+sgnvec(x-mu)*epsi)))**a)
      end function
   end subroutine fitaep

   subroutine regaep(y, x, fit, max_iter, tol)
      real(dp),intent(in)::y(:),x(:,:)
      type(aep_reg_result),intent(out)::fit
      integer,intent(in),optional::max_iter
      real(dp),intent(in),optional::tol
      integer::n,q,p,it,mit,i,j
      real(dp)::alpha,sigma,epsi,mu,ttol,conv,old_alpha,old_sigma,old_eps,st,se
      real(dp),allocatable::z(:,:),beta(:),old_beta(:),r(:),e(:),wgt(:),a(:,:),rhs(:),d(:,:)
      n=size(y);q=size(x,2);p=q+1
      if(size(x,1)/=n.or.n<=p+3)error stop "regaep: invalid dimensions"
      mit=2000;if(present(max_iter))mit=max_iter;ttol=1.0e-4_dp;if(present(tol))ttol=tol
      allocate(z(n,p),beta(p),old_beta(p),r(n),e(n),wgt(n),a(p,p),rhs(p),d(n,p+3))
      z(:,1)=1.0_dp;if(q>0)z(:,2:)=x
      call least_squares(z,y,beta)
      r=y-matmul(z,beta)
      alpha=1.0_dp;mu=sample_median(r);epsi=1.0_dp-2.0_dp*real(count(r<mu),dp)/n
      epsi=max(-0.999_dp,min(0.999_dp,epsi));sigma=sqrt(max(tiny(1.0_dp),sample_var(r)*gamma(1.0_dp/alpha)/gamma(3.0_dp/alpha)))
      fit%converged=.false.
      do it=1,mit
         old_beta=beta;old_alpha=alpha;old_sigma=sigma;old_eps=epsi
         r=y-matmul(z,beta);mu=0.0_dp
         do i=1,n
            if(abs(r(i)-mu)<=1.0e-8_dp)then;e(i)=1.0e-12_dp
            else;e(i)=alpha/2.0_dp*(abs(r(i)-mu)/sigma)**(alpha-2.0_dp)* &
                 abs(1.0_dp+sign(1.0_dp,r(i)-mu)*epsi)**(2.0_dp-alpha);end if
         end do
         sigma=sqrt(max(tiny(1.0_dp),2.0_dp/n*sum((r-mu)**2*e/(1.0_dp+sgnvec(r-mu)*epsi)**2)))
         epsi=golden_min(reps,-0.999_dp,0.999_dp,1.0e-9_dp)
         alpha=golden_min(ralpha,0.01_dp,2.0_dp,1.0e-9_dp)
         wgt=e/(1.0_dp+sgnvec(r-mu)*epsi)**2
         a=0.0_dp;rhs=0.0_dp
         do i=1,n
            do j=1,p
               rhs(j)=rhs(j)+z(i,j)*(y(i)-mu)*wgt(i)
            end do
            a=a+wgt(i)*outer(z(i,:),z(i,:))
         end do
         call solve_linear(a,rhs,beta)
         conv=sum(abs(beta-old_beta))+abs(alpha-old_alpha)+abs(sigma-old_sigma)+abs(epsi-old_eps)
         if(conv<ttol)then;fit%converged=.true.;exit;end if
      end do
      r=y-matmul(z,beta);se=sum(r*r);st=sum((y-sum(y)/n)**2)
      allocate(fit%beta(p),fit%residuals(n),fit%inv_ofim(p+3,p+3));fit%beta=beta;fit%residuals=r
      fit%alpha=alpha;fit%sigma=sigma;fit%epsilon=epsi;fit%iterations=it
      fit%r2=1.0_dp-se/st;fit%adjusted_r2=1.0_dp-real(n-1,dp)/real(n-p,dp)*se/st
      if(p>1)fit%f_stat=((st-se)/real(p-1,dp))/(se/real(n-p,dp))
      call regression_score(z,r,alpha,sigma,epsi,d);call inv_xtx(d,fit%inv_ofim)
   contains
      real(dp) function reps(zeta) result(v)
         real(dp),intent(in)::zeta
         v=sum((r-mu)**2*e/(sigma*sigma*(1.0_dp+sgnvec(r-mu)*zeta)**2))
      end function
      real(dp) function ralpha(aa) result(v)
         real(dp),intent(in)::aa
         v=n*log_gamma(1.0_dp+1.0_dp/aa)+sum((abs(r-mu)/(sigma*(1.0_dp+sgnvec(r-mu)*epsi)))**aa)
      end function
   end subroutine regaep

   subroutine score_matrix(x,alpha,sigma,mu,epsi,d)
      real(dp),intent(in)::x(:),alpha,sigma,mu,epsi;real(dp),intent(out)::d(:,:)
      integer::i;real(dp)::s,z,lz
      do i=1,size(x)
         s=sign(1.0_dp,x(i)-mu);if(x(i)==mu)s=0.0_dp
         z=abs(x(i)-mu)/(sigma*(1.0_dp+s*epsi));lz=log(max(z,tiny(1.0_dp)))
         d(i,1)=digamma_aep(1.0_dp+1.0_dp/alpha)/alpha**2-z**alpha*lz
         d(i,2)=-1.0_dp/sigma+alpha/sigma*z**alpha
         d(i,3)=alpha*s/(sigma*(1.0_dp+s*epsi))*z**(alpha-1.0_dp)
         d(i,4)=alpha*s/(1.0_dp+s*epsi)*z**alpha
      end do
   end subroutine

   subroutine regression_score(z,r,alpha,sigma,epsi,d)
      real(dp),intent(in)::z(:,:),r(:),alpha,sigma,epsi;real(dp),intent(out)::d(:,:)
      integer::i,p;real(dp)::s,t,lz,common
      p=size(z,2)
      do i=1,size(r)
         s=sign(1.0_dp,r(i));if(r(i)==0.0_dp)s=0.0_dp
         t=abs(r(i))/(sigma*(1.0_dp+s*epsi));lz=log(max(t,tiny(1.0_dp)))
         common=alpha*s/(sigma*(1.0_dp+s*epsi))*t**(alpha-1.0_dp)
         d(i,1:p)=z(i,:)*common
         d(i,p+1)=digamma_aep(1.0_dp+1.0_dp/alpha)/alpha**2-t**alpha*lz
         d(i,p+2)=-1.0_dp/sigma+alpha/sigma*t**alpha
         d(i,p+3)=alpha*s/(1.0_dp+s*epsi)*t**alpha
      end do
   end subroutine

   subroutine inv_xtx(d,inv)
      real(dp),intent(in)::d(:,:);real(dp),intent(out)::inv(:,:)
      real(dp)::a(size(d,2),size(d,2)),e(size(d,2)),sol(size(d,2));integer::j
      a=matmul(transpose(d),d);inv=0.0_dp
      do j=1,size(d,2);e=0.0_dp;e(j)=1.0_dp;call solve_linear(a,e,sol);inv(:,j)=sol;end do
   end subroutine

   subroutine least_squares(x,y,b)
      real(dp),intent(in)::x(:,:),y(:);real(dp),intent(out)::b(:)
      real(dp)::a(size(b),size(b)),rhs(size(b));a=matmul(transpose(x),x);rhs=matmul(transpose(x),y);call solve_linear(a,rhs,b)
   end subroutine

   subroutine solve_linear(ain,bin,x)
      real(dp),intent(in)::ain(:,:),bin(:);real(dp),intent(out)::x(:)
      real(dp)::a(size(bin),size(bin)),b(size(bin)),tmp,factor;integer::n,i,j,k,piv
      n=size(bin);a=ain;b=bin
      do k=1,n-1
         piv=k;do i=k+1,n;if(abs(a(i,k))>abs(a(piv,k)))piv=i;end do
         if(abs(a(piv,k))<1.0e-14_dp)then;x=0.0_dp;return;end if
         if(piv/=k)then
            do j=k,n;tmp=a(k,j);a(k,j)=a(piv,j);a(piv,j)=tmp;end do
            tmp=b(k);b(k)=b(piv);b(piv)=tmp
         end if
         do i=k+1,n;factor=a(i,k)/a(k,k);a(i,k:n)=a(i,k:n)-factor*a(k,k:n);b(i)=b(i)-factor*b(k);end do
      end do
      if(abs(a(n,n))<1.0e-14_dp)then;x=0.0_dp;return;end if
      x(n)=b(n)/a(n,n);do i=n-1,1,-1;x(i)=(b(i)-sum(a(i,i+1:n)*x(i+1:n)))/a(i,i);end do
   end subroutine

   pure function outer(a,b) result(c)
      real(dp),intent(in)::a(:),b(:);real(dp)::c(size(a),size(b));integer::i
      do i=1,size(a);c(i,:)=a(i)*b;end do
   end function
   pure function sgnvec(x) result(s)
      real(dp),intent(in)::x(:);real(dp)::s(size(x));integer::i
      do i=1,size(x);if(x(i)>0)s(i)=1; if(x(i)<0)s(i)=-1; if(x(i)==0)s(i)=0;end do
   end function
   pure real(dp) function sample_var(x) result(v)
      real(dp),intent(in)::x(:);real(dp)::m;m=sum(x)/size(x);v=sum((x-m)**2)/real(size(x)-1,dp)
   end function
   real(dp) function sample_median(x) result(v)
      real(dp),intent(in)::x(:);real(dp),allocatable::z(:);integer::n
      z=x;call sort_real(z);n=size(z);if(mod(n,2)==1)then;v=z((n+1)/2);else;v=0.5_dp*(z(n/2)+z(n/2+1));end if
   end function
   subroutine sort_real(x)
      real(dp),intent(inout)::x(:);integer::i,j;real(dp)::t
      do i=2,size(x);t=x(i);j=i-1;do while(j>=1);if(x(j)<=t)exit;x(j+1)=x(j);j=j-1;end do;x(j+1)=t;end do
   end subroutine
end module aep_fit
