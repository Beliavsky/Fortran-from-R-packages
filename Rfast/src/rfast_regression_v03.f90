module rfast_regression_v03
   use rfast_special, only : dp, pi, normal_pdf, normal_cdf, chisq_cdf, f_cdf
   use rfast_linalg, only : solve_linear, inverse_matrix
   use rfast_regression, only : regression_result, glm_poisson
   use rfast_regression_v02, only : gamma_regression, invgauss_regression, quasipoisson_regression, &
                                    proportion_regression, multinomial_regression, multinomial_result
   use rfast_extra_mle, only : weibull_mle
   use rfast_mle, only : mle_result
   implicit none
   private

   integer, parameter, public :: ORDINAL_LOGIT = 1, ORDINAL_PROBIT = 2

   type, public :: tobit_result
      real(dp) :: location = 0.0_dp
      real(dp) :: scale = 1.0_dp
      real(dp) :: loglik = -huge(1.0_dp)
      integer :: iterations = 0
      integer :: status = 0
   end type tobit_result

   type, public :: ordinal_result
      real(dp), allocatable :: threshold(:)
      real(dp) :: loglik = -huge(1.0_dp)
      integer :: status = 0
   end type ordinal_result

   type, public :: weibull_regression_result
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: fitted(:)
      real(dp) :: shape = 1.0_dp
      real(dp) :: loglik = -huge(1.0_dp)
      integer :: iterations = 0
      integer :: status = 0
   end type weibull_regression_result

   public :: normlog_regression, weibull_regression, tobit_mle, ordinal_mle
   public :: normlog_regs, gamma_regs, invgauss_regs, quasipoisson_regs, proportion_regs
   public :: multinomial_regs, poisson_regs, univglms

contains

   subroutine add_intercept(x,xx)
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable, intent(out) :: xx(:,:)
      allocate(xx(size(x,1),size(x,2)+1))
      xx(:,1)=1.0_dp
      xx(:,2:)=x
   end subroutine add_intercept

   function normlog_regression(y,x,tol,maxiter) result(res)
      real(dp), intent(in) :: y(:),x(:,:)
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxiter
      type(regression_result) :: res
      real(dp), allocatable :: xx(:,:),ly(:),b1(:),b2(:),yhat(:),w1(:),w2(:),der(:),h(:,:),step(:)
      real(dp) :: eps,deviance
      integer :: n,p,it,mi,j,info
      if(size(y)==0.or.size(x,1)/=size(y).or.any(y<0.0_dp))then
         res%status=1
         return
      end if
      call add_intercept(x,xx)
      n=size(y);p=size(xx,2);eps=1.0e-7_dp;if(present(tol))eps=tol
      mi=100;if(present(maxiter))mi=maxiter
      allocate(ly(n),b1(p),b2(p),yhat(n),w1(n),w2(n),der(p),h(p,p),step(p))
      ly=log(y+0.1_dp)
      h=matmul(transpose(xx),xx)
      der=matmul(transpose(xx),ly)
      call solve_linear(h,der,b1,info)
      if(info/=0)then;res%status=info;return;end if
      do it=1,mi
         yhat=exp(max(-700.0_dp,min(700.0_dp,matmul(xx,b1))))
         w1=yhat*yhat
         w2=yhat*y
         der=matmul(transpose(xx),w1-w2)
         h=0.0_dp
         do j=1,p
            h(j,:)=matmul(2.0_dp*w1-w2,spread(xx(:,j),2,p)*xx)
         end do
         call solve_linear(h,der,step,info)
         if(info/=0)then;res%status=info;exit;end if
         b2=b1-step
         if(sum(abs(b2-b1))<=eps)then;b1=b2;exit;end if
         b1=b2
      end do
      yhat=exp(max(-700.0_dp,min(700.0_dp,matmul(xx,b1))))
      deviance=sum((y-yhat)**2)
      allocate(res%beta(p),res%fitted(n),res%residuals(n),res%covariance(p,p))
      res%beta=b1;res%fitted=yhat;res%residuals=y-yhat;res%deviance=deviance;res%iterations=it
      res%loglik=-0.5_dp*real(n,dp)*(log(max(tiny(1.0_dp),deviance/real(n,dp)))+log(2.0_dp*pi)+1.0_dp)
      res%dispersion=deviance/real(max(1,n-p),dp)
      ! Rfast does not return a covariance matrix for this routine; provide a local Hessian inverse when nonsingular.
      w1=yhat*yhat;w2=yhat*y;h=0.0_dp
      do j=1,p;h(j,:)=matmul(2.0_dp*w1-w2,spread(xx(:,j),2,p)*xx);end do
      call inverse_matrix(h,res%covariance,info)
      if(info/=0)res%covariance=0.0_dp
   end function normlog_regression

   function weibull_regression(y,x,tol,maxiter) result(res)
      real(dp), intent(in) :: y(:),x(:,:)
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxiter
      type(weibull_regression_result) :: res
      real(dp), allocatable :: xx(:,:),be(:),ben(:),lam(:),com(:),logcom(:),xcom(:,:),derb(:),h(:,:),step(:)
      real(dp) :: eps,sly,ek,k,derk,derk2,lik1,lik2
      integer :: n,p,it,mi,j,info
      type(mle_result) :: ini
      if(size(y)==0.or.size(x,1)/=size(y).or.any(y<=0.0_dp))then;res%status=1;return;end if
      call add_intercept(x,xx);n=size(y);p=size(xx,2);eps=1.0e-7_dp;if(present(tol))eps=tol
      mi=100;if(present(maxiter))mi=maxiter
      ini=weibull_mle(y,eps,mi)
      if(ini%status/=0)then;res%status=ini%status;return;end if
      allocate(be(p),ben(p),lam(n),com(n),logcom(n),xcom(n,p),derb(p),h(p,p),step(p))
      be=0.0_dp;be(1)=log(ini%param(2));ek=ini%param(1);k=log(ek);sly=sum(log(y));lik1=ini%loglik
      do it=1,mi
         lam=-matmul(xx,be)
         com=exp(max(-700.0_dp,min(700.0_dp,ek*(log(y)+lam))))
         logcom=log(max(tiny(1.0_dp),com))
         derk=real(n,dp)+ek*(sly+sum(lam))-sum(com*logcom)
         derk2=derk-real(n,dp)-sum(com*logcom*logcom)
         do j=1,p;xcom(:,j)=com*xx(:,j);end do
         derb=sum(xcom,dim=1)-sum(xx,dim=1)
         h=-ek*matmul(transpose(xx),xcom)
         call solve_linear(h,derb,step,info)
         if(info/=0)then;res%status=info;return;end if
         if(abs(derk2)<=tiny(1.0_dp))then;res%status=2;return;end if
         k=k-derk/derk2
         ben=be-step
         ek=exp(max(-50.0_dp,min(50.0_dp,k)))
         lam=-matmul(xx,ben)
         com=exp(max(-700.0_dp,min(700.0_dp,ek*(log(y)+lam))))
         lik2=real(n,dp)*k+(ek-1.0_dp)*sly+ek*sum(lam)-sum(com)
         if(lik2-lik1<=eps.and.maxval(abs(ben-be))<=sqrt(eps))then;be=ben;exit;end if
         be=ben;lik1=lik2
      end do
      ek=exp(max(-50.0_dp,min(50.0_dp,k)));lam=-matmul(xx,be)
      com=exp(max(-700.0_dp,min(700.0_dp,ek*(log(y)+lam))))
      allocate(res%beta(p),res%fitted(n));res%beta=be;res%shape=ek;res%fitted=exp(-lam);res%iterations=it
      res%loglik=real(n,dp)*k+(ek-1.0_dp)*sly+ek*sum(lam)-sum(com)
   end function weibull_regression

   function tobit_mle(y,tol,maxiter) result(res)
      real(dp), intent(in) :: y(:)
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxiter
      type(tobit_result) :: res
      real(dp), allocatable :: y1(:)
      real(dp) :: eps,m,s,sy1,sy12,com,cdf,derm,derm2,ders,ders2,derms,det
      real(dp) :: aold(2),anew(2)
      integer :: n,n1,n0,it,mi
      n=size(y);if(n==0.or.any(y<0.0_dp))then;res%status=1;return;end if
      y1=pack(y,y>0.0_dp);n1=size(y1);n0=n-n1
      if(n1<2)then;res%status=2;return;end if
      sy1=sum(y1);sy12=sum(y1*y1);m=sum(y)/real(n,dp)
      s=sqrt(max(1.0e-12_dp,sy12/real(n,dp)-m*m))
      eps=1.0e-9_dp;if(present(tol))eps=tol;mi=100;if(present(maxiter))mi=maxiter
      aold=[m,log(s)];anew=aold
      do it=1,mi
         m=aold(1);s=exp(aold(2));cdf=max(tiny(1.0_dp),normal_cdf(-m/s));com=normal_pdf(m/s)/(s*cdf)
         derm=(sy1-real(n1,dp)*m)/(s*s)-real(n0,dp)*com
         derm2=-real(n1,dp)/(s*s)-real(n0,dp)*(-m/(s*s)*com+com*com)
         ders=-real(n1,dp)+(sy12-2.0_dp*m*sy1+real(n1,dp)*m*m)/(s*s)+real(n0,dp)*m*com
         ders2=-2.0_dp*(sy12-2.0_dp*m*sy1+real(n1,dp)*m*m)/(s*s) &
               +real(n0,dp)*m*(-com+m*m/(s*s)*com-com*com*m)
         derms=-2.0_dp*(sy1-real(n1,dp)*m)/(s*s) &
                +real(n0,dp)*(com-m*m/(s*s)*com+com*com*m)
         det=derm2*ders2-derms*derms
         if(abs(det)<=tiny(1.0_dp))then;res%status=3;return;end if
         anew=aold-[ders2*derm-derms*ders,-derms*derm+derm2*ders]/det
         if(sum(abs(anew-aold))<=eps)exit
         aold=anew
      end do
      m=anew(1);s=exp(anew(2));cdf=max(tiny(1.0_dp),normal_cdf(-m/s))
      res%location=m;res%scale=s;res%iterations=it
      res%loglik=-0.5_dp*real(n1,dp)*log(2.0_dp*pi*s*s) &
                 -0.5_dp*(sy12-2.0_dp*m*sy1+real(n1,dp)*m*m)/(s*s)+real(n0,dp)*log(cdf)
   end function tobit_mle

   function ordinal_mle(y,link) result(res)
      integer, intent(in) :: y(:)
      integer, intent(in), optional :: link
      type(ordinal_result) :: res
      integer, allocatable :: counts(:)
      real(dp), allocatable :: cum(:),prob(:)
      integer :: k,j,l,n
      real(dp) :: p
      n=size(y);if(n==0.or.minval(y)<1)then;res%status=1;return;end if
      k=maxval(y);if(k<2)then;res%status=2;return;end if
      l=ORDINAL_LOGIT;if(present(link))l=link
      allocate(counts(k),cum(k),prob(k),res%threshold(k-1));counts=0
      do j=1,n;counts(y(j))=counts(y(j))+1;end do
      cum=0.0_dp;do j=1,k;cum(j)=real(sum(counts(1:j)),dp)/real(n,dp);end do
      prob(1)=cum(1);do j=2,k;prob(j)=cum(j)-cum(j-1);end do
      if(any(prob<=0.0_dp))then;res%status=3;return;end if
      do j=1,k-1
         p=min(1.0_dp-1.0e-15_dp,max(1.0e-15_dp,cum(j)))
         select case(l)
         case(ORDINAL_LOGIT);res%threshold(j)=log(p/(1.0_dp-p))
         case(ORDINAL_PROBIT);res%threshold(j)=normal_quantile_local(p)
         case default;res%status=4;return
         end select
      end do
      res%loglik=sum(real(counts,dp)*log(prob))
   end function ordinal_mle

   real(dp) function normal_quantile_local(p) result(q)
      use rfast_special, only : normal_quantile
      real(dp), intent(in) :: p
      q=normal_quantile(p)
   end function normal_quantile_local

   function normlog_regs(y,x,tol,maxiter,logged) result(out)
      real(dp),intent(in)::y(:),x(:,:)
      real(dp),intent(in),optional::tol
      integer,intent(in),optional::maxiter
      logical,intent(in),optional::logged
      real(dp)::out(size(x,2),2)
      type(regression_result)::fit
      real(dp)::null_dev,stat,pv
      real(dp),allocatable::one(:,:)
      integer::j,n
      logical::lg
      n=size(y);lg=.false.;if(present(logged))lg=logged
      null_dev=sum((y-sum(y)/real(n,dp))**2)
      allocate(one(n,1))
      do j=1,size(x,2)
         one(:,1)=x(:,j);fit=normlog_regression(y,one,tol,maxiter)
         stat=real(n,dp)*log(max(tiny(1.0_dp),null_dev/max(tiny(1.0_dp),fit%deviance)))
         pv=max(tiny(1.0_dp),1.0_dp-chisq_cdf(max(0.0_dp,stat),1.0_dp))
         out(j,1)=stat;out(j,2)=merge(log(pv),pv,lg)
      end do
   end function normlog_regs

   function gamma_regs(y,x,tol,maxiter,logged) result(out)
      real(dp),intent(in)::y(:),x(:,:)
      real(dp),intent(in),optional::tol
      integer,intent(in),optional::maxiter
      logical,intent(in),optional::logged
      real(dp)::out(size(x,2),2),null_dev,stat,pv
      type(regression_result)::fit
      real(dp),allocatable::one(:,:)
      integer::j,n
      logical::lg
      n=size(y);lg=.false.;if(present(logged))lg=logged
      null_dev=2.0_dp*sum(y/(sum(y)/real(n,dp))-1.0_dp-log(y/(sum(y)/real(n,dp))))
      allocate(one(n,1))
      do j=1,size(x,2)
         one(:,1)=x(:,j);fit=gamma_regression(y,one,tol,maxiter)
         stat=(null_dev-fit%deviance)/max(tiny(1.0_dp),fit%dispersion)
         pv=max(tiny(1.0_dp),1.0_dp-f_cdf(max(0.0_dp,stat),1.0_dp,real(max(1,n-2),dp)))
         out(j,:)=[stat,merge(log(pv),pv,lg)]
      end do
   end function gamma_regs

   function invgauss_regs(y,x,tol,maxiter,logged) result(out)
      real(dp),intent(in)::y(:),x(:,:)
      real(dp),intent(in),optional::tol
      integer,intent(in),optional::maxiter
      logical,intent(in),optional::logged
      real(dp)::out(size(x,2),2),stat,pv,null_obj
      type(regression_result)::fit
      real(dp),allocatable::one(:,:)
      integer::j,n
      logical::lg
      n=size(y);lg=.false.;if(present(logged))lg=logged
      null_obj=sum(1.0_dp/y)/real(n,dp)-1.0_dp/(sum(y)/real(n,dp))
      allocate(one(n,1))
      do j=1,size(x,2)
         one(:,1)=x(:,j);fit=invgauss_regression(y,one,tol,maxiter)
         stat=(real(n,dp)*null_obj-fit%deviance)/max(tiny(1.0_dp),fit%dispersion)
         pv=max(tiny(1.0_dp),1.0_dp-f_cdf(max(0.0_dp,stat),1.0_dp,real(max(1,n-2),dp)))
         out(j,:)=[stat,merge(log(pv),pv,lg)]
      end do
   end function invgauss_regs

   function quasipoisson_regs(y,x,tol,maxiter,logged) result(out)
      real(dp),intent(in)::y(:),x(:,:)
      real(dp),intent(in),optional::tol
      integer,intent(in),optional::maxiter
      logical,intent(in),optional::logged
      real(dp)::out(size(x,2),2),null_dev,stat,pv,m
      type(regression_result)::fit
      real(dp),allocatable::one(:,:)
      integer::j,n,i
      logical::lg
      n=size(y);m=sum(y)/real(n,dp);null_dev=0.0_dp;lg=.false.;if(present(logged))lg=logged
      do i=1,n
         if(y(i)>0.0_dp)then;null_dev=null_dev+2.0_dp*(y(i)*log(y(i)/m)-(y(i)-m));else;null_dev=null_dev+2.0_dp*m;end if
      end do
      allocate(one(n,1))
      do j=1,size(x,2)
         one(:,1)=x(:,j);fit=quasipoisson_regression(y,one,tol,maxiter)
         stat=(null_dev-fit%deviance)/max(tiny(1.0_dp),fit%dispersion)
         pv=max(tiny(1.0_dp),1.0_dp-chisq_cdf(max(0.0_dp,stat),1.0_dp))
         out(j,:)=[stat,merge(log(pv),pv,lg)]
      end do
   end function quasipoisson_regs

   function proportion_regs(y,x,quasi,tol,maxiter,logged) result(out)
      real(dp),intent(in)::y(:),x(:,:)
      logical,intent(in),optional::quasi,logged
      real(dp),intent(in),optional::tol
      integer,intent(in),optional::maxiter
      real(dp)::out(size(x,2),2),null_dev,stat,pv,m
      type(regression_result)::fit
      real(dp),allocatable::one(:,:)
      integer::j,n,i
      logical::q,lg
      n=size(y);m=sum(y)/real(n,dp);q=.true.;if(present(quasi))q=quasi;lg=.false.;if(present(logged))lg=logged
      null_dev=0.0_dp
      do i=1,n
         if(y(i)>0.0_dp)null_dev=null_dev+2.0_dp*y(i)*log(y(i)/max(tiny(1.0_dp),m))
         if(y(i)<1.0_dp)null_dev=null_dev+2.0_dp*(1.0_dp-y(i))*log((1.0_dp-y(i))/max(tiny(1.0_dp),1.0_dp-m))
      end do
      allocate(one(n,1))
      do j=1,size(x,2)
         one(:,1)=x(:,j);fit=proportion_regression(y,one,q,tol,maxiter)
         stat=(null_dev-fit%deviance)/merge(max(tiny(1.0_dp),fit%dispersion),1.0_dp,q)
         pv=max(tiny(1.0_dp),1.0_dp-chisq_cdf(max(0.0_dp,stat),1.0_dp))
         out(j,:)=[stat,merge(log(pv),pv,lg)]
      end do
   end function proportion_regs


   function poisson_regs(y,x,tol,maxiter,logged) result(out)
      real(dp),intent(in)::y(:),x(:,:)
      real(dp),intent(in),optional::tol
      integer,intent(in),optional::maxiter
      logical,intent(in),optional::logged
      real(dp)::out(size(x,2),2),null_dev,stat,pv,m
      type(regression_result)::fit
      real(dp),allocatable::xx(:,:)
      integer::j,n,i
      logical::lg
      n=size(y);m=sum(y)/real(n,dp);null_dev=0.0_dp;lg=.false.;if(present(logged))lg=logged
      do i=1,n
         if(y(i)>0.0_dp)then
            null_dev=null_dev+2.0_dp*(y(i)*log(y(i)/max(tiny(1.0_dp),m))-(y(i)-m))
         else
            null_dev=null_dev+2.0_dp*m
         end if
      end do
      allocate(xx(n,2));xx(:,1)=1.0_dp
      do j=1,size(x,2)
         xx(:,2)=x(:,j);fit=glm_poisson(xx,y,tol,maxiter)
         stat=max(0.0_dp,null_dev-fit%deviance)
         pv=max(tiny(1.0_dp),1.0_dp-chisq_cdf(stat,1.0_dp))
         out(j,:)=[stat,merge(log(pv),pv,lg)]
      end do
   end function poisson_regs

   function multinomial_regs(y,x,tol,maxiter,logged) result(out)
      integer,intent(in)::y(:)
      real(dp),intent(in)::x(:,:)
      real(dp),intent(in),optional::tol
      integer,intent(in),optional::maxiter
      logical,intent(in),optional::logged
      real(dp)::out(size(x,2),2),null_ll,stat,pv
      type(multinomial_result)::fit
      real(dp),allocatable::one(:,:)
      integer::j,n,k,c
      logical::lg
      n=size(y);k=maxval(y);lg=.false.;if(present(logged))lg=logged;null_ll=0.0_dp
      do c=1,k
         if(count(y==c)>0)null_ll=null_ll+real(count(y==c),dp)*log(real(count(y==c),dp)/real(n,dp))
      end do
      allocate(one(n,1))
      do j=1,size(x,2)
         one(:,1)=x(:,j);fit=multinomial_regression(y,one,tol,maxiter)
         stat=2.0_dp*(fit%loglik-null_ll)
         pv=max(tiny(1.0_dp),1.0_dp-chisq_cdf(max(0.0_dp,stat),real(max(1,k-1),dp)))
         out(j,:)=[stat,merge(log(pv),pv,lg)]
      end do
   end function multinomial_regs

   function univglms(y,x,family,logged) result(out)
      real(dp),intent(in)::y(:),x(:,:)
      integer,intent(in),optional::family
      logical,intent(in),optional::logged
      real(dp)::out(size(x,2),2)
      integer::fam,j,n
      real(dp)::m,rho,sy,sx,stat,pv
      logical::lg
      fam=0;if(present(family))fam=family;lg=.false.;if(present(logged))lg=logged;n=size(y)
      if(fam==0)then
         if(all(abs(y-nint(y))<1.0e-12_dp).and.minval(y)>=0.0_dp)then
            if(maxval(y)<=1.0_dp)then;fam=1;else;fam=2;end if
         else;fam=3;end if
      end if
      select case(fam)
      case(1)
         out=proportion_regs(y,x,.false.,logged=lg)
      case(2)
         out=poisson_regs(y,x,logged=lg)
      case(3)
         m=sum(y)/real(n,dp);sy=sqrt(sum((y-m)**2))
         do j=1,size(x,2)
            sx=sqrt(sum((x(:,j)-sum(x(:,j))/real(n,dp))**2))
            rho=dot_product(y-m,x(:,j)-sum(x(:,j))/real(n,dp))/max(tiny(1.0_dp),sy*sx)
            rho=max(-1.0_dp+1e-14_dp,min(1.0_dp-1e-14_dp,rho))
            stat=rho*sqrt(real(max(1,n-2),dp))/sqrt(max(tiny(1.0_dp),1.0_dp-rho*rho))
            pv=max(tiny(1.0_dp),2.0_dp*(1.0_dp-student_t_cdf_local(abs(stat),real(max(1,n-2),dp))))
            out(j,:)=[stat,merge(log(pv),pv,lg)]
         end do
      case default
         out=huge(1.0_dp)
      end select
   end function univglms

   real(dp) function student_t_cdf_local(x,df) result(p)
      use rfast_special, only : student_t_cdf
      real(dp),intent(in)::x,df
      p=student_t_cdf(x,df)
   end function student_t_cdf_local

end module rfast_regression_v03
