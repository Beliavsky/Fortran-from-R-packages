module new_dist_distributions
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use new_dist_kinds, only : dp
   use new_dist_numerics, only : pi, target_probability, finish_probability, &
      rand_uniform, sign0, nan_dp, log1p_stable
   use pracma_special, only : lambertWn
   use vgam_special, only : regularized_gamma_p, gamma_quantile, normal_cdf, lerch_phi
   use expint_mod, only : gamma_inc
   implicit none
   private

   public :: dEPd,pEPd,qEPd,rEPd
   public :: dLd,pLd,qLd,rLd
   public :: dRA,pRA,qRA,rRA
   public :: dbwd,pbwd,qbwd,rbwd
   public :: ddLd1,pdLd1,qdLd1,rdLd1
   public :: ddLd2,pdLd2,qdLd2,rdLd2
   public :: dgld,pgld,qgld,rgld
   public :: dkd,pkd,qkd,rkd
   public :: dmd,pmd,qmd,rmd
   public :: domd,pomd,qomd,romd
   public :: dpldd,ppldd,qpldd,rpldd
   public :: dsgrd,psgrd,qsgrd,rsgrd
   public :: dsod,psod,qsod,rsod
   public :: dtpmd,ptpmd,qtpmd,rtpmd
   public :: dtprd,ptprd,qtprd,rtprd
   public :: dugd,pugd,qugd,rugd
   public :: duigd,puigd,quigd,ruigd
   public :: dwgd,pwgd,qwgd,rwgd

contains

   ! EP distribution ---------------------------------------------------------
   real(dp) function dEPd(x,lambda,beta,log_p) result(v)
      real(dp), intent(in) :: x,lambda,beta
      logical, intent(in), optional :: log_p
      logical :: lp
      real(dp) :: lv
      lp=.false.; if(present(log_p)) lp=log_p
      if(lambda<=0.0_dp .or. beta<=0.0_dp) then
         v=nan_dp(); return
      end if
      if(x<=0.0_dp) then
         v=merge(-huge(1.0_dp),0.0_dp,lp); return
      end if
      lv=log(lambda)+log(beta)-log(1.0_dp-exp(-lambda))-lambda-beta*x+ &
         lambda*exp(-beta*x)
      v=merge(lv,exp(lv),lp)
   end function dEPd

   real(dp) function pEPd(q,lambda,beta,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,lambda,beta
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: cdf
      if(lambda<=0.0_dp .or. beta<=0.0_dp) then
         v=nan_dp(); return
      end if
      if(q<=0.0_dp) then
         cdf=0.0_dp
      else
         cdf=(exp(lambda*exp(-beta*q))-exp(lambda))/(1.0_dp-exp(lambda))
      end if
      v=finish_probability(cdf,lower_tail,log_p)
   end function pEPd

   real(dp) function qEPd(p,lambda,beta,lower_tail) result(x)
      real(dp), intent(in) :: p,lambda,beta
      logical, intent(in), optional :: lower_tail
      real(dp) :: t,a
      if(p<0.0_dp .or. p>1.0_dp .or. lambda<=0.0_dp .or. beta<=0.0_dp) then
         x=nan_dp(); return
      end if
      t=target_probability(p,lower_tail)
      if(t<=0.0_dp) then; x=0.0_dp; return; end if
      if(t>=1.0_dp) then; x=huge(1.0_dp); return; end if
      a=exp(lambda)+t-t*exp(lambda)
      x=-log(log(a)/lambda)/beta
   end function qEPd

   subroutine rEPd(x,lambda,beta)
      real(dp), intent(out) :: x(:)
      real(dp), intent(in) :: lambda,beta
      integer :: i
      do i=1,size(x); x(i)=qEPd(rand_uniform(),lambda,beta); end do
   end subroutine rEPd

   ! Lindley -----------------------------------------------------------------
   real(dp) function dLd(x,theta,log_p) result(v)
      real(dp), intent(in) :: x,theta
      logical, intent(in), optional :: log_p
      logical :: lp
      real(dp) :: lv
      lp=.false.; if(present(log_p)) lp=log_p
      if(theta<=0.0_dp) then; v=nan_dp(); return; end if
      if(x<=0.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lp); return; end if
      lv=2.0_dp*log(theta)-log1p_stable(theta)+log1p_stable(x)-theta*x
      v=merge(lv,exp(lv),lp)
   end function dLd

   real(dp) function pLd(q,theta,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,theta
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: cdf
      if(theta<=0.0_dp) then; v=nan_dp(); return; end if
      if(q<=0.0_dp) then
         cdf=0.0_dp
      else
         cdf=1.0_dp-(1.0_dp+theta*q/(1.0_dp+theta))*exp(-theta*q)
      end if
      v=finish_probability(cdf,lower_tail,log_p)
   end function pLd

   real(dp) function qLd(p,theta,lower_tail) result(x)
      real(dp), intent(in) :: p,theta
      logical, intent(in), optional :: lower_tail
      real(dp) :: t,arg,w
      integer :: status
      if(p<0.0_dp .or. p>1.0_dp .or. theta<=0.0_dp) then; x=nan_dp(); return; end if
      t=target_probability(p,lower_tail)
      if(t<=0.0_dp) then; x=0.0_dp; return; end if
      if(t>=1.0_dp) then; x=huge(1.0_dp); return; end if
      arg=-(1.0_dp+theta)*(1.0_dp-t)*exp(-1.0_dp-theta)
      w=lambertWn(arg,status)
      if(status/=0) then; x=nan_dp(); else; x=-(w+1.0_dp+theta)/theta; end if
   end function qLd

   subroutine rLd(x,theta)
      real(dp), intent(out) :: x(:); real(dp), intent(in) :: theta
      integer :: i
      do i=1,size(x); x(i)=qLd(rand_uniform(),theta); end do
   end subroutine rLd

   ! Ram Awadh ---------------------------------------------------------------
   real(dp) function dRA(x,theta,log_p) result(v)
      real(dp), intent(in) :: x,theta
      logical, intent(in), optional :: log_p
      logical :: lp
      real(dp) :: lv
      lp=.false.; if(present(log_p)) lp=log_p
      if(theta<=0.0_dp) then; v=nan_dp(); return; end if
      if(x<=0.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lp); return; end if
      lv=6.0_dp*log(theta)-log(theta**6+120.0_dp)+log(theta+x**5)-theta*x
      v=merge(lv,exp(lv),lp)
   end function dRA

   real(dp) function pRA(q,theta,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,theta
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: cdf,poly
      if(theta<=0.0_dp) then; v=nan_dp(); return; end if
      if(q<=0.0_dp) then
         cdf=0.0_dp
      else
         poly=theta**4*q**4+5.0_dp*theta**3*q**3+20.0_dp*theta**2*q**2+ &
              60.0_dp*theta*q+120.0_dp
         cdf=1.0_dp-(1.0_dp+theta*q*poly/(theta**6+120.0_dp))*exp(-theta*q)
      end if
      v=finish_probability(cdf,lower_tail,log_p)
   end function pRA

   real(dp) function qRA(p,theta,lower_tail) result(x)
      real(dp), intent(in) :: p,theta
      logical, intent(in), optional :: lower_tail
      real(dp) :: t,lo,hi,mid
      integer :: it
      if(p<0.0_dp .or. p>1.0_dp .or. theta<=0.0_dp) then; x=nan_dp(); return; end if
      t=target_probability(p,lower_tail)
      if(t<=0.0_dp) then; x=0.0_dp; return; end if
      if(t>=1.0_dp) then; x=huge(1.0_dp); return; end if
      lo=0.0_dp; hi=max(1.0_dp,1.0_dp/theta)
      do while(pRA(hi,theta)<t); hi=2.0_dp*hi; if(hi>huge(hi)/4) exit; end do
      do it=1,100
         mid=0.5_dp*(lo+hi)
         if(pRA(mid,theta)<t) then; lo=mid; else; hi=mid; end if
      end do
      x=0.5_dp*(lo+hi)
   end function qRA

   subroutine rRA(x,theta)
      real(dp), intent(out) :: x(:); real(dp), intent(in) :: theta
      integer :: i
      do i=1,size(x); x(i)=qRA(rand_uniform(),theta); end do
   end subroutine rRA

   ! Bimodal Weibull ---------------------------------------------------------
   real(dp) function bwd_z(alpha,beta,sigma) result(z)
      real(dp), intent(in) :: alpha,beta,sigma
      z=2.0_dp+sigma*sigma*beta*beta*gamma(1.0_dp+2.0_dp/alpha)- &
        2.0_dp*sigma*beta*gamma(1.0_dp+1.0_dp/alpha)
   end function bwd_z

   real(dp) function dbwd(x,alpha,beta,sigma,log_p) result(v)
      real(dp), intent(in) :: x,alpha,beta,sigma
      logical, intent(in), optional :: log_p
      logical :: lp
      real(dp) :: z,lv
      lp=.false.; if(present(log_p)) lp=log_p
      if(alpha<=0.0_dp .or. beta<=0.0_dp) then; v=nan_dp(); return; end if
      if(x<0.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lp); return; end if
      z=bwd_z(alpha,beta,sigma)
      if(z<=0.0_dp) then; v=nan_dp(); return; end if
      if(x<=0.0_dp .and. alpha<1.0_dp) then
         v=merge(huge(1.0_dp),huge(1.0_dp),lp); return
      end if
      lv=log(alpha)-log(beta)-log(z)+log(1.0_dp+(1.0_dp-sigma*x)**2)+ &
         (alpha-1.0_dp)*log(max(x/beta,tiny(1.0_dp)))-(x/beta)**alpha
      v=merge(lv,exp(lv),lp)
   end function dbwd

   real(dp) function pbwd(q,alpha,beta,sigma,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,alpha,beta,sigma
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: cdf,z,y,g1,g2
      if(alpha<=0.0_dp .or. beta<=0.0_dp) then; v=nan_dp(); return; end if
      if(q<0.0_dp) then
         cdf=0.0_dp
      else
         z=bwd_z(alpha,beta,sigma)
         y=(q/beta)**alpha
         g1=gamma(1.0_dp/alpha)*regularized_gamma_p(1.0_dp/alpha,y)
         g2=gamma(2.0_dp/alpha)*regularized_gamma_p(2.0_dp/alpha,y)
         cdf=(2.0_dp-(1.0_dp+(1.0_dp-sigma*q)**2)*exp(-y)- &
              (2.0_dp*sigma*beta/alpha)*(g1-sigma*beta*g2))/z
      end if
      v=finish_probability(cdf,lower_tail,log_p)
   end function pbwd

   real(dp) function qbwd(p,alpha,beta,sigma,lower_tail) result(x)
      real(dp), intent(in) :: p,alpha,beta,sigma
      logical, intent(in), optional :: lower_tail
      real(dp) :: t,lo,hi,mid,mean0,z
      integer :: it
      if(p<0.0_dp .or. p>1.0_dp .or. alpha<=0.0_dp .or. beta<=0.0_dp) then
         x=nan_dp(); return
      end if
      t=target_probability(p,lower_tail)
      if(t<=0.0_dp) then; x=0.0_dp; return; end if
      if(t>=1.0_dp) then; x=huge(1.0_dp); return; end if
      z=bwd_z(alpha,beta,sigma)
      mean0=beta*(2.0_dp*gamma(1.0_dp+1.0_dp/alpha)+sigma*sigma*beta*beta* &
         gamma(1.0_dp+3.0_dp/alpha)-2.0_dp*sigma*beta*gamma(1.0_dp+2.0_dp/alpha))/z
      lo=0.0_dp; hi=max(beta,max(1.0_dp,2.0_dp*abs(mean0)))
      do while(pbwd(hi,alpha,beta,sigma)<t); hi=2.0_dp*hi; if(hi>huge(hi)/4) exit; end do
      do it=1,100
         mid=0.5_dp*(lo+hi)
         if(pbwd(mid,alpha,beta,sigma)<t) then; lo=mid; else; hi=mid; end if
      end do
      x=0.5_dp*(lo+hi)
   end function qbwd

   subroutine rbwd(x,alpha,beta,sigma)
      real(dp), intent(out) :: x(:); real(dp), intent(in) :: alpha,beta,sigma
      integer :: i
      do i=1,size(x); x(i)=qbwd(rand_uniform(),alpha,beta,sigma); end do
   end subroutine rbwd

   ! Discrete Lindley I ------------------------------------------------------
   real(dp) function ddLd1(x,theta,log_p) result(v)
      real(dp), intent(in) :: x,theta
      logical, intent(in), optional :: log_p
      logical :: lp
      integer :: k
      real(dp) :: lambda,pv
      lp=.false.; if(present(log_p)) lp=log_p
      if(theta<=0.0_dp) then; v=nan_dp(); return; end if
      k=floor(x); lambda=exp(-theta)
      if(k<0) then; pv=0.0_dp
      else
         pv=lambda**k/(1.0_dp-log(lambda))*(lambda*log(lambda)+ &
            (1.0_dp-lambda)*(1.0_dp-log(lambda**(k+1))))
      end if
      if(lp) then; if(pv<=0.0_dp) then; v=-huge(1.0_dp); else; v=log(pv); end if
      else; v=pv; end if
   end function ddLd1

   real(dp) function pdLd1(q,theta,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,theta
      logical, intent(in), optional :: lower_tail,log_p
      integer :: k
      real(dp) :: lambda,cdf
      if(theta<=0.0_dp) then; v=nan_dp(); return; end if
      k=floor(q); lambda=exp(-theta)
      if(k<0) then; cdf=0.0_dp
      else
         cdf=(1.0_dp-lambda**(k+1)+((real(k+2,dp)*lambda**(k+1))-1.0_dp)* &
              log(lambda))/(1.0_dp-log(lambda))
      end if
      v=finish_probability(cdf,lower_tail,log_p)
   end function pdLd1

   integer function qdLd1(p,theta,lower_tail) result(k)
      real(dp), intent(in) :: p,theta
      logical, intent(in), optional :: lower_tail
      real(dp) :: t
      integer :: lo,hi,mid
      if(p<0.0_dp .or. p>1.0_dp .or. theta<=0.0_dp) then; k=-huge(1); return; end if
      t=target_probability(p,lower_tail)
      if(t<=0.0_dp) then; k=0; return; end if
      if(t>=1.0_dp) then; k=huge(1); return; end if
      lo=0; hi=1
      do while(pdLd1(real(hi,dp),theta)<t); hi=2*hi+1; if(hi>500000000) exit; end do
      do while(lo<hi)
         mid=lo+(hi-lo)/2
         if(pdLd1(real(mid,dp),theta)>=t) then; hi=mid; else; lo=mid+1; end if
      end do
      k=lo
   end function qdLd1

   subroutine rdLd1(x,theta)
      integer, intent(out) :: x(:); real(dp), intent(in) :: theta
      integer :: i
      do i=1,size(x); x(i)=qdLd1(rand_uniform(),theta); end do
   end subroutine rdLd1

   ! Discrete Lindley II -----------------------------------------------------
   real(dp) function ddLd2(x,theta,log_p) result(v)
      real(dp), intent(in) :: x,theta
      logical, intent(in), optional :: log_p
      logical :: lp
      integer :: k
      real(dp) :: lambda,pv
      lp=.false.; if(present(log_p)) lp=log_p
      if(theta<=0.0_dp) then; v=nan_dp(); return; end if
      k=floor(x); lambda=exp(-theta)
      if(k<0) then; pv=0.0_dp
      else
         pv=lambda**k/(1.0_dp+theta)*(theta*(1.0_dp-2.0_dp*lambda)+ &
            (1.0_dp-lambda)*(1.0_dp+theta*real(k,dp)))
      end if
      if(lp) then; if(pv<=0.0_dp) then; v=-huge(1.0_dp); else; v=log(pv); end if
      else; v=pv; end if
   end function ddLd2

   real(dp) function pdLd2(q,theta,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,theta
      logical, intent(in), optional :: lower_tail,log_p
      integer :: k
      real(dp) :: lambda,cdf
      if(theta<=0.0_dp) then; v=nan_dp(); return; end if
      k=floor(q)+1; lambda=exp(-theta)
      if(k<0) then; cdf=0.0_dp
      else; cdf=1.0_dp-(1.0_dp+theta+theta*real(k,dp))/(1.0_dp+theta)*lambda**k; end if
      v=finish_probability(cdf,lower_tail,log_p)
   end function pdLd2

   integer function qdLd2(p,theta,lower_tail) result(k)
      real(dp), intent(in) :: p,theta
      logical, intent(in), optional :: lower_tail
      real(dp) :: t
      integer :: lo,hi,mid
      if(p<0.0_dp .or. p>1.0_dp .or. theta<=0.0_dp) then; k=-huge(1); return; end if
      t=target_probability(p,lower_tail)
      if(t<=0.0_dp) then; k=0; return; end if
      if(t>=1.0_dp) then; k=huge(1); return; end if
      lo=0; hi=1
      do while(pdLd2(real(hi,dp),theta)<t); hi=2*hi+1; if(hi>500000000) exit; end do
      do while(lo<hi)
         mid=lo+(hi-lo)/2
         if(pdLd2(real(mid,dp),theta)>=t) then; hi=mid; else; lo=mid+1; end if
      end do
      k=lo
   end function qdLd2

   subroutine rdLd2(x,theta)
      integer, intent(out) :: x(:); real(dp), intent(in) :: theta
      integer :: i
      do i=1,size(x); x(i)=qdLd2(rand_uniform(),theta); end do
   end subroutine rdLd2

   ! Gamma-Lomax -------------------------------------------------------------
   real(dp) function dgld(x,a,alpha,beta,log_p) result(v)
      real(dp), intent(in) :: x,a,alpha,beta
      logical, intent(in), optional :: log_p
      logical :: lp
      real(dp) :: t,lv
      lp=.false.; if(present(log_p)) lp=log_p
      if(a<=0.0_dp .or. alpha<=0.0_dp .or. beta<=0.0_dp) then; v=nan_dp(); return; end if
      if(x<=0.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lp); return; end if
      t=-alpha*log(beta/(beta+x))
      lv=log(alpha)+alpha*log(beta)-log_gamma(a)-(alpha+1.0_dp)*log(beta+x)+ &
         (a-1.0_dp)*log(t)
      v=merge(lv,exp(lv),lp)
   end function dgld

   real(dp) function pgld(q,a,alpha,beta,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,a,alpha,beta
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: t,cdf
      if(a<=0.0_dp .or. alpha<=0.0_dp .or. beta<=0.0_dp) then; v=nan_dp(); return; end if
      if(q<=0.0_dp) then; cdf=0.0_dp
      else
         t=-alpha*log(beta/(beta+q))
         cdf=(gamma(a)-gamma_inc(a,t))/gamma(a)
      end if
      v=finish_probability(cdf,lower_tail,log_p)
   end function pgld

   real(dp) function qgld(p,a,alpha,beta,lower_tail) result(x)
      real(dp), intent(in) :: p,a,alpha,beta
      logical, intent(in), optional :: lower_tail
      real(dp) :: t,y
      if(p<0.0_dp .or. p>1.0_dp .or. a<=0.0_dp .or. alpha<=0.0_dp .or. beta<=0.0_dp) then
         x=nan_dp(); return
      end if
      t=target_probability(p,lower_tail)
      if(t<=0.0_dp) then; x=0.0_dp; return; end if
      if(t>=1.0_dp) then; x=huge(1.0_dp); return; end if
      y=gamma_quantile(t,a,1.0_dp)
      x=beta*(exp(y/alpha)-1.0_dp)
   end function qgld

   subroutine rgld(x,a,alpha,beta)
      real(dp), intent(out) :: x(:); real(dp), intent(in) :: a,alpha,beta
      integer :: i
      do i=1,size(x); x(i)=qgld(rand_uniform(),a,alpha,beta); end do
   end subroutine rgld

   ! Kumaraswamy -------------------------------------------------------------
   real(dp) function dkd(x,lambda,alpha,log_p) result(v)
      real(dp), intent(in) :: x,lambda,alpha
      logical, intent(in), optional :: log_p
      logical :: lp
      real(dp) :: lv
      lp=.false.; if(present(log_p)) lp=log_p
      if(lambda<=0.0_dp .or. alpha<=0.0_dp) then; v=nan_dp(); return; end if
      if(x<=0.0_dp .or. x>=1.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lp); return; end if
      lv=log(alpha)+log(lambda)+(lambda-1.0_dp)*log(x)+(alpha-1.0_dp)*log(1.0_dp-x**lambda)
      v=merge(lv,exp(lv),lp)
   end function dkd

   real(dp) function pkd(q,lambda,alpha,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,lambda,alpha
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: cdf
      if(lambda<=0.0_dp .or. alpha<=0.0_dp) then; v=nan_dp(); return; end if
      if(q<=0.0_dp) then; cdf=0.0_dp
      else if(q>=1.0_dp) then; cdf=1.0_dp
      else; cdf=1.0_dp-(1.0_dp-q**lambda)**alpha; end if
      v=finish_probability(cdf,lower_tail,log_p)
   end function pkd

   real(dp) function qkd(p,lambda,alpha,lower_tail) result(x)
      real(dp), intent(in) :: p,lambda,alpha
      logical, intent(in), optional :: lower_tail
      real(dp) :: t
      if(p<0.0_dp .or. p>1.0_dp .or. lambda<=0.0_dp .or. alpha<=0.0_dp) then
         x=nan_dp(); return
      end if
      t=target_probability(p,lower_tail)
      if(t<=0.0_dp) then; x=0.0_dp
      else if(t>=1.0_dp) then; x=1.0_dp
      else; x=(1.0_dp-(1.0_dp-t)**(1.0_dp/alpha))**(1.0_dp/lambda); end if
   end function qkd

   subroutine rkd(x,lambda,alpha)
      real(dp), intent(out) :: x(:); real(dp), intent(in) :: lambda,alpha
      integer :: i
      do i=1,size(x); x(i)=qkd(rand_uniform(),lambda,alpha); end do
   end subroutine rkd

   ! Maxwell -----------------------------------------------------------------
   real(dp) function dmd(x,theta,log_p) result(v)
      real(dp), intent(in) :: x,theta
      logical, intent(in), optional :: log_p
      logical :: lp
      real(dp) :: lv
      lp=.false.; if(present(log_p)) lp=log_p
      if(theta<=0.0_dp) then; v=nan_dp(); return; end if
      if(x<0.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lp); return; end if
      if(x<=0.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lp); return; end if
      lv=log(4.0_dp)-0.5_dp*log(pi)-1.5_dp*log(theta)+2.0_dp*log(x)-x*x/theta
      v=merge(lv,exp(lv),lp)
   end function dmd

   real(dp) function pmd(q,theta,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,theta
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: cdf
      if(theta<=0.0_dp) then; v=nan_dp(); return; end if
      if(q<0.0_dp) then; cdf=0.0_dp
      else; cdf=regularized_gamma_p(1.5_dp,q*q/theta); end if
      v=finish_probability(cdf,lower_tail,log_p)
   end function pmd

   real(dp) function qmd(p,theta,lower_tail) result(x)
      real(dp), intent(in) :: p,theta
      logical, intent(in), optional :: lower_tail
      real(dp) :: t,y
      if(p<0.0_dp .or. p>1.0_dp .or. theta<=0.0_dp) then; x=nan_dp(); return; end if
      t=target_probability(p,lower_tail)
      if(t<=0.0_dp) then; x=0.0_dp; return; end if
      if(t>=1.0_dp) then; x=huge(1.0_dp); return; end if
      y=gamma_quantile(t,1.5_dp,1.0_dp)
      x=sqrt(theta*y)
   end function qmd

   subroutine rmd(x,theta)
      real(dp), intent(out) :: x(:); real(dp), intent(in) :: theta
      integer :: i
      do i=1,size(x); x(i)=qmd(rand_uniform(),theta); end do
   end subroutine rmd

   ! Muth --------------------------------------------------------------------
   real(dp) function domd(x,alpha,log_p) result(v)
      real(dp), intent(in) :: x,alpha
      logical, intent(in), optional :: log_p
      logical :: lp
      real(dp) :: ea,lv
      lp=.false.; if(present(log_p)) lp=log_p
      if(alpha<=0.0_dp .or. alpha>1.0_dp) then; v=nan_dp(); return; end if
      if(x<=0.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lp); return; end if
      ea=exp(alpha*x)
      lv=log(ea-alpha)+alpha*x-(ea-1.0_dp)/alpha
      v=merge(lv,exp(lv),lp)
   end function domd

   real(dp) function pomd(q,alpha,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,alpha
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: cdf
      if(alpha<=0.0_dp .or. alpha>1.0_dp) then; v=nan_dp(); return; end if
      if(q<=0.0_dp) then; cdf=0.0_dp
      else; cdf=1.0_dp-exp(alpha*q-(exp(alpha*q)-1.0_dp)/alpha); end if
      v=finish_probability(cdf,lower_tail,log_p)
   end function pomd

   real(dp) function qomd(p,alpha,lower_tail) result(x)
      real(dp), intent(in) :: p,alpha
      logical, intent(in), optional :: lower_tail
      real(dp) :: t,l1,arg,w
      integer :: status
      if(p<0.0_dp .or. p>1.0_dp .or. alpha<=0.0_dp .or. alpha>1.0_dp) then
         x=nan_dp(); return
      end if
      t=target_probability(p,lower_tail)
      if(t<=0.0_dp) then; x=0.0_dp; return; end if
      if(t>=1.0_dp) then; x=huge(1.0_dp); return; end if
      l1=log(1.0_dp-t)
      arg=-exp((-1.0_dp+alpha*l1)/alpha)/alpha
      w=lambertWn(arg,status)
      if(status/=0) then; x=nan_dp(); else; x=(-alpha*w-1.0_dp+alpha*l1)/(alpha*alpha); end if
   end function qomd

   subroutine romd(x,alpha)
      real(dp), intent(out) :: x(:); real(dp), intent(in) :: alpha
      integer :: i
      do i=1,size(x); x(i)=qomd(rand_uniform(),alpha); end do
   end subroutine romd

   ! Power log-Dagum ---------------------------------------------------------
   real(dp) function pldd_g(x,beta,theta) result(g)
      real(dp), intent(in) :: x,beta,theta
      if(abs(beta)<=tiny(1.0_dp)) then; g=nan_dp(); return; end if
      g=beta*x+sign0(x)*(theta/beta)*abs(x)**beta
   end function pldd_g

   real(dp) function dpldd(x,alpha,beta,theta,log_p) result(v)
      real(dp), intent(in) :: x,alpha,beta,theta
      logical, intent(in), optional :: log_p
      logical :: lp
      real(dp) :: g,fac,lv
      lp=.false.; if(present(log_p)) lp=log_p
      if(alpha<=0.0_dp .or. theta<0.0_dp .or. abs(beta)<=tiny(1.0_dp)) then; v=nan_dp(); return; end if
      g=pldd_g(x,beta,theta)
      if(abs(x)<=tiny(1.0_dp)) then
         if(beta>1.0_dp .or. theta<=0.0_dp) then; fac=beta
         else if(abs(beta-1.0_dp)<=epsilon(1.0_dp)) then; fac=beta+theta
         else; v=huge(1.0_dp); return; end if
      else; fac=beta+theta*abs(x)**(beta-1.0_dp); end if
      if(fac<=0.0_dp) then; v=nan_dp(); return; end if
      lv=log(alpha)+log(fac)-g-(alpha+1.0_dp)*log(1.0_dp+exp(-g))
      v=merge(lv,exp(lv),lp)
   end function dpldd

   real(dp) function ppldd(q,alpha,beta,theta,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,alpha,beta,theta
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: g,cdf
      if(alpha<=0.0_dp .or. theta<0.0_dp .or. abs(beta)<=tiny(1.0_dp)) then; v=nan_dp(); return; end if
      g=pldd_g(q,beta,theta)
      if(g>=0.0_dp) then; cdf=(1.0_dp+exp(-g))**(-alpha)
      else; cdf=(exp(g)/(1.0_dp+exp(g)))**alpha; end if
      v=finish_probability(cdf,lower_tail,log_p)
   end function ppldd

   real(dp) function qpldd(p,alpha,beta,theta,lower_tail) result(x)
      real(dp), intent(in) :: p,alpha,beta,theta
      logical, intent(in), optional :: lower_tail
      real(dp) :: t,lo,hi,mid,fl,fh
      integer :: it
      if(p<0.0_dp .or. p>1.0_dp .or. alpha<=0.0_dp .or. theta<0.0_dp .or. &
         abs(beta)<=tiny(1.0_dp)) then; x=nan_dp(); return; end if
      t=target_probability(p,lower_tail)
      if(t<=0.0_dp) then; x=-huge(1.0_dp); return; end if
      if(t>=1.0_dp) then; x=huge(1.0_dp); return; end if
      lo=-10.0_dp; hi=10.0_dp; fl=ppldd(lo,alpha,beta,theta); fh=ppldd(hi,alpha,beta,theta)
      do while(fl>t .and. abs(lo)<1.0e6_dp); lo=2.0_dp*lo; fl=ppldd(lo,alpha,beta,theta); end do
      do while(fh<t .and. abs(hi)<1.0e6_dp); hi=2.0_dp*hi; fh=ppldd(hi,alpha,beta,theta); end do
      if(fl>t .or. fh<t) then; x=nan_dp(); return; end if
      do it=1,120
         mid=0.5_dp*(lo+hi)
         if(ppldd(mid,alpha,beta,theta)<t) then; lo=mid; else; hi=mid; end if
      end do
      x=0.5_dp*(lo+hi)
   end function qpldd

   subroutine rpldd(x,alpha,beta,theta)
      real(dp), intent(out) :: x(:); real(dp), intent(in) :: alpha,beta,theta
      integer :: i
      do i=1,size(x); x(i)=qpldd(rand_uniform(),alpha,beta,theta); end do
   end subroutine rpldd

   ! Slashed generalized Rayleigh -------------------------------------------
   real(dp) function dsgrd(x,theta,alpha,beta,log_p) result(v)
      real(dp), intent(in) :: x,theta,alpha,beta
      logical, intent(in), optional :: log_p
      logical :: lp
      real(dp) :: a,reg,lv
      lp=.false.; if(present(log_p)) lp=log_p
      if(alpha<=-1.0_dp .or. theta<=0.0_dp .or. beta<=2.0_dp) then; v=nan_dp(); return; end if
      if(x<=0.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lp); return; end if
      a=(2.0_dp*alpha+beta+2.0_dp)/2.0_dp
      reg=regularized_gamma_p(a,theta*x*x)
      if(reg<=0.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lp); return; end if
      lv=log(beta)-(beta+1.0_dp)*log(x)-log_gamma(alpha+1.0_dp)- &
         0.5_dp*beta*log(theta)+log_gamma(a)+log(reg)
      v=merge(lv,exp(lv),lp)
   end function dsgrd

   real(dp) function psgrd(q,theta,alpha,beta,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,theta,alpha,beta
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: a,y,low1,lowa,cdf
      if(alpha<=-1.0_dp .or. theta<=0.0_dp .or. beta<=2.0_dp) then; v=nan_dp(); return; end if
      if(q<=0.0_dp) then; cdf=0.0_dp
      else
         a=(2.0_dp*alpha+beta+2.0_dp)/2.0_dp; y=theta*q*q
         low1=gamma(alpha+1.0_dp)*regularized_gamma_p(alpha+1.0_dp,y)
         lowa=gamma(a)*regularized_gamma_p(a,y)
         cdf=(low1-lowa*y**(-0.5_dp*beta))/gamma(alpha+1.0_dp)
      end if
      v=finish_probability(cdf,lower_tail,log_p)
   end function psgrd

   real(dp) function qsgrd(p,theta,alpha,beta,lower_tail) result(x)
      real(dp), intent(in) :: p,theta,alpha,beta
      logical, intent(in), optional :: lower_tail
      real(dp) :: t,lo,hi,mid,mean0
      integer :: it
      if(p<0.0_dp .or. p>1.0_dp .or. alpha<=-1.0_dp .or. theta<=0.0_dp .or. beta<=2.0_dp) then
         x=nan_dp(); return
      end if
      t=target_probability(p,lower_tail)
      if(t<=0.0_dp) then; x=0.0_dp; return; end if
      if(t>=1.0_dp) then; x=huge(1.0_dp); return; end if
      mean0=gamma(alpha+1.5_dp)*beta/(gamma(alpha+1.0_dp)*sqrt(theta)*(beta-1.0_dp))
      lo=0.0_dp; hi=max(1.0_dp,2.0_dp*mean0)
      do while(psgrd(hi,theta,alpha,beta)<t); hi=2.0_dp*hi; if(hi>1.0e10_dp) exit; end do
      do it=1,100
         mid=0.5_dp*(lo+hi)
         if(psgrd(mid,theta,alpha,beta)<t) then; lo=mid; else; hi=mid; end if
      end do
      x=0.5_dp*(lo+hi)
   end function qsgrd

   subroutine rsgrd(x,theta,alpha,beta)
      real(dp), intent(out) :: x(:); real(dp), intent(in) :: theta,alpha,beta
      integer :: i
      do i=1,size(x); x(i)=qsgrd(rand_uniform(),theta,alpha,beta); end do
   end subroutine rsgrd

   ! Standard Omega ----------------------------------------------------------
   real(dp) function dsod(x,alpha,beta,log_p) result(v)
      real(dp), intent(in) :: x,alpha,beta
      logical, intent(in), optional :: log_p
      logical :: lp
      real(dp) :: xb,lv
      lp=.false.; if(present(log_p)) lp=log_p
      if(alpha<=0.0_dp .or. beta<=0.0_dp) then; v=nan_dp(); return; end if
      if(x<=0.0_dp .or. x>=1.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lp); return; end if
      xb=x**beta
      lv=log(alpha)+log(beta)+(beta-1.0_dp)*log(x)-log(1.0_dp-xb*xb)- &
         0.5_dp*alpha*log((1.0_dp+xb)/(1.0_dp-xb))
      v=merge(lv,exp(lv),lp)
   end function dsod

   real(dp) function psod(q,alpha,beta,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,alpha,beta
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: cdf,xb
      if(alpha<=0.0_dp .or. beta<=0.0_dp) then; v=nan_dp(); return; end if
      if(q<=0.0_dp) then; cdf=0.0_dp
      else if(q>=1.0_dp) then; cdf=1.0_dp
      else; xb=q**beta; cdf=1.0_dp-((1.0_dp+xb)/(1.0_dp-xb))**(-0.5_dp*alpha); end if
      v=finish_probability(cdf,lower_tail,log_p)
   end function psod

   real(dp) function qsod(p,alpha,beta,lower_tail) result(x)
      real(dp), intent(in) :: p,alpha,beta
      logical, intent(in), optional :: lower_tail
      real(dp) :: t,a,y
      if(p<0.0_dp .or. p>1.0_dp .or. alpha<=0.0_dp .or. beta<=0.0_dp) then; x=nan_dp(); return; end if
      t=target_probability(p,lower_tail)
      if(t<=0.0_dp) then; x=0.0_dp; return; end if
      if(t>=1.0_dp) then; x=1.0_dp; return; end if
      a=(1.0_dp-t)**(-2.0_dp/alpha); y=(a-1.0_dp)/(a+1.0_dp)
      x=y**(1.0_dp/beta)
   end function qsod

   subroutine rsod(x,alpha,beta)
      real(dp), intent(out) :: x(:); real(dp), intent(in) :: alpha,beta
      integer :: i
      do i=1,size(x); x(i)=qsod(rand_uniform(),alpha,beta); end do
   end subroutine rsod

   ! Two-parameter Muth ------------------------------------------------------
   real(dp) function dtpmd(x,beta,alpha,log_p) result(v)
      real(dp), intent(in) :: x,beta,alpha
      logical, intent(in), optional :: log_p
      logical :: lp
      real(dp) :: u,eu,lv
      lp=.false.; if(present(log_p)) lp=log_p
      if(beta<=0.0_dp .or. alpha<=0.0_dp) then; v=nan_dp(); return; end if
      if(x<=0.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lp); return; end if
      u=(x/beta)**alpha; eu=exp(u)
      lv=log(alpha)-alpha*log(beta)+(alpha-1.0_dp)*log(x)+log(eu-1.0_dp)+u-(eu-1.0_dp)
      v=merge(lv,exp(lv),lp)
   end function dtpmd

   real(dp) function ptpmd(q,beta,alpha,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,beta,alpha
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: u,cdf
      if(beta<=0.0_dp .or. alpha<=0.0_dp) then; v=nan_dp(); return; end if
      if(q<=0.0_dp) then; cdf=0.0_dp
      else; u=(q/beta)**alpha; cdf=1.0_dp-exp(u-(exp(u)-1.0_dp)); end if
      v=finish_probability(cdf,lower_tail,log_p)
   end function ptpmd

   real(dp) function qtpmd(p,beta,alpha,lower_tail) result(x)
      real(dp), intent(in) :: p,beta,alpha
      logical, intent(in), optional :: lower_tail
      real(dp) :: t,l1,arg,w,z
      integer :: status
      if(p<0.0_dp .or. p>1.0_dp .or. beta<=0.0_dp .or. alpha<=0.0_dp) then; x=nan_dp(); return; end if
      t=target_probability(p,lower_tail)
      if(t<=0.0_dp) then; x=0.0_dp; return; end if
      if(t>=1.0_dp) then; x=huge(1.0_dp); return; end if
      l1=log(1.0_dp-t); arg=-(1.0_dp-t)/exp(1.0_dp)
      w=lambertWn(arg,status)
      if(status/=0) then; x=nan_dp(); return; end if
      z=-w-1.0_dp+l1
      x=beta*z**(1.0_dp/alpha)
   end function qtpmd

   subroutine rtpmd(x,beta,alpha)
      real(dp), intent(out) :: x(:); real(dp), intent(in) :: beta,alpha
      integer :: i
      do i=1,size(x); x(i)=qtpmd(rand_uniform(),beta,alpha); end do
   end subroutine rtpmd

   ! Two-parameter Rayleigh --------------------------------------------------
   real(dp) function dtprd(x,lambda,mu,log_p) result(v)
      real(dp), intent(in) :: x,lambda,mu
      logical, intent(in), optional :: log_p
      logical :: lp
      real(dp) :: y,lv
      lp=.false.; if(present(log_p)) lp=log_p
      if(lambda<=0.0_dp) then; v=nan_dp(); return; end if
      if(x<=mu) then; v=merge(-huge(1.0_dp),0.0_dp,lp); return; end if
      y=x-mu; lv=log(2.0_dp*lambda)+log(y)-lambda*y*y
      v=merge(lv,exp(lv),lp)
   end function dtprd

   real(dp) function ptprd(q,lambda,mu,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,lambda,mu
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: cdf
      if(lambda<=0.0_dp) then; v=nan_dp(); return; end if
      if(q<=mu) then; cdf=0.0_dp; else; cdf=1.0_dp-exp(-lambda*(q-mu)**2); end if
      v=finish_probability(cdf,lower_tail,log_p)
   end function ptprd

   real(dp) function qtprd(p,lambda,mu,lower_tail) result(x)
      real(dp), intent(in) :: p,lambda,mu
      logical, intent(in), optional :: lower_tail
      real(dp) :: t
      if(p<0.0_dp .or. p>1.0_dp .or. lambda<=0.0_dp) then; x=nan_dp(); return; end if
      t=target_probability(p,lower_tail)
      if(t<=0.0_dp) then; x=mu; return; end if
      if(t>=1.0_dp) then; x=huge(1.0_dp); return; end if
      x=mu+sqrt(-log(1.0_dp-t)/lambda)
   end function qtprd

   subroutine rtprd(x,lambda,mu)
      real(dp), intent(out) :: x(:); real(dp), intent(in) :: lambda,mu
      integer :: i
      do i=1,size(x); x(i)=qtprd(rand_uniform(),lambda,mu); end do
   end subroutine rtprd

   ! Uniform-geometric -------------------------------------------------------
   real(dp) function dugd(x,theta,log_p) result(v)
      real(dp), intent(in) :: x,theta
      logical, intent(in), optional :: log_p
      logical :: lp
      integer :: k
      real(dp) :: pv
      lp=.false.; if(present(log_p)) lp=log_p
      if(theta<=0.0_dp .or. theta>=1.0_dp) then; v=nan_dp(); return; end if
      k=floor(x)
      if(k<=0) then; pv=0.0_dp
      else; pv=theta*(1.0_dp-theta)**(k-1)*lerch_phi(1.0_dp-theta,1.0_dp,real(k,dp)); end if
      if(lp) then; if(pv<=0.0_dp) then; v=-huge(1.0_dp); else; v=log(pv); end if
      else; v=pv; end if
   end function dugd

   real(dp) function pugd(q,theta,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,theta
      logical, intent(in), optional :: lower_tail,log_p
      integer :: k
      real(dp) :: cdf
      if(theta<=0.0_dp .or. theta>=1.0_dp) then; v=nan_dp(); return; end if
      k=floor(q)
      if(k<1) then; cdf=0.0_dp
      else
         cdf=1.0_dp-theta*(1.0_dp-theta)**k*(1.0_dp/theta-real(k,dp)* &
            lerch_phi(1.0_dp-theta,1.0_dp,real(k+1,dp)))
      end if
      v=finish_probability(cdf,lower_tail,log_p)
   end function pugd

   integer function qugd(p,theta,lower_tail) result(k)
      real(dp), intent(in) :: p,theta
      logical, intent(in), optional :: lower_tail
      real(dp) :: t
      integer :: lo,hi,mid
      if(p<0.0_dp .or. p>1.0_dp .or. theta<=0.0_dp .or. theta>=1.0_dp) then; k=-huge(1); return; end if
      t=target_probability(p,lower_tail)
      if(t<=0.0_dp) then; k=1; return; end if
      if(t>=1.0_dp) then; k=huge(1); return; end if
      lo=1; hi=2
      do while(pugd(real(hi,dp),theta)<t); hi=2*hi; if(hi>500000000) exit; end do
      do while(lo<hi)
         mid=lo+(hi-lo)/2
         if(pugd(real(mid,dp),theta)>=t) then; hi=mid; else; lo=mid+1; end if
      end do
      k=lo
   end function qugd

   subroutine rugd(x,theta)
      integer, intent(out) :: x(:); real(dp), intent(in) :: theta
      integer :: i
      do i=1,size(x); x(i)=qugd(rand_uniform(),theta); end do
   end subroutine rugd

   ! Inverse Gaussian (named Unit IG upstream) -------------------------------
   pure real(dp) function normal_logcdf(z) result(lv)
      real(dp), intent(in) :: z
      real(dp) :: zz,inv2,series
      if(z>-10.0_dp) then
         lv=log(max(normal_cdf(z),tiny(1.0_dp)))
      else
         zz=-z; inv2=1.0_dp/(zz*zz)
         series=1.0_dp-inv2+3.0_dp*inv2*inv2-15.0_dp*inv2**3
         lv=-0.5_dp*zz*zz-log(zz)-0.5_dp*log(2.0_dp*pi)+log(series)
      end if
   end function normal_logcdf

   real(dp) function duigd(x,mu,lambda,log_p) result(v)
      real(dp), intent(in) :: x,mu,lambda
      logical, intent(in), optional :: log_p
      logical :: lp
      real(dp) :: lv
      lp=.false.; if(present(log_p)) lp=log_p
      if(mu<=0.0_dp .or. lambda<=0.0_dp) then; v=nan_dp(); return; end if
      if(x<=0.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lp); return; end if
      lv=0.5_dp*(log(lambda)-log(2.0_dp*pi))-1.5_dp*log(x)- &
         lambda*(x-mu)**2/(2.0_dp*mu*mu*x)
      v=merge(lv,exp(lv),lp)
   end function duigd

   real(dp) function puigd(q,mu,lambda,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,mu,lambda
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: cdf,z1,z2,l2
      if(mu<=0.0_dp .or. lambda<=0.0_dp) then; v=nan_dp(); return; end if
      if(q<=0.0_dp) then; cdf=0.0_dp
      else
         z1=sqrt(lambda/q)*(q/mu-1.0_dp)
         z2=-sqrt(lambda/q)*(q/mu+1.0_dp)
         l2=2.0_dp*lambda/mu+normal_logcdf(z2)
         cdf=normal_cdf(z1)+exp(l2)
      end if
      v=finish_probability(cdf,lower_tail,log_p)
   end function puigd

   real(dp) function quigd(p,mu,lambda,lower_tail) result(x)
      real(dp), intent(in) :: p,mu,lambda
      logical, intent(in), optional :: lower_tail
      real(dp) :: t,lo,hi,mid
      integer :: it
      if(p<0.0_dp .or. p>1.0_dp .or. mu<=0.0_dp .or. lambda<=0.0_dp) then; x=nan_dp(); return; end if
      t=target_probability(p,lower_tail)
      if(t<=0.0_dp) then; x=0.0_dp; return; end if
      if(t>=1.0_dp) then; x=huge(1.0_dp); return; end if
      lo=0.0_dp; hi=max(1.0_dp,mu)
      do while(puigd(hi,mu,lambda)<t); hi=2.0_dp*hi; if(hi>1.0e12_dp) exit; end do
      do it=1,120
         mid=0.5_dp*(lo+hi)
         if(puigd(mid,mu,lambda)<t) then; lo=mid; else; hi=mid; end if
      end do
      x=0.5_dp*(lo+hi)
   end function quigd

   subroutine ruigd(x,mu,lambda)
      real(dp), intent(out) :: x(:); real(dp), intent(in) :: mu,lambda
      integer :: i
      do i=1,size(x); x(i)=quigd(rand_uniform(),mu,lambda); end do
   end subroutine ruigd

   ! Weighted geometric ------------------------------------------------------
   real(dp) function dwgd(x,alpha,lambda,log_p) result(v)
      real(dp), intent(in) :: x,alpha,lambda
      logical, intent(in), optional :: log_p
      logical :: lp
      integer :: k
      real(dp) :: pv
      lp=.false.; if(present(log_p)) lp=log_p
      if(alpha<=0.0_dp .or. alpha>=1.0_dp .or. lambda<=0.0_dp) then; v=nan_dp(); return; end if
      k=floor(x)
      if(k<=0) then; pv=0.0_dp
      else
         pv=(1.0_dp-alpha)*(1.0_dp-alpha**(lambda+1.0_dp))/(1.0_dp-alpha**lambda)* &
            alpha**(k-1)*(1.0_dp-alpha**(lambda*real(k,dp)))
      end if
      if(lp) then; if(pv<=0.0_dp) then; v=-huge(1.0_dp); else; v=log(pv); end if
      else; v=pv; end if
   end function dwgd

   real(dp) function pwgd(q,alpha,lambda,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,alpha,lambda
      logical, intent(in), optional :: lower_tail,log_p
      integer :: k
      real(dp) :: cdf
      if(alpha<=0.0_dp .or. alpha>=1.0_dp .or. lambda<=0.0_dp) then; v=nan_dp(); return; end if
      k=floor(q)
      if(k<=0) then; cdf=0.0_dp
      else
         cdf=1.0_dp-((1.0_dp-alpha**(lambda+1.0_dp)-alpha**(lambda*real(k+1,dp))* &
            (1.0_dp-alpha))/(1.0_dp-alpha**lambda))*alpha**k
      end if
      v=finish_probability(cdf,lower_tail,log_p)
   end function pwgd

   integer function qwgd(p,alpha,lambda,lower_tail) result(k)
      real(dp), intent(in) :: p,alpha,lambda
      logical, intent(in), optional :: lower_tail
      real(dp) :: t
      integer :: lo,hi,mid
      if(p<0.0_dp .or. p>1.0_dp .or. alpha<=0.0_dp .or. alpha>=1.0_dp .or. lambda<=0.0_dp) then
         k=-huge(1); return
      end if
      t=target_probability(p,lower_tail)
      if(t<=0.0_dp) then; k=1; return; end if
      if(t>=1.0_dp) then; k=huge(1); return; end if
      lo=1; hi=2
      do while(pwgd(real(hi,dp),alpha,lambda)<t); hi=2*hi; if(hi>500000000) exit; end do
      do while(lo<hi)
         mid=lo+(hi-lo)/2
         if(pwgd(real(mid,dp),alpha,lambda)>=t) then; hi=mid; else; lo=mid+1; end if
      end do
      k=lo
   end function qwgd

   subroutine rwgd(x,alpha,lambda)
      integer, intent(out) :: x(:); real(dp), intent(in) :: alpha,lambda
      integer :: i
      do i=1,size(x); x(i)=qwgd(rand_uniform(),alpha,lambda); end do
   end subroutine rwgd

end module new_dist_distributions
