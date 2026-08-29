! SPDX-License-Identifier: GPL-3.0-only
module matrixdist_iph
   use r_compat, only: dp
   use matrixdist_linalg, only: matrix_exponential
   use matrixdist_ph, only: ph_exit_rates
   implicit none
   private
   public :: iph_density, iph_cdf, iph_inverse_transform, transform_time_jac
   public :: mweibull_density, mweibull_cdf, mpareto_density, mpareto_cdf
   public :: mlognormal_density, mlognormal_cdf, mloglogistic_density, mloglogistic_cdf
   public :: mgompertz_density, mgompertz_cdf, mgev_density, mgev_cdf

contains

   subroutine transform_time_jac(x, kind, beta, t, jac, valid)
      real(dp), intent(in) :: x, beta(:)
      character(len=*), intent(in) :: kind
      real(dp), intent(out) :: t, jac
      logical, intent(out) :: valid
      real(dp) :: z
      valid=.true.
      t=0.0_dp
      jac=0.0_dp
      select case(trim(kind))
      case('identity')
         if(x<0.0_dp) then
         valid=.false.
         return
         end if
         t=x
         jac=1.0_dp
      case('weibull')
         if(x<0.0_dp .or. beta(1)<=0.0_dp) then
         valid=.false.
         return
         end if
         t=x**beta(1)
         if(x==0.0_dp .and. beta(1)<1.0_dp) then
            jac=huge(1.0_dp)
         else if(x==0.0_dp .and. beta(1)>1.0_dp) then
            jac=0.0_dp
         else
            jac=beta(1)*x**(beta(1)-1.0_dp)
         end if
      case('pareto')
         if(x<0.0_dp .or. beta(1)<=0.0_dp) then
         valid=.false.
         return
         end if
         t=log(1.0_dp+x/beta(1))
         jac=1.0_dp/(x+beta(1))
      case('lognormal')
         if(x<0.0_dp .or. beta(1)<=0.0_dp) then
         valid=.false.
         return
         end if
         z=log(1.0_dp+x)
         t=z**beta(1)
         if(z==0.0_dp .and. beta(1)<1.0_dp) then
            jac=huge(1.0_dp)
         else if(z==0.0_dp .and. beta(1)>1.0_dp) then
            jac=0.0_dp
         else
            jac=beta(1)*z**(beta(1)-1.0_dp)/(x+1.0_dp)
         end if
      case('loglogistic')
         if(x<0.0_dp .or. size(beta)<2 .or. beta(1)<=0.0_dp .or. beta(2)<=0.0_dp) then
            valid=.false.
            return
         end if
         z=(x/beta(1))**beta(2)
         t=log(1.0_dp+z)
         if(x==0.0_dp .and. beta(2)<1.0_dp) then
            jac=huge(1.0_dp)
         else if(x==0.0_dp .and. beta(2)>1.0_dp) then
            jac=0.0_dp
         else
            jac=(x/beta(1))**(beta(2)-1.0_dp)*beta(2)/beta(1)/(1.0_dp+z)
         end if
      case('gompertz')
         if(x<0.0_dp .or. beta(1)==0.0_dp) then
            if(x<0.0_dp) then
            valid=.false.
            return
            end if
            t=x
            jac=1.0_dp
         else
            t=(exp(beta(1)*x)-1.0_dp)/beta(1)
            jac=exp(beta(1)*x)
         end if
      case('gev')
         if(size(beta)<3 .or. beta(2)<=0.0_dp) then
         valid=.false.
         return
         end if
         if(beta(3)==0.0_dp) then
            t=exp(-(x-beta(1))/beta(2))
            jac=t/beta(2)
         else
            z=1.0_dp+(beta(3)/beta(2))*(x-beta(1))
            if(z<=0.0_dp) then
            valid=.false.
            return
            end if
            t=z**(-1.0_dp/beta(3))
            jac=z**(-(1.0_dp+beta(3))/beta(3))/beta(2)
         end if
      case default
         valid=.false.
      end select
   end subroutine transform_time_jac

   function iph_density(x,alpha,s,kind,beta,scale) result(f)
      real(dp),intent(in)::x,alpha(:),s(:,:),beta(:)
      character(len=*),intent(in)::kind
      real(dp),intent(in),optional::scale
      real(dp)::f,t,jac,sc,xx
      logical::ok
      real(dp),allocatable::m(:,:),exitv(:)
      sc=1.0_dp
      if(present(scale))sc=scale
      if(sc<=0.0_dp)error stop 'iph_density: scale must be positive'
      xx=x/sc
      if(trim(kind)/='gev' .and. xx==0.0_dp) then
         f=max(0.0_dp,1.0_dp-sum(alpha))
         return
      end if
      call transform_time_jac(xx,kind,beta,t,jac,ok)
      if(.not.ok) then
      f=0.0_dp
      return
      end if
      exitv=ph_exit_rates(s)
      m=matrix_exponential(s*t)
      f=dot_product(alpha,matmul(m,exitv))*jac/sc
   end function iph_density

   function iph_cdf(x,alpha,s,kind,beta,lower_tail,scale) result(pv)
      real(dp),intent(in)::x,alpha(:),s(:,:),beta(:)
      character(len=*),intent(in)::kind
      logical,intent(in),optional::lower_tail
      real(dp),intent(in),optional::scale
      real(dp)::pv,t,jac,sv,sc,xx
      logical::ok,lower
      real(dp),allocatable::m(:,:),e(:)
      lower=.true.
      if(present(lower_tail))lower=lower_tail
      sc=1.0_dp
      if(present(scale))sc=scale
      if(sc<=0.0_dp)error stop 'iph_cdf: scale must be positive'
      xx=x/sc
      if(trim(kind)/='gev' .and. xx<0.0_dp) then
         pv=merge(0.0_dp,1.0_dp,lower)
         return
      end if
      call transform_time_jac(xx,kind,beta,t,jac,ok)
      if(.not.ok) then
         if(trim(kind)=='gev' .and. size(beta)>=3) then
            ! The GEV time transformation is decreasing.  Below a finite
            ! lower endpoint (xi>0) the CDF is zero; above a finite upper
            ! endpoint (xi<0) it is one.
            if(beta(3)>0.0_dp .and. xx<=beta(1)-beta(2)/beta(3)) then
               pv=merge(0.0_dp,1.0_dp,lower)
            else if(beta(3)<0.0_dp .and. xx>=beta(1)-beta(2)/beta(3)) then
               pv=merge(1.0_dp,0.0_dp,lower)
            else
               pv=merge(0.0_dp,1.0_dp,lower)
            end if
         else
            pv=merge(0.0_dp,1.0_dp,lower)
         end if
         return
      end if
      allocate(e(size(alpha)))
      e=1.0_dp
      m=matrix_exponential(s*t)
      sv=dot_product(alpha,matmul(m,e))
      if(trim(kind)=='gev') then
         ! Since t=g^{-1}(x) decreases with x, F_X(x)=P(T>=t).
         ! This corrects the sign/orientation inconsistency in upstream
         ! mgevcdf(), whose returned value is incompatible with mgevden()
         ! and rmatrixgev() even in the one-state GEV special case.
         pv=sv
      else
         pv=1.0_dp-sv
      end if
      if(.not.lower)pv=1.0_dp-pv
      pv=min(1.0_dp,max(0.0_dp,pv))
   end function iph_cdf

   function iph_inverse_transform(t,kind,beta) result(x)
      real(dp),intent(in)::t,beta(:)
      character(len=*),intent(in)::kind
      real(dp)::x
      select case(trim(kind))
      case('pareto'); x=beta(1)*(exp(t)-1.0_dp)
      case('weibull'); x=t**(1.0_dp/beta(1))
      case('lognormal'); x=exp(t**(1.0_dp/beta(1)))-1.0_dp
      case('loglogistic'); x=beta(1)*(exp(t)-1.0_dp)**(1.0_dp/beta(2))
      case('gompertz')
         if(beta(1)==0.0_dp) then
         x=t
         else
         x=log(1.0_dp+beta(1)*t)/beta(1)
         end if
      case('gev')
         if(beta(3)==0.0_dp) then
            x=beta(1)-beta(2)*log(t)
         else
            x=beta(1)+beta(2)*(t**(-beta(3))-1.0_dp)/beta(3)
         end if
      case default; x=t
      end select
   end function iph_inverse_transform

   function mweibull_density(x,a,s,beta) result(v)
      real(dp),intent(in)::x,a(:),s(:,:),beta
      real(dp)::v
      v=iph_density(x,a,s,'weibull',[beta])
   end function
   function mweibull_cdf(x,a,s,beta,lower_tail) result(v)
      real(dp),intent(in)::x,a(:),s(:,:),beta
      logical,intent(in),optional::lower_tail
      real(dp)::v
      if(present(lower_tail))then
      v=iph_cdf(x,a,s,'weibull',[beta],lower_tail)
      else
      v=iph_cdf(x,a,s,'weibull',[beta])
      end if
   end function
   function mpareto_density(x,a,s,beta) result(v)
      real(dp),intent(in)::x,a(:),s(:,:),beta
      real(dp)::v
      v=iph_density(x,a,s,'pareto',[beta])
   end function
   function mpareto_cdf(x,a,s,beta,lower_tail) result(v)
      real(dp),intent(in)::x,a(:),s(:,:),beta
      logical,intent(in),optional::lower_tail
      real(dp)::v
      if(present(lower_tail))then
      v=iph_cdf(x,a,s,'pareto',[beta],lower_tail)
      else
      v=iph_cdf(x,a,s,'pareto',[beta])
      end if
   end function
   function mlognormal_density(x,a,s,beta) result(v)
      real(dp),intent(in)::x,a(:),s(:,:),beta
      real(dp)::v
      v=iph_density(x,a,s,'lognormal',[beta])
   end function
   function mlognormal_cdf(x,a,s,beta,lower_tail) result(v)
      real(dp),intent(in)::x,a(:),s(:,:),beta
      logical,intent(in),optional::lower_tail
      real(dp)::v
      if(present(lower_tail))then
      v=iph_cdf(x,a,s,'lognormal',[beta],lower_tail)
      else
      v=iph_cdf(x,a,s,'lognormal',[beta])
      end if
   end function
   function mloglogistic_density(x,a,s,beta) result(v)
      real(dp),intent(in)::x,a(:),s(:,:),beta(:)
      real(dp)::v
      v=iph_density(x,a,s,'loglogistic',beta)
   end function
   function mloglogistic_cdf(x,a,s,beta,lower_tail) result(v)
      real(dp),intent(in)::x,a(:),s(:,:),beta(:)
      logical,intent(in),optional::lower_tail
      real(dp)::v
      if(present(lower_tail))then
      v=iph_cdf(x,a,s,'loglogistic',beta,lower_tail)
      else
      v=iph_cdf(x,a,s,'loglogistic',beta)
      end if
   end function
   function mgompertz_density(x,a,s,beta) result(v)
      real(dp),intent(in)::x,a(:),s(:,:),beta
      real(dp)::v
      v=iph_density(x,a,s,'gompertz',[beta])
   end function
   function mgompertz_cdf(x,a,s,beta,lower_tail) result(v)
      real(dp),intent(in)::x,a(:),s(:,:),beta
      logical,intent(in),optional::lower_tail
      real(dp)::v
      if(present(lower_tail))then
      v=iph_cdf(x,a,s,'gompertz',[beta],lower_tail)
      else
      v=iph_cdf(x,a,s,'gompertz',[beta])
      end if
   end function
   function mgev_density(x,a,s,beta) result(v)
      real(dp),intent(in)::x,a(:),s(:,:),beta(:)
      real(dp)::v
      v=iph_density(x,a,s,'gev',beta)
   end function
   function mgev_cdf(x,a,s,beta,lower_tail) result(v)
      real(dp),intent(in)::x,a(:),s(:,:),beta(:)
      logical,intent(in),optional::lower_tail
      real(dp)::v
      if(present(lower_tail))then
      v=iph_cdf(x,a,s,'gev',beta,lower_tail)
      else
      v=iph_cdf(x,a,s,'gev',beta)
      end if
   end function

end module matrixdist_iph
