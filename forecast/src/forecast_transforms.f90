module forecast_transforms
   use forecast_kinds, only : dp
   use forecast_stats, only : mean_value, variance_value
   use forecast_optimize, only : golden_minimize
   implicit none
   private
   public :: boxcox, inv_boxcox, boxcox_lambda_guerrero, boxcox_lambda_loglik
   real(dp), allocatable, save :: lambda_x(:)
   integer, save :: lambda_m=1
contains
   pure function boxcox(x,lambda) result(y)
      real(dp),intent(in) :: x(:),lambda
      real(dp),allocatable :: y(:)
      allocate(y(size(x)))
      if(abs(lambda)<1.0e-12_dp) then
      y=log(x)
      else
      y=(x**lambda-1.0_dp)/lambda
      end if
   end function
   pure function inv_boxcox(x,lambda,biasadj,fvar) result(y)
      real(dp),intent(in) :: x(:),lambda
      logical,intent(in),optional :: biasadj
      real(dp),intent(in),optional :: fvar(:)
      real(dp),allocatable :: y(:),fv(:)
      logical :: ba
      allocate(y(size(x)))
      ba=.false.
      if(present(biasadj)) ba=biasadj
      if(abs(lambda)<1.0e-12_dp) then
         y=exp(x)
      else
         y=max(0.0_dp,lambda*x+1.0_dp)**(1.0_dp/lambda)
      end if
      if(ba .and. present(fvar)) then
         allocate(fv(size(x)))
         fv=fvar
         if(abs(lambda)<1.0e-12_dp) then
            y=y*(1.0_dp+0.5_dp*fv)
         else
            y=y*(1.0_dp+0.5_dp*fv*(1.0_dp-lambda)/(max(1.0e-12_dp,lambda*x+1.0_dp)**2))
         end if
      end if
   end function
   function boxcox_lambda_guerrero(x,m,lower,upper) result(lambda)
      real(dp),intent(in) :: x(:)
      integer,intent(in) :: m
      real(dp),intent(in),optional :: lower,upper
      real(dp) :: lambda,l,u
      lambda_x=x
      lambda_m=max(1,m)
      l=-1.0_dp
      u=2.0_dp
      if(present(lower)) l=lower
      if(present(upper)) u=upper
      lambda=golden_minimize(guerrero_objective,l,u,1.0e-5_dp)
   end function
   function guerrero_objective(lambda) result(v)
      real(dp),intent(in)::lambda
      real(dp)::v
      integer :: n,nj,j,k,start
      real(dp),allocatable :: ratios(:),z(:)
      real(dp)::mu,sd
      n=size(lambda_x)
      nj=n/lambda_m
      if(nj<2) then
      v=huge(1.0_dp)
      return
      end if
      allocate(ratios(nj))
      start=n-nj*lambda_m+1
      do j=1,nj
         k=start+(j-1)*lambda_m
         z=lambda_x(k:k+lambda_m-1)
         mu=mean_value(z)
         if(lambda_m>1) then
         sd=sqrt(variance_value(z,.true.))
         else
         sd=0.0_dp
         end if
         ratios(j)=sd/max(abs(mu)**(1.0_dp-lambda),1.0e-12_dp)
      end do
      v=sqrt(variance_value(ratios,.true.))/max(mean_value(ratios),1.0e-12_dp)
   end function
   function boxcox_lambda_loglik(x,lower,upper) result(lambda)
      real(dp),intent(in) :: x(:)
      real(dp),intent(in),optional :: lower,upper
      real(dp)::lambda,l,u
      lambda_x=x
      l=-1.0_dp
      u=2.0_dp
      if(present(lower))l=lower
      if(present(upper))u=upper
      lambda=golden_minimize(loglik_objective,l,u,1.0e-5_dp)
   end function
   function loglik_objective(lambda) result(v)
      real(dp),intent(in)::lambda
      real(dp)::v
      real(dp),allocatable::z(:)
      real(dp)::gm
      gm=exp(sum(log(lambda_x))/real(size(lambda_x),dp))
      z=boxcox(lambda_x/gm,lambda)
      v=real(size(z),dp)*log(max(sum((z-mean_value(z))**2)/real(size(z),dp),tiny(1.0_dp)))
   end function
end module forecast_transforms
