module forecast_spline
   use forecast_kinds, only : dp
   use forecast_types, only : spline_model_t, forecast_result
   use forecast_linalg, only : inverse_matrix, logdet_spd
   use forecast_optimize, only : golden_minimize
   use forecast_stats, only : normal_quantile
   implicit none
   private
   public :: spline_model_fit, spline_forecast, spline_sigma_matrix, spline_gcv
   real(dp), allocatable, save :: objective_y(:)
contains
   function spline_sigma_matrix(n, n0) result(sigma)
      integer, intent(in) :: n
      integer, intent(in), optional :: n0
      real(dp), allocatable :: sigma(:,:)
      integer :: nn, extra, i, j
      extra = 0
      if (present(n0)) extra = max(0,n0)
      nn = n + extra
      allocate(sigma(nn,nn))
      do i = 1, nn
         do j = 1, nn
            sigma(i,j) = real(min(i,j),dp)**2 * &
               real(3*max(i,j)-min(i,j),dp)/6.0_dp/real(n,dp)**3
         end do
      end do
   end function spline_sigma_matrix

   subroutine spline_omega(n,beta,n0,omega)
      integer,intent(in)::n,n0
      real(dp),intent(in)::beta
      real(dp),allocatable,intent(out)::omega(:,:)
      real(dp),allocatable::sigma(:,:),s(:,:)
      integer::nn,i
      nn=n+n0
      sigma=spline_sigma_matrix(n,n0)
      allocate(s(nn,2),omega(nn,nn))
      s(:,1)=1.0_dp
      do i=1,nn
      s(i,2)=real(i,dp)/real(n,dp)
      end do
      omega=100.0_dp*matmul(s,transpose(s))+sigma/max(beta,1.0e-15_dp)
      do i=1,nn
      omega(i,i)=omega(i,i)+1.0_dp
      end do
   end subroutine spline_omega

   function spline_objective(logbeta) result(value)
      real(dp),intent(in)::logbeta
      real(dp)::value,ld,q
      real(dp),allocatable::omega(:,:),inv(:,:)
      integer::info,n
      n=size(objective_y)
      call spline_omega(n,exp(logbeta),0,omega)
      call logdet_spd(omega,ld,info)
      if(info/=0)then
      value=huge(1.0_dp)
      return
      end if
      call inverse_matrix(omega,inv,info)
      if(info/=0)then
      value=huge(1.0_dp)
      return
      end if
      q=dot_product(objective_y,matmul(inv,objective_y))
      value=0.5_dp*ld+0.5_dp*real(n,dp)*log(max(q/real(n,dp),1.0e-300_dp))
   end function spline_objective


   subroutine spline_universal_smoother(y,beta,fitted,resid,trace_residual,info)
      real(dp),intent(in)::y(:),beta
      real(dp),allocatable,intent(out)::fitted(:),resid(:)
      real(dp),intent(out)::trace_residual
      integer,intent(out)::info
      real(dp),allocatable::sigma(:,:),m(:,:),minv(:,:),x(:,:),mx(:,:),a(:,:),ainv(:,:),proj(:,:)
      integer::n,i
      n=size(y)
      info=0
      sigma=spline_sigma_matrix(n)
      allocate(m(n,n),x(n,2))
      m=sigma/max(beta,1.0e-15_dp)
      do i=1,n
         m(i,i)=m(i,i)+1.0_dp
         x(i,1)=1.0_dp
         x(i,2)=real(i,dp)/real(n,dp)
      end do
      call inverse_matrix(m,minv,info)
      if(info/=0)return
      mx=matmul(minv,x)
      a=matmul(transpose(x),mx)
      call inverse_matrix(a,ainv,info)
      if(info/=0)return
      proj=minv-matmul(matmul(mx,ainv),transpose(mx))
      resid=matmul(proj,y)
      fitted=y-resid
      trace_residual=0.0_dp
      do i=1,n
         trace_residual=trace_residual+proj(i,i)
      end do
   end subroutine spline_universal_smoother

   function spline_gcv(logbeta) result(value)
      real(dp),intent(in)::logbeta
      real(dp)::value,trp,rss
      real(dp),allocatable::fitted(:),resid(:)
      integer::info,n
      n=size(objective_y)
      call spline_universal_smoother(objective_y,exp(logbeta),fitted,resid,trp,info)
      if(info/=0.or.trp<=tiny(1.0_dp))then
         value=huge(1.0_dp)
         return
      end if
      rss=sum(resid**2)
      value=(rss/real(n,dp))/max((trp/real(n,dp))**2,tiny(1.0_dp))
   end function spline_gcv

   function spline_model_fit(y,beta,method) result(model)
      real(dp),intent(in)::y(:)
      real(dp),intent(in),optional::beta
      character(len=*),intent(in),optional::method
      type(spline_model_t)::model
      real(dp),allocatable::omega(:,:),inv(:,:),oi(:,:),u(:),smooth_resid(:)
      real(dp)::b,sd,trp
      integer::n,i,info
      character(len=8)::meth
      n=size(y)
      if(n<4)error stop 'spline_model_fit: at least four observations required'
      meth='gcv'
      if(present(method))meth=adjustl(method)
      if(present(beta))then
         b=max(beta/real(n,dp)**3,1.0e-12_dp)
      else
         if(trim(meth)=='mle' .or. trim(meth)=='MLE')then
            objective_y=y(max(1,n-99):n)
            b=exp(golden_minimize(spline_objective,log(1.0e-12_dp),log(10.0_dp),tol=1.0e-5_dp))
         else
            objective_y=y
            b=exp(golden_minimize(spline_gcv,log(1.0e-12_dp),log(10.0_dp),tol=1.0e-5_dp))
         end if
      end if
      model%x=[(real(i,dp),i=1,n)]
      model%y=y
      model%beta=b*real(n,dp)**3
      call spline_omega(n,b,0,omega)
      call inverse_matrix(omega,inv,info)
      if(info/=0)error stop 'spline_model_fit: singular covariance matrix'
      call spline_universal_smoother(y,b,model%fitted,smooth_resid,trp,info)
      if(info/=0)error stop 'spline_model_fit: universal smoother failed'
      allocate(model%onestep(n),model%residuals(n))
      model%onestep=0.0_dp
      model%residuals=0.0_dp
      do i=2,n
         call spline_omega(i-1,b,1,oi)
         u=oi(1:i-1,i)
         call inverse_matrix(oi(1:i-1,1:i-1),inv,info)
         if(info/=0)cycle
         model%onestep(i)=dot_product(u,matmul(inv,y(1:i-1)))
         sd=sqrt(max(oi(i,i)-dot_product(u,matmul(inv,u)),1.0e-14_dp))
         model%residuals(i)=(y(i)-model%onestep(i))/sd
      end do
      model%sigma2=sum(model%residuals(2:)**2)/real(n-1,dp)
   end function spline_model_fit

   function spline_forecast(model,h,levels) result(fc)
      type(spline_model_t),intent(in)::model
      integer,intent(in)::h
      real(dp),intent(in),optional::levels(:)
      type(forecast_result)::fc
      real(dp),allocatable::omega(:,:),newomega(:,:),inv(:,:),u(:,:),covf(:,:)
      real(dp)::b,z
      integer::n,info,i,j
      n=size(model%y)
      b=model%beta/real(n,dp)**3
      call spline_omega(n,b,0,omega)
      call inverse_matrix(omega,inv,info)
      if(info/=0)error stop 'spline_forecast: matrix inversion failed'
      call spline_omega(n,b,h,newomega)
      u=newomega(1:n,n+1:n+h)
      allocate(fc%mean(h),fc%se(h))
      fc%mean=matmul(transpose(u),matmul(inv,model%y))
      covf=newomega(n+1:n+h,n+1:n+h)-matmul(transpose(u),matmul(inv,u))
      do i=1,h
      fc%se(i)=sqrt(max(model%sigma2*covf(i,i),0.0_dp))
      end do
      if(present(levels))then
         fc%level=levels
         allocate(fc%lower(h,size(levels)),fc%upper(h,size(levels)))
         do j=1,size(levels)
         z=normal_quantile(0.5_dp+0.005_dp*levels(j))
         fc%lower(:,j)=fc%mean-z*fc%se
         fc%upper(:,j)=fc%mean+z*fc%se
         end do
      end if
   end function spline_forecast
end module forecast_spline
