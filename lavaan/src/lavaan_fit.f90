module lavaan_fit
   use lavaan_kinds, only : dp
   use lavaan_ram, only : ram_model, ram_free_map, ram_sigma, ram_mu, ram_set_free, ram_get_free
   use lavaan_linalg, only : sample_mean_cov, vech, inverse_general, logdet_spd
   use lavaan_objectives, only : objective_ml, objective_gls, objective_uls, objective_wls, objective_dwls
   use lavaan_objectives, only : mvn_loglik_complete, mvn_loglik_missing
   use lavaan_optimizer, only : bfgs_minimize
   use numderiv, only : hessian, nd_success
   implicit none
   private
   type, public :: sem_fit_result
      real(dp), allocatable :: par(:), se(:), vcov(:, :), sigma(:, :), mu(:)
      real(dp) :: objective=huge(1.0_dp), loglik=-huge(1.0_dp)
      real(dp) :: chisq=huge(1.0_dp), df=0.0_dp, aic=huge(1.0_dp), bic=huge(1.0_dp)
      real(dp) :: rmsea=huge(1.0_dp), cfi=0.0_dp, tli=0.0_dp, srmr=huge(1.0_dp)
      logical :: converged=.false.
      integer :: iterations=0, status=0
   end type sem_fit_result
   public :: fit_ram_cov, fit_ram_data, fit_ram_fiml, standardized_ram, residual_covariance
contains
   subroutine fit_ram_cov(template,map,data_cov,data_mean,nobs,result,estimator,wls_v,wls_vd)
      type(ram_model),intent(in)::template
      type(ram_free_map),intent(in)::map
      real(dp),intent(in)::data_cov(:,:),data_mean(:)
      integer,intent(in)::nobs
      type(sem_fit_result),intent(out)::result
      character(len=*),intent(in),optional::estimator
      real(dp),intent(in),optional::wls_v(:,:),wls_vd(:)
      character(len=8)::est
      real(dp),allocatable::x(:),hess(:,:),hi(:,:),diagcov(:,:)
      real(dp)::fval,fbase,chi_base,df_base,ld,q
      integer::info,p,k,status,i
      logical::meanstructure
      type(ram_model)::work
      est='ML'
      if(present(estimator)) est=adjustl(estimator)
      work=template
      meanstructure=allocated(template%m)
      x=ram_get_free(template,map)
      k=size(x)
      p=size(data_cov,1)
      call bfgs_minimize(criterion,x,fval,result%converged,result%iterations,maxiter=1000,tol=1e-7_dp)
      work=template
      call ram_set_free(work,map,x)
      call ram_sigma(work,result%sigma,info)
      call ram_mu(work,result%mu,info)
      result%par=x
      result%objective=fval
      result%status=info
      result%df=real(p*(p+1)/2 + merge(p,0,meanstructure) - k,dp)
      if(trim(est)=='ML') then
         result%chisq=real(nobs,dp)*fval
         ld=logdet_spd(result%sigma,info)
         if(info==0) then
            q=sum(data_cov*spd_inverse(result%sigma,info))
            if(meanstructure) q=q+dot_product(data_mean-result%mu, &
               matmul(spd_inverse(result%sigma,info),data_mean-result%mu))
            result%loglik=-0.5_dp*real(nobs,dp)*(real(p,dp)*log(2*acos(-1.0_dp))+ld+q)
         end if
      else
         result%chisq=real(nobs,dp)*fval
         result%loglik=-0.5_dp*result%chisq
      end if
      result%aic=-2*result%loglik+2*real(k,dp)
      result%bic=-2*result%loglik+log(real(nobs,dp))*real(k,dp)
      allocate(diagcov(p,p))
      diagcov=0
      do i=1,p
      diagcov(i,i)=data_cov(i,i)
      end do
      fbase=objective_ml(diagcov,data_mean,data_cov,data_mean,meanstructure,info)
      chi_base=real(nobs,dp)*fbase
      df_base=real(p*(p-1)/2,dp)
      if(result%df>0) result%rmsea=sqrt(max((result%chisq-result%df)/(result%df*real(nobs-1,dp)),0.0_dp))
      result%cfi=1.0_dp-max(result%chisq-result%df,0.0_dp)/ &
         max(max(chi_base-df_base,result%chisq-result%df),1.0e-12_dp)
      if(result%df>0 .and. df_base>0 .and. chi_base>1e-12_dp) then
         result%tli=(chi_base/df_base-result%chisq/result%df)/(chi_base/df_base-1.0_dp)
      end if
      result%srmr=compute_srmr(data_cov,result%sigma)
      call hessian(total_criterion,x,hess,status=status)
      allocate(result%vcov(k,k),result%se(k))
      result%vcov=0
      result%se=huge(1.0_dp)
      if(status==nd_success) then
         call inverse_general(hess,hi,info)
         if(info==0) then
            result%vcov=hi
            do i=1,k
            if(hi(i,i)>=0) result%se(i)=sqrt(hi(i,i))
            end do
         end if
      end if
   contains
      function criterion(z) result(v)
         real(dp),intent(in)::z(:)
         real(dp)::v
         real(dp),allocatable::sg(:,:),mm(:),ev(:),ov(:)
         integer::istat
         work=template
         call ram_set_free(work,map,z)
         call ram_sigma(work,sg,istat)
         if(istat/=0 .or. any([(sg(i,i)<=0.0_dp,i=1,size(sg,1))])) then
         v=huge(1.0_dp)/100
         return
         end if
         call ram_mu(work,mm,istat)
         select case(trim(est))
         case('ML'); v=objective_ml(sg,mm,data_cov,data_mean,meanstructure,istat)
         case('GLS'); v=objective_gls(sg,mm,data_cov,data_mean,meanstructure,istat)
         case default
            ev=vech(sg)
            ov=vech(data_cov)
            if(meanstructure) then
            ev=[ev,mm]
            ov=[ov,data_mean]
            end if
            select case(trim(est))
            case('ULS'); v=objective_uls(ev,ov)
            case('WLS')
            if(present(wls_v)) then
            v=objective_wls(ev,ov,wls_v)
            else
            v=objective_uls(ev,ov)
            end if
            case('DWLS')
            if(present(wls_vd)) then
            v=objective_dwls(ev,ov,wls_vd)
            else
            v=objective_uls(ev,ov)
            end if
            case default; v=huge(1.0_dp)/100
            end select
         end select
         if(.not.(v<huge(1.0_dp)/10)) v=huge(1.0_dp)/100
      end function criterion
      function total_criterion(z) result(v)
         real(dp),intent(in)::z(:)
         real(dp)::v
         v=0.5_dp*real(nobs,dp)*criterion(z)
      end function total_criterion
   end subroutine fit_ram_cov

   subroutine fit_ram_data(template,map,data,result,estimator)
      type(ram_model),intent(in)::template
      type(ram_free_map),intent(in)::map
      real(dp),intent(in)::data(:,:)
      type(sem_fit_result),intent(out)::result
      character(len=*),intent(in),optional::estimator
      real(dp),allocatable::mean(:),cov(:,:)
      call sample_mean_cov(data,mean,cov)
      if(present(estimator)) then
         call fit_ram_cov(template,map,cov,mean,size(data,1),result,estimator)
      else
         call fit_ram_cov(template,map,cov,mean,size(data,1),result)
      end if
   end subroutine fit_ram_data

   subroutine fit_ram_fiml(template,map,data,result)
      use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
      type(ram_model),intent(in)::template
      type(ram_free_map),intent(in)::map
      real(dp),intent(in)::data(:,:)
      type(sem_fit_result),intent(out)::result
      type(ram_model)::work
      real(dp),allocatable::x(:),hess(:,:),hi(:,:)
      real(dp)::fval
      integer::info,status,i,k,p,nobs
      logical::hasmiss
      work=template
      x=ram_get_free(template,map)
      k=size(x)
      p=size(data,2)
      nobs=size(data,1)
      hasmiss=any(ieee_is_nan(data))
      call bfgs_minimize(nll,x,fval,result%converged,result%iterations,maxiter=1200,tol=1e-7_dp)
      work=template
      call ram_set_free(work,map,x)
      call ram_sigma(work,result%sigma,info)
      call ram_mu(work,result%mu,info)
      result%par=x
      result%loglik=-fval
      result%objective=2*fval/real(nobs,dp)
      result%df=real(p*(p+1)/2 + merge(p,0,allocated(template%m))-k,dp)
      result%aic=2*fval+2*real(k,dp)
      result%bic=2*fval+log(real(nobs,dp))*real(k,dp)
      call hessian(nll,x,hess,status=status)
      allocate(result%vcov(k,k),result%se(k))
      result%vcov=0
      result%se=huge(1.0_dp)
      if(status==nd_success) then
      call inverse_general(hess,hi,info)
      if(info==0) then
      result%vcov=hi
      do i=1,k
         if(hi(i,i)>=0) result%se(i)=sqrt(hi(i,i))
         end do
         end if
         end if
      result%status=info
   contains
      function nll(z) result(v)
         real(dp),intent(in)::z(:)
         real(dp)::v
         real(dp),allocatable::sg(:,:),mm(:)
         integer::istat
         work=template
         call ram_set_free(work,map,z)
         call ram_sigma(work,sg,istat)
         if(istat/=0) then
         v=huge(1.0_dp)/100
         return
         end if
         call ram_mu(work,mm,istat)
         if(hasmiss) then
         v=-mvn_loglik_missing(data,mm,sg,istat)
         else
         v=-mvn_loglik_complete(data,mm,sg,istat)
         end if
         if(istat/=0) v=huge(1.0_dp)/100
      end function nll
   end subroutine fit_ram_fiml

   subroutine standardized_ram(model,a_std,s_std,info)
      type(ram_model),intent(in)::model
      real(dp),allocatable,intent(out)::a_std(:,:),s_std(:,:)
      integer,intent(out)::info
      real(dp),allocatable::allcov(:,:),ia(:,:),inv(:,:),sd(:)
      integer::n,i,j
      n=size(model%a,1)
      allocate(ia(n,n))
      ia=-model%a
      do i=1,n
      ia(i,i)=ia(i,i)+1
      end do
      call inverse_general(ia,inv,info)
      allocate(a_std(n,n),s_std(n,n),sd(n))
      if(info/=0) return
      allcov=matmul(inv,matmul(model%s,transpose(inv)))
      do i=1,n
      sd(i)=sqrt(max(allcov(i,i),tiny(1.0_dp)))
      end do
      do j=1,n
      do i=1,n
      a_std(i,j)=model%a(i,j)*sd(j)/sd(i)
      s_std(i,j)=model%s(i,j)/(sd(i)*sd(j))
      end do
      end do
   end subroutine standardized_ram

   function residual_covariance(observed,implied,standardized) result(resid)
      real(dp),intent(in)::observed(:,:),implied(:,:)
      logical,intent(in),optional::standardized
      real(dp),allocatable::resid(:,:)
      logical::std
      real(dp),allocatable::sd(:)
      integer::i,j,p
      std=.false.
      if(present(standardized)) std=standardized
      resid=observed-implied
      if(std) then
      p=size(observed,1)
      allocate(sd(p))
      do i=1,p
      sd(i)=sqrt(max(observed(i,i),tiny(1.0_dp)))
      end do
         do j=1,p
         do i=1,p
         resid(i,j)=resid(i,j)/(sd(i)*sd(j))
         end do
         end do
      end if
   end function residual_covariance

   function spd_inverse(a,info) result(inv)
      use lavaan_linalg, only : inverse_spd
      real(dp),intent(in)::a(:,:)
      integer,intent(out)::info
      real(dp),allocatable::inv(:,:)
      call inverse_spd(a,inv,info)
   end function spd_inverse

   function compute_srmr(obs,imp) result(v)
      real(dp),intent(in)::obs(:,:),imp(:,:)
      real(dp)::v,s
      integer::p,i,j,n
      real(dp),allocatable::sd(:)
      p=size(obs,1)
      allocate(sd(p))
      do i=1,p
      sd(i)=sqrt(max(obs(i,i),tiny(1.0_dp)))
      end do
      s=0
      n=0
      do j=1,p
      do i=j,p
      s=s+((obs(i,j)-imp(i,j))/(sd(i)*sd(j)))**2
      n=n+1
      end do
      end do
      v=sqrt(s/real(n,dp))
   end function compute_srmr
end module lavaan_fit
