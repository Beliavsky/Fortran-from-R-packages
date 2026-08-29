module lavaan_sam
   use lavaan_kinds, only : dp
   use lavaan_ram, only : ram_model, ram_free_map, ram_set_free
   use lavaan_fit, only : sem_fit_result, fit_ram_cov, fit_ram_data
   use lavaan_linalg, only : inverse_general
   use lavaan_robust_tests, only : scaled_test_result, scaled_tests_from_ugamma
   implicit none
   private

   type, public :: sam_result
      type(sem_fit_result) :: measurement
      type(sem_fit_result) :: structural
      type(ram_model) :: structural_model
      real(dp), allocatable :: stage1_jacobian(:, :)
      real(dp), allocatable :: structural_vcov_corrected(:, :)
      real(dp), allocatable :: structural_se_corrected(:)
      integer :: status=0
      logical :: uncertainty_propagated=.false.
   end type sam_result

   public :: sam_fit_cov, sam_fit_data, sam_fix_measurement, sam_propagate_uncertainty
   public :: sam_yuan_chan_test

contains

   subroutine sam_fit_cov(measurement_template,measurement_map,structural_template,structural_map, &
                          data_cov,data_mean,nobs,result,estimator)
      type(ram_model),intent(in)::measurement_template,structural_template
      type(ram_free_map),intent(in)::measurement_map,structural_map
      real(dp),intent(in)::data_cov(:,:),data_mean(:)
      integer,intent(in)::nobs
      type(sam_result),intent(out)::result
      character(len=*),intent(in),optional::estimator
      type(ram_model)::fixed
      if(present(estimator)) then
         call fit_ram_cov(measurement_template,measurement_map,data_cov,data_mean,nobs,result%measurement,estimator)
      else
         call fit_ram_cov(measurement_template,measurement_map,data_cov,data_mean,nobs,result%measurement)
      end if
      if(.not.result%measurement%converged) then
      result%status=1
      return
      end if
      call sam_fix_measurement(structural_template,measurement_map,result%measurement%par,fixed)
      result%structural_model=fixed
      if(present(estimator)) then
         call fit_ram_cov(fixed,structural_map,data_cov,data_mean,nobs,result%structural,estimator)
      else
         call fit_ram_cov(fixed,structural_map,data_cov,data_mean,nobs,result%structural)
      end if
      result%status=merge(0,2,result%structural%converged)
      if(result%status==0) call sam_propagate_uncertainty(structural_template,measurement_map,structural_map, &
         data_cov,data_mean,nobs,result,estimator)
   end subroutine sam_fit_cov

   subroutine sam_fit_data(measurement_template,measurement_map,structural_template,structural_map,data,result,estimator)
      use lavaan_linalg, only : sample_mean_cov
      type(ram_model),intent(in)::measurement_template,structural_template
      type(ram_free_map),intent(in)::measurement_map,structural_map
      real(dp),intent(in)::data(:,:)
      type(sam_result),intent(out)::result
      character(len=*),intent(in),optional::estimator
      real(dp),allocatable::mean(:),cov(:,:)
      call sample_mean_cov(data,mean,cov)
      if(present(estimator)) then
         call sam_fit_cov(measurement_template,measurement_map,structural_template,structural_map, &
                          cov,mean,size(data,1),result,estimator)
      else
         call sam_fit_cov(measurement_template,measurement_map,structural_template,structural_map, &
                          cov,mean,size(data,1),result)
      end if
   end subroutine sam_fit_data

   subroutine sam_propagate_uncertainty(structural_template,measurement_map,structural_map, &
                                        data_cov,data_mean,nobs,result,estimator)
      type(ram_model),intent(in)::structural_template
      type(ram_free_map),intent(in)::measurement_map,structural_map
      real(dp),intent(in)::data_cov(:,:),data_mean(:)
      integer,intent(in)::nobs
      type(sam_result),intent(inout)::result
      character(len=*),intent(in),optional::estimator
      type(ram_model)::fixed
      type(sem_fit_result)::plusfit,minusfit
      real(dp),allocatable::pplus(:),pminus(:),jac(:,:),add(:,:)
      real(dp)::h
      integer::j,k1,k2,i
      if(.not.result%measurement%converged .or. .not.result%structural%converged) return
      if(.not.allocated(result%measurement%vcov) .or. .not.allocated(result%structural%vcov)) return
      k1=size(result%measurement%par)
      k2=size(result%structural%par)
      if(k1==0 .or. k2==0) return
      if(any(shape(result%measurement%vcov)/=[k1,k1])) return
      allocate(jac(k2,k1))
      jac=0.0_dp
      do j=1,k1
         h=sqrt(epsilon(1.0_dp))*max(1.0_dp,abs(result%measurement%par(j)))
         h=max(h,1.0e-3_dp)
         pplus=result%measurement%par
         pminus=result%measurement%par
         pplus(j)=pplus(j)+h
         pminus(j)=pminus(j)-h
         call sam_fix_measurement(structural_template,measurement_map,pplus,fixed)
         if(present(estimator)) then
            call fit_ram_cov(fixed,structural_map,data_cov,data_mean,nobs,plusfit,estimator)
         else
            call fit_ram_cov(fixed,structural_map,data_cov,data_mean,nobs,plusfit)
         end if
         call sam_fix_measurement(structural_template,measurement_map,pminus,fixed)
         if(present(estimator)) then
            call fit_ram_cov(fixed,structural_map,data_cov,data_mean,nobs,minusfit,estimator)
         else
            call fit_ram_cov(fixed,structural_map,data_cov,data_mean,nobs,minusfit)
         end if
         if(plusfit%converged .and. minusfit%converged) then
            jac(:,j)=(plusfit%par-minusfit%par)/(2.0_dp*h)
         else if(plusfit%converged) then
            jac(:,j)=(plusfit%par-result%structural%par)/h
         else if(minusfit%converged) then
            jac(:,j)=(result%structural%par-minusfit%par)/h
         end if
      end do
      add=matmul(jac,matmul(result%measurement%vcov,transpose(jac)))
      result%stage1_jacobian=jac
      result%structural_vcov_corrected=result%structural%vcov+add
      allocate(result%structural_se_corrected(k2))
      do i=1,k2
         result%structural_se_corrected(i)=sqrt(max(0.0_dp,result%structural_vcov_corrected(i,i)))
      end do
      result%uncertainty_propagated=.true.
   end subroutine sam_propagate_uncertainty

   subroutine sam_yuan_chan_test(chisq,df,delta_theta,delta_gamma,p_jac,gamma,w,result,status)
      real(dp),intent(in)::chisq,df,delta_theta(:,:),delta_gamma(:,:),p_jac(:,:),gamma(:,:),w(:,:)
      type(scaled_test_result),intent(out)::result
      integer,intent(out)::status
      real(dp),allocatable::dmat(:,:),gtilde(:,:),e(:,:),ei(:,:),u(:,:),ug(:,:)
      integer::nstat,ntheta,ngamma,info,i
      nstat=size(delta_theta,1)
      ntheta=size(delta_theta,2)
      ngamma=size(delta_gamma,2)
      if(df<=0.0_dp .or. size(delta_gamma,1)/=nstat .or. size(p_jac,1)/=ngamma .or. &
         size(p_jac,2)/=nstat .or. any(shape(gamma)/=[nstat,nstat]) .or. &
         any(shape(w)/=[nstat,nstat]) .or. ntheta<1) then
         status=-1
         result%status=-1
         return
      end if
      allocate(dmat(nstat,nstat))
      dmat=-matmul(delta_gamma,p_jac)
      do i=1,nstat
      dmat(i,i)=dmat(i,i)+1.0_dp
      end do
      gtilde=matmul(dmat,matmul(gamma,transpose(dmat)))
      e=matmul(transpose(delta_theta),matmul(w,delta_theta))
      call inverse_general(e,ei,info)
      if(info/=0) then
      status=info
      result%status=info
      return
      end if
      u=w-matmul(w,matmul(delta_theta,matmul(ei,matmul(transpose(delta_theta),w))))
      ug=matmul(u,gtilde)
      call scaled_tests_from_ugamma(chisq,df,ug,result)
      if(result%status==0) then
         result%yb_scaling=result%sb_scaling
         result%chisq_yb=result%chisq_sb
      end if
      status=result%status
   end subroutine sam_yuan_chan_test

   subroutine sam_fix_measurement(structural_template,measurement_map,measurement_par,fixed)
      type(ram_model),intent(in)::structural_template
      type(ram_free_map),intent(in)::measurement_map
      real(dp),intent(in)::measurement_par(:)
      type(ram_model),intent(out)::fixed
      fixed=structural_template
      call ram_set_free(fixed,measurement_map,measurement_par)
   end subroutine sam_fix_measurement
end module lavaan_sam
