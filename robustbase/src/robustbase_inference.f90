! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of robustbase.
! It may be redistributed and/or modified under GPL version 2 or later.
module robustbase_inference
   use robustbase_kinds, only: dp
   use robustbase_linalg, only: invert_symmetric, symmetric_eigen
   use robustbase_probability, only: normal_quantile, chi_square_quantile, chi_square_cdf
   use robustbase_lmrob, only: robust_m_scale
   use robustbase_sort, only: median
   implicit none
   private
   public :: robust_prediction_result, robust_test_result, robust_outlier_stats_result, &
             robust_linear_predict, robust_wald_test, robust_deviance_test, robust_r_squared, &
             robust_outlier_stats, tolerance_ellipse_points

   type :: robust_prediction_result
      real(dp),allocatable :: fit(:)
      real(dp),allocatable :: standard_error(:)
      real(dp),allocatable :: lower(:)
      real(dp),allocatable :: upper(:)
      character(len=12) :: interval='confidence'
      real(dp) :: level=0.95_dp
   end type robust_prediction_result

   type :: robust_test_result
      real(dp) :: statistic=0.0_dp
      real(dp) :: p_value=1.0_dp
      integer :: degrees_freedom=0
      character(len=12) :: test='Wald'
   end type robust_test_result

   type :: robust_outlier_stats_result
      real(dp),allocatable :: leverage(:)
      real(dp),allocatable :: standardized_residuals(:)
      logical,allocatable :: high_leverage(:)
      logical,allocatable :: outlier(:)
      integer :: n_high_leverage=0
      integer :: n_outliers=0
      real(dp) :: leverage_cutoff=0.0_dp
      real(dp) :: residual_cutoff=2.5_dp
   end type robust_outlier_stats_result
contains
   subroutine robust_linear_predict(x_new,coefficients,covariance,scale,result,level,interval,df,prediction_variance)
      real(dp),intent(in)::x_new(:,:),coefficients(:),covariance(:,:),scale
      type(robust_prediction_result),intent(out)::result
      real(dp),intent(in),optional::level,prediction_variance
      character(len=*),intent(in),optional::interval
      integer,intent(in),optional::df
      integer::n,p,i,nu
      real(dp)::lev,z,predvar
      n=size(x_new,1);p=size(x_new,2)
      if(size(coefficients)/=p .or. any(shape(covariance)/=[p,p]))error stop 'robust_linear_predict: size mismatch'
      result%level=0.95_dp;if(present(level))result%level=level
      result%interval='confidence';if(present(interval))result%interval=adjustl(interval)
      if(result%level<=0.0_dp .or. result%level>=1.0_dp)error stop 'robust_linear_predict: level must be in (0,1)'
      nu=1000000;if(present(df))nu=max(1,df)
      z=student_t_quantile_approx(0.5_dp+0.5_dp*result%level,real(nu,dp))
      predvar=scale*scale;if(present(prediction_variance))predvar=max(0.0_dp,prediction_variance)
      allocate(result%fit(n),result%standard_error(n),result%lower(n),result%upper(n))
      result%fit=matmul(x_new,coefficients)
      do i=1,n
         lev=max(0.0_dp,dot_product(x_new(i,:),matmul(covariance,x_new(i,:))))
         if(trim(result%interval)=='prediction')lev=lev+predvar
         result%standard_error(i)=sqrt(lev)
      end do
      result%lower=result%fit-z*result%standard_error;result%upper=result%fit+z*result%standard_error
   end subroutine robust_linear_predict

   subroutine robust_wald_test(coefficients,covariance,indices,result,null_values)
      real(dp),intent(in)::coefficients(:),covariance(:,:)
      integer,intent(in)::indices(:)
      type(robust_test_result),intent(out)::result
      real(dp),intent(in),optional::null_values(:)
      real(dp),allocatable::diff(:),subcov(:,:),inv(:,:)
      integer::k,i,j,info
      k=size(indices)
      if(k<1)error stop 'robust_wald_test: empty index set'
      allocate(diff(k),subcov(k,k),inv(k,k))
      do i=1,k
         if(indices(i)<1 .or. indices(i)>size(coefficients))error stop 'robust_wald_test: invalid index'
         diff(i)=coefficients(indices(i))
      end do
      if(present(null_values))then
         if(size(null_values)/=k)error stop 'robust_wald_test: null size mismatch'
         diff=diff-null_values
      end if
      do i=1,k;do j=1,k;subcov(i,j)=covariance(indices(i),indices(j));end do;end do
      call invert_symmetric(subcov,inv,info,ridge=1.0e-12_dp)
      if(info/=0)error stop 'robust_wald_test: covariance inversion failed'
      result%statistic=max(0.0_dp,dot_product(diff,matmul(inv,diff)))
      result%degrees_freedom=k;result%p_value=1.0_dp-chi_square_cdf(result%statistic,real(k,dp));result%test='Wald'
   end subroutine robust_wald_test

   subroutine robust_deviance_test(full_residuals,reduced_residuals,scale,result,tuning,degrees_freedom)
      real(dp),intent(in)::full_residuals(:),reduced_residuals(:),scale
      type(robust_test_result),intent(out)::result
      real(dp),intent(in),optional::tuning
      integer,intent(in),optional::degrees_freedom
      real(dp)::c,tau_star,mean_deriv,mean_psi2,full_rho,reduced_rho
      integer::df
      if(size(full_residuals)/=size(reduced_residuals))error stop 'robust_deviance_test: size mismatch'
      c=4.685061_dp;if(present(tuning))c=tuning
      mean_deriv=sum(tukey_derivative(full_residuals/max(scale,1.0e-14_dp),c))/real(size(full_residuals),dp)
      mean_psi2=sum(tukey_psi_local(full_residuals/max(scale,1.0e-14_dp),c)**2)/real(size(full_residuals),dp)
      tau_star=mean_deriv/max(mean_psi2,1.0e-14_dp)
      full_rho=sum(tukey_rho_local(full_residuals/max(scale,1.0e-14_dp),c))
      reduced_rho=sum(tukey_rho_local(reduced_residuals/max(scale,1.0e-14_dp),c))
      result%statistic=max(0.0_dp,2.0_dp*tau_star*(reduced_rho-full_rho))
      df=1;if(present(degrees_freedom))df=max(1,degrees_freedom);result%degrees_freedom=df;result%p_value=1.0_dp-chi_square_cdf(result%statistic,real(df,dp));result%test='Deviance'
   end subroutine robust_deviance_test

   function robust_r_squared(y,residuals,tuning) result(r2)
      real(dp),intent(in)::y(:),residuals(:)
      real(dp),intent(in),optional::tuning
      real(dp)::r2,c,sres,stot,med
      if(size(y)/=size(residuals))error stop 'robust_r_squared: size mismatch'
      c=1.54764_dp;if(present(tuning))c=tuning
      med=median(y);sres=robust_m_scale(residuals,tuning=c,b=0.5_dp);stot=robust_m_scale(y-med,tuning=c,b=0.5_dp)
      if(stot<=1.0e-14_dp)then;r2=0.0_dp;else;r2=1.0_dp-(sres/stot)**2;end if
   end function robust_r_squared

   subroutine robust_outlier_stats(x,residuals,scale,result,leverage_multiplier,residual_cutoff)
      real(dp),intent(in)::x(:,:),residuals(:),scale
      type(robust_outlier_stats_result),intent(out)::result
      real(dp),intent(in),optional::leverage_multiplier,residual_cutoff
      real(dp),allocatable::xtx(:,:),inv(:,:)
      real(dp)::lm
      integer::n,p,i,info
      n=size(x,1);p=size(x,2)
      if(size(residuals)/=n)error stop 'robust_outlier_stats: size mismatch'
      lm=2.0_dp;if(present(leverage_multiplier))lm=leverage_multiplier
      result%residual_cutoff=2.5_dp;if(present(residual_cutoff))result%residual_cutoff=residual_cutoff
      result%leverage_cutoff=lm*real(p,dp)/real(n,dp)
      allocate(result%leverage(n),result%standardized_residuals(n),result%high_leverage(n),result%outlier(n),xtx(p,p),inv(p,p))
      xtx=matmul(transpose(x),x);call invert_symmetric(xtx,inv,info,ridge=1.0e-12_dp)
      if(info/=0)error stop 'robust_outlier_stats: inverse failed'
      do i=1,n
         result%leverage(i)=dot_product(x(i,:),matmul(inv,x(i,:)))
      end do
      result%standardized_residuals=residuals/(max(scale,1.0e-14_dp)*sqrt(max(1.0_dp-result%leverage,1.0e-8_dp)))
      result%high_leverage=result%leverage>result%leverage_cutoff
      result%outlier=abs(result%standardized_residuals)>result%residual_cutoff
      result%n_high_leverage=count(result%high_leverage);result%n_outliers=count(result%outlier)
   end subroutine robust_outlier_stats

   subroutine tolerance_ellipse_points(center,covariance,points,probability,n_points)
      real(dp),intent(in)::center(:),covariance(:,:)
      real(dp),allocatable,intent(out)::points(:,:)
      real(dp),intent(in),optional::probability
      integer,intent(in),optional::n_points
      real(dp)::prob,q,angle
      real(dp)::values(2),vectors(2,2),radius(2),unit(2)
      integer::np,i,info
      if(size(center)/=2 .or. any(shape(covariance)/=[2,2]))error stop 'tolerance_ellipse_points: requires two dimensions'
      prob=0.975_dp;if(present(probability))prob=probability
      np=629;if(present(n_points))np=max(4,n_points)
      call symmetric_eigen(covariance,values,vectors,info)
      if(info/=0 .or. any(values<0.0_dp))error stop 'tolerance_ellipse_points: covariance is not PSD'
      q=chi_square_quantile(prob,2.0_dp);radius=sqrt(q*max(values,0.0_dp))
      allocate(points(np,2))
      do i=1,np
         angle=2.0_dp*acos(-1.0_dp)*real(i-1,dp)/real(np-1,dp)
         unit=[radius(1)*cos(angle),radius(2)*sin(angle)]
         points(i,:)=center+matmul(vectors,unit)
      end do
   end subroutine tolerance_ellipse_points

   elemental function student_t_quantile_approx(p,df) result(q)
      real(dp),intent(in)::p,df
      real(dp)::q,z,z2
      z=normal_quantile(p);z2=z*z
      q=z+(z*z2+z)/(4.0_dp*df)+(5.0_dp*z**5+16.0_dp*z2*z+3.0_dp*z)/(96.0_dp*df*df) &
         +(3.0_dp*z**7+19.0_dp*z**5+17.0_dp*z**3-15.0_dp*z)/(384.0_dp*df**3)
   end function student_t_quantile_approx

   elemental function tukey_psi_local(x,c) result(value)
      real(dp),intent(in)::x,c
      real(dp)::value,u
      u=x/c;if(abs(u)<1.0_dp)then;value=x*(1.0_dp-u*u)**2;else;value=0.0_dp;end if
   end function tukey_psi_local
   elemental function tukey_derivative(x,c) result(value)
      real(dp),intent(in)::x,c
      real(dp)::value,u
      u=x/c;if(abs(u)<1.0_dp)then;value=(1.0_dp-u*u)*(1.0_dp-5.0_dp*u*u);else;value=0.0_dp;end if
   end function tukey_derivative
   elemental function tukey_rho_local(x,c) result(value)
      real(dp),intent(in)::x,c
      real(dp)::value,u
      u=x/c;if(abs(u)<1.0_dp)then;value=c*c/6.0_dp*(1.0_dp-(1.0_dp-u*u)**3);else;value=c*c/6.0_dp;end if
   end function tukey_rho_local
end module robustbase_inference
