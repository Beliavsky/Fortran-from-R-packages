program parity_driver
   use grpreg
   implicit none
   integer, parameter :: n = 40, p = 4
   real(dp) :: x(n,p), yg(n), yb(n), yp(n), time(n), event(n)
   integer :: group(p), i, u
   type(grpreg_fit_type) :: fit
   type(spline_expansion_type) :: spl
   do i = 1, n
      x(i,1) = (real(i,dp)-20.5_dp)/12.0_dp
      x(i,2) = sin(0.37_dp*real(i,dp))
      x(i,3) = cos(0.23_dp*real(i,dp))
      x(i,4) = sin(0.11_dp*real(i*i,dp))
      yb(i) = merge(1.0_dp,0.0_dp,mod(i,5) == 0 .or. mod(i,7) == 0 .or. x(i,2) > 0.65_dp)
      yp(i) = real(mod(i,4),dp)
      time(i) = real(mod(17*i,41)+1,dp)+0.01_dp*real(i,dp)
      event(i) = merge(1.0_dp,0.0_dp,mod(5*i+2,7) < 4)
   end do
   group = [1,1,2,2]
   yg = 1.25_dp+1.7_dp*x(:,1)-0.9_dp*x(:,2)+0.6_dp*x(:,3)-0.3_dp*x(:,4)
   open(newunit=u,file='fortran_parity.csv',status='replace',action='write')
   write(u,'(a)') 'case,lambda,intercept,beta1,beta2,beta3,beta4,deviance'
   call grpreg_fit(x,yg,group,fit,family='gaussian',penalty='grLasso',lambda=[0.15_dp],eps=1.0e-11_dp)
   call write_fit(u,'gaussian_group_lasso',fit,1)
   call grpreg_fit(x,yb,group,fit,family='binomial',penalty='grLasso',lambda=[1.0e-10_dp],eps=1.0e-10_dp)
   call write_fit(u,'binomial_near_mle',fit,1)
   call grpreg_fit(x,yp,group,fit,family='poisson',penalty='grLasso',lambda=[1.0e-10_dp],eps=1.0e-10_dp)
   call write_fit(u,'poisson_near_mle',fit,1)
   call grpsurv_fit(x,time,event,group,fit,penalty='grLasso',lambda=[1.0e-10_dp],eps=1.0e-10_dp)
   call write_fit(u,'cox_near_mle',fit,1)
   call expand_spline(x(:,1:1),spl,df=4,degree=3,spline_type='bs')
   write(u,'(a)') 'spline_row,0,0,'//trim(real_text(spl%x(13,1)))//','//trim(real_text(spl%x(13,2)))//','// &
      trim(real_text(spl%x(13,3)))//','//trim(real_text(spl%x(13,4)))//',0'
   close(u)
contains
   subroutine write_fit(unit, name, object, index)
      integer, intent(in) :: unit !! Open output unit receiving one parity row.
      character(len=*), intent(in) :: name !! Case identifier written as the first CSV field.
      type(grpreg_fit_type), intent(in) :: object !! Fitted path supplying coefficients and deviance.
      integer, intent(in) :: index !! One-based fitted-lambda index to serialize.
      write(unit,'(a,",",es24.16,",",es24.16,5(",",es24.16))') trim(name),object%lambda(index), &
         object%intercept(index),object%beta(1,index),object%beta(2,index),object%beta(3,index), &
         object%beta(4,index),object%deviance(index)
   end subroutine write_fit
   function real_text(x) result(text)
      real(dp), intent(in) :: x !! Scalar real formatted for the parity CSV helper row.
      character(len=32) :: text
      write(text,'(es24.16)') x
   end function real_text
end program parity_driver
