! SPDX-License-Identifier: GPL-3.0-only
module ufrisk_loggarch
   use kind_mod, only : dp
   use distribution_mod, only : standardized_log_density
   implicit none
   private
   public :: arfilt_coefficients, log_variance_forecasts, estimate_student_df
contains
   pure function arfilt_coefficients(ar,ma,d,k) result(coefficients)
      real(dp), intent(in) :: ar(:),ma(:),d
      integer, intent(in) :: k
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: ar_polynomial(:), arma_polynomial(:), fractional(:), combined(:)
      integer :: p,q,i,j,upper
      p = count(abs(ar)>0.0_dp)
      q = count(abs(ma)>0.0_dp)
      allocate(coefficients(max(0,k)))
      if (k <= 0) return
      allocate(ar_polynomial(0:k),arma_polynomial(0:k),fractional(0:k),combined(0:k))
      ar_polynomial = 0.0_dp; arma_polynomial = 0.0_dp
      fractional = 0.0_dp; combined = 0.0_dp
      ar_polynomial(0) = 1.0_dp
      if (p > 0) ar_polynomial(1:min(p,k)) = ar(1:min(p,k))
      arma_polynomial(0) = 1.0_dp
      if (p > 0 .or. q > 0) then
         do i = 1,k
            upper = min(q,i)
            if (upper > 0) then
               do j = 1,upper
                  arma_polynomial(i) = arma_polynomial(i)+ma(j)*arma_polynomial(i-j)
               end do
            end if
            arma_polynomial(i) = arma_polynomial(i)-ar_polynomial(i)
         end do
      end if
      fractional(0) = 1.0_dp
      do i = 1,k
         fractional(i) = -fractional(i-1)*(d-real(i-1,dp))/real(i,dp)
      end do
      do i = 0,k
         do j = 0,i
            combined(i) = combined(i)+fractional(j)*arma_polynomial(i-j)
         end do
      end do
      coefficients = -combined(1:k)
   end function arfilt_coefficients

   pure subroutine log_variance_forecasts(observed_in,observed_out,coefficients,center,forecast)
      real(dp), intent(in) :: observed_in(:),observed_out(:),coefficients(:),center
      real(dp), allocatable, intent(out) :: forecast(:)
      real(dp), allocatable :: history(:)
      integer :: n_in,n_out,k,t,j,index
      n_in = size(observed_in); n_out = size(observed_out); k = size(coefficients)
      allocate(forecast(n_out),history(k+n_in+n_out)); history = 0.0_dp
      history(k+1:k+n_in) = observed_in-center
      history(k+n_in+1:k+n_in+n_out) = observed_out-center
      do t = 1,n_out
         index = k+n_in+t
         forecast(t) = 0.0_dp
         do j = 1,k
            forecast(t) = forecast(t)+coefficients(j)*history(index-j)
         end do
      end do
   end subroutine log_variance_forecasts

   pure real(dp) function estimate_student_df(x,lower,upper) result(degrees_freedom)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in), optional :: lower,upper
      real(dp) :: a,b,c,d,fc,fd,lo,hi,ratio
      integer :: iteration
      lo = 2.05_dp
      hi = 100.0_dp
      if (present(lower)) lo = max(2.0001_dp,lower)
      if (present(upper)) hi = max(lo+0.1_dp,upper)
      ratio = 0.5_dp*(sqrt(5.0_dp)-1.0_dp)
      a = lo; b = hi
      c = b-ratio*(b-a); d = a+ratio*(b-a)
      fc = negative_log_likelihood(c)
      fd = negative_log_likelihood(d)
      do iteration = 1,100
         if (fc < fd) then
            b = d; d = c; fd = fc; c = b-ratio*(b-a); fc = negative_log_likelihood(c)
         else
            a = c; c = d; fc = fd; d = a+ratio*(b-a); fd = negative_log_likelihood(d)
         end if
         if (b-a < 1.0e-6_dp*max(1.0_dp,0.5_dp*(a+b))) exit
      end do
      degrees_freedom = 0.5_dp*(a+b)
   contains
      pure real(dp) function negative_log_likelihood(nu) result(value)
         real(dp), intent(in) :: nu
         integer :: i
         value = 0.0_dp
         do i = 1,size(x)
            value = value-standardized_log_density(x(i),2,nu,1.0_dp,-0.5_dp)
         end do
      end function negative_log_likelihood
   end function estimate_student_df
end module ufrisk_loggarch
