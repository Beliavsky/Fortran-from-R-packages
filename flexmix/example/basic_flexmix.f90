program basic_flexmix
   use flexmix, only : dp, flexmix_result, flexmix_control, flexmix_gaussian
   implicit none
   integer, parameter :: n = 20
   real(dp) :: x(n,2), y(n)
   type(flexmix_result) :: fit
   type(flexmix_control) :: control
   integer :: i

   control%minprior = 0.01_dp
   do i = 1, n
      x(i,1) = 1.0_dp
      x(i,2) = real(mod(i-1,10),dp) / 9.0_dp
      if (i <= 10) then
         y(i) = 1.0_dp + 2.0_dp * x(i,2) + 0.02_dp * real((-1)**i,dp)
      else
         y(i) = 8.0_dp - 1.5_dp * x(i,2) + 0.02_dp * real((-1)**i,dp)
      end if
   end do

   call flexmix_gaussian(x, y, 2, fit, control)
   write(*,'(a,l1)') 'converged: ', fit%converged
   write(*,'(a,i0)') 'components: ', fit%k
   write(*,'(a,f12.6)') 'log likelihood: ', fit%loglik
   write(*,'(a,*(1x,f10.6))') 'priors:', fit%prior
   write(*,'(a)') 'component coefficients (columns):'
   do i = 1, size(fit%beta,1)
      write(*,'(*(1x,f12.6))') fit%beta(i,:)
   end do
end program basic_flexmix
