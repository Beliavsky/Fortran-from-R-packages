! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of robustbase.
! It may be redistributed and/or modified under GPL version 2 or later.
module nonlinear_example_model
   use robustbase_kinds, only: dp
   implicit none
contains
   subroutine exponential_model(theta,x,yhat)
      real(dp),intent(in)::theta(:),x(:,:)
      real(dp),intent(out)::yhat(:)
      yhat=theta(1)*exp(theta(2)*x(:,1))
   end subroutine exponential_model
end module nonlinear_example_model
program nonlinear_example
   use robustbase
   use nonlinear_example_model, only: exponential_model
   implicit none
   real(dp)::x(50,1),y(50),start(2)
   type(robust_nls_result)::fit
   integer::i
   do i=1,50
      x(i,1)=2.0_dp*real(i-1,dp)/49.0_dp
      y(i)=1.5_dp*exp(0.7_dp*x(i,1))+0.02_dp*sin(real(i,dp))
   end do
   y(47:50)=y(47:50)+8.0_dp
   start=[1.0_dp,0.3_dp]
   call robust_nls_fit(exponential_model,x,y,start,fit)
   write(*,'(a,*(1x,f12.6))')'parameters',fit%parameters
   write(*,'(a,1x,l1)')'converged',fit%converged
end program nonlinear_example
