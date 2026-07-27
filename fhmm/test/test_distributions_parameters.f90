! SPDX-License-Identifier: GPL-3.0-only
program test_distributions_parameters
   use fhmm
   implicit none
   type(hmm_parameters) :: par, par2
   real(dp), allocatable :: u(:)
   real(dp) :: expected

   call assert_close(distribution_pdf(dist_normal,0.0_dp,0.0_dp,1.0_dp,10.0_dp), &
      0.3989422804014327_dp,1.0e-13_dp,'normal density')
   call assert_close(distribution_cdf(dist_normal,0.0_dp,0.0_dp,1.0_dp,10.0_dp),0.5_dp,1.0e-15_dp,'normal cdf')
   expected=gamma(2.5_dp)/(sqrt(4.0_dp*pi)*gamma(2.0_dp))
   call assert_close(distribution_pdf(dist_student_t,0.0_dp,0.0_dp,1.0_dp,4.0_dp),expected,1.0e-13_dp,'t density')
   call assert_close(distribution_pdf(dist_gamma,2.0_dp,2.0_dp,1.0_dp,10.0_dp), &
      exp(3.0_dp*log(2.0_dp)-4.0_dp-log_gamma(4.0_dp)-4.0_dp*log(0.5_dp)), &
      1.0e-13_dp,'gamma density')
   call assert_close(distribution_pdf(dist_poisson,3.0_dp,2.0_dp,1.0_dp,10.0_dp), &
      exp(3.0_dp*log(2.0_dp)-2.0_dp-log_gamma(4.0_dp)),1.0e-13_dp,'poisson density')

   par%distribution=dist_student_t
   allocate(par%gamma(2,2),par%mu(2),par%sigma(2),par%df(2))
   par%gamma=reshape([0.9_dp,0.2_dp,0.1_dp,0.8_dp],[2,2])
   par%mu=[-1.0_dp,1.0_dp];par%sigma=[0.5_dp,1.5_dp];par%df=[5.0_dp,9.0_dp]
   u=pack_hmm_parameters(par)
   par2=unpack_hmm_parameters(u,2,dist_student_t)
   call assert_array_close(par2%gamma,par%gamma,2.0e-14_dp,'transition roundtrip')
   call assert_vec_close(par2%mu,par%mu,2.0e-14_dp,'mu roundtrip')
   call assert_vec_close(par2%sigma,par%sigma,2.0e-14_dp,'sigma roundtrip')
   call assert_vec_close(par2%df,par%df,2.0e-14_dp,'df roundtrip')
   if(.not.validate_hmm_parameters(par2))error stop 1
   print '(a)','test_distributions_parameters: PASS'
contains
   subroutine assert_close(x,y,tol,label)
      real(dp),intent(in)::x,y,tol;character(len=*),intent(in)::label
      if(abs(x-y)>tol)then;write(*,*)trim(label),x,y,abs(x-y);error stop 1;end if
   end subroutine
   subroutine assert_vec_close(x,y,tol,label)
      real(dp),intent(in)::x(:),y(:),tol;character(len=*),intent(in)::label
      if(maxval(abs(x-y))>tol)then;write(*,*)trim(label),maxval(abs(x-y));error stop 1;end if
   end subroutine
   subroutine assert_array_close(x,y,tol,label)
      real(dp),intent(in)::x(:,:),y(:,:),tol;character(len=*),intent(in)::label
      if(maxval(abs(x-y))>tol)then;write(*,*)trim(label),maxval(abs(x-y));error stop 1;end if
   end subroutine
end program test_distributions_parameters
