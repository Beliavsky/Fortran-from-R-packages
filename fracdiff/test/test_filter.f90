! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran port of fracdiff; see NOTICE.md for attribution.

program test_filter
   use fracdiff, only : dp, fd_ok, haslett_raftery_filter, arma_residual_jacobian
   use test_support
   implicit none

   real(dp) :: x(10), y(10), expected_y(10), sum_log_v, estimated_mean
   real(dp) :: ar(1), ma(1), residuals(9), jacobian(9,2)
   real(dp) :: residual_plus(9), residual_minus(9), numerical(9)
   real(dp) :: ar_work(1), ma_work(1), step
   integer :: status, j

   x = [0.3_dp,-0.2_dp,1.1_dp,0.7_dp,-0.4_dp,0.9_dp,1.5_dp,-0.6_dp,0.2_dp,0.8_dp]
   expected_y = [-0.2332098827149459_dp,-0.6883135612841985_dp, &
      0.7329786487868368_dp,0.0733218427924248_dp,-1.0148027646578506_dp, &
      0.5208192372608322_dp,0.9129588070380624_dp,-1.3861457409121550_dp, &
      -0.1941674518906057_dp,0.3533699606521097_dp]

   call haslett_raftery_filter(x,0.23_dp,5,.true.,y,sum_log_v,estimated_mean,status)
   call assert_true(status==fd_ok,"Haslett-Raftery filter status")
   call assert_close(estimated_mean,0.44107591886122954_dp,2.0e-14_dp,"filtered mean")
   call assert_close(sum_log_v,0.22488117121652798_dp,2.0e-14_dp,"sum log innovation variance")
   call assert_vector_close(y,expected_y,3.0e-14_dp,"Haslett-Raftery filtered series")

   ar=0.2_dp
   ma=-0.35_dp
   call arma_residual_jacobian(y,ar,ma,residuals,jacobian,status)
   call assert_true(status==fd_ok,"ARMA Jacobian status")
   step=1.0e-6_dp
   do j=1,2
      ar_work=ar
      ma_work=ma
      if(j==1) then
         ma_work(1)=ma_work(1)+step
      else
         ar_work(1)=ar_work(1)+step
      end if
      call residual_only(y,ar_work,ma_work,residual_plus,status)
      ar_work=ar
      ma_work=ma
      if(j==1) then
         ma_work(1)=ma_work(1)-step
      else
         ar_work(1)=ar_work(1)-step
      end if
      call residual_only(y,ar_work,ma_work,residual_minus,status)
      numerical=(residual_plus-residual_minus)/(2.0_dp*step)
      call assert_vector_close(jacobian(:,j),numerical,2.0e-8_dp,"analytical ARMA Jacobian")
   end do

   write(*,'(a)') "test_filter: PASS"

contains

   subroutine residual_only(values,ar_values,ma_values,res,status_out)
      use fracdiff, only : arma_residuals
      real(dp),intent(in)::values(:),ar_values(:),ma_values(:)
      real(dp),intent(out)::res(:)
      integer,intent(out)::status_out
      call arma_residuals(values,ar_values,ma_values,res,status_out)
   end subroutine residual_only

end program test_filter
