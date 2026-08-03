! SPDX-License-Identifier: GPL-3.0-or-later
program test_cla_queries
   use kind_mod, only: dp
   use cla, only: cla_result_t, cla_path_query_t, critical_line, find_sigma, find_mu
   implicit none
   real(dp) :: mu(3), covar(3,3), lower(3), upper(3)
   type(cla_result_t) :: result
   type(cla_path_query_t) :: sig_result, mu_result

   mu = [0.0408_dp,0.102_dp,-0.023_dp]
   covar = reshape([0.00648_dp,0.00792_dp,0.00473_dp, &
                    0.00792_dp,0.0334_dp,0.0121_dp, &
                    0.00473_dp,0.0121_dp,0.0793_dp],[3,3])
   lower = 0.0_dp
   upper = 1.0_dp
   result = critical_line(mu,covar,lower,upper)

   sig_result = find_sigma([0.04_dp],result,covar)
   call check(sig_result%info == 0,'find_sigma status')
   call close(sig_result%value(1),0.080300141561153_dp,2.0e-10_dp,'interpolated sigma')
   call close(maxval(abs(sig_result%weights(:,1)- &
      [0.987460815047022_dp,0.0_dp,0.012539184952978_dp])),0.0_dp,2.0e-10_dp, &
      'interpolated weights')

   mu_result = find_mu(sig_result%value,result,covar,tolerance=1.0e-11_dp)
   call check(mu_result%info == 0,'find_mu status')
   call close(mu_result%value(1),0.04_dp,2.0e-9_dp,'inverted mean')
   call close(sum(mu_result%weights(:,1)),1.0_dp,1.0e-12_dp,'query budget')
   print '(a)', 'test_cla_queries: PASS'

contains
   subroutine close(actual,expected,tolerance,label)
      real(dp),intent(in)::actual,expected,tolerance
      character(len=*),intent(in)::label
      if(abs(actual-expected)>tolerance)then
         write(*,'(a,2es24.14)')trim(label)//' failed: ',actual,expected
         error stop 1
      end if
   end subroutine close
   subroutine check(condition,label)
      logical,intent(in)::condition
      character(len=*),intent(in)::label
      if(.not.condition)then
         write(*,'(a)')trim(label)//' failed'
         error stop 1
      end if
   end subroutine check
end program test_cla_queries
