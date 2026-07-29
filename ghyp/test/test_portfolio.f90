! SPDX-License-Identifier: GPL-2.0-or-later
program test_portfolio
   use ghyp
   implicit none
   type(ghyp_model_type) :: mv
   type(portfolio_result) :: fit, target
   real(dp) :: sc(2,2)
   sc=reshape([1.0_dp,0.3_dp,0.3_dp,0.8_dp],[2,2])
   mv=ghyp_mv(0.7_dp,1.4_dp,2.3_dp,[0.1_dp,-0.2_dp],sc,[0.2_dp,-0.1_dp])
   fit=portfolio_optimize(mv,'sd','minimum.risk')
   call assert_true(fit%ok,'minimum-risk portfolio')
   call assert_close(fit%weights(1),0.41297879_dp,2.0e-7_dp,'minimum-risk weight')
   call assert_close(sum(fit%weights),1.0_dp,2.0e-13_dp,'weight sum')
   call assert_close(fit%standard_deviation,0.8845168136201007_dp,3.0e-8_dp,'portfolio sd')
   target=portfolio_optimize(mv,'sd','target.return',target_return=0.0_dp)
   call assert_true(target%ok,'target-return portfolio')
   call assert_close(target%expected_return,0.0_dp,3.0e-10_dp,'target return')
   print '(a)', 'test_portfolio: PASS'
contains
   subroutine assert_close(actual,expected,tol,label)
      real(dp),intent(in)::actual,expected,tol
      character(len=*),intent(in)::label
      if(abs(actual-expected)>tol*(1.0_dp+abs(expected)))then
         write(*,'(a,3es24.16)')trim(label)//' mismatch: ',actual,expected,abs(actual-expected)
         error stop 1
      end if
   end subroutine assert_close
   subroutine assert_true(condition,label)
      logical,intent(in)::condition
      character(len=*),intent(in)::label
      if(.not.condition)then;write(*,'(a)')trim(label)//' failed';error stop 1;end if
   end subroutine assert_true
end program test_portfolio
