! SPDX-License-Identifier: GPL-2.0-or-later
program randomness_demo
   use randtests, only: dp, randomness_result, rank_test, runs_test, bartels_rank_test, pruns
   implicit none
   real(dp), parameter :: x(8) = [4.0_dp, 1.0_dp, 8.0_dp, 2.0_dp, 6.0_dp, 3.0_dp, 7.0_dp, 5.0_dp]
   type(randomness_result) :: result
   result = rank_test(x)
   call show('Rank test', result)
   result = runs_test(x, pvalue='exact')
   call show('Exact runs test', result)
   result = bartels_rank_test(x, pvalue='exact')
   call show('Exact Bartels test', result)
   print '(a,f12.8)', 'P(R <= 3 | n1=2,n2=3) = ', pruns(3.0_dp, 2, 3)
contains
   subroutine show(label, result)
      character(len=*), intent(in) :: label !! Human-readable name of the statistical test being displayed.
      type(randomness_result), intent(in) :: result !! Completed result, whose status must be checked before use.
      if (result % status /= 0) then
         print '(a)', trim(result % message)
         error stop 1
      end if
      print '(a,a,f12.8,a,f12.8)', label, ': statistic = ', result % statistic, ', p = ', result % p_value
   end subroutine show
end program randomness_demo
