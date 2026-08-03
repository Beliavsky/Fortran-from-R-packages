program journal_and_pl
   use pmwr, only : dp, journal_type, pl_summary_result, pl_path_result, &
                    make_journal, positions_at, pl_summary, pl_path
   implicit none
   type(journal_type) :: journal
   type(pl_summary_result) :: summary
   type(pl_path_result) :: path
   real(dp), allocatable :: position(:,:)

   call make_journal(timestamp=[1.0_dp, 2.0_dp, 3.0_dp], &
                     amount=[10.0_dp, -4.0_dp, -6.0_dp], &
                     price=[100.0_dp, 110.0_dp, 120.0_dp], &
                     instrument=[1, 1, 1], journal=journal)
   call positions_at(journal, [1.0_dp, 2.0_dp, 3.0_dp], 1, position)
   call pl_summary(journal%amount, journal%price, summary)
   call pl_path(journal%amount, journal%price, path)

   print '(a,*(f8.2,1x))', 'Position: ', position(:,1)
   print '(a,f10.2)', 'Closed-trade P/L: ', summary%total_pl
   print '(a,*(f8.2,1x))', 'Realized path: ', path%realized
end program journal_and_pl
