program roommate_ttc_example
   use matchingr
   implicit none
   integer :: rp(3,4), tp(4,4)
   integer, allocatable :: tm(:)
   type(roommate_result_t) :: rr
   rp=reshape([2,3,4,1,3,4,1,2,4,1,2,3],[3,4])
   rr=stable_roommates_preferences(rp)
   print '(a,*(i0,1x))','roommates: ',rr%matching
   tp=reshape([2,3,4,1,4,3,2,1,3,4,2,1,4,2,1,3],[4,4])
   tm=top_trading_cycles_preferences(tp)
   print '(a,*(i0,1x))','TTC: ',tm
end program
