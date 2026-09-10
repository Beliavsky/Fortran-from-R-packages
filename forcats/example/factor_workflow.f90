! SPDX-License-Identifier: MIT
! SPDX-FileComment: Example for the Fortran forcats translation.
program factor_workflow
   !! Demonstrates deterministic factor recoding, lumping, and counting.
   use forcats
   implicit none

   type(factor_type) :: sector, compact
   type(factor_count_type) :: counts
   character(len=:), allocatable :: labels(:)
   integer :: i

   sector = fct([character(len=10) :: "technology", "finance", "technology", &
      "energy", "technology", "utilities"])
   compact = fct_lump_min(sector, 2.0_dp, other_level="other")
   counts = fct_count(compact, sort=.true., prop=.true.)
   labels = counts%factor%values()

   print '(a)', "level count proportion"
   do i = 1, size(counts%count)
      print '(a,1x,i0,1x,f6.3)', trim(labels(i)), counts%count(i), counts%proportion(i)
   end do
end program factor_workflow
