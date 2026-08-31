program test_rank
   use ordinal, only : dp, independent_columns, drop_rank_deficient_columns
   implicit none
   real(dp) :: x(6, 4)
   real(dp), allocatable :: reduced(:, :)
   logical :: keep0(4)
   logical, allocatable :: keep(:)
   integer :: i, rank, status

   do i = 1, 6
      x(i, 1) = 1.0_dp
      x(i, 2) = real(i, dp)
      x(i, 3) = 2.0_dp*x(i, 2)
      x(i, 4) = real((-1)**i, dp)
   end do
   call independent_columns(x, keep0, rank, status=status)
   if (status /= 0 .or. rank /= 3) error stop 'rank detection failed'
   if (count(keep0) /= 3) error stop 'rank keep mask failed'
   call drop_rank_deficient_columns(x, reduced, keep, rank, status=status)
   if (status /= 0 .or. size(reduced, 2) /= 3) error stop 'rank-deficient column drop failed'
   if (count(keep) /= 3) error stop 'drop keep mask failed'

   print *, 'test_rank: PASS'
end program test_rank
