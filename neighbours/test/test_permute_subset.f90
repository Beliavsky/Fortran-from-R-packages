program test_permute_subset
   use neighbours, only: rng_state, rng_seed, permute_neighbour, next_subset, &
      neighbours_ok
   use neighbours_kinds, only: dp, i8
   implicit none
   type(rng_state) :: rng
   integer :: xi(5), a(2), st, nseen
   real(dp) :: xr(6), sum_before
   character(len=1) :: xc(5)
   logical :: has_next

   call rng_seed(rng, 34567_i8)
   xi = [1,2,3,4,5]
   call permute_neighbour(xi, rng, 2, st)
   call check(st == neighbours_ok, 'integer permute status')
   call check(sum(xi) == 15, 'integer permutation content')
   call check(count(xi == [1,2,3,4,5]) == 3, 'two-position swap')

   xr = [1,2,3,4,5,6]
   sum_before = sum(xr)
   call permute_neighbour(xr, rng, 5, st)
   call check(abs(sum(xr)-sum_before) < 1.0e-15_dp, 'real permutation')

   xc = ['a','b','c','d','e']
   call permute_neighbour(xc, rng, 5, st)
   call check(all([(count(xc == achar(iachar('a')+nseen)) == 1, nseen=0,4)]), &
      'character permutation')

   a = [1,2]
   nseen = 1
   do
      call next_subset(a, 4, 2, has_next, st)
      call check(st == neighbours_ok, 'next_subset status')
      if (.not. has_next) exit
      nseen = nseen + 1
   end do
   call check(nseen == 6, 'all 4 choose 2 subsets')
   call check(all(a == [3,4]), 'last subset retained')

   print *, 'test_permute_subset: PASS'
contains
   subroutine check(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         print *, 'FAIL: ', trim(message)
         error stop 1
      end if
   end subroutine check
end program test_permute_subset
