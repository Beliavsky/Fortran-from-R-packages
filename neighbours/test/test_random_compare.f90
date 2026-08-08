program test_random_compare
   use neighbours, only: rng_state, rng_seed, random_logical_vector, &
      random_logical_vectors, random_numeric_vector, random_numeric_vectors, &
      compare_logical_vectors, neighbours_ok
   use neighbours_kinds, only: dp, i8
   implicit none
   type(rng_state) :: rng
   logical :: lx(20), lm(20,30), v(5,3)
   real(dp) :: x(20), xm(20,20)
   integer, allocatable :: d(:)
   logical, allocatable :: mask(:,:)
   integer :: st, j

   call rng_seed(rng, 45678_i8)
   call random_logical_vector(lx, rng, 3, 7, st)
   call check(st == neighbours_ok .and. count(lx) >= 3 .and. count(lx) <= 7, &
      'random logical')
   call random_logical_vectors(lm, rng, 5, 5, st)
   do j = 1, size(lm,2)
      call check(count(lm(:,j)) == 5, 'random logical matrix')
   end do

   call random_numeric_vector(x, rng, -2.0_dp, 3.0_dp, 4, 8, st)
   call check(all(x >= -2.0_dp) .and. all(x <= 3.0_dp), 'random numeric bounds')
   call check(count(abs(x) > 0.0_dp) >= 4 .and. count(abs(x) > 0.0_dp) <= 8, &
      'random numeric cardinality')
   call random_numeric_vectors(xm, rng, 0.0_dp, 1.0_dp, 6, 6, st)
   do j = 1, size(xm,2)
      call check(count(abs(xm(:,j)) > 0.0_dp) == 6, 'numeric matrix cardinality')
   end do

   v(:,1) = [.false.,.true.,.false.,.true.,.false.]
   v(:,2) = [.true., .true.,.false.,.true.,.false.]
   v(:,3) = [.true., .false.,.false.,.false.,.false.]
   call compare_logical_vectors(v, d, mask, st)
   call check(all(d == [1,2]), 'logical distances')
   call check(count(mask(:,1)) == 1 .and. count(mask(:,2)) == 2, 'difference mask')

   print *, 'test_random_compare: PASS'
contains
   subroutine check(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         print *, 'FAIL: ', trim(message)
         error stop 1
      end if
   end subroutine check
end program test_random_compare
