program test_logical
   use neighbours, only: logical_neighbour_config, rng_state, rng_seed, &
      init_logical_neighbour, logical_neighbour, neighbours_ok
   use neighbours_kinds, only: i8
   implicit none
   type(logical_neighbour_config) :: cfg
   type(rng_state) :: rng
   logical :: x(10), xn(10)
   integer :: i, st

   call rng_seed(rng, 23456_i8)
   x = .false.
   call init_logical_neighbour(cfg, 10, stepsize=3, status=st)
   do i = 1, 500
      call logical_neighbour(cfg, x, xn, rng, status=st)
      call check(st == neighbours_ok, 'toggle status')
      call check(count(x .neqv. xn) == 3, 'toggle three')
      x = xn
   end do

   x = .false.
   x(1:3) = .true.
   call init_logical_neighbour(cfg, 10, stepsize=1, kmin=3, kmax=3, status=st)
   do i = 1, 500
      call logical_neighbour(cfg, x, xn, rng, status=st)
      call check(count(xn) == 3, 'constant cardinality')
      x = xn
   end do

   x = .false.
   call init_logical_neighbour(cfg, 10, stepsize=1, kmin=0, kmax=5, status=st)
   do i = 1, 500
      call logical_neighbour(cfg, x, xn, rng, status=st)
      call check(count(xn) <= 5, 'upper cardinality')
      x = xn
   end do

   x = .false.
   x(1:2) = .true.
   call init_logical_neighbour(cfg, 10, stepsize=1, kmin=2, kmax=2, &
      active=[1,2,3,4,5], status=st)
   do i = 1, 100
      call logical_neighbour(cfg, x, xn, rng, status=st)
      call check(all(xn(6:) .eqv. x(6:)), 'logical active')
      call check(count(xn(1:5)) == 2, 'active fixed cardinality')
      x = xn
   end do

   print *, 'test_logical: PASS'
contains
   subroutine check(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         print *, 'FAIL: ', trim(message)
         error stop 1
      end if
   end subroutine check
end program test_logical
