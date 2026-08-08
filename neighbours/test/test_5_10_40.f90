program test_5_10_40
   use neighbours, only: rng_state, rng_seed, portfolio_5_10_40_neighbour, &
      neighbours_ok
   use neighbours_kinds, only: dp, i8
   implicit none
   type(rng_state) :: rng
   real(dp) :: x(20), xn(20), s
   integer :: i, st

   call rng_seed(rng, 56789_i8)
   x = 0.05_dp
   s = sum(x)
   do i = 1, 1000
      call portfolio_5_10_40_neighbour(x, xn, rng, status=st)
      call check(st == neighbours_ok, '5/10/40 status')
      call check(abs(sum(xn)-s) < 5.0e-14_dp, '5/10/40 sum')
      call check(all(xn >= 0.0_dp) .and. all(xn <= 0.1_dp + 1.0e-15_dp), &
         '5/10/40 individual bounds')
      call check(sum(pack(xn, xn > 0.05_dp)) <= 0.4_dp + 1.0e-12_dp, &
         '5/10/40 concentration')
      x = xn
   end do
   print *, 'test_5_10_40: PASS'
contains
   subroutine check(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         print *, 'FAIL: ', trim(message)
         error stop 1
      end if
   end subroutine check
end program test_5_10_40
