program test_numeric
   use neighbours, only: numeric_neighbour_config, rng_state, rng_seed, &
      init_numeric_neighbour, numeric_neighbour, neighbours_ok, budget_range, budget_none
   use neighbours_kinds, only: dp, i8
   implicit none
   type(numeric_neighbour_config) :: cfg
   type(rng_state) :: rng
   real(dp) :: x(25), xn(25), old_sum, a(6,25), ax(6), ax_exact(6)
   integer :: i, st

   call rng_seed(rng, 12345_i8)
   x = 1.0_dp / 25.0_dp
   call init_numeric_neighbour(cfg, 25, 0.015_dp, lower=[-0.05_dp], &
      upper=[0.05_dp], status=st)
   call check(st == neighbours_ok, 'numeric config')
   old_sum = sum(x)
   do i = 1, 1000
      call numeric_neighbour(cfg, x, xn, rng, status=st)
      call check(st == neighbours_ok, 'numeric neighbour')
      call check(all(xn >= -0.05_dp) .and. all(xn <= 0.05_dp), 'bounds')
      call check(abs(sum(xn) - old_sum) < 5.0e-14_dp, 'zero-sum')
      x = xn
   end do

   call init_numeric_neighbour(cfg, 25, 0.01_dp, lower=[-0.2_dp], &
      upper=[0.2_dp], random_step=.false., budget_mode_in=budget_range, &
      budget=[0.8_dp, 1.2_dp], status=st)
   x = 1.0_dp / 25.0_dp
   do i = 1, 1000
      call numeric_neighbour(cfg, x, xn, rng, status=st)
      call check(st == neighbours_ok, 'range neighbour')
      call check(sum(xn) >= 0.8_dp - 1.0e-12_dp, 'budget lower')
      call check(sum(xn) <= 1.2_dp + 1.0e-12_dp, 'budget upper')
      call check(count(abs(xn - x) > 1.0e-14_dp) <= 1, 'one coordinate')
      x = xn
   end do

   call init_numeric_neighbour(cfg, 25, 0.01_dp, lower=[-0.2_dp], &
      upper=[0.2_dp], budget_mode_in=budget_none, active=[1,2,3,4,5], status=st)
   x = 1.0_dp / 25.0_dp
   do i = 1, 200
      call numeric_neighbour(cfg, x, xn, rng, status=st)
      call check(all(abs(xn(6:) - x(6:)) < 1.0e-15_dp), 'active indices')
      x = xn
   end do

   do i = 1, 25
      a(:,i) = [real(i,dp), real(i*i,dp), 1.0_dp, -real(i,dp), &
                 0.5_dp*real(i,dp), 2.0_dp]
   end do
   call init_numeric_neighbour(cfg, 25, 0.02_dp, lower=[-0.1_dp], &
      upper=[0.1_dp], a=a, status=st)
   x = 1.0_dp / 25.0_dp
   ax = matmul(a, x)
   do i = 1, 100
      call numeric_neighbour(cfg, x, xn, rng, ax=ax, status=st)
      call check(st == neighbours_ok, 'Ax neighbour')
      x = xn
   end do
   ax_exact = matmul(a, x)
   call check(maxval(abs(ax - ax_exact)) < 2.0e-12_dp, 'Ax update')

   print *, 'test_numeric: PASS'
contains
   subroutine check(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         print *, 'FAIL: ', trim(message)
         error stop 1
      end if
   end subroutine check
end program test_numeric
