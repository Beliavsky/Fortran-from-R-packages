program numeric_portfolio
   use neighbours, only: numeric_neighbour_config, rng_state, rng_seed, &
      init_numeric_neighbour, numeric_neighbour, neighbours_ok
   use neighbours_kinds, only: dp, i8
   implicit none
   type(numeric_neighbour_config) :: config
   type(rng_state) :: rng
   real(dp) :: x(5), xn(5)
   integer :: i, status

   x = 0.20_dp
   call rng_seed(rng, 101_i8)
   call init_numeric_neighbour(config, 5, 0.03_dp, lower=[0.0_dp], &
      upper=[0.40_dp], status=status)
   if (status /= neighbours_ok) error stop 1

   do i = 1, 10
      call numeric_neighbour(config, x, xn, rng, status=status)
      if (status /= neighbours_ok) error stop 1
      x = xn
   end do
   print '(a,5f10.6)', 'weights = ', x
   print '(a,f10.6)', 'sum     = ', sum(x)
end program numeric_portfolio
