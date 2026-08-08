program logical_fixed_k
   use neighbours, only: logical_neighbour_config, rng_state, rng_seed, &
      init_logical_neighbour, logical_neighbour, neighbours_ok
   use neighbours_kinds, only: i8
   implicit none
   type(logical_neighbour_config) :: config
   type(rng_state) :: rng
   logical :: x(8), xn(8)
   integer :: i, status

   x = .false.
   x(1:3) = .true.
   call rng_seed(rng, 202_i8)
   call init_logical_neighbour(config, 8, kmin=3, kmax=3, status=status)
   if (status /= neighbours_ok) error stop 1

   do i = 1, 5
      call logical_neighbour(config, x, xn, rng, status=status)
      if (status /= neighbours_ok) error stop 1
      x = xn
      write (*,'(8l2)') x
   end do
end program logical_fixed_k
