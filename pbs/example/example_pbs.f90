program example_pbs
   use pbs_mod, only : dp, pbs_basis, pbs, PBS_OK
   implicit none
   type(pbs_basis) :: basis
   real(dp) :: x(9)
   integer :: i

   x = [(real(i - 1, dp) * acos(-1.0_dp) / 4.0_dp, i=1, 9)]
   basis = pbs(x, df=5, degree=3, intercept=.false., &
      boundary_knots=[0.0_dp, 2.0_dp * acos(-1.0_dp)], periodic=.true.)
   if (basis%status /= PBS_OK) error stop trim(basis%message)

   print '(a,i0,a,i0)', 'Periodic B-spline basis: ', size(basis%values, 1), ' x ', size(basis%values, 2)
   do i = 1, size(x)
      write (*, '(f8.4,1x,*(f10.6,1x))') x(i), basis%values(i, :)
   end do
end program example_pbs
