program parity_driver
   use pbs_mod, only : dp, pbs_basis, pbs, PBS_OK
   implicit none
   type(pbs_basis) :: periodic_fit
   type(pbs_basis) :: ordinary_fit
   real(dp), parameter :: xp(7) = [0.0_dp, 0.5_dp, 1.0_dp, 2.0_dp, 3.0_dp, 5.5_dp, 6.0_dp]
   real(dp), parameter :: xo(6) = [-1.0_dp, 0.0_dp, 0.5_dp, 1.0_dp, 2.0_dp, 3.0_dp]
   integer :: i
   integer :: j
   integer :: unit

   periodic_fit = pbs(xp, knots=[1.0_dp, 2.0_dp, 3.0_dp], degree=3, intercept=.true., &
      boundary_knots=[0.0_dp, 6.0_dp], periodic=.true.)
   ordinary_fit = pbs(xo, knots=[1.0_dp, 2.0_dp], degree=2, intercept=.true., &
      boundary_knots=[0.0_dp, 2.5_dp], periodic=.false.)
   if (periodic_fit%status /= PBS_OK .or. ordinary_fit%status /= PBS_OK) error stop 1

   open(newunit=unit, file='fortran_parity.csv', status='replace', action='write')
   write(unit, '(a)') 'case,row,x,column,value'
   do i = 1, size(xp)
      do j = 1, size(periodic_fit%values, 2)
         write(unit, '(a,i0,a,es24.16,a,i0,a,es24.16)') 'periodic,', i, ',', xp(i), ',', j, ',', &
            periodic_fit%values(i, j)
      end do
   end do
   do i = 1, size(xo)
      do j = 1, size(ordinary_fit%values, 2)
         write(unit, '(a,i0,a,es24.16,a,i0,a,es24.16)') 'ordinary,', i, ',', xo(i), ',', j, ',', &
            ordinary_fit%values(i, j)
      end do
   end do
   close(unit)
end program parity_driver
