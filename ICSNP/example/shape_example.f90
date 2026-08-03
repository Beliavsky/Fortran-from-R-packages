! SPDX-License-Identifier: GPL-2.0-or-later
program shape_example
   use icsnp, only : dp, icsnp_ok, tyler_shape, duembgen_shape, symm_huber
   implicit none
   integer, parameter :: n = 100, p = 3
   real(dp) :: x(n, p)
   real(dp), allocatable :: shape(:,:)
   integer :: i, status, iterations

   do i = 1, n
      x(i, 1) = sin(0.17_dp * real(i, dp))
      x(i, 2) = 0.6_dp * x(i, 1) + cos(0.23_dp * real(i, dp))
      x(i, 3) = -0.25_dp * x(i, 1) + sin(0.31_dp * real(i, dp))
   end do
   x(n, :) = [8.0_dp, -6.0_dp, 5.0_dp]

   call tyler_shape(x, shape, status, iterations)
   if (status /= icsnp_ok) error stop 'tyler_shape failed'
   call print_matrix('Tyler shape', shape)

   call duembgen_shape(x, shape, status, iterations)
   if (status /= icsnp_ok) error stop 'duembgen_shape failed'
   call print_matrix('Duembgen shape', shape)

   call symm_huber(x, shape, status, iterations)
   if (status /= icsnp_ok) error stop 'symm_huber failed'
   call print_matrix('Symmetrized Huber scatter', shape)
contains
   subroutine print_matrix(label, a)
      character(len=*), intent(in) :: label
      real(dp), intent(in) :: a(:,:)
      integer :: j
      write(*, '(a)') trim(label)
      do j = 1, size(a, 1)
         write(*, '(*(f12.6,1x))') a(j, :)
      end do
   end subroutine print_matrix
end program shape_example
