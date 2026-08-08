! SPDX-License-Identifier: LGPL-3.0-only
module pso_random
   use pso_kinds, only : dp
   implicit none
   private
   public :: seed_random, random_box_matrix, random_uniform_vector
   public :: random_r_sphere_like, random_r_spheres_like, random_permutation

contains

   subroutine seed_random(seed)
      integer, intent(in) :: seed
      integer :: n, i
      integer, allocatable :: put(:)

      call random_seed(size=n)
      allocate(put(n))
      do i = 1, n
         put(i) = modulo(seed + 104729 * (i - 1), huge(0) - 1)
         if (put(i) <= 0) put(i) = i
      end do
      call random_seed(put=put)
   end subroutine seed_random

   subroutine random_box_matrix(x, lower, upper)
      real(dp), intent(out) :: x(:,:)
      real(dp), intent(in) :: lower(:), upper(:)
      integer :: j

      call random_number(x)
      do j = 1, size(x, 2)
         x(:,j) = lower + x(:,j) * (upper - lower)
      end do
   end subroutine random_box_matrix

   subroutine random_uniform_vector(x, lower, upper)
      real(dp), intent(out) :: x(:)
      real(dp), intent(in) :: lower, upper
      call random_number(x)
      x = lower + x * (upper - lower)
   end subroutine random_uniform_vector

   subroutine random_r_sphere_like(x, radius)
      real(dp), intent(out) :: x(:)
      real(dp), intent(in) :: radius
      real(dp) :: nrm, u

      call random_number(x)
      nrm = sqrt(sum(x * x))
      call random_number(u)
      if (nrm > 0.0_dp) then
         x = (radius * u / nrm) * x
      else
         x = 0.0_dp
      end if
   end subroutine random_r_sphere_like

   subroutine random_r_spheres_like(x, radius)
      real(dp), intent(out) :: x(:,:)
      real(dp), intent(in) :: radius(:)
      real(dp) :: nrm, u
      integer :: j

      call random_number(x)
      do j = 1, size(x, 2)
         nrm = sqrt(sum(x(:,j) * x(:,j)))
         call random_number(u)
         if (nrm > 0.0_dp) then
            x(:,j) = (radius(j) * u / nrm) * x(:,j)
         else
            x(:,j) = 0.0_dp
         end if
      end do
   end subroutine random_r_spheres_like

   subroutine random_permutation(index)
      integer, intent(out) :: index(:)
      integer :: i, j, tmp
      real(dp) :: u

      do i = 1, size(index)
         index(i) = i
      end do
      do i = size(index), 2, -1
         call random_number(u)
         j = 1 + int(u * real(i, dp))
         if (j > i) j = i
         tmp = index(i)
         index(i) = index(j)
         index(j) = tmp
      end do
   end subroutine random_permutation

end module pso_random
