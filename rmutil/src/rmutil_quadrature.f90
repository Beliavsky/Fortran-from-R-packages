! rmutil computational translation
! Copyright (C) 1998-2001 J.K. Lindsey
! Copyright (C) 2026 OpenAI (modern Fortran translation)
! SPDX-License-Identifier: GPL-2.0-or-later
module rmutil_quadrature
   use rmutil_kinds, only : dp, pi
   implicit none
   private
   public :: gauss_hermite
contains
   subroutine gauss_hermite(points, nodes, weights, iterlim, status)
      integer, intent(in) :: points
      real(dp), allocatable, intent(out) :: nodes(:), weights(:)
      integer, intent(in), optional :: iterlim
      integer, intent(out), optional :: status
      integer :: i, j, m, limit
      real(dp) :: z, z1, p, pp
      allocate(nodes(points), weights(points))
      nodes = 0.0_dp
      weights = 0.0_dp
      limit = 20
      if (present(iterlim)) limit = iterlim
      if (present(status)) status = 0
      m = (points+1)/2
      z = 0.0_dp
      do i = 1, m
         if (i == 1) then
            z = sqrt(2.0_dp*real(points,dp)+1.0_dp) - &
               2.0_dp*(2.0_dp*real(points,dp)+1.0_dp)**(-1.0_dp/6.0_dp)
         else if (i == 2) then
            z = z - sqrt(real(points,dp))/z
         else if (i == 3 .or. i == 4) then
            z = 1.9_dp*z - 0.9_dp*nodes(i-2)
         else
            z = 2.0_dp*z - nodes(i-2)
         end if
         do j = 1, limit
            z1 = z
            call hermite_value_derivative(points,z,p,pp)
            z = z1 - p/pp
            if (abs(z-z1) <= 1.0e-15_dp) exit
         end do
         if (j > limit .and. present(status)) status = 1
         nodes(i) = z
         nodes(points+1-i) = -z
         weights(i) = 2.0_dp/(pp*pp)
         weights(points+1-i) = weights(i)
      end do
      nodes = nodes*sqrt(2.0_dp)
      weights = weights/sum(weights)
   end subroutine gauss_hermite

   subroutine hermite_value_derivative(points,z,p,pp)
      integer, intent(in) :: points
      real(dp), intent(in) :: z
      real(dp), intent(out) :: p, pp
      integer :: j
      real(dp) :: p1, p2, p3
      p1 = 1.0_dp/pi**0.4_dp
      p2 = 0.0_dp
      do j = 1, points
         p3 = p2
         p2 = p1
         p1 = z*sqrt(2.0_dp/real(j,dp))*p2 - &
            sqrt(real(j-1,dp)/real(j,dp))*p3
      end do
      p = p1
      pp = sqrt(2.0_dp*real(points,dp))*p2
   end subroutine hermite_value_derivative
end module rmutil_quadrature
