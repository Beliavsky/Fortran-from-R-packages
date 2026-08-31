! Gaussian quadrature utilities used by ordinal mixed-model likelihoods.
! Copyright (C) 2011-2026 R. H. B. Christensen; translation (C) 2026.
! Distributed under GPL-2.0-or-later.
module ordinal_quadrature
   use ordinal_kinds, only : dp
   implicit none
   private
   real(dp), parameter :: pi = acos(-1.0_dp)
   public :: gauss_hermite_rule
contains
   pure subroutine gauss_hermite_rule(n, nodes, weights, status)
      integer, intent(in) :: n !! Number of physicists' Gauss-Hermite points; must be positive.
      real(dp), allocatable, intent(out) :: nodes(:) !! Increasing abscissae for integration with weight exp(-x**2).
      real(dp), allocatable, intent(out) :: weights(:) !! Positive weights corresponding to nodes.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid order or root-iteration failure.
      real(dp) :: z, zold, p1, p2, p3, pp
      integer :: i, j, its, m
      if (n < 1) then
         allocate(nodes(0), weights(0))
         status = 1
         return
      end if
      allocate(nodes(n), weights(n))
      nodes = 0.0_dp
      weights = 0.0_dp
      m = (n + 1)/2
      do i = 1, m
         select case (i)
         case (1)
            z = sqrt(2.0_dp*real(n, dp) + 1.0_dp) &
                 - 1.85575_dp*(2.0_dp*real(n, dp) + 1.0_dp)**(-1.0_dp/6.0_dp)
         case (2)
            z = z - 1.14_dp*real(n, dp)**0.426_dp/z
         case (3)
            z = 1.86_dp*z - 0.86_dp*nodes(n)
         case (4)
            z = 1.91_dp*z - 0.91_dp*nodes(n - 1)
         case default
            z = 2.0_dp*z - nodes(n - i + 3)
         end select
         do its = 1, 30
            p1 = pi**(-0.25_dp)
            p2 = 0.0_dp
            do j = 1, n
               p3 = p2
               p2 = p1
               p1 = z*sqrt(2.0_dp/real(j, dp))*p2 &
                    - sqrt(real(j - 1, dp)/real(j, dp))*p3
            end do
            pp = sqrt(2.0_dp*real(n, dp))*p2
            zold = z
            z = zold - p1/pp
            if (abs(z - zold) <= 20.0_dp*epsilon(1.0_dp)*max(1.0_dp, abs(z))) exit
         end do
         if (its > 30) then
            status = 2
            return
         end if
         nodes(i) = -z
         nodes(n + 1 - i) = z
         weights(i) = 2.0_dp/(pp*pp)
         weights(n + 1 - i) = weights(i)
      end do
      status = 0
   end subroutine gauss_hermite_rule
end module ordinal_quadrature
