! SPDX-License-Identifier: GPL-3.0-only
module fingraph_rng
   use fingraph_kinds, only : dp
   use fingraph_linalg, only : symmetric_eigen_jacobi
   implicit none
   private
   public :: rng_state, seed_rng, uniform_random, normal_random
   public :: gamma_random, random_mvn, random_mvt

   type :: rng_state
      integer(kind=8) :: state = 104729_8
      logical :: has_spare = .false.
      real(dp) :: spare = 0.0_dp
   end type rng_state
contains
   subroutine seed_rng(rng, seed)
      type(rng_state), intent(inout) :: rng
      integer, intent(in) :: seed
      rng%state = max(1_8, int(abs(seed),8))
      rng%has_spare = .false.
   end subroutine seed_rng

   function uniform_random(rng) result(u)
      type(rng_state), intent(inout) :: rng
      real(dp) :: u
      integer(kind=8), parameter :: a = 48271_8, m = 2147483647_8
      rng%state = modulo(a*rng%state, m)
      u = real(rng%state,dp)/real(m,dp)
      u = max(tiny(1.0_dp), min(1.0_dp-epsilon(1.0_dp), u))
   end function uniform_random

   function normal_random(rng) result(z)
      type(rng_state), intent(inout) :: rng
      real(dp) :: z, r, theta
      real(dp), parameter :: two_pi = 6.283185307179586476925286766559_dp
      if (rng%has_spare) then
         z = rng%spare
         rng%has_spare = .false.
      else
         r = sqrt(-2.0_dp*log(uniform_random(rng)))
         theta = two_pi*uniform_random(rng)
         z = r*cos(theta)
         rng%spare = r*sin(theta)
         rng%has_spare = .true.
      end if
   end function normal_random

   recursive function gamma_random(rng, shape) result(g)
      type(rng_state), intent(inout) :: rng
      real(dp), intent(in) :: shape
      real(dp) :: g, d, c, x, v, u
      if (shape <= 0.0_dp) then
         g = 0.0_dp
      else if (shape < 1.0_dp) then
         g = gamma_random(rng,shape+1.0_dp)*uniform_random(rng)**(1.0_dp/shape)
      else
         d = shape - 1.0_dp/3.0_dp
         c = 1.0_dp/sqrt(9.0_dp*d)
         do
            x = normal_random(rng)
            v = 1.0_dp + c*x
            if (v <= 0.0_dp) cycle
            v = v*v*v
            u = uniform_random(rng)
            if (u < 1.0_dp-0.0331_dp*x**4) exit
            if (log(u) < 0.5_dp*x*x+d*(1.0_dp-v+log(v))) exit
         end do
         g = d*v
      end if
   end function gamma_random

   subroutine random_mvn(covariance, nobs, x, seed)
      real(dp), intent(in) :: covariance(:,:)
      integer, intent(in) :: nobs
      real(dp), allocatable, intent(out) :: x(:,:)
      integer, intent(in), optional :: seed
      type(rng_state) :: rng
      real(dp), allocatable :: values(:), vectors(:,:), root(:,:), z(:)
      integer :: i, j, p, info
      p = size(covariance,1)
      allocate(x(nobs,p),z(p))
      call seed_rng(rng,12345)
      if (present(seed)) call seed_rng(rng,seed)
      call symmetric_eigen_jacobi(covariance,values,vectors,info)
      values = sqrt(max(values,0.0_dp))
      root = vectors
      do j = 1,p
         root(:,j) = root(:,j)*values(j)
      end do
      do i = 1,nobs
         do j = 1,p
            z(j) = normal_random(rng)
         end do
         x(i,:) = matmul(root,z)
      end do
   end subroutine random_mvn

   subroutine random_mvt(scatter, nu, nobs, x, seed)
      real(dp), intent(in) :: scatter(:,:), nu
      integer, intent(in) :: nobs
      real(dp), allocatable, intent(out) :: x(:,:)
      integer, intent(in), optional :: seed
      type(rng_state) :: rng
      real(dp), allocatable :: values(:), vectors(:,:), root(:,:), z(:)
      real(dp) :: chi
      integer :: i, j, p, info
      p = size(scatter,1)
      allocate(x(nobs,p),z(p))
      call seed_rng(rng,12345)
      if (present(seed)) call seed_rng(rng,seed)
      call symmetric_eigen_jacobi(scatter,values,vectors,info)
      values = sqrt(max(values,0.0_dp))
      root = vectors
      do j = 1,p
         root(:,j) = root(:,j)*values(j)
      end do
      do i = 1,nobs
         chi = 2.0_dp*gamma_random(rng,0.5_dp*nu)
         do j = 1,p
            z(j) = normal_random(rng)
         end do
         x(i,:) = matmul(root,z)*sqrt(nu/max(chi,tiny(1.0_dp)))
      end do
   end subroutine random_mvt
end module fingraph_rng
