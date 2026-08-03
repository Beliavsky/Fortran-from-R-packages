! SPDX-License-Identifier: GPL-2.0-only
module gensa_rng
   use gensa_kinds, only : dp, i8
   implicit none
   private

   type, public :: ran2_state
      integer(i8) :: idum = -100377_i8
      integer(i8) :: idum2 = 123456789_i8
      integer(i8) :: iv(32) = 0_i8
      integer(i8) :: iy = 0_i8
      logical :: initialized = .false.
      logical :: has_spare_normal = .false.
      real(dp) :: spare_normal = 0.0_dp
   contains
      procedure :: seed => ran2_seed
      procedure :: uniform => ran2_uniform
      procedure :: normal => ran2_normal
   end type ran2_state

   public :: gensa_visit

contains

   subroutine ran2_seed(self, seed_value)
      class(ran2_state), intent(inout) :: self
      integer(i8), intent(in) :: seed_value

      self%idum = -max(1_i8, abs(seed_value))
      self%idum2 = 123456789_i8
      self%iv = 0_i8
      self%iy = 0_i8
      self%initialized = .false.
      self%has_spare_normal = .false.
      self%spare_normal = 0.0_dp
   end subroutine ran2_seed

   function ran2_uniform(self) result(u)
      class(ran2_state), intent(inout) :: self
      real(dp) :: u
      integer(i8) :: j, k
      integer :: jj
      integer(i8), parameter :: im1 = 2147483563_i8
      integer(i8), parameter :: im2 = 2147483399_i8
      integer(i8), parameter :: imm1 = im1 - 1_i8
      real(dp), parameter :: am = 1.0_dp / real(im1, dp)
      real(dp), parameter :: rnmx = 0.99999988_dp

      if (.not. self%initialized .or. self%idum <= 0_i8) then
         self%idum = max(1_i8, -self%idum)
         self%idum2 = self%idum
         do jj = 40, 33, -1
            k = self%idum / 53668_i8
            self%idum = (self%idum - k * 53668_i8) * 40014_i8 - k * 12211_i8
            if (self%idum < 0_i8) self%idum = self%idum + im1
         end do
         do jj = 32, 1, -1
            k = self%idum / 53668_i8
            self%idum = (self%idum - k * 53668_i8) * 40014_i8 - k * 12211_i8
            if (self%idum < 0_i8) self%idum = self%idum + im1
            self%iv(jj) = self%idum
         end do
         self%iy = self%iv(1)
         self%initialized = .true.
      end if

      k = self%idum / 53668_i8
      self%idum = (self%idum - k * 53668_i8) * 40014_i8 - k * 12211_i8
      if (self%idum < 0_i8) self%idum = self%idum + im1

      k = self%idum2 / 52774_i8
      self%idum2 = (self%idum2 - k * 52774_i8) * 40692_i8 - k * 3791_i8
      if (self%idum2 < 0_i8) self%idum2 = self%idum2 + im2

      j = self%iy / 67108862_i8 + 1_i8
      self%iy = self%iv(int(j)) - self%idum2
      self%iv(int(j)) = self%idum
      if (self%iy < 1_i8) self%iy = self%iy + imm1

      u = min(am * real(self%iy, dp), rnmx)
   end function ran2_uniform

   function ran2_normal(self) result(z)
      class(ran2_state), intent(inout) :: self
      real(dp) :: z
      real(dp) :: x, y, s, scale

      if (self%has_spare_normal) then
         z = self%spare_normal
         self%has_spare_normal = .false.
         return
      end if

      do
         x = 2.0_dp * self%uniform() - 1.0_dp
         y = 2.0_dp * self%uniform() - 1.0_dp
         s = x * x + y * y
         if (s > 0.0_dp .and. s < 1.0_dp) exit
      end do

      scale = sqrt(-2.0_dp * log(s) / s)
      z = y * scale
      self%spare_normal = x * scale
      self%has_spare_normal = .true.
   end function ran2_normal

   function gensa_visit(qv, temperature, rng) result(step)
      real(dp), intent(in) :: qv, temperature
      type(ran2_state), intent(inout) :: rng
      real(dp) :: step
      real(dp) :: pi, factor1, factor2, factor3, factor4
      real(dp) :: factor5, factor6, sigma_x, denominator, y

      pi = acos(-1.0_dp)
      factor1 = exp(log(temperature) / (qv - 1.0_dp))
      factor2 = exp((4.0_dp - qv) * log(qv - 1.0_dp))
      factor3 = exp((2.0_dp - qv) * log(2.0_dp) / (qv - 1.0_dp))
      factor4 = sqrt(pi) * factor1 * factor2 / (factor3 * (3.0_dp - qv))
      factor5 = 1.0_dp / (qv - 1.0_dp) - 0.5_dp
      factor6 = pi * (1.0_dp - factor5) / sin(pi * (1.0_dp - factor5)) &
         / exp(log_gamma(2.0_dp - factor5))
      sigma_x = exp(-(qv - 1.0_dp) * log(factor6 / factor4) / (3.0_dp - qv))

      y = rng%normal()
      denominator = exp((qv - 1.0_dp) * log(max(abs(y), tiny(1.0_dp))) / (3.0_dp - qv))
      step = sigma_x * rng%normal() / denominator
   end function gensa_visit

end module gensa_rng
