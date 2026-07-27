! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2008-2021 David Ardia
! Modern Fortran translation of computational routines from bayesGARCH.
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
module bayesgarch_garch
   use bayesgarch_kinds, only : dp
   use bayesgarch_math, only : student_log_density
   use bayesgarch_rng, only : random_standardized_student
   implicit none
   private

   public :: garch11_filter
   public :: garch_filter
   public :: threshold_garch_filter
   public :: simulate_garch11_student
   public :: simulate_garch
   public :: simulate_threshold_garch
   public :: garch11_student_loglik
   public :: filter_alpha
   public :: filter_alpha_asymmetric
   public :: filter_w
   public :: filter_w_asymmetric
   public :: quasi_difference

contains

   subroutine garch11_filter(y, theta, h)
      real(dp), intent(in) :: y(:)
      real(dp), intent(in) :: theta(3)
      real(dp), intent(out) :: h(size(y) + 1)
      integer :: t

      if (any(theta <= 0.0_dp)) error stop "garch11_filter: parameters must be positive"
      h(1) = theta(1)
      do t = 2, size(y) + 1
         h(t) = theta(1) + theta(2) * y(t - 1)**2 + theta(3) * h(t - 1)
      end do
   end subroutine garch11_filter

   subroutine garch_filter(y, omega, alpha, beta, h)
      real(dp), intent(in) :: y(:)
      real(dp), intent(in) :: omega
      real(dp), intent(in) :: alpha(:)
      real(dp), intent(in) :: beta(:)
      real(dp), intent(out) :: h(size(y) + 1)

      call threshold_garch_filter(y, omega, alpha, beta, alpha, h)
   end subroutine garch_filter

   subroutine threshold_garch_filter(y, omega, alpha_positive, beta, alpha_negative, h)
      real(dp), intent(in) :: y(:)
      real(dp), intent(in) :: omega
      real(dp), intent(in) :: alpha_positive(:)
      real(dp), intent(in) :: beta(:)
      real(dp), intent(in) :: alpha_negative(:)
      real(dp), intent(out) :: h(size(y) + 1)
      integer :: i
      integer :: j
      integer :: p
      integer :: q
      integer :: s
      real(dp) :: value

      q = size(alpha_positive)
      p = size(beta)
      if (size(alpha_negative) /= q) error stop "threshold_garch_filter: ARCH arrays differ in size"
      if (omega <= 0.0_dp) error stop "threshold_garch_filter: omega must be positive"
      if (any(alpha_positive < 0.0_dp) .or. any(alpha_negative < 0.0_dp) .or. any(beta < 0.0_dp)) &
         error stop "threshold_garch_filter: coefficients must be nonnegative"
      s = max(p, q)
      if (s > size(y)) error stop "threshold_garch_filter: series is shorter than model order"

      h = omega
      do i = s + 1, size(y) + 1
         value = omega
         do j = 1, q
            if (y(i - j) >= 0.0_dp) then
               value = value + alpha_positive(j) * y(i - j)**2
            else
               value = value + alpha_negative(j) * y(i - j)**2
            end if
         end do
         do j = 1, p
            value = value + beta(j) * h(i - j)
         end do
         h(i) = value
      end do
   end subroutine threshold_garch_filter


   subroutine simulate_garch(innovations, omega, alpha, beta, y, h)
      real(dp), intent(in) :: innovations(:)
      real(dp), intent(in) :: omega
      real(dp), intent(in) :: alpha(:)
      real(dp), intent(in) :: beta(:)
      real(dp), intent(out) :: y(size(innovations))
      real(dp), intent(out) :: h(size(innovations) + 1)

      call simulate_threshold_garch(innovations, omega, alpha, beta, alpha, y, h)
   end subroutine simulate_garch

   subroutine simulate_threshold_garch(innovations, omega, alpha_positive, beta, alpha_negative, y, h)
      real(dp), intent(in) :: innovations(:)
      real(dp), intent(in) :: omega
      real(dp), intent(in) :: alpha_positive(:)
      real(dp), intent(in) :: beta(:)
      real(dp), intent(in) :: alpha_negative(:)
      real(dp), intent(out) :: y(size(innovations))
      real(dp), intent(out) :: h(size(innovations) + 1)
      integer :: i
      integer :: j
      integer :: n
      integer :: p
      integer :: q
      integer :: s
      real(dp) :: value

      n = size(innovations)
      q = size(alpha_positive)
      p = size(beta)
      if (size(alpha_negative) /= q) error stop "simulate_threshold_garch: ARCH arrays differ in size"
      if (omega <= 0.0_dp) error stop "simulate_threshold_garch: omega must be positive"
      if (any(alpha_positive < 0.0_dp) .or. any(alpha_negative < 0.0_dp) .or. any(beta < 0.0_dp)) &
         error stop "simulate_threshold_garch: coefficients must be nonnegative"
      s = max(p, q)
      if (n < max(1, s)) error stop "simulate_threshold_garch: innovations are shorter than model order"

      h = omega
      y = 0.0_dp
      if (s > 0) then
         do i = 1, s
            y(i) = sqrt(h(i)) * innovations(i)
         end do
      end if
      do i = s + 1, n
         value = omega
         do j = 1, q
            if (y(i - j) >= 0.0_dp) then
               value = value + alpha_positive(j) * y(i - j)**2
            else
               value = value + alpha_negative(j) * y(i - j)**2
            end if
         end do
         do j = 1, p
            value = value + beta(j) * h(i - j)
         end do
         h(i) = value
         y(i) = sqrt(value) * innovations(i)
      end do

      value = omega
      do j = 1, q
         if (y(n + 1 - j) >= 0.0_dp) then
            value = value + alpha_positive(j) * y(n + 1 - j)**2
         else
            value = value + alpha_negative(j) * y(n + 1 - j)**2
         end if
      end do
      do j = 1, p
         value = value + beta(j) * h(n + 1 - j)
      end do
      h(n + 1) = value
   end subroutine simulate_threshold_garch

   subroutine simulate_garch11_student(n, theta, nu, y, h, burn)
      integer, intent(in) :: n
      real(dp), intent(in) :: theta(3)
      real(dp), intent(in) :: nu
      real(dp), intent(out) :: y(n)
      real(dp), intent(out) :: h(n + 1)
      integer, intent(in), optional :: burn
      integer :: b
      integer :: i
      integer :: total
      real(dp), allocatable :: all_y(:)
      real(dp), allocatable :: all_h(:)

      if (n < 1) error stop "simulate_garch11_student: n must be positive"
      if (nu <= 2.0_dp) error stop "simulate_garch11_student: nu must exceed 2"
      if (any(theta <= 0.0_dp)) error stop "simulate_garch11_student: parameters must be positive"
      b = 0
      if (present(burn)) b = max(0, burn)
      total = n + b
      allocate(all_y(total), all_h(total + 1))
      all_h(1) = theta(1)
      do i = 1, total
         all_y(i) = sqrt(all_h(i)) * random_standardized_student(nu)
         all_h(i + 1) = theta(1) + theta(2) * all_y(i)**2 + theta(3) * all_h(i)
      end do
      y = all_y(b + 1:b + n)
      h = all_h(b + 1:b + n + 1)
   end subroutine simulate_garch11_student

   function garch11_student_loglik(y, theta, nu) result(value)
      real(dp), intent(in) :: y(:)
      real(dp), intent(in) :: theta(3)
      real(dp), intent(in) :: nu
      real(dp) :: value
      real(dp) :: h(size(y) + 1)
      integer :: i

      call garch11_filter(y, theta, h)
      value = 0.0_dp
      do i = 1, size(y)
         value = value + student_log_density(y(i), h(i), nu)
      end do
   end function garch11_student_loglik

   subroutine filter_alpha(u, beta, x)
      real(dp), intent(in) :: u(:)
      real(dp), intent(in) :: beta
      real(dp), intent(out) :: x(size(u), 2)
      real(dp) :: vstar(size(u))
      real(dp) :: lstar(size(u))
      integer :: i

      vstar(1) = u(1)**2
      lstar(1) = 1.0_dp
      do i = 2, size(u)
         vstar(i) = u(i)**2 + beta * vstar(i - 1)
         lstar(i) = 1.0_dp + beta * lstar(i - 1)
      end do
      x(:, 1) = lstar
      x(1, 2) = 0.0_dp
      if (size(u) > 1) x(2:, 2) = vstar(:size(u) - 1)
   end subroutine filter_alpha

   subroutine filter_alpha_asymmetric(u, beta, x)
      real(dp), intent(in) :: u(:)
      real(dp), intent(in) :: beta
      real(dp), intent(out) :: x(size(u), 3)
      real(dp) :: vpos(size(u))
      real(dp) :: vneg(size(u))
      real(dp) :: lstar(size(u))
      integer :: i

      vpos(1) = u(1)**2
      vneg(1) = u(1)**2
      lstar(1) = 0.0_dp
      do i = 2, size(u)
         vpos(i) = beta * vpos(i - 1)
         if (u(i) >= 0.0_dp) vpos(i) = vpos(i) + u(i)**2
         vneg(i) = beta * vneg(i - 1)
         if (u(i) < 0.0_dp) vneg(i) = vneg(i) + u(i)**2
         lstar(i) = 1.0_dp + beta * lstar(i - 1)
      end do
      x(:, 1) = lstar
      x(1, 2:3) = 0.0_dp
      if (size(u) > 1) then
         x(2:, 2) = vpos(:size(u) - 1)
         x(2:, 3) = vneg(:size(u) - 1)
      end if
   end subroutine filter_alpha_asymmetric

   subroutine filter_w(u, alpha, beta, w)
      real(dp), intent(in) :: u(:)
      real(dp), intent(in) :: alpha(2)
      real(dp), intent(in) :: beta
      real(dp), intent(out) :: w(size(u))
      integer :: i

      w(1) = u(1)**2 - alpha(1)
      do i = 2, size(u)
         w(i) = u(i)**2 - alpha(1) - (alpha(2) + beta) * u(i - 1)**2 + beta * w(i - 1)
      end do
   end subroutine filter_w

   subroutine filter_w_asymmetric(u, alpha, beta, w)
      real(dp), intent(in) :: u(:)
      real(dp), intent(in) :: alpha(3)
      real(dp), intent(in) :: beta
      real(dp), intent(out) :: w(size(u))
      integer :: i

      w(1) = u(1)**2 - alpha(1)
      do i = 2, size(u)
         w(i) = u(i)**2 - alpha(1) + beta * w(i - 1)
         if (u(i - 1) >= 0.0_dp) then
            w(i) = w(i) - (alpha(2) + beta) * u(i - 1)**2
         else
            w(i) = w(i) - (alpha(3) + beta) * u(i - 1)**2
         end if
      end do
   end subroutine filter_w_asymmetric

   subroutine quasi_difference(u, theta, w, w0)
      real(dp), intent(in) :: u(:)
      real(dp), intent(in) :: theta
      real(dp), intent(out) :: w(size(u))
      real(dp), intent(in), optional :: w0
      integer :: i

      w(1) = 0.0_dp
      if (present(w0)) w(1) = w0
      do i = 2, size(u)
         w(i) = u(i - 1) - theta * w(i - 1)
      end do
   end subroutine quasi_difference

end module bayesgarch_garch
