! SPDX-License-Identifier: GPL-3.0-or-later
module cec2013_functions
   use cec2013_kinds, only : dp
   use cec2013_data, only : cec2013_context, CEC2013_OK, CEC2013_BAD_PROBLEM, CEC2013_BAD_SHAPE
   implicit none
   private

   real(dp), parameter :: PI = 3.1415926535897932384626433832795029_dp
   real(dp), parameter :: E_CONST = 2.7182818284590452353602874713526625_dp
   real(dp), parameter :: INF_WEIGHT = 1.0e99_dp

   public :: cec2013_evaluate, cec2013_evaluate_batch, cec2013_optimum_value

contains

   function cec2013_evaluate(ctx, problem, x, status) result(f)
      type(cec2013_context), intent(in) :: ctx
      integer, intent(in) :: problem
      real(dp), intent(in) :: x(:)
      integer, intent(out), optional :: status
      real(dp) :: f
      real(dp), allocatable :: y(:), z(:)

      if (present(status)) status = CEC2013_OK
      f = huge(1.0_dp)
      if (size(x) /= ctx%n .or. ctx%n <= 0) then
         if (present(status)) status = CEC2013_BAD_SHAPE
         return
      end if
      if (problem < 1 .or. problem > 28) then
         if (present(status)) status = CEC2013_BAD_PROBLEM
         return
      end if

      allocate(y(ctx%n), z(ctx%n))
      y = 0.0_dp
      z = 0.0_dp
      call evaluate_core(problem, x, f, ctx%n, ctx%shift, ctx%rotation, y, z)
   end function cec2013_evaluate

   subroutine cec2013_evaluate_batch(ctx, problem, x, f, status)
      type(cec2013_context), intent(in) :: ctx
      integer, intent(in) :: problem
      real(dp), intent(in) :: x(:, :)
      real(dp), intent(out) :: f(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: y(:), z(:)
      integer :: i

      if (present(status)) status = CEC2013_OK
      if (size(x, 1) /= ctx%n .or. size(x, 2) /= size(f) .or. ctx%n <= 0) then
         if (present(status)) status = CEC2013_BAD_SHAPE
         f = huge(1.0_dp)
         return
      end if
      if (problem < 1 .or. problem > 28) then
         if (present(status)) status = CEC2013_BAD_PROBLEM
         f = huge(1.0_dp)
         return
      end if

      allocate(y(ctx%n), z(ctx%n))
      y = 0.0_dp
      z = 0.0_dp
      do i = 1, size(f)
         call evaluate_core(problem, x(:, i), f(i), ctx%n, ctx%shift, ctx%rotation, y, z)
      end do
   end subroutine cec2013_evaluate_batch

   pure real(dp) function cec2013_optimum_value(problem) result(fopt)
      integer, intent(in) :: problem
      if (problem >= 1 .and. problem <= 14) then
         fopt = -1500.0_dp + 100.0_dp*real(problem, dp)
      else if (problem >= 15 .and. problem <= 28) then
         fopt = -1400.0_dp + 100.0_dp*real(problem, dp)
      else
         fopt = huge(1.0_dp)
      end if
   end function cec2013_optimum_value

   subroutine evaluate_core(problem, x, f, n, os, mr, y, z)
      integer, intent(in) :: problem, n
      real(dp), intent(in) :: x(n), os(:), mr(:)
      real(dp), intent(inout) :: y(n), z(n)
      real(dp), intent(out) :: f

      select case (problem)
      case (1)
         call sphere_func(x, f, n, os, 0, mr, 0, .false., y, z)
         f = f - 1400.0_dp
      case (2)
         call ellips_func(x, f, n, os, 0, mr, 0, .true., y, z)
         f = f - 1300.0_dp
      case (3)
         call bent_cigar_func(x, f, n, os, 0, mr, 0, .true., y, z)
         f = f - 1200.0_dp
      case (4)
         call discus_func(x, f, n, os, 0, mr, 0, .true., y, z)
         f = f - 1100.0_dp
      case (5)
         call dif_powers_func(x, f, n, os, 0, mr, 0, .false., y, z)
         f = f - 1000.0_dp
      case (6)
         call rosenbrock_func(x, f, n, os, 0, mr, 0, .true., y, z)
         f = f - 900.0_dp
      case (7)
         call schaffer_f7_func(x, f, n, os, 0, mr, 0, .true., y, z)
         f = f - 800.0_dp
      case (8)
         call ackley_func(x, f, n, os, 0, mr, 0, .true., y, z)
         f = f - 700.0_dp
      case (9)
         call weierstrass_func(x, f, n, os, 0, mr, 0, .true., y, z)
         f = f - 600.0_dp
      case (10)
         call griewank_func(x, f, n, os, 0, mr, 0, .true., y, z)
         f = f - 500.0_dp
      case (11)
         call rastrigin_func(x, f, n, os, 0, mr, 0, .false., y, z)
         f = f - 400.0_dp
      case (12)
         call rastrigin_func(x, f, n, os, 0, mr, 0, .true., y, z)
         f = f - 300.0_dp
      case (13)
         call step_rastrigin_func(x, f, n, os, 0, mr, 0, .true., y, z)
         f = f - 200.0_dp
      case (14)
         call schwefel_func(x, f, n, os, 0, mr, 0, .false., y, z)
         f = f - 100.0_dp
      case (15)
         call schwefel_func(x, f, n, os, 0, mr, 0, .true., y, z)
         f = f + 100.0_dp
      case (16)
         call katsuura_func(x, f, n, os, 0, mr, 0, .true., y, z)
         f = f + 200.0_dp
      case (17)
         call bi_rastrigin_func(x, f, n, os, 0, mr, 0, .false., y, z)
         f = f + 300.0_dp
      case (18)
         call bi_rastrigin_func(x, f, n, os, 0, mr, 0, .true., y, z)
         f = f + 400.0_dp
      case (19)
         call grie_rosen_func(x, f, n, os, 0, mr, 0, .true., y, z)
         f = f + 500.0_dp
      case (20)
         call escaffer6_func(x, f, n, os, 0, mr, 0, .true., y, z)
         f = f + 600.0_dp
      case (21)
         call cf01(x, f, n, os, mr, .true., y, z)
         f = f + 700.0_dp
      case (22)
         call cf02(x, f, n, os, mr, .false., y, z)
         f = f + 800.0_dp
      case (23)
         call cf03(x, f, n, os, mr, .true., y, z)
         f = f + 900.0_dp
      case (24)
         call cf04(x, f, n, os, mr, .true., y, z)
         f = f + 1000.0_dp
      case (25)
         call cf05(x, f, n, os, mr, .true., y, z)
         f = f + 1100.0_dp
      case (26)
         call cf06(x, f, n, os, mr, .true., y, z)
         f = f + 1200.0_dp
      case (27)
         call cf07(x, f, n, os, mr, .true., y, z)
         f = f + 1300.0_dp
      case (28)
         call cf08(x, f, n, os, mr, .true., y, z)
         f = f + 1400.0_dp
      end select
   end subroutine evaluate_core

   subroutine sphere_func(x, f, n, os, osoff, mr, moff, rotate, y, z)
      integer, intent(in) :: n, osoff, moff
      logical, intent(in) :: rotate
      real(dp), intent(in) :: x(n), os(:), mr(:)
      real(dp), intent(inout) :: y(n), z(n)
      real(dp), intent(out) :: f
      call shiftfunc(x, y, n, os, osoff)
      if (rotate) then
         call rotatefunc(y, z, n, mr, moff)
      else
         z = y
      end if
      f = sum(z*z)
   end subroutine sphere_func

   subroutine ellips_func(x, f, n, os, osoff, mr, moff, rotate, y, z)
      integer, intent(in) :: n, osoff, moff
      logical, intent(in) :: rotate
      real(dp), intent(in) :: x(n), os(:), mr(:)
      real(dp), intent(inout) :: y(n), z(n)
      real(dp), intent(out) :: f
      integer :: i
      call shiftfunc(x, y, n, os, osoff)
      if (rotate) then
         call rotatefunc(y, z, n, mr, moff)
      else
         z = y
      end if
      call oszfunc(z, y, n)
      f = 0.0_dp
      do i = 1, n
         f = f + 10.0_dp**(6.0_dp*real(i-1,dp)/real(n-1,dp))*y(i)*y(i)
      end do
   end subroutine ellips_func

   subroutine bent_cigar_func(x, f, n, os, osoff, mr, moff, rotate, y, z)
      integer, intent(in) :: n, osoff, moff
      logical, intent(in) :: rotate
      real(dp), intent(in) :: x(n), os(:), mr(:)
      real(dp), intent(inout) :: y(n), z(n)
      real(dp), intent(out) :: f
      call shiftfunc(x, y, n, os, osoff)
      if (rotate) then
         call rotatefunc(y, z, n, mr, moff)
      else
         z = y
      end if
      call asyfunc(z, y, n, 0.5_dp)
      if (rotate) then
         call rotatefunc(y, z, n, mr, moff+n*n)
      else
         z = y
      end if
      f = z(1)*z(1) + 1.0e6_dp*sum(z(2:n)*z(2:n))
   end subroutine bent_cigar_func

   subroutine discus_func(x, f, n, os, osoff, mr, moff, rotate, y, z)
      integer, intent(in) :: n, osoff, moff
      logical, intent(in) :: rotate
      real(dp), intent(in) :: x(n), os(:), mr(:)
      real(dp), intent(inout) :: y(n), z(n)
      real(dp), intent(out) :: f
      call shiftfunc(x, y, n, os, osoff)
      if (rotate) then
         call rotatefunc(y, z, n, mr, moff)
      else
         z = y
      end if
      call oszfunc(z, y, n)
      f = 1.0e6_dp*y(1)*y(1) + sum(y(2:n)*y(2:n))
   end subroutine discus_func

   subroutine dif_powers_func(x, f, n, os, osoff, mr, moff, rotate, y, z)
      integer, intent(in) :: n, osoff, moff
      logical, intent(in) :: rotate
      real(dp), intent(in) :: x(n), os(:), mr(:)
      real(dp), intent(inout) :: y(n), z(n)
      real(dp), intent(out) :: f
      integer :: i
      call shiftfunc(x, y, n, os, osoff)
      if (rotate) then
         call rotatefunc(y, z, n, mr, moff)
      else
         z = y
      end if
      f = 0.0_dp
      do i = 1, n
         f = f + abs(z(i))**real(2 + (4*(i-1))/(n-1), dp)
      end do
      f = sqrt(f)
   end subroutine dif_powers_func

   subroutine rosenbrock_func(x, f, n, os, osoff, mr, moff, rotate, y, z)
      integer, intent(in) :: n, osoff, moff
      logical, intent(in) :: rotate
      real(dp), intent(in) :: x(n), os(:), mr(:)
      real(dp), intent(inout) :: y(n), z(n)
      real(dp), intent(out) :: f
      integer :: i
      real(dp) :: tmp1, tmp2
      call shiftfunc(x, y, n, os, osoff)
      y = y*2.048_dp/100.0_dp
      if (rotate) then
         call rotatefunc(y, z, n, mr, moff)
      else
         z = y
      end if
      z = z + 1.0_dp
      f = 0.0_dp
      do i = 1, n-1
         tmp1 = z(i)*z(i) - z(i+1)
         tmp2 = z(i) - 1.0_dp
         f = f + 100.0_dp*tmp1*tmp1 + tmp2*tmp2
      end do
   end subroutine rosenbrock_func

   subroutine schaffer_f7_func(x, f, n, os, osoff, mr, moff, rotate, y, z)
      integer, intent(in) :: n, osoff, moff
      logical, intent(in) :: rotate
      real(dp), intent(in) :: x(n), os(:), mr(:)
      real(dp), intent(inout) :: y(n), z(n)
      real(dp), intent(out) :: f
      integer :: i
      real(dp) :: tmp
      call shiftfunc(x, y, n, os, osoff)
      if (rotate) then
         call rotatefunc(y, z, n, mr, moff)
      else
         z = y
      end if
      call asyfunc(z, y, n, 0.5_dp)
      do i = 1, n
         z(i) = y(i)*10.0_dp**(real(i-1,dp)/real(n-1,dp)/2.0_dp)
      end do
      if (rotate) then
         call rotatefunc(z, y, n, mr, moff+n*n)
      else
         y = z
      end if
      do i = 1, n-1
         z(i) = sqrt(y(i)*y(i) + y(i+1)*y(i+1))
      end do
      f = 0.0_dp
      do i = 1, n-1
         tmp = sin(50.0_dp*z(i)**0.2_dp)
         f = f + sqrt(z(i))*(1.0_dp + tmp*tmp)
      end do
      f = f*f/real((n-1)*(n-1), dp)
   end subroutine schaffer_f7_func

   subroutine ackley_func(x, f, n, os, osoff, mr, moff, rotate, y, z)
      integer, intent(in) :: n, osoff, moff
      logical, intent(in) :: rotate
      real(dp), intent(in) :: x(n), os(:), mr(:)
      real(dp), intent(inout) :: y(n), z(n)
      real(dp), intent(out) :: f
      integer :: i
      real(dp) :: sum1, sum2
      call shiftfunc(x, y, n, os, osoff)
      if (rotate) then
         call rotatefunc(y, z, n, mr, moff)
      else
         z = y
      end if
      call asyfunc(z, y, n, 0.5_dp)
      do i = 1, n
         z(i) = y(i)*10.0_dp**(real(i-1,dp)/real(n-1,dp)/2.0_dp)
      end do
      if (rotate) then
         call rotatefunc(z, y, n, mr, moff+n*n)
      else
         y = z
      end if
      sum1 = sum(y*y)
      sum2 = sum(cos(2.0_dp*PI*y))
      sum1 = -0.2_dp*sqrt(sum1/real(n,dp))
      sum2 = sum2/real(n,dp)
      f = E_CONST - 20.0_dp*exp(sum1) - exp(sum2) + 20.0_dp
   end subroutine ackley_func

   subroutine weierstrass_func(x, f, n, os, osoff, mr, moff, rotate, y, z)
      integer, intent(in) :: n, osoff, moff
      logical, intent(in) :: rotate
      real(dp), intent(in) :: x(n), os(:), mr(:)
      real(dp), intent(inout) :: y(n), z(n)
      real(dp), intent(out) :: f
      integer :: i, j
      real(dp) :: sum1, sum2, a, b
      call shiftfunc(x, y, n, os, osoff)
      y = y*0.5_dp/100.0_dp
      if (rotate) then
         call rotatefunc(y, z, n, mr, moff)
      else
         z = y
      end if
      call asyfunc(z, y, n, 0.5_dp)
      do i = 1, n
         z(i) = y(i)*10.0_dp**(real(i-1,dp)/real(n-1,dp)/2.0_dp)
      end do
      if (rotate) then
         call rotatefunc(z, y, n, mr, moff+n*n)
      else
         y = z
      end if
      a = 0.5_dp
      b = 3.0_dp
      f = 0.0_dp
      sum2 = 0.0_dp
      do i = 1, n
         sum1 = 0.0_dp
         sum2 = 0.0_dp
         do j = 0, 20
            sum1 = sum1 + a**j*cos(2.0_dp*PI*b**j*(y(i)+0.5_dp))
            sum2 = sum2 + a**j*cos(2.0_dp*PI*b**j*0.5_dp)
         end do
         f = f + sum1
      end do
      f = f - real(n,dp)*sum2
   end subroutine weierstrass_func

   subroutine griewank_func(x, f, n, os, osoff, mr, moff, rotate, y, z)
      integer, intent(in) :: n, osoff, moff
      logical, intent(in) :: rotate
      real(dp), intent(in) :: x(n), os(:), mr(:)
      real(dp), intent(inout) :: y(n), z(n)
      real(dp), intent(out) :: f
      integer :: i
      real(dp) :: s, p
      call shiftfunc(x, y, n, os, osoff)
      y = y*6.0_dp
      if (rotate) then
         call rotatefunc(y, z, n, mr, moff)
      else
         z = y
      end if
      do i = 1, n
         z(i) = z(i)*100.0_dp**(real(i-1,dp)/real(n-1,dp)/2.0_dp)
      end do
      s = sum(z*z)
      p = 1.0_dp
      do i = 1, n
         p = p*cos(z(i)/sqrt(real(i,dp)))
      end do
      f = 1.0_dp + s/4000.0_dp - p
   end subroutine griewank_func

   subroutine rastrigin_func(x, f, n, os, osoff, mr, moff, rotate, y, z)
      integer, intent(in) :: n, osoff, moff
      logical, intent(in) :: rotate
      real(dp), intent(in) :: x(n), os(:), mr(:)
      real(dp), intent(inout) :: y(n), z(n)
      real(dp), intent(out) :: f
      integer :: i
      call shiftfunc(x, y, n, os, osoff)
      y = y*5.12_dp/100.0_dp
      if (rotate) then
         call rotatefunc(y, z, n, mr, moff)
      else
         z = y
      end if
      call oszfunc(z, y, n)
      call asyfunc(y, z, n, 0.2_dp)
      if (rotate) then
         call rotatefunc(z, y, n, mr, moff+n*n)
      else
         y = z
      end if
      do i = 1, n
         y(i) = y(i)*10.0_dp**(real(i-1,dp)/real(n-1,dp)/2.0_dp)
      end do
      if (rotate) then
         call rotatefunc(y, z, n, mr, moff)
      else
         z = y
      end if
      f = sum(z*z - 10.0_dp*cos(2.0_dp*PI*z) + 10.0_dp)
   end subroutine rastrigin_func

   subroutine step_rastrigin_func(x, f, n, os, osoff, mr, moff, rotate, y, z)
      integer, intent(in) :: n, osoff, moff
      logical, intent(in) :: rotate
      real(dp), intent(in) :: x(n), os(:), mr(:)
      real(dp), intent(inout) :: y(n), z(n)
      real(dp), intent(out) :: f
      integer :: i
      call shiftfunc(x, y, n, os, osoff)
      y = y*5.12_dp/100.0_dp
      if (rotate) then
         call rotatefunc(y, z, n, mr, moff)
      else
         z = y
      end if
      do i = 1, n
         if (abs(z(i)) > 0.5_dp) z(i) = floor_as_real(2.0_dp*z(i)+0.5_dp)/2.0_dp
      end do
      call oszfunc(z, y, n)
      call asyfunc(y, z, n, 0.2_dp)
      if (rotate) then
         call rotatefunc(z, y, n, mr, moff+n*n)
      else
         y = z
      end if
      do i = 1, n
         y(i) = y(i)*10.0_dp**(real(i-1,dp)/real(n-1,dp)/2.0_dp)
      end do
      if (rotate) then
         call rotatefunc(y, z, n, mr, moff)
      else
         z = y
      end if
      f = sum(z*z - 10.0_dp*cos(2.0_dp*PI*z) + 10.0_dp)
   end subroutine step_rastrigin_func

   subroutine schwefel_func(x, f, n, os, osoff, mr, moff, rotate, y, z)
      integer, intent(in) :: n, osoff, moff
      logical, intent(in) :: rotate
      real(dp), intent(in) :: x(n), os(:), mr(:)
      real(dp), intent(inout) :: y(n), z(n)
      real(dp), intent(out) :: f
      integer :: i
      real(dp) :: tmp
      call shiftfunc(x, y, n, os, osoff)
      y = y*10.0_dp
      if (rotate) then
         call rotatefunc(y, z, n, mr, moff)
      else
         z = y
      end if
      do i = 1, n
         y(i) = z(i)*10.0_dp**(real(i-1,dp)/real(n-1,dp)/2.0_dp)
      end do
      z = y + 420.9687462275036_dp
      f = 0.0_dp
      do i = 1, n
         if (z(i) > 500.0_dp) then
            f = f - (500.0_dp-mod(z(i),500.0_dp))*sin(sqrt(500.0_dp-mod(z(i),500.0_dp)))
            tmp = (z(i)-500.0_dp)/100.0_dp
            f = f + tmp*tmp/real(n,dp)
         else if (z(i) < -500.0_dp) then
            f = f - (-500.0_dp+mod(abs(z(i)),500.0_dp))*sin(sqrt(500.0_dp-mod(abs(z(i)),500.0_dp)))
            tmp = (z(i)+500.0_dp)/100.0_dp
            f = f + tmp*tmp/real(n,dp)
         else
            f = f - z(i)*sin(sqrt(abs(z(i))))
         end if
      end do
      f = 418.9828872724338_dp*real(n,dp) + f
   end subroutine schwefel_func

   subroutine katsuura_func(x, f, n, os, osoff, mr, moff, rotate, y, z)
      integer, intent(in) :: n, osoff, moff
      logical, intent(in) :: rotate
      real(dp), intent(in) :: x(n), os(:), mr(:)
      real(dp), intent(inout) :: y(n), z(n)
      real(dp), intent(out) :: f
      integer :: i, j
      real(dp) :: temp, tmp1, tmp2, tmp3
      tmp3 = real(n,dp)**1.2_dp
      call shiftfunc(x, y, n, os, osoff)
      y = y*0.05_dp
      if (rotate) then
         call rotatefunc(y, z, n, mr, moff)
      else
         z = y
      end if
      do i = 1, n
         z(i) = z(i)*100.0_dp**(real(i-1,dp)/real(n-1,dp)/2.0_dp)
      end do
      if (rotate) then
         call rotatefunc(z, y, n, mr, moff+n*n)
      else
         y = z
      end if
      f = 1.0_dp
      do i = 1, n
         temp = 0.0_dp
         do j = 1, 32
            tmp1 = 2.0_dp**j
            tmp2 = tmp1*y(i)
            temp = temp + abs(tmp2-floor_as_real(tmp2+0.5_dp))/tmp1
         end do
         f = f*(1.0_dp + real(i,dp)*temp)**(10.0_dp/tmp3)
      end do
      tmp1 = 10.0_dp/real(n*n,dp)
      f = f*tmp1 - tmp1
   end subroutine katsuura_func

   subroutine bi_rastrigin_func(x, f, n, os, osoff, mr, moff, rotate, y, z)
      integer, intent(in) :: n, osoff, moff
      logical, intent(in) :: rotate
      real(dp), intent(in) :: x(n), os(:), mr(:)
      real(dp), intent(inout) :: y(n), z(n)
      real(dp), intent(out) :: f
      real(dp), allocatable :: tmpx(:)
      real(dp) :: mu0, d, s, mu1, tmp, tmp1, tmp2
      integer :: i
      allocate(tmpx(n))
      mu0 = 2.5_dp
      d = 1.0_dp
      s = 1.0_dp - 1.0_dp/(2.0_dp*sqrt(real(n,dp)+20.0_dp)-8.2_dp)
      mu1 = -sqrt((mu0*mu0-d)/s)
      call shiftfunc(x, y, n, os, osoff)
      y = y*0.1_dp
      do i = 1, n
         tmpx(i) = 2.0_dp*y(i)
         if (os(osoff+i) < 0.0_dp) tmpx(i) = -tmpx(i)
      end do
      z = tmpx
      tmpx = tmpx + mu0
      if (rotate) then
         call rotatefunc(z, y, n, mr, moff)
      else
         y = z
      end if
      do i = 1, n
         y(i) = y(i)*100.0_dp**(real(i-1,dp)/real(n-1,dp)/2.0_dp)
      end do
      if (rotate) then
         call rotatefunc(y, z, n, mr, moff+n*n)
      else
         z = y
      end if
      tmp1 = sum((tmpx-mu0)**2)
      tmp2 = s*sum((tmpx-mu1)**2) + d*real(n,dp)
      tmp = sum(cos(2.0_dp*PI*z))
      f = min(tmp1,tmp2) + 10.0_dp*(real(n,dp)-tmp)
   end subroutine bi_rastrigin_func

   subroutine grie_rosen_func(x, f, n, os, osoff, mr, moff, rotate, y, z)
      integer, intent(in) :: n, osoff, moff
      logical, intent(in) :: rotate
      real(dp), intent(in) :: x(n), os(:), mr(:)
      real(dp), intent(inout) :: y(n), z(n)
      real(dp), intent(out) :: f
      integer :: i
      real(dp) :: temp, tmp1, tmp2
      call shiftfunc(x, y, n, os, osoff)
      y = y*0.05_dp
      if (rotate) then
         call rotatefunc(y, z, n, mr, moff)
      else
         z = y
      end if
      ! Preserve the upstream C source: it overwrites z with y+1 here,
      ! so the just-computed rotation does not affect this function.
      z = y + 1.0_dp
      f = 0.0_dp
      do i = 1, n-1
         tmp1 = z(i)*z(i) - z(i+1)
         tmp2 = z(i) - 1.0_dp
         temp = 100.0_dp*tmp1*tmp1 + tmp2*tmp2
         f = f + temp*temp/4000.0_dp - cos(temp) + 1.0_dp
      end do
      tmp1 = z(n)*z(n) - z(1)
      tmp2 = z(n) - 1.0_dp
      temp = 100.0_dp*tmp1*tmp1 + tmp2*tmp2
      f = f + temp*temp/4000.0_dp - cos(temp) + 1.0_dp
   end subroutine grie_rosen_func

   subroutine escaffer6_func(x, f, n, os, osoff, mr, moff, rotate, y, z)
      integer, intent(in) :: n, osoff, moff
      logical, intent(in) :: rotate
      real(dp), intent(in) :: x(n), os(:), mr(:)
      real(dp), intent(inout) :: y(n), z(n)
      real(dp), intent(out) :: f
      integer :: i
      real(dp) :: temp1, temp2
      call shiftfunc(x, y, n, os, osoff)
      if (rotate) then
         call rotatefunc(y, z, n, mr, moff)
      else
         z = y
      end if
      call asyfunc(z, y, n, 0.5_dp)
      if (rotate) then
         call rotatefunc(y, z, n, mr, moff+n*n)
      else
         z = y
      end if
      f = 0.0_dp
      do i = 1, n-1
         temp1 = sin(sqrt(z(i)*z(i)+z(i+1)*z(i+1)))**2
         temp2 = 1.0_dp + 0.001_dp*(z(i)*z(i)+z(i+1)*z(i+1))
         f = f + 0.5_dp + (temp1-0.5_dp)/(temp2*temp2)
      end do
      temp1 = sin(sqrt(z(n)*z(n)+z(1)*z(1)))**2
      temp2 = 1.0_dp + 0.001_dp*(z(n)*z(n)+z(1)*z(1))
      f = f + 0.5_dp + (temp1-0.5_dp)/(temp2*temp2)
   end subroutine escaffer6_func

   subroutine cf01(x, f, n, os, mr, rotate, y, z)
      integer, intent(in) :: n
      logical, intent(in) :: rotate
      real(dp), intent(in) :: x(n), os(:), mr(:)
      real(dp), intent(inout) :: y(n), z(n)
      real(dp), intent(out) :: f
      real(dp) :: fit(5)
      real(dp), parameter :: delta(5) = [10.0_dp,20.0_dp,30.0_dp,40.0_dp,50.0_dp]
      real(dp), parameter :: bias(5) = [0.0_dp,100.0_dp,200.0_dp,300.0_dp,400.0_dp]
      call rosenbrock_func(x, fit(1), n, os, 0*n, mr, 0*n*n, rotate, y, z)
      fit(1) = 10000.0_dp*fit(1)/1.0e4_dp
      call dif_powers_func(x, fit(2), n, os, 1*n, mr, 1*n*n, rotate, y, z)
      fit(2) = 10000.0_dp*fit(2)/1.0e10_dp
      call bent_cigar_func(x, fit(3), n, os, 2*n, mr, 2*n*n, rotate, y, z)
      fit(3) = 10000.0_dp*fit(3)/1.0e30_dp
      call discus_func(x, fit(4), n, os, 3*n, mr, 3*n*n, rotate, y, z)
      fit(4) = 10000.0_dp*fit(4)/1.0e10_dp
      call sphere_func(x, fit(5), n, os, 4*n, mr, 4*n*n, .false., y, z)
      fit(5) = 10000.0_dp*fit(5)/1.0e5_dp
      call cf_cal(x, f, n, os, delta, bias, fit)
   end subroutine cf01

   subroutine cf02(x, f, n, os, mr, rotate, y, z)
      integer, intent(in) :: n
      logical, intent(in) :: rotate
      real(dp), intent(in) :: x(n), os(:), mr(:)
      real(dp), intent(inout) :: y(n), z(n)
      real(dp), intent(out) :: f
      real(dp) :: fit(3)
      real(dp), parameter :: delta(3) = [20.0_dp,20.0_dp,20.0_dp]
      real(dp), parameter :: bias(3) = [0.0_dp,100.0_dp,200.0_dp]
      integer :: i
      do i = 1, 3
         call schwefel_func(x, fit(i), n, os, (i-1)*n, mr, (i-1)*n*n, rotate, y, z)
      end do
      call cf_cal(x, f, n, os, delta, bias, fit)
   end subroutine cf02

   subroutine cf03(x, f, n, os, mr, rotate, y, z)
      integer, intent(in) :: n
      logical, intent(in) :: rotate
      real(dp), intent(in) :: x(n), os(:), mr(:)
      real(dp), intent(inout) :: y(n), z(n)
      real(dp), intent(out) :: f
      real(dp) :: fit(3)
      real(dp), parameter :: delta(3) = [20.0_dp,20.0_dp,20.0_dp]
      real(dp), parameter :: bias(3) = [0.0_dp,100.0_dp,200.0_dp]
      integer :: i
      do i = 1, 3
         call schwefel_func(x, fit(i), n, os, (i-1)*n, mr, (i-1)*n*n, rotate, y, z)
      end do
      call cf_cal(x, f, n, os, delta, bias, fit)
   end subroutine cf03

   subroutine cf04(x, f, n, os, mr, rotate, y, z)
      integer, intent(in) :: n
      logical, intent(in) :: rotate
      real(dp), intent(in) :: x(n), os(:), mr(:)
      real(dp), intent(inout) :: y(n), z(n)
      real(dp), intent(out) :: f
      real(dp) :: fit(3)
      real(dp), parameter :: delta(3) = [20.0_dp,20.0_dp,20.0_dp]
      real(dp), parameter :: bias(3) = [0.0_dp,100.0_dp,200.0_dp]
      call schwefel_func(x, fit(1), n, os, 0*n, mr, 0*n*n, rotate, y, z)
      fit(1) = 1000.0_dp*fit(1)/4.0e3_dp
      call rastrigin_func(x, fit(2), n, os, 1*n, mr, 1*n*n, rotate, y, z)
      fit(2) = 1000.0_dp*fit(2)/1.0e3_dp
      call weierstrass_func(x, fit(3), n, os, 2*n, mr, 2*n*n, rotate, y, z)
      fit(3) = 1000.0_dp*fit(3)/400.0_dp
      call cf_cal(x, f, n, os, delta, bias, fit)
   end subroutine cf04

   subroutine cf05(x, f, n, os, mr, rotate, y, z)
      integer, intent(in) :: n
      logical, intent(in) :: rotate
      real(dp), intent(in) :: x(n), os(:), mr(:)
      real(dp), intent(inout) :: y(n), z(n)
      real(dp), intent(out) :: f
      real(dp) :: fit(3)
      real(dp), parameter :: delta(3) = [10.0_dp,30.0_dp,50.0_dp]
      real(dp), parameter :: bias(3) = [0.0_dp,100.0_dp,200.0_dp]
      call schwefel_func(x, fit(1), n, os, 0*n, mr, 0*n*n, rotate, y, z)
      fit(1) = 1000.0_dp*fit(1)/4.0e3_dp
      call rastrigin_func(x, fit(2), n, os, 1*n, mr, 1*n*n, rotate, y, z)
      fit(2) = 1000.0_dp*fit(2)/1.0e3_dp
      call weierstrass_func(x, fit(3), n, os, 2*n, mr, 2*n*n, rotate, y, z)
      fit(3) = 1000.0_dp*fit(3)/400.0_dp
      call cf_cal(x, f, n, os, delta, bias, fit)
   end subroutine cf05

   subroutine cf06(x, f, n, os, mr, rotate, y, z)
      integer, intent(in) :: n
      logical, intent(in) :: rotate
      real(dp), intent(in) :: x(n), os(:), mr(:)
      real(dp), intent(inout) :: y(n), z(n)
      real(dp), intent(out) :: f
      real(dp) :: fit(5)
      real(dp), parameter :: delta(5) = [10.0_dp,10.0_dp,10.0_dp,10.0_dp,10.0_dp]
      real(dp), parameter :: bias(5) = [0.0_dp,100.0_dp,200.0_dp,300.0_dp,400.0_dp]
      call schwefel_func(x, fit(1), n, os, 0*n, mr, 0*n*n, rotate, y, z)
      fit(1) = 1000.0_dp*fit(1)/4.0e3_dp
      call rastrigin_func(x, fit(2), n, os, 1*n, mr, 1*n*n, rotate, y, z)
      fit(2) = 1000.0_dp*fit(2)/1.0e3_dp
      call ellips_func(x, fit(3), n, os, 2*n, mr, 2*n*n, rotate, y, z)
      fit(3) = 1000.0_dp*fit(3)/1.0e10_dp
      call weierstrass_func(x, fit(4), n, os, 3*n, mr, 3*n*n, rotate, y, z)
      fit(4) = 1000.0_dp*fit(4)/400.0_dp
      call griewank_func(x, fit(5), n, os, 4*n, mr, 4*n*n, rotate, y, z)
      fit(5) = 1000.0_dp*fit(5)/100.0_dp
      call cf_cal(x, f, n, os, delta, bias, fit)
   end subroutine cf06

   subroutine cf07(x, f, n, os, mr, rotate, y, z)
      integer, intent(in) :: n
      logical, intent(in) :: rotate
      real(dp), intent(in) :: x(n), os(:), mr(:)
      real(dp), intent(inout) :: y(n), z(n)
      real(dp), intent(out) :: f
      real(dp) :: fit(5)
      real(dp), parameter :: delta(5) = [10.0_dp,10.0_dp,10.0_dp,20.0_dp,20.0_dp]
      real(dp), parameter :: bias(5) = [0.0_dp,100.0_dp,200.0_dp,300.0_dp,400.0_dp]
      call griewank_func(x, fit(1), n, os, 0*n, mr, 0*n*n, rotate, y, z)
      fit(1) = 10000.0_dp*fit(1)/100.0_dp
      call rastrigin_func(x, fit(2), n, os, 1*n, mr, 1*n*n, rotate, y, z)
      fit(2) = 10000.0_dp*fit(2)/1.0e3_dp
      call schwefel_func(x, fit(3), n, os, 2*n, mr, 2*n*n, rotate, y, z)
      fit(3) = 10000.0_dp*fit(3)/4.0e3_dp
      call weierstrass_func(x, fit(4), n, os, 3*n, mr, 3*n*n, rotate, y, z)
      fit(4) = 10000.0_dp*fit(4)/400.0_dp
      call sphere_func(x, fit(5), n, os, 4*n, mr, 4*n*n, .false., y, z)
      fit(5) = 10000.0_dp*fit(5)/1.0e5_dp
      call cf_cal(x, f, n, os, delta, bias, fit)
   end subroutine cf07

   subroutine cf08(x, f, n, os, mr, rotate, y, z)
      integer, intent(in) :: n
      logical, intent(in) :: rotate
      real(dp), intent(in) :: x(n), os(:), mr(:)
      real(dp), intent(inout) :: y(n), z(n)
      real(dp), intent(out) :: f
      real(dp) :: fit(5)
      real(dp), parameter :: delta(5) = [10.0_dp,20.0_dp,30.0_dp,40.0_dp,50.0_dp]
      real(dp), parameter :: bias(5) = [0.0_dp,100.0_dp,200.0_dp,300.0_dp,400.0_dp]
      call grie_rosen_func(x, fit(1), n, os, 0*n, mr, 0*n*n, rotate, y, z)
      fit(1) = 10000.0_dp*fit(1)/4.0e3_dp
      call schaffer_f7_func(x, fit(2), n, os, 1*n, mr, 1*n*n, rotate, y, z)
      fit(2) = 10000.0_dp*fit(2)/4.0e6_dp
      call schwefel_func(x, fit(3), n, os, 2*n, mr, 2*n*n, rotate, y, z)
      fit(3) = 10000.0_dp*fit(3)/4.0e3_dp
      call escaffer6_func(x, fit(4), n, os, 3*n, mr, 3*n*n, rotate, y, z)
      fit(4) = 10000.0_dp*fit(4)/2.0e7_dp
      call sphere_func(x, fit(5), n, os, 4*n, mr, 4*n*n, .false., y, z)
      fit(5) = 10000.0_dp*fit(5)/1.0e5_dp
      call cf_cal(x, f, n, os, delta, bias, fit)
   end subroutine cf08

   pure real(dp) function floor_as_real(x) result(y)
      use iso_fortran_env, only : int64
      real(dp), intent(in) :: x
      y = real(floor(x, kind=int64), dp)
   end function floor_as_real

   subroutine shiftfunc(x, xshift, n, os, osoff)
      integer, intent(in) :: n, osoff
      real(dp), intent(in) :: x(n), os(:)
      real(dp), intent(out) :: xshift(n)
      xshift = x - os(osoff+1:osoff+n)
   end subroutine shiftfunc

   subroutine rotatefunc(x, xrot, n, mr, moff)
      integer, intent(in) :: n, moff
      real(dp), intent(in) :: x(n), mr(:)
      real(dp), intent(out) :: xrot(n)
      integer :: i, j, idx
      do i = 1, n
         xrot(i) = 0.0_dp
         idx = moff + (i-1)*n
         do j = 1, n
            xrot(i) = xrot(i) + x(j)*mr(idx+j)
         end do
      end do
   end subroutine rotatefunc

   subroutine asyfunc(x, xasy, n, beta)
      integer, intent(in) :: n
      real(dp), intent(in) :: x(n), beta
      real(dp), intent(inout) :: xasy(n)
      integer :: i
      ! Preserve the upstream C implementation: for x <= 0, the destination
      ! entry is left unchanged rather than copied from x.
      do i = 1, n
         if (x(i) > 0.0_dp) then
            xasy(i) = x(i)**(1.0_dp + beta*real(i-1,dp)/real(n-1,dp)*sqrt(x(i)))
         end if
      end do
   end subroutine asyfunc

   subroutine oszfunc(x, xosz, n)
      integer, intent(in) :: n
      real(dp), intent(in) :: x(n)
      real(dp), intent(out) :: xosz(n)
      integer :: i, sx
      real(dp) :: c1, c2, xx
      xx = 0.0_dp
      do i = 1, n
         if (i == 1 .or. i == n) then
            if (abs(x(i)) > 0.0_dp) xx = log(abs(x(i)))
            if (x(i) > 0.0_dp) then
               c1 = 10.0_dp
               c2 = 7.9_dp
               sx = 1
            else if (x(i) < 0.0_dp) then
               c1 = 5.5_dp
               c2 = 3.1_dp
               sx = -1
            else
               c1 = 5.5_dp
               c2 = 3.1_dp
               sx = 0
            end if
            xosz(i) = real(sx,dp)*exp(xx + 0.049_dp*(sin(c1*xx)+sin(c2*xx)))
         else
            xosz(i) = x(i)
         end if
      end do
   end subroutine oszfunc

   subroutine cf_cal(x, f, n, os, delta, bias, fit)
      integer, intent(in) :: n
      real(dp), intent(in) :: x(n), os(:), delta(:), bias(:)
      real(dp), intent(inout) :: fit(:)
      real(dp), intent(out) :: f
      real(dp) :: w(size(fit)), dist2, wmax, wsum
      integer :: i
      do i = 1, size(fit)
         fit(i) = fit(i) + bias(i)
         dist2 = sum((x-os((i-1)*n+1:i*n))**2)
         if (dist2 > 0.0_dp) then
            w(i) = sqrt(1.0_dp/dist2)*exp(-dist2/(2.0_dp*real(n,dp)*delta(i)**2))
         else
            w(i) = INF_WEIGHT
         end if
      end do
      wmax = maxval(w)
      wsum = sum(w)
      if (wmax <= 0.0_dp) then
         w = 1.0_dp
         wsum = real(size(w),dp)
      end if
      f = sum(w/wsum*fit)
   end subroutine cf_cal

end module cec2013_functions
