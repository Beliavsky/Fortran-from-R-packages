! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran computational translation of tvgarch 2.4.3.
module tvgarch_transition
   use garchx_kinds, only : dp
   implicit none
   private

   type, public :: tv_spec
      integer, allocatable :: orders(:)
      real(dp) :: intercept = 1.0_dp
      real(dp), allocatable :: sizes(:)
      real(dp), allocatable :: speeds(:)
      real(dp), allocatable :: locations(:)
      integer :: speed_option = 2
   end type tv_spec

   public :: make_tv_spec, tv_transition, tv_component, tv_parameter_count
   public :: pack_tv_parameters, unpack_tv_parameters, combinations_binary
contains
   subroutine make_tv_spec(spec, orders, intercept, sizes, speeds, locations, speed_option)
      type(tv_spec), intent(out) :: spec
      integer, intent(in), optional :: orders(:), speed_option
      real(dp), intent(in), optional :: intercept, sizes(:), speeds(:), locations(:)
      integer :: s, nloc

      if (present(orders)) then
         allocate(spec%orders(size(orders)))
         spec%orders = orders
      else
         allocate(spec%orders(0))
      end if
      s = size(spec%orders)
      nloc = sum(spec%orders)
      allocate(spec%sizes(s), spec%speeds(s), spec%locations(nloc))
      spec%sizes = 0.1_dp
      spec%speeds = 10.0_dp
      spec%locations = 0.5_dp
      if (present(intercept)) spec%intercept = intercept
      if (present(sizes)) then
         if (size(sizes) == s) spec%sizes = sizes
      end if
      if (present(speeds)) then
         if (size(speeds) == s) spec%speeds = speeds
      end if
      if (present(locations)) then
         if (size(locations) == nloc) spec%locations = locations
      end if
      if (present(speed_option)) spec%speed_option = speed_option
   end subroutine make_tv_spec

   pure integer function tv_parameter_count(orders) result(npar)
      integer, intent(in) :: orders(:)
      npar = 1 + 2*size(orders) + sum(orders)
   end function tv_parameter_count

   pure real(dp) function sample_sd(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: xbar
      if (size(x) <= 1) then
         value = 1.0_dp
      else
         xbar = sum(x)/real(size(x), dp)
         value = sqrt(sum((x-xbar)**2)/real(size(x)-1, dp))
         if (value <= epsilon(1.0_dp)) value = 1.0_dp
      end if
   end function sample_sd

   pure elemental real(dp) function logistic_stable(x) result(value)
      real(dp), intent(in) :: x
      if (x >= 0.0_dp) then
         value = 1.0_dp/(1.0_dp+exp(-min(x, 700.0_dp)))
      else
         value = exp(max(x, -700.0_dp))/(1.0_dp+exp(max(x, -700.0_dp)))
      end if
   end function logistic_stable

   subroutine tv_transition(speed, location, xtv, values, status, speed_option)
      real(dp), intent(in) :: speed, location(:), xtv(:)
      real(dp), allocatable, intent(out) :: values(:)
      integer, intent(out) :: status
      integer, intent(in), optional :: speed_option
      integer :: i, j, option
      real(dp) :: adjusted, product_term

      status = 0
      option = 0
      if (present(speed_option)) option = speed_option
      if (size(location) < 1 .or. size(xtv) < 1 .or. option < 0 .or. option > 2) then
         status = 1
         allocate(values(0))
         return
      end if
      adjusted = speed
      if (option == 1) adjusted = speed/sample_sd(xtv)
      if (option == 2) adjusted = exp(min(speed, log(huge(1.0_dp))*0.25_dp))
      allocate(values(size(xtv)))
      do i = 1, size(xtv)
         product_term = 1.0_dp
         do j = 1, size(location)
            product_term = product_term*(xtv(i)-location(j))
         end do
         values(i) = logistic_stable(adjusted*product_term)
      end do
   end subroutine tv_transition

   subroutine tv_component(spec, xtv, g, status)
      type(tv_spec), intent(in) :: spec
      real(dp), intent(in) :: xtv(:)
      real(dp), allocatable, intent(out) :: g(:)
      integer, intent(out) :: status
      integer :: s, j, first, last, local_status
      real(dp), allocatable :: transition(:)

      status = 0
      s = size(spec%orders)
      if (size(spec%sizes) /= s .or. size(spec%speeds) /= s .or. &
          size(spec%locations) /= sum(spec%orders) .or. any(spec%orders < 1)) then
         status = 2
         allocate(g(0))
         return
      end if
      allocate(g(size(xtv)))
      g = spec%intercept
      first = 1
      do j = 1, s
         last = first + spec%orders(j)-1
         call tv_transition(spec%speeds(j), spec%locations(first:last), xtv, transition, &
                            local_status, spec%speed_option)
         if (local_status /= 0) then
            status = local_status
            g = 0.0_dp
            return
         end if
         g = g + spec%sizes(j)*transition
         first = last+1
      end do
      if (any(g <= 0.0_dp)) status = 3
   end subroutine tv_component

   subroutine pack_tv_parameters(spec, par)
      type(tv_spec), intent(in) :: spec
      real(dp), allocatable, intent(out) :: par(:)
      integer :: s
      s = size(spec%orders)
      allocate(par(tv_parameter_count(spec%orders)))
      par(1) = spec%intercept
      if (s > 0) then
         par(2:1+s) = spec%sizes
         par(2+s:1+2*s) = spec%speeds
      end if
      if (sum(spec%orders) > 0) par(2+2*s:) = spec%locations
   end subroutine pack_tv_parameters

   subroutine unpack_tv_parameters(orders, par, speed_option, spec, status)
      integer, intent(in) :: orders(:), speed_option
      real(dp), intent(in) :: par(:)
      type(tv_spec), intent(out) :: spec
      integer, intent(out) :: status
      integer :: s
      status = 0
      s = size(orders)
      if (size(par) /= tv_parameter_count(orders) .or. any(orders < 1)) then
         status = 1
         call make_tv_spec(spec)
         return
      end if
      call make_tv_spec(spec, orders=orders, speed_option=speed_option)
      spec%intercept = par(1)
      if (s > 0) then
         spec%sizes = par(2:1+s)
         spec%speeds = par(2+s:1+2*s)
      end if
      if (sum(orders) > 0) spec%locations = par(2+2*s:)
   end subroutine unpack_tv_parameters

   subroutine combinations_binary(n, binary, status)
      integer, intent(in) :: n
      integer, allocatable, intent(out) :: binary(:, :)
      integer, intent(out) :: status
      integer :: rows, i, j
      status = 0
      if (n < 1 .or. n > 30) then
         status = 1
         allocate(binary(0, 0))
         return
      end if
      rows = 2**n-1
      allocate(binary(rows, n))
      do i = 1, rows
         do j = 1, n
            binary(i, j) = merge(1, 0, btest(i, j-1))
         end do
      end do
   end subroutine combinations_binary
end module tvgarch_transition
