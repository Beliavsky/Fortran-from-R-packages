module ycevo_prediction
   use ycevo_kinds, only : dp
   use ycevo_status, only : ycevo_success, ycevo_err_input
   use ycevo_linalg, only : weighted_quadratic_fit
   use ycevo_curve, only : discount_to_yield
   use ycevo_types, only : yield_surface_t
   implicit none
   private

   public :: linear_interpolate, bilinear_interpolate
   public :: loess_predict, predict_discount, predict_yield

contains

   real(dp) function linear_interpolate(x, y, xout) result(value)
      real(dp), intent(in) :: x(:), y(:), xout
      integer :: i, n
      real(dp) :: t

      n = size(x)
      if (n == 0 .or. size(y) /= n) then
         value = huge(1.0_dp)
         return
      end if
      if (n == 1 .or. xout <= x(1)) then
         value = y(1)
         return
      end if
      if (xout >= x(n)) then
         value = y(n)
         return
      end if
      do i = 1, n - 1
         if (xout >= x(i) .and. xout <= x(i+1)) then
            t = (xout - x(i))/(x(i+1) - x(i))
            value = (1.0_dp - t)*y(i) + t*y(i+1)
            return
         end if
      end do
      value = y(n)
   end function linear_interpolate

   real(dp) function bilinear_interpolate(x, z, values, xout, zout) result(value)
      real(dp), intent(in) :: x(:), z(:), values(:, :), xout, zout
      real(dp), allocatable :: temp(:)
      integer :: j

      if (size(values,1) /= size(x) .or. size(values,2) /= size(z)) then
         value = huge(1.0_dp)
         return
      end if
      allocate(temp(size(z)))
      do j = 1, size(z)
         temp(j) = linear_interpolate(x, values(:,j), xout)
      end do
      value = linear_interpolate(z, temp, zout)
   end function bilinear_interpolate

   subroutine loess_predict(x, y, xout, yout, status, span)
      real(dp), intent(in) :: x(:), y(:), xout(:)
      real(dp), allocatable, intent(out) :: yout(:)
      integer, intent(out) :: status
      real(dp), intent(in), optional :: span
      real(dp), allocatable :: distance(:), sorted(:), w(:), xc(:)
      real(dp) :: fraction, scale, coef(3)
      integer :: i, n, neighbours, fit_status

      status = ycevo_err_input
      n = size(x)
      allocate(yout(size(xout)))
      yout = huge(1.0_dp)
      if (n < 3 .or. size(y) /= n) return
      fraction = 0.75_dp
      if (present(span)) fraction = span
      if (fraction <= 0.0_dp .or. fraction > 1.0_dp) return
      neighbours = max(3, min(n, ceiling(fraction*real(n,dp))))
      allocate(distance(n), sorted(n), w(n), xc(n))

      do i = 1, size(xout)
         distance = abs(x - xout(i))
         sorted = distance
         call sort_real(sorted)
         scale = sorted(neighbours)
         if (scale <= epsilon(1.0_dp)) then
            yout(i) = sum(y, mask=distance <= epsilon(1.0_dp))/ &
               real(count(distance <= epsilon(1.0_dp)),dp)
            cycle
         end if
         xc = x - xout(i)
         w = 0.0_dp
         where (distance < scale)
            w = (1.0_dp - (distance/scale)**3)**3
         end where
         call weighted_quadratic_fit(xc, y, w, coef, fit_status)
         if (fit_status == ycevo_success) then
            yout(i) = coef(1)
         else
            yout(i) = linear_interpolate(x, y, xout(i))
         end if
      end do
      status = ycevo_success
   end subroutine loess_predict

   subroutine predict_discount(surface, xout, tauout, discount, status, loess_span)
      type(yield_surface_t), intent(in) :: surface
      real(dp), intent(in) :: xout, tauout
      real(dp), intent(out) :: discount
      integer, intent(out) :: status
      real(dp), intent(in), optional :: loess_span
      real(dp), allocatable :: logd(:), smooth(:), at_time(:)
      integer :: j, nx

      status = ycevo_err_input
      discount = huge(1.0_dp)
      nx = size(surface%xgrid)
      if (.not. allocated(surface%tau) .or. .not. allocated(surface%discount)) return
      if (size(surface%discount,1) /= size(surface%tau) .or. &
          size(surface%discount,2) /= nx .or. any(surface%discount <= 0.0_dp)) return
      allocate(at_time(nx), logd(size(surface%tau)))
      do j = 1, nx
         logd = log(surface%discount(:,j))
         call loess_predict(surface%tau, logd, [tauout], smooth, status, loess_span)
         if (status /= ycevo_success) return
         at_time(j) = smooth(1)
      end do
      discount = exp(linear_interpolate(surface%xgrid, at_time, xout))
      status = ycevo_success
   end subroutine predict_discount

   subroutine predict_yield(surface, xout, tauout, yield_value, status, loess_span)
      type(yield_surface_t), intent(in) :: surface
      real(dp), intent(in) :: xout, tauout
      real(dp), intent(out) :: yield_value
      integer, intent(out) :: status
      real(dp), intent(in), optional :: loess_span
      real(dp) :: discount

      call predict_discount(surface, xout, tauout, discount, status, loess_span)
      if (status == ycevo_success) yield_value = discount_to_yield(discount, tauout)
   end subroutine predict_yield

   subroutine sort_real(x)
      real(dp), intent(inout) :: x(:)
      real(dp) :: key
      integer :: i, j

      do i = 2, size(x)
         key = x(i)
         j = i - 1
         do while (j >= 1)
            if (x(j) <= key) exit
            x(j+1) = x(j)
            j = j - 1
         end do
         x(j+1) = key
      end do
   end subroutine sort_real

end module ycevo_prediction
