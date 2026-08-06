module ycevo_curve
   use ycevo_kinds, only : dp
   implicit none
   private

   public :: nelson_siegel, get_yield_at, generate_yield
   public :: discount_to_yield, default_tau_grid, bandwidth_from_tau

contains

   elemental real(dp) function nelson_siegel(maturity, b0, b1, b2, t1, t2) result(yield_value)
      real(dp), intent(in) :: maturity
      real(dp), intent(in), optional :: b0, b1, b2, t1, t2
      real(dp) :: bb0, bb1, bb2, tt1, tt2, z1, z2

      bb0 = 0.0_dp
      bb1 = 0.05_dp
      bb2 = 2.0_dp
      tt1 = 0.75_dp
      tt2 = 125.0_dp
      if (present(b0)) bb0 = b0
      if (present(b1)) bb1 = b1
      if (present(b2)) bb2 = b2
      if (present(t1)) tt1 = t1
      if (present(t2)) tt2 = t2

      if (maturity <= epsilon(1.0_dp)) then
         yield_value = bb0 + bb1
         return
      end if
      z1 = maturity / tt1
      z2 = maturity / tt2
      yield_value = bb0 + bb1 * (one_minus_exp_neg(z1) / z1) + &
         bb2 * (one_minus_exp_neg(z2) / z2 - exp(-z2))
   end function nelson_siegel

   elemental real(dp) function get_yield_at(time, maturity, b0, b1, b2, t1, t2, &
                                             linear, quadratic, cubic) result(yield_value)
      real(dp), intent(in) :: time, maturity
      real(dp), intent(in), optional :: b0, b1, b2, t1, t2
      real(dp), intent(in), optional :: linear, quadratic, cubic
      real(dp) :: ll, qq, cc

      ll = -0.55_dp
      qq = 0.55_dp
      cc = -0.55_dp
      if (present(linear)) ll = linear
      if (present(quadratic)) qq = quadratic
      if (present(cubic)) cc = cubic

      yield_value = nelson_siegel(maturity, b0, b1, b2, t1, t2) * &
         (1.0_dp + ll*time + qq*time*time + cc*time*time*time)
   end function get_yield_at

   subroutine generate_yield(yields, n_qdate, periods, b0, b1, b2, t1, t2, &
                             linear, quadratic, cubic)
      real(dp), allocatable, intent(out) :: yields(:, :)
      integer, intent(in), optional :: n_qdate, periods
      real(dp), intent(in), optional :: b0, b1, b2, t1, t2
      real(dp), intent(in), optional :: linear, quadratic, cubic
      integer :: nq, np, i, j
      real(dp) :: time, maturity

      nq = 12
      np = 36
      if (present(n_qdate)) nq = n_qdate
      if (present(periods)) np = periods
      allocate(yields(np, nq))
      do j = 1, nq
         time = real(j, dp) / real(nq, dp)
         do i = 1, np
            maturity = real(i, dp) / (real(np, dp) / 10.0_dp)
            yields(i, j) = get_yield_at(time, maturity, b0, b1, b2, t1, t2, &
                                        linear, quadratic, cubic)
         end do
      end do
   end subroutine generate_yield

   elemental real(dp) function discount_to_yield(discount, tau) result(yield_value)
      real(dp), intent(in) :: discount, tau

      if (discount <= 0.0_dp .or. tau <= 0.0_dp) then
         yield_value = huge(1.0_dp)
      else
         yield_value = -log(discount) / tau
      end if
   end function discount_to_yield

   subroutine default_tau_grid(max_tau, tau)
      real(dp), intent(in) :: max_tau
      real(dp), allocatable, intent(out) :: tau(:)
      real(dp), allocatable :: work(:)
      integer :: n, day

      allocate(work(200))
      n = 0
      do day = 30, 180, 30
         call append_day(day)
      end do
      do day = 240, 730, 60
         call append_day(day)
      end do
      do day = 810, 2190, 90
         call append_day(day)
      end do
      do day = 2280, 7300, 120
         call append_day(day)
      end do
      do day = 7482, int(30.6_dp*365.0_dp), 182
         call append_day(day)
      end do
      allocate(tau(n))
      tau = work(:n)

   contains

      subroutine append_day(day_value)
         integer, intent(in) :: day_value
         real(dp) :: value

         value = real(day_value, dp) / 365.0_dp
         if (value < max_tau) then
            n = n + 1
            work(n) = value
         end if
      end subroutine append_day
   end subroutine default_tau_grid

   subroutine bandwidth_from_tau(tau, bandwidth)
      real(dp), intent(in) :: tau(:)
      real(dp), allocatable, intent(out) :: bandwidth(:)
      integer :: i, n

      n = size(tau)
      allocate(bandwidth(n))
      if (n == 1) then
         bandwidth(1) = max(tau(1), 1.0_dp/365.0_dp)
         return
      end if
      bandwidth(1) = tau(2) - tau(1)
      bandwidth(n) = tau(n) - tau(n-1)
      do i = 2, n - 1
         bandwidth(i) = max(tau(i) - tau(i-1), tau(i+1) - tau(i))
      end do
   end subroutine bandwidth_from_tau

   elemental real(dp) function one_minus_exp_neg(x) result(value)
      real(dp), intent(in) :: x

      if (abs(x) < 1.0e-5_dp) then
         value = x*(1.0_dp - x*(0.5_dp - x*(1.0_dp/6.0_dp - &
            x*(1.0_dp/24.0_dp - x/120.0_dp))))
      else
         value = 1.0_dp - exp(-x)
      end if
   end function one_minus_exp_neg

end module ycevo_curve
