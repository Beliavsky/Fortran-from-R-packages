module ycevo_simulation
   use ycevo_kinds, only : dp
   use ycevo_kernel, only : epanechnikov_quantile
   use ycevo_curve, only : get_yield_at
   use ycevo_types, only : bond_panel_t
   implicit none
   private

   public :: simulate_bond_panel, seed_random_number

contains

   subroutine seed_random_number(seed)
      integer, intent(in) :: seed
      integer, allocatable :: values(:)
      integer :: i, n

      call random_seed(size=n)
      allocate(values(n))
      do i = 1, n
         values(i) = modulo(seed + 104729*i + 8191*i*i, huge(1) - 1)
         if (values(i) == 0) values(i) = i
      end do
      call random_seed(put=values)
   end subroutine seed_random_number

   subroutine simulate_bond_panel(panel, nday, n_bonds, seed, max_maturity_years)
      type(bond_panel_t), intent(out) :: panel
      integer, intent(in), optional :: nday, n_bonds, seed
      real(dp), intent(in), optional :: max_maturity_years
      integer :: nd, nb, i, b, d, count_rows, pos, payment_abs
      integer, allocatable :: maturity(:)
      real(dp), allocatable :: coupon(:), u(:), prices(:), cashflows(:), maturities(:)
      real(dp) :: max_years, time, y, discount

      nd = 60
      nb = 40
      max_years = 10.0_dp
      if (present(nday)) nd = nday
      if (present(n_bonds)) nb = n_bonds
      if (present(max_maturity_years)) max_years = max_maturity_years
      if (present(seed)) call seed_random_number(seed)

      allocate(maturity(nb), coupon(nb), u(nb))
      call random_number(u)
      do b = 1, nb
         maturity(b) = nd + nint((0.5_dp + u(b)*(max_years - 0.5_dp))*365.0_dp)
      end do
      call random_number(u)
      coupon = epanechnikov_quantile(u, 4.0_dp, sqrt(7.0_dp))
      coupon = max(coupon, 0.0_dp)

      count_rows = 0
      do d = 1, nd
         do b = 1, nb
            payment_abs = maturity(b)
            do while (payment_abs > d)
               count_rows = count_rows + 1
               payment_abs = payment_abs - 182
            end do
         end do
      end do

      allocate(panel%day(count_rows), panel%id(count_rows), panel%tupq(count_rows), &
               panel%price(count_rows), panel%cashflow(count_rows))
      panel%nday = nd
      pos = 0
      do d = 1, nd
         time = merge(real(d-1,dp)/real(nd-1,dp), 0.0_dp, nd > 1)
         do b = 1, nb
            i = 0
            payment_abs = maturity(b)
            do while (payment_abs > d)
               i = i + 1
               payment_abs = payment_abs - 182
            end do
            allocate(prices(i), cashflows(i), maturities(i))
            payment_abs = maturity(b)
            do i = 1, size(prices)
               maturities(i) = real(payment_abs - d,dp)/365.0_dp
               cashflows(i) = coupon(b)/2.0_dp
               if (i == 1) cashflows(i) = cashflows(i) + 100.0_dp
               y = get_yield_at(time, maturities(i))
               discount = exp(-maturities(i)*y)
               prices(i) = cashflows(i)*discount
               payment_abs = payment_abs - 182
            end do
            do i = 1, size(prices)
               pos = pos + 1
               panel%day(pos) = d
               panel%id(pos) = b
               panel%tupq(pos) = nint(maturities(i)*365.0_dp)
               panel%cashflow(pos) = cashflows(i)
               panel%price(pos) = sum(prices)
            end do
            deallocate(prices, cashflows, maturities)
         end do
      end do
      call panel%sort()
   end subroutine simulate_bond_panel

end module ycevo_simulation
