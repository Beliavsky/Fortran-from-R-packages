module ycevo_estimation
   use ycevo_kinds, only : dp
   use ycevo_status, only : ycevo_success, ycevo_err_input, ycevo_err_singular
   use ycevo_kernel, only : epanechnikov
   use ycevo_curve, only : discount_to_yield
   use ycevo_types, only : bond_panel_t, yield_curve_t, yield_surface_t
   use ycevo_linalg, only : solve_linear_system
   implicit none
   private

   public :: calc_dbar, calc_hhat_numerator, interpolation_weights
   public :: estimate_yield, estimate_yield_surface, count_maturing_bonds
   public :: adjust_tau_grid

contains

   subroutine calc_dbar(panel, xgrid, hx, tau, ht, numerator, denominator, &
                        interest, rgrid, hr, status)
      type(bond_panel_t), intent(in) :: panel
      real(dp), intent(in) :: xgrid, hx, tau(:), ht(:)
      real(dp), allocatable, intent(out) :: numerator(:), denominator(:)
      real(dp), intent(in), optional :: interest(:), rgrid, hr
      integer, intent(out) :: status
      real(dp) :: wd, wt, maturity
      integer :: i, j, n

      status = ycevo_err_input
      n = size(tau)
      allocate(numerator(n), denominator(n))
      numerator = 0.0_dp
      denominator = 0.0_dp
      if (size(ht) /= n .or. hx <= 0.0_dp .or. any(ht <= 0.0_dp)) return
      if ((present(interest) .neqv. present(rgrid)) .or. &
          (present(interest) .neqv. present(hr))) return
      if (present(interest)) then
         if (size(interest) /= panel%nday .or. hr <= 0.0_dp) return
      end if

      do i = 1, panel%size()
         wd = epanechnikov((xgrid - real(panel%day(i), dp)/real(panel%nday, dp))/hx)
         if (present(interest)) then
            wd = wd*epanechnikov((rgrid - interest(panel%day(i)))/hr)
         end if
         if (wd <= 0.0_dp) cycle
         maturity = real(panel%tupq(i), dp)/365.0_dp
         do j = 1, n
            wt = epanechnikov((tau(j) - maturity)/ht(j))
            if (wt <= 0.0_dp) cycle
            numerator(j) = numerator(j) + panel%price(i)*panel%cashflow(i)*wt*wd
            denominator(j) = denominator(j) + panel%cashflow(i)**2*wt*wd
         end do
      end do
      status = ycevo_success
   end subroutine calc_dbar

   subroutine calc_hhat_numerator(panel_in, xgrid, hx, tau, ht, tau_p, htp, hhat, &
                                  interest, rgrid, hr, status)
      type(bond_panel_t), intent(in) :: panel_in
      real(dp), intent(in) :: xgrid, hx, tau(:), ht(:), tau_p(:), htp(:)
      real(dp), allocatable, intent(out) :: hhat(:, :)
      real(dp), intent(in), optional :: interest(:), rgrid, hr
      integer, intent(out) :: status
      type(bond_panel_t) :: panel
      real(dp), allocatable :: sum_p(:), row_wp(:)
      real(dp) :: wd, maturity, wt
      integer :: g1, g2, i, j, p, nt, np

      status = ycevo_err_input
      nt = size(tau)
      np = size(tau_p)
      allocate(hhat(nt, np))
      hhat = 0.0_dp
      if (size(ht) /= nt .or. size(htp) /= np .or. hx <= 0.0_dp .or. &
          any(ht <= 0.0_dp) .or. any(htp <= 0.0_dp)) return
      if ((present(interest) .neqv. present(rgrid)) .or. &
          (present(interest) .neqv. present(hr))) return
      if (present(interest)) then
         if (size(interest) /= panel_in%nday .or. hr <= 0.0_dp) return
      end if

      panel = panel_in
      call panel%sort()
      allocate(sum_p(np), row_wp(np))
      g1 = 1
      do while (g1 <= panel%size())
         g2 = g1
         do while (g2 < panel%size())
            if (panel%day(g2+1) /= panel%day(g1) .or. panel%id(g2+1) /= panel%id(g1)) exit
            g2 = g2 + 1
         end do

         wd = epanechnikov((xgrid - real(panel%day(g1), dp)/real(panel%nday, dp))/hx)
         if (present(interest)) then
            wd = wd*epanechnikov((rgrid - interest(panel%day(g1)))/hr)
         end if
         if (wd > 0.0_dp) then
            sum_p = 0.0_dp
            do i = g1, g2
               maturity = real(panel%tupq(i), dp)/365.0_dp
               do p = 1, np
                  sum_p(p) = sum_p(p) + panel%cashflow(i)* &
                     epanechnikov((tau_p(p) - maturity)/htp(p))
               end do
            end do
            do i = g1, g2
               maturity = real(panel%tupq(i), dp)/365.0_dp
               do p = 1, np
                  row_wp(p) = epanechnikov((tau_p(p) - maturity)/htp(p))
               end do
               do j = 1, nt
                  wt = epanechnikov((tau(j) - maturity)/ht(j))
                  if (wt <= 0.0_dp) cycle
                  do p = 1, np
                     hhat(j,p) = hhat(j,p) + &
                        (sum_p(p) - panel%cashflow(i)*row_wp(p))* &
                        panel%cashflow(i)*wt*wd
                  end do
               end do
            end do
         end if
         g1 = g2 + 1
      end do
      status = ycevo_success
   end subroutine calc_hhat_numerator

   subroutine interpolation_weights(tau, tau_p, weights, status)
      real(dp), intent(in) :: tau(:), tau_p(:)
      real(dp), allocatable, intent(out) :: weights(:, :)
      integer, intent(out) :: status
      integer :: j, k, lower, upper, nt, np
      real(dp) :: distance, tol

      nt = size(tau)
      np = size(tau_p)
      allocate(weights(nt, np))
      weights = 0.0_dp
      status = ycevo_err_input
      if (nt == 0 .or. np == 0 .or. any(tau(2:) <= tau(:nt-1))) return
      if (minval(tau_p) < tau(1) .or. maxval(tau_p) > tau(nt)) return
      tol = 50.0_dp*epsilon(1.0_dp)

      do j = 1, np
         k = minloc(abs(tau - tau_p(j)), dim=1)
         if (abs(tau(k) - tau_p(j)) <= tol*max(1.0_dp, abs(tau_p(j)))) then
            weights(k,j) = 1.0_dp
         else
            lower = maxloc(tau, dim=1, mask=tau < tau_p(j))
            upper = minloc(tau, dim=1, mask=tau > tau_p(j))
            distance = tau(upper) - tau(lower)
            ! These are the weights used by the R implementation.
            weights(lower,j) = (tau_p(j) - tau(lower))/distance
            weights(upper,j) = (tau(upper) - tau_p(j))/distance
         end if
      end do
      status = ycevo_success
   end subroutine interpolation_weights

   subroutine estimate_yield(panel, xgrid, hx, tau, ht, curve, status, message, &
                             tau_p, htp, interest, rgrid, hr)
      type(bond_panel_t), intent(in) :: panel
      real(dp), intent(in) :: xgrid, hx, tau(:), ht(:)
      type(yield_curve_t), intent(out) :: curve
      integer, intent(out) :: status
      character(len=*), intent(out), optional :: message
      real(dp), intent(in), optional :: tau_p(:), htp(:), interest(:), rgrid, hr
      real(dp), allocatable :: tp(:), hp(:), num(:), den(:), dbar(:)
      real(dp), allocatable :: hnum(:, :), h(:, :), iw(:, :), matrix(:, :), solution(:)
      integer :: i, n, np, stat_panel
      character(len=256) :: panel_message

      status = ycevo_err_input
      if (present(message)) message = ''
      call panel%validate(stat_panel, panel_message)
      if (stat_panel /= ycevo_success) then
         if (present(message)) message = trim(panel_message)
         return
      end if
      n = size(tau)
      if (n == 0 .or. size(ht) /= n .or. any(tau <= 0.0_dp) .or. &
          any(tau(2:) <= tau(:n-1))) then
         if (present(message)) message = 'tau must be positive and strictly increasing.'
         return
      end if
      if (present(tau_p) .neqv. present(htp)) then
         if (present(message)) message = 'tau_p and htp must be supplied together.'
         return
      end if
      if (present(tau_p)) then
         allocate(tp(size(tau_p)), hp(size(htp)))
         tp = tau_p
         hp = htp
      else
         allocate(tp(n), hp(n))
         tp = tau
         hp = ht
      end if
      np = size(tp)
      if (size(hp) /= np .or. any(tp(2:) <= tp(:np-1)) .or. &
          minval(tp) < tau(1) .or. maxval(tp) > tau(n)) then
         if (present(message)) message = 'Invalid tau_p or htp grid.'
         return
      end if

      call calc_dbar(panel, xgrid, hx, tau, ht, num, den, interest, rgrid, hr, status)
      if (status /= ycevo_success) then
         if (present(message)) message = 'Invalid kernel or covariate arguments.'
         return
      end if
      if (any(den <= tiny(1.0_dp))) then
         if (present(message)) message = 'At least one tau window has no usable cash flows.'
         status = ycevo_err_input
         return
      end if
      dbar = num/den

      call calc_hhat_numerator(panel, xgrid, hx, tau, ht, tp, hp, hnum, &
                                interest, rgrid, hr, status)
      if (status /= ycevo_success) then
         if (present(message)) message = 'Unable to calculate cash-flow cross-products.'
         return
      end if
      allocate(h(n, np))
      do i = 1, n
         h(i,:) = hnum(i,:)/den(i)
      end do
      call interpolation_weights(tau, tp, iw, status)
      if (status /= ycevo_success) then
         if (present(message)) message = 'Unable to form interpolation weights.'
         return
      end if
      allocate(matrix(n,n))
      matrix = matmul(h, transpose(iw))
      do i = 1, n
         matrix(i,i) = matrix(i,i) + 1.0_dp
      end do
      call solve_linear_system(matrix, dbar, solution, status)
      if (status /= ycevo_success) then
         if (present(message)) message = 'The discount-function linear system is singular.'
         status = ycevo_err_singular
         return
      end if

      curve%xgrid = xgrid
      allocate(curve%tau(n), curve%discount(n), curve%yield(n))
      curve%tau = tau
      curve%discount = solution
      curve%yield = discount_to_yield(solution, tau)
      status = ycevo_success
   end subroutine estimate_yield

   subroutine estimate_yield_surface(panel, xgrid, hx, tau, ht, surface, status, message)
      type(bond_panel_t), intent(in) :: panel
      real(dp), intent(in) :: xgrid(:), hx(:), tau(:), ht(:)
      type(yield_surface_t), intent(out) :: surface
      integer, intent(out) :: status
      character(len=*), intent(out), optional :: message
      type(yield_curve_t) :: curve
      integer :: j, nx, nt
      character(len=256) :: local_message

      status = ycevo_err_input
      if (present(message)) message = ''
      nx = size(xgrid)
      nt = size(tau)
      if (size(hx) /= nx) then
         if (present(message)) message = 'xgrid and hx must have equal length.'
         return
      end if
      allocate(surface%xgrid(nx), surface%tau(nt), surface%discount(nt,nx), surface%yield(nt,nx))
      surface%xgrid = xgrid
      surface%tau = tau
      do j = 1, nx
         call estimate_yield(panel, xgrid(j), hx(j), tau, ht, curve, status, local_message)
         if (status /= ycevo_success) then
            if (present(message)) write(message, '(a,i0,a,a)') 'Curve ', j, ': ', trim(local_message)
            return
         end if
         surface%discount(:,j) = curve%discount
         surface%yield(:,j) = curve%yield
      end do
      status = ycevo_success
   end subroutine estimate_yield_surface

   subroutine count_maturing_bonds(panel, xgrid, hx, tau, ht, counts, min_cashflow)
      type(bond_panel_t), intent(in) :: panel
      real(dp), intent(in) :: xgrid, hx, tau(:), ht(:)
      integer, allocatable, intent(out) :: counts(:)
      real(dp), intent(in), optional :: min_cashflow
      real(dp) :: cutoff, wd, wt, maturity
      integer :: i, j

      cutoff = 100.0_dp
      if (present(min_cashflow)) cutoff = min_cashflow
      allocate(counts(size(tau)))
      counts = 0
      do i = 1, panel%size()
         if (panel%cashflow(i) < cutoff) cycle
         wd = epanechnikov((xgrid - real(panel%day(i),dp)/real(panel%nday,dp))/hx)
         if (wd <= 0.0_dp) cycle
         maturity = real(panel%tupq(i),dp)/365.0_dp
         do j = 1, size(tau)
            wt = epanechnikov((tau(j) - maturity)/ht(j))
            if (wt > 0.01_dp) counts(j) = counts(j) + 1
         end do
      end do
   end subroutine count_maturing_bonds


   subroutine adjust_tau_grid(panel, xgrid, hx, tau, ht, tau_out, ht_out, min_points)
      type(bond_panel_t), intent(in) :: panel
      real(dp), intent(in) :: xgrid, hx, tau(:), ht(:)
      real(dp), allocatable, intent(out) :: tau_out(:), ht_out(:)
      integer, intent(in), optional :: min_points
      integer, allocatable :: counts(:)
      logical, allocatable :: keep(:)
      real(dp), allocatable :: adjusted(:)
      integer :: minimum, n, last_good, i, left, right, nout
      real(dp) :: gap

      minimum = 5
      if (present(min_points)) minimum = min_points
      n = size(tau)
      if (n == 0 .or. size(ht) /= n) then
         allocate(tau_out(0), ht_out(0))
         return
      end if
      call count_maturing_bonds(panel, xgrid, hx, tau, ht, counts)
      last_good = n
      do while (last_good > 0 .and. counts(last_good) < minimum)
         last_good = last_good - 1
      end do
      if (last_good == 0) then
         allocate(tau_out(0), ht_out(0))
         return
      end if
      allocate(keep(last_good), adjusted(last_good))
      keep = counts(:last_good) >= minimum
      adjusted = ht(:last_good)

      i = 1
      do while (i <= last_good)
         if (keep(i)) then
            i = i + 1
            cycle
         end if
         left = i - 1
         do while (i <= last_good .and. .not. keep(i))
            i = i + 1
         end do
         right = i
         if (left >= 1 .and. right <= last_good) then
            gap = 0.5_dp*(tau(right) - tau(left))
            adjusted(left) = max(adjusted(left), gap)
            adjusted(right) = max(adjusted(right), gap)
         end if
      end do

      nout = count(keep)
      allocate(tau_out(nout), ht_out(nout))
      tau_out = pack(tau(:last_good), keep)
      ht_out = pack(adjusted, keep)
   end subroutine adjust_tau_grid

end module ycevo_estimation
