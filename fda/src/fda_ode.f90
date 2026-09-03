! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from the computational code of R package fda 6.3.0.
module fda_ode
   use r_kinds, only : dp
   use fda_fd, only : fd_type, eval_fd
   implicit none
   private

   type, public :: ode_solution_type
      real(dp), allocatable :: t(:)
      real(dp), allocatable :: y(:, :, :)
   end type ode_solution_type

   public :: linear_ode_rhs
   public :: odesolv

contains

   subroutine linear_ode_rhs(tnow, y, weights, dydt, info)
      real(dp), intent(in) :: tnow !! Current independent-variable value at which the coefficient functions are evaluated.
      real(dp), intent(in) :: y(:, :) !! Current state matrix with derivative order in rows and independent solutions in columns.
      type(fd_type), intent(in) :: weights(:) !! Coefficients `w_1,...,w_m` of the order-`m` homogeneous operator.
      real(dp), allocatable, intent(out) :: dydt(:, :) !! Allocated first-order-system derivative matrix matching the shape of `y`.
      integer, intent(out) :: info !! Zero on success; nonzero for inconsistent shape or coefficient evaluation.
      real(dp), allocatable :: value(:, :), wmat(:, :)
      real(dp) :: point(1)
      integer :: j, m

      m = size(weights)
      info = 0
      if (m < 1 .or. size(y, 1) /= m) then
         allocate(dydt(0, 0))
         info = 1
         return
      end if
      allocate(wmat(m, m))
      wmat = 0.0_dp
      do j = 1, m - 1
         wmat(j, j + 1) = 1.0_dp
      end do
      point(1) = tnow
      do j = 1, m
         call eval_fd(point, weights(j), 0, value, info)
         if (info /= 0 .or. size(value, 2) /= 1) then
            if (info == 0) info = 2
            allocate(dydt(0, 0))
            return
         end if
         wmat(m, j) = -value(1, 1)
      end do
      allocate(dydt(size(y, 1), size(y, 2)))
      dydt = matmul(wmat, y)
   end subroutine linear_ode_rhs

   subroutine odesolv(weights, ystart, solution, info, h0, hmin, hmax, eps, maxstp)
      type(fd_type), intent(in) :: weights(:) !! Coefficient functions of an order-`m` homogeneous linear differential operator.
      real(dp), intent(in) :: ystart(:, :) !! Initial derivatives `0:m-1` by independent solution.
      type(ode_solution_type), intent(out) :: solution !! Adaptive time grid and Cash-Karp state history.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid shape/range, step failure, or step limit.
      real(dp), intent(in), optional :: h0 !! Initial positive trial step; defaults to one hundredth of the integration width.
      real(dp), intent(in), optional :: hmin !! Minimum accepted positive step; defaults to `width*1e-10`.
      real(dp), intent(in), optional :: hmax !! Maximum positive step; defaults to half the integration width.
      real(dp), intent(in), optional :: eps !! Positive relative local-error tolerance; defaults to `1e-4`.
      integer, intent(in), optional :: maxstp !! Maximum adaptive Runge-Kutta steps; defaults to 1000.
      real(dp), allocatable :: dydt(:, :), tstore(:), y(:, :), ystore(:, :, :)
      real(dp) :: eps_use, h, h_max, h_min, h_next, h_try, tbeg, tend, tnow, width
      integer :: m, nsol, nstep, nstored, step_limit

      info = 0
      m = size(weights)
      if (m < 1 .or. size(ystart, 1) /= m .or. size(ystart, 2) < 1) then
         info = 1
         return
      end if
      tbeg = weights(1)%basis%rangeval(1)
      tend = weights(1)%basis%rangeval(2)
      width = tend - tbeg
      if (width <= 0.0_dp) then
         info = 2
         return
      end if
      do nstep = 2, m
         if (maxval(abs(weights(nstep)%basis%rangeval - weights(1)%basis%rangeval)) > &
             64.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(width))) then
            info = 3
            return
         end if
      end do
      eps_use = 1.0e-4_dp
      if (present(eps)) eps_use = eps
      h_min = width * 1.0e-10_dp
      if (present(hmin)) h_min = hmin
      h_max = 0.5_dp * width
      if (present(hmax)) h_max = hmax
      h = width / 100.0_dp
      if (present(h0)) h = h0
      step_limit = 1000
      if (present(maxstp)) step_limit = maxstp
      if (eps_use <= 0.0_dp .or. h_min <= 0.0_dp .or. h_max <= 0.0_dp .or. h <= 0.0_dp .or. step_limit < 1) then
         info = 4
         return
      end if
      h = min(h, h_max)
      nsol = size(ystart, 2)
      allocate(y(m, nsol), tstore(step_limit + 1), ystore(m, nsol, step_limit + 1))
      y = ystart
      tnow = tbeg
      nstored = 1
      tstore(1) = tnow
      ystore(:, :, 1) = y

      do nstep = 1, step_limit
         call linear_ode_rhs(tnow, y, weights, dydt, info)
         if (info /= 0) return
         h_try = min(h, tend - tnow)
         call rkqs(y, dydt, tnow, h_try, weights, eps_use, h_min, h_max, h_next, info)
         if (info /= 0) return
         h = h_next
         nstored = nstored + 1
         tstore(nstored) = tnow
         ystore(:, :, nstored) = y
         if (tnow >= tend - 16.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(tend))) then
            allocate(solution%t(nstored), solution%y(m, nsol, nstored))
            solution%t = tstore(1:nstored)
            solution%y = ystore(:, :, 1:nstored)
            info = 0
            return
         end if
      end do
      allocate(solution%t(nstored), solution%y(m, nsol, nstored))
      solution%t = tstore(1:nstored)
      solution%y = ystore(:, :, 1:nstored)
      info = 5
   end subroutine odesolv

   subroutine rkqs(y, dydt, tnow, htry, weights, eps, hmin, hmax, hnext, info)
      real(dp), intent(inout) :: y(:, :) !! State matrix, replaced by the accepted state at the end of the step.
      real(dp), intent(in) :: dydt(:, :) !! State derivative evaluated at the beginning of the proposed step.
      real(dp), intent(inout) :: tnow !! Independent-variable value, advanced by the accepted step.
      real(dp), intent(in) :: htry !! Positive trial step length.
      type(fd_type), intent(in) :: weights(:) !! Differential-operator coefficient functions passed to intermediate RHS evaluations.
      real(dp), intent(in) :: eps !! Positive relative local-error tolerance.
      real(dp), intent(in) :: hmin !! Positive minimum permitted step length.
      real(dp), intent(in) :: hmax !! Positive maximum proposed next step length.
      real(dp), intent(out) :: hnext !! Proposed step length for the subsequent integration step.
      integer, intent(out) :: info !! Zero on success; nonzero for RHS failure or inability to meet tolerance above `hmin`.
      real(dp), allocatable :: yerr(:, :), yscal(:, :), ytemp(:, :)
      real(dp) :: errmax, h, htemp
      integer :: attempt

      info = 0
      h = htry
      allocate(yscal(size(y, 1), size(y, 2)))
      yscal = abs(y) + abs(h * dydt) + tiny(1.0_dp)
      do attempt = 1, 100
         call rkck(y, dydt, tnow, h, weights, ytemp, yerr, info)
         if (info /= 0) return
         errmax = maxval(abs(yerr / yscal)) / eps
         if (errmax <= 1.0_dp) exit
         htemp = 0.9_dp * h * errmax**(-0.25_dp)
         h = max(0.1_dp * h, htemp)
         if (h < hmin) then
            info = 1
            return
         end if
      end do
      if (errmax > 1.89e-4_dp) then
         hnext = 0.9_dp * h * errmax**(-0.2_dp)
      else
         hnext = 5.0_dp * h
      end if
      hnext = min(hmax, max(hmin, hnext))
      tnow = tnow + h
      y = ytemp
   end subroutine rkqs

   subroutine rkck(y, dydt, tnow, h, weights, yout, yerr, info)
      real(dp), intent(in) :: y(:, :) !! State matrix at the beginning of the Cash-Karp step.
      real(dp), intent(in) :: dydt(:, :) !! State derivative at the beginning of the Cash-Karp step.
      real(dp), intent(in) :: tnow !! Independent-variable value at the beginning of the step.
      real(dp), intent(in) :: h !! Positive step length for the Cash-Karp evaluation.
      type(fd_type), intent(in) :: weights(:) !! Differential-operator coefficient functions for intermediate RHS evaluations.
      real(dp), allocatable, intent(out) :: yout(:, :) !! Fifth-order Cash-Karp state estimate at `tnow+h`.
      real(dp), allocatable, intent(out) :: yerr(:, :) !! Difference between embedded fifth- and fourth-order estimates.
      integer, intent(out) :: info !! Zero on success; nonzero when an intermediate RHS evaluation fails.
      real(dp), allocatable :: ak2(:, :), ak3(:, :), ak4(:, :), ak5(:, :), ak6(:, :), ytemp(:, :)
      real(dp), parameter :: c1 = 37.0_dp / 378.0_dp
      real(dp), parameter :: c3 = 250.0_dp / 621.0_dp
      real(dp), parameter :: c4 = 125.0_dp / 594.0_dp
      real(dp), parameter :: c6 = 512.0_dp / 1771.0_dp
      real(dp), parameter :: dc1 = c1 - 2825.0_dp / 27648.0_dp
      real(dp), parameter :: dc3 = c3 - 18575.0_dp / 48384.0_dp
      real(dp), parameter :: dc4 = c4 - 13525.0_dp / 55296.0_dp
      real(dp), parameter :: dc5 = -277.0_dp / 14336.0_dp
      real(dp), parameter :: dc6 = c6 - 0.25_dp

      allocate(ytemp(size(y, 1), size(y, 2)))
      ytemp = y + h * 0.2_dp * dydt
      call linear_ode_rhs(tnow + 0.2_dp * h, ytemp, weights, ak2, info)
      if (info /= 0) return
      ytemp = y + h * (0.075_dp * dydt + 0.225_dp * ak2)
      call linear_ode_rhs(tnow + 0.3_dp * h, ytemp, weights, ak3, info)
      if (info /= 0) return
      ytemp = y + h * (0.3_dp * dydt - 0.9_dp * ak2 + 1.2_dp * ak3)
      call linear_ode_rhs(tnow + 0.6_dp * h, ytemp, weights, ak4, info)
      if (info /= 0) return
      ytemp = y + h * (-11.0_dp * dydt / 54.0_dp + 2.5_dp * ak2 + (-70.0_dp * ak3 + 35.0_dp * ak4) / 27.0_dp)
      call linear_ode_rhs(tnow + h, ytemp, weights, ak5, info)
      if (info /= 0) return
      ytemp = y + h * (1631.0_dp * dydt / 55296.0_dp + 175.0_dp * ak2 / 512.0_dp + &
         575.0_dp * ak3 / 13824.0_dp + 44275.0_dp * ak4 / 110592.0_dp + 253.0_dp * ak5 / 4096.0_dp)
      call linear_ode_rhs(tnow + 0.875_dp * h, ytemp, weights, ak6, info)
      if (info /= 0) return
      allocate(yout(size(y, 1), size(y, 2)), yerr(size(y, 1), size(y, 2)))
      yout = y + h * (c1 * dydt + c3 * ak3 + c4 * ak4 + c6 * ak6)
      yerr = h * (dc1 * dydt + dc3 * ak3 + dc4 * ak4 + dc5 * ak5 + dc6 * ak6)
   end subroutine rkck

end module fda_ode
