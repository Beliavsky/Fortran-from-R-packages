! SPDX-License-Identifier: GPL-3.0-only
!
! Experimental modern Fortran translation of computational methods from the
! R package rmgarch, copyright (C) 2008-2025 Alexios Galanos.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 3 only.
module rmgarch_scenario
   use rmgarch_kinds, only : dp
   use rmgarch_math, only : outer_product, normalize_covariance, make_positive_definite
   use rmgarch_distributions, only : random_multivariate
   use rmgarch_types, only : dcc_fit_result
   implicit none
   private

   public :: simulate_dcc_scenarios, scenario_moments, portfolio_scenarios

contains

   subroutine simulate_dcc_scenarios(fit, horizon, nsim, scenarios, correlations, qpaths, valid)
      type(dcc_fit_result), intent(in) :: fit
      integer, intent(in) :: horizon, nsim
      real(dp), intent(out) :: scenarios(horizon,size(fit%qbar,1),nsim)
      real(dp), intent(out) :: correlations(size(fit%qbar,1),size(fit%qbar,2),horizon,nsim)
      real(dp), intent(out), optional :: qpaths(size(fit%qbar,1),size(fit%qbar,2),horizon,nsim)
      logical, intent(out), optional :: valid
      real(dp), allocatable :: zhistory(:,:), qhistory(:,:,:)
      real(dp) :: qnext(size(fit%qbar,1),size(fit%qbar,2)), rnext(size(fit%qbar,1),size(fit%qbar,2))
      real(dp) :: nvec(size(fit%qbar,1)), zero(size(fit%qbar,1)), draw(size(fit%qbar,1))
      real(dp) :: intercept(size(fit%qbar,1),size(fit%qbar,2))
      integer :: simulation, h, i, source, nz, nq, m
      logical :: ok, draw_ok

      scenarios = 0.0_dp
      correlations = 0.0_dp
      if (present(qpaths)) qpaths = 0.0_dp
      ok = horizon > 0 .and. nsim > 0 .and. allocated(fit%qbar) .and. &
         allocated(fit%nbar) .and. allocated(fit%q) .and. &
         allocated(fit%standardized_residuals) .and. allocated(fit%spec%alpha) .and. &
         allocated(fit%spec%beta) .and. allocated(fit%spec%gamma)
      if (.not. ok) then
         if (present(valid)) valid = .false.
         return
      end if
      m = size(fit%qbar,1)
      nz = size(fit%q,3)
      nq = size(fit%q,3)
      zero = 0.0_dp
      intercept = (1.0_dp-sum(fit%spec%alpha)-sum(fit%spec%beta))*fit%qbar- &
         sum(fit%spec%gamma)*fit%nbar

      do simulation = 1, nsim
         allocate(zhistory(nz+horizon,m),qhistory(m,m,nq+horizon))
         zhistory = 0.0_dp
         qhistory = 0.0_dp
         zhistory(1:nz,:) = fit%standardized_residuals
         qhistory(:,:,1:nq) = fit%q
         do h = 1, horizon
            qnext = intercept
            do i = 1, fit%spec%p
               source = nz+h-i
               if (source > nz) then
                  qnext = qnext+fit%spec%alpha(i)*outer_product(zhistory(source,:),zhistory(source,:))
               else
                  qnext = qnext+fit%spec%alpha(i)*fit%qbar
               end if
            end do
            do i = 1, fit%spec%g
               source = nz+h-i
               if (source > nz) then
                  where (zhistory(source,:) < 0.0_dp)
                     nvec = zhistory(source,:)
                  elsewhere
                     nvec = 0.0_dp
                  end where
                  qnext = qnext+fit%spec%gamma(i)*outer_product(nvec,nvec)
               else
                  qnext = qnext+fit%spec%gamma(i)*fit%nbar
               end if
            end do
            do i = 1, fit%spec%q
               source = nq+h-i
               if (source >= 1) then
                  qnext = qnext+fit%spec%beta(i)*qhistory(:,:,source)
               else
                  qnext = qnext+fit%spec%beta(i)*fit%qbar
               end if
            end do
            qnext = make_positive_definite(qnext,1.0e-10_dp)
            rnext = normalize_covariance(qnext)
            call random_multivariate(fit%spec%distribution,zero,rnext,draw, &
               fit%spec%shape,draw_ok)
            if (draw_ok) zhistory(nz+h,:) = draw
            ok = ok .and. draw_ok
            qhistory(:,:,nq+h) = qnext
            scenarios(h,:,simulation) = zhistory(nz+h,:)
            correlations(:,:,h,simulation) = rnext
            if (present(qpaths)) qpaths(:,:,h,simulation) = qnext
         end do
         deallocate(zhistory,qhistory)
      end do
      if (present(valid)) valid = ok
   end subroutine simulate_dcc_scenarios

   subroutine scenario_moments(scenarios, mean, covariance, coskew, cokurt)
      real(dp), intent(in) :: scenarios(:,:,:)
      real(dp), intent(out) :: mean(size(scenarios,1),size(scenarios,2))
      real(dp), intent(out) :: covariance(size(scenarios,2),size(scenarios,2),size(scenarios,1))
      real(dp), intent(out), optional :: coskew(size(scenarios,2),size(scenarios,2), &
         size(scenarios,2),size(scenarios,1))
      real(dp), intent(out), optional :: cokurt(size(scenarios,2),size(scenarios,2), &
         size(scenarios,2),size(scenarios,2),size(scenarios,1))
      real(dp) :: centered(size(scenarios,2))
      integer :: h, s, i, j, k, l, nsim, m

      nsim = size(scenarios,3)
      m = size(scenarios,2)
      mean = 0.0_dp
      covariance = 0.0_dp
      if (present(coskew)) coskew = 0.0_dp
      if (present(cokurt)) cokurt = 0.0_dp
      if (nsim < 1) return
      do h = 1, size(scenarios,1)
         mean(h,:) = sum(scenarios(h,:,:),dim=2)/real(nsim,dp)
         do s = 1, nsim
            centered = scenarios(h,:,s)-mean(h,:)
            do j = 1, m
               do i = 1, m
                  covariance(i,j,h) = covariance(i,j,h)+centered(i)*centered(j)
               end do
            end do
            if (present(coskew)) then
               do k = 1, m
                  do j = 1, m
                     do i = 1, m
                        coskew(i,j,k,h) = coskew(i,j,k,h)+centered(i)*centered(j)*centered(k)
                     end do
                  end do
               end do
            end if
            if (present(cokurt)) then
               do l = 1, m
                  do k = 1, m
                     do j = 1, m
                        do i = 1, m
                           cokurt(i,j,k,l,h) = cokurt(i,j,k,l,h)+ &
                              centered(i)*centered(j)*centered(k)*centered(l)
                        end do
                     end do
                  end do
               end do
            end if
         end do
         if (nsim > 1) covariance(:,:,h) = covariance(:,:,h)/real(nsim-1,dp)
         if (present(coskew)) coskew(:,:,:,h) = coskew(:,:,:,h)/real(nsim,dp)
         if (present(cokurt)) cokurt(:,:,:,:,h) = cokurt(:,:,:,:,h)/real(nsim,dp)
      end do
   end subroutine scenario_moments

   subroutine portfolio_scenarios(scenarios, weights, portfolio)
      real(dp), intent(in) :: scenarios(:,:,:), weights(:,:)
      real(dp), intent(out) :: portfolio(size(scenarios,1),size(scenarios,3))
      integer :: h, s

      portfolio = 0.0_dp
      if (size(weights,2) /= size(scenarios,2)) return
      if (size(weights,1) /= 1 .and. size(weights,1) /= size(scenarios,1)) return
      do h = 1, size(scenarios,1)
         do s = 1, size(scenarios,3)
            if (size(weights,1) == 1) then
               portfolio(h,s) = dot_product(weights(1,:),scenarios(h,:,s))
            else
               portfolio(h,s) = dot_product(weights(h,:),scenarios(h,:,s))
            end if
         end do
      end do
   end subroutine portfolio_scenarios

end module rmgarch_scenario
