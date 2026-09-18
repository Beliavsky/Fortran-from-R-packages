module grpreg_mfdr
   use grpreg_kinds, only : dp
   use grpreg_types, only : grpreg_fit_type, mfdr_result_type
   use grpreg_preprocess, only : grpreg_preprocess_type, prepare_design
   use grpreg_math, only : logistic, chi_square_cdf
   implicit none
   private
   public :: mfdr_grpreg

contains

   subroutine mfdr_grpreg(fit, result, x)
      type(grpreg_fit_type), intent(in) :: fit !! Fitted grouped model whose marginal false-discovery rates are estimated.
      type(mfdr_result_type), intent(out) :: result !! Expected false groups, selected groups, and marginal FDR by lambda.
      real(dp), optional, intent(in) :: x(:,:) !! Raw design matrix required for binomial and Cox curvature calculations.
      type(grpreg_preprocess_type) :: prep
      real(dp), allocatable :: xw(:,:), w(:), hazard(:), risk(:)
      real(dp) :: tau_sq, threshold, tail
      integer :: l, g, qg, j, i
      if (.not. (trim(fit%penalty) == 'grLasso' .or. trim(fit%penalty) == 'grMCP' .or. &
         trim(fit%penalty) == 'grSCAD')) error stop 'mfdr: supported only for grLasso, grMCP, and grSCAD'
      if ((trim(fit%family) == 'binomial' .or. trim(fit%family) == 'cox') .and. .not. present(x)) &
         error stop 'mfdr: X is required for binomial and Cox fits'
      allocate(result%ef(fit%nlambda),result%selected(fit%nlambda),result%mfdr(fit%nlambda))
      allocate(w(fit%n))
      w = 0.0_dp
      result%ef = 0.0_dp
      do l = 1, fit%nlambda
         result%selected(l) = 0
         do g = 1, fit%ngroups
            if (g > size(fit%group_multiplier)) cycle
            if (fit%group_multiplier(g) <= tiny(1.0_dp)) cycle
            if (any(abs(fit%beta(:,l)) > 1.0e-12_dp .and. fit%group == g)) result%selected(l) = result%selected(l)+1
         end do
      end do
      if (trim(fit%family) /= 'gaussian') then
         call prepare_design(x,fit%group,.false.,xw,prep,fit%group_multiplier)
         if (trim(fit%family) == 'cox') allocate(hazard(fit%n),risk(fit%n))
      end if
      do l = 1, fit%nlambda
         if (trim(fit%family) == 'gaussian') then
            tau_sq = fit%deviance(l)/max(real(fit%n,dp)-fit%df(l),1.0_dp)
            do g = 1, fit%ngroups
               qg = count(fit%group == g)
               if (qg == 0 .or. g > size(fit%group_multiplier)) cycle
               if (fit%group_multiplier(g) <= tiny(1.0_dp)) cycle
               threshold = real(fit%n,dp)*(fit%lambda(l)*fit%group_multiplier(g))**2*fit%alpha/max(tau_sq,tiny(1.0_dp))
               tail = max(0.0_dp,1.0_dp-chi_square_cdf(threshold,real(qg,dp)))
               result%ef(l) = result%ef(l)+tail
            end do
         else if (trim(fit%family) == 'binomial') then
            w = logistic(fit%eta(:,l))*(1.0_dp-logistic(fit%eta(:,l)))
            do g = 1, fit%ngroups
               qg = count(prep%group_reduced == g)
               if (qg == 0 .or. g > size(fit%group_multiplier)) cycle
               if (fit%group_multiplier(g) <= tiny(1.0_dp)) cycle
               tau_sq = 0.0_dp
               do j = 1, size(xw,2)
                  if (prep%group_reduced(j) == g) tau_sq = tau_sq+sum(w*xw(:,j)**2)/real(fit%n,dp)
               end do
               threshold = real(fit%n,dp)*(fit%lambda(l)*fit%group_multiplier(g)**2)**2*fit%alpha/ &
                  max(tau_sq,tiny(1.0_dp))
               tail = max(0.0_dp,1.0_dp-chi_square_cdf(threshold,real(qg,dp)))
               result%ef(l) = result%ef(l)+tail
            end do
         else
            hazard = exp(min(fit%eta(:,l),40.0_dp))
            risk(fit%n) = hazard(fit%n)
            do i = fit%n-1, 1, -1
               risk(i) = risk(i+1)+hazard(i)
            end do
            do j = 1, fit%n
               w(j) = 0.0_dp
               do i = 1, j
                  if (fit%fail(i) > 0.5_dp) &
                     w(j) = w(j)+fit%fail(i)*hazard(j)/risk(i)*(1.0_dp-hazard(j)/risk(i))
               end do
            end do
            do g = 1, fit%ngroups
               qg = count(prep%group_reduced == g)
               if (qg == 0 .or. g > size(fit%group_multiplier)) cycle
               if (fit%group_multiplier(g) <= tiny(1.0_dp)) cycle
               tau_sq = 0.0_dp
               do j = 1, size(xw,2)
                  if (prep%group_reduced(j) == g) tau_sq = tau_sq+sum(w*xw(:,j)**2)/real(fit%n,dp)
               end do
               threshold = real(fit%n,dp)*(fit%lambda(l)*fit%group_multiplier(g)**2)**2*fit%alpha/ &
                  max(tau_sq,tiny(1.0_dp))
               tail = max(0.0_dp,1.0_dp-chi_square_cdf(threshold,real(qg,dp)))
               result%ef(l) = result%ef(l)+tail
            end do
         end if
         result%ef(l) = min(result%ef(l),real(result%selected(l),dp))
         if (result%selected(l) > 0) then
            result%mfdr(l) = result%ef(l)/real(result%selected(l),dp)
         else
            result%mfdr(l) = 0.0_dp
         end if
      end do
   end subroutine mfdr_grpreg

end module grpreg_mfdr
