! SPDX-License-Identifier: GPL-2.0-or-later
module flexmix_em_utils
   use flexmix_kinds, only : dp
   use flexmix_types, only : flexmix_class_weighted, flexmix_class_hard
   use flexmix_numeric, only : normalize_log_probabilities, solve_linear_system
   implicit none
   private
   public :: initialize_posteriors
   public :: classify_posteriors
   public :: group_log_densities
   public :: group_first_mask
   public :: fit_constant_prior
   public :: fit_multinomial_prior
   public :: multinomial_prior
   public :: select_active_components
   public :: weighted_loglikelihood
   public :: argmax_rows

contains

   pure subroutine initialize_posteriors(signal, k, posterior, initial_cluster, initial_posterior, info)
      real(dp), intent(in) :: signal(:) !! Scalar initialization signal, commonly the response or first clustering variable.
      integer, intent(in) :: k !! Requested number of mixture components, at least one.
      real(dp), intent(out) :: posterior(:,:) !! Initialized membership weights, shape `(n, k)`.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based initial component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional nonnegative initial membership matrix, shape `(n.
      integer, intent(out) :: info !! Zero on success; nonzero for inconsistent initialization input.
      integer, allocatable :: order(:)
      real(dp) :: rowsum
      integer :: i, j, n, pos, tmp, target
      n = size(signal)
      posterior = 0.0_dp
      info = 0
      if (size(posterior, 1) /= n .or. size(posterior, 2) /= k .or. k < 1) then
         info = -1
         return
      end if
      if (present(initial_posterior)) then
         if (size(initial_posterior, 1) /= n .or. size(initial_posterior, 2) /= k) then
            info = -2
            return
         end if
         if (any(initial_posterior < 0.0_dp)) then
            info = -3
            return
         end if
         posterior = initial_posterior
         do i = 1, n
            rowsum = sum(posterior(i,:))
            if (rowsum <= 0.0_dp) then
               posterior(i,:) = 1.0_dp / real(k, dp)
            else
               posterior(i,:) = posterior(i,:) / rowsum
            end if
         end do
         return
      end if
      if (present(initial_cluster)) then
         if (size(initial_cluster) /= n) then
            info = -4
            return
         end if
         if (any(initial_cluster < 1) .or. any(initial_cluster > k)) then
            info = -5
            return
         end if
         posterior = 0.1_dp
         do i = 1, n
            posterior(i, initial_cluster(i)) = 0.9_dp
            posterior(i,:) = posterior(i,:) / sum(posterior(i,:))
         end do
         return
      end if
      allocate(order(n))
      order = [(i, i = 1, n)]
      do i = 2, n
         tmp = order(i)
         pos = i - 1
         do while (pos >= 1)
            if (signal(order(pos)) <= signal(tmp)) exit
            order(pos + 1) = order(pos)
            pos = pos - 1
         end do
         order(pos + 1) = tmp
      end do
      posterior = 0.1_dp
      do j = 1, n
         target = min(k, 1 + (j - 1) * k / max(1, n))
         posterior(order(j), target) = 0.9_dp
      end do
      do i = 1, n
         posterior(i,:) = posterior(i,:) / sum(posterior(i,:))
      end do
   end subroutine initialize_posteriors

   pure subroutine classify_posteriors(posterior, classify, classified)
      real(dp), intent(in) :: posterior(:,:) !! Current row-normalized posterior matrix, shape `(n, k)`.
      integer, intent(in) :: classify !! Classification mode; weighted or hard/CEM in this translation.
      real(dp), intent(out) :: classified(:,:) !! Posterior weights used in the M step, shape `(n, k)`.
      integer :: i, winner
      classified = posterior
      if (classify == flexmix_class_hard) then
         classified = 0.0_dp
         do i = 1, size(posterior, 1)
            winner = maxloc(posterior(i,:), dim=1)
            classified(i,winner) = 1.0_dp
         end do
      else if (classify /= flexmix_class_weighted) then
         classified = posterior
      end if
   end subroutine classify_posteriors

   pure subroutine group_first_mask(group, first)
      integer, intent(in) :: group(:) !! Integer group labels, size `n`; equal labels define one grouped observation.
      logical, intent(out) :: first(:) !! True at the first occurrence of each group label, size `n`.
      integer :: i, j
      first = .true.
      do i = 1, size(group)
         do j = 1, i - 1
            if (group(j) == group(i)) then
               first(i) = .false.
               exit
            end if
         end do
      end do
   end subroutine group_first_mask

   pure subroutine group_log_densities(log_density, group, grouped)
      real(dp), intent(in) :: log_density(:,:) !! Per-observation component log densities, shape `(n, k)`.
      integer, intent(in) :: group(:) !! Integer group labels, size `n`; observations in one group share a posterior.
      real(dp), intent(out) :: grouped(:,:) !! Group-summed component log densities repeated on group rows, shape `(n, k)`.
      integer :: i, j
      grouped = 0.0_dp
      do i = 1, size(log_density, 1)
         do j = 1, size(log_density, 1)
            if (group(j) == group(i)) grouped(i,:) = grouped(i,:) + log_density(j,:)
         end do
      end do
   end subroutine group_log_densities

   pure subroutine fit_constant_prior(classified, case_weights, first, active, prior)
      real(dp), intent(in) :: classified(:,:) !! M-step component weights, shape `(n, k0)`.
      real(dp), intent(in) :: case_weights(:) !! Nonnegative integer-like case weights represented as reals, size `n`.
      logical, intent(in) :: first(:) !! Group-first mask; priors are estimated once per grouped observation.
      logical, intent(in) :: active(:) !! Active-component mask, size `k0`.
      real(dp), intent(out) :: prior(:) !! Fitted constant component probabilities, size `k0`.
      real(dp) :: sw, total
      integer :: i, k
      prior = 0.0_dp
      sw = 0.0_dp
      do i = 1, size(classified, 1)
         if (.not. first(i)) cycle
         sw = sw + case_weights(i)
         do k = 1, size(classified, 2)
            if (active(k)) prior(k) = prior(k) + case_weights(i) * classified(i,k)
         end do
      end do
      if (sw <= tiny(1.0_dp)) then
         total = real(count(active), dp)
         if (total > 0.0_dp) then
            do k = 1, size(prior)
               if (active(k)) prior(k) = 1.0_dp / total
            end do
         end if
      else
         prior = prior / sw
         total = sum(prior)
         if (total > 0.0_dp) prior = prior / total
      end if
   end subroutine fit_constant_prior

   pure subroutine multinomial_prior(x, coef, active, prior)
      real(dp), intent(in) :: x(:,:) !! Concomitant-model design matrix, shape `(n, q)`.
      real(dp), intent(in) :: coef(:,:) !! Softmax coefficients, shape `(q, k0)`, with a zero reference column.
      logical, intent(in) :: active(:) !! Active-component mask, size `k0`.
      real(dp), intent(out) :: prior(:,:) !! Row-wise concomitant component probabilities, shape `(n, k0)`.
      real(dp) :: m, s
      integer :: i, k
      prior = 0.0_dp
      do i = 1, size(x, 1)
         m = -huge(1.0_dp)
         do k = 1, size(coef, 2)
            if (active(k)) m = max(m, dot_product(x(i,:), coef(:,k)))
         end do
         s = 0.0_dp
         do k = 1, size(coef, 2)
            if (active(k)) then
               prior(i,k) = exp(dot_product(x(i,:), coef(:,k)) - m)
               s = s + prior(i,k)
            end if
         end do
         if (s > 0.0_dp) then
            prior(i,:) = prior(i,:) / s
         else
            do k = 1, size(coef, 2)
               if (active(k)) prior(i,k) = 1.0_dp / real(count(active), dp)
            end do
         end if
      end do
   end subroutine multinomial_prior

   pure subroutine fit_multinomial_prior(x, classified, case_weights, first, active, coef, prior, info)
      real(dp), intent(in) :: x(:,:) !! Concomitant-model design matrix, shape `(n, q)` and constant within any group.
      real(dp), intent(in) :: classified(:,:) !! M-step component targets, shape `(n, k0)`.
      real(dp), intent(in) :: case_weights(:) !! Nonnegative case weights, size `n`.
      logical, intent(in) :: first(:) !! Group-first mask selecting one row per grouped observation.
      logical, intent(in) :: active(:) !! Active-component mask, size `k0`.
      real(dp), intent(inout) :: coef(:,:) !! Softmax coefficients, shape `(q, k0)`; updated in place.
      real(dp), intent(out) :: prior(:,:) !! Fitted row-wise concomitant priors, shape `(n, k0)`.
      integer, intent(out) :: info !! Zero on convergence; nonzero if a Newton system cannot be solved.
      integer, allocatable :: classes(:), row_index(:)
      real(dp), allocatable :: grad(:), hess(:,:), delta(:), prob(:,:), old_coef(:,:)
      real(dp) :: wi, block_weight, change, ridge
      integer :: a, b, c, d, i, iter, k, m, q, r, s, nrows, solve_info
      q = size(x, 2)
      nrows = count(first)
      allocate(row_index(nrows))
      r = 0
      do i = 1, size(first)
         if (first(i)) then
            r = r + 1
            row_index(r) = i
         end if
      end do
      allocate(classes(count(active)))
      r = 0
      do k = 1, size(active)
         if (active(k)) then
            r = r + 1
            classes(r) = k
         end if
      end do
      m = size(classes)
      if (m == 0) then
         coef = 0.0_dp
         prior = 0.0_dp
         info = 1
         return
      end if
      coef(:, classes(1)) = 0.0_dp
      if (m == 1) then
         coef = 0.0_dp
         call multinomial_prior(x, coef, active, prior)
         info = 0
         return
      end if
      allocate(grad(q * (m - 1)), hess(q * (m - 1), q * (m - 1)))
      allocate(delta(q * (m - 1)), prob(size(x,1), size(classified,2)), old_coef(size(coef,1), size(coef,2)))
      info = 0
      do iter = 1, 40
         old_coef = coef
         call multinomial_prior(x, coef, active, prob)
         grad = 0.0_dp
         hess = 0.0_dp
         do r = 1, nrows
            i = row_index(r)
            wi = case_weights(i)
            if (wi <= 0.0_dp) cycle
            do a = 2, m
               c = classes(a)
               do s = 1, q
                  grad((a - 2) * q + s) = grad((a - 2) * q + s) + &
                     wi * x(i,s) * (classified(i,c) - prob(i,c))
               end do
               do b = 2, m
                  d = classes(b)
                  if (c == d) then
                     block_weight = prob(i,c) * (1.0_dp - prob(i,c))
                  else
                     block_weight = -prob(i,c) * prob(i,d)
                  end if
                  do s = 1, q
                     do k = 1, q
                        hess((a - 2) * q + s, (b - 2) * q + k) = &
                           hess((a - 2) * q + s, (b - 2) * q + k) + wi * block_weight * x(i,s) * x(i,k)
                     end do
                  end do
               end do
            end do
         end do
         ridge = 1.0e-8_dp * max(1.0_dp, maxval(abs(hess)))
         do s = 1, size(hess, 1)
            hess(s,s) = hess(s,s) + ridge
         end do
         call solve_linear_system(hess, grad, delta, solve_info)
         if (solve_info /= 0) then
            info = solve_info
            exit
         end if
         do a = 2, m
            c = classes(a)
            coef(:,c) = coef(:,c) + delta((a - 2) * q + 1:(a - 1) * q)
         end do
         coef(:, classes(1)) = 0.0_dp
         change = maxval(abs(coef - old_coef))
         if (change <= 1.0e-8_dp * (1.0_dp + maxval(abs(coef)))) exit
      end do
      call multinomial_prior(x, coef, active, prior)
   end subroutine fit_multinomial_prior

   pure subroutine select_active_components(prior, minprior, active, changed)
      real(dp), intent(in) :: prior(:) !! Current component prior masses, size `k0`.
      real(dp), intent(in) :: minprior !! Minimum retained marginal prior probability.
      logical, intent(inout) :: active(:) !! Active-component mask updated by removing too-small components.
      logical, intent(out) :: changed !! True if one or more active components were removed.
      integer :: k, best
      changed = .false.
      do k = 1, size(active)
         if (active(k) .and. prior(k) < minprior) then
            active(k) = .false.
            changed = .true.
         end if
      end do
      if (count(active) == 0) then
         best = maxloc(prior, dim=1)
         active(best) = .true.
         changed = .true.
      end if
   end subroutine select_active_components

   pure function weighted_loglikelihood(log_row_sum, case_weights, first) result(value)
      real(dp), intent(in) :: log_row_sum(:) !! Log mixture likelihood for each observation/group row.
      real(dp), intent(in) :: case_weights(:) !! Nonnegative case weights, size `n`.
      logical, intent(in) :: first(:) !! Group-first mask selecting one contribution per grouped observation.
      real(dp) :: value
      integer :: i
      value = 0.0_dp
      do i = 1, size(log_row_sum)
         if (first(i)) value = value + case_weights(i) * log_row_sum(i)
      end do
   end function weighted_loglikelihood

   pure subroutine argmax_rows(x, index)
      real(dp), intent(in) :: x(:,:) !! Row-wise scores or probabilities, shape `(n, k)`.
      integer, intent(out) :: index(:) !! One-based location of each row maximum, size `n`.
      integer :: i
      do i = 1, size(x, 1)
         index(i) = maxloc(x(i,:), dim=1)
      end do
   end subroutine argmax_rows

end module flexmix_em_utils
