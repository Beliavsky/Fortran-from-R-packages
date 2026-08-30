! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of gbm3 3.0.3; see NOTICE.md and upstream/.
module gbm3_distributions
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use gbm3_kinds, only : dp
   use gbm3_constants
   use gbm3_math, only : sigmoid, softplus, log_add_exp, exp_diff, weighted_quantile, location_m_tdist, quiet_nan
   use gbm3_types, only : gbm_options, gbm_tree
   implicit none
   private
   public :: dist_init_f, dist_working_response, dist_deviance
   public :: dist_fit_best_constant, dist_bag_improvement

contains

   real(dp) function dist_init_f(y, offset, weight, options) result(init_f)
      real(dp), intent(in) :: y(:), offset(:), weight(:)
      type(gbm_options), intent(in) :: options
      real(dp), allocatable :: v(:)
      real(dp) :: sw, step, num, den, p, log_num, log_den, log_w
      real(dp) :: max_off, min_off
      integer :: i, it

      if (size(y) /= size(offset) .or. size(y) /= size(weight)) error stop "dist_init_f: shape mismatch"
      sw = sum(weight)
      select case (options%distribution)
      case (GBM_GAUSSIAN)
         init_f = sum(weight * (y - offset)) / sw
      case (GBM_BERNOULLI)
         init_f = 0.0_dp
         step = 1.0_dp
         do it = 1, 6
            if (abs(step) <= 0.001_dp) exit
            num = 0.0_dp
            den = 0.0_dp
            do i = 1, size(y)
               p = sigmoid(offset(i) + init_f)
               num = num + weight(i) * (y(i) - p)
               den = den + weight(i) * p * (1.0_dp - p)
            end do
            if (abs(den) <= tiny(1.0_dp)) exit
            step = num / den
            init_f = init_f + step
         end do
      case (GBM_POISSON)
         log_num = -huge(1.0_dp)
         log_den = -huge(1.0_dp)
         max_off = -huge(1.0_dp)
         min_off = huge(1.0_dp)
         do i = 1, size(y)
            if (weight(i) > 0.0_dp) then
               log_w = log(weight(i))
               log_den = log_add_exp(log_den, log_w + offset(i))
               if (y(i) > 0.0_dp) log_num = log_add_exp(log_num, log_w + log(y(i)))
               max_off = max(max_off, offset(i))
               min_off = min(min_off, offset(i))
            end if
         end do
         if (log_num <= -huge(1.0_dp)) then
            init_f = -19.0_dp
         else
            init_f = log_num - log_den
            if (max_off + init_f > 19.0_dp) init_f = 19.0_dp - max_off
            if (min_off + init_f < -19.0_dp) init_f = -19.0_dp - min_off
         end if
      case (GBM_GAMMA)
         log_num = -huge(1.0_dp)
         log_den = -huge(1.0_dp)
         max_off = -huge(1.0_dp)
         min_off = huge(1.0_dp)
         do i = 1, size(y)
            if (weight(i) > 0.0_dp) then
               log_w = log(weight(i))
               log_den = log_add_exp(log_den, log_w)
               if (y(i) > 0.0_dp) log_num = log_add_exp(log_num, log_w + log(y(i)) - offset(i))
               max_off = max(max_off, offset(i))
               min_off = min(min_off, offset(i))
            end if
         end do
         if (log_num <= -huge(1.0_dp)) then
            init_f = -19.0_dp
         else
            init_f = log_num - log_den
         end if
         if (max_off + init_f > 19.0_dp) init_f = 19.0_dp - max_off
         if (min_off + init_f < -19.0_dp) init_f = -19.0_dp - min_off
      case (GBM_LAPLACE)
         allocate(v(size(y)))
         v = y - offset
         init_f = weighted_quantile(v, weight, 0.5_dp)
      case (GBM_TDIST)
         allocate(v(size(y)))
         v = y - offset
         init_f = location_m_tdist(v, weight, options%t_df, 0.5_dp)
      case (GBM_QUANTILE)
         allocate(v(size(y)))
         v = y - offset
         init_f = weighted_quantile(v, weight, options%quantile_alpha)
      case (GBM_ADABOOST)
         log_num = -huge(1.0_dp)
         log_den = -huge(1.0_dp)
         do i = 1, size(y)
            if (weight(i) <= 0.0_dp) cycle
            log_w = log(weight(i))
            if (abs(y(i) - 1.0_dp) <= tiny(1.0_dp)) then
               log_num = log_add_exp(log_num, log_w - offset(i))
            else
               log_den = log_add_exp(log_den, log_w + offset(i))
            end if
         end do
         init_f = 0.5_dp * (log_num - log_den)
      case (GBM_HUBERIZED)
         num = sum(weight, mask=(abs(y - 1.0_dp) <= tiny(1.0_dp)))
         den = sum(weight, mask=(abs(y - 1.0_dp) > tiny(1.0_dp)))
         if (den <= tiny(1.0_dp)) then
            init_f = huge(1.0_dp)
         else
            init_f = num / den
         end if
      case (GBM_TWEEDIE)
         log_num = -huge(1.0_dp)
         log_den = -huge(1.0_dp)
         max_off = -huge(1.0_dp)
         min_off = huge(1.0_dp)
         do i = 1, size(y)
            if (weight(i) <= 0.0_dp) cycle
            log_w = log(weight(i))
            log_den = log_add_exp(log_den, log_w + offset(i) * (2.0_dp - options%tweedie_power))
            if (y(i) > 0.0_dp) then
               log_num = log_add_exp(log_num, log_w + log(y(i)) + offset(i) * (1.0_dp - options%tweedie_power))
            end if
            max_off = max(max_off, offset(i))
            min_off = min(min_off, offset(i))
         end do
         if (log_num <= -huge(1.0_dp)) then
            init_f = -19.0_dp
         else
            init_f = log_num - log_den
         end if
         if (max_off + init_f > 19.0_dp) init_f = 19.0_dp - max_off
         if (min_off + init_f < -19.0_dp) init_f = -19.0_dp - min_off
      case default
         error stop "dist_init_f: distribution handled by another module or unknown"
      end select
   end function dist_init_f

   subroutine dist_working_response(y, offset, f, options, residual)
      real(dp), intent(in) :: y(:), offset(:), f(:)
      type(gbm_options), intent(in) :: options
      real(dp), intent(out) :: residual(size(y))
      integer :: i
      real(dp) :: eta, u, outcome

      select case (options%distribution)
      case (GBM_GAUSSIAN)
         residual = y - offset - f
      case (GBM_BERNOULLI)
         do i = 1, size(y)
            eta = f(i) + offset(i)
            residual(i) = y(i) - sigmoid(eta)
         end do
      case (GBM_POISSON)
         do i = 1, size(y)
            eta = f(i) + offset(i)
            residual(i) = y(i) - exp(eta)
         end do
      case (GBM_GAMMA)
         do i = 1, size(y)
            eta = f(i) + offset(i)
            residual(i) = y(i) * exp(-eta) - 1.0_dp
         end do
      case (GBM_LAPLACE)
         do i = 1, size(y)
            if (y(i) - offset(i) - f(i) > 0.0_dp) then
               residual(i) = 1.0_dp
            else
               residual(i) = -1.0_dp
            end if
         end do
      case (GBM_TDIST)
         do i = 1, size(y)
            u = y(i) - offset(i) - f(i)
            residual(i) = 2.0_dp * u / (options%t_df + u * u)
         end do
      case (GBM_QUANTILE)
         do i = 1, size(y)
            if (y(i) > f(i) + offset(i)) then
               residual(i) = options%quantile_alpha
            else
               residual(i) = -(1.0_dp - options%quantile_alpha)
            end if
         end do
      case (GBM_ADABOOST)
         do i = 1, size(y)
            outcome = 2.0_dp * y(i) - 1.0_dp
            residual(i) = -outcome * exp(-outcome * (offset(i) + f(i)))
         end do
      case (GBM_HUBERIZED)
         do i = 1, size(y)
            outcome = 2.0_dp * y(i) - 1.0_dp
            eta = f(i) + offset(i)
            if (outcome * eta < -1.0_dp) then
               residual(i) = -4.0_dp * outcome
            else if (1.0_dp - outcome * eta < 0.0_dp) then
               residual(i) = 0.0_dp
            else
               residual(i) = -2.0_dp * outcome * (1.0_dp - outcome * eta)
            end if
         end do
      case (GBM_TWEEDIE)
         do i = 1, size(y)
            eta = f(i) + offset(i)
            residual(i) = y(i) * exp(eta * (1.0_dp - options%tweedie_power)) - &
                          exp(eta * (2.0_dp - options%tweedie_power))
         end do
      case default
         error stop "dist_working_response: distribution handled by another module or unknown"
      end select
   end subroutine dist_working_response

   real(dp) function dist_deviance(y, offset, weight, f, options) result(dev)
      real(dp), intent(in) :: y(:), offset(:), weight(:), f(:)
      type(gbm_options), intent(in) :: options
      real(dp) :: loss, sw, eta, u, outcome, pwr
      integer :: i

      loss = 0.0_dp
      sw = sum(weight)
      select case (options%distribution)
      case (GBM_GAUSSIAN)
         loss = sum(weight * (y - offset - f) ** 2)
         dev = loss_over_weight(loss, sw, 1.0_dp)
      case (GBM_BERNOULLI)
         do i = 1, size(y)
            eta = f(i) + offset(i)
            loss = loss + weight(i) * (y(i) * eta - softplus(eta))
         end do
         dev = loss_over_weight(-loss, sw, 2.0_dp)
      case (GBM_POISSON)
         do i = 1, size(y)
            eta = f(i) + offset(i)
            loss = loss + weight(i) * (y(i) * eta - exp(eta))
         end do
         dev = loss_over_weight(-loss, sw, 2.0_dp)
      case (GBM_GAMMA)
         do i = 1, size(y)
            eta = f(i) + offset(i)
            loss = loss + weight(i) * (y(i) * exp(-eta) + eta)
         end do
         dev = loss_over_weight(loss, sw, 2.0_dp)
      case (GBM_LAPLACE)
         loss = sum(weight * abs(y - offset - f))
         dev = loss_over_weight(loss, sw, 1.0_dp)
      case (GBM_TDIST)
         do i = 1, size(y)
            u = y(i) - offset(i) - f(i)
            loss = loss + weight(i) * log(options%t_df + u * u)
         end do
         dev = loss_over_weight(loss, sw, 1.0_dp)
      case (GBM_QUANTILE)
         do i = 1, size(y)
            eta = f(i) + offset(i)
            if (y(i) > eta) then
               loss = loss + weight(i) * options%quantile_alpha * (y(i) - eta)
            else
               loss = loss + weight(i) * (1.0_dp - options%quantile_alpha) * (eta - y(i))
            end if
         end do
         dev = loss_over_weight(loss, sw, 1.0_dp)
      case (GBM_ADABOOST)
         do i = 1, size(y)
            outcome = 2.0_dp * y(i) - 1.0_dp
            loss = loss + weight(i) * exp(-outcome * (offset(i) + f(i)))
         end do
         dev = loss_over_weight(loss, sw, 1.0_dp)
      case (GBM_HUBERIZED)
         do i = 1, size(y)
            outcome = 2.0_dp * y(i) - 1.0_dp
            eta = offset(i) + f(i)
            if (outcome * eta < -1.0_dp) then
               loss = loss - weight(i) * 4.0_dp * outcome * eta
            else if (1.0_dp - outcome * eta < 0.0_dp) then
               loss = loss
            else
               loss = loss + weight(i) * (1.0_dp - outcome * eta) ** 2
            end if
         end do
         dev = loss_over_weight(loss, sw, 1.0_dp)
      case (GBM_TWEEDIE)
         pwr = options%tweedie_power
         do i = 1, size(y)
            eta = f(i) + offset(i)
            loss = loss + weight(i) * (y(i) ** (2.0_dp - pwr) / ((1.0_dp - pwr) * (2.0_dp - pwr)) - &
                   y(i) * exp(eta * (1.0_dp - pwr)) / (1.0_dp - pwr) + &
                   exp(eta * (2.0_dp - pwr)) / (2.0_dp - pwr))
         end do
         dev = loss_over_weight(loss, sw, 2.0_dp)
      case default
         error stop "dist_deviance: distribution handled by another module or unknown"
      end select
   end function dist_deviance

   pure real(dp) function loss_over_weight(loss, sw, factor) result(v)
      real(dp), intent(in) :: loss, sw, factor
      if (abs(sw) <= tiny(1.0_dp)) then
         if (abs(loss) <= tiny(1.0_dp)) then
            v = ieee_value(0.0_dp, ieee_quiet_nan)
         else
            v = sign(huge(1.0_dp), loss)
         end if
      else
         v = factor * loss / sw
      end if
   end function loss_over_weight

   subroutine dist_fit_best_constant(y, offset, weight, f, residual, in_bag, assignment, &
                                     min_obs, options, tree)
      real(dp), intent(in) :: y(:), offset(:), weight(:), f(:), residual(:)
      logical, intent(in) :: in_bag(:)
      integer, intent(in) :: assignment(:), min_obs
      type(gbm_options), intent(in) :: options
      type(gbm_tree), intent(inout) :: tree

      real(dp), allocatable :: num(:), den(:), maxeta(:), mineta(:), scale(:)
      real(dp), allocatable :: v(:), wv(:)
      logical, allocatable :: mask(:)
      real(dp) :: eta, temp, outcome, lw, scaled, rescale
      integer :: i, node

      if (options%distribution == GBM_GAUSSIAN) return
      allocate(num(tree%n_nodes), den(tree%n_nodes))
      num = 0.0_dp
      den = 0.0_dp

      select case (options%distribution)
      case (GBM_BERNOULLI)
         do i = 1, size(y)
            if (.not. in_bag(i)) cycle
            node = assignment(i)
            num(node) = num(node) + weight(i) * residual(i)
            den(node) = den(node) + weight(i) * (y(i) - residual(i)) * (1.0_dp - y(i) + residual(i))
         end do
         do node = 1, tree%n_nodes
            if (.not. tree%nodes(node)%is_terminal) cycle
            if (abs(den(node)) <= tiny(1.0_dp)) then
               tree%nodes(node)%prediction = 0.0_dp
            else
               temp = num(node) / den(node)
               tree%nodes(node)%prediction = min(10.0_dp, max(-10.0_dp, temp))
            end if
         end do

      case (GBM_POISSON, GBM_GAMMA, GBM_TWEEDIE)
         allocate(maxeta(tree%n_nodes), mineta(tree%n_nodes))
         maxeta = -huge(1.0_dp)
         mineta = huge(1.0_dp)
         do i = 1, size(y)
            if (.not. in_bag(i)) cycle
            node = assignment(i)
            eta = f(i) + offset(i)
            select case (options%distribution)
            case (GBM_POISSON)
               num(node) = num(node) + weight(i) * y(i)
               den(node) = den(node) + weight(i) * exp(eta)
            case (GBM_GAMMA)
               num(node) = num(node) + weight(i) * y(i) * exp(-eta)
               den(node) = den(node) + weight(i)
            case (GBM_TWEEDIE)
               num(node) = num(node) + weight(i) * y(i) * exp(eta * (1.0_dp - options%tweedie_power))
               den(node) = den(node) + weight(i) * exp(eta * (2.0_dp - options%tweedie_power))
            end select
            maxeta(node) = max(maxeta(node), eta)
            mineta(node) = min(mineta(node), eta)
         end do
         do node = 1, tree%n_nodes
            if (.not. tree%nodes(node)%is_terminal) cycle
            if (abs(num(node)) <= tiny(1.0_dp)) then
               tree%nodes(node)%prediction = -19.0_dp
            else if (abs(den(node)) <= tiny(1.0_dp)) then
               tree%nodes(node)%prediction = 0.0_dp
            else
               tree%nodes(node)%prediction = log(num(node) / den(node))
            end if
            if (maxeta(node) > -huge(1.0_dp)) then
               tree%nodes(node)%prediction = min(tree%nodes(node)%prediction, 19.0_dp - maxeta(node))
               tree%nodes(node)%prediction = max(tree%nodes(node)%prediction, -19.0_dp - mineta(node))
            end if
         end do

      case (GBM_LAPLACE, GBM_TDIST, GBM_QUANTILE)
         allocate(mask(size(y)))
         do node = 1, tree%n_nodes
            if (.not. tree%nodes(node)%is_terminal) cycle
            if (tree%nodes(node)%num_obs < min_obs) cycle
            mask = in_bag .and. (assignment == node)
            if (count(mask) == 0) cycle
            v = pack(y - offset - f, mask)
            wv = pack(weight, mask)
            select case (options%distribution)
            case (GBM_LAPLACE)
               tree%nodes(node)%prediction = weighted_quantile(v, wv, 0.5_dp)
            case (GBM_TDIST)
               tree%nodes(node)%prediction = location_m_tdist(v, wv, options%t_df, 0.5_dp)
            case (GBM_QUANTILE)
               tree%nodes(node)%prediction = weighted_quantile(v, wv, options%quantile_alpha)
            end select
         end do

      case (GBM_ADABOOST)
         allocate(scale(tree%n_nodes))
         scale = -huge(1.0_dp)
         do i = 1, size(y)
            if (.not. in_bag(i) .or. weight(i) <= 0.0_dp) cycle
            node = assignment(i)
            outcome = 2.0_dp * y(i) - 1.0_dp
            eta = f(i) + offset(i)
            lw = log(weight(i)) - outcome * eta
            scaled = 1.0_dp
            if (lw > scale(node)) then
               if (scale(node) > -huge(1.0_dp)) then
                  rescale = exp(scale(node) - lw)
                  num(node) = num(node) * rescale
                  den(node) = den(node) * rescale
               end if
               scale(node) = lw
            else
               scaled = exp(lw - scale(node))
            end if
            num(node) = num(node) + outcome * scaled
            den(node) = den(node) + scaled
         end do
         do node = 1, tree%n_nodes
            if (.not. tree%nodes(node)%is_terminal) cycle
            if (abs(den(node)) <= tiny(1.0_dp)) then
               tree%nodes(node)%prediction = 0.0_dp
            else
               tree%nodes(node)%prediction = num(node) / den(node)
            end if
         end do

      case (GBM_HUBERIZED)
         do i = 1, size(y)
            if (.not. in_bag(i)) cycle
            node = assignment(i)
            outcome = 2.0_dp * y(i) - 1.0_dp
            eta = f(i) + offset(i)
            if (outcome * eta < -1.0_dp) then
               num(node) = num(node) + weight(i) * 4.0_dp * outcome
               den(node) = den(node) - weight(i) * 4.0_dp * outcome * eta
            else if (1.0_dp - outcome * eta < 0.0_dp) then
               continue
            else
               num(node) = num(node) + weight(i) * 2.0_dp * outcome * (1.0_dp - outcome * eta)
               den(node) = den(node) + weight(i) * (1.0_dp - outcome * eta) ** 2
            end if
         end do
         do node = 1, tree%n_nodes
            if (.not. tree%nodes(node)%is_terminal) cycle
            if (abs(den(node)) <= tiny(1.0_dp)) then
               tree%nodes(node)%prediction = 0.0_dp
            else
               tree%nodes(node)%prediction = num(node) / den(node)
            end if
         end do
      case default
         error stop "dist_fit_best_constant: distribution handled by another module or unknown"
      end select
   end subroutine dist_fit_best_constant

   real(dp) function dist_bag_improvement(y, offset, weight, f, delta, in_bag, options) result(improvement)
      real(dp), intent(in) :: y(:), offset(:), weight(:), f(:), delta(:)
      logical, intent(in) :: in_bag(:)
      type(gbm_options), intent(in) :: options
      real(dp) :: val, sw, eta, outcome, u, v, pwr, old_loss, new_loss
      integer :: i

      val = 0.0_dp
      sw = 0.0_dp
      do i = 1, size(y)
         if (in_bag(i)) cycle
         eta = f(i) + offset(i)
         select case (options%distribution)
         case (GBM_GAUSSIAN)
            val = val + weight(i) * options%shrinkage * delta(i) * &
                  (2.0_dp * (y(i) - eta) - options%shrinkage * delta(i))
         case (GBM_BERNOULLI)
            if (abs(y(i) - 1.0_dp) <= tiny(1.0_dp)) val = val + weight(i) * options%shrinkage * delta(i)
            val = val + weight(i) * (softplus(eta) - softplus(eta + options%shrinkage * delta(i)))
         case (GBM_POISSON)
            val = val + weight(i) * (y(i) * options%shrinkage * delta(i) - &
                  exp(eta + options%shrinkage * delta(i)) + exp(eta))
         case (GBM_GAMMA)
            val = val + weight(i) * (y(i) * exp(-eta) * (1.0_dp - exp(-options%shrinkage * delta(i))) - &
                  options%shrinkage * delta(i))
         case (GBM_LAPLACE)
            val = val + weight(i) * (abs(y(i) - eta) - abs(y(i) - eta - options%shrinkage * delta(i)))
         case (GBM_TDIST)
            u = y(i) - eta
            v = y(i) - eta - options%shrinkage * delta(i)
            val = val + weight(i) * (log(options%t_df + u * u) - log(options%t_df + v * v))
         case (GBM_QUANTILE)
            old_loss = quantile_loss(y(i), eta, options%quantile_alpha)
            new_loss = quantile_loss(y(i), eta + options%shrinkage * delta(i), options%quantile_alpha)
            val = val + weight(i) * (old_loss - new_loss)
         case (GBM_ADABOOST)
            outcome = 2.0_dp * y(i) - 1.0_dp
            val = val + weight(i) * exp_diff(-outcome * eta, &
                  -outcome * (eta + options%shrinkage * delta(i)))
         case (GBM_HUBERIZED)
            old_loss = huberized_loss(y(i), eta)
            new_loss = huberized_loss(y(i), eta + options%shrinkage * delta(i))
            val = val + weight(i) * (old_loss - new_loss)
         case (GBM_TWEEDIE)
            pwr = options%tweedie_power
            val = val + weight(i) * (exp(eta * (1.0_dp - pwr)) * y(i) / (1.0_dp - pwr) * &
                  (exp(options%shrinkage * delta(i) * (1.0_dp - pwr)) - 1.0_dp) + &
                  exp(eta * (2.0_dp - pwr)) / (2.0_dp - pwr) * &
                  (1.0_dp - exp(options%shrinkage * delta(i) * (2.0_dp - pwr))))
         case default
            error stop "dist_bag_improvement: distribution handled by another module or unknown"
         end select
         ! Upstream huberized omits weight accumulation in one branch; preserve
         ! that behavior explicitly.
         if (options%distribution == GBM_HUBERIZED) then
            outcome = 2.0_dp * y(i) - 1.0_dp
            if (outcome * eta < -1.0_dp .or. 1.0_dp - outcome * eta < 0.0_dp) sw = sw + weight(i)
         else
            sw = sw + weight(i)
         end if
      end do

      if (abs(sw) <= tiny(1.0_dp)) then
         if (abs(val) <= tiny(1.0_dp)) then
            improvement = quiet_nan()
         else
            improvement = sign(huge(1.0_dp), val)
         end if
      else
         select case (options%distribution)
         case (GBM_GAMMA, GBM_TWEEDIE)
            improvement = 2.0_dp * val / sw
         case default
            improvement = val / sw
         end select
      end if
   end function dist_bag_improvement

   pure real(dp) function quantile_loss(y, eta, alpha) result(v)
      real(dp), intent(in) :: y, eta, alpha
      if (y > eta) then
         v = alpha * (y - eta)
      else
         v = (1.0_dp - alpha) * (eta - y)
      end if
   end function quantile_loss

   pure real(dp) function huberized_loss(y, eta) result(v)
      real(dp), intent(in) :: y, eta
      real(dp) :: outcome
      outcome = 2.0_dp * y - 1.0_dp
      if (outcome * eta < -1.0_dp) then
         v = -4.0_dp * outcome * eta
      else if (1.0_dp - outcome * eta < 0.0_dp) then
         v = 0.0_dp
      else
         v = (1.0_dp - outcome * eta) ** 2
      end if
   end function huberized_loss

end module gbm3_distributions
