! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2012-2014 Alexios Galanos and Bernhard Pfaff.
! Modern Fortran translation Copyright (C) 2026 OpenAI.
! This derivative work is distributed under GPL-3.0-or-later.
! Derived from parma 1.7, Copyright (C) 2012-2014 Alexios Galanos and Bernhard Pfaff.
module parma_portfolio
   use parma_kinds, only: dp
   use parma_types
   use parma_linalg, only: covariance_matrix, mean_columns, project_box_budget, &
      project_linear_equalities, vector_norm
   use parma_risk, only: risk_value, quadratic_variance, benchmark_variance
   use parma_constraints, only: constraint_penalty, constraints_feasible
   use parma_cmaes, only: cmaes_minimize
   use parma_utility, only: cara2_value, cara4_value
   implicit none
   private
   public :: parmaspec, parmasolve, parmafrontier, parmautility2, parmautility4
   public :: validate_spec, repair_weights, reward_value, spec_risk_value

contains

   subroutine parmaspec(spec, data, risk, objective, mu, covariance, lb, ub, &
      budget, target, alpha, moment, threshold, benchmark, benchmark_covariance, &
      initial, leverage, max_positions, info)
      type(parma_spec), intent(out) :: spec
      real(dp), intent(in), optional :: data(:,:), mu(:), covariance(:,:)
      integer, intent(in), optional :: risk, objective, max_positions
      real(dp), intent(in), optional :: lb(:), ub(:), budget, target, alpha, moment, threshold
      real(dp), intent(in), optional :: benchmark(:), benchmark_covariance(:), initial(:), leverage
      integer, intent(out), optional :: info
      integer :: m, ierr

      ierr = 0
      if (present(data)) then
         m = size(data,2)
         allocate(spec%data(size(data,1),m))
         spec%data = data
      else if (present(mu)) then
         m = size(mu)
      else if (present(covariance)) then
         m = size(covariance,1)
      else
         m = 0
         ierr = -1
      end if
      if (m > 0) then
         allocate(spec%mu(m),spec%cov(m,m),spec%lb(m),spec%ub(m),spec%initial(m))
         if (present(mu)) then
            spec%mu = mu
         else if (present(data)) then
            spec%mu = mean_columns(data)
         else
            spec%mu = 0.0_dp
         end if
         if (present(covariance)) then
            spec%cov = covariance
         else if (present(data)) then
            spec%cov = covariance_matrix(data)
         else
            spec%cov = 0.0_dp
         end if
         spec%lb = 0.0_dp
         spec%ub = 1.0_dp
         if (present(lb)) spec%lb = lb
         if (present(ub)) spec%ub = ub
         spec%initial = 1.0_dp/real(m,dp)
         if (present(initial)) spec%initial = initial
      end if
      if (present(risk)) spec%risk = risk
      if (present(objective)) spec%objective = objective
      if (present(budget)) spec%budget = budget
      if (present(leverage)) spec%leverage = leverage
      if (present(max_positions)) spec%max_positions = max_positions
      if (present(target)) spec%target = target
      if (present(alpha)) spec%alpha = alpha
      if (present(moment)) spec%moment = moment
      if (present(threshold)) spec%threshold = threshold
      if (present(benchmark)) then
         allocate(spec%benchmark(size(benchmark)))
         spec%benchmark = benchmark
      end if
      if (present(benchmark_covariance)) then
         allocate(spec%benchmark_cov(size(benchmark_covariance)))
         spec%benchmark_cov = benchmark_covariance
      end if
      if (present(info)) info = ierr
   end subroutine parmaspec

   subroutine validate_spec(spec, info, message)
      type(parma_spec), intent(in) :: spec
      integer, intent(out) :: info
      character(len=*), intent(out), optional :: message
      integer :: m

      info = 0
      if (.not. allocated(spec%mu)) then
         info = 1
         if (present(message)) message = 'mu is not allocated'
         return
      end if
      m = size(spec%mu)
      if (.not. allocated(spec%lb) .or. .not. allocated(spec%ub)) then
         info = 2
         if (present(message)) message = 'bounds are not allocated'
         return
      end if
      if (size(spec%lb) /= m .or. size(spec%ub) /= m) then
         info = 3
         if (present(message)) message = 'bounds have incorrect size'
         return
      end if
      if (any(spec%ub < spec%lb)) then
         info = 4
         if (present(message)) message = 'upper bound is below lower bound'
         return
      end if
      if (spec%leverage <= 0.0_dp) then
         if (sum(spec%lb) > spec%budget+1.0e-12_dp .or. sum(spec%ub) < spec%budget-1.0e-12_dp) then
            info = 5
            if (present(message)) message = 'budget is infeasible under bounds'
            return
         end if
      else if (spec%leverage > sum(max(abs(spec%lb),abs(spec%ub)))+1.0e-12_dp) then
         info = 5
         if (present(message)) message = 'leverage is infeasible under bounds'
         return
      end if
      if (spec%risk /= risk_ev .and. .not. allocated(spec%data)) then
         info = 6
         if (present(message)) message = 'scenario data are required for selected risk measure'
         return
      end if
      if (spec%risk == risk_ev .and. .not. allocated(spec%cov) .and. .not. allocated(spec%data)) then
         info = 7
         if (present(message)) message = 'covariance or scenario data are required'
         return
      end if
      if (allocated(spec%benchmark)) then
         if (.not. allocated(spec%data) .or. size(spec%benchmark) /= size(spec%data,1)) then
            info = 8
            if (present(message)) message = 'benchmark length must equal number of scenarios'
            return
         end if
      end if
      if (allocated(spec%benchmark_cov)) then
         if (size(spec%benchmark_cov) /= m+1) then
            info = 9
            if (present(message)) message = 'benchmark covariance must have m+1 entries'
            return
         end if
      end if
      if (present(message)) message = 'ok'
   end subroutine validate_spec

   function reward_value(weights, spec) result(value)
      real(dp), intent(in) :: weights(:)
      type(parma_spec), intent(in) :: spec
      real(dp) :: value
      value = dot_product(weights,spec%mu)
   end function reward_value

   subroutine repair_weights(x, spec, weights, info)
      real(dp), intent(in) :: x(:)
      type(parma_spec), intent(in) :: spec
      real(dp), intent(out) :: weights(:)
      integer, intent(out) :: info
      real(dp), allocatable :: temp(:), adjusted(:), direction(:)
      real(dp) :: delta, denom, reward, scale
      integer :: iter, ierr

      allocate(temp(size(x)),adjusted(size(x)),direction(size(x)))
      temp = x
      info = 0
      do iter = 1, 20
         if (spec%leverage > 0.0_dp) then
            adjusted = min(max(temp,spec%lb),spec%ub)
            if (sum(abs(adjusted)) > tiny(1.0_dp)) then
               scale = spec%leverage/sum(abs(adjusted))
               adjusted = min(max(scale*adjusted,spec%lb),spec%ub)
            end if
            ierr = 0
         else
            call project_box_budget(temp,spec%lb,spec%ub,spec%budget,adjusted,ierr)
         end if
         if (ierr /= 0) then
            info = ierr
            weights = adjusted
            return
         end if
         temp = adjusted
         if (allocated(spec%eq_a) .and. allocated(spec%eq_b)) then
            call project_linear_equalities(temp,spec%eq_a,spec%eq_b,adjusted,ierr)
            if (ierr == 0) temp = adjusted
         end if
         if (allocated(spec%mu)) then
            reward = dot_product(temp,spec%mu)
            direction = spec%mu-sum(spec%mu)/real(size(spec%mu),dp)
            denom = dot_product(direction,spec%mu)
            if (abs(denom) > tiny(1.0_dp)) then
               if (spec%target_is_equality) then
                  delta = (spec%target-reward)/denom
                  temp = temp+delta*direction
               else if (abs(spec%target) > tiny(1.0_dp) .and. reward < spec%target) then
                  delta = (spec%target-reward)/denom
                  temp = temp+delta*direction
               end if
            end if
         end if
      end do
      if (spec%max_positions > 0 .and. spec%max_positions < size(temp)) then
         call enforce_cardinality(temp,spec,weights,ierr)
      else if (spec%leverage > 0.0_dp) then
         weights = min(max(temp,spec%lb),spec%ub)
         if (sum(abs(weights)) > tiny(1.0_dp)) weights = spec%leverage*weights/sum(abs(weights))
         weights = min(max(weights,spec%lb),spec%ub)
         ierr = 0
      else
         call project_box_budget(temp,spec%lb,spec%ub,spec%budget,weights,ierr)
      end if
      if (ierr /= 0) info = ierr
   end subroutine repair_weights

   subroutine enforce_cardinality(x,spec,weights,info)
      real(dp), intent(in) :: x(:)
      type(parma_spec), intent(in) :: spec
      real(dp), intent(out) :: weights(:)
      integer, intent(out) :: info
      real(dp), allocatable :: score(:),lbx(:),ubx(:),trial(:)
      logical, allocatable :: active(:)
      real(dp) :: scale
      integer :: i,idx

      allocate(score(size(x)),lbx(size(x)),ubx(size(x)),trial(size(x)),active(size(x)))
      score = abs(x)
      active = .false.
      do i = 1,min(spec%max_positions,size(x))
         idx = maxloc(score,dim=1)
         active(idx) = .true.
         score(idx) = -1.0_dp
      end do
      lbx = merge(spec%lb,0.0_dp,active)
      ubx = merge(spec%ub,0.0_dp,active)
      trial = merge(x,0.0_dp,active)
      if (spec%leverage > 0.0_dp) then
         weights = min(max(trial,lbx),ubx)
         if (sum(abs(weights)) > tiny(1.0_dp)) then
            scale = spec%leverage/sum(abs(weights))
            weights = min(max(scale*weights,lbx),ubx)
         end if
         info = 0
      else
         call project_box_budget(trial,lbx,ubx,spec%budget,weights,info)
      end if
   end subroutine enforce_cardinality

   function spec_risk_value(weights,spec,var_level,dar_level) result(value)
      real(dp), intent(in) :: weights(:)
      type(parma_spec), intent(in) :: spec
      real(dp), intent(out), optional :: var_level,dar_level
      real(dp) :: value,varx,darx

      varx = 0.0_dp
      darx = 0.0_dp
      if (spec%risk == risk_ev .and. allocated(spec%benchmark_cov)) then
         value = abs(benchmark_variance(weights,spec%cov,spec%benchmark_cov))
      else if (allocated(spec%data)) then
         if (allocated(spec%benchmark)) then
            value = risk_value(weights,spec%data,spec%risk,spec%alpha,spec%moment, &
               spec%threshold,benchmark=spec%benchmark,covariance=spec%cov, &
               lpm_legacy=spec%lpm_legacy,var_level=varx,dar_level=darx)
         else
            value = risk_value(weights,spec%data,spec%risk,spec%alpha,spec%moment, &
               spec%threshold,covariance=spec%cov,lpm_legacy=spec%lpm_legacy, &
               var_level=varx,dar_level=darx)
         end if
      else
         value = risk_value(weights,risk=spec%risk,alpha=spec%alpha,moment=spec%moment, &
            threshold=spec%threshold,covariance=spec%cov,lpm_legacy=spec%lpm_legacy, &
            var_level=varx,dar_level=darx)
      end if
      if (present(var_level)) var_level = varx
      if (present(dar_level)) dar_level = darx
   end function spec_risk_value

   subroutine parmasolve(spec, solution, options)
      type(parma_spec), intent(in) :: spec
      type(parma_port), intent(out) :: solution
      type(parma_options), intent(in), optional :: options
      type(parma_options) :: opt
      type(cmaes_result) :: cres
      real(dp), allocatable :: x0(:), repaired(:)
      real(dp) :: varx, darx
      integer :: info, m
      character(len=160) :: message

      opt = parma_options()
      if (present(options)) opt = options
      call validate_spec(spec,info,message)
      if (info /= 0) then
         solution%status = info
         solution%message = message
         return
      end if
      m = size(spec%mu)
      allocate(x0(m),repaired(m),solution%weights(m))
      if (allocated(spec%initial)) then
         x0 = spec%initial
      else
         x0 = spec%budget/real(m,dp)
      end if
      call repair_weights(x0,spec,repaired,info)
      x0 = repaired

      if (spec%risk == risk_ev .and. spec%objective == solve_min_risk .and. &
          .not. allocated(spec%ineq_a) .and. .not. allocated(spec%benchmark_cov) .and. &
          spec%leverage <= 0.0_dp .and. spec%max_positions == 0 .and. &
          spec%turnover_limit >= huge(1.0_dp)/2.0_dp .and. &
          spec%variance_limit >= huge(1.0_dp)/2.0_dp) then
         call solve_variance_projected(spec,x0,repaired,solution%iterations,opt)
         cres%status = 0
         cres%value = penalized_objective(repaired)
      else
         call cmaes_minimize(penalized_objective,x0,cres,opt,spec%lb,spec%ub)
         call repair_weights(cres%x,spec,repaired,info)
         if (opt%polish) call coordinate_polish(repaired,spec,opt)
         solution%iterations = cres%iterations
      end if
      solution%weights = repaired
      solution%reward = reward_value(repaired,spec)
      solution%risk = spec_risk_value(repaired,spec,varx,darx)
      solution%objective = raw_objective(repaired)
      solution%var_level = varx
      solution%dar_level = darx
      if (constraints_feasible(repaired,spec,max(1.0e-5_dp,sqrt(opt%tol)))) then
         solution%status = 0
         solution%message = 'converged'
      else
         solution%status = 2
         solution%message = 'candidate returned with residual constraint violation'
      end if

   contains

      function raw_objective(weights) result(value)
         real(dp), intent(in) :: weights(:)
         real(dp) :: value, rsk, rew, excess

         rew = reward_value(weights,spec)
         rsk = spec_risk_value(weights,spec)
         select case (spec%objective)
         case (solve_min_risk)
            value = rsk
         case (solve_max_reward)
            value = -rew
         case (solve_max_ratio)
            excess = rew-spec%ratio_rf
            if (excess <= 1.0e-12_dp) then
               value = 1.0e6_dp+abs(excess)*1.0e6_dp+rsk
            else
               value = rsk/excess
            end if
         case (solve_utility)
            value = spec%risk_aversion*rsk-rew
         case default
            value = rsk
         end select
      end function raw_objective

      function penalized_objective(weights) result(value)
         real(dp), intent(in) :: weights(:)
         real(dp) :: value
         value = raw_objective(weights)+opt%penalty*constraint_penalty(weights,spec)
      end function penalized_objective

   end subroutine parmasolve

   subroutine solve_variance_projected(spec, x0, weights, iterations, options)
      type(parma_spec), intent(in) :: spec
      real(dp), intent(in) :: x0(:)
      real(dp), intent(out) :: weights(:)
      integer, intent(out) :: iterations
      type(parma_options), intent(in) :: options
      real(dp), allocatable :: candidate(:), gradient(:), repaired(:)
      real(dp) :: step, f, fnew, reward, shortfall
      integer :: iter, info

      allocate(candidate(size(x0)),gradient(size(x0)),repaired(size(x0)))
      weights = x0
      step = options%step0/max(1.0_dp,maxval(abs(spec%cov)))
      f = quadratic_variance(weights,spec%cov)
      do iter = 1, options%max_iter
         gradient = 2.0_dp*matmul(spec%cov,weights)
         if (abs(spec%target) > tiny(1.0_dp)) then
            reward = dot_product(weights,spec%mu)
            if (spec%target_is_equality) then
               gradient = gradient+2.0_dp*options%penalty*(reward-spec%target)*spec%mu
            else
               shortfall = max(spec%target-reward,0.0_dp)
               if (shortfall > 0.0_dp) gradient = gradient-2.0_dp*options%penalty*shortfall*spec%mu
            end if
         end if
         candidate = weights-step*gradient
         call repair_weights(candidate,spec,repaired,info)
         fnew = quadratic_variance(repaired,spec%cov)+options%penalty*constraint_penalty(repaired,spec)
         if (fnew <= f) then
            if (vector_norm(repaired-weights) <= options%tol*(1.0_dp+vector_norm(weights))) then
               weights = repaired
               exit
            end if
            weights = repaired
            f = fnew
            step = min(step*1.05_dp,1.0_dp)
         else
            step = step*0.5_dp
            if (step < 1.0e-14_dp) exit
         end if
      end do
      iterations = min(iter,options%max_iter)
   end subroutine solve_variance_projected

   subroutine coordinate_polish(weights,spec,options)
      real(dp), intent(inout) :: weights(:)
      type(parma_spec), intent(in) :: spec
      type(parma_options), intent(in) :: options
      real(dp), allocatable :: candidate(:), repaired(:)
      real(dp) :: step, best, value
      integer :: i,j,iter,info

      allocate(candidate(size(weights)),repaired(size(weights)))
      best = portfolio_objective(weights,spec)+options%penalty*constraint_penalty(weights,spec)
      step = 0.05_dp*max(1.0_dp,maxval(spec%ub-spec%lb))
      do iter = 1, 100
         do i = 1, size(weights)
            do j = 1, 2
               candidate = weights
               if (j == 1) then
                  candidate(i) = candidate(i)+step
               else
                  candidate(i) = candidate(i)-step
               end if
               call repair_weights(candidate,spec,repaired,info)
               value = portfolio_objective(repaired,spec)+options%penalty*constraint_penalty(repaired,spec)
               if (value < best) then
                  weights = repaired
                  best = value
               end if
            end do
         end do
         step = step*0.7_dp
         if (step < options%tol) exit
      end do
   end subroutine coordinate_polish

   function portfolio_objective(weights,spec) result(value)
      real(dp), intent(in) :: weights(:)
      type(parma_spec), intent(in) :: spec
      real(dp) :: value, rsk, rew, excess

      rew = reward_value(weights,spec)
      rsk = spec_risk_value(weights,spec)
      select case (spec%objective)
      case (solve_min_risk)
         value = rsk
      case (solve_max_reward)
         value = -rew
      case (solve_max_ratio)
         excess = rew-spec%ratio_rf
         if (excess <= 1.0e-12_dp) then
            value = 1.0e6_dp+abs(excess)*1.0e6_dp+rsk
         else
            value = rsk/excess
         end if
      case (solve_utility)
         value = spec%risk_aversion*rsk-rew
      case default
         value = rsk
      end select
   end function portfolio_objective

   subroutine parmafrontier(spec, targets, frontier, options)
      type(parma_spec), intent(in) :: spec
      real(dp), intent(in) :: targets(:)
      type(parma_port), allocatable, intent(out) :: frontier(:)
      type(parma_options), intent(in), optional :: options
      type(parma_spec) :: work
      integer :: i

      allocate(frontier(size(targets)))
      do i = 1, size(targets)
         work = spec
         work%target = targets(i)
         if (present(options)) then
            call parmasolve(work,frontier(i),options)
         else
            call parmasolve(work,frontier(i))
         end if
      end do
   end subroutine parmafrontier

   subroutine parmautility2(m1,m2,risk_aversion,budget,lb,ub,solution,options)
      real(dp), intent(in) :: m1(:),m2(:,:),risk_aversion,budget,lb(:),ub(:)
      type(parma_port), intent(out) :: solution
      type(parma_options), intent(in), optional :: options
      type(parma_options) :: opt
      type(cmaes_result) :: cres
      real(dp), allocatable :: x0(:),x(:)
      integer :: info

      opt = parma_options()
      if (present(options)) opt = options
      allocate(x0(size(m1)),x(size(m1)),solution%weights(size(m1)))
      x0 = budget/real(size(m1),dp)
      call cmaes_minimize(objective,x0,cres,opt,lb,ub)
      call project_box_budget(cres%x,lb,ub,budget,x,info)
      solution%weights = x
      solution%reward = dot_product(x,m1)
      solution%objective = cara2_value(x,risk_aversion,m1,m2)
      solution%risk = dot_product(x,matmul(m2,x))
      solution%iterations = cres%iterations
      solution%status = cres%status
      solution%message = 'CARA second-moment optimization'
   contains
      function objective(w) result(value)
         real(dp), intent(in) :: w(:)
         real(dp) :: value
         value = cara2_value(w,risk_aversion,m1,m2)+opt%penalty*(sum(w)-budget)**2
      end function objective
   end subroutine parmautility2

   subroutine parmautility4(m1,m2,m3,m4,risk_aversion,budget,lb,ub,solution,options)
      real(dp), intent(in) :: m1(:),m2(:,:),m3(:,:),m4(:,:),risk_aversion,budget,lb(:),ub(:)
      type(parma_port), intent(out) :: solution
      type(parma_options), intent(in), optional :: options
      type(parma_options) :: opt
      type(cmaes_result) :: cres
      real(dp), allocatable :: x0(:),x(:)
      integer :: info

      opt = parma_options()
      if (present(options)) opt = options
      allocate(x0(size(m1)),x(size(m1)),solution%weights(size(m1)))
      x0 = budget/real(size(m1),dp)
      call cmaes_minimize(objective,x0,cres,opt,lb,ub)
      call project_box_budget(cres%x,lb,ub,budget,x,info)
      solution%weights = x
      solution%reward = dot_product(x,m1)
      solution%objective = cara4_value(x,risk_aversion,m1,m2,m3,m4)
      solution%risk = dot_product(x,matmul(m2,x))
      solution%iterations = cres%iterations
      solution%status = cres%status
      solution%message = 'CARA fourth-moment optimization'
   contains
      function objective(w) result(value)
         real(dp), intent(in) :: w(:)
         real(dp) :: value
         value = cara4_value(w,risk_aversion,m1,m2,m3,m4)+opt%penalty*(sum(w)-budget)**2
      end function objective
   end subroutine parmautility4

end module parma_portfolio
