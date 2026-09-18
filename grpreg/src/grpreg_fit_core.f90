module grpreg_fit_core
   use grpreg_kinds, only : dp
   use grpreg_math, only : soft_threshold, firm_threshold, scad_threshold, mcp_penalty, dmcp, logistic, log1pexp
   implicit none
   private
   public :: fit_path_gaussian, fit_path_glm, fit_path_cox, make_lambda_sequence

contains

   subroutine make_lambda_sequence(x, score, group, multiplier, alpha, bilevel, nlambda, lambda_min, log_lambda, lambda)
      real(dp), intent(in) :: x(:,:) !! Working design matrix with observations in rows.
      real(dp), intent(in) :: score(:) !! Null-model score or residual vector used to determine lambda maximum.
      integer, intent(in) :: group(:) !! Working predictor groups, zero for unpenalized predictors.
      real(dp), intent(in) :: multiplier(:) !! Penalty multipliers for positive groups.
      real(dp), intent(in) :: alpha !! Fraction of regularization assigned to the sparsity penalty.
      logical, intent(in) :: bilevel !! True uses maximum coordinate gradient; false uses group gradient norm.
      integer, intent(in) :: nlambda !! Number of requested regularization values, at least two.
      real(dp), intent(in) :: lambda_min !! Smallest lambda as a fraction of lambda maximum; zero appends exact zero.
      logical, intent(in) :: log_lambda !! True creates a log-spaced path; false creates a linear path.
      real(dp), allocatable, intent(out) :: lambda(:) !! Generated regularization sequence from largest to smallest.
      real(dp) :: z, zmax, lo, t
      real(dp), allocatable :: zg(:)
      integer :: g, j, l, ng, k
      ng = maxval(group)
      zmax = 0.0_dp
      do g = 1, ng
         if (bilevel) then
            do j = 1, size(group)
               if (group(j) /= g) cycle
               z = abs(dot_product(x(:,j), score))/multiplier(g)/real(size(x,1),dp)
               zmax = max(zmax, z)
            end do
         else
            allocate(zg(count(group == g)))
            k = 0
            do j = 1, size(group)
               if (group(j) /= g) cycle
               k = k + 1
               zg(k) = dot_product(x(:,j), score)/real(size(x,1),dp)
            end do
            z = sqrt(sum(zg*zg))/multiplier(g)
            zmax = max(zmax, z)
            deallocate(zg)
         end if
      end do
      zmax = zmax/max(alpha, tiny(1.0_dp))
      if (zmax <= 0.0_dp) zmax = 1.0_dp
      allocate(lambda(nlambda))
      if (lambda_min <= tiny(1.0_dp)) then
         lo = 0.001_dp*zmax
         do l = 1, nlambda - 1
            t = real(l - 1, dp)/real(max(1, nlambda - 2), dp)
            if (log_lambda) then
               lambda(l) = exp(log(zmax) + t*(log(lo) - log(zmax)))
            else
               lambda(l) = zmax + t*(lo - zmax)
            end if
         end do
         lambda(nlambda) = 0.0_dp
      else
         lo = lambda_min*zmax
         do l = 1, nlambda
            t = real(l - 1, dp)/real(max(1, nlambda - 1), dp)
            if (log_lambda) then
               lambda(l) = exp(log(zmax) + t*(log(lo) - log(zmax)))
            else
               lambda(l) = zmax + t*(lo - zmax)
            end if
         end do
      end if
   end subroutine make_lambda_sequence

   subroutine fit_path_gaussian(x, y, group, multiplier, penalty, lambda, alpha, gamma, tau, bridge_delta, eps, max_iter, &
      beta, intercept, deviance, df, iter, eta)
      real(dp), intent(in) :: x(:,:) !! Preprocessed design matrix with unit feature norms.
      real(dp), intent(in) :: y(:) !! Centered Gaussian response vector.
      integer, intent(in) :: group(:) !! Working group labels, with zero denoting unpenalized predictors.
      real(dp), intent(in) :: multiplier(:) !! Penalty multiplier for every positive group.
      character(len=*), intent(in) :: penalty !! One of grLasso, grMCP, grSCAD, gel, cMCP, or gBridge.
      real(dp), intent(in) :: lambda(:) !! Regularization path in fitting order.
      real(dp), intent(in) :: alpha !! Sparsity-versus-ridge mixing proportion in (0,1].
      real(dp), intent(in) :: gamma !! Nonconvex penalty concavity parameter.
      real(dp), intent(in) :: tau !! Group exponential lasso tuning parameter.
      real(dp), intent(in) :: bridge_delta !! Positive lower bound used to delete numerically tiny group-bridge groups.
      real(dp), intent(in) :: eps !! Coefficient-change convergence threshold.
      integer, intent(in) :: max_iter !! Maximum coordinate sweeps over the full path.
      real(dp), allocatable, intent(out) :: beta(:,:) !! Working-scale coefficient path, predictors by lambda.
      real(dp), allocatable, intent(out) :: intercept(:) !! Working-scale intercepts, zero for centered Gaussian fits.
      real(dp), allocatable, intent(out) :: deviance(:) !! Residual sum of squares at each lambda.
      real(dp), allocatable, intent(out) :: df(:) !! Upstream-style effective degrees-of-freedom approximation.
      integer, allocatable, intent(out) :: iter(:) !! Number of coordinate sweeps for each lambda.
      real(dp), allocatable, intent(out) :: eta(:,:) !! Linear predictor path on the working scale.
      real(dp), allocatable :: b(:), r(:), zvec(:)
      real(dp) :: max_change, old, newb, z, znorm, len, l1, l2, sdy, sgroup, ljk
      integer :: n, p, nlp, l, j, g, k, ng, total_iter
      logical :: bilevel
      n = size(x,1)
      p = size(x,2)
      nlp = size(lambda)
      ng = maxval(group)
      bilevel = trim(penalty) == 'gel' .or. trim(penalty) == 'cMCP' .or. trim(penalty) == 'gBridge'
      allocate(beta(p,nlp), intercept(nlp), deviance(nlp), df(nlp), iter(nlp), eta(n,nlp), b(p), r(n))
      beta = 0.0_dp
      intercept = 0.0_dp
      deviance = 0.0_dp
      df = 0.0_dp
      iter = 0
      eta = 0.0_dp
      b = 0.0_dp
      r = y
      sdy = sqrt(sum(y*y)/real(n,dp))
      total_iter = 0
      if (trim(penalty) == 'gBridge') call initialize_gaussian_ls(x, y, b, group)
      r = y - matmul(x,b)
      do l = 1, nlp
         if (l > 1) b = beta(:,l-1)
         r = y - matmul(x,b)
         do while (total_iter < max_iter)
            iter(l) = iter(l) + 1
            total_iter = total_iter + 1
            max_change = 0.0_dp
            df(l) = 0.0_dp
            do j = 1, p
               if (group(j) /= 0) cycle
               old = b(j)
               z = dot_product(x(:,j),r)/real(n,dp) + old
               b(j) = z
               r = r - x(:,j)*(b(j)-old)
               max_change = max(max_change, abs(b(j)-old))
               df(l) = df(l) + 1.0_dp
            end do
            if (.not. bilevel) then
               do g = 1, ng
                  allocate(zvec(count(group == g)))
                  k = 0
                  do j = 1, p
                     if (group(j) /= g) cycle
                     k = k + 1
                     zvec(k) = dot_product(x(:,j),r)/real(n,dp) + b(j)
                  end do
                  znorm = sqrt(sum(zvec*zvec))
                  l1 = lambda(l)*multiplier(g)*alpha
                  l2 = lambda(l)*multiplier(g)*(1.0_dp-alpha)
                  if (trim(penalty) == 'grLasso') then
                     len = max(soft_threshold(znorm,l1),0.0_dp)/(1.0_dp+l2)
                  else if (trim(penalty) == 'grMCP') then
                     len = firm_threshold(znorm,l1,l2,gamma)
                  else
                     len = scad_threshold(znorm,l1,l2,gamma)
                  end if
                  if (znorm > 0.0_dp) then
                     k = 0
                     do j = 1, p
                        if (group(j) /= g) cycle
                        k = k + 1
                        old = b(j)
                        b(j) = len*zvec(k)/znorm
                        r = r - x(:,j)*(b(j)-old)
                        max_change = max(max_change,abs(b(j)-old))
                     end do
                     if (len > 0.0_dp) df(l) = df(l) + real(size(zvec),dp)*len/znorm
                  end if
                  deallocate(zvec)
               end do
            else
               do g = 1, ng
                  l1 = lambda(l)*multiplier(g)*alpha
                  l2 = lambda(l)*multiplier(g)*(1.0_dp-alpha)
                  sgroup = bilevel_group_measure(b, group, g, penalty, l1, gamma)
                  if (trim(penalty) == 'gBridge' .and. sgroup <= tiny(1.0_dp)) cycle
                  if (trim(penalty) == 'gBridge' .and. sgroup < bridge_delta) then
                     do j = 1, p
                        if (group(j) /= g) cycle
                        old = b(j)
                        b(j) = 0.0_dp
                        r = r+x(:,j)*old
                        max_change = max(max_change,abs(old))
                     end do
                     cycle
                  end if
                  do j = 1, p
                     if (group(j) /= g) cycle
                     old = b(j)
                     z = dot_product(x(:,j),r)/real(n,dp) + old
                     ljk = bilevel_threshold(old, sgroup, count(group == g), penalty, l1, gamma, tau)
                     newb = soft_threshold(z,ljk)/(1.0_dp+l2)
                     b(j) = newb
                     r = r - x(:,j)*(newb-old)
                     max_change = max(max_change,abs(newb-old))
                     if (abs(z) > tiny(1.0_dp)) df(l) = df(l) + abs(newb)/abs(z)
                     sgroup = update_bilevel_measure(sgroup, old, newb, penalty, l1, gamma)
                  end do
               end do
            end if
            if (max_change <= eps*max(sdy,1.0_dp)) exit
         end do
         beta(:,l) = b
         eta(:,l) = matmul(x,b)
         deviance(l) = sum((y-eta(:,l))**2)
      end do
   end subroutine fit_path_gaussian

   subroutine fit_path_glm(x, y, group, multiplier, family, penalty, lambda, alpha, gamma, tau, bridge_delta, eps, max_iter, &
      beta, intercept, deviance, df, iter, eta)
      real(dp), intent(in) :: x(:,:) !! Preprocessed design matrix with unit feature norms.
      real(dp), intent(in) :: y(:) !! Binomial zero-one outcomes or nonnegative Poisson counts.
      integer, intent(in) :: group(:) !! Working group labels, with zero denoting unpenalized predictors.
      real(dp), intent(in) :: multiplier(:) !! Penalty multiplier for every positive group.
      character(len=*), intent(in) :: family !! GLM family, either binomial or poisson.
      character(len=*), intent(in) :: penalty !! Group or bi-level penalty name.
      real(dp), intent(in) :: lambda(:) !! Regularization path in fitting order.
      real(dp), intent(in) :: alpha !! Sparsity-versus-ridge mixing proportion in (0,1].
      real(dp), intent(in) :: gamma !! Nonconvex penalty concavity parameter.
      real(dp), intent(in) :: tau !! Group exponential lasso tuning parameter.
      real(dp), intent(in) :: bridge_delta !! Positive lower bound used to delete numerically tiny group-bridge groups.
      real(dp), intent(in) :: eps !! Maximum coefficient-change convergence threshold.
      integer, intent(in) :: max_iter !! Maximum coordinate sweeps over the full path.
      real(dp), allocatable, intent(out) :: beta(:,:) !! Working-scale coefficient path.
      real(dp), allocatable, intent(out) :: intercept(:) !! Intercept path.
      real(dp), allocatable, intent(out) :: deviance(:) !! Model deviance at each lambda.
      real(dp), allocatable, intent(out) :: df(:) !! Effective degrees-of-freedom approximation.
      integer, allocatable, intent(out) :: iter(:) !! Number of coordinate sweeps for each lambda.
      real(dp), allocatable, intent(out) :: eta(:,:) !! Linear predictor path.
      real(dp), allocatable :: b(:), r(:), mu(:), zvec(:)
      real(dp) :: b0, ybar, v, max_change, old, z, znorm, len, l1, l2, sgroup, ljk, newb, shift
      integer :: n, p, nlp, l, j, g, k, ng, total_iter
      logical :: bilevel
      n = size(x,1)
      p = size(x,2)
      nlp = size(lambda)
      ng = maxval(group)
      bilevel = trim(penalty) == 'gel' .or. trim(penalty) == 'cMCP' .or. trim(penalty) == 'gBridge'
      allocate(beta(p,nlp), intercept(nlp), deviance(nlp), df(nlp), iter(nlp), eta(n,nlp), b(p), r(n), mu(n))
      beta = 0.0_dp
      intercept = 0.0_dp
      deviance = 0.0_dp
      df = 0.0_dp
      iter = 0
      eta = 0.0_dp
      b = 0.0_dp
      ybar = sum(y)/real(n,dp)
      if (trim(family) == 'binomial') then
         ybar = min(max(ybar,1.0e-8_dp),1.0_dp-1.0e-8_dp)
         b0 = log(ybar/(1.0_dp-ybar))
      else
         b0 = log(max(ybar,1.0e-12_dp))
      end if
      total_iter = 0
      do l = 1, nlp
         if (l > 1) then
            b = beta(:,l-1)
            b0 = intercept(l-1)
         end if
         eta(:,l) = b0 + matmul(x,b)
         do while (total_iter < max_iter)
            iter(l) = iter(l) + 1
            total_iter = total_iter + 1
            if (trim(family) == 'binomial') then
               mu = logistic(eta(:,l))
               v = 0.25_dp
            else
               mu = exp(min(eta(:,l),40.0_dp))
               v = max(maxval(mu),1.0e-8_dp)
            end if
            r = (y-mu)/v
            max_change = 0.0_dp
            shift = sum(r)/real(n,dp)
            b0 = b0 + shift
            r = r - shift
            eta(:,l) = eta(:,l) + shift
            max_change = abs(shift)
            df(l) = 1.0_dp
            do j = 1, p
               if (group(j) /= 0) cycle
               old = b(j)
               z = dot_product(x(:,j),r)/real(n,dp) + old
               b(j) = z
               shift = b(j)-old
               r = r - x(:,j)*shift
               eta(:,l) = eta(:,l) + x(:,j)*shift
               max_change = max(max_change,abs(shift))
               df(l) = df(l) + 1.0_dp
            end do
            if (.not. bilevel) then
               do g = 1, ng
                  allocate(zvec(count(group == g)))
                  k = 0
                  do j = 1, p
                     if (group(j) /= g) cycle
                     k = k + 1
                     zvec(k) = dot_product(x(:,j),r)/real(n,dp) + b(j)
                  end do
                  znorm = sqrt(sum(zvec*zvec))
                  l1 = lambda(l)*multiplier(g)*alpha
                  l2 = lambda(l)*multiplier(g)*(1.0_dp-alpha)
                  if (trim(penalty) == 'grLasso') then
                     len = max(soft_threshold(v*znorm,l1),0.0_dp)/(v*(1.0_dp+l2))
                  else if (trim(penalty) == 'grMCP') then
                     len = firm_threshold(v*znorm,l1,l2,gamma)/v
                  else
                     len = scad_threshold(v*znorm,l1,l2,gamma)/v
                  end if
                  if (znorm > 0.0_dp) then
                     k = 0
                     do j = 1, p
                        if (group(j) /= g) cycle
                        k = k + 1
                        old = b(j)
                        b(j) = len*zvec(k)/znorm
                        shift = b(j)-old
                        r = r - x(:,j)*shift
                        eta(:,l) = eta(:,l) + x(:,j)*shift
                        max_change = max(max_change,abs(shift))
                     end do
                     if (len > 0.0_dp) df(l) = df(l) + real(size(zvec),dp)*len/znorm
                  end if
                  deallocate(zvec)
               end do
            else
               do g = 1, ng
                  l1 = lambda(l)*multiplier(g)*alpha
                  l2 = lambda(l)*multiplier(g)*(1.0_dp-alpha)
                  sgroup = bilevel_group_measure(b/v, group, g, penalty, l1, gamma)
                  if (trim(penalty) == 'gBridge' .and. sgroup <= tiny(1.0_dp)) cycle
                  if (trim(penalty) == 'gBridge' .and. sgroup < bridge_delta) then
                     do j = 1, p
                        if (group(j) /= g) cycle
                        old = b(j)
                        b(j) = 0.0_dp
                        shift = -old
                        r = r-x(:,j)*shift
                        eta(:,l) = eta(:,l)+x(:,j)*shift
                        max_change = max(max_change,abs(old))
                     end do
                     cycle
                  end if
                  do j = 1, p
                     if (group(j) /= g) cycle
                     old = b(j)
                     z = dot_product(x(:,j),r)/real(n,dp) + old
                     ljk = bilevel_threshold(old, sgroup, count(group == g), penalty, l1, gamma, tau)
                     newb = soft_threshold(v*z,ljk)/(v*(1.0_dp+l2))
                     b(j) = newb
                     shift = newb-old
                     r = r - x(:,j)*shift
                     eta(:,l) = eta(:,l) + x(:,j)*shift
                     max_change = max(max_change,abs(shift))
                     if (abs(z) > tiny(1.0_dp)) df(l) = df(l) + abs(newb)/abs(z)
                     sgroup = update_bilevel_measure(sgroup, old, newb, penalty, l1, gamma)
                  end do
               end do
            end if
            if (max_change < eps) exit
         end do
         beta(:,l) = b
         intercept(l) = b0
         if (trim(family) == 'binomial') then
            deviance(l) = binomial_deviance(y,eta(:,l))
         else
            deviance(l) = poisson_deviance(y,eta(:,l))
         end if
      end do
   end subroutine fit_path_glm

   subroutine fit_path_cox(x, event, group, multiplier, penalty, lambda, alpha, gamma, tau, eps, max_iter, &
      beta, deviance, df, iter, eta)
      real(dp), intent(in) :: x(:,:) !! Preprocessed design matrix ordered by increasing follow-up time.
      real(dp), intent(in) :: event(:) !! Event indicators, zero for censoring and one for failures.
      integer, intent(in) :: group(:) !! Working group labels, zero for unpenalized predictors.
      real(dp), intent(in) :: multiplier(:) !! Penalty multiplier for every positive group.
      character(len=*), intent(in) :: penalty !! Group or bi-level penalty name.
      real(dp), intent(in) :: lambda(:) !! Regularization path in fitting order.
      real(dp), intent(in) :: alpha !! Sparsity-versus-ridge mixing proportion.
      real(dp), intent(in) :: gamma !! Nonconvex penalty concavity parameter.
      real(dp), intent(in) :: tau !! Group exponential lasso tuning parameter.
      real(dp), intent(in) :: eps !! Maximum coefficient-change convergence threshold.
      integer, intent(in) :: max_iter !! Maximum coordinate sweeps over the full path.
      real(dp), allocatable, intent(out) :: beta(:,:) !! Working-scale Cox coefficient path.
      real(dp), allocatable, intent(out) :: deviance(:) !! Minus twice the Cox partial log likelihood.
      real(dp), allocatable, intent(out) :: df(:) !! Effective degrees-of-freedom approximation.
      integer, allocatable, intent(out) :: iter(:) !! Number of coordinate sweeps for each lambda.
      real(dp), allocatable, intent(out) :: eta(:,:) !! Cox linear predictor path.
      real(dp), allocatable :: b(:), haz(:), risk(:), h(:), r(:), zvec(:), vj(:)
      real(dp) :: max_change, old, shift, l1, l2, znorm, len, z, sgroup, ljk, newb, xwr
      integer :: n, p, nlp, l, j, g, k, i, ng, total_iter
      logical :: bilevel
      n = size(x,1)
      p = size(x,2)
      nlp = size(lambda)
      ng = maxval(group)
      bilevel = trim(penalty) == 'gel' .or. trim(penalty) == 'cMCP' .or. trim(penalty) == 'gBridge'
      allocate(beta(p,nlp), deviance(nlp), df(nlp), iter(nlp), eta(n,nlp), b(p), haz(n), risk(n), h(n), r(n), vj(p))
      beta = 0.0_dp
      deviance = 0.0_dp
      df = 0.0_dp
      iter = 0
      eta = 0.0_dp
      b = 0.0_dp
      total_iter = 0
      do l = 1, nlp
         if (l > 1) b = beta(:,l-1)
         eta(:,l) = matmul(x,b)
         do while (total_iter < max_iter)
            iter(l) = iter(l) + 1
            total_iter = total_iter + 1
            haz = exp(min(eta(:,l),40.0_dp))
            risk(n) = haz(n)
            do i = n - 1, 1, -1
               risk(i) = risk(i+1) + haz(i)
            end do
            h(1) = event(1)/max(risk(1),tiny(1.0_dp))
            do i = 2, n
               h(i) = h(i-1) + event(i)/max(risk(i),tiny(1.0_dp))
            end do
            h = h*haz
            r = event-h
            max_change = 0.0_dp
            df(l) = 0.0_dp
            do j = 1, p
               if (group(j) /= 0) cycle
               old = b(j)
               z = dot_product(x(:,j),r)/real(n,dp) + old
               b(j) = z
               shift = b(j)-old
               r = r - x(:,j)*shift
               eta(:,l) = eta(:,l) + x(:,j)*shift
               max_change = max(max_change,abs(shift))
               df(l) = df(l) + 1.0_dp
            end do
            if (.not. bilevel) then
               do g = 1, ng
                  allocate(zvec(count(group == g)))
                  k = 0
                  do j = 1, p
                     if (group(j) /= g) cycle
                     k = k + 1
                     zvec(k) = dot_product(x(:,j),r)/real(n,dp) + b(j)
                  end do
                  znorm = sqrt(sum(zvec*zvec))
                  l1 = lambda(l)*multiplier(g)*alpha
                  l2 = lambda(l)*multiplier(g)*(1.0_dp-alpha)
                  if (trim(penalty) == 'grLasso') then
                     len = max(soft_threshold(znorm,l1),0.0_dp)/(1.0_dp+l2)
                  else if (trim(penalty) == 'grMCP') then
                     len = firm_threshold(znorm,l1,l2,gamma)
                  else
                     len = scad_threshold(znorm,l1,l2,gamma)
                  end if
                  if (znorm > 0.0_dp) then
                     k = 0
                     do j = 1, p
                        if (group(j) /= g) cycle
                        k = k + 1
                        old = b(j)
                        b(j) = len*zvec(k)/znorm
                        shift = b(j)-old
                        r = r - x(:,j)*shift
                        eta(:,l) = eta(:,l) + x(:,j)*shift
                        max_change = max(max_change,abs(shift))
                     end do
                     if (len > 0.0_dp) df(l) = df(l) + real(size(zvec),dp)*len/znorm
                  end if
                  deallocate(zvec)
               end do
            else
               do j = 1, p
                  vj(j) = max(sum(h*x(:,j)*x(:,j))/real(n,dp),1.0e-8_dp)
               end do
               do g = 1, ng
                  l1 = lambda(l)*multiplier(g)*alpha
                  l2 = lambda(l)*multiplier(g)*(1.0_dp-alpha)
                  sgroup = 0.0_dp
                  do j = 1, p
                     if (group(j) /= g) cycle
                     if (trim(penalty) == 'cMCP') then
                        sgroup = sgroup + mcp_penalty(b(j)/vj(j),sqrt(max(l1,0.0_dp)),gamma)
                     else
                        sgroup = sgroup + abs(b(j))/vj(j)
                     end if
                  end do
                  do j = 1, p
                     if (group(j) /= g) cycle
                     old = b(j)
                     xwr = sum(x(:,j)*r*h)/real(n,dp)
                     z = xwr + vj(j)*old
                     ljk = bilevel_threshold(old,sgroup,count(group == g),penalty,l1,gamma,tau)
                     newb = soft_threshold(z,ljk)/(vj(j)*(1.0_dp+l2))
                     b(j) = newb
                     shift = newb-old
                     r = r - x(:,j)*shift
                     eta(:,l) = eta(:,l) + x(:,j)*shift
                     max_change = max(max_change,abs(shift))
                     if (abs(z) > tiny(1.0_dp)) df(l) = df(l) + abs(newb)/abs(z)
                     sgroup = update_bilevel_measure(sgroup,old,newb,penalty,l1,gamma)
                  end do
               end do
            end if
            if (max_change < eps) exit
         end do
         beta(:,l) = b
         deviance(l) = -2.0_dp*cox_loglik(event,eta(:,l))
      end do
   end subroutine fit_path_cox

   subroutine initialize_gaussian_ls(x, y, beta, group)
      real(dp), intent(in) :: x(:,:) !! Standardized design matrix used for group-bridge warm start.
      real(dp), intent(in) :: y(:) !! Centered response vector used for group-bridge warm start.
      real(dp), intent(inout) :: beta(:) !! Coefficient vector overwritten by one coordinate least-squares pass.
      integer, intent(in) :: group(:) !! Working group labels used only to preserve API symmetry.
      real(dp), allocatable :: r(:)
      real(dp) :: old
      integer :: j
      allocate(r(size(y)))
      r = y
      beta = 0.0_dp
      do j = 1, size(beta)
         old = dot_product(x(:,j),r)/real(size(y),dp)
         beta(j) = old
         r = r - old*x(:,j)
      end do
      if (size(group) < 0) error stop 'unreachable'
   end subroutine initialize_gaussian_ls

   pure real(dp) function bilevel_group_measure(beta, group, g, penalty, lambda1, gamma) result(value)
      real(dp), intent(in) :: beta(:) !! Current coefficient vector on the working scale.
      integer, intent(in) :: group(:) !! Working group labels for beta.
      integer, intent(in) :: g !! Positive group whose inner penalty is summarized.
      character(len=*), intent(in) :: penalty !! Bi-level penalty name.
      real(dp), intent(in) :: lambda1 !! Group-specific sparsity tuning parameter.
      real(dp), intent(in) :: gamma !! Concavity parameter for MCP or bridge penalty.
      integer :: j
      real(dp) :: lam
      value = 0.0_dp
      lam = lambda1
      if (trim(penalty) == 'cMCP') lam = sqrt(max(lambda1,0.0_dp))
      do j = 1, size(beta)
         if (group(j) /= g) cycle
         if (trim(penalty) == 'cMCP') then
            value = value + mcp_penalty(beta(j),lam,gamma)
         else
            value = value + abs(beta(j))
         end if
      end do
   end function bilevel_group_measure

   pure real(dp) function bilevel_threshold(beta_j, sgroup, kgroup, penalty, lambda1, gamma, tau) result(value)
      real(dp), intent(in) :: beta_j !! Current coefficient within the active group.
      real(dp), intent(in) :: sgroup !! Current group-level sum of inner penalties.
      integer, intent(in) :: kgroup !! Number of coefficients in the current group.
      character(len=*), intent(in) :: penalty !! Bi-level penalty name.
      real(dp), intent(in) :: lambda1 !! Group-specific sparsity tuning parameter.
      real(dp), intent(in) :: gamma !! Concavity parameter for MCP or bridge penalty.
      real(dp), intent(in) :: tau !! Group exponential lasso decay parameter.
      real(dp) :: lam, outer_gamma
      value = 0.0_dp
      if (lambda1 <= 0.0_dp) return
      select case (trim(penalty))
      case ('gel')
         value = lambda1*exp(-tau/lambda1*sgroup)
      case ('cMCP')
         lam = sqrt(lambda1)
         outer_gamma = real(kgroup,dp)*gamma*lam*lam/(2.0_dp*lam)
         value = dmcp(sgroup,lam,outer_gamma)*dmcp(beta_j,lam,gamma)
      case ('gBridge')
         if (sgroup > 0.0_dp) value = lambda1*gamma*sgroup**(gamma-1.0_dp)
      end select
   end function bilevel_threshold

   pure real(dp) function update_bilevel_measure(sgroup, old, new, penalty, lambda1, gamma) result(value)
      real(dp), intent(in) :: sgroup !! Group-level inner-penalty sum before a coordinate update.
      real(dp), intent(in) :: old !! Previous coefficient value.
      real(dp), intent(in) :: new !! Updated coefficient value.
      character(len=*), intent(in) :: penalty !! Bi-level penalty name.
      real(dp), intent(in) :: lambda1 !! Group-specific sparsity tuning parameter.
      real(dp), intent(in) :: gamma !! Concavity parameter used by composite MCP.
      real(dp) :: lam
      value = sgroup
      if (trim(penalty) == 'cMCP') then
         lam = sqrt(max(lambda1,0.0_dp))
         value = sgroup + mcp_penalty(new,lam,gamma) - mcp_penalty(old,lam,gamma)
      else
         value = sgroup + abs(new) - abs(old)
      end if
   end function update_bilevel_measure

   pure real(dp) function binomial_deviance(y, eta) result(value)
      real(dp), intent(in) :: y(:) !! Zero-one outcomes.
      real(dp), intent(in) :: eta(:) !! Logistic linear predictors corresponding to y.
      integer :: i
      value = 0.0_dp
      do i = 1, size(y)
         value = value + 2.0_dp*(log1pexp(eta(i)) - y(i)*eta(i))
      end do
   end function binomial_deviance

   pure real(dp) function poisson_deviance(y, eta) result(value)
      real(dp), intent(in) :: y(:) !! Nonnegative Poisson outcomes.
      real(dp), intent(in) :: eta(:) !! Log-mean linear predictors corresponding to y.
      real(dp) :: mu
      integer :: i
      value = 0.0_dp
      do i = 1, size(y)
         mu = exp(min(eta(i),40.0_dp))
         if (y(i) > 0.0_dp) then
            value = value + 2.0_dp*(y(i)*log(y(i)/mu) + mu - y(i))
         else
            value = value + 2.0_dp*mu
         end if
      end do
   end function poisson_deviance

   pure real(dp) function cox_loglik(event, eta) result(value)
      real(dp), intent(in) :: event(:) !! Event indicators ordered by increasing follow-up time.
      real(dp), intent(in) :: eta(:) !! Cox linear predictors in the same order.
      real(dp) :: risk(size(event))
      integer :: i, n
      n = size(event)
      risk(n) = exp(min(eta(n),40.0_dp))
      do i = n - 1, 1, -1
         risk(i) = risk(i+1) + exp(min(eta(i),40.0_dp))
      end do
      value = sum(event*(eta-log(max(risk,tiny(1.0_dp)))))
   end function cox_loglik

end module grpreg_fit_core
