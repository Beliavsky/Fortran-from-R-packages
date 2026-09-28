module psych_correlations
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use r_kinds, only : dp
   use r_linalg, only : symmetric_eigen, thin_svd, determinant
   use psych_stats, only : finite_mean, finite_sd, finite_variance, finite_covariance, finite_correlation
   use psych_stats, only : covariance_matrix, correlation_matrix, trace_matrix
   use psych_distributions, only : chi_square_cdf, student_t_cdf, normal_cdf
      use psych_types, only : kmo_result, correlation_test_result
   implicit none
   private
   public :: is_correlation_matrix, is_covariance_matrix, correlation_to_distance
   public :: correct_correlation, range_correction, weighted_correlation_matrix
   public :: pseudoinverse, squared_multiple_correlations, kmo_statistic, partial_correlation
   public :: congruence_matrix, profile_distance_matrix, factor_fit_value, factor_model_matrix
   public :: bartlett_correlation_test, correlation_significance, paired_correlation_test
   public :: smooth_correlation, rv_coefficient, cohen_kappa_value, weighted_kappa_value
   public :: cohen_d_value, pooled_sd_value, mahalanobis_squared, mardia_statistics
   public :: auto_correlation, biserial_correlation, point_biserial_correlation

contains

   pure logical function is_correlation_matrix(r, tolerance) result(ok)
      real(dp), intent(in) :: r(:, :) !! Candidate square correlation matrix.
      real(dp), intent(in), optional :: tolerance !! Absolute symmetry/diagonal tolerance; default `1e-8`.
      real(dp) :: tol
      integer :: i
      tol = 1.0e-8_dp
      if (present(tolerance)) tol = tolerance
      ok = size(r,1) == size(r,2)
      if (.not. ok) return
      ok = maxval(abs(r-transpose(r))) <= tol
      if (.not. ok) return
      do i = 1, size(r,1)
         if (abs(r(i,i)-1.0_dp) > tol) then
            ok = .false.
            return
         end if
      end do
      ok = maxval(abs(r)) <= 1.0_dp+tol
   end function is_correlation_matrix

   pure logical function is_covariance_matrix(c, tolerance) result(ok)
      real(dp), intent(in) :: c(:, :) !! Candidate covariance matrix.
      real(dp), intent(in), optional :: tolerance !! Absolute symmetry tolerance; default `1e-8`.
      real(dp) :: tol
      integer :: i
      tol = 1.0e-8_dp
      if (present(tolerance)) tol = tolerance
      ok = size(c,1) == size(c,2)
      if (.not. ok) return
      ok = maxval(abs(c-transpose(c))) <= tol
      if (.not. ok) return
      do i = 1, size(c,1)
         if (c(i,i) < -tol) then
            ok = .false.
            return
         end if
      end do
   end function is_covariance_matrix

   pure subroutine correlation_to_distance(r, distance)
      real(dp), intent(in) :: r(:, :) !! Correlation matrix with entries in `[-1,1]`.
      real(dp), intent(out) :: distance(:, :) !! Chord distance matrix `sqrt(2*(1-r))`.
      distance = sqrt(max(0.0_dp,2.0_dp*(1.0_dp-r)))
   end subroutine correlation_to_distance

   pure subroutine correct_correlation(r, reliability, corrected)
      real(dp), intent(in) :: r(:, :) !! Observed correlation matrix.
      real(dp), intent(in) :: reliability(:) !! Reliability coefficients for each variable.
      real(dp), intent(out) :: corrected(:, :) !! Disattenuated correlation matrix.
      integer :: i, j
      corrected = r
      do i = 1, size(r,1)
         corrected(i,i) = reliability(i)
      end do
      do j = 2, size(r,2)
         do i = 1, j-1
            if (reliability(i) > 0.0_dp .and. reliability(j) > 0.0_dp) then
               corrected(i,j) = r(i,j)/sqrt(reliability(i)*reliability(j))
               corrected(j,i) = corrected(i,j)
            end if
         end do
      end do
   end subroutine correct_correlation

   pure real(dp) function range_correction(r, sdu, sdr, sdxu, sdxr, correction_case) result(value)
      real(dp), intent(in) :: r !! Restricted-range observed correlation.
      real(dp), intent(in) :: sdu !! Unrestricted standard deviation for the criterion or first variable.
      real(dp), intent(in) :: sdr !! Restricted standard deviation for the criterion or first variable.
      real(dp), intent(in), optional :: sdxu !! Unrestricted standard deviation for the second variable in case four.
      real(dp), intent(in), optional :: sdxr !! Restricted standard deviation for the second variable in case four.
      integer, intent(in), optional :: correction_case !! Thorndike-style case selector; default is two.
      integer :: mode
      mode = 2
      if (present(correction_case)) mode = correction_case
      select case (mode)
      case (1)
         value = sqrt(max(0.0_dp,1.0_dp-(sdr*sdr/(sdu*sdu))*(1.0_dp-r*r)))
      case (4)
         if (present(sdxu) .and. present(sdxr)) then
            value = r*(sdr/sdu)*(sdxr/sdxu)+sqrt(max(0.0_dp,(1.0_dp-(sdr/sdu)**2)*(1.0_dp-(sdxr/sdxu)**2)))
         else
            value = r
         end if
      case default
         value = r*sdu/(sdr*sqrt(max(tiny(1.0_dp),1.0_dp-r*r+r*r*(sdu*sdu/(sdr*sdr)))))
      end select
   end function range_correction

   pure subroutine weighted_correlation_matrix(x, weights, r)
      real(dp), intent(in) :: x(:, :) !! Observation-by-variable matrix without missing-data imputation.
      real(dp), intent(in) :: weights(:) !! Nonnegative observation weights with length equal to rows of `x`.
      real(dp), intent(out) :: r(:, :) !! Weighted Pearson correlation matrix.
      real(dp), allocatable :: means(:), sds(:)
      real(dp) :: sw, covij
      integer :: i, j, k, p
      p = size(x,2)
      allocate(means(p),sds(p))
      sw = sum(weights)
      if (sw <= 0.0_dp) then
         r = 0.0_dp
         return
      end if
      do j = 1, p
         means(j) = sum(weights*x(:,j))/sw
         sds(j) = sqrt(max(0.0_dp,sum(weights*(x(:,j)-means(j))**2)/sw))
      end do
      do j = 1, p
         do i = 1, p
            covij = 0.0_dp
            do k = 1, size(x,1)
               covij = covij+weights(k)*(x(k,i)-means(i))*(x(k,j)-means(j))
            end do
            covij = covij/sw
            if (sds(i) > 0.0_dp .and. sds(j) > 0.0_dp) then
               r(i,j) = covij/(sds(i)*sds(j))
            else
               r(i,j) = 0.0_dp
            end if
         end do
      end do
   end subroutine weighted_correlation_matrix

   subroutine pseudoinverse(a, inverse, tolerance, info)
      real(dp), intent(in) :: a(:, :) !! Matrix to pseudo-invert.
      real(dp), allocatable, intent(out) :: inverse(:, :) !! Moore-Penrose pseudoinverse with shape `(n,m)`.
      real(dp), intent(in), optional :: tolerance !! Singular-value cutoff; defaults to a size-scaled machine threshold.
      integer, intent(out) :: info !! Zero on successful SVD, otherwise the linear-algebra backend status.
      real(dp), allocatable :: u(:, :), singular(:), vt(:, :), scaled_v(:, :)
      real(dp) :: cutoff
      integer :: i, m, n, k
      m = size(a,1)
      n = size(a,2)
      call thin_svd(a,u,singular,vt,info)
      allocate(inverse(n,m))
      inverse = 0.0_dp
      if (info /= 0 .or. size(singular) == 0) return
      cutoff = real(max(m,n),dp)*epsilon(1.0_dp)*max(1.0_dp,singular(1))
      if (present(tolerance)) cutoff = max(0.0_dp,tolerance)
      k = size(singular)
      allocate(scaled_v(n,k))
      scaled_v = transpose(vt)
      do i = 1, k
         if (singular(i) > cutoff) then
            scaled_v(:,i) = scaled_v(:,i)/singular(i)
         else
            scaled_v(:,i) = 0.0_dp
         end if
      end do
      inverse = matmul(scaled_v,transpose(u))
   end subroutine pseudoinverse

   subroutine squared_multiple_correlations(r, smc, info)
      real(dp), intent(in) :: r(:, :) !! Correlation or covariance matrix.
      real(dp), intent(out) :: smc(:) !! Squared multiple correlation for each variable, clipped to `[0,1]`.
      integer, intent(out) :: info !! Zero when the pseudoinverse calculation succeeds.
      real(dp), allocatable :: rinv(:, :)
      integer :: i
      call pseudoinverse(r,rinv,info=info)
      if (info /= 0) then
         smc = 1.0_dp
         return
      end if
      do i = 1, size(smc)
         if (abs(rinv(i,i)) > tiny(1.0_dp)) then
            smc(i) = min(1.0_dp,max(0.0_dp,1.0_dp-1.0_dp/rinv(i,i)))
         else
            smc(i) = 1.0_dp
         end if
      end do
   end subroutine squared_multiple_correlations

   subroutine kmo_statistic(r, result, info)
      real(dp), intent(in) :: r(:, :) !! Correlation matrix.
      type(kmo_result), intent(out) :: result !! Overall and item-level Kaiser-Meyer-Olkin sampling adequacy.
      integer, intent(out) :: info !! Zero on successful matrix inversion/pseudoinversion.
      real(dp), allocatable :: q(:, :), work_r(:, :), work_q(:, :)
      real(dp) :: sum_r2, sum_q2
      integer :: i, j, p
      p = size(r,1)
      call pseudoinverse(r,q,info=info)
      allocate(result%item_msa(p),result%image(p,p),work_r(p,p),work_q(p,p))
      if (info /= 0) then
         result%msa = 0.0_dp
         result%item_msa = 0.0_dp
         result%image = 0.0_dp
         return
      end if
      do j = 1, p
         do i = 1, p
            if (q(i,i) > 0.0_dp .and. q(j,j) > 0.0_dp) then
               result%image(i,j) = q(i,j)/sqrt(q(i,i)*q(j,j))
            else
               result%image(i,j) = 0.0_dp
            end if
         end do
      end do
      work_r = r
      work_q = result%image
      do i = 1, p
         work_r(i,i) = 0.0_dp
         work_q(i,i) = 0.0_dp
      end do
      sum_r2 = sum(work_r*work_r)
      sum_q2 = sum(work_q*work_q)
      if (sum_r2+sum_q2 > 0.0_dp) then
         result%msa = sum_r2/(sum_r2+sum_q2)
      else
         result%msa = 0.0_dp
      end if
      do j = 1, p
         sum_r2 = sum(work_r(:,j)*work_r(:,j))
         sum_q2 = sum(work_q(:,j)*work_q(:,j))
         if (sum_r2+sum_q2 > 0.0_dp) then
            result%item_msa(j) = sum_r2/(sum_r2+sum_q2)
         else
            result%item_msa(j) = 0.0_dp
         end if
      end do
   end subroutine kmo_statistic

   subroutine partial_correlation(r, partial, info)
      real(dp), intent(in) :: r(:, :) !! Correlation matrix whose inverse defines partial correlations.
      real(dp), intent(out) :: partial(:, :) !! Full-order partial-correlation matrix.
      integer, intent(out) :: info !! Zero on successful pseudo-inversion.
      real(dp), allocatable :: q(:, :)
      integer :: i, j, p
      p = size(r,1)
      call pseudoinverse(r,q,info=info)
      if (info /= 0) then
         partial = 0.0_dp
         return
      end if
      do j = 1, p
         do i = 1, p
            if (i == j) then
               partial(i,j) = 1.0_dp
            else if (q(i,i) > 0.0_dp .and. q(j,j) > 0.0_dp) then
               partial(i,j) = -q(i,j)/sqrt(q(i,i)*q(j,j))
            else
               partial(i,j) = 0.0_dp
            end if
         end do
      end do
   end subroutine partial_correlation

   pure subroutine congruence_matrix(x, y, congruence)
      real(dp), intent(in) :: x(:, :) !! First variable-by-loading matrix.
      real(dp), intent(in) :: y(:, :) !! Second variable-by-loading matrix with matching rows.
      real(dp), intent(out) :: congruence(:, :) !! Tucker congruence coefficients, rows for `x`, columns for `y`.
      real(dp) :: denom
      integer :: i, j
      do i = 1, size(x,2)
         do j = 1, size(y,2)
            denom = sqrt(sum(x(:,i)*x(:,i))*sum(y(:,j)*y(:,j)))
            if (denom > 0.0_dp) then
               congruence(i,j) = sum(x(:,i)*y(:,j))/denom
            else
               congruence(i,j) = 0.0_dp
            end if
         end do
      end do
   end subroutine congruence_matrix

   pure subroutine profile_distance_matrix(x, y, power, distance)
      real(dp), intent(in) :: x(:, :) !! First matrix, with profiles in columns.
      real(dp), intent(in) :: y(:, :) !! Second matrix, with profiles in columns.
      real(dp), intent(in), optional :: power !! Minkowski order; default is two.
      real(dp), intent(out) :: distance(:, :) !! Pairwise profile distances with shape `(ny,nx)`.
      real(dp) :: p
      integer :: i, j
      p = 2.0_dp
      if (present(power)) p = power
      do i = 1, size(x,2)
         do j = 1, size(y,2)
            distance(j,i) = sum(abs(x(:,i)-y(:,j))**p)**(1.0_dp/p)
         end do
      end do
   end subroutine profile_distance_matrix

   pure real(dp) function factor_fit_value(r, loadings, uniqueness) result(fit)
      real(dp), intent(in) :: r(:, :) !! Observed correlation matrix.
      real(dp), intent(in) :: loadings(:, :) !! Factor loading matrix.
      real(dp), intent(in), optional :: uniqueness(:) !! Optional uniquenesses; defaults to observed diagonal minus communality.
      real(dp), allocatable :: model(:, :), residual(:, :)
      integer :: i, p
      p = size(r,1)
      allocate(model(p,p),residual(p,p))
      model = matmul(loadings,transpose(loadings))
      do i = 1, p
         if (present(uniqueness)) then
            model(i,i) = model(i,i)+uniqueness(i)
         else
            model(i,i) = r(i,i)
         end if
      end do
      residual = r-model
      do i = 1, p
         residual(i,i) = 0.0_dp
      end do
      fit = 1.0_dp-sum(residual*residual)/max(tiny(1.0_dp),sum((r-diagonal_only(r))**2))
   end function factor_fit_value

   pure function diagonal_only(a) result(d)
      real(dp), intent(in) :: a(:, :) !! Square matrix whose diagonal is retained.
      real(dp) :: d(size(a,1),size(a,2))
      integer :: i
      d = 0.0_dp
      do i = 1, min(size(a,1),size(a,2))
         d(i,i) = a(i,i)
      end do
   end function diagonal_only

   pure subroutine factor_model_matrix(loadings, uniqueness, model)
      real(dp), intent(in) :: loadings(:, :) !! Factor loading matrix.
      real(dp), intent(in) :: uniqueness(:) !! Uniqueness vector with one entry per observed variable.
      real(dp), intent(out) :: model(:, :) !! Reproduced covariance/correlation matrix.
      integer :: i
      model = matmul(loadings,transpose(loadings))
      do i = 1, min(size(model,1),size(uniqueness))
         model(i,i) = model(i,i)+uniqueness(i)
      end do
   end subroutine factor_model_matrix

   subroutine bartlett_correlation_test(r, n, result)
      real(dp), intent(in) :: r(:, :) !! Correlation matrix to test against the identity.
      integer, intent(in) :: n !! Sample size.
      type(correlation_test_result), intent(out) :: result !! Bartlett chi-square statistic, degrees of freedom, and p-value.
      real(dp) :: det_r
      integer :: info, p
      p = size(r,1)
      call determinant(r,det_r,info)
      if (info /= 0 .or. det_r <= 0.0_dp) then
         result%statistic = huge(1.0_dp)
         result%df = real(p*(p-1)/2,dp)
         result%p_value = 0.0_dp
         return
      end if
      result%statistic = -log(det_r)*(real(n-1,dp)-real(2*p+5,dp)/6.0_dp)
      result%df = real(p*(p-1)/2,dp)
      result%p_value = 1.0_dp-chi_square_cdf(result%statistic,result%df)
   end subroutine bartlett_correlation_test

   pure subroutine correlation_significance(r, n, p_values)
      real(dp), intent(in) :: r(:, :) !! Correlation matrix.
      integer, intent(in) :: n !! Sample size applied to each correlation.
      real(dp), intent(out) :: p_values(:, :) !! Two-sided Student-t p-values.
      real(dp) :: t, tail
      integer :: i, j
      do j = 1, size(r,2)
         do i = 1, size(r,1)
            if (i == j) then
               p_values(i,j) = 0.0_dp
            else
               t = abs(r(i,j))*sqrt(real(max(1,n-2),dp)/max(tiny(1.0_dp),1.0_dp-r(i,j)*r(i,j)))
               tail = 1.0_dp-student_t_cdf(t,real(max(1,n-2),dp))
               p_values(i,j) = min(1.0_dp,2.0_dp*tail)
            end if
         end do
      end do
   end subroutine correlation_significance

   pure subroutine paired_correlation_test(rxy, rxz, ryz, n, result)
      real(dp), intent(in) :: rxy !! First correlation sharing one variable with `rxz`.
      real(dp), intent(in) :: rxz !! Second correlation sharing one variable with `rxy`.
      real(dp), intent(in) :: ryz !! Correlation between the two nonshared variables.
      integer, intent(in) :: n !! Common sample size.
      type(correlation_test_result), intent(out) :: result !! Williams/Steiger-style t statistic and p-value.
      real(dp) :: diff, determinant3, average, cube, denominator
      diff = rxy-rxz
      determinant3 = 1.0_dp-rxy*rxy-rxz*rxz-ryz*ryz+2.0_dp*rxy*rxz*ryz
      average = 0.5_dp*(rxy+rxz)
      cube = (1.0_dp-ryz)**3
      denominator = (2.0_dp*real(n-1,dp)/real(max(1,n-3),dp))*determinant3+average*average*cube
      if (denominator <= 0.0_dp) then
         result%statistic = 0.0_dp
      else
         result%statistic = diff*sqrt(real(n-1,dp)*(1.0_dp+ryz)/denominator)
      end if
      result%df = real(max(1,n-3),dp)
      result%p_value = min(1.0_dp,2.0_dp*(1.0_dp-student_t_cdf(abs(result%statistic),result%df)))
   end subroutine paired_correlation_test

   subroutine smooth_correlation(r, smoothed, info, epsilon_value)
      real(dp), intent(in) :: r(:, :) !! Symmetric correlation matrix, possibly not positive definite.
      real(dp), intent(out) :: smoothed(:, :) !! Positive-semidefinite correlation matrix with unit diagonal.
      integer, intent(out) :: info !! Eigen-solver status.
      real(dp), intent(in), optional :: epsilon_value !! Minimum eigenvalue, default `1e-10`.
      real(dp), allocatable :: values(:), vectors(:, :), reconstructed(:, :), scales(:)
      real(dp) :: floor_value
      integer :: i, j, p
      p = size(r,1)
      floor_value = 1.0e-10_dp
      if (present(epsilon_value)) floor_value = epsilon_value
      call symmetric_eigen(r,values,vectors,info,descending=.false.)
      if (info /= 0) then
         smoothed = r
         return
      end if
      values = max(values,floor_value)
      reconstructed = matmul(vectors*spread(values,1,p),transpose(vectors))
      allocate(scales(p))
      do i = 1, p
         scales(i) = sqrt(max(tiny(1.0_dp),reconstructed(i,i)))
      end do
      do j = 1, p
         do i = 1, p
            smoothed(i,j) = reconstructed(i,j)/(scales(i)*scales(j))
         end do
      end do
   end subroutine smooth_correlation

   pure real(dp) function rv_coefficient(x, y, xy) result(rv)
      real(dp), intent(in) :: x(:, :) !! Within-set covariance/correlation matrix for the first variable set.
      real(dp), intent(in) :: y(:, :) !! Within-set covariance/correlation matrix for the second variable set.
      real(dp), intent(in) :: xy(:, :) !! Cross-covariance/correlation matrix from first to second set.
      real(dp) :: numerator, denominator
      numerator = trace_matrix(matmul(xy,transpose(xy)))
      denominator = sqrt(trace_matrix(matmul(x,transpose(x)))*trace_matrix(matmul(y,transpose(y))))
      if (denominator > 0.0_dp) then
         rv = numerator/denominator
      else
         rv = 0.0_dp
      end if
   end function rv_coefficient

   pure real(dp) function cohen_kappa_value(table) result(kappa)
      real(dp), intent(in) :: table(:, :) !! Square agreement-count table.
      real(dp) :: total, observed, expected
      real(dp), allocatable :: rowp(:), colp(:)
      integer :: i, p
      p = min(size(table,1),size(table,2))
      total = sum(table)
      if (total <= 0.0_dp) then
         kappa = 0.0_dp
         return
      end if
      allocate(rowp(size(table,1)),colp(size(table,2)))
      rowp = sum(table,dim=2)/total
      colp = sum(table,dim=1)/total
      observed = 0.0_dp
      expected = 0.0_dp
      do i = 1, p
         observed = observed+table(i,i)/total
         expected = expected+rowp(i)*colp(i)
      end do
      if (1.0_dp-expected > tiny(1.0_dp)) then
         kappa = (observed-expected)/(1.0_dp-expected)
      else
         kappa = 1.0_dp
      end if
   end function cohen_kappa_value

   pure real(dp) function weighted_kappa_value(table, power) result(kappa)
      real(dp), intent(in) :: table(:, :) !! Square agreement-count table.
      real(dp), intent(in), optional :: power !! Distance exponent in the default weights; default two.
      real(dp), allocatable :: weights(:, :), expected(:, :), probs(:, :), rowp(:), colp(:)
      real(dp) :: pwr, total, wpo, wpc
      integer :: i, j, p
      p = min(size(table,1),size(table,2))
      pwr = 2.0_dp
      if (present(power)) pwr = power
      total = sum(table)
      if (total <= 0.0_dp .or. p <= 1) then
         kappa = 0.0_dp
         return
      end if
      allocate(weights(p,p),expected(p,p),probs(p,p),rowp(p),colp(p))
      probs = table(1:p,1:p)/total
      rowp = sum(probs,dim=2)
      colp = sum(probs,dim=1)
      do j = 1, p
         do i = 1, p
            weights(i,j) = 1.0_dp-(abs(real(i-j,dp))/real(p-1,dp))**pwr
            expected(i,j) = rowp(i)*colp(j)
         end do
      end do
      wpo = sum(weights*probs)
      wpc = sum(weights*expected)
      if (1.0_dp-wpc > tiny(1.0_dp)) then
         kappa = (wpo-wpc)/(1.0_dp-wpc)
      else
         kappa = 1.0_dp
      end if
   end function weighted_kappa_value

   pure real(dp) function pooled_sd_value(x, groups) result(sp)
      real(dp), intent(in) :: x(:) !! Numeric response values.
      integer, intent(in) :: groups(:) !! Two-group labels coded with distinct integer values.
      real(dp) :: m1, m2, s1, s2
      integer :: g1, g2, i, n1, n2
      g1 = groups(1)
      g2 = g1
      do i = 1, size(groups)
         if (groups(i) /= g1) then
            g2 = groups(i)
            exit
         end if
      end do
      m1 = 0.0_dp
      m2 = 0.0_dp
      n1 = 0
      n2 = 0
      do i = 1, size(x)
         if (groups(i) == g1 .and. .not. ieee_is_nan(x(i))) then
            m1 = m1+x(i)
            n1 = n1+1
         else if (groups(i) == g2 .and. .not. ieee_is_nan(x(i))) then
            m2 = m2+x(i)
            n2 = n2+1
         end if
      end do
      if (n1 > 0) m1 = m1/real(n1,dp)
      if (n2 > 0) m2 = m2/real(n2,dp)
      s1 = 0.0_dp
      s2 = 0.0_dp
      do i = 1, size(x)
         if (groups(i) == g1 .and. .not. ieee_is_nan(x(i))) s1 = s1+(x(i)-m1)**2
         if (groups(i) == g2 .and. .not. ieee_is_nan(x(i))) s2 = s2+(x(i)-m2)**2
      end do
      if (n1+n2 > 2) then
         sp = sqrt((s1+s2)/real(n1+n2-2,dp))
      else
         sp = 0.0_dp
      end if
   end function pooled_sd_value

   pure real(dp) function cohen_d_value(x, groups) result(d)
      real(dp), intent(in) :: x(:) !! Numeric response values.
      integer, intent(in) :: groups(:) !! Two-group labels coded with distinct integers.
      real(dp) :: mean1, mean2, sp
      integer :: g1, g2, i, n1, n2
      g1 = groups(1)
      g2 = g1
      do i = 1, size(groups)
         if (groups(i) /= g1) then
            g2 = groups(i)
            exit
         end if
      end do
      mean1 = 0.0_dp
      mean2 = 0.0_dp
      n1 = 0
      n2 = 0
      do i = 1, size(x)
         if (groups(i) == g1 .and. .not. ieee_is_nan(x(i))) then
            mean1 = mean1+x(i)
            n1 = n1+1
         else if (groups(i) == g2 .and. .not. ieee_is_nan(x(i))) then
            mean2 = mean2+x(i)
            n2 = n2+1
         end if
      end do
      if (n1 > 0) mean1 = mean1/real(n1,dp)
      if (n2 > 0) mean2 = mean2/real(n2,dp)
      sp = pooled_sd_value(x,groups)
      if (sp > 0.0_dp) then
         d = (mean2-mean1)/sp
      else
         d = 0.0_dp
      end if
   end function cohen_d_value

   subroutine mahalanobis_squared(x, distances, info)
      real(dp), intent(in) :: x(:, :) !! Observation-by-variable data matrix without missing rows.
      real(dp), intent(out) :: distances(:) !! Squared Mahalanobis distance for each row.
      integer, intent(out) :: info !! Matrix-inversion status.
      real(dp), allocatable :: cov(:, :), covinv(:, :), means(:), d(:)
      integer :: i, p
      p = size(x,2)
      allocate(cov(p,p),means(p),d(p))
      call covariance_matrix(x,cov)
      call pseudoinverse(cov,covinv,info=info)
      if (info /= 0) then
         distances = 0.0_dp
         return
      end if
      do i = 1, p
         means(i) = finite_mean(x(:,i))
      end do
      do i = 1, size(x,1)
         d = x(i,:)-means
         distances(i) = dot_product(d,matmul(covinv,d))
      end do
   end subroutine mahalanobis_squared

   subroutine mardia_statistics(x, skewness_value, kurtosis_value, info)
      real(dp), intent(in) :: x(:, :) !! Complete observation-by-variable matrix.
      real(dp), intent(out) :: skewness_value !! Mardia multivariate skewness coefficient `b1p`.
      real(dp), intent(out) :: kurtosis_value !! Mardia multivariate kurtosis coefficient `b2p`.
      integer, intent(out) :: info !! Matrix-inversion status.
      real(dp), allocatable :: cov(:, :), covinv(:, :), centered(:, :), means(:), d(:, :)
      integer :: i, n, p
      n = size(x,1)
      p = size(x,2)
      allocate(cov(p,p),centered(n,p),means(p))
      do i = 1, p
         means(i) = finite_mean(x(:,i))
         centered(:,i) = x(:,i)-means(i)
      end do
      call covariance_matrix(x,cov)
      call pseudoinverse(cov,covinv,info=info)
      if (info /= 0) then
         skewness_value = 0.0_dp
         kurtosis_value = 0.0_dp
         return
      end if
      d = matmul(matmul(centered,covinv),transpose(centered))
      skewness_value = sum(d**3)/real(n*n,dp)
      kurtosis_value = trace_matrix(matmul(d,d))/real(n,dp)
   end subroutine mardia_statistics

   pure real(dp) function auto_correlation(x, lag) result(r)
      real(dp), intent(in) :: x(:) !! Ordered series.
      integer, intent(in), optional :: lag !! Positive lag; default one.
      integer :: k
      k = 1
      if (present(lag)) k = max(1,lag)
      if (size(x) <= k) then
         r = 0.0_dp
      else
         r = finite_correlation(x(1:size(x)-k),x(1+k:size(x)))
      end if
   end function auto_correlation

   pure real(dp) function point_biserial_correlation(x, group) result(r)
      real(dp), intent(in) :: x(:) !! Continuous variable.
      integer, intent(in) :: group(:) !! Binary group indicator using two distinct integer codes.
      real(dp), allocatable :: g(:)
      integer :: i, base
      allocate(g(size(group)))
      base = group(1)
      do i = 1, size(group)
         g(i) = merge(0.0_dp,1.0_dp,group(i)==base)
      end do
      r = finite_correlation(x,g)
   end function point_biserial_correlation

   pure real(dp) function biserial_correlation(x, group) result(rbis)
      real(dp), intent(in) :: x(:) !! Continuous variable.
      integer, intent(in) :: group(:) !! Binary group indicator using two distinct integer codes.
      real(dp) :: rpb, p, q, threshold, density
      integer :: i, base, n0
      base = group(1)
      n0 = 0
      do i = 1, size(group)
         if (group(i) == base) n0 = n0+1
      end do
      p = real(n0,dp)/real(size(group),dp)
      q = 1.0_dp-p
      rpb = point_biserial_correlation(x,group)
      if (p <= 0.0_dp .or. q <= 0.0_dp) then
         rbis = 0.0_dp
         return
      end if
      threshold = inverse_normal_local(p)
      density = exp(-0.5_dp*threshold*threshold)/sqrt(2.0_dp*acos(-1.0_dp))
      if (density > tiny(1.0_dp)) then
         rbis = rpb*sqrt(p*q)/density
      else
         rbis = rpb
      end if
   end function biserial_correlation

   pure real(dp) function inverse_normal_local(p) result(x)
      real(dp), intent(in) :: p !! Probability in `(0,1)` used by the biserial correction.
      real(dp) :: lo, hi, mid
      integer :: iter
      lo = -8.0_dp
      hi = 8.0_dp
      do iter = 1, 100
         mid = 0.5_dp*(lo+hi)
         if (normal_cdf(mid) < p) then
            lo = mid
         else
            hi = mid
         end if
      end do
      x = 0.5_dp*(lo+hi)
   end function inverse_normal_local

end module psych_correlations
