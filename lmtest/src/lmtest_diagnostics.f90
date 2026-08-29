module lmtest_diagnostics
   use lmtest_kinds, only : dp
   use lmtest_types, only : test_result, bg_test_result, lm_result
   use lmtest_regression, only : lm_fit, lm_wfit, recursive_residuals
   use lmtest_distributions, only : chi_square_sf, f_cdf, f_sf, &
      student_t_sf, normal_cdf, normal_sf
   use lmtest_linalg, only : invert_spd, symmetric_eigenvalues, &
      sort_order, reorder_rows, reorder_vector, covariance_matrix
   use lmtest_pan, only : pan_probability
   implicit none
   private

   public :: breusch_godfrey_test
   public :: breusch_pagan_test
   public :: durbin_watson_test
   public :: goldfeld_quandt_test
   public :: harvey_collier_test
   public :: harrison_mccabe_test
   public :: rainbow_test
   public :: reset_test

contains

   function breusch_godfrey_test(xin, yin, order, use_f, fill, &
      order_by) result(out)
      real(dp), intent(in) :: xin(:,:), yin(:)
      integer, intent(in), optional :: order
      logical, intent(in), optional :: use_f
      real(dp), intent(in), optional :: fill
      real(dp), intent(in), optional :: order_by(:)
      type(bg_test_result) :: out

      type(lm_result) :: base, aux
      real(dp), allocatable :: x(:,:), y(:), z(:,:), xz(:,:)
      integer, allocatable :: idx(:)
      real(dp) :: fill_value
      integer :: n, k, m, lag, i
      logical :: as_f

      x = xin
      y = yin
      if (present(order_by)) then
         idx = sort_order(order_by)
         x = reorder_rows(xin, idx)
         y = reorder_vector(yin, idx)
      end if

      n = size(x, 1)
      k = size(x, 2)
      m = 1
      if (present(order)) m = order
      fill_value = 0.0_dp
      if (present(fill)) fill_value = fill
      as_f = .false.
      if (present(use_f)) as_f = use_f

      base = lm_fit(x, y)
      if (base%info /= 0 .or. m < 1 .or. n <= k + m) return

      allocate(z(n, m), xz(n, k + m))
      z = fill_value
      do lag = 1, m
         do i = lag + 1, n
            z(i, lag) = base%residuals(i - lag)
         end do
      end do
      xz(:, 1:k) = x
      xz(:, k + 1:k + m) = z

      aux = lm_fit(xz, base%residuals)
      if (aux%info /= 0) return

      allocate(out%coefficients(k + m), out%vcov(k + m, k + m))
      out%coefficients = aux%beta
      out%vcov = aux%vcov

      if (as_f) then
         out%test%statistic = &
            ((base%rss - aux%rss) / real(m, dp)) / &
            (aux%rss / real(n - k - m, dp))
         out%test%df1 = real(m, dp)
         out%test%df2 = real(n - k - m, dp)
         out%test%p_value = f_sf(out%test%statistic, &
            out%test%df1, out%test%df2)
      else
         out%test%statistic = real(n, dp) * sum(aux%fitted**2) / &
            base%rss
         out%test%df1 = real(m, dp)
         out%test%p_value = chi_square_sf(out%test%statistic, &
            out%test%df1)
      end if
   end function breusch_godfrey_test

   function breusch_pagan_test(x, y, z, studentize, weights) result(out)
      real(dp), intent(in) :: x(:,:), y(:), z(:,:)
      logical, intent(in), optional :: studentize
      real(dp), intent(in), optional :: weights(:)
      type(test_result) :: out

      type(lm_result) :: base, aux
      real(dp), allocatable :: wts(:), w(:)
      real(dp) :: sigma2
      integer :: n, n_positive
      logical :: stud

      n = size(x, 1)
      allocate(wts(n))
      wts = 1.0_dp
      if (present(weights)) wts = weights
      stud = .true.
      if (present(studentize)) stud = studentize
      if (size(z, 1) /= n .or. size(z, 2) < 2) return

      base = lm_wfit(x, y, wts)
      if (base%info /= 0) return
      n_positive = count(wts > 0.0_dp)
      sigma2 = sum(wts * base%residuals**2) / real(n_positive, dp)
      allocate(w(n))

      if (stud) then
         w = base%residuals**2 - sigma2
         aux = lm_wfit(z, w, wts)
         out%statistic = real(n_positive, dp) * &
            sum(wts * aux%fitted**2) / sum(wts * w**2)
      else
         w = base%residuals**2 / sigma2 - 1.0_dp
         aux = lm_wfit(z, w, wts)
         out%statistic = 0.5_dp * sum(wts * aux%fitted**2)
      end if

      out%df1 = real(aux%rank - 1, dp)
      out%p_value = chi_square_sf(out%statistic, out%df1)
   end function breusch_pagan_test

   function durbin_watson_test(xin, yin, alternative, exact, &
      iterations, tol, order_by) result(out)
      real(dp), intent(in) :: xin(:,:), yin(:)
      character(len=*), intent(in), optional :: alternative
      logical, intent(in), optional :: exact
      integer, intent(in), optional :: iterations
      real(dp), intent(in), optional :: tol
      real(dp), intent(in), optional :: order_by(:)
      type(test_result) :: out

      type(lm_result) :: fit
      real(dp), allocatable :: x(:,:), y(:)
      real(dp), allocatable :: xtx(:,:), q1(:,:), mmat(:,:)
      real(dp), allocatable :: amat(:,:), bmat(:,:), all_eigen(:), ev(:)
      real(dp), allocatable :: ax(:,:), xaxq(:,:), temp(:,:)
      integer, allocatable :: idx(:)
      real(dp) :: tolerance, pdw, p, qval, dmean, dvar, zval
      integer :: n, k, i, j, info, nev, nterms
      logical :: do_exact
      character(len=16) :: alt

      x = xin
      y = yin
      if (present(order_by)) then
         idx = sort_order(order_by)
         x = reorder_rows(xin, idx)
         y = reorder_vector(yin, idx)
      end if

      n = size(x, 1)
      k = size(x, 2)
      tolerance = 1.0e-10_dp
      if (present(tol)) tolerance = tol
      nterms = 15
      if (present(iterations)) nterms = iterations
      do_exact = n < 100
      if (present(exact)) do_exact = exact
      alt = 'greater'
      if (present(alternative)) alt = adjustl(alternative)

      fit = lm_fit(x, y)
      if (fit%info /= 0 .or. fit%rss <= 0.0_dp) return
      out%statistic = sum((fit%residuals(2:n) - &
         fit%residuals(1:n - 1))**2) / fit%rss
      if (n < 3) then
         out%p_value = 1.0_dp
         return
      end if

      xtx = matmul(transpose(x), x)
      call invert_spd(xtx, q1, info)
      if (info /= 0) return

      if (do_exact) then
         allocate(mmat(n, n), amat(n, n), bmat(n, n))
         mmat = 0.0_dp
         do i = 1, n
            mmat(i, i) = 1.0_dp
         end do
         mmat = mmat - matmul(x, matmul(q1, transpose(x)))

         amat = 0.0_dp
         do i = 1, n
            if (i == 1 .or. i == n) then
               amat(i, i) = 1.0_dp
            else
               amat(i, i) = 2.0_dp
            end if
         end do
         do i = 1, n - 1
            amat(i, i + 1) = -1.0_dp
            amat(i + 1, i) = -1.0_dp
         end do

         bmat = matmul(mmat, matmul(amat, mmat))
         bmat = 0.5_dp * (bmat + transpose(bmat))
         call symmetric_eigenvalues(bmat, all_eigen, info)
         if (info == 0) then
            nev = count(all_eigen > tolerance)
            allocate(ev(nev))
            j = 0
            do i = 1, size(all_eigen)
               if (all_eigen(i) > tolerance) then
                  j = j + 1
                  ev(j) = all_eigen(i)
               end if
            end do
            if (nev > 0) then
               pdw = pan_probability(ev, out%statistic, 0.0_dp, nterms)
               select case (trim(alt))
               case ('two.sided')
                  out%p_value = 2.0_dp * min(pdw, 1.0_dp - pdw)
               case ('less')
                  out%p_value = 1.0_dp - pdw
               case default
                  out%p_value = pdw
               end select
               if (out%p_value >= 0.0_dp .and. &
                  out%p_value <= 1.0_dp) return
            end if
         end if
      end if

      if (n < max(5, k)) then
         out%p_value = 1.0_dp
         return
      end if

      allocate(ax(n, k))
      do j = 1, k
         ax(1, j) = x(1, j) - x(2, j)
         do i = 2, n - 1
            ax(i, j) = 2.0_dp * x(i, j) - x(i - 1, j) - x(i + 1, j)
         end do
         ax(n, j) = x(n, j) - x(n - 1, j)
      end do

      xaxq = matmul(transpose(x), matmul(ax, q1))
      p = 2.0_dp * real(n - 1, dp) - trace_matrix(xaxq)
      temp = matmul(matmul(transpose(ax), ax), q1)
      qval = 2.0_dp * real(3 * n - 4, dp) - &
         2.0_dp * trace_matrix(temp) + &
         trace_matrix(matmul(xaxq, xaxq))
      dmean = p / real(n - k, dp)
      dvar = 2.0_dp / (real(n - k, dp) * real(n - k + 2, dp)) * &
         (qval - p * dmean)
      if (dvar <= 0.0_dp) then
         out%p_value = 1.0_dp
         return
      end if

      zval = (out%statistic - dmean) / sqrt(dvar)
      select case (trim(alt))
      case ('two.sided')
         out%p_value = 2.0_dp * normal_sf(abs(zval))
      case ('less')
         out%p_value = normal_sf(zval)
      case default
         out%p_value = normal_cdf(zval)
      end select
   end function durbin_watson_test

   function goldfeld_quandt_test(xin, yin, point, fraction, &
      alternative, order_by) result(out)
      real(dp), intent(in) :: xin(:,:), yin(:)
      real(dp), intent(in), optional :: point, fraction
      character(len=*), intent(in), optional :: alternative
      real(dp), intent(in), optional :: order_by(:)
      type(test_result) :: out

      real(dp), allocatable :: x(:,:), y(:)
      integer, allocatable :: idx(:)
      type(lm_result) :: fit1, fit2
      real(dp) :: split, frac, mss1, mss2, lower, upper
      integer :: n, k, point1, point2
      character(len=16) :: alt

      x = xin
      y = yin
      n = size(x, 1)
      k = size(x, 2)
      split = 0.5_dp
      if (present(point)) split = point
      frac = 0.0_dp
      if (present(fraction)) frac = fraction
      alt = 'greater'
      if (present(alternative)) alt = adjustl(alternative)

      if (present(order_by)) then
         idx = sort_order(order_by)
         x = reorder_rows(xin, idx)
         y = reorder_vector(yin, idx)
      end if

      if (split > 1.0_dp) then
         if (frac < 1.0_dp) frac = floor(frac * real(n, dp))
         point1 = int(split - ceiling(frac / 2.0_dp))
         point2 = int(split + ceiling(frac / 2.0_dp + 0.01_dp))
      else
         if (frac >= 1.0_dp) frac = frac / real(n, dp)
         point1 = int(floor((split - frac / 2.0_dp) * real(n, dp)))
         point2 = int(ceiling((split + frac / 2.0_dp) * &
            real(n, dp) + 0.01_dp))
      end if
      if (point2 > n - k + 1 .or. point1 < k) return

      fit1 = lm_fit(x(1:point1, :), y(1:point1))
      fit2 = lm_fit(x(point2:n, :), y(point2:n))
      mss1 = fit1%rss / real(point1 - k, dp)
      mss2 = fit2%rss / real(n - point2 + 1 - k, dp)

      out%statistic = mss2 / mss1
      out%df1 = real(n - point2 + 1 - k, dp)
      out%df2 = real(point1 - k, dp)
      lower = f_cdf(out%statistic, out%df1, out%df2)
      upper = f_sf(out%statistic, out%df1, out%df2)
      select case (trim(alt))
      case ('two.sided')
         out%p_value = 2.0_dp * min(lower, upper)
      case ('less')
         out%p_value = lower
      case default
         out%p_value = upper
      end select
   end function goldfeld_quandt_test

   function harvey_collier_test(xin, yin, order_by) result(out)
      real(dp), intent(in) :: xin(:,:), yin(:)
      real(dp), intent(in), optional :: order_by(:)
      type(test_result) :: out

      real(dp), allocatable :: x(:,:), y(:), w(:), scaled(:)
      integer, allocatable :: idx(:)
      integer :: n, k, info, m
      real(dp) :: meanw, varw, sigma, var_scaled

      x = xin
      y = yin
      n = size(x, 1)
      k = size(x, 2)
      if (present(order_by)) then
         idx = sort_order(order_by)
         x = reorder_rows(xin, idx)
         y = reorder_vector(yin, idx)
      end if

      call recursive_residuals(x, y, w, info)
      if (info /= 0) return
      m = size(w)
      meanw = sum(w) / real(m, dp)
      varw = sum((w - meanw)**2) / real(max(1, m - 1), dp)
      sigma = sqrt(varw * real(m - 1, dp) / real(n - k - 1, dp))
      scaled = w / sigma
      meanw = sum(scaled) / real(m, dp)
      var_scaled = sum((scaled - meanw)**2) / real(max(1, m - 1), dp)

      out%statistic = abs(sum(scaled) / sqrt(real(n - k, dp))) / &
         sqrt(var_scaled)
      out%df1 = real(n - k - 1, dp)
      out%p_value = 2.0_dp * student_t_sf(out%statistic, out%df1)
   end function harvey_collier_test

   function harrison_mccabe_test(xin, yin, point, order_by, &
      nsim, seed) result(out)
      real(dp), intent(in) :: xin(:,:), yin(:)
      real(dp), intent(in), optional :: point, order_by(:)
      integer, intent(in), optional :: nsim, seed
      type(test_result) :: out

      real(dp), allocatable :: x(:,:), y(:), z(:)
      integer, allocatable :: idx(:), seed_values(:)
      type(lm_result) :: fit
      real(dp) :: split, u1, u2, meanz, stat
      integer :: n, k, breakpoint, nsim_use
      integer :: i, j, nseed, count_less

      x = xin
      y = yin
      n = size(x, 1)
      k = size(x, 2)
      split = 0.5_dp
      if (present(point)) split = point
      if (split < 1.0_dp) then
         breakpoint = int(floor(split * real(n, dp)))
      else
         breakpoint = int(split)
      end if
      if (breakpoint > n - k .or. breakpoint < k) return

      if (present(order_by)) then
         idx = sort_order(order_by)
         x = reorder_rows(xin, idx)
         y = reorder_vector(yin, idx)
      end if

      fit = lm_fit(x, y)
      if (fit%info /= 0) return
      out%statistic = sum(fit%residuals(1:breakpoint)**2) / fit%rss

      nsim_use = 1000
      if (present(nsim)) nsim_use = nsim
      if (nsim_use <= 0) then
         out%p_value = -1.0_dp
         return
      end if

      if (present(seed)) then
         call random_seed(size=nseed)
         allocate(seed_values(nseed))
         do i = 1, nseed
            seed_values(i) = seed + 37 * (i - 1)
         end do
         call random_seed(put=seed_values)
      end if

      allocate(z(n))
      count_less = 0
      do i = 1, nsim_use
         j = 1
         do while (j <= n)
            call random_number(u1)
            call random_number(u2)
            u1 = max(u1, tiny(1.0_dp))
            z(j) = sqrt(-2.0_dp * log(u1)) * &
               cos(2.0_dp * acos(-1.0_dp) * u2)
            if (j + 1 <= n) then
               z(j + 1) = sqrt(-2.0_dp * log(u1)) * &
                  sin(2.0_dp * acos(-1.0_dp) * u2)
            end if
            j = j + 2
         end do
         meanz = sum(z) / real(n, dp)
         z = z - meanz
         stat = sum(z(1:breakpoint)**2) / sum(z**2)
         if (stat <= out%statistic) count_less = count_less + 1
      end do
      out%p_value = real(count_less, dp) / real(nsim_use, dp)
   end function harrison_mccabe_test

   function rainbow_test(xin, yin, fraction, center, order_by, &
      mahalanobis, mahalanobis_center) result(out)
      real(dp), intent(in) :: xin(:,:), yin(:)
      real(dp), intent(in), optional :: fraction, center, order_by(:)
      logical, intent(in), optional :: mahalanobis
      real(dp), intent(in), optional :: mahalanobis_center(:)
      type(test_result) :: out

      real(dp), allocatable :: x(:,:), y(:), distances(:)
      real(dp), allocatable :: xtx(:,:), inv(:,:), ctr(:), delta(:)
      integer, allocatable :: idx(:)
      type(lm_result) :: allfit, subfit
      real(dp) :: frac, center_use, quant
      integer :: n, k, first, last, n1, i, info
      logical :: use_mahalanobis

      x = xin
      y = yin
      n = size(x, 1)
      k = size(x, 2)
      frac = 0.5_dp
      if (present(fraction)) frac = fraction
      use_mahalanobis = .false.
      if (present(mahalanobis)) use_mahalanobis = mahalanobis

      if (use_mahalanobis) then
         allocate(ctr(k))
         if (present(mahalanobis_center)) then
            if (size(mahalanobis_center) /= k) return
            ctr = mahalanobis_center
         else
            ctr = sum(x, dim=1) / real(n, dp)
         end if
         xtx = matmul(transpose(x), x)
         call invert_spd(xtx, inv, info)
         if (info /= 0) return
         allocate(distances(n), delta(k))
         do i = 1, n
            delta = x(i, :) - ctr
            distances(i) = dot_product(delta, matmul(inv, delta))
         end do
         idx = sort_order(distances)
         x = reorder_rows(x, idx)
         y = reorder_vector(y, idx)
         first = 1
         last = int(floor(frac * real(n, dp)))
      else
         if (present(order_by)) then
            idx = sort_order(order_by)
            x = reorder_rows(x, idx)
            y = reorder_vector(y, idx)
         end if
         center_use = 0.5_dp
         if (present(center)) center_use = center
         if (center_use > 1.0_dp) center_use = center_use / real(n, dp)
         quant = 1.0_dp + (real(n, dp) - 1.0_dp) * &
            (center_use - frac / 2.0_dp)
         first = int(ceiling(quant))
         last = first + int(floor(frac * real(n, dp))) - 1
      end if

      n1 = last - first + 1
      if (n1 < k .or. first < 1 .or. last > n) return
      allfit = lm_fit(x, y)
      subfit = lm_fit(x(first:last, :), y(first:last))
      out%statistic = ((allfit%rss - subfit%rss) / real(n - n1, dp)) / &
         (subfit%rss / real(n1 - k, dp))
      out%df1 = real(n - n1, dp)
      out%df2 = real(n1 - k, dp)
      out%p_value = f_sf(out%statistic, out%df1, out%df2)
   end function rainbow_test

   function reset_test(x, y, powers, kind) result(out)
      real(dp), intent(in) :: x(:,:), y(:)
      integer, intent(in), optional :: powers(:)
      character(len=*), intent(in), optional :: kind
      type(test_result) :: out

      integer, allocatable :: power_use(:)
      real(dp), allocatable :: z(:,:), xz(:,:), basevars(:,:)
      real(dp), allocatable :: pc(:), eigenvector(:), cov(:,:)
      type(lm_result) :: fit1, fit2
      character(len=16) :: kind_use
      integer :: n, k, q, i, j, nv, col
      real(dp) :: normv

      n = size(x, 1)
      k = size(x, 2)
      kind_use = 'fitted'
      if (present(kind)) kind_use = adjustl(kind)
      if (present(powers)) then
         allocate(power_use(size(powers)))
         power_use = powers
      else
         allocate(power_use(2))
         power_use = [2, 3]
      end if

      select case (trim(kind_use))
      case ('regressor')
         nv = 0
         do j = 1, k
            if (maxval(abs(x(:, j) - x(1, j))) > 1.0e-12_dp) nv = nv + 1
         end do
         if (nv < 1) return
         allocate(basevars(n, nv))
         col = 0
         do j = 1, k
            if (maxval(abs(x(:, j) - x(1, j))) > 1.0e-12_dp) then
               col = col + 1
               basevars(:, col) = x(:, j)
            end if
         end do
         allocate(z(n, nv * size(power_use)))
         col = 0
         do j = 1, nv
            do i = 1, size(power_use)
               col = col + 1
               z(:, col) = basevars(:, j)**power_use(i)
            end do
         end do

      case ('princomp')
         nv = 0
         do j = 1, k
            if (maxval(abs(x(:, j) - x(1, j))) > 1.0e-12_dp) nv = nv + 1
         end do
         if (nv < 1) return
         allocate(basevars(n, nv))
         col = 0
         do j = 1, k
            if (maxval(abs(x(:, j) - x(1, j))) > 1.0e-12_dp) then
               col = col + 1
               basevars(:, col) = x(:, j)
            end if
         end do
         cov = covariance_matrix(basevars)
         allocate(eigenvector(nv))
         eigenvector = 1.0_dp / sqrt(real(nv, dp))
         do i = 1, 100
            eigenvector = matmul(cov, eigenvector)
            normv = sqrt(sum(eigenvector**2))
            if (normv > 0.0_dp) eigenvector = eigenvector / normv
         end do
         pc = matmul(basevars, eigenvector)
         allocate(z(n, size(power_use)))
         do i = 1, size(power_use)
            z(:, i) = pc**power_use(i)
         end do

      case default
         fit1 = lm_fit(x, y)
         if (fit1%info /= 0) return
         allocate(z(n, size(power_use)))
         do i = 1, size(power_use)
            z(:, i) = fit1%fitted**power_use(i)
         end do
      end select

      q = size(z, 2)
      allocate(xz(n, k + q))
      xz(:, 1:k) = x
      xz(:, k + 1:k + q) = z
      fit1 = lm_fit(x, y)
      fit2 = lm_fit(xz, y)
      if (fit1%info /= 0 .or. fit2%info /= 0) return

      out%df1 = real(q, dp)
      out%df2 = real(n - k - q, dp)
      out%statistic = (out%df2 / out%df1) * &
         ((fit1%rss - fit2%rss) / fit2%rss)
      out%p_value = f_sf(out%statistic, out%df1, out%df2)
   end function reset_test

   pure real(dp) function trace_matrix(a) result(v)
      real(dp), intent(in) :: a(:,:)
      integer :: i
      v = 0.0_dp
      do i = 1, min(size(a, 1), size(a, 2))
         v = v + a(i, i)
      end do
   end function trace_matrix

end module lmtest_diagnostics
