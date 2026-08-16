module rfast_tests
   use rfast_special, only : dp, student_t_cdf, f_cdf, chisq_cdf, normal_cdf, nan_r
   use rfast_arrays, only : mean_r, variance_r, rank_average
   implicit none
   private

   type, public :: test_result
      real(dp) :: statistic = 0.0_dp
      real(dp) :: pvalue = 1.0_dp
      real(dp) :: df1 = 0.0_dp
      real(dp) :: df2 = 0.0_dp
      real(dp) :: estimate = 0.0_dp
   end type test_result

   public :: ttest1, ttest2, paired_ttest, ftest_variance, vartest_chisq
   public :: chi2_test, g2_test, mcnemar_test, proportion_test
   public :: odds_ratio_2x2, relative_risk_2x2, auc_score, kruskal_test
   public :: poisson_dispersion_test

contains

   function ttest1(x, mu, alternative) result(res)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in), optional :: mu
      integer, intent(in), optional :: alternative
      type(test_result) :: res
      real(dp) :: m0, m, se, c
      integer :: alt, n

      m0 = 0.0_dp
      if (present(mu)) m0 = mu
      alt = 0
      if (present(alternative)) alt = alternative
      n = size(x)
      m = mean_r(x)
      se = sqrt(variance_r(x) / real(n, dp))
      res%statistic = (m - m0) / se
      res%df1 = real(n - 1, dp)
      res%estimate = m
      c = student_t_cdf(res%statistic, res%df1)
      select case (alt)
      case (-1)
         res%pvalue = c
      case (1)
         res%pvalue = 1.0_dp - c
      case default
         res%pvalue = 2.0_dp * min(c, 1.0_dp - c)
      end select
   end function ttest1

   function ttest2(x, y, equal_var, alternative) result(res)
      real(dp), intent(in) :: x(:), y(:)
      logical, intent(in), optional :: equal_var
      integer, intent(in), optional :: alternative
      type(test_result) :: res
      logical :: eq
      integer :: nx, ny, alt
      real(dp) :: mx, my, vx, vy, se, df, c, sp

      eq = .false.
      if (present(equal_var)) eq = equal_var
      alt = 0
      if (present(alternative)) alt = alternative
      nx = size(x)
      ny = size(y)
      mx = mean_r(x)
      my = mean_r(y)
      vx = variance_r(x)
      vy = variance_r(y)
      if (eq) then
         sp = (real(nx - 1, dp) * vx + real(ny - 1, dp) * vy) / real(nx + ny - 2, dp)
         se = sqrt(sp * (1.0_dp / real(nx, dp) + 1.0_dp / real(ny, dp)))
         df = real(nx + ny - 2, dp)
      else
         se = sqrt(vx / real(nx, dp) + vy / real(ny, dp))
         df = (vx / real(nx, dp) + vy / real(ny, dp))**2
         df = df / ((vx / real(nx, dp))**2 / real(nx - 1, dp) + &
                    (vy / real(ny, dp))**2 / real(ny - 1, dp))
      end if
      res%statistic = (mx - my) / se
      res%df1 = df
      res%estimate = mx - my
      c = student_t_cdf(res%statistic, df)
      select case (alt)
      case (-1)
         res%pvalue = c
      case (1)
         res%pvalue = 1.0_dp - c
      case default
         res%pvalue = 2.0_dp * min(c, 1.0_dp - c)
      end select
   end function ttest2

   function paired_ttest(x, y, alternative) result(res)
      real(dp), intent(in) :: x(:), y(:)
      integer, intent(in), optional :: alternative
      type(test_result) :: res
      if (present(alternative)) then
         res = ttest1(x - y, 0.0_dp, alternative)
      else
         res = ttest1(x - y)
      end if
   end function paired_ttest

   function ftest_variance(x, y, alternative) result(res)
      real(dp), intent(in) :: x(:), y(:)
      integer, intent(in), optional :: alternative
      type(test_result) :: res
      integer :: alt
      real(dp) :: c

      alt = 0
      if (present(alternative)) alt = alternative
      res%statistic = variance_r(x) / variance_r(y)
      res%df1 = real(size(x) - 1, dp)
      res%df2 = real(size(y) - 1, dp)
      c = f_cdf(res%statistic, res%df1, res%df2)
      select case (alt)
      case (-1)
         res%pvalue = c
      case (1)
         res%pvalue = 1.0_dp - c
      case default
         res%pvalue = 2.0_dp * min(c, 1.0_dp - c)
      end select
   end function ftest_variance

   function vartest_chisq(x, var0, alternative) result(res)
      real(dp), intent(in) :: x(:), var0
      integer, intent(in), optional :: alternative
      type(test_result) :: res
      integer :: alt
      real(dp) :: c

      alt = 0
      if (present(alternative)) alt = alternative
      res%df1 = real(size(x) - 1, dp)
      res%statistic = res%df1 * variance_r(x) / var0
      c = chisq_cdf(res%statistic, res%df1)
      select case (alt)
      case (-1)
         res%pvalue = c
      case (1)
         res%pvalue = 1.0_dp - c
      case default
         res%pvalue = 2.0_dp * min(c, 1.0_dp - c)
      end select
   end function vartest_chisq

   function chi2_test(table, yates) result(res)
      real(dp), intent(in) :: table(:,:)
      logical, intent(in), optional :: yates
      type(test_result) :: res
      real(dp) :: rt(size(table,1)), ct(size(table,2)), tot, e, delta
      integer :: i, j
      logical :: yc

      yc = .false.
      if (present(yates)) yc = yates
      rt = sum(table, dim=2)
      ct = sum(table, dim=1)
      tot = sum(table)
      res%statistic = 0.0_dp
      do j = 1, size(table,2)
         do i = 1, size(table,1)
            e = rt(i) * ct(j) / tot
            if (e <= 0.0_dp) cycle
            delta = abs(table(i,j) - e)
            if (yc .and. size(table,1) == 2 .and. size(table,2) == 2) then
               delta = max(0.0_dp, delta - 0.5_dp)
            end if
            res%statistic = res%statistic + delta * delta / e
         end do
      end do
      res%df1 = real((size(table,1) - 1) * (size(table,2) - 1), dp)
      res%pvalue = 1.0_dp - chisq_cdf(res%statistic, res%df1)
   end function chi2_test

   function g2_test(table) result(res)
      real(dp), intent(in) :: table(:,:)
      type(test_result) :: res
      real(dp) :: rt(size(table,1)), ct(size(table,2)), tot, e
      integer :: i, j

      rt = sum(table, dim=2)
      ct = sum(table, dim=1)
      tot = sum(table)
      res%statistic = 0.0_dp
      do j = 1, size(table,2)
         do i = 1, size(table,1)
            e = rt(i) * ct(j) / tot
            if (table(i,j) > 0.0_dp .and. e > 0.0_dp) then
               res%statistic = res%statistic + 2.0_dp * table(i,j) * log(table(i,j) / e)
            end if
         end do
      end do
      res%df1 = real((size(table,1) - 1) * (size(table,2) - 1), dp)
      res%pvalue = 1.0_dp - chisq_cdf(res%statistic, res%df1)
   end function g2_test

   function mcnemar_test(table, correct) result(res)
      real(dp), intent(in) :: table(2,2)
      logical, intent(in), optional :: correct
      type(test_result) :: res
      real(dp) :: b, c, delta
      logical :: cor

      cor = .true.
      if (present(correct)) cor = correct
      b = table(1,2)
      c = table(2,1)
      delta = abs(b - c)
      if (cor) delta = max(0.0_dp, delta - 1.0_dp)
      if (b + c > 0.0_dp) res%statistic = delta * delta / (b + c)
      res%df1 = 1.0_dp
      res%pvalue = 1.0_dp - chisq_cdf(res%statistic, 1.0_dp)
   end function mcnemar_test

   function proportion_test(success, n, p0, alternative, correct) result(res)
      integer, intent(in) :: success, n
      real(dp), intent(in), optional :: p0
      integer, intent(in), optional :: alternative
      logical, intent(in), optional :: correct
      type(test_result) :: res
      real(dp) :: p, pp, z, c, adj
      integer :: alt
      logical :: cor

      pp = 0.5_dp
      if (present(p0)) pp = p0
      alt = 0
      if (present(alternative)) alt = alternative
      cor = .false.
      if (present(correct)) cor = correct
      p = real(success, dp) / real(n, dp)
      adj = 0.0_dp
      if (cor) adj = 0.5_dp / real(n, dp)
      z = (abs(p - pp) - adj) / sqrt(pp * (1.0_dp - pp) / real(n, dp))
      if (p < pp) z = -z
      res%statistic = z
      res%estimate = p
      c = normal_cdf(z)
      select case (alt)
      case (-1)
         res%pvalue = c
      case (1)
         res%pvalue = 1.0_dp - c
      case default
         res%pvalue = 2.0_dp * min(c, 1.0_dp - c)
      end select
   end function proportion_test

   pure real(dp) function odds_ratio_2x2(table, correction) result(oratio)
      real(dp), intent(in) :: table(2,2)
      real(dp), intent(in), optional :: correction
      real(dp) :: c
      c = 0.0_dp
      if (present(correction)) c = correction
      oratio = ((table(1,1) + c) * (table(2,2) + c)) / &
               ((table(1,2) + c) * (table(2,1) + c))
   end function odds_ratio_2x2

   pure real(dp) function relative_risk_2x2(table, correction) result(rr)
      real(dp), intent(in) :: table(2,2)
      real(dp), intent(in), optional :: correction
      real(dp) :: c, r1, r2
      c = 0.0_dp
      if (present(correction)) c = correction
      r1 = (table(1,1) + c) / (sum(table(1,:)) + 2.0_dp * c)
      r2 = (table(2,1) + c) / (sum(table(2,:)) + 2.0_dp * c)
      rr = r1 / r2
   end function relative_risk_2x2

   function auc_score(score, label) result(auc)
      real(dp), intent(in) :: score(:)
      integer, intent(in) :: label(:)
      real(dp) :: auc
      real(dp), allocatable :: r(:)
      integer :: n1, n0

      r = rank_average(score)
      n1 = count(label /= 0)
      n0 = size(label) - n1
      if (n1 == 0 .or. n0 == 0) then
         auc = nan_r()
         return
      end if
      auc = (sum(r, mask=label /= 0) - real(n1 * (n1 + 1), dp) / 2.0_dp) / real(n1 * n0, dp)
   end function auc_score

   function kruskal_test(x, group) result(res)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: group(:)
      type(test_result) :: res
      real(dp), allocatable :: r(:)
      integer, allocatable :: ug(:)
      integer :: k, g, ng
      real(dp) :: s

      r = rank_average(x)
      ug = unique_groups_local(group)
      s = 0.0_dp
      do k = 1, size(ug)
         g = ug(k)
         ng = count(group == g)
         if (ng > 0) s = s + sum(r, mask=group == g)**2 / real(ng, dp)
      end do
      res%statistic = 12.0_dp * s / (real(size(x),dp) * real(size(x) + 1,dp)) &
                    - 3.0_dp * real(size(x) + 1,dp)
      res%df1 = real(size(ug) - 1, dp)
      res%pvalue = 1.0_dp - chisq_cdf(res%statistic, res%df1)
   end function kruskal_test

   function poisson_dispersion_test(x) result(res)
      integer, intent(in) :: x(:)
      type(test_result) :: res
      real(dp) :: m

      m = sum(real(x,dp)) / real(size(x),dp)
      res%statistic = sum((real(x,dp) - m)**2 / max(m, tiny(1.0_dp)))
      res%df1 = real(size(x) - 1,dp)
      res%pvalue = 1.0_dp - chisq_cdf(res%statistic, res%df1)
      res%estimate = res%statistic / res%df1
   end function poisson_dispersion_test

   function unique_groups_local(g) result(u)
      integer, intent(in) :: g(:)
      integer, allocatable :: u(:), tmp(:)
      integer :: i, j, n, key

      tmp = g
      do i = 2, size(tmp)
         key = tmp(i)
         j = i - 1
         do while (j >= 1)
            if (tmp(j) <= key) exit
            tmp(j+1) = tmp(j)
            j = j - 1
         end do
         tmp(j+1) = key
      end do
      n = 0
      do i = 1, size(tmp)
         if (i == 1) then
            n = n + 1
         else if (tmp(i) /= tmp(i-1)) then
            n = n + 1
         end if
      end do
      allocate(u(n))
      n = 0
      do i = 1, size(tmp)
         if (i == 1) then
            n = n + 1
            u(n) = tmp(i)
         else if (tmp(i) /= tmp(i-1)) then
            n = n + 1
            u(n) = tmp(i)
         end if
      end do
   end function unique_groups_local

end module rfast_tests
