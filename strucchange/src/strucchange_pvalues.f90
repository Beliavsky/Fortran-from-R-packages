! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Derived from the R package strucchange 1.6-0. See NOTICE.md and UPSTREAM.md.
module strucchange_pvalues
   use r_kinds, only : dp
   use r_distributions, only : r_pchisq, r_pnorm
   use strucchange_tables, only : sc_beta_sup_1, sc_beta_sup_2, sc_beta_sup_3
   use strucchange_tables, only : sc_beta_ave_1, sc_beta_ave_2, sc_beta_ave_3
   use strucchange_tables, only : sc_beta_ave_4, sc_beta_ave_5
   use strucchange_tables, only : sc_beta_exp_1, sc_beta_exp_2, sc_beta_exp_3
   use strucchange_tables, only : sc_beta_exp_4, sc_beta_exp_5
   use strucchange_tables, only : sc_me, sc_mean_l2, sc_max_l2
   use strucchange_tables, only : mean_l2_probabilities, max_l2_probabilities
   use strucchange_utils, only : linear_interp
   implicit none
   private
   public :: efp_pvalue
   public :: fstats_pvalue
contains
   real(dp) function fstats_pvalue(x, test_type, k, lambda) result(p)
      real(dp), intent(in) :: x, lambda
      character(len = *), intent(in) :: test_type
      integer, intent(in) :: k
      real(dp) :: dummy(25), pp(25), tau, taua
      integer :: i, k_use, row0, tau1

      k_use = min(max(k, 1), 40)
      select case (trim(test_type))
      case ("supF", "supf")
         row0 = (k_use - 1) * 25
         do i = 1, 25
            dummy(i) = sc_beta_sup_1(row0 + i) + &
               sc_beta_sup_2(row0 + i) * x
            dummy(i) = max(dummy(i), 0.0_dp)
            pp(i) = r_pchisq(dummy(i), sc_beta_sup_3(row0 + i), &
               lower_tail = .false.)
         end do
      case ("aveF", "avef")
         row0 = (k_use - 1) * 25
         do i = 1, 25
            dummy(i) = sc_beta_ave_1(row0 + i) + &
               sc_beta_ave_2(row0 + i) * x + sc_beta_ave_3(row0 + i) * x ** 2 + &
               sc_beta_ave_4(row0 + i) * x ** 3
            dummy(i) = max(dummy(i), 0.0_dp)
            pp(i) = r_pchisq(dummy(i), sc_beta_ave_5(row0 + i), &
               lower_tail = .false.)
         end do
      case ("expF", "expf")
         row0 = (k_use - 1) * 25
         do i = 1, 25
            dummy(i) = sc_beta_exp_1(row0 + i) + &
               sc_beta_exp_2(row0 + i) * x + sc_beta_exp_3(row0 + i) * x ** 2 + &
               sc_beta_exp_4(row0 + i) * x ** 3
            dummy(i) = max(dummy(i), 0.0_dp)
            pp(i) = r_pchisq(dummy(i), sc_beta_exp_5(row0 + i), &
               lower_tail = .false.)
         end do
      case default
         p = 1.0_dp
         return
      end select

      if (lambda < 1.0_dp) then
         tau = lambda
      else
         tau = 1.0_dp / (1.0_dp + sqrt(lambda))
      end if
      if (abs(tau - 0.5_dp) <= 8.0_dp * epsilon(1.0_dp)) then
         p = r_pchisq(x, real(k_use, dp), lower_tail = .false.)
      else if (tau <= 0.01_dp) then
         p = pp(25)
      else if (tau >= 0.49_dp) then
         p = ((0.5_dp - tau) * pp(1) + (tau - 0.49_dp) * &
            r_pchisq(x, real(k_use, dp), lower_tail = .false.)) * 100.0_dp
      else
         taua = (0.51_dp - tau) * 50.0_dp
         tau1 = floor(taua)
         tau1 = min(max(tau1, 1), 24)
         p = (real(tau1 + 1, dp) - taua) * pp(tau1) + &
            (taua - real(tau1, dp)) * pp(tau1 + 1)
      end if
      p = min(max(p, 0.0_dp), 1.0_dp)
   end function fstats_pvalue

   real(dp) function efp_pvalue(x, limiting_process, functional, k, h, &
      alternative_boundary) result(p)
      real(dp), intent(in) :: x
      character(len = *), intent(in) :: limiting_process, functional
      integer, intent(in) :: k
      real(dp), intent(in), optional :: h
      logical, intent(in), optional :: alternative_boundary
      real(dp) :: h_use
      logical :: alternative

      alternative = .false.
      if (present(alternative_boundary)) alternative = alternative_boundary
      h_use = 0.15_dp
      if (present(h)) h_use = h
      select case (trim(limiting_process))
      case ("Brownian motion")
         p = brownian_motion_pvalue(x, k, alternative, functional)
      case ("Brownian bridge")
         p = brownian_bridge_pvalue(x, k, alternative, functional)
      case ("Brownian motion increments")
         p = brownian_motion_increments_pvalue(x, h_use, k, alternative, &
            functional)
      case ("Brownian bridge increments")
         p = brownian_bridge_increments_pvalue(x, h_use, k, functional)
      case default
         p = 1.0_dp
      end select
      p = min(max(p, 0.0_dp), 1.0_dp)
   end function efp_pvalue

   real(dp) function brownian_motion_pvalue(x, k, alternative, functional) &
      result(p)
      real(dp), intent(in) :: x
      integer, intent(in) :: k
      logical, intent(in) :: alternative
      character(len = *), intent(in) :: functional
      real(dp), parameter :: pval(35) = [ &
         1.0_dp, 0.997_dp, 0.99_dp, 0.975_dp, 0.949_dp, 0.912_dp, &
         0.864_dp, 0.806_dp, 0.739_dp, 0.666_dp, 0.589_dp, 0.512_dp, &
         0.437_dp, 0.368_dp, 0.307_dp, 0.253_dp, 0.205_dp, 0.163_dp, &
         0.129_dp, 0.100_dp, 0.077_dp, 0.058_dp, 0.043_dp, 0.032_dp, &
         0.024_dp, 0.018_dp, 0.012_dp, 0.009_dp, 0.006_dp, 0.004_dp, &
         0.003_dp, 0.002_dp, 0.001_dp, 0.001_dp, 0.001_dp ]
      real(dp) :: crit(35), base
      integer :: i

      if (trim(functional) /= "max") then
         p = 1.0_dp
         return
      end if
      if (alternative) then
         do i = 1, 35
            crit(i) = real(i + 9, dp) / 10.0_dp
         end do
         base = linear_interp(crit, pval, x)
      else if (x < 0.3_dp) then
         base = 1.0_dp - 0.1465_dp * x
      else
         base = 2.0_dp * (1.0_dp - r_pnorm(3.0_dp * x) + &
            exp(-4.0_dp * x * x) * (r_pnorm(x) + r_pnorm(5.0_dp * x) - 1.0_dp) - &
            exp(-16.0_dp * x * x) * (1.0_dp - r_pnorm(x)))
      end if
      p = 1.0_dp - (1.0_dp - base) ** max(k, 1)
   end function brownian_motion_pvalue

   real(dp) function brownian_bridge_pvalue(x, k, alternative, functional) &
      result(p)
      real(dp), intent(in) :: x
      integer, intent(in) :: k
      logical, intent(in) :: alternative
      character(len = *), intent(in) :: functional
      real(dp), parameter :: pval(35) = [ &
         1.0_dp, 1.0_dp, 0.997_dp, 0.99_dp, 0.977_dp, 0.954_dp, &
         0.919_dp, 0.871_dp, 0.812_dp, 0.743_dp, 0.666_dp, 0.585_dp, &
         0.504_dp, 0.426_dp, 0.353_dp, 0.288_dp, 0.230_dp, 0.182_dp, &
         0.142_dp, 0.109_dp, 0.082_dp, 0.062_dp, 0.046_dp, 0.034_dp, &
         0.025_dp, 0.017_dp, 0.011_dp, 0.008_dp, 0.005_dp, 0.004_dp, &
         0.003_dp, 0.002_dp, 0.001_dp, 0.001_dp, 0.0001_dp ]
      real(dp) :: base, crit(35), series, term, xgrid(24), ygrid(24)
      integer :: i, k_use

      k_use = max(k, 1)
      select case (trim(functional))
      case ("max")
         if (alternative) then
            do i = 1, 35
               crit(i) = real(i + 11, dp) / 10.0_dp
            end do
            base = linear_interp(crit, pval, x)
            p = 1.0_dp - (1.0_dp - base) ** k_use
         else if (x < 0.1_dp) then
            p = 1.0_dp
         else
            series = 0.0_dp
            do i = 1, 100
               term = exp(-2.0_dp * real(i * i, dp) * x * x)
               if (mod(i, 2) == 1) term = -term
               series = series + term
            end do
            p = 1.0_dp - (1.0_dp + 2.0_dp * series) ** k_use
         end if
      case ("range")
         if (x < 0.4_dp) then
            p = 1.0_dp
         else
            series = 0.0_dp
            do i = 1, 10
               series = series + (4.0_dp * real(i * i, dp) * x * x - 1.0_dp) * &
                  exp(-2.0_dp * real(i * i, dp) * x * x)
            end do
            p = 1.0_dp - (1.0_dp - 2.0_dp * series) ** k_use
         end if
      case ("maxL2")
         k_use = min(k_use, 25)
         xgrid(1) = 0.0_dp
         ygrid(1) = 1.0_dp
         do i = 1, 23
            xgrid(i + 1) = sc_max_l2(k_use, i)
            ygrid(i + 1) = 1.0_dp - max_l2_probabilities(i)
         end do
         p = linear_interp(xgrid, ygrid, x)
      case ("meanL2")
         k_use = min(k_use, 25)
         p = mean_l2_pvalue(x, k_use)
      case default
         p = 1.0_dp
      end select
   end function brownian_bridge_pvalue

   real(dp) function mean_l2_pvalue(x, k) result(p)
      real(dp), intent(in) :: x
      integer, intent(in) :: k
      real(dp) :: xgrid(16), ygrid(16)
      integer :: i

      xgrid(1) = 0.0_dp
      ygrid(1) = 1.0_dp
      do i = 1, 15
         xgrid(i + 1) = sc_mean_l2(k, i)
         ygrid(i + 1) = 1.0_dp - mean_l2_probabilities(i)
      end do
      p = linear_interp(xgrid, ygrid, x)
   end function mean_l2_pvalue

   real(dp) function brownian_motion_increments_pvalue(x, h, k, &
      alternative, functional) result(p)
      real(dp), intent(in) :: x, h
      integer, intent(in) :: k
      logical, intent(in) :: alternative
      character(len = *), intent(in) :: functional
      real(dp), parameter :: table(10, 6) = reshape([ &
         3.2165_dp, 2.9795_dp, 2.8289_dp, 2.7099_dp, 2.6061_dp, &
         2.5111_dp, 2.4283_dp, 2.3464_dp, 2.2686_dp, 2.2255_dp, &
         3.3185_dp, 3.0894_dp, 2.9479_dp, 2.8303_dp, 2.7325_dp, &
         2.6418_dp, 2.5609_dp, 2.4840_dp, 2.4083_dp, 2.3668_dp, &
         3.4554_dp, 3.2368_dp, 3.1028_dp, 2.9874_dp, 2.8985_dp, &
         2.8134_dp, 2.7327_dp, 2.6605_dp, 2.5899_dp, 2.5505_dp, &
         3.6622_dp, 3.4681_dp, 3.3382_dp, 3.2351_dp, 3.1531_dp, &
         3.0728_dp, 3.0043_dp, 2.9333_dp, 2.8743_dp, 2.8334_dp, &
         3.8632_dp, 3.6707_dp, 3.5598_dp, 3.4604_dp, 3.3845_dp, &
         3.3102_dp, 3.2461_dp, 3.1823_dp, 3.1229_dp, 3.0737_dp, &
         4.1009_dp, 3.9397_dp, 3.8143_dp, 3.7337_dp, 3.6626_dp, &
         3.5907_dp, 3.5333_dp, 3.4895_dp, 3.4123_dp, 3.3912_dp &
      ], [10, 6])
      real(dp), parameter :: table_p(6) = [ &
         0.2_dp, 0.15_dp, 0.1_dp, 0.05_dp, 0.025_dp, 0.01_dp ]
      real(dp) :: hgrid(10), interpolated(6), scaled(10), xgrid(7)
      real(dp) :: ygrid(7), base
      integer :: i, j

      if (alternative .or. trim(functional) /= "max") then
         p = 1.0_dp
         return
      end if
      do i = 1, 10
         hgrid(i) = 0.05_dp * real(i, dp)
      end do
      do j = 1, 6
         do i = 1, 10
            scaled(i) = table(i, j) * sqrt(hgrid(i))
         end do
         interpolated(j) = linear_interp(hgrid, scaled, h)
      end do
      xgrid(1) = 0.0_dp
      xgrid(2:7) = interpolated
      ygrid(1) = 1.0_dp
      ygrid(2:7) = table_p
      base = linear_interp(xgrid, ygrid, x)
      p = 1.0_dp - (1.0_dp - base) ** max(k, 1)
   end function brownian_motion_increments_pvalue

   real(dp) function brownian_bridge_increments_pvalue(x, h, k, functional) &
      result(p)
      real(dp), intent(in) :: x, h
      integer, intent(in) :: k
      character(len = *), intent(in) :: functional
      real(dp), parameter :: table_p(4) = [ &
         0.1_dp, 0.05_dp, 0.025_dp, 0.01_dp ]
      real(dp) :: hgrid(10), interpolated(4), xgrid(5), ygrid(5), series
      integer :: i, j, k_use, row0

      k_use = min(max(k, 1), 6)
      select case (trim(functional))
      case ("max")
         row0 = (k_use - 1) * 10
         do i = 1, 10
            hgrid(i) = 0.05_dp * real(i, dp)
         end do
         do j = 1, 4
            interpolated(j) = linear_interp(hgrid, &
               sc_me(row0 + 1:row0 + 10, j), h)
         end do
         xgrid(1) = 0.0_dp
         xgrid(2:5) = interpolated
         ygrid(1) = 1.0_dp
         ygrid(2:5) = table_p
         p = linear_interp(xgrid, ygrid, x)
      case ("range")
         if (x < 0.53_dp) then
            p = 1.0_dp
         else
            series = 0.0_dp
            do i = 1, 10
               if (mod(i, 2) == 1) then
                  series = series + real(i, dp) * r_pnorm(-real(i, dp) * x)
               else
                  series = series - real(i, dp) * r_pnorm(-real(i, dp) * x)
               end if
            end do
            p = 1.0_dp - (1.0_dp - 8.0_dp * series) ** k_use
         end if
      case default
         p = 1.0_dp
      end select
   end function brownian_bridge_increments_pvalue
end module strucchange_pvalues
