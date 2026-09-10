! SPDX-License-Identifier: GPL-2.0-or-later
! Additional multiple-testing threshold rules from wavethresh 4.7.3.
module wavethresh_threshold_extra
   use wavethresh_types, only : dp, wd_t
   use wavethresh_stats, only : shrink_soft, kolsmi_chi2
   implicit none
   private

   public :: to_getthrda1, to_one_by_one1
   public :: to_getthrda2, to_one_by_one2
   public :: to_threshda1, to_threshda1_thresholds
   public :: to_threshda2, to_threshda2_thresholds

   real(dp), parameter :: to2_alpha_table(172) = [ &
      0.99999899999999997_dp, 0.999996_dp, 0.99999099999999996_dp, 0.99997899999999995_dp, &
      0.99995400000000001_dp, 0.99990900000000005_dp, 0.99982899999999997_dp, 0.99969699999999995_dp, &
      0.99948899999999996_dp, 0.99917400000000001_dp, 0.99871500000000002_dp, 0.99807100000000004_dp, &
      0.99719199999999997_dp, 0.99602800000000002_dp, 0.99452399999999996_dp, 0.99262300000000003_dp, &
      0.99026999999999998_dp, 0.98741000000000001_dp, 0.98399499999999995_dp, 0.97997800000000002_dp, &
      0.97531800000000002_dp, 0.96998300000000004_dp, 0.96394500000000005_dp, 0.95718599999999998_dp, &
      0.94969400000000004_dp, 0.94146600000000003_dp, 0.93250299999999997_dp, 0.922817_dp, &
      0.91242299999999998_dp, 0.90134400000000003_dp, 0.88960499999999998_dp, 0.87724000000000002_dp, &
      0.86428199999999999_dp, 0.85077100000000005_dp, 0.83677500000000005_dp, 0.82224699999999995_dp, &
      0.80732300000000001_dp, 0.79201299999999997_dp, 0.77636300000000003_dp, 0.76041800000000004_dp, &
      0.74421999999999999_dp, 0.72781099999999999_dp, 0.71123499999999995_dp, 0.69452899999999995_dp, &
      0.67773499999999998_dp, 0.660887_dp, 0.64401900000000001_dp, 0.62716700000000003_dp, &
      0.61036000000000001_dp, 0.59362800000000004_dp, 0.57699800000000001_dp, 0.56049499999999997_dp, &
      0.54414300000000004_dp, 0.52795899999999996_dp, 0.51197000000000004_dp, 0.49619200000000002_dp, &
      0.48063400000000001_dp, 0.46531800000000001_dp, 0.45025599999999999_dp, 0.43545400000000001_dp, &
      0.42093000000000003_dp, 0.40668399999999999_dp, 0.39273000000000002_dp, 0.37907200000000002_dp, &
      0.36571399999999998_dp, 0.35266199999999998_dp, 0.339918_dp, 0.327484_dp, &
      0.31536399999999998_dp, 0.30355599999999999_dp, 0.29205999999999999_dp, 0.28087400000000001_dp, &
      0.27000000000000002_dp, 0.259434_dp, 0.24917400000000001_dp, 0.23921999999999999_dp, &
      0.22956599999999999_dp, 0.22020600000000001_dp, 0.21113999999999999_dp, 0.20236399999999999_dp, &
      0.19387199999999999_dp, 0.18565799999999999_dp, 0.17771799999999999_dp, 0.17005000000000001_dp, &
      0.16264400000000001_dp, 0.155498_dp, 0.14860599999999999_dp, 0.141962_dp, &
      0.13555800000000001_dp, 0.129388_dp, 0.12345200000000001_dp, 0.117742_dp, &
      0.11225_dp, 0.10697_dp, 0.101896_dp, 0.097028000000000003_dp, &
      0.092352000000000004_dp, 0.087868000000000002_dp, 0.083568000000000003_dp, 0.079444000000000001_dp, &
      0.075495000000000007_dp, 0.071711999999999998_dp, 0.068092_dp, 0.064630000000000007_dp, &
      0.061317999999999998_dp, 0.058152000000000002_dp, 0.055128000000000003_dp, 0.052243999999999999_dp, &
      0.049487999999999997_dp, 0.046857999999999997_dp, 0.044350000000000001_dp, 0.041959999999999997_dp, &
      0.039682000000000002_dp, 0.037513999999999999_dp, 0.035448_dp, 0.033484_dp, &
      0.031618_dp, 0.029842_dp, 0.028153999999999998_dp, 0.026551999999999999_dp, &
      0.02503_dp, 0.023588000000000001_dp, 0.022218000000000002_dp, 0.019689999999999999_dp, &
      0.017422_dp, 0.015389999999999999_dp, 0.013573999999999999_dp, 0.011952000000000001_dp, &
      0.010508_dp, 0.0092230000000000003_dp, 0.0080829999999999999_dp, 0.0070720000000000002_dp, &
      0.0061770000000000002_dp, 0.0053880000000000004_dp, 0.0046909999999999999_dp, 0.004078_dp, &
      0.0035400000000000002_dp, 0.003068_dp, 0.0026540000000000001_dp, 0.0022929999999999999_dp, &
      0.001977_dp, 0.0017030000000000001_dp, 0.001464_dp, 0.001256_dp, &
      0.0010759999999999999_dp, 0.00092100000000000005_dp, 0.00078700000000000005_dp, 0.00067100000000000005_dp, &
      0.00057200000000000003_dp, 0.000484_dp, 0.00041199999999999999_dp, 0.00035_dp, &
      0.00029500000000000001_dp, 0.00025000000000000001_dp, 0.00021000000000000001_dp, 0.00017799999999999999_dp, &
      0.00014799999999999999_dp, 0.000126_dp, 0.00010399999999999999_dp, 8.7999999999999998e-05_dp, &
      7.3999999999999996e-05_dp, 6.0000000000000002e-05_dp, 5.1e-05_dp, 4.1999999999999998e-05_dp, &
      3.4999999999999997e-05_dp, 3.0000000000000001e-05_dp, 2.4000000000000001e-05_dp, 2.0000000000000002e-05_dp, &
      1.5999999999999999e-05_dp, 1.2999999999999999e-05_dp, 1.1e-05_dp, 9.0000000000000002e-06_dp &
   ]

contains

   pure function to_one_by_one1(dat, alpha) result(cc)
      real(dp), intent(in) :: dat(:) !! Ascending squared coefficients used by the sequential chi-square test.
      real(dp), intent(in) :: alpha !! Sequential-test significance probability in the open interval (0, 1).
      real(dp), allocatable :: cc(:)
      real(dp), allocatable :: work(:)
      integer :: i
      integer :: count_values

      if (size(dat) == 0) then
         allocate(cc(0))
         return
      end if
      allocate(work(size(dat)))
      i = size(dat)
      count_values = 1
      work(count_values) = 1.0_dp - pchisq1(dat(i))**i
      do while (work(count_values) < alpha .and. i > 1)
         i = i - 1
         count_values = count_values + 1
         work(count_values) = 1.0_dp - pchisq1(dat(i))**i
      end do
      cc = work(1:count_values)
   end function to_one_by_one1

   pure function to_getthrda1(dat, alpha) result(threshold)
      real(dp), intent(in) :: dat(:) !! Wavelet coefficients whose squared magnitudes determine the level threshold.
      real(dp), intent(in) :: alpha !! Sequential-test significance probability in the open interval (0, 1).
      real(dp) :: threshold
      real(dp), allocatable :: squared(:)
      real(dp), allocatable :: cc(:)
      integer :: index

      if (size(dat) == 0 .or. alpha <= 0.0_dp .or. alpha >= 1.0_dp) then
         threshold = 0.0_dp
         return
      end if
      squared = dat**2
      call sort_ascending(squared)
      cc = to_one_by_one1(squared, alpha)
      if (size(cc) == size(squared)) then
         if (1.0_dp - pchisq1(squared(1)) < alpha) then
            threshold = 0.0_dp
         else
            threshold = sqrt(max(squared(1), 0.0_dp))
         end if
      else
         index = size(squared) - size(cc) + 1
         threshold = sqrt(max(squared(index), 0.0_dp))
      end if
   end function to_getthrda1

   pure function to_threshda1(object, alpha) result(thresholded)
      type(wd_t), intent(in) :: object !! Decimated wavelet object whose detail levels are soft-thresholded independently.
      real(dp), intent(in), optional :: alpha !! Sequential-test significance probability; default is 0.05.
      type(wd_t) :: thresholded
      real(dp) :: probability
      real(dp) :: threshold
      integer :: level

      thresholded = object
      if (.not. object%ok) then
         thresholded%message = "TOthreshda1 requires a valid wd object"
         thresholded%ok = .false.
         return
      end if
      probability = 0.05_dp
      if (present(alpha)) probability = alpha
      if (probability <= 0.0_dp .or. probability >= 1.0_dp) then
         thresholded%message = "alpha must lie strictly between zero and one"
         thresholded%ok = .false.
         return
      end if
      do level = object%nlevels - 1, 0, -1
         threshold = to_getthrda1(object%detail(level)%values, probability)
         thresholded%detail(level)%values = shrink_soft(object%detail(level)%values, threshold)
      end do
      thresholded%ok = .true.
      thresholded%message = "ok"
   end function to_threshda1

   pure function to_threshda1_thresholds(object, alpha) result(thresholds)
      type(wd_t), intent(in) :: object !! Decimated wavelet object whose per-level TOthreshda1 thresholds are requested.
      real(dp), intent(in), optional :: alpha !! Sequential-test significance probability; default is 0.05.
      real(dp), allocatable :: thresholds(:)
      real(dp) :: probability
      integer :: level

      if (.not. object%ok) then
         allocate(thresholds(0))
         return
      end if
      probability = 0.05_dp
      if (present(alpha)) probability = alpha
      if (probability <= 0.0_dp .or. probability >= 1.0_dp) then
         allocate(thresholds(0))
         return
      end if
      allocate(thresholds(object%nlevels))
      do level = 0, object%nlevels - 1
         thresholds(level + 1) = to_getthrda1(object%detail(level)%values, probability)
      end do
   end function to_threshda1_thresholds


   pure function to_one_by_one2(dat, alpha) result(cc)
      real(dp), intent(in) :: dat(:) !! Squared coefficients in original order for the sequential goodness-of-fit test.
      real(dp), intent(in) :: alpha !! Significance probability; must lie between the extrema of the upstream empirical table.
      real(dp), allocatable :: cc(:)
      real(dp), allocatable :: work(:)
      real(dp), allocatable :: subset(:)
      real(dp) :: critical_value
      integer :: i
      integer :: ind
      integer :: count_values

      if (size(dat) == 0 .or. alpha <= to2_alpha_table(172) .or. alpha >= to2_alpha_table(1)) then
         allocate(cc(0))
         return
      end if
      ind = 0
      do i = 2, size(to2_alpha_table)
         if (alpha > to2_alpha_table(i)) then
            ind = i
            exit
         end if
      end do
      if (ind == 0) then
         allocate(cc(0))
         return
      end if
      critical_value = to2_critical(ind - 1) + &
         (to2_alpha_table(ind - 1) - alpha) * (to2_critical(ind) - to2_critical(ind - 1)) / &
         (to2_alpha_table(ind - 1) - to2_alpha_table(ind))

      allocate(work(size(dat)))
      i = size(dat)
      count_values = 1
      work(count_values) = kolsmi_chi2(dat)
      do while (work(count_values) > critical_value .and. i > 1)
         i = i - 1
         count_values = count_values + 1
         subset = smallest_subset_original_order(dat, i)
         work(count_values) = kolsmi_chi2(subset)
      end do
      cc = work(1:count_values)
   end function to_one_by_one2

   pure function to_getthrda2(dat, alpha) result(threshold)
      real(dp), intent(in) :: dat(:) !! Squared wavelet coefficients used by the second sequential threshold rule.
      real(dp), intent(in) :: alpha !! Significance probability within the upstream empirical calibration range.
      real(dp) :: threshold
      real(dp), allocatable :: cc(:)
      real(dp), allocatable :: sorted(:)
      integer :: rank

      if (size(dat) == 0 .or. alpha <= to2_alpha_table(172) .or. alpha >= to2_alpha_table(1)) then
         threshold = 0.0_dp
         return
      end if
      cc = to_one_by_one2(dat, alpha)
      sorted = dat
      call sort_ascending(sorted)
      if (size(cc) == size(dat)) then
         if (1.0_dp - pchisq1(sorted(1)) < alpha) then
            threshold = 0.0_dp
         else
            threshold = sqrt(max(sorted(1), 0.0_dp))
         end if
      else
         rank = size(dat) - size(cc) + 1
         threshold = sqrt(max(sorted(rank), 0.0_dp))
      end if
   end function to_getthrda2

   pure function to_threshda2(object, alpha) result(thresholded)
      type(wd_t), intent(in) :: object !! Decimated wavelet object whose detail levels are soft-thresholded by TOthreshda2.
      real(dp), intent(in), optional :: alpha !! Empirical-test significance probability; default is 0.05.
      type(wd_t) :: thresholded
      real(dp) :: probability
      real(dp) :: threshold
      integer :: level

      thresholded = object
      if (.not. object%ok) then
         thresholded%message = "TOthreshda2 requires a valid wd object"
         thresholded%ok = .false.
         return
      end if
      probability = 0.05_dp
      if (present(alpha)) probability = alpha
      if (probability <= to2_alpha_table(172) .or. probability >= to2_alpha_table(1)) then
         thresholded%message = "alpha is outside the TOthreshda2 calibration range"
         thresholded%ok = .false.
         return
      end if
      do level = object%nlevels - 1, 0, -1
         threshold = to_getthrda2(object%detail(level)%values**2, probability)
         thresholded%detail(level)%values = shrink_soft(object%detail(level)%values, threshold)
      end do
      thresholded%ok = .true.
      thresholded%message = "ok"
   end function to_threshda2

   pure function to_threshda2_thresholds(object, alpha) result(thresholds)
      type(wd_t), intent(in) :: object !! Decimated wavelet object whose per-level TOthreshda2 thresholds are requested.
      real(dp), intent(in), optional :: alpha !! Empirical-test significance probability; default is 0.05.
      real(dp), allocatable :: thresholds(:)
      real(dp) :: probability
      integer :: level

      if (.not. object%ok) then
         allocate(thresholds(0))
         return
      end if
      probability = 0.05_dp
      if (present(alpha)) probability = alpha
      if (probability <= to2_alpha_table(172) .or. probability >= to2_alpha_table(1)) then
         allocate(thresholds(0))
         return
      end if
      allocate(thresholds(object%nlevels))
      do level = 0, object%nlevels - 1
         thresholds(level + 1) = to_getthrda2(object%detail(level)%values**2, probability)
      end do
   end function to_threshda2_thresholds

   pure elemental function to2_critical(index) result(value)
      integer, intent(in) :: index !! One-based index into the 172-entry upstream TOonebyone2 critical-value grid.
      real(dp) :: value

      if (index <= 122) then
         value = 0.28_dp + 0.01_dp * real(index - 1, dp)
      else
         value = 1.5_dp + 0.02_dp * real(index - 123, dp)
      end if
   end function to2_critical

   pure function smallest_subset_original_order(values, nkeep) result(subset)
      real(dp), intent(in) :: values(:) !! Values from which the nkeep smallest entries are selected.
      integer, intent(in) :: nkeep !! Number of smallest entries to retain while restoring their original order.
      real(dp), allocatable :: subset(:)
      integer, allocatable :: order(:)
      integer, allocatable :: kept(:)
      integer :: i
      integer :: j
      integer :: key

      if (nkeep <= 0) then
         allocate(subset(0))
         return
      end if
      allocate(order(size(values)))
      order = [(i, i=1, size(values))]
      do i = 2, size(order)
         key = order(i)
         j = i - 1
         do while (j >= 1)
            if (values(order(j)) <= values(key)) exit
            order(j + 1) = order(j)
            j = j - 1
         end do
         order(j + 1) = key
      end do
      kept = order(1:min(nkeep, size(order)))
      call sort_integer_ascending(kept)
      allocate(subset(size(kept)))
      do i = 1, size(kept)
         subset(i) = values(kept(i))
      end do
   end function smallest_subset_original_order

   pure subroutine sort_integer_ascending(values)
      integer, intent(inout) :: values(:) !! Integer indices sorted into ascending order in place.
      integer :: i
      integer :: j
      integer :: key

      do i = 2, size(values)
         key = values(i)
         j = i - 1
         do while (j >= 1)
            if (values(j) <= key) exit
            values(j + 1) = values(j)
            j = j - 1
         end do
         values(j + 1) = key
      end do
   end subroutine sort_integer_ascending

   pure elemental function pchisq1(x) result(probability)
      real(dp), intent(in) :: x !! Chi-square variate with one degree of freedom; negative values map to probability zero.
      real(dp) :: probability

      if (x <= 0.0_dp) then
         probability = 0.0_dp
      else
         probability = erf(sqrt(0.5_dp * x))
      end if
   end function pchisq1

   pure subroutine sort_ascending(values)
      real(dp), intent(inout) :: values(:) !! Real vector sorted in ascending order in place.
      real(dp) :: key
      integer :: i
      integer :: j

      do i = 2, size(values)
         key = values(i)
         j = i - 1
         do while (j >= 1)
            if (values(j) <= key) exit
            values(j + 1) = values(j)
            j = j - 1
         end do
         values(j + 1) = key
      end do
   end subroutine sort_ascending

end module wavethresh_threshold_extra
