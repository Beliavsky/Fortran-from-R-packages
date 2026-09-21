module ksamples_ad
   use r_kinds, only : dp
   use ksamples_types, only : ad_result, sample_block
   use ksamples_utils, only : multinomial_count, shuffle_real, convolve_values
   use ksamples_utils, only : initial_group_labels, next_group_labels, grouped_values
   implicit none
   private

   public :: ad_pval, ad_statistic, ad_test, ad_test_combined

   real(dp), parameter :: ad_s_grid(8) = [ &
      1.0_dp, 0.7071067811865475_dp, 0.5773502691896258_dp, 0.5_dp, &
      0.4082482904638630_dp, 0.3535533905932738_dp, 0.3162277660168379_dp, 0.0_dp ]
   real(dp), parameter :: ad_tail(35) = [ &
      1e-05_dp, 5e-05_dp, 0.0001_dp, 0.0005_dp, 0.001_dp, &
      0.005_dp, 0.01_dp, 0.025_dp, 0.05_dp, 0.075_dp, &
      0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp, &
      0.6_dp, 0.7_dp, 0.8_dp, 0.9_dp, 0.925_dp, &
      0.95_dp, 0.975_dp, 0.99_dp, 0.9925_dp, 0.995_dp, &
      0.9975_dp, 0.999_dp, 0.99925_dp, 0.9995_dp, 0.99975_dp, &
      0.9999_dp, 0.999925_dp, 0.99995_dp, 0.999975_dp, 0.99999_dp &
   ]
   real(dp), parameter :: ad_table1(8, 35) = reshape([ &
      -1.1954_dp, -1.5806_dp, -1.8172_dp, -2.0032_dp, &
      -2.2526_dp, -2.4204_dp, -2.5283_dp, -4.2649_dp, &
      -1.1786_dp, -1.5394_dp, -1.7728_dp, -1.9426_dp, &
      -2.1685_dp, -2.3288_dp, -2.4374_dp, -3.8906_dp, &
      -1.166_dp, -1.5193_dp, -1.7462_dp, -1.9067_dp, &
      -2.126_dp, -2.2818_dp, -2.3926_dp, -3.719_dp, &
      -1.1407_dp, -1.4659_dp, -1.671_dp, -1.8105_dp, &
      -2.0048_dp, -2.1356_dp, -2.2348_dp, -3.2905_dp, &
      -1.1253_dp, -1.4371_dp, -1.6314_dp, -1.7619_dp, &
      -1.9396_dp, -2.0637_dp, -2.1521_dp, -3.0902_dp, &
      -1.0777_dp, -1.3503_dp, -1.5102_dp, -1.6177_dp, &
      -1.761_dp, -1.8537_dp, -1.9178_dp, -2.5758_dp, &
      -1.0489_dp, -1.2984_dp, -1.4415_dp, -1.5355_dp, &
      -1.6625_dp, -1.738_dp, -1.7936_dp, -2.3263_dp, &
      -0.9978_dp, -1.2098_dp, -1.3251_dp, -1.4007_dp, &
      -1.4977_dp, -1.5555_dp, -1.5941_dp, -1.96_dp, &
      -0.9417_dp, -1.1187_dp, -1.209_dp, -1.2671_dp, &
      -1.3382_dp, -1.379_dp, -1.405_dp, -1.6449_dp, &
      -0.8981_dp, -1.0491_dp, -1.1235_dp, -1.1692_dp, &
      -1.2249_dp, -1.2552_dp, -1.2755_dp, -1.4395_dp, &
      -0.8598_dp, -0.9904_dp, -1.0513_dp, -1.0879_dp, &
      -1.1317_dp, -1.155_dp, -1.1694_dp, -1.2816_dp, &
      -0.7258_dp, -0.7938_dp, -0.8188_dp, -0.8312_dp, &
      -0.8435_dp, -0.8471_dp, -0.8496_dp, -0.8416_dp, &
      -0.5966_dp, -0.617_dp, -0.6177_dp, -0.6139_dp, &
      -0.6073_dp, -0.5987_dp, -0.5941_dp, -0.5244_dp, &
      -0.4572_dp, -0.4383_dp, -0.419_dp, -0.4033_dp, &
      -0.3834_dp, -0.3676_dp, -0.3587_dp, -0.2533_dp, &
      -0.2966_dp, -0.2428_dp, -0.2078_dp, -0.1844_dp, &
      -0.1548_dp, -0.1346_dp, -0.1224_dp, 0.0_dp, &
      -0.1009_dp, -0.0169_dp, 0.0304_dp, 0.0596_dp, &
      0.0933_dp, 0.1156_dp, 0.1294_dp, 0.2533_dp, &
      0.1571_dp, 0.2635_dp, 0.3169_dp, 0.348_dp, &
      0.3823_dp, 0.4038_dp, 0.4166_dp, 0.5244_dp, &
      0.5357_dp, 0.6496_dp, 0.6992_dp, 0.7246_dp, &
      0.7528_dp, 0.7683_dp, 0.7771_dp, 0.8416_dp, &
      1.2255_dp, 1.2989_dp, 1.3202_dp, 1.3254_dp, &
      1.3305_dp, 1.3286_dp, 1.3257_dp, 1.2816_dp, &
      1.5262_dp, 1.5677_dp, 1.5709_dp, 1.5663_dp, &
      1.5561_dp, 1.5449_dp, 1.5356_dp, 1.4395_dp, &
      1.9633_dp, 1.943_dp, 1.919_dp, 1.8975_dp, &
      1.8641_dp, 1.8389_dp, 1.8212_dp, 1.6449_dp, &
      2.7314_dp, 2.5899_dp, 2.5_dp, 2.4451_dp, &
      2.3664_dp, 2.3155_dp, 2.2823_dp, 1.96_dp, &
      3.7825_dp, 3.4425_dp, 3.2582_dp, 3.1423_dp, &
      3.0036_dp, 2.9101_dp, 2.8579_dp, 2.3263_dp, &
      4.1241_dp, 3.716_dp, 3.4984_dp, 3.3651_dp, &
      3.2003_dp, 3.0928_dp, 3.0311_dp, 2.4324_dp, &
      4.6044_dp, 4.0847_dp, 3.8348_dp, 3.6714_dp, &
      3.4721_dp, 3.3453_dp, 3.2777_dp, 2.5758_dp, &
      5.409_dp, 4.7223_dp, 4.4022_dp, 4.1791_dp, &
      3.9357_dp, 3.7809_dp, 3.6963_dp, 2.807_dp, &
      6.4954_dp, 5.5823_dp, 5.1456_dp, 4.8657_dp, &
      4.5506_dp, 4.3275_dp, 4.2228_dp, 3.0902_dp, &
      6.8279_dp, 5.8282_dp, 5.3658_dp, 5.0749_dp, &
      4.7318_dp, 4.4923_dp, 4.3642_dp, 3.1747_dp, &
      7.2755_dp, 6.197_dp, 5.6715_dp, 5.3642_dp, &
      4.9991_dp, 4.7135_dp, 4.5945_dp, 3.2905_dp, &
      8.1885_dp, 6.8537_dp, 6.2077_dp, 5.8499_dp, &
      5.4246_dp, 5.1137_dp, 4.9555_dp, 3.4808_dp, &
      9.3061_dp, 7.6592_dp, 6.85_dp, 6.4806_dp, &
      5.9919_dp, 5.6122_dp, 5.5136_dp, 3.719_dp, &
      9.6132_dp, 7.9234_dp, 7.1025_dp, 6.6731_dp, &
      6.1549_dp, 5.8217_dp, 5.7345_dp, 3.7911_dp, &
      10.0989_dp, 8.2395_dp, 7.4326_dp, 6.9567_dp, &
      6.3908_dp, 6.011_dp, 5.9566_dp, 3.8906_dp, &
      10.8825_dp, 8.8994_dp, 7.8934_dp, 7.4501_dp, &
      6.9009_dp, 6.4538_dp, 6.2705_dp, 4.0556_dp, &
      11.8537_dp, 9.5482_dp, 8.5568_dp, 8.0283_dp, &
      7.4418_dp, 6.9524_dp, 6.6195_dp, 4.2649_dp &
   ], [8, 35])
   real(dp), parameter :: ad_table2(8, 35) = reshape([ &
      -1.1976_dp, -1.5824_dp, -1.8195_dp, -2.005_dp, &
      -2.2546_dp, -2.422_dp, -2.5292_dp, -4.2649_dp, &
      -1.1806_dp, -1.5416_dp, -1.7747_dp, -1.9434_dp, &
      -2.1687_dp, -2.3301_dp, -2.438_dp, -3.8906_dp, &
      -1.1681_dp, -1.5212_dp, -1.7479_dp, -1.9078_dp, &
      -2.1268_dp, -2.2827_dp, -2.3937_dp, -3.719_dp, &
      -1.1427_dp, -1.4677_dp, -1.6724_dp, -1.8115_dp, &
      -2.0059_dp, -2.1363_dp, -2.2359_dp, -3.2905_dp, &
      -1.1272_dp, -1.4387_dp, -1.6325_dp, -1.7629_dp, &
      -1.9405_dp, -2.0649_dp, -2.1527_dp, -3.0902_dp, &
      -1.0794_dp, -1.3518_dp, -1.5112_dp, -1.6187_dp, &
      -1.7617_dp, -1.8545_dp, -1.9182_dp, -2.5758_dp, &
      -1.0504_dp, -1.2997_dp, -1.4425_dp, -1.5362_dp, &
      -1.6632_dp, -1.7387_dp, -1.7943_dp, -2.3263_dp, &
      -0.999_dp, -1.2109_dp, -1.3259_dp, -1.4014_dp, &
      -1.4981_dp, -1.5561_dp, -1.5945_dp, -1.96_dp, &
      -0.9428_dp, -1.1196_dp, -1.2098_dp, -1.2677_dp, &
      -1.3386_dp, -1.3795_dp, -1.4054_dp, -1.6449_dp, &
      -0.8991_dp, -1.05_dp, -1.1241_dp, -1.1697_dp, &
      -1.2253_dp, -1.2557_dp, -1.2758_dp, -1.4395_dp, &
      -0.8607_dp, -0.9911_dp, -1.0518_dp, -1.0883_dp, &
      -1.1321_dp, -1.1555_dp, -1.1698_dp, -1.2816_dp, &
      -0.7264_dp, -0.7944_dp, -0.8192_dp, -0.8315_dp, &
      -0.8437_dp, -0.8473_dp, -0.8498_dp, -0.8416_dp, &
      -0.597_dp, -0.6173_dp, -0.6179_dp, -0.6141_dp, &
      -0.6074_dp, -0.5989_dp, -0.5942_dp, -0.5244_dp, &
      -0.4574_dp, -0.4385_dp, -0.4191_dp, -0.4034_dp, &
      -0.3835_dp, -0.3677_dp, -0.3588_dp, -0.2533_dp, &
      -0.2966_dp, -0.2427_dp, -0.2078_dp, -0.1844_dp, &
      -0.1548_dp, -0.1347_dp, -0.1224_dp, 0.0_dp, &
      -0.1007_dp, -0.0168_dp, 0.0305_dp, 0.0596_dp, &
      0.0934_dp, 0.1157_dp, 0.1295_dp, 0.2533_dp, &
      0.1573_dp, 0.2638_dp, 0.3171_dp, 0.3482_dp, &
      0.3825_dp, 0.404_dp, 0.4168_dp, 0.5244_dp, &
      0.5363_dp, 0.6501_dp, 0.6996_dp, 0.7249_dp, &
      0.753_dp, 0.7685_dp, 0.7773_dp, 0.8416_dp, &
      1.2263_dp, 1.2997_dp, 1.3209_dp, 1.3258_dp, &
      1.3309_dp, 1.329_dp, 1.326_dp, 1.2816_dp, &
      1.5274_dp, 1.5686_dp, 1.5716_dp, 1.5667_dp, &
      1.5565_dp, 1.5453_dp, 1.536_dp, 1.4395_dp, &
      1.9644_dp, 1.944_dp, 1.92_dp, 1.8983_dp, &
      1.8647_dp, 1.8396_dp, 1.8216_dp, 1.6449_dp, &
      2.7334_dp, 2.5915_dp, 2.5013_dp, 2.4457_dp, &
      2.3671_dp, 2.3162_dp, 2.2827_dp, 1.96_dp, &
      3.7851_dp, 3.4443_dp, 3.2595_dp, 3.1436_dp, &
      3.0046_dp, 2.9111_dp, 2.8585_dp, 2.3263_dp, &
      4.1255_dp, 3.7175_dp, 3.4997_dp, 3.3661_dp, &
      3.2011_dp, 3.0939_dp, 3.0318_dp, 2.4324_dp, &
      4.6067_dp, 4.0869_dp, 3.8363_dp, 3.6724_dp, &
      3.4729_dp, 3.3463_dp, 3.278_dp, 2.5758_dp, &
      5.4121_dp, 4.7248_dp, 4.4032_dp, 4.1812_dp, &
      3.9369_dp, 3.7819_dp, 3.6977_dp, 2.807_dp, &
      6.5_dp, 5.5856_dp, 5.1469_dp, 4.8683_dp, &
      4.552_dp, 4.3284_dp, 4.2229_dp, 3.0902_dp, &
      6.8324_dp, 5.8302_dp, 5.3678_dp, 5.0769_dp, &
      4.7332_dp, 4.4933_dp, 4.3654_dp, 3.1747_dp, &
      7.278_dp, 6.1999_dp, 5.674_dp, 5.3661_dp, &
      5.0001_dp, 4.7147_dp, 4.5956_dp, 3.2905_dp, &
      8.1926_dp, 6.8586_dp, 6.2082_dp, 5.8524_dp, &
      5.4265_dp, 5.115_dp, 4.9571_dp, 3.4808_dp, &
      9.3096_dp, 7.6673_dp, 6.8522_dp, 6.4825_dp, &
      5.9934_dp, 5.6135_dp, 5.5147_dp, 3.719_dp, &
      9.6207_dp, 7.929_dp, 7.1051_dp, 6.6763_dp, &
      6.1548_dp, 5.8229_dp, 5.7358_dp, 3.7911_dp, &
      10.1076_dp, 8.2437_dp, 7.4349_dp, 6.9593_dp, &
      6.3923_dp, 6.0136_dp, 5.9573_dp, 3.8906_dp, &
      10.8874_dp, 8.9034_dp, 7.8991_dp, 7.4543_dp, &
      6.9017_dp, 6.4568_dp, 6.2723_dp, 4.0556_dp, &
      11.8602_dp, 9.5499_dp, 8.5596_dp, 8.0315_dp, &
      7.4425_dp, 6.9537_dp, 6.6213_dp, 4.2649_dp &
   ], [8, 35])

contains

   pure real(dp) function ad_pval(tx, m, version) result(pvalue)
      real(dp), intent(in) :: tx !! Standardized Anderson-Darling statistic T_m.
      real(dp), intent(in) :: m !! Positive index m, conventionally the number of samples minus one.
      integer, intent(in), optional :: version !! Statistic version: 1 uses the first table; other values use the second.
      real(dp) :: s, w, tm(35), lp(35), y
      integer :: i, j, v
      v = 1
      if (present(version)) v = version
      if (m <= 0.0_dp) then
         pvalue = 1.0_dp
         return
      end if
      s = 1.0_dp/sqrt(m)
      s = min(1.0_dp, max(0.0_dp, s))
      do j = 1, 35
         if (s >= ad_s_grid(1)) then
            tm(j) = merge(ad_table1(1,j), ad_table2(1,j), v == 1)
         else if (s <= ad_s_grid(8)) then
            tm(j) = merge(ad_table1(8,j), ad_table2(8,j), v == 1)
         else
            do i = 1, 7
               if (s <= ad_s_grid(i) .and. s >= ad_s_grid(i + 1)) then
                  w = (s - ad_s_grid(i + 1))/(ad_s_grid(i) - ad_s_grid(i + 1))
                  if (v == 1) then
                     tm(j) = w*ad_table1(i,j) + (1.0_dp - w)*ad_table1(i + 1,j)
                  else
                     tm(j) = w*ad_table2(i,j) + (1.0_dp - w)*ad_table2(i + 1,j)
                  end if
                  exit
               end if
            end do
         end if
         lp(j) = log((1.0_dp - ad_tail(j))/ad_tail(j))
      end do
      y = lp(1)
      if (tx <= tm(1)) then
         y = lp(1) + (tx - tm(1))*(lp(2) - lp(1))/(tm(2) - tm(1))
      else if (tx >= tm(35)) then
         y = lp(34) + (tx - tm(34))*(lp(35) - lp(34))/(tm(35) - tm(34))
      else
         do j = 1, 34
            if (tx >= tm(j) .and. tx <= tm(j + 1)) then
               w = (tx - tm(j))/(tm(j + 1) - tm(j))
               y = (1.0_dp - w)*lp(j) + w*lp(j + 1)
               exit
            end if
         end do
      end if
      if (y >= 0.0_dp) then
         pvalue = 1.0_dp/(1.0_dp + exp(-y))
      else
         pvalue = exp(y)/(1.0_dp + exp(y))
      end if
   end function ad_pval

   pure subroutine ad_statistic(x, ns, statistic)
      real(dp), intent(in) :: x(:) !! Samples concatenated in group order.
      integer, intent(in) :: ns(:) !! Positive group sizes whose sum must equal size(x).
      real(dp), intent(out) :: statistic(2) !! Version-1 and version-2 Anderson-Darling statistics.
      real(dp), allocatable :: z(:)
      integer, allocatable :: fij(:, :), mult(:)
      integer :: i, j, k, l, n, offset
      real(dp) :: mij, maij, bj, baj, tmp, inner_sum, alt_sum, denom
      if (sum(ns) /= size(x) .or. size(ns) < 2 .or. any(ns <= 0)) &
         error stop 'ad_statistic: invalid group sizes'
      n = size(x)
      k = size(ns)
      call unique_values(x, z)
      l = size(z)
      allocate(fij(k,l), mult(l))
      fij = 0
      offset = 0
      do i = 1, k
         do j = 1, ns(i)
            call add_count(x(offset + j), z, fij(i,:))
         end do
         offset = offset + ns(i)
      end do
      do j = 1, l
         mult(j) = sum(fij(:,j))
      end do
      statistic = 0.0_dp
      do i = 1, k
         mij = 0.0_dp
         inner_sum = 0.0_dp
         alt_sum = 0.0_dp
         bj = 0.0_dp
         do j = 1, l
            mij = mij + real(fij(i,j), dp)
            maij = mij - 0.5_dp*real(fij(i,j), dp)
            bj = bj + real(mult(j), dp)
            baj = bj - 0.5_dp*real(mult(j), dp)
            if (j < l) then
               tmp = real(n,dp)*mij - real(ns(i),dp)*bj
               denom = bj*(real(n,dp) - bj)
               if (denom > 0.0_dp) inner_sum = inner_sum + real(mult(j),dp)*tmp*tmp/denom
            end if
            tmp = real(n,dp)*maij - real(ns(i),dp)*baj
            denom = baj*(real(n,dp) - baj) - 0.25_dp*real(n*mult(j),dp)
            if (denom > 0.0_dp) alt_sum = alt_sum + real(mult(j),dp)*tmp*tmp/denom
         end do
         statistic(1) = statistic(1) + inner_sum/real(ns(i),dp)
         statistic(2) = statistic(2) + alt_sum/real(ns(i),dp)
      end do
      statistic(1) = statistic(1)/real(n,dp)
      statistic(2) = real(n - 1,dp)*statistic(2)/real(n*n,dp)
   contains
      pure subroutine unique_values(a, values)
         real(dp), intent(in) :: a(:) !! Unsorted pooled observations.
         real(dp), allocatable, intent(out) :: values(:) !! Sorted distinct observations.
         real(dp), allocatable :: tmp(:)
         integer :: ii, jj, nn
         real(dp) :: yy
         tmp = a
         do ii = 2, size(tmp)
            yy = tmp(ii)
            jj = ii - 1
            do while (jj >= 1)
               if (tmp(jj) <= yy) exit
               tmp(jj + 1) = tmp(jj)
               jj = jj - 1
            end do
            tmp(jj + 1) = yy
         end do
         nn = 0
         do ii = 1, size(tmp)
            if (ii == 1) then
               nn = nn + 1
            else if (tmp(ii) /= tmp(ii - 1)) then
               nn = nn + 1
            end if
         end do
         allocate(values(nn))
         nn = 0
         do ii = 1, size(tmp)
            if (ii == 1) then
               nn = nn + 1
               values(nn) = tmp(ii)
            else if (tmp(ii) /= tmp(ii - 1)) then
               nn = nn + 1
               values(nn) = tmp(ii)
            end if
         end do
      end subroutine unique_values
      pure subroutine add_count(value, values, counts)
         real(dp), intent(in) :: value !! Observation to locate in the distinct-value support.
         real(dp), intent(in) :: values(:) !! Sorted distinct pooled observations.
         integer, intent(inout) :: counts(:) !! Per-support counts for one sample.
         integer :: jj
         do jj = 1, size(values)
            if (value == values(jj)) then
               counts(jj) = counts(jj) + 1
               return
            end if
         end do
      end subroutine add_count
   end subroutine ad_statistic

   pure real(dp) function ad_sigma(ns) result(sigma)
      integer, intent(in) :: ns(:) !! Positive group sizes used by the k-sample test.
      integer :: i, n, k
      real(dp) :: hsum, h, g, a, b, c, d, sig2
      n = sum(ns)
      k = size(ns)
      if (n == 3 .and. k == 2) then
         sigma = 0.3535534_dp
         return
      end if
      if (n <= 3) then
         sigma = 0.0_dp
         return
      end if
      hsum = sum(1.0_dp/real(ns,dp))
      h = 0.0_dp
      do i = 1, n - 1
         h = h + 1.0_dp/real(i,dp)
      end do
      g = 0.0_dp
      do i = 1, n - 2
         g = g + (h - harmonic(i))/real(n - i,dp)
      end do
      a = (4.0_dp*g - 6.0_dp)*real(k - 1,dp) + (10.0_dp - 6.0_dp*g)*hsum
      b = (2.0_dp*g - 4.0_dp)*real(k*k,dp) + 8.0_dp*h*real(k,dp) + &
         (2.0_dp*g - 14.0_dp*h - 4.0_dp)*hsum - 8.0_dp*h + 4.0_dp*g - 6.0_dp
      c = (6.0_dp*h + 2.0_dp*g - 2.0_dp)*real(k*k,dp) + &
         (4.0_dp*h - 4.0_dp*g + 6.0_dp)*real(k,dp) + (2.0_dp*h - 6.0_dp)*hsum + 4.0_dp*h
      d = (2.0_dp*h + 6.0_dp)*real(k*k,dp) - 4.0_dp*h*real(k,dp)
      sig2 = (a*real(n,dp)**3 + b*real(n,dp)**2 + c*real(n,dp) + d)/ &
         real((n - 1)*(n - 2)*(n - 3),dp)
      sigma = sqrt(max(0.0_dp, sig2))
   contains
      pure real(dp) function harmonic(nterm) result(value)
         integer, intent(in) :: nterm !! Number of reciprocal terms to sum.
         integer :: ii
         value = 0.0_dp
         do ii = 1, nterm
            value = value + 1.0_dp/real(ii,dp)
         end do
      end function harmonic
   end function ad_sigma

   subroutine ad_test(x, ns, method, nsim, keep_null, result)
      real(dp), intent(in) :: x(:) !! Samples concatenated in group order.
      integer, intent(in) :: ns(:) !! Positive group sizes whose sum equals size(x).
      character(len=*), intent(in), optional :: method !! Requested method: asymptotic, simulated, or exact.
      integer, intent(in), optional :: nsim !! Simulation cap and exact-enumeration threshold; default 10000.
      logical, intent(in), optional :: keep_null !! Retain the generated null statistics when true.
      type(ad_result), intent(out) :: result !! Test statistics, asymptotic and randomization p-values, and optional null sample.
      character(len=10) :: chosen
      integer :: simulations, ncount, nsave, stored
      integer, allocatable :: labels(:)
      logical :: keep, more
      real(dp) :: ncomb
      real(dp), allocatable :: work(:), null(:, :), grouped(:)
      real(dp) :: stat(2)
      call ad_statistic(x, ns, result%statistic)
      result%sigma = ad_sigma(ns)
      if (result%sigma > 0.0_dp) then
         result%standardized = (result%statistic - real(size(ns) - 1,dp))/result%sigma
         result%asymptotic_p(1) = ad_pval(result%standardized(1), real(size(ns) - 1,dp), 1)
         result%asymptotic_p(2) = ad_pval(result%standardized(2), real(size(ns) - 1,dp), 2)
      else
         result%standardized = 0.0_dp
         result%asymptotic_p = 1.0_dp
      end if
      chosen = 'asymptotic'
      if (present(method)) chosen = adjustl(method)
      simulations = 10000
      if (present(nsim)) simulations = max(1, nsim)
      keep = .false.
      if (present(keep_null)) keep = keep_null
      if (trim(chosen) == 'asymptotic') then
         result%method = 'asymptotic'
         return
      end if
      ncomb = multinomial_count(ns)
      if (trim(chosen) == 'exact' .and. ncomb <= real(simulations,dp) .and. ncomb <= real(huge(1),dp)) then
         nsave = nint(ncomb)
         allocate(null(nsave,2))
         call initial_group_labels(ns, labels)
         stored = 0
         more = .true.
         do while (more)
            call grouped_values(x, labels, size(ns), grouped)
            stored = stored + 1
            call ad_statistic(grouped, ns, stat)
            null(stored,:) = stat
            call next_group_labels(labels, more)
         end do
         result%randomization_p(1) = real(count(null(:,1) >= result%statistic(1)),dp)/real(nsave,dp)
         result%randomization_p(2) = real(count(null(:,2) >= result%statistic(2)),dp)/real(nsave,dp)
         result%method = 'exact'
      else
         allocate(null(simulations,2), work(size(x)))
         work = x
         do ncount = 1, simulations
            call shuffle_real(work)
            call ad_statistic(work, ns, null(ncount,:))
         end do
         result%randomization_p(1) = real(count(null(:,1) >= result%statistic(1)),dp)/real(simulations,dp)
         result%randomization_p(2) = real(count(null(:,2) >= result%statistic(2)),dp)/real(simulations,dp)
         result%method = 'simulated'
      end if
      if (keep) then
         allocate(result%null_dist(size(null,1),2))
         result%null_dist = null
      end if
   end subroutine ad_test

   subroutine ad_test_combined(blocks, method, nsim, keep_null, result)
      type(sample_block), intent(in) :: blocks(:) !! Independent blocks, each containing pooled values and group sizes.
      character(len=*), intent(in), optional :: method !! Requested method: asymptotic, simulated, or exact.
      integer, intent(in), optional :: nsim !! Simulation cap and exact joint-enumeration threshold; default 10000.
      logical, intent(in), optional :: keep_null !! Retain the generated combined null statistics when true.
      type(ad_result), intent(out) :: result !! Combined Anderson-Darling result for both statistic versions.
      character(len=10) :: chosen
      integer :: b, i, simulations, nsave
      logical :: keep
      real(dp) :: mu, joint_count, s2
      real(dp) :: st(2)
      real(dp), allocatable :: work(:), null(:, :), acc1(:), acc2(:), next(:)
      type(ad_result) :: one
      one = ad_result()
      result%statistic = 0.0_dp
      mu = 0.0_dp
      s2 = 0.0_dp
      joint_count = 1.0_dp
      do b = 1, size(blocks)
         call ad_statistic(blocks(b)%x, blocks(b)%ns, st)
         result%statistic = result%statistic + st
         mu = mu + real(size(blocks(b)%ns) - 1,dp)
         s2 = s2 + ad_sigma(blocks(b)%ns)**2
         joint_count = min(huge(1.0_dp), joint_count*multinomial_count(blocks(b)%ns))
      end do
      result%sigma = sqrt(s2)
      if (result%sigma > 0.0_dp) then
         result%standardized = (result%statistic - mu)/result%sigma
         result%asymptotic_p(1) = ad_pval(result%standardized(1), mu, 1)
         result%asymptotic_p(2) = ad_pval(result%standardized(2), mu, 2)
      end if
      chosen = 'asymptotic'
      if (present(method)) chosen = adjustl(method)
      simulations = 10000
      if (present(nsim)) simulations = max(1, nsim)
      keep = .false.
      if (present(keep_null)) keep = keep_null
      if (trim(chosen) == 'asymptotic') then
         result%method = 'asymptotic'
         return
      end if
      if (trim(chosen) == 'exact' .and. joint_count <= real(simulations,dp) .and. joint_count <= real(huge(1),dp)) then
         do b = 1, size(blocks)
            call ad_test(blocks(b)%x, blocks(b)%ns,'exact',simulations,.true.,one)
            if (b == 1) then
               acc1 = one%null_dist(:,1)
               acc2 = one%null_dist(:,2)
            else
               call convolve_values(acc1, one%null_dist(:,1), next)
               call move_alloc(next, acc1)
               call convolve_values(acc2, one%null_dist(:,2), next)
               call move_alloc(next, acc2)
            end if
         end do
         nsave = size(acc1)
         allocate(null(nsave,2))
         null(:,1) = acc1
         null(:,2) = acc2
         result%method = 'exact'
      else
         allocate(null(simulations,2))
         null = 0.0_dp
         do i = 1, simulations
            do b = 1, size(blocks)
               work = blocks(b)%x
               call shuffle_real(work)
               call ad_statistic(work, blocks(b)%ns, st)
               null(i,:) = null(i,:) + st
            end do
         end do
         result%method = 'simulated'
      end if
      result%randomization_p(1) = real(count(null(:,1) >= result%statistic(1)),dp)/real(size(null,1),dp)
      result%randomization_p(2) = real(count(null(:,2) >= result%statistic(2)),dp)/real(size(null,1),dp)
      if (keep) then
         allocate(result%null_dist(size(null,1),2))
         result%null_dist = null
      end if
   end subroutine ad_test_combined

end module ksamples_ad
