! SPDX-License-Identifier: GPL-2.0-or-later
! Computational translation derived from R package spdep 1.4-2.
! Translated/modified 2026-08-30 for modern free-form Fortran.
module spdep_statistics
   use spdep_kinds, only : dp
   use spdep_types, only : neighbor_list, spatial_weights, spatial_test_result, &
      local_stat_result, eb_result, weights_constants
   use spdep_math, only : mean_dp, normal_two_sided_p, chi_square_sf, poisson_cdf, &
      deterministic_seed, random_permutation, safe_nan
   use spdep_graph, only : include_self
   use spdep_weights, only : nb2listw, lag_listw, listw2mat, spweights_constants
   implicit none
   private

   public :: moran
   public :: moran_test
   public :: moran_mc
   public :: moran_bv
   public :: local_moran
   public :: local_moran_bv
   public :: geary
   public :: geary_test
   public :: geary_mc
   public :: global_g_test
   public :: local_g
   public :: lee
   public :: local_lee
   public :: losh
   public :: joincount_test
   public :: ebest
   public :: eblocal
   public :: choynowski
   public :: ebi_moran

contains

   pure real(dp) function moran(x, listw) result(value)
      real(dp), intent(in) :: x(:) !! Numeric observations with one value per spatial region.
      type(spatial_weights), intent(in) :: listw !! Spatial weights defining the lag used by Moran's I.
      real(dp), allocatable :: z(:)
      real(dp), allocatable :: lz(:)
      real(dp) :: zz
      real(dp) :: s0
      integer :: n
      integer :: i

      n = listw%size()
      if (size(x) /= n .or. n < 2) then
         value = safe_nan()
         return
      end if
      allocate(z(n))
      z = x - mean_dp(x)
      zz = sum(z ** 2)
      if (zz <= 0.0_dp) then
         value = safe_nan()
         return
      end if
      s0 = 0.0_dp
      do i = 1, n
         s0 = s0 + sum(listw%weights(i)%values)
      end do
      if (s0 == 0.0_dp) then
         value = safe_nan()
         return
      end if
      lz = lag_listw(listw, z)
      value = (real(n, dp) / s0) * dot_product(z, lz) / zz
   end function moran

   pure function moran_test(x, listw, randomization) result(res)
      real(dp), intent(in) :: x(:) !! Numeric observations with one value per spatial region.
      type(spatial_weights), intent(in) :: listw !! Spatial weights used to calculate Moran's I and its moments.
      logical, intent(in), optional :: randomization !! If true, use randomization moments; otherwise use normality moments.
      type(spatial_test_result) :: res
      type(weights_constants) :: c
      logical :: use_randomization
      integer :: n
      real(dp) :: nn
      real(dp) :: k
      real(dp) :: vi
      real(dp) :: tmp

      n = listw%size()
      use_randomization = .true.
      if (present(randomization)) use_randomization = randomization
      res%statistic = moran(x, listw)
      k = centered_kurtosis_ratio(x)
      res%kurtosis = k
      if (n < 4) then
         res%expectation = safe_nan()
         res%variance = safe_nan()
         res%z_score = safe_nan()
         res%p_value = safe_nan()
         return
      end if
      c = spweights_constants(listw)
      if (c%s0 == 0.0_dp) then
         res%expectation = safe_nan()
         res%variance = safe_nan()
         res%z_score = safe_nan()
         res%p_value = safe_nan()
         return
      end if
      nn = real(n, dp)
      res%expectation = -1.0_dp / real(n - 1, dp)
      if (use_randomization) then
         vi = nn * (c%s1 * (nn * nn - 3.0_dp * nn + 3.0_dp) - nn * c%s2 &
            + 3.0_dp * c%s0 ** 2)
         tmp = k * (c%s1 * (nn * nn - nn) - 2.0_dp * nn * c%s2 &
            + 6.0_dp * c%s0 ** 2)
         vi = (vi - tmp) / (real((n - 1) * (n - 2) * (n - 3), dp) * c%s0 ** 2)
         vi = vi - res%expectation ** 2
      else
         vi = (nn * nn * c%s1 - nn * c%s2 + 3.0_dp * c%s0 ** 2) &
            / (c%s0 ** 2 * (nn * nn - 1.0_dp)) - res%expectation ** 2
      end if
      res%variance = vi
      if (vi > 0.0_dp) then
         res%z_score = (res%statistic - res%expectation) / sqrt(vi)
         res%p_value = normal_two_sided_p(res%z_score)
      else
         res%z_score = safe_nan()
         res%p_value = safe_nan()
      end if
   end function moran_test

   function moran_mc(x, listw, nsim, seed) result(res)
      real(dp), intent(in) :: x(:) !! Numeric observations permuted across regions under the Monte Carlo null.
      type(spatial_weights), intent(in) :: listw !! Spatial weights held fixed during permutations.
      integer, intent(in) :: nsim !! Number of random permutations; must be positive.
      integer, intent(in), optional :: seed !! Deterministic random seed; default is 12345.
      type(spatial_test_result) :: res
      integer, allocatable :: p(:)
      real(dp), allocatable :: xp(:)
      real(dp), allocatable :: sim(:)
      integer :: use_seed
      integer :: i
      integer :: more_extreme
      real(dp) :: center

      if (nsim <= 0 .or. size(x) /= listw%size()) then
         res%statistic = safe_nan()
         res%expectation = safe_nan()
         res%variance = safe_nan()
         res%z_score = safe_nan()
         res%p_value = safe_nan()
         return
      end if
      use_seed = 12345
      if (present(seed)) use_seed = seed
      call deterministic_seed(use_seed)
      allocate(p(size(x)), xp(size(x)), sim(nsim))
      res%statistic = moran(x, listw)
      do i = 1, nsim
         call random_permutation(p)
         xp = x(p)
         sim(i) = moran(xp, listw)
      end do
      center = sum(sim) / real(nsim, dp)
      res%expectation = center
      if (nsim > 1) then
         res%variance = sum((sim - center) ** 2) / real(nsim - 1, dp)
      else
         res%variance = 0.0_dp
      end if
      if (res%variance > 0.0_dp) then
         res%z_score = (res%statistic - center) / sqrt(res%variance)
      else
         res%z_score = safe_nan()
      end if
      more_extreme = count(abs(sim - center) >= abs(res%statistic - center))
      res%p_value = real(more_extreme + 1, dp) / real(nsim + 1, dp)
   end function moran_mc

   pure real(dp) function moran_bv(x, y, listw) result(value)
      real(dp), intent(in) :: x(:) !! First numeric variable with one observation per region.
      real(dp), intent(in) :: y(:) !! Second numeric variable whose spatial lag is cross-correlated with x.
      type(spatial_weights), intent(in) :: listw !! Spatial weights defining the lag of y.
      real(dp), allocatable :: zx(:)
      real(dp), allocatable :: zy(:)
      real(dp), allocatable :: lzy(:)
      real(dp) :: s0
      integer :: i
      integer :: n

      n = listw%size()
      if (size(x) /= n .or. size(y) /= n .or. n < 2) then
         value = safe_nan()
         return
      end if
      zx = x - mean_dp(x)
      zy = y - mean_dp(y)
      if (sum(zx ** 2) <= 0.0_dp) then
         value = safe_nan()
         return
      end if
      lzy = lag_listw(listw, zy)
      s0 = 0.0_dp
      do i = 1, n
         s0 = s0 + sum(listw%weights(i)%values)
      end do
      if (s0 == 0.0_dp) then
         value = safe_nan()
      else
         value = real(n, dp) * dot_product(zx, lzy) / (s0 * sum(zx ** 2))
      end if
   end function moran_bv

   pure function local_moran(x, listw, conditional, mlvar) result(res)
      real(dp), intent(in) :: x(:) !! Numeric observations with one value per region.
      type(spatial_weights), intent(in) :: listw !! Spatial weights excluding self-neighbors for local Moran statistics.
      logical, intent(in), optional :: conditional !! If true, use Sokal conditional moments; default is true.
      logical, intent(in), optional :: mlvar !! If true, use divisor n in the local-I scale; default is true.
      type(local_stat_result) :: res
      logical :: use_conditional
      logical :: use_mlvar
      real(dp), allocatable :: z(:)
      real(dp), allocatable :: lz(:)
      real(dp), allocatable :: wi(:)
      real(dp), allocatable :: wi2(:)
      real(dp) :: s2
      real(dp) :: m2
      real(dp) :: b2
      real(dp) :: aa
      real(dp) :: bb
      integer :: i
      integer :: n

      n = listw%size()
      allocate(res%statistic(n), res%expectation(n), res%variance(n), &
         res%z_score(n), res%p_value(n), z(n), wi(n), wi2(n))
      if (size(x) /= n .or. n < 3) then
         res%statistic = safe_nan()
         res%expectation = safe_nan()
         res%variance = safe_nan()
         res%z_score = safe_nan()
         res%p_value = safe_nan()
         return
      end if
      use_conditional = .true.
      if (present(conditional)) use_conditional = conditional
      use_mlvar = .true.
      if (present(mlvar)) use_mlvar = mlvar
      z = x - mean_dp(x)
      m2 = sum(z ** 2) / real(n, dp)
      if (use_mlvar) then
         s2 = m2
         b2 = (sum(z ** 4) / real(n, dp)) / (s2 ** 2)
      else
         s2 = sum(z ** 2) / real(n - 1, dp)
         b2 = (sum(z ** 4) / real(n - 1, dp)) / (s2 ** 2)
      end if
      lz = lag_listw(listw, z)
      res%statistic = (z / s2) * lz
      do i = 1, n
         wi(i) = sum(listw%weights(i)%values)
         wi2(i) = sum(listw%weights(i)%values ** 2)
      end do
      if (use_conditional) then
         res%expectation = -(z ** 2 * wi) / (real(n - 1, dp) * m2)
      else
         res%expectation = -wi / real(n - 1, dp)
      end if
      aa = (real(n, dp) - b2) / real(n - 1, dp)
      bb = (2.0_dp * b2 - real(n, dp)) / real((n - 1) * (n - 2), dp)
      if (use_conditional) then
         res%variance = (z / m2) ** 2 * (real(n, dp) / real(n - 2, dp)) &
            * (wi2 - wi ** 2 / real(n - 1, dp)) &
            * (m2 - z ** 2 / real(n - 1, dp))
      else
         res%variance = aa * wi2 + bb * (wi ** 2 - wi2) - res%expectation ** 2
      end if
      do i = 1, n
         if (res%variance(i) > 0.0_dp) then
            res%z_score(i) = (res%statistic(i) - res%expectation(i)) / sqrt(res%variance(i))
            res%p_value(i) = normal_two_sided_p(res%z_score(i))
         else
            res%z_score(i) = safe_nan()
            res%p_value(i) = safe_nan()
         end if
      end do
   end function local_moran

   pure function local_moran_bv(x, y, listw) result(statistic)
      real(dp), intent(in) :: x(:) !! First variable centered and scaled in the local bivariate Moran product.
      real(dp), intent(in) :: y(:) !! Second variable whose centered values are spatially lagged.
      type(spatial_weights), intent(in) :: listw !! Spatial weights defining the lag of y.
      real(dp), allocatable :: statistic(:)
      real(dp), allocatable :: zx(:)
      real(dp), allocatable :: zy(:)
      real(dp), allocatable :: lzy(:)
      real(dp) :: sx
      integer :: n

      n = listw%size()
      allocate(statistic(n))
      if (size(x) /= n .or. size(y) /= n .or. n < 2) then
         statistic = safe_nan()
         return
      end if
      zx = x - mean_dp(x)
      zy = y - mean_dp(y)
      sx = sqrt(sum(zx ** 2) / real(n, dp))
      if (sx <= 0.0_dp) then
         statistic = safe_nan()
         return
      end if
      lzy = lag_listw(listw, zy)
      statistic = (zx / (sx * sx)) * lzy
   end function local_moran_bv

   pure real(dp) function geary(x, listw) result(value)
      real(dp), intent(in) :: x(:) !! Numeric observations with one value per region.
      type(spatial_weights), intent(in) :: listw !! Spatial weights used in pairwise squared differences.
      real(dp), allocatable :: z(:)
      real(dp) :: zz
      real(dp) :: numerator
      real(dp) :: s0
      integer :: i
      integer :: j
      integer :: k
      integer :: n

      n = listw%size()
      if (size(x) /= n .or. n < 2) then
         value = safe_nan()
         return
      end if
      z = x - mean_dp(x)
      zz = sum(z ** 2)
      if (zz <= 0.0_dp) then
         value = safe_nan()
         return
      end if
      numerator = 0.0_dp
      s0 = 0.0_dp
      do i = 1, n
         do k = 1, size(listw%nb%neighbors(i)%values)
            j = listw%nb%neighbors(i)%values(k)
            numerator = numerator + listw%weights(i)%values(k) * (z(i) - z(j)) ** 2
            s0 = s0 + listw%weights(i)%values(k)
         end do
      end do
      if (s0 == 0.0_dp) then
         value = safe_nan()
      else
         value = real(n - 1, dp) * numerator / (2.0_dp * s0 * zz)
      end if
   end function geary

   pure function geary_test(x, listw, randomization) result(res)
      real(dp), intent(in) :: x(:) !! Numeric observations with one value per region.
      type(spatial_weights), intent(in) :: listw !! Spatial weights used to calculate Geary's C and its moments.
      logical, intent(in), optional :: randomization !! If true, use randomization moments; otherwise use normality moments.
      type(spatial_test_result) :: res
      type(weights_constants) :: c
      logical :: use_randomization
      integer :: n
      real(dp) :: nn
      real(dp) :: k
      real(dp) :: vc

      n = listw%size()
      use_randomization = .true.
      if (present(randomization)) use_randomization = randomization
      res%statistic = geary(x, listw)
      k = centered_kurtosis_ratio(x)
      res%kurtosis = k
      res%expectation = 1.0_dp
      if (n < 4) then
         res%variance = safe_nan()
         res%z_score = safe_nan()
         res%p_value = safe_nan()
         return
      end if
      c = spweights_constants(listw)
      nn = real(n, dp)
      if (c%s0 == 0.0_dp) then
         res%variance = safe_nan()
         res%z_score = safe_nan()
         res%p_value = safe_nan()
         return
      end if
      if (use_randomization) then
         vc = real(n - 1, dp) * c%s1 * (nn * nn - 3.0_dp * nn + 3.0_dp - k * real(n - 1, dp))
         vc = vc - 0.25_dp * real(n - 1, dp) * c%s2 &
            * (nn * nn + 3.0_dp * nn - 6.0_dp - k * (nn * nn - nn + 2.0_dp))
         vc = vc + c%s0 ** 2 * (nn * nn - 3.0_dp - k * real((n - 1) ** 2, dp))
         vc = vc / (nn * real((n - 2) * (n - 3), dp) * c%s0 ** 2)
      else
         vc = ((2.0_dp * c%s1 + c%s2) * real(n - 1, dp) - 4.0_dp * c%s0 ** 2) &
            / (2.0_dp * real(n + 1, dp) * c%s0 ** 2)
      end if
      res%variance = vc
      if (vc > 0.0_dp) then
         res%z_score = (res%expectation - res%statistic) / sqrt(vc)
         res%p_value = normal_two_sided_p(res%z_score)
      else
         res%z_score = safe_nan()
         res%p_value = safe_nan()
      end if
   end function geary_test

   function geary_mc(x, listw, nsim, seed) result(res)
      real(dp), intent(in) :: x(:) !! Numeric observations permuted across regions under the Monte Carlo null.
      type(spatial_weights), intent(in) :: listw !! Spatial weights held fixed during permutations.
      integer, intent(in) :: nsim !! Number of random permutations; must be positive.
      integer, intent(in), optional :: seed !! Deterministic random seed; default is 12345.
      type(spatial_test_result) :: res
      integer, allocatable :: p(:)
      real(dp), allocatable :: xp(:)
      real(dp), allocatable :: sim(:)
      integer :: use_seed
      integer :: i
      integer :: more_extreme
      real(dp) :: center

      if (nsim <= 0 .or. size(x) /= listw%size()) then
         res%statistic = safe_nan()
         res%expectation = safe_nan()
         res%variance = safe_nan()
         res%z_score = safe_nan()
         res%p_value = safe_nan()
         return
      end if
      use_seed = 12345
      if (present(seed)) use_seed = seed
      call deterministic_seed(use_seed)
      allocate(p(size(x)), xp(size(x)), sim(nsim))
      res%statistic = geary(x, listw)
      do i = 1, nsim
         call random_permutation(p)
         xp = x(p)
         sim(i) = geary(xp, listw)
      end do
      center = sum(sim) / real(nsim, dp)
      res%expectation = center
      if (nsim > 1) then
         res%variance = sum((sim - center) ** 2) / real(nsim - 1, dp)
      else
         res%variance = 0.0_dp
      end if
      if (res%variance > 0.0_dp) then
         res%z_score = (center - res%statistic) / sqrt(res%variance)
      else
         res%z_score = safe_nan()
      end if
      more_extreme = count(abs(sim - center) >= abs(res%statistic - center))
      res%p_value = real(more_extreme + 1, dp) / real(nsim + 1, dp)
   end function geary_mc

   pure function global_g_test(x, listw, b1_correct) result(res)
      real(dp), intent(in) :: x(:) !! Nonnegative numeric observations used in Getis-Ord global G.
      type(spatial_weights), intent(in) :: listw !! Spatial weights used by global G and its randomization moments.
      logical, intent(in), optional :: b1_correct !! If true, use the corrected 6*S0^2 B1 term; default is true.
      type(spatial_test_result) :: res
      type(weights_constants) :: c
      logical :: correct_b1
      real(dp), allocatable :: lx(:)
      real(dp) :: sx
      real(dp) :: sx2
      real(dp) :: sx3
      real(dp) :: sx4
      real(dp) :: denom
      real(dp) :: b0
      real(dp) :: b1
      real(dp) :: b2
      real(dp) :: b3
      real(dp) :: b4
      real(dp) :: nn
      real(dp) :: numerator_var
      integer :: n

      n = listw%size()
      correct_b1 = .true.
      if (present(b1_correct)) correct_b1 = b1_correct
      if (size(x) /= n .or. n < 4 .or. any(x < 0.0_dp)) then
         res%statistic = safe_nan()
         res%expectation = safe_nan()
         res%variance = safe_nan()
         res%z_score = safe_nan()
         res%p_value = safe_nan()
         return
      end if
      lx = lag_listw(listw, x)
      sx = sum(x)
      sx2 = sum(x ** 2)
      sx3 = sum(x ** 3)
      sx4 = sum(x ** 4)
      denom = sx * sx - sx2
      if (denom == 0.0_dp) then
         res%statistic = safe_nan()
         res%expectation = safe_nan()
         res%variance = safe_nan()
         res%z_score = safe_nan()
         res%p_value = safe_nan()
         return
      end if
      res%statistic = dot_product(x, lx) / denom
      c = spweights_constants(listw)
      nn = real(n, dp)
      res%expectation = c%s0 / (nn * real(n - 1, dp))
      b0 = (nn * nn - 3.0_dp * nn + 3.0_dp) * c%s1 - nn * c%s2 + 3.0_dp * c%s0 ** 2
      if (correct_b1) then
         b1 = -((nn * nn - nn) * c%s1 - 2.0_dp * nn * c%s2 + 6.0_dp * c%s0 ** 2)
      else
         b1 = -((nn * nn - nn) * c%s1 - 2.0_dp * nn * c%s2 + 3.0_dp * c%s0 ** 2)
      end if
      b2 = -(2.0_dp * nn * c%s1 - (nn + 3.0_dp) * c%s2 + 6.0_dp * c%s0 ** 2)
      b3 = 4.0_dp * real(n - 1, dp) * c%s1 - 2.0_dp * real(n + 1, dp) * c%s2 &
         + 8.0_dp * c%s0 ** 2
      b4 = c%s1 - c%s2 + c%s0 ** 2
      numerator_var = b0 * sx2 ** 2 + b1 * sx4 + b2 * sx ** 2 * sx2 &
         + b3 * sx * sx3 + b4 * sx ** 4
      res%variance = numerator_var / (denom ** 2 * nn * real((n - 1) * (n - 2) * (n - 3), dp)) &
         - res%expectation ** 2
      if (res%variance > 0.0_dp) then
         res%z_score = (res%statistic - res%expectation) / sqrt(res%variance)
         res%p_value = normal_two_sided_p(res%z_score)
      else
         res%z_score = safe_nan()
         res%p_value = safe_nan()
      end if
   end function global_g_test

   pure function local_g(x, listw) result(res)
      real(dp), intent(in) :: x(:) !! Numeric observations with one value per region.
      type(spatial_weights), intent(in) :: listw !! Spatial weights; self-inclusion selects G-star rather than G.
      type(local_stat_result) :: res
      real(dp), allocatable :: lx(:)
      real(dp), allocatable :: wi(:)
      real(dp), allocatable :: s1i(:)
      real(dp), allocatable :: xibar(:)
      real(dp), allocatable :: si2(:)
      real(dp) :: xstar
      real(dp) :: mu
      integer :: i
      integer :: n
      logical :: gstar

      n = listw%size()
      allocate(res%statistic(n), res%expectation(n), res%variance(n), &
         res%z_score(n), res%p_value(n), wi(n), s1i(n), xibar(n), si2(n))
      if (size(x) /= n .or. n < 3) then
         res%statistic = safe_nan()
         res%expectation = safe_nan()
         res%variance = safe_nan()
         res%z_score = safe_nan()
         res%p_value = safe_nan()
         return
      end if
      gstar = listw%nb%self_included
      lx = lag_listw(listw, x)
      xstar = sum(x)
      do i = 1, n
         wi(i) = sum(listw%weights(i)%values)
         s1i(i) = sum(listw%weights(i)%values ** 2)
      end do
      if (gstar) then
         mu = mean_dp(x)
         xibar = mu
         si2 = sum((x - mu) ** 2) / real(n, dp)
         res%expectation = wi * xibar
         res%variance = si2 * (real(n, dp) * s1i - wi ** 2) / real(n - 1, dp)
         do i = 1, n
            if (xstar /= 0.0_dp) then
               res%statistic(i) = lx(i) / xstar
               res%expectation(i) = res%expectation(i) / xstar
               res%variance(i) = res%variance(i) / (xstar ** 2)
            else
               res%statistic(i) = safe_nan()
               res%expectation(i) = safe_nan()
               res%variance(i) = safe_nan()
            end if
         end do
      else
         xibar = (xstar - x) / real(n - 1, dp)
         si2 = (sum(x ** 2) - x ** 2) / real(n - 1, dp) - xibar ** 2
         res%expectation = wi * xibar
         res%variance = si2 * (real(n - 1, dp) * s1i - wi ** 2) / real(n - 2, dp)
         do i = 1, n
            if (xstar - x(i) /= 0.0_dp) then
               res%statistic(i) = lx(i) / (xstar - x(i))
               res%expectation(i) = res%expectation(i) / (xstar - x(i))
               res%variance(i) = res%variance(i) / ((xstar - x(i)) ** 2)
            else
               res%statistic(i) = safe_nan()
               res%expectation(i) = safe_nan()
               res%variance(i) = safe_nan()
            end if
         end do
      end if
      do i = 1, n
         if (res%variance(i) > 0.0_dp) then
            res%z_score(i) = (res%statistic(i) - res%expectation(i)) / sqrt(res%variance(i))
            res%p_value(i) = normal_two_sided_p(res%z_score(i))
         else
            res%z_score(i) = safe_nan()
            res%p_value(i) = safe_nan()
         end if
      end do
   end function local_g

   pure real(dp) function lee(x, y, listw) result(value)
      real(dp), intent(in) :: x(:) !! First numeric variable used in Lee's global spatial association statistic.
      real(dp), intent(in) :: y(:) !! Second numeric variable used in Lee's global spatial association statistic.
      type(spatial_weights), intent(in) :: listw !! Spatial weights applied to both centered variables.
      real(dp), allocatable :: zx(:)
      real(dp), allocatable :: zy(:)
      real(dp), allocatable :: lzx(:)
      real(dp), allocatable :: lzy(:)
      real(dp) :: s2
      real(dp) :: denom
      real(dp) :: wi
      integer :: i
      integer :: n

      n = listw%size()
      if (size(x) /= n .or. size(y) /= n .or. n < 2) then
         value = safe_nan()
         return
      end if
      zx = x - mean_dp(x)
      zy = y - mean_dp(y)
      denom = sqrt(sum(zx ** 2)) * sqrt(sum(zy ** 2))
      if (denom <= 0.0_dp) then
         value = safe_nan()
         return
      end if
      lzx = lag_listw(listw, zx)
      lzy = lag_listw(listw, zy)
      s2 = 0.0_dp
      do i = 1, n
         wi = sum(listw%weights(i)%values)
         s2 = s2 + wi * wi
      end do
      if (s2 == 0.0_dp) then
         value = safe_nan()
      else
         value = real(n, dp) * dot_product(lzx, lzy) / (s2 * denom)
      end if
   end function lee

   pure function local_lee(x, y, listw) result(value)
      real(dp), intent(in) :: x(:) !! First numeric variable used in Lee's local spatial association statistic.
      real(dp), intent(in) :: y(:) !! Second numeric variable used in Lee's local spatial association statistic.
      type(spatial_weights), intent(in) :: listw !! Spatial weights applied to both centered variables.
      real(dp), allocatable :: value(:)
      real(dp), allocatable :: zx(:)
      real(dp), allocatable :: zy(:)
      real(dp), allocatable :: lzx(:)
      real(dp), allocatable :: lzy(:)
      real(dp) :: denom
      integer :: n

      n = listw%size()
      allocate(value(n))
      if (size(x) /= n .or. size(y) /= n .or. n < 2) then
         value = safe_nan()
         return
      end if
      zx = x - mean_dp(x)
      zy = y - mean_dp(y)
      denom = sqrt(sum(zx ** 2)) * sqrt(sum(zy ** 2))
      if (denom <= 0.0_dp) then
         value = safe_nan()
         return
      end if
      lzx = lag_listw(listw, zx)
      lzy = lag_listw(listw, zy)
      value = real(n, dp) * lzx * lzy / denom
   end function local_lee

   pure function losh(x, listw, a) result(res)
      real(dp), intent(in) :: x(:) !! Numeric observations used by the local spatial heteroscedasticity statistic.
      type(spatial_weights), intent(in) :: listw !! Spatial weights used for local means and dispersion lags.
      real(dp), intent(in), optional :: a !! Positive absolute-deviation exponent; default 2 enables chi-square inference.
      type(local_stat_result) :: res
      real(dp), allocatable :: wi(:)
      real(dp), allocatable :: wi2(:)
      real(dp), allocatable :: xbar(:)
      real(dp), allocatable :: ei(:)
      real(dp), allocatable :: denom(:)
      real(dp), allocatable :: lei(:)
      real(dp) :: aa
      real(dp) :: mean_ei
      real(dp) :: var_ei
      integer :: i
      integer :: n

      n = listw%size()
      allocate(res%statistic(n), res%expectation(n), res%variance(n), &
         res%z_score(n), res%p_value(n), wi(n), wi2(n), xbar(n), ei(n), denom(n))
      if (size(x) /= n .or. n < 2) then
         res%statistic = safe_nan()
         res%expectation = safe_nan()
         res%variance = safe_nan()
         res%z_score = safe_nan()
         res%p_value = safe_nan()
         return
      end if
      aa = 2.0_dp
      if (present(a)) aa = a
      if (aa <= 0.0_dp) then
         res%statistic = safe_nan()
         res%expectation = safe_nan()
         res%variance = safe_nan()
         res%z_score = safe_nan()
         res%p_value = safe_nan()
         return
      end if
      do i = 1, n
         wi(i) = sum(listw%weights(i)%values)
         wi2(i) = sum(listw%weights(i)%values ** 2)
      end do
      xbar = lag_listw(listw, x)
      do i = 1, n
         if (wi(i) /= 0.0_dp) then
            xbar(i) = xbar(i) / wi(i)
         else
            xbar(i) = safe_nan()
         end if
      end do
      ei = abs(x - xbar) ** aa
      mean_ei = sum(ei) / real(n, dp)
      denom = mean_ei * wi
      lei = lag_listw(listw, ei)
      do i = 1, n
         if (denom(i) /= 0.0_dp) then
            res%statistic(i) = lei(i) / denom(i)
         else
            res%statistic(i) = safe_nan()
         end if
      end do
      var_ei = sum(ei ** 2) / real(n, dp) - mean_ei ** 2
      res%expectation = 1.0_dp
      do i = 1, n
         if (denom(i) /= 0.0_dp) then
            res%variance(i) = var_ei * (real(n, dp) * wi2(i) - wi(i) ** 2) &
               / (real(n - 1, dp) * denom(i) ** 2)
         else
            res%variance(i) = safe_nan()
         end if
         if (res%variance(i) > 0.0_dp) then
            res%z_score(i) = 2.0_dp * res%statistic(i) / res%variance(i)
            if (abs(aa - 2.0_dp) <= epsilon(1.0_dp)) then
               res%p_value(i) = chi_square_sf(res%z_score(i), 2.0_dp / res%variance(i))
            else
               res%p_value(i) = safe_nan()
            end if
         else
            res%z_score(i) = safe_nan()
            res%p_value(i) = safe_nan()
         end if
      end do
   end function losh

   pure function joincount_test(x, listw, free_sampling) result(res)
      integer, intent(in) :: x(:) !! Binary indicator vector; entries must be zero or one.
      type(spatial_weights), intent(in) :: listw !! Spatial weights used to count same-class one-one joins.
      logical, intent(in), optional :: free_sampling !! If true, use Bernoulli/free moments rather than fixed-count moments.
      type(spatial_test_result) :: res
      type(weights_constants) :: c
      real(dp), allocatable :: xr(:)
      real(dp), allocatable :: lx(:)
      logical :: use_free
      integer :: n
      integer :: m
      real(dp) :: p
      real(dp) :: nn
      real(dp) :: mm

      n = listw%size()
      use_free = .false.
      if (present(free_sampling)) use_free = free_sampling
      if (size(x) /= n .or. any(x < 0) .or. any(x > 1) .or. n < 4) then
         res%statistic = safe_nan()
         res%expectation = safe_nan()
         res%variance = safe_nan()
         res%z_score = safe_nan()
         res%p_value = safe_nan()
         return
      end if
      allocate(xr(n))
      xr = real(x, dp)
      lx = lag_listw(listw, xr)
      res%statistic = 0.5_dp * dot_product(xr, lx)
      c = spweights_constants(listw)
      m = sum(x)
      nn = real(n, dp)
      mm = real(m, dp)
      if (use_free) then
         p = mm / nn
         res%expectation = 0.5_dp * c%s0 * p ** 2
         res%variance = 0.25_dp * (c%s1 * p ** 2 + (c%s2 - 2.0_dp * c%s1) * p ** 3 &
            + (c%s1 - c%s2 + c%s0 ** 2) * p ** 4)
      else
         res%expectation = c%s0 * mm * real(m - 1, dp) / (2.0_dp * nn * real(n - 1, dp))
         res%variance = c%s1 * mm * real(m - 1, dp) / (nn * real(n - 1, dp))
         res%variance = res%variance + (c%s2 - 2.0_dp * c%s1) * mm * real((m - 1) * (m - 2), dp) &
            / (nn * real((n - 1) * (n - 2), dp))
         res%variance = res%variance + (c%s0 ** 2 + c%s1 - c%s2) * mm &
            * real((m - 1) * (m - 2) * (m - 3), dp) &
            / (nn * real((n - 1) * (n - 2) * (n - 3), dp))
         res%variance = 0.25_dp * res%variance - res%expectation ** 2
      end if
      if (res%variance > 0.0_dp) then
         res%z_score = (res%statistic - res%expectation) / sqrt(res%variance)
         res%p_value = normal_two_sided_p(res%z_score)
      else
         res%z_score = safe_nan()
         res%p_value = safe_nan()
      end if
   end function joincount_test

   pure function ebest(cases, population, binomial) result(res)
      real(dp), intent(in) :: cases(:) !! Nonnegative observed case counts, one per region.
      real(dp), intent(in) :: population(:) !! Strictly positive population at risk, one per region.
      logical, intent(in), optional :: binomial !! If true, use the upstream binomial shrinkage formula; default is Poisson.
      type(eb_result) :: res
      logical :: use_binomial
      real(dp), allocatable :: rho(:)
      real(dp) :: s2
      real(dp) :: b
      real(dp) :: a
      real(dp) :: mean_pop
      real(dp) :: pop_sum
      integer :: n

      n = size(cases)
      allocate(res%raw(n), res%estimate(n))
      use_binomial = .false.
      if (present(binomial)) use_binomial = binomial
      if (size(population) /= n .or. n == 0 .or. any(population <= 0.0_dp) .or. any(cases < 0.0_dp)) then
         res%raw = safe_nan()
         res%estimate = safe_nan()
         res%global_rate = safe_nan()
         res%dispersion = safe_nan()
         return
      end if
      res%raw = cases / population
      pop_sum = sum(population)
      b = sum(cases) / pop_sum
      s2 = sum(population * (res%raw - b) ** 2) / pop_sum
      res%global_rate = b
      if (use_binomial) then
         mean_pop = mean_dp(population)
         allocate(rho(n))
         rho = (population * s2 - (population / mean_pop) * b * (1.0_dp - b)) &
            / ((population - 1.0_dp) * s2 + ((mean_pop - population) / mean_pop) * b * (1.0_dp - b))
         res%estimate = rho * res%raw + (1.0_dp - rho) * b
         res%dispersion = s2
      else
         a = s2 - b / (pop_sum / real(n, dp))
         a = max(0.0_dp, a)
         res%estimate = b + a * (res%raw - b) / (a + b / population)
         res%dispersion = a
      end if
   end function ebest

   pure function eblocal(cases, population, nb) result(res)
      real(dp), intent(in) :: cases(:) !! Nonnegative observed case counts, one per region.
      real(dp), intent(in) :: population(:) !! Strictly positive population at risk, one per region.
      type(neighbor_list), intent(in) :: nb !! Neighbor graph defining each local empirical-Bayes borrowing set.
      type(eb_result) :: res
      type(neighbor_list) :: self_nb
      type(spatial_weights) :: binary_w
      type(spatial_weights) :: row_w
      real(dp), allocatable :: ri(:)
      real(dp), allocatable :: ni(:)
      real(dp), allocatable :: nbar(:)
      real(dp), allocatable :: mi(:)
      real(dp), allocatable :: ci(:)
      real(dp), allocatable :: ai(:)
      integer :: n

      n = nb%size()
      allocate(res%raw(n), res%estimate(n))
      if (size(cases) /= n .or. size(population) /= n .or. any(population <= 0.0_dp) &
         .or. any(cases < 0.0_dp)) then
         res%raw = safe_nan()
         res%estimate = safe_nan()
         res%global_rate = safe_nan()
         res%dispersion = safe_nan()
         return
      end if
      self_nb = include_self(nb)
      binary_w = nb2listw(self_nb, "B")
      row_w = nb2listw(self_nb, "W")
      res%raw = cases / population
      ri = lag_listw(binary_w, cases)
      ni = lag_listw(binary_w, population)
      nbar = lag_listw(row_w, population)
      allocate(mi(n), ci(n), ai(n))
      mi = ri / ni
      ci = lag_listw(binary_w, population * (res%raw - mi) ** 2)
      ai = ci / ni - mi / nbar
      ai = max(ai, 0.0_dp)
      res%estimate = mi + (res%raw - mi) * ai / (ai + mi / population)
      res%global_rate = sum(cases) / sum(population)
      res%dispersion = sum(ai) / real(max(1, n), dp)
   end function eblocal

   pure function choynowski(cases, population) result(probability)
      integer, intent(in) :: cases(:) !! Nonnegative integer case counts, one per region.
      real(dp), intent(in) :: population(:) !! Positive population at risk, one per region.
      real(dp), allocatable :: probability(:)
      real(dp) :: global_rate
      real(dp) :: expected
      integer :: i
      integer :: n

      n = size(cases)
      allocate(probability(n))
      if (size(population) /= n .or. any(population <= 0.0_dp) .or. any(cases < 0)) then
         probability = safe_nan()
         return
      end if
      global_rate = real(sum(cases), dp) / sum(population)
      if (global_rate > 1.0_dp) then
         probability = safe_nan()
         return
      end if
      do i = 1, n
         expected = population(i) * global_rate
         if (real(cases(i), dp) < expected) then
            probability(i) = poisson_cdf(cases(i), expected)
         else
            probability(i) = 1.0_dp - poisson_cdf(cases(i) - 1, expected)
         end if
      end do
   end function choynowski

   pure real(dp) function ebi_moran(rate, listw, subtract_mean) result(value)
      real(dp), intent(in) :: rate(:) !! Empirical-Bayes adjusted or raw rates with one value per region.
      type(spatial_weights), intent(in) :: listw !! Spatial weights used in the rate autocorrelation statistic.
      logical, intent(in), optional :: subtract_mean !! If true, center both focal and lag variables; default is true.
      logical :: use_center
      real(dp), allocatable :: z(:)
      real(dp), allocatable :: lz(:)
      real(dp) :: zz
      real(dp) :: s0
      integer :: i
      integer :: n

      n = listw%size()
      use_center = .true.
      if (present(subtract_mean)) use_center = subtract_mean
      if (size(rate) /= n .or. n < 2) then
         value = safe_nan()
         return
      end if
      allocate(z(n))
      if (use_center) then
         z = rate - mean_dp(rate)
      else
         z = rate
      end if
      zz = sum(z ** 2)
      lz = lag_listw(listw, z)
      s0 = 0.0_dp
      do i = 1, n
         s0 = s0 + sum(listw%weights(i)%values)
      end do
      if (zz <= 0.0_dp .or. s0 == 0.0_dp) then
         value = safe_nan()
      else
         value = real(n, dp) * dot_product(z, lz) / (s0 * zz)
      end if
   end function ebi_moran

   pure real(dp) function centered_kurtosis_ratio(x) result(k)
      real(dp), intent(in) :: x(:) !! Numeric observations used to compute n*sum(z^4)/(sum(z^2)^2).
      real(dp), allocatable :: z(:)
      real(dp) :: zz

      if (size(x) == 0) then
         k = safe_nan()
         return
      end if
      z = x - mean_dp(x)
      zz = sum(z ** 2)
      if (zz <= 0.0_dp) then
         k = safe_nan()
      else
         k = real(size(x), dp) * sum(z ** 4) / (zz ** 2)
      end if
   end function centered_kurtosis_ratio

end module spdep_statistics
