! SPDX-License-Identifier: GPL-2.0-or-later
! Computational translation derived from R package spdep 1.4-2.
! Translated/modified 2026-08-30 for modern free-form Fortran.
module spdep_delta
   use spdep_kinds, only : dp
   use spdep_types, only : spatial_delta_result
   use spdep_math, only : matrix_trace, symmetric_max_eigenvalue, normal_cdf, &
      normal_two_sided_p, safe_nan
   implicit none
   private

   public :: spatialdelta
   public :: linearised_diffusive_weights
   public :: metropolis_hastings_weights
   public :: iterative_proportional_fitting_weights
   public :: graph_distance_weights
   public :: localdelta
   public :: cornish_fisher

contains

   pure function spatialdelta(dissimilarity_matrix, adjusted_spatial_weights, &
      regional_weights, alternative) result(res)
      real(dp), intent(in) :: dissimilarity_matrix(:, :) !! Square pairwise dissimilarity matrix for regional feature profiles.
      real(dp), intent(in) :: adjusted_spatial_weights(:, :) !! Square adjusted spatial transition-weight matrix.
      real(dp), intent(in) :: regional_weights(:) !! Strictly positive regional masses; internally normalized to sum to one.
      character(len=*), intent(in), optional :: alternative !! Tail alternative: greater, less, or two.sided; default greater.
      type(spatial_delta_result) :: res
      real(dp), allocatable :: rw(:)
      real(dp), allocatable :: h(:, :)
      real(dp), allocatable :: b(:, :)
      real(dp), allocatable :: kx(:, :)
      real(dp), allocatable :: kw(:, :)
      real(dp), allocatable :: kw2(:, :)
      real(dp), allocatable :: kx2(:, :)
      real(dp), allocatable :: tmp(:, :)
      real(dp), allocatable :: sf(:)
      real(dp) :: tr_kx
      real(dp) :: tr_kx2
      real(dp) :: tr_kw_kx
      real(dp) :: tr_kw2
      real(dp) :: tr_kw
      real(dp) :: tr_kw3
      real(dp) :: tr_kx3
      real(dp) :: tr_kw4
      real(dp) :: tr_kx4
      real(dp) :: mubar2
      real(dp) :: lbar2
      real(dp) :: mubar3
      real(dp) :: lbar3
      real(dp) :: mubar4
      real(dp) :: lbar4
      real(dp) :: vd1
      real(dp) :: vd2
      real(dp) :: vd3
      real(dp) :: vi
      real(dp) :: vx
      real(dp) :: skew1
      real(dp) :: alphamu
      real(dp) :: alphalambda
      real(dp) :: gammamu
      real(dp) :: gammalambda
      real(dp) :: kurt1
      real(dp) :: kurt2
      real(dp) :: kurt3
      real(dp) :: nn
      character(len=12) :: alt
      integer :: n
      integer :: i
      integer :: j

      n = size(regional_weights)
      alt = "greater"
      if (present(alternative)) alt = adjustl(alternative)
      if (n < 5 .or. size(dissimilarity_matrix, 1) /= n .or. &
         size(dissimilarity_matrix, 2) /= n .or. &
         size(adjusted_spatial_weights, 1) /= n .or. &
         size(adjusted_spatial_weights, 2) /= n .or. any(regional_weights <= 0.0_dp)) then
         call set_delta_nan(res)
         return
      end if
      allocate(rw(n), h(n, n), b(n, n), kx(n, n), kw(n, n), sf(n))
      rw = regional_weights / sum(regional_weights)
      h = 0.0_dp
      do i = 1, n
         h(i, i) = 1.0_dp
         do j = 1, n
            h(i, j) = h(i, j) - rw(j)
         end do
      end do
      b = -0.5_dp * matmul(matmul(h, dissimilarity_matrix), transpose(h))
      sf = sqrt(rw)
      do i = 1, n
         do j = 1, n
            kx(i, j) = sf(i) * b(i, j) * sf(j)
            kw(i, j) = sf(i) * adjusted_spatial_weights(i, j) / sf(j) - sf(i) * sf(j)
         end do
      end do
      kx2 = matmul(kx, kx)
      kw2 = matmul(kw, kw)
      tr_kx = matrix_trace(kx)
      tr_kx2 = matrix_trace(kx2)
      tr_kw_kx = matrix_trace(matmul(kw, kx))
      tr_kw2 = matrix_trace(kw2)
      tr_kw = matrix_trace(kw)
      tmp = matmul(kw2, kw)
      tr_kw3 = matrix_trace(tmp)
      tr_kw4 = matrix_trace(matmul(tmp, kw))
      tmp = matmul(kx2, kx)
      tr_kx3 = matrix_trace(tmp)
      tr_kx4 = matrix_trace(matmul(tmp, kx))
      nn = real(n, dp)
      mubar2 = tr_kw2 / real(n - 1, dp) - (tr_kw / real(n - 1, dp)) ** 2
      lbar2 = tr_kx2 / real(n - 1, dp) - (tr_kx / real(n - 1, dp)) ** 2
      mubar3 = tr_kw3 / real(n - 1, dp) &
         - 3.0_dp * tr_kw2 * tr_kw / real((n - 1) ** 2, dp) &
         + 2.0_dp * (tr_kw / real(n - 1, dp)) ** 3
      lbar3 = tr_kx3 / real(n - 1, dp) &
         - 3.0_dp * tr_kx2 * tr_kx / real((n - 1) ** 2, dp) &
         + 2.0_dp * (tr_kx / real(n - 1, dp)) ** 3
      mubar4 = tr_kw4 / real(n - 1, dp) &
         - 4.0_dp * tr_kw3 * tr_kw / real((n - 1) ** 2, dp) &
         + 6.0_dp * tr_kw2 * (tr_kw / real(n - 1, dp)) ** 2 / real(n - 1, dp) &
         - 3.0_dp * (tr_kw / real(n - 1, dp)) ** 4
      lbar4 = tr_kx4 / real(n - 1, dp) &
         - 4.0_dp * tr_kx3 * tr_kx / real((n - 1) ** 2, dp) &
         + 6.0_dp * tr_kx2 * (tr_kx / real(n - 1, dp)) ** 2 / real(n - 1, dp) &
         - 3.0_dp * (tr_kx / real(n - 1, dp)) ** 4
      if (tr_kx == 0.0_dp .or. tr_kw2 == 0.0_dp .or. tr_kx2 == 0.0_dp &
         .or. mubar2 <= 0.0_dp .or. lbar2 <= 0.0_dp) then
         call set_delta_nan(res)
         return
      end if
      res%delta = tr_kw_kx / tr_kx
      res%rv = res%delta * tr_kx / (tr_kw2 * tr_kx2)
      res%expectation = tr_kw / real(n - 1, dp)
      vd1 = 2.0_dp / (real(n - 2, dp) * real((n - 1) ** 2, dp) * real(n + 1, dp))
      vd2 = real(n - 1, dp) * tr_kw2 - tr_kw ** 2
      vd3 = (real(n - 1, dp) * tr_kx2 - tr_kx ** 2) / (tr_kx ** 2)
      res%variance = vd1 * vd2 * vd3
      vi = 2.0_dp / (nn ** 2 - 1.0_dp) &
         * (matrix_trace(matmul(adjusted_spatial_weights, adjusted_spatial_weights)) &
         - 1.0_dp - (matrix_trace(adjusted_spatial_weights) - 1.0_dp) ** 2 / real(n - 1, dp))
      vx = (real(n - 1, dp) / (tr_kx ** 2 / tr_kx2) - 1.0_dp) / real(n - 2, dp)
      res%variance_product = vi * vx
      alphamu = mubar3 / (mubar2 ** 1.5_dp)
      alphalambda = lbar3 / (lbar2 ** 1.5_dp)
      skew1 = sqrt(8.0_dp * real((n - 2) * (n + 1), dp)) &
         / real((n - 3) * (n + 3), dp)
      res%skewness = skew1 * alphamu * alphalambda
      gammamu = mubar4 / (mubar2 ** 2) - 3.0_dp
      gammalambda = lbar4 / (lbar2 ** 2) - 3.0_dp
      kurt1 = 3.0_dp * real((n - 2) * (n + 1), dp) &
         / real((n - 4) * (n - 3) * (n - 1) * n * (n + 3) * (n + 5), dp)
      kurt2 = 4.0_dp * real(n * n - n + 2, dp) * gammamu * gammalambda &
         + real(4 * n * n - 8 * n + 52, dp) * (gammamu + gammalambda)
      kurt3 = 4.0_dp * real(5 * n ** 3 - 57 * n ** 2 + 25 * n + 169, dp) &
         / real((n - 2) * (n - 1), dp)
      res%excess_kurtosis = kurt1 * (kurt2 - kurt3)
      if (res%variance > 0.0_dp) then
         res%z_score = (res%delta - res%expectation) / sqrt(res%variance)
         select case (trim(alt))
         case ("two.sided")
            res%p_value = normal_two_sided_p(res%z_score)
         case ("less")
            res%p_value = normal_cdf(res%z_score)
         case default
            res%p_value = 1.0_dp - normal_cdf(res%z_score)
         end select
      else
         res%z_score = safe_nan()
         res%p_value = safe_nan()
      end if
   end function spatialdelta

   pure function linearised_diffusive_weights(adjacency_matrix, regional_weights, &
      t_choice) result(w)
      real(dp), intent(in) :: adjacency_matrix(:, :) !! Square adjacency matrix, normally symmetric with zero diagonal.
      real(dp), intent(in) :: regional_weights(:) !! Strictly positive regional masses; internally normalized to sum to one.
      integer, intent(in), optional :: t_choice !! Step rule: 1 uses min(r_i/degree_i); 2 uses inverse largest eigenvalue.
      real(dp), allocatable :: w(:, :)
      real(dp), allocatable :: rw(:)
      real(dp), allocatable :: row_sum(:)
      real(dp), allocatable :: lap(:, :)
      real(dp), allocatable :: la(:, :)
      real(dp) :: t
      real(dp) :: lambda_max
      integer :: choice
      integer :: n
      integer :: i
      integer :: j

      n = size(regional_weights)
      allocate(w(n, n))
      w = safe_nan()
      if (n == 0 .or. size(adjacency_matrix, 1) /= n .or. &
         size(adjacency_matrix, 2) /= n .or. any(regional_weights <= 0.0_dp)) then
         return
      end if
      rw = regional_weights / sum(regional_weights)
      row_sum = sum(adjacency_matrix, dim = 2)
      allocate(lap(n, n), la(n, n))
      lap = -adjacency_matrix
      do i = 1, n
         lap(i, i) = lap(i, i) + row_sum(i)
      end do
      do i = 1, n
         do j = 1, n
            la(i, j) = lap(i, j) / sqrt(rw(i) * rw(j))
         end do
      end do
      choice = 2
      if (present(t_choice)) choice = t_choice
      if (choice == 1) then
         t = huge(1.0_dp)
         do i = 1, n
            if (row_sum(i) > 0.0_dp) t = min(t, rw(i) / row_sum(i))
         end do
         if (t == huge(1.0_dp)) t = 0.0_dp
      else
         lambda_max = symmetric_max_eigenvalue(la)
         if (lambda_max > 0.0_dp) then
            t = 1.0_dp / lambda_max
         else
            t = 0.0_dp
         end if
      end if
      w = 0.0_dp
      do i = 1, n
         w(i, i) = 1.0_dp
         w(i, :) = w(i, :) - t * lap(i, :) / rw(i)
      end do
   end function linearised_diffusive_weights

   pure function metropolis_hastings_weights(adjacency_matrix, regional_weights) result(w)
      real(dp), intent(in) :: adjacency_matrix(:, :) !! Square nonnegative adjacency matrix used as the proposal graph.
      real(dp), intent(in) :: regional_weights(:) !! Positive target stationary masses; normalized internally to sum to one.
      real(dp), allocatable :: w(:, :)
      real(dp), allocatable :: rw(:)
      real(dp), allocatable :: row_sum(:)
      real(dp), allocatable :: p(:, :)
      real(dp), allocatable :: fp(:, :)
      real(dp), allocatable :: g(:, :)
      real(dp) :: offsum
      integer :: n
      integer :: i
      integer :: j

      n = size(regional_weights)
      allocate(w(n, n))
      w = safe_nan()
      if (n == 0 .or. size(adjacency_matrix, 1) /= n .or. &
         size(adjacency_matrix, 2) /= n .or. any(regional_weights <= 0.0_dp)) return
      rw = regional_weights / sum(regional_weights)
      row_sum = sum(adjacency_matrix, dim = 2)
      allocate(p(n, n), fp(n, n), g(n, n))
      p = 0.0_dp
      do i = 1, n
         if (row_sum(i) > 0.0_dp) p(i, :) = adjacency_matrix(i, :) / row_sum(i)
      end do
      do i = 1, n
         fp(i, :) = rw(i) * p(i, :)
      end do
      do i = 1, n
         do j = 1, n
            g(i, j) = min(fp(i, j), fp(j, i))
         end do
      end do
      w = 0.0_dp
      do i = 1, n
         do j = 1, n
            w(i, j) = g(i, j) / rw(i)
         end do
         offsum = sum(g(i, :))
         w(i, i) = w(i, i) + (rw(i) - offsum) / rw(i)
      end do
   end function metropolis_hastings_weights

   pure function iterative_proportional_fitting_weights(adjacency_matrix, &
      regional_weights, g, max_iter, tol) result(w)
      real(dp), intent(in) :: adjacency_matrix(:, :) !! Square adjacency seed; a small positive g is added to every cell.
      real(dp), intent(in) :: regional_weights(:) !! Positive target row and column margins; internally normalized to one.
      real(dp), intent(in), optional :: g !! Positive additive seed smoothing constant; default is 0.001.
      integer, intent(in), optional :: max_iter !! Maximum alternating row/column scaling iterations; default is 1000.
      real(dp), intent(in), optional :: tol !! Maximum target-margin error required for convergence; default is 1e-10.
      real(dp), allocatable :: w(:, :)
      real(dp), allocatable :: rw(:)
      real(dp), allocatable :: fitted(:, :)
      real(dp), allocatable :: rs(:)
      real(dp), allocatable :: cs(:)
      real(dp) :: gg
      real(dp) :: eps
      integer :: nit
      integer :: iter
      integer :: i
      integer :: j
      integer :: n

      n = size(regional_weights)
      allocate(w(n, n))
      w = safe_nan()
      if (n == 0 .or. size(adjacency_matrix, 1) /= n .or. &
         size(adjacency_matrix, 2) /= n .or. any(regional_weights <= 0.0_dp)) return
      gg = 0.001_dp
      if (present(g)) gg = max(0.0_dp, g)
      nit = 1000
      if (present(max_iter)) nit = max(1, max_iter)
      eps = 1.0e-10_dp
      if (present(tol)) eps = max(0.0_dp, tol)
      rw = regional_weights / sum(regional_weights)
      allocate(fitted(n, n), rs(n), cs(n))
      fitted = adjacency_matrix + gg
      do iter = 1, nit
         rs = sum(fitted, dim = 2)
         do i = 1, n
            if (rs(i) > 0.0_dp) fitted(i, :) = fitted(i, :) * rw(i) / rs(i)
         end do
         cs = sum(fitted, dim = 1)
         do j = 1, n
            if (cs(j) > 0.0_dp) fitted(:, j) = fitted(:, j) * rw(j) / cs(j)
         end do
         rs = sum(fitted, dim = 2)
         cs = sum(fitted, dim = 1)
         if (max(maxval(abs(rs - rw)), maxval(abs(cs - rw))) <= eps) exit
      end do
      do i = 1, n
         w(i, :) = fitted(i, :) / rw(i)
      end do
   end function iterative_proportional_fitting_weights

   pure function graph_distance_weights(adjacency_matrix, regional_weights, c) result(w)
      real(dp), intent(in) :: adjacency_matrix(:, :) !! Square adjacency matrix interpreted as an unweighted graph.
      real(dp), intent(in) :: regional_weights(:) !! Strictly positive regional masses; internally normalized to sum to one.
      real(dp), intent(in), optional :: c !! Positive scaling constant not exceeding -1/min(B); default is that upper bound.
      real(dp), allocatable :: w(:, :)
      real(dp), allocatable :: rw(:)
      real(dp), allocatable :: d(:, :)
      real(dp), allocatable :: h(:, :)
      real(dp), allocatable :: b(:, :)
      real(dp) :: cc
      real(dp) :: cmax
      integer :: n
      integer :: i
      integer :: j
      integer :: k

      n = size(regional_weights)
      allocate(w(n, n))
      w = safe_nan()
      if (n == 0 .or. size(adjacency_matrix, 1) /= n .or. &
         size(adjacency_matrix, 2) /= n .or. any(regional_weights <= 0.0_dp)) return
      rw = regional_weights / sum(regional_weights)
      allocate(d(n, n))
      d = huge(1.0_dp)
      do i = 1, n
         d(i, i) = 0.0_dp
         do j = 1, n
            if (adjacency_matrix(i, j) /= 0.0_dp) d(i, j) = 1.0_dp
         end do
      end do
      do k = 1, n
         do i = 1, n
            do j = 1, n
               if (d(i, k) < huge(1.0_dp) / 4.0_dp .and. &
                  d(k, j) < huge(1.0_dp) / 4.0_dp) then
                  d(i, j) = min(d(i, j), d(i, k) + d(k, j))
               end if
            end do
         end do
      end do
      if (any(d >= huge(1.0_dp) / 4.0_dp)) then
         w = metropolis_hastings_weights(adjacency_matrix, rw)
         return
      end if
      allocate(h(n, n), b(n, n))
      h = 0.0_dp
      do i = 1, n
         h(i, i) = 1.0_dp
         do j = 1, n
            h(i, j) = h(i, j) - rw(j)
         end do
      end do
      b = -0.5_dp * matmul(matmul(h, d), transpose(h))
      if (minval(b) >= 0.0_dp) then
         w = metropolis_hastings_weights(adjacency_matrix, rw)
         return
      end if
      cmax = -1.0_dp / minval(b)
      cc = cmax
      if (present(c)) cc = c
      if (cc <= 0.0_dp .or. cc > cmax) then
         w = safe_nan()
         return
      end if
      do i = 1, n
         do j = 1, n
            w(i, j) = (1.0_dp + cc * b(i, j)) * rw(j)
         end do
      end do
   end function graph_distance_weights

   pure function localdelta(dissimilarity_matrix, adjusted_spatial_weights, &
      regional_weights) result(di)
      real(dp), intent(in) :: dissimilarity_matrix(:, :) !! Square pairwise dissimilarity matrix for regional features.
      real(dp), intent(in) :: adjusted_spatial_weights(:, :) !! Square adjusted spatial transition-weight matrix.
      real(dp), intent(in) :: regional_weights(:) !! Positive regional masses used for weighted double centering.
      real(dp), allocatable :: di(:)
      real(dp), allocatable :: rw(:)
      real(dp), allocatable :: h(:, :)
      real(dp), allocatable :: b(:, :)
      real(dp), allocatable :: product(:, :)
      real(dp) :: trace_kx
      integer :: n
      integer :: i
      integer :: j

      n = size(regional_weights)
      allocate(di(n))
      di = safe_nan()
      if (n == 0 .or. size(dissimilarity_matrix, 1) /= n .or. &
         size(dissimilarity_matrix, 2) /= n .or. &
         size(adjusted_spatial_weights, 1) /= n .or. &
         size(adjusted_spatial_weights, 2) /= n .or. any(regional_weights <= 0.0_dp)) return
      rw = regional_weights / sum(regional_weights)
      allocate(h(n, n), b(n, n))
      h = 0.0_dp
      do i = 1, n
         h(i, i) = 1.0_dp
         do j = 1, n
            h(i, j) = h(i, j) - rw(j)
         end do
      end do
      b = -0.5_dp * matmul(matmul(h, dissimilarity_matrix), transpose(h))
      trace_kx = 0.0_dp
      do i = 1, n
         trace_kx = trace_kx + rw(i) * b(i, i)
      end do
      if (trace_kx == 0.0_dp) return
      product = matmul(adjusted_spatial_weights, b)
      do i = 1, n
         di(i) = product(i, i) / trace_kx
      end do
   end function localdelta

   pure function cornish_fisher(input, alternative) result(output)
      type(spatial_delta_result), intent(in) :: input !! Spatial-delta result whose standardized deviate receives the correction.
      character(len=*), intent(in), optional :: alternative !! Tail alternative: greater, less, or two.sided; default greater.
      type(spatial_delta_result) :: output
      character(len=12) :: alt
      real(dp) :: s
      real(dp) :: k
      real(dp) :: s2
      real(dp) :: k8
      real(dp) :: domain
      real(dp) :: z
      real(dp) :: correction

      output = input
      alt = "greater"
      if (present(alternative)) alt = adjustl(alternative)
      s = input%skewness
      k = input%excess_kurtosis
      s2 = s * s
      k8 = k / 8.0_dp
      domain = s2 / 9.0_dp - 4.0_dp * (k8 - s2 / 6.0_dp) &
         * (1.0_dp - k8 - 5.0_dp * s2 / 36.0_dp)
      if (domain > 0.0_dp) return
      z = input%z_score
      correction = (s / 6.0_dp) * (z * z - 1.0_dp) &
         + (k / 24.0_dp) * (z ** 3 - 3.0_dp * z) &
         - (s2 / 36.0_dp) * (2.0_dp * z ** 3 - 5.0_dp * z)
      if (z <= 0.0_dp) then
         output%z_score = z + correction
      else
         output%z_score = z - correction
      end if
      select case (trim(alt))
      case ("two.sided")
         output%p_value = normal_two_sided_p(output%z_score)
      case ("less")
         output%p_value = normal_cdf(output%z_score)
      case default
         output%p_value = 1.0_dp - normal_cdf(output%z_score)
      end select
   end function cornish_fisher

   pure subroutine set_delta_nan(res)
      type(spatial_delta_result), intent(out) :: res !! Result object filled with quiet NaNs for an invalid delta calculation.

      res%delta = safe_nan()
      res%rv = safe_nan()
      res%expectation = safe_nan()
      res%variance = safe_nan()
      res%variance_product = safe_nan()
      res%z_score = safe_nan()
      res%p_value = safe_nan()
      res%skewness = safe_nan()
      res%excess_kurtosis = safe_nan()
   end subroutine set_delta_nan

end module spdep_delta
