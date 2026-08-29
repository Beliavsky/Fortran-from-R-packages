module lavaan_muthen1984
   use lavaan_kinds, only : dp
   use lavaan_categorical, only : categorical_stats_result
   use lavaan_ordinal, only : ordinal_thresholds, polychoric_matrix, bvn_rectangle
   use lavaan_linalg, only : inverse_general
   implicit none
   private

   type, public :: muthen1984_result
      type(categorical_stats_result) :: categorical
      real(dp), allocatable :: score(:, :), inner(:, :)
      real(dp), allocatable :: a11(:, :), a21(:, :), a22(:, :), b_inv(:, :)
      integer :: status = 0
   end type muthen1984_result

   public :: muthen1984_ordinal, categorical_wls_statistics_muthen

contains

   subroutine categorical_wls_statistics_muthen(data, result)
      integer, intent(in) :: data(:, :)
      type(categorical_stats_result), intent(out) :: result
      type(muthen1984_result) :: fit
      call muthen1984_ordinal(data, fit)
      result = fit%categorical
      result%status = fit%status
   end subroutine categorical_wls_statistics_muthen

   subroutine muthen1984_ordinal(data, result)
      integer, intent(in) :: data(:, :)
      type(muthen1984_result), intent(out) :: result
      real(dp), allocatable :: sc_th(:, :), sc_cor(:, :), pair_th(:, :)
      real(dp), allocatable :: binv11(:, :), binv22(:, :), tmp(:, :), sandwich(:, :)
      real(dp), allocatable :: ridge(:, :), winv(:, :)
      integer :: n, p, nth, ncor, q, i, j, k, info, off1, off2

      n = size(data, 1)
      p = size(data, 2)
      if (n < 5 .or. p < 2 .or. minval(data) < 1) then
         result%status = -1
         result%categorical%status = -1
         return
      end if
      call setup_statistics(data, result%categorical, info)
      if (info /= 0) then
         result%status = info
         result%categorical%status = info
         return
      end if
      nth = size(result%categorical%thresholds)
      ncor = p * (p - 1) / 2
      q = nth + ncor
      allocate(sc_th(n, nth), sc_cor(n, ncor))
      call threshold_scores(data, result%categorical, sc_th, info)
      if (info /= 0) then
         result%status = info
         result%categorical%status = info
         return
      end if

      allocate(result%a11(nth, nth))
      result%a11 = 0.0_dp
      do j = 1, p
         off1 = result%categorical%threshold_offset(j)
         off2 = result%categorical%threshold_offset(j + 1) - 1
         result%a11(off1:off2, off1:off2) = &
            matmul(transpose(sc_th(:, off1:off2)), sc_th(:, off1:off2))
      end do

      allocate(result%a21(ncor, nth), result%a22(ncor, ncor))
      result%a21 = 0.0_dp
      result%a22 = 0.0_dp
      sc_cor = 0.0_dp
      k = 0
      do j = 1, p - 1
         do i = j + 1, p
            k = k + 1
            allocate(pair_th(n, nth))
            pair_th = 0.0_dp
            call pair_scores(data(:, j), data(:, i), j, i, result%categorical, &
                             result%categorical%correlation(i, j), sc_cor(:, k), pair_th, info)
            if (info /= 0) then
               result%status = 1000 + info
               result%categorical%status = result%status
               return
            end if
            result%a21(k, :) = matmul(sc_cor(:, k), pair_th)
            result%a22(k, k) = dot_product(sc_cor(:, k), sc_cor(:, k))
            deallocate(pair_th)
         end do
      end do

      call inverse_general(result%a11, binv11, info)
      if (info /= 0) then
         result%status = 2000 + info
         result%categorical%status = result%status
         return
      end if
      allocate(binv22(ncor, ncor))
      binv22 = 0.0_dp
      do k = 1, ncor
         if (result%a22(k, k) <= tiny(1.0_dp)) then
            result%status = 3000 + k
            result%categorical%status = result%status
            return
         end if
         binv22(k, k) = 1.0_dp / result%a22(k, k)
      end do
      allocate(result%b_inv(q, q))
      result%b_inv = 0.0_dp
      result%b_inv(1:nth, 1:nth) = binv11
      result%b_inv(nth + 1:q, nth + 1:q) = binv22
      tmp = -matmul(binv22, matmul(result%a21, binv11))
      result%b_inv(nth + 1:q, 1:nth) = tmp

      allocate(result%score(n, q))
      result%score(:, 1:nth) = sc_th
      result%score(:, nth + 1:q) = sc_cor
      result%inner = matmul(transpose(result%score), result%score)
      sandwich = matmul(result%b_inv, matmul(result%inner, transpose(result%b_inv)))
      allocate(result%categorical%gamma(q, q))
      result%categorical%gamma = real(n, dp) * 0.5_dp * (sandwich + transpose(sandwich))

      ridge = result%categorical%gamma
      do k = 1, q
         ridge(k, k) = ridge(k, k) + 1.0e-10_dp * max(1.0_dp, abs(ridge(k, k)))
      end do
      call inverse_general(ridge, winv, info)
      allocate(result%categorical%weight(q, q), result%categorical%dwls_weight(q))
      result%categorical%weight = 0.0_dp
      if (info == 0) then
         result%categorical%weight = winv
      else
         do k = 1, q
            if (result%categorical%gamma(k, k) > tiny(1.0_dp)) &
               result%categorical%weight(k, k) = 1.0_dp / result%categorical%gamma(k, k)
         end do
      end if
      do k = 1, q
         result%categorical%dwls_weight(k) = result%categorical%weight(k, k)
      end do
      result%categorical%n_jackknife = 0
      result%categorical%status = 0
      result%status = 0
   end subroutine muthen1984_ordinal

   subroutine setup_statistics(data, result, info)
      integer, intent(in) :: data(:, :)
      type(categorical_stats_result), intent(inout) :: result
      integer, intent(out) :: info
      integer :: p, n, j, i, nth, pos, ncor, k
      integer, allocatable :: counts(:)
      real(dp), allocatable :: th(:)
      p = size(data, 2)
      n = size(data, 1)
      info = 0
      allocate(result%ncat(p), result%threshold_offset(p + 1))
      result%threshold_offset(1) = 1
      nth = 0
      do j = 1, p
         result%ncat(j) = maxval(data(:, j))
         if (result%ncat(j) < 2 .or. minval(data(:, j)) < 1) then
            info = j
            return
         end if
         nth = nth + result%ncat(j) - 1
         result%threshold_offset(j + 1) = nth + 1
      end do
      allocate(result%thresholds(nth))
      pos = 1
      do j = 1, p
         allocate(counts(result%ncat(j)))
         counts = 0
         do i = 1, n
            counts(data(i, j)) = counts(data(i, j)) + 1
         end do
         if (any(counts == 0)) then
            info = 100 + j
            return
         end if
         th = ordinal_thresholds(counts)
         result%thresholds(pos:pos + size(th) - 1) = th
         pos = pos + size(th)
         deallocate(counts)
      end do
      call polychoric_matrix(data, result%correlation, info)
      if (info /= 0) return
      ncor = p * (p - 1) / 2
      allocate(result%stats(nth + ncor))
      result%stats(1:nth) = result%thresholds
      k = nth
      do j = 1, p - 1
         do i = j + 1, p
            k = k + 1
            result%stats(k) = result%correlation(i, j)
         end do
      end do
   end subroutine setup_statistics

   subroutine threshold_scores(data, result, score, info)
      integer, intent(in) :: data(:, :)
      type(categorical_stats_result), intent(in) :: result
      real(dp), intent(out) :: score(:, :)
      integer, intent(out) :: info
      integer :: n, p, r, j, c, base, idx
      real(dp) :: lo, hi, prob
      n = size(data, 1)
      p = size(data, 2)
      score = 0.0_dp
      info = 0
      do j = 1, p
         base = result%threshold_offset(j)
         do r = 1, n
            c = data(r, j)
            call marginal_bounds(j, c, result, lo, hi)
            prob = max(normal_cdf(hi) - normal_cdf(lo), 1.0e-300_dp)
            if (c < result%ncat(j)) then
               idx = base + c - 1
               score(r, idx) = score(r, idx) + normal_pdf(hi) / prob
            end if
            if (c > 1) then
               idx = base + c - 2
               score(r, idx) = score(r, idx) - normal_pdf(lo) / prob
            end if
         end do
      end do
   end subroutine threshold_scores

   subroutine pair_scores(y1, y2, j1, j2, result, rho_in, score_rho, score_th, info)
      integer, intent(in) :: y1(:), y2(:), j1, j2
      type(categorical_stats_result), intent(in) :: result
      real(dp), intent(in) :: rho_in
      real(dp), intent(out) :: score_rho(:), score_th(:, :)
      integer, intent(out) :: info
      integer :: r, c1, c2, idx, base1, base2
      real(dp) :: l1, u1, l2, u2, rho, prob, drho, dth, s
      if (size(y1) /= size(y2) .or. size(score_rho) /= size(y1) .or. &
          size(score_th, 1) /= size(y1)) then
         info = -1
         return
      end if
      rho = max(-0.999_dp, min(0.999_dp, rho_in))
      s = sqrt(max(1.0e-12_dp, 1.0_dp - rho * rho))
      score_rho = 0.0_dp
      score_th = 0.0_dp
      info = 0
      base1 = result%threshold_offset(j1)
      base2 = result%threshold_offset(j2)
      do r = 1, size(y1)
         c1 = y1(r)
         c2 = y2(r)
         call marginal_bounds(j1, c1, result, l1, u1)
         call marginal_bounds(j2, c2, result, l2, u2)
         prob = bvn_rectangle(l1, u1, l2, u2, rho)
         drho = bvn_pdf(u1, u2, rho) - bvn_pdf(l1, u2, rho) - &
                bvn_pdf(u1, l2, rho) + bvn_pdf(l1, l2, rho)
         score_rho(r) = drho / prob
         if (c1 < result%ncat(j1)) then
            idx = base1 + c1 - 1
            dth = normal_pdf(u1) * (normal_cdf((u2 - rho * u1) / s) - &
                                    normal_cdf((l2 - rho * u1) / s))
            score_th(r, idx) = score_th(r, idx) + dth / prob
         end if
         if (c1 > 1) then
            idx = base1 + c1 - 2
            dth = -normal_pdf(l1) * (normal_cdf((u2 - rho * l1) / s) - &
                                     normal_cdf((l2 - rho * l1) / s))
            score_th(r, idx) = score_th(r, idx) + dth / prob
         end if
         if (c2 < result%ncat(j2)) then
            idx = base2 + c2 - 1
            dth = normal_pdf(u2) * (normal_cdf((u1 - rho * u2) / s) - &
                                    normal_cdf((l1 - rho * u2) / s))
            score_th(r, idx) = score_th(r, idx) + dth / prob
         end if
         if (c2 > 1) then
            idx = base2 + c2 - 2
            dth = -normal_pdf(l2) * (normal_cdf((u1 - rho * l2) / s) - &
                                     normal_cdf((l1 - rho * l2) / s))
            score_th(r, idx) = score_th(r, idx) + dth / prob
         end if
      end do
   end subroutine pair_scores

   subroutine marginal_bounds(j, c, result, lo, hi)
      integer, intent(in) :: j, c
      type(categorical_stats_result), intent(in) :: result
      real(dp), intent(out) :: lo, hi
      integer :: base
      base = result%threshold_offset(j)
      if (c == 1) then
         lo = -huge(1.0_dp)
      else
         lo = result%thresholds(base + c - 2)
      end if
      if (c == result%ncat(j)) then
         hi = huge(1.0_dp)
      else
         hi = result%thresholds(base + c - 1)
      end if
   end subroutine marginal_bounds

   pure function normal_cdf(x) result(p)
      real(dp), intent(in) :: x
      real(dp) :: p
      if (x > 8.0_dp) then
         p = 1.0_dp
      else if (x < -8.0_dp) then
         p = 0.0_dp
      else
         p = 0.5_dp * erfc(-x / sqrt(2.0_dp))
      end if
   end function normal_cdf

   pure function normal_pdf(x) result(v)
      real(dp), intent(in) :: x
      real(dp) :: v
      if (abs(x) > 30.0_dp) then
         v = 0.0_dp
      else
         v = exp(-0.5_dp * x * x) / sqrt(2.0_dp * acos(-1.0_dp))
      end if
   end function normal_pdf

   pure function bvn_pdf(x, y, rho) result(v)
      real(dp), intent(in) :: x, y, rho
      real(dp) :: v, omr2
      if (abs(x) > 30.0_dp .or. abs(y) > 30.0_dp) then
         v = 0.0_dp
         return
      end if
      omr2 = max(1.0e-12_dp, 1.0_dp - rho * rho)
      v = exp(-(x * x - 2.0_dp * rho * x * y + y * y) / (2.0_dp * omr2)) / &
          (2.0_dp * acos(-1.0_dp) * sqrt(omr2))
   end function bvn_pdf

end module lavaan_muthen1984
