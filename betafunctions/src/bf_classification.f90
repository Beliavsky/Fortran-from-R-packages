! SPDX-License-Identifier: CC0-1.0
module bf_classification
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
  use bf_kinds, only: dp
  use bf_special, only: binomial_prob_range, chi_square_sf
  use bf_quadrature, only: integrate_gk
  use bf_distributions, only: beta4_pdf, compound_binomial_pmf
  use bf_moments, only: beta_params, beta_true_score_fit, hb_beta_true_score_fit
  implicit none
  private

  public :: accuracy_stats, consistency_stats, classification_result, roc_result, omega_result, model_fit_result
  public :: etl, reliability_from_etl, confmat, ca_stats, cc_stats, auc
  public :: cronbach_alpha, mcdonald_omega, lords_k
  public :: ll_classify, ll_classify_params, hb_classify, hb_classify_params
  public :: ll_roc, ll_roc_params, hb_roc, hb_roc_params
  public :: ll_model_fit, hb_model_fit

  type :: accuracy_stats
    real(dp) :: sensitivity = 0.0_dp
    real(dp) :: specificity = 0.0_dp
    real(dp) :: ppv = 0.0_dp
    real(dp) :: npv = 0.0_dp
    real(dp) :: youden_j = 0.0_dp
    real(dp) :: accuracy = 0.0_dp
  end type accuracy_stats

  type :: consistency_stats
    real(dp) :: p = 0.0_dp
    real(dp) :: p_c = 0.0_dp
    real(dp) :: p_c_pos = 0.0_dp
    real(dp) :: p_c_neg = 0.0_dp
    real(dp) :: kappa = 0.0_dp
  end type consistency_stats

  type :: classification_result
    type(beta_params) :: parameters
    real(dp), allocatable :: accuracy_matrix(:,:)
    real(dp), allocatable :: consistency_matrix(:,:)
    type(accuracy_stats), allocatable :: category_accuracy(:)
    type(consistency_stats), allocatable :: category_consistency(:)
    real(dp) :: overall_accuracy = 0.0_dp
    type(consistency_stats) :: overall_consistency
  end type classification_result

  type :: roc_result
    ! Columns: FPR, TPR, Youden.J, Cutoff, Accuracy, PPV, NPV
    real(dp), allocatable :: table(:,:)
    real(dp) :: area = 0.0_dp
    integer :: max_youden_index = 0
    integer :: max_accuracy_index = 0
  end type roc_result

  type :: omega_result
    real(dp) :: omega = 0.0_dp
    real(dp) :: gfi = 0.0_dp
    real(dp), allocatable :: loadings(:)
    real(dp), allocatable :: error_variances(:)
    real(dp), allocatable :: observed(:,:)
    real(dp), allocatable :: fitted(:,:)
    real(dp), allocatable :: discrepancy(:,:)
  end type omega_result

  type :: model_fit_result
    real(dp), allocatable :: expected(:)
    real(dp), allocatable :: observed(:)
    real(dp) :: chi_square = 0.0_dp
    integer :: df = 0
    real(dp) :: p_value = 0.0_dp
  end type model_fit_result

  type :: class_integrand_ctx
    type(beta_params) :: par
    integer :: n = 0
    integer :: obs_lo = 0
    integer :: obs_hi = 0
    integer :: obs2_lo = 0
    integer :: obs2_hi = 0
    logical :: use_hb = .false.
    logical :: product_two = .false.
  end type class_integrand_ctx

  type :: fit_integrand_ctx
    type(beta_params) :: par
    integer :: n = 0
    integer :: lo = 0
    integer :: hi = 0
    logical :: use_hb = .false.
  end type fit_integrand_ctx

contains

  pure real(dp) function safe_ratio(a, b) result(v)
    real(dp), intent(in) :: a, b
    if (abs(b) <= tiny(1.0_dp)) then
      v = ieee_value(0.0_dp, ieee_quiet_nan)
    else
      v = a / b
    end if
  end function safe_ratio

  pure real(dp) function etl(mean_value, variance, min_value, max_value, reliability) result(n)
    real(dp), intent(in) :: mean_value, variance, min_value, max_value, reliability
    n = ((mean_value - min_value) * (max_value - mean_value) - reliability * variance) / &
        (variance * (1.0_dp - reliability))
  end function etl

  pure real(dp) function reliability_from_etl(mean_value, variance, min_value, max_value, n_eff) result(r)
    real(dp), intent(in) :: mean_value, variance, min_value, max_value, n_eff
    r = (n_eff * variance + min_value * (max_value - mean_value) + mean_value * (mean_value - max_value)) / &
        ((n_eff - 1.0_dp) * variance)
  end function reliability_from_etl

  subroutine confmat(tp, tn, fp, fn, mat, proportions)
    real(dp), intent(in) :: tp, tn, fp, fn
    real(dp), intent(out) :: mat(3,3)
    logical, intent(in), optional :: proportions
    logical :: prop
    real(dp) :: total

    prop = .false.
    if (present(proportions)) prop = proportions
    mat = 0.0_dp
    mat(1,1) = tp
    mat(1,2) = tn
    mat(2,1) = fp
    mat(2,2) = fn
    mat(1,3) = sum(mat(1,1:2))
    mat(2,3) = sum(mat(2,1:2))
    mat(3,1) = sum(mat(1:2,1))
    mat(3,2) = sum(mat(1:2,2))
    mat(3,3) = sum(mat(1:2,1:2))
    if (prop) then
      total = tp + tn + fp + fn
      if (abs(total) > tiny(1.0_dp)) mat = mat / total
    end if
  end subroutine confmat

  pure function ca_stats(tp, tn, fp, fn) result(s)
    real(dp), intent(in) :: tp, tn, fp, fn
    type(accuracy_stats) :: s
    s%sensitivity = safe_ratio(tp, tp + fn)
    s%specificity = safe_ratio(tn, tn + fp)
    s%ppv = safe_ratio(tp, tp + fp)
    s%npv = safe_ratio(tn, tn + fn)
    s%accuracy = safe_ratio(tp + tn, tp + tn + fp + fn)
    s%youden_j = s%sensitivity + s%specificity - 1.0_dp
  end function ca_stats

  pure function cc_stats(ii, ij, ji, jj) result(s)
    real(dp), intent(in) :: ii, ij, ji, jj
    type(consistency_stats) :: s
    real(dp) :: total
    total = ii + ij + ji + jj
    s%p = safe_ratio(ii + jj, total)
    s%p_c_pos = (ii + ij) * (ii + ji)
    s%p_c_neg = (ij + jj) * (ji + jj)
    s%p_c = s%p_c_pos + s%p_c_neg
    s%kappa = safe_ratio(s%p - s%p_c, 1.0_dp - s%p_c)
  end function cc_stats

  pure real(dp) function auc(fpr, tpr) result(a)
    real(dp), intent(in) :: fpr(:), tpr(:)
    integer :: i, n
    n = min(size(fpr), size(tpr))
    a = 0.0_dp
    do i = 1, n - 1
      a = a + 0.5_dp * (tpr(i) + tpr(i + 1)) * (fpr(i + 1) - fpr(i))
    end do
  end function auc

  subroutine sample_covariance(x, cov)
    real(dp), intent(in) :: x(:,:)
    real(dp), intent(out) :: cov(size(x,2), size(x,2))
    real(dp), allocatable :: means(:)
    integer :: i, j, n
    n = size(x,1)
    allocate(means(size(x,2)))
    do j = 1, size(x,2)
      means(j) = sum(x(:,j)) / real(n, dp)
    end do
    do i = 1, size(x,2)
      do j = 1, size(x,2)
        cov(i,j) = sum((x(:,i) - means(i)) * (x(:,j) - means(j))) / real(max(n - 1, 1), dp)
      end do
    end do
  end subroutine sample_covariance

  real(dp) function cronbach_alpha(x) result(alpha)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable :: cov(:,:)
    integer :: p, i
    real(dp) :: tracev
    p = size(x,2)
    allocate(cov(p,p))
    call sample_covariance(x, cov)
    tracev = 0.0_dp
    do i = 1, p
      tracev = tracev + cov(i,i)
    end do
    alpha = real(p, dp) / real(p - 1, dp) * (1.0_dp - tracev / sum(cov))
  end function cronbach_alpha

  subroutine mcdonald_omega(x, out)
    real(dp), intent(in) :: x(:,:)
    type(omega_result), intent(out) :: out
    integer :: p, i, j, k, count
    real(dp) :: s, denom

    p = size(x,2)
    allocate(out%observed(p,p), out%fitted(p,p), out%discrepancy(p,p))
    allocate(out%loadings(p), out%error_variances(p))
    call sample_covariance(x, out%observed)

    do i = 1, p
      do j = 1, p
        if (i /= j .and. out%observed(i,j) <= 0.0_dp) then
          out%omega = ieee_value(0.0_dp, ieee_quiet_nan)
          out%gfi = out%omega
          out%loadings = out%omega
          out%error_variances = out%omega
          out%fitted = out%omega
          out%discrepancy = out%omega
          return
        end if
      end do
    end do

    do i = 1, p
      s = 0.0_dp
      count = 0
      do j = 1, p - 1
        if (j == i) cycle
        do k = j + 1, p
          if (k == i) cycle
          if (out%observed(j,k) > 0.0_dp) then
            s = s + sqrt(out%observed(i,j) * out%observed(i,k) / out%observed(j,k))
            count = count + 1
          end if
        end do
      end do
      if (count > 0) then
        out%loadings(i) = s / real(count, dp)
      else
        out%loadings(i) = sqrt(max(out%observed(i,i), 0.0_dp))
      end if
      out%error_variances(i) = out%observed(i,i) - out%loadings(i)**2
    end do

    denom = sum(out%error_variances) + sum(out%loadings)**2
    out%omega = safe_ratio(sum(out%loadings)**2, denom)
    do i = 1, p
      do j = 1, p
        if (i == j) then
          out%fitted(i,j) = out%loadings(i)**2 + out%error_variances(i)
        else
          out%fitted(i,j) = out%loadings(i) * out%loadings(j)
        end if
      end do
    end do
    out%discrepancy = out%observed - out%fitted
    out%gfi = 1.0_dp - sum(out%discrepancy**2) / real(p*p, dp) / &
              (sum(out%observed) / real(p*p, dp))
  end subroutine mcdonald_omega

  real(dp) function lords_k(x, n, reliability) result(k)
    real(dp), intent(in) :: x(:), reliability
    integer, intent(in) :: n
    real(dp) :: mu, sigma2, sigma2e, num, den
    mu = sum(x) / real(size(x), dp)
    sigma2 = sum((x - mu)**2) / real(max(size(x) - 1, 1), dp)
    sigma2e = sigma2 * (1.0_dp - reliability)
    num = real(n, dp) * (real(n - 1, dp) * (sigma2 - sigma2e) - real(n, dp) * sigma2 + mu * (real(n,dp) - mu))
    den = 2.0_dp * (mu * (real(n,dp) - mu) - (sigma2 - sigma2e))
    k = num / den
  end function lords_k

  pure real(dp) function observed_category_prob(x, ctx, second) result(p)
    real(dp), intent(in) :: x
    type(class_integrand_ctx), intent(in) :: ctx
    logical, intent(in) :: second
    integer :: lo, hi, j

    if (second) then
      lo = ctx%obs2_lo
      hi = ctx%obs2_hi
    else
      lo = ctx%obs_lo
      hi = ctx%obs_hi
    end if
    if (.not. ctx%use_hb) then
      p = binomial_prob_range(lo, hi, ctx%n, x)
    else
      p = 0.0_dp
      do j = lo, hi
        p = p + compound_binomial_pmf(j, ctx%n, ctx%par%k, x)
      end do
    end if
  end function observed_category_prob

  real(dp) function class_integrand(x, anyctx) result(v)
    real(dp), intent(in) :: x
    class(*), intent(in) :: anyctx
    real(dp) :: p1, p2
    select type (ctx => anyctx)
    type is (class_integrand_ctx)
      p1 = observed_category_prob(x, ctx, .false.)
      if (ctx%product_two) then
        p2 = observed_category_prob(x, ctx, .true.)
      else
        p2 = 1.0_dp
      end if
      v = beta4_pdf(x, ctx%par%l, ctx%par%u, ctx%par%alpha, ctx%par%beta) * p1 * p2
    class default
      v = 0.0_dp
    end select
  end function class_integrand

  subroutine make_category_ranges(cuts, n, lo, hi)
    real(dp), intent(in) :: cuts(:)
    integer, intent(in) :: n
    integer, intent(out) :: lo(size(cuts)+1), hi(size(cuts)+1)
    integer :: i, nc
    integer, allocatable :: icut(:)
    nc = size(cuts) + 1
    allocate(icut(size(cuts)))
    do i = 1, size(cuts)
      icut(i) = nint(cuts(i) * real(n, dp))
    end do
    lo(1) = 0
    if (size(cuts) > 0) hi(1) = icut(1) - 1
    do i = 2, nc - 1
      lo(i) = icut(i - 1)
      hi(i) = icut(i) - 1
    end do
    if (nc > 1) then
      lo(nc) = icut(size(cuts))
      hi(nc) = n
    else
      hi(1) = n
    end if
  end subroutine make_category_ranges

  subroutine make_hb_ranges(cuts, n, lo, hi)
    real(dp), intent(in) :: cuts(:)
    integer, intent(in) :: n
    integer, intent(out) :: lo(size(cuts)+1), hi(size(cuts)+1)
    integer :: i, nc
    integer, allocatable :: icut(:)
    nc = size(cuts) + 1
    allocate(icut(size(cuts)))
    do i = 1, size(cuts)
      icut(i) = nint(cuts(i))
    end do
    lo(1) = 0
    if (size(cuts) > 0) hi(1) = icut(1) - 1
    do i = 2, nc - 1
      lo(i) = icut(i - 1)
      hi(i) = icut(i) - 1
    end do
    if (nc > 1) then
      lo(nc) = icut(size(cuts))
      hi(nc) = n
    else
      hi(1) = n
    end if
  end subroutine make_hb_ranges

  subroutine finish_classification(camat, ccmat, par, out)
    real(dp), intent(in) :: camat(:,:), ccmat(:,:)
    type(beta_params), intent(in) :: par
    type(classification_result), intent(out) :: out
    integer :: nc, i
    real(dp) :: total_a, total_c, tp, tn, fp, fn, rowp, colp

    nc = size(camat,1)
    out%parameters = par
    allocate(out%accuracy_matrix(nc,nc), out%consistency_matrix(nc,nc))
    allocate(out%category_accuracy(nc), out%category_consistency(nc))
    total_a = sum(camat)
    total_c = sum(ccmat)
    if (total_a > 0.0_dp) then
      out%accuracy_matrix = camat / total_a
    else
      out%accuracy_matrix = camat
    end if
    if (total_c > 0.0_dp) then
      out%consistency_matrix = ccmat / total_c
    else
      out%consistency_matrix = ccmat
    end if
    out%overall_accuracy = 0.0_dp
    do i = 1, nc
      out%overall_accuracy = out%overall_accuracy + out%accuracy_matrix(i,i)
      tp = out%accuracy_matrix(i,i)
      fn = sum(out%accuracy_matrix(:,i)) - tp
      fp = sum(out%accuracy_matrix(i,:)) - tp
      tn = 1.0_dp - tp - fn - fp
      out%category_accuracy(i) = ca_stats(tp, tn, fp, fn)
    end do

    out%overall_consistency%p = 0.0_dp
    out%overall_consistency%p_c = 0.0_dp
    do i = 1, nc
      out%overall_consistency%p = out%overall_consistency%p + out%consistency_matrix(i,i)
      colp = sum(out%consistency_matrix(:,i))
      out%overall_consistency%p_c = out%overall_consistency%p_c + colp**2
    end do
    out%overall_consistency%kappa = safe_ratio(out%overall_consistency%p - out%overall_consistency%p_c, &
                                               1.0_dp - out%overall_consistency%p_c)
    do i = 1, nc
      rowp = sum(out%consistency_matrix(i,:))
      out%category_consistency(i)%p = out%consistency_matrix(i,i)
      out%category_consistency(i)%p_c = rowp**2
      out%category_consistency(i)%kappa = safe_ratio(out%category_consistency(i)%p - out%category_consistency(i)%p_c, &
                                                     1.0_dp - out%category_consistency(i)%p_c)
    end do
  end subroutine finish_classification

  subroutine classify_params(par, cuts, true_cuts, use_hb, out)
    type(beta_params), intent(in) :: par
    real(dp), intent(in) :: cuts(:), true_cuts(:)
    logical, intent(in) :: use_hb
    type(classification_result), intent(out) :: out
    integer :: n, nc, i, j
    integer, allocatable :: lo(:), hi(:)
    real(dp), allocatable :: true_bounds(:), camat(:,:), ccmat(:,:)
    type(class_integrand_ctx) :: ctx

    n = nint(merge(par%n, par%etl, use_hb))
    nc = size(cuts) + 1
    allocate(lo(nc), hi(nc), true_bounds(nc+1), camat(nc,nc), ccmat(nc,nc))
    if (use_hb) then
      call make_hb_ranges(cuts, n, lo, hi)
      true_bounds(1) = 0.0_dp
      do i = 1, size(true_cuts)
        true_bounds(i+1) = true_cuts(i) / real(n, dp)
      end do
      true_bounds(nc+1) = 1.0_dp
    else
      call make_category_ranges(cuts, n, lo, hi)
      true_bounds(1) = 0.0_dp
      do i = 1, size(true_cuts)
        true_bounds(i+1) = true_cuts(i)
      end do
      true_bounds(nc+1) = 1.0_dp
    end if

    camat = 0.0_dp
    ccmat = 0.0_dp
    ctx%par = par
    ctx%n = n
    ctx%use_hb = use_hb

    do i = 1, nc
      do j = 1, nc
        ctx%obs_lo = lo(j)
        ctx%obs_hi = hi(j)
        ctx%product_two = .false.
        camat(j,i) = integrate_gk(class_integrand, ctx, true_bounds(i), true_bounds(i+1), 1.0e-10_dp, 2.0e-8_dp)
      end do
    end do

    do i = 1, nc
      do j = 1, nc
        ctx%obs_lo = lo(j)
        ctx%obs_hi = hi(j)
        ctx%obs2_lo = lo(i)
        ctx%obs2_hi = hi(i)
        ctx%product_two = .true.
        ccmat(j,i) = integrate_gk(class_integrand, ctx, 0.0_dp, 1.0_dp, 1.0e-10_dp, 2.0e-8_dp)
      end do
    end do
    call finish_classification(camat, ccmat, par, out)
  end subroutine classify_params

  subroutine ll_classify_params(par, cuts, out, true_cuts)
    type(beta_params), intent(in) :: par
    real(dp), intent(in) :: cuts(:)
    type(classification_result), intent(out) :: out
    real(dp), intent(in), optional :: true_cuts(:)
    if (present(true_cuts)) then
      call classify_params(par, cuts, true_cuts, .false., out)
    else
      call classify_params(par, cuts, cuts, .false., out)
    end if
  end subroutine ll_classify_params

  subroutine ll_classify(scores, reliability, cuts, min_value, max_value, four_parameter, failsafe, l, u, out, true_cuts)
    real(dp), intent(in) :: scores(:), reliability, cuts(:), min_value, max_value, l, u
    logical, intent(in) :: four_parameter, failsafe
    type(classification_result), intent(out) :: out
    real(dp), intent(in), optional :: true_cuts(:)
    type(beta_params) :: par
    real(dp), allocatable :: scaled_cuts(:), scaled_true(:)

    par = beta_true_score_fit(scores, min_value, max_value, 0.0_dp, reliability, .false., &
                              four_parameter, failsafe, l, u)
    allocate(scaled_cuts(size(cuts)))
    scaled_cuts = (cuts - min_value) / (max_value - min_value)
    if (present(true_cuts)) then
      allocate(scaled_true(size(true_cuts)))
      scaled_true = (true_cuts - min_value) / (max_value - min_value)
      call ll_classify_params(par, scaled_cuts, out, scaled_true)
    else
      call ll_classify_params(par, scaled_cuts, out)
    end if
  end subroutine ll_classify

  subroutine hb_classify_params(par, cuts, out, true_cuts)
    type(beta_params), intent(in) :: par
    real(dp), intent(in) :: cuts(:)
    type(classification_result), intent(out) :: out
    real(dp), intent(in), optional :: true_cuts(:)
    if (present(true_cuts)) then
      call classify_params(par, cuts, true_cuts, .true., out)
    else
      call classify_params(par, cuts, cuts, .true., out)
    end if
  end subroutine hb_classify_params

  subroutine hb_classify(scores, reliability, cuts, testlength, four_parameter, failsafe, l, u, out, true_cuts)
    real(dp), intent(in) :: scores(:), reliability, cuts(:), l, u
    integer, intent(in) :: testlength
    logical, intent(in) :: four_parameter, failsafe
    type(classification_result), intent(out) :: out
    real(dp), intent(in), optional :: true_cuts(:)
    type(beta_params) :: par
    real(dp) :: k
    k = lords_k(scores, testlength, reliability)
    par = hb_beta_true_score_fit(scores, real(testlength,dp), k, four_parameter, failsafe, l, u)
    if (present(true_cuts)) then
      call hb_classify_params(par, cuts, out, true_cuts)
    else
      call hb_classify_params(par, cuts, out)
    end if
  end subroutine hb_classify

  subroutine ll_roc_params(par, min_value, max_value, truecut, grainsize, out)
    type(beta_params), intent(in) :: par
    real(dp), intent(in) :: min_value, max_value, truecut
    integer, intent(in) :: grainsize
    type(roc_result), intent(out) :: out
    type(classification_result) :: cr
    real(dp) :: cut(1), tcut(1), scaled_cut(1), scaled_true(1)
    integer :: i

    allocate(out%table(grainsize+1,7))
    scaled_true(1) = (truecut - min_value) / (max_value - min_value)
    do i = 1, grainsize + 1
      cut(1) = min_value + real(i - 1, dp) * (max_value - min_value) / real(grainsize, dp)
      scaled_cut(1) = (cut(1) - min_value) / (max_value - min_value)
      tcut(1) = scaled_true(1)
      call ll_classify_params(par, scaled_cut, cr, tcut)
      out%table(i,1) = 1.0_dp - cr%category_accuracy(1)%specificity
      out%table(i,2) = cr%category_accuracy(1)%sensitivity
      out%table(i,3) = cr%category_accuracy(1)%youden_j
      out%table(i,4) = cut(1)
      out%table(i,5) = cr%category_accuracy(1)%accuracy
      out%table(i,6) = cr%category_accuracy(1)%ppv
      out%table(i,7) = cr%category_accuracy(1)%npv
      if (ieee_is_nan(out%table(i,1))) out%table(i,1) = 0.0_dp
      if (ieee_is_nan(out%table(i,2))) out%table(i,2) = 1.0_dp
      if (ieee_is_nan(out%table(i,6))) out%table(i,6) = 1.0_dp
      if (ieee_is_nan(out%table(i,7))) out%table(i,7) = 1.0_dp
    end do
    out%area = auc(out%table(:,1), out%table(:,2))
    out%max_youden_index = maxloc(out%table(:,3), dim=1)
    out%max_accuracy_index = maxloc(out%table(:,5), dim=1)
  end subroutine ll_roc_params

  subroutine ll_roc(scores, reliability, min_value, max_value, truecut, four_parameter, failsafe, l, u, grainsize, out)
    real(dp), intent(in) :: scores(:), reliability, min_value, max_value, truecut, l, u
    logical, intent(in) :: four_parameter, failsafe
    integer, intent(in) :: grainsize
    type(roc_result), intent(out) :: out
    type(beta_params) :: par
    par = beta_true_score_fit(scores, min_value, max_value, 0.0_dp, reliability, .false., &
                              four_parameter, failsafe, l, u)
    call ll_roc_params(par, min_value, max_value, truecut, grainsize, out)
  end subroutine ll_roc

  subroutine hb_roc_params(par, truecut, grainsize, out)
    type(beta_params), intent(in) :: par
    real(dp), intent(in) :: truecut
    integer, intent(in) :: grainsize
    type(roc_result), intent(out) :: out
    type(classification_result) :: cr
    real(dp) :: cut(1), tcut(1)
    integer :: i, n
    n = nint(par%n)
    allocate(out%table(grainsize+1,7))
    tcut(1) = truecut
    do i = 1, grainsize + 1
      cut(1) = real(i - 1, dp) * real(n, dp) / real(grainsize, dp)
      call hb_classify_params(par, cut, cr, tcut)
      out%table(i,1) = 1.0_dp - cr%category_accuracy(1)%specificity
      out%table(i,2) = cr%category_accuracy(1)%sensitivity
      out%table(i,3) = cr%category_accuracy(1)%youden_j
      out%table(i,4) = cut(1)
      out%table(i,5) = cr%category_accuracy(1)%accuracy
      out%table(i,6) = cr%category_accuracy(1)%ppv
      out%table(i,7) = cr%category_accuracy(1)%npv
      if (ieee_is_nan(out%table(i,1))) out%table(i,1) = 0.0_dp
      if (ieee_is_nan(out%table(i,2))) out%table(i,2) = 1.0_dp
      if (ieee_is_nan(out%table(i,6))) out%table(i,6) = 1.0_dp
      if (ieee_is_nan(out%table(i,7))) out%table(i,7) = 1.0_dp
    end do
    out%area = auc(out%table(:,1), out%table(:,2))
    out%max_youden_index = maxloc(out%table(:,3), dim=1)
    out%max_accuracy_index = maxloc(out%table(:,5), dim=1)
  end subroutine hb_roc_params

  subroutine hb_roc(scores, reliability, testlength, truecut, four_parameter, failsafe, l, u, grainsize, out)
    real(dp), intent(in) :: scores(:), reliability, truecut, l, u
    integer, intent(in) :: testlength, grainsize
    logical, intent(in) :: four_parameter, failsafe
    type(roc_result), intent(out) :: out
    type(beta_params) :: par
    real(dp) :: k
    k = lords_k(scores, testlength, reliability)
    par = hb_beta_true_score_fit(scores, real(testlength,dp), k, four_parameter, failsafe, l, u)
    call hb_roc_params(par, truecut, grainsize, out)
  end subroutine hb_roc

  real(dp) function fit_integrand(x, anyctx) result(v)
    real(dp), intent(in) :: x
    class(*), intent(in) :: anyctx
    integer :: j
    real(dp) :: p
    select type (ctx => anyctx)
    type is (fit_integrand_ctx)
      if (.not. ctx%use_hb) then
        p = binomial_prob_range(ctx%lo, ctx%hi, ctx%n, x)
      else
        p = 0.0_dp
        do j = ctx%lo, ctx%hi
          p = p + compound_binomial_pmf(j, ctx%n, ctx%par%k, x)
        end do
      end if
      v = beta4_pdf(x, ctx%par%l, ctx%par%u, ctx%par%alpha, ctx%par%beta) * p
    class default
      v = 0.0_dp
    end select
  end function fit_integrand

  subroutine merge_small_bins(expected, observed, min_expected)
    real(dp), allocatable, intent(inout) :: expected(:), observed(:)
    real(dp), intent(in) :: min_expected
    real(dp), allocatable :: e2(:), o2(:)
    integer :: i, n

    n = size(expected)
    i = 1
    do while (i < n)
      if (expected(i) < min_expected) then
        expected(i+1) = expected(i+1) + expected(i)
        observed(i+1) = observed(i+1) + observed(i)
        if (n - 1 > 0) then
          allocate(e2(n-1), o2(n-1))
          if (i > 1) then
            e2(1:i-1) = expected(1:i-1)
            o2(1:i-1) = observed(1:i-1)
          end if
          e2(i:n-1) = expected(i+1:n)
          o2(i:n-1) = observed(i+1:n)
          call move_alloc(e2, expected)
          call move_alloc(o2, observed)
          n = n - 1
        end if
      else
        i = i + 1
      end if
    end do
    if (n > 1 .and. expected(n) < min_expected) then
      expected(n-1) = expected(n-1) + expected(n)
      observed(n-1) = observed(n-1) + observed(n)
      allocate(e2(n-1), o2(n-1))
      e2 = expected(1:n-1)
      o2 = observed(1:n-1)
      call move_alloc(e2, expected)
      call move_alloc(o2, observed)
    end if
  end subroutine merge_small_bins

  subroutine ll_model_fit(scores, par, initial_bins, min_expected, four_parameter, out)
    real(dp), intent(in) :: scores(:)
    type(beta_params), intent(in) :: par
    integer, intent(in) :: initial_bins
    real(dp), intent(in) :: min_expected
    logical, intent(in) :: four_parameter
    type(model_fit_result), intent(out) :: out
    integer :: n, b, i, lo, hi, idx
    real(dp) :: width, scaled
    real(dp), allocatable :: expected(:), observed(:)
    type(fit_integrand_ctx) :: ctx

    n = nint(par%etl)
    b = max(initial_bins, 1)
    allocate(expected(b), observed(b))
    expected = 0.0_dp
    observed = 0.0_dp
    width = real(n, dp) / real(b, dp)
    ctx%par = par
    ctx%n = n
    ctx%use_hb = .false.
    do i = 1, b
      lo = ceiling(real(i-1,dp) * width)
      if (i == b) then
        hi = n
      else
        hi = ceiling(real(i,dp) * width) - 1
      end if
      ctx%lo = lo
      ctx%hi = hi
      expected(i) = integrate_gk(fit_integrand, ctx, 0.0_dp, 1.0_dp, 1.0e-10_dp, 2.0e-8_dp) * real(size(scores),dp)
    end do
    do i = 1, size(scores)
      scaled = scores(i)
      idx = min(b, max(1, int(floor(scaled / max(width, tiny(1.0_dp)))) + 1))
      observed(idx) = observed(idx) + 1.0_dp
    end do
    call merge_small_bins(expected, observed, min_expected)
    out%chi_square = 0.0_dp
    do i = 1, size(expected)
      if (expected(i) > 0.0_dp) out%chi_square = out%chi_square + (observed(i)-expected(i))**2 / expected(i)
    end do
    out%df = size(expected) - merge(4, 2, four_parameter .and. .not. par%used_failsafe)
    out%p_value = chi_square_sf(out%chi_square, out%df)
    call move_alloc(expected, out%expected)
    call move_alloc(observed, out%observed)
  end subroutine ll_model_fit

  subroutine hb_model_fit(scores, par, min_expected, four_parameter, out)
    real(dp), intent(in) :: scores(:)
    type(beta_params), intent(in) :: par
    real(dp), intent(in) :: min_expected
    logical, intent(in) :: four_parameter
    type(model_fit_result), intent(out) :: out
    integer :: n, i, idx
    real(dp), allocatable :: expected(:), observed(:)
    type(fit_integrand_ctx) :: ctx

    n = nint(par%n)
    allocate(expected(n+1), observed(n+1))
    expected = 0.0_dp
    observed = 0.0_dp
    ctx%par = par
    ctx%n = n
    ctx%use_hb = .true.
    do i = 0, n
      ctx%lo = i
      ctx%hi = i
      expected(i+1) = integrate_gk(fit_integrand, ctx, 0.0_dp, 1.0_dp, 1.0e-10_dp, 2.0e-8_dp) * real(size(scores),dp)
    end do
    do i = 1, size(scores)
      idx = nint(scores(i)) + 1
      if (idx >= 1 .and. idx <= n+1) observed(idx) = observed(idx) + 1.0_dp
    end do
    call merge_small_bins(expected, observed, min_expected)
    out%chi_square = 0.0_dp
    do i = 1, size(expected)
      if (expected(i) > 0.0_dp) out%chi_square = out%chi_square + (observed(i)-expected(i))**2 / expected(i)
    end do
    out%df = size(expected) - merge(4, 2, four_parameter .and. .not. par%used_failsafe)
    out%p_value = chi_square_sf(out%chi_square, out%df)
    call move_alloc(expected, out%expected)
    call move_alloc(observed, out%observed)
  end subroutine hb_model_fit

end module bf_classification
