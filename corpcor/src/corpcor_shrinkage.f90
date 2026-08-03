! SPDX-License-Identifier: GPL-3.0-or-later
module corpcor_shrinkage
  use corpcor_kinds, only : dp
  use corpcor_types, only : matrix_shrinkage_result, vector_shrinkage_result, svd_result, scale_result, &
    corpcor_success, corpcor_invalid_argument, corpcor_dimension_error
  use corpcor_weighted, only : normalized_weights, weighted_moments, weighted_scale, median_value
  use corpcor_linalg, only : matrix_power, pseudoinverse, covariance_to_correlation, fast_svd, identity_matrix
  use corpcor_matrix_tools, only : rebuild_covariance
  implicit none
  private
  public :: estimate_lambda, estimate_lambda_var
  public :: variance_shrinkage, correlation_power_shrinkage
  public :: correlation_shrinkage, inverse_correlation_shrinkage
  public :: covariance_shrinkage, inverse_covariance_shrinkage
  public :: partial_correlation_shrinkage, partial_variance_shrinkage
  public :: correlation_to_partial, partial_to_correlation
  public :: crossprod_correlation_power_shrinkage

contains

  function estimate_lambda_var(x, w, status) result(lambda_var)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(in), optional :: w(:)
    integer, intent(out), optional :: status
    real(dp) :: lambda_var
    real(dp), allocatable :: wn(:), xc(:, :), v(:), q1(:), q2(:)
    real(dp) :: w2, h1, h1w2, target, numerator, denominator
    type_scale: block
      use corpcor_types, only : scale_result
      type(scale_result) :: sc
      integer :: n, p, j, istat
      n = size(x, 1)
      p = size(x, 2)
      lambda_var = 0.0_dp
      if (n < 3 .or. p < 1) then
        if (present(status)) status = corpcor_invalid_argument
        exit type_scale
      end if
      wn = normalized_weights(n, w, istat)
      if (istat /= corpcor_success) then
        if (present(status)) status = istat
        exit type_scale
      end if
      w2 = sum(wn * wn)
      if (w2 >= 1.0_dp) then
        if (present(status)) status = corpcor_invalid_argument
        exit type_scale
      end if
      h1 = 1.0_dp / (1.0_dp - w2)
      h1w2 = w2 / (1.0_dp - w2)
      sc = weighted_scale(x, wn, center=.true., scale=.false.)
      xc = sc%x
      allocate(v(p), q1(p), q2(p))
      do j = 1, p
        v(j) = h1 * sum(wn * xc(:, j) ** 2)
        q1(j) = sum(wn * xc(:, j) ** 2)
        q2(j) = sum(wn * xc(:, j) ** 4) - q1(j) ** 2
      end do
      target = median_value(v)
      numerator = sum(q2)
      denominator = sum((q1 - target / h1) ** 2)
      if (abs(denominator) <= epsilon(1.0_dp)) then
        lambda_var = 1.0_dp
      else
        lambda_var = max(0.0_dp, min(1.0_dp, numerator / denominator * h1w2))
      end if
      if (present(status)) status = corpcor_success
    end block type_scale
  end function estimate_lambda_var

  function estimate_lambda(x, w, status) result(lambda)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(in), optional :: w(:)
    integer, intent(out), optional :: status
    real(dp) :: lambda
    real(dp), allocatable :: wn(:), xsw(:, :), xs2w(:, :), gram(:, :), column_norm2(:)
    real(dp) :: w2, h1w2, s_e2r, s_er2, numerator, denominator
    integer :: n, p, i, istat
    type(scale_result) :: sc

    n = size(x, 1)
    p = size(x, 2)
    lambda = 0.0_dp
    if (p == 1) then
      lambda = 1.0_dp
      if (present(status)) status = corpcor_success
      return
    end if
    if (n < 3 .or. p < 1) then
      if (present(status)) status = corpcor_invalid_argument
      return
    end if
    wn = normalized_weights(n, w, istat)
    if (istat /= corpcor_success) then
      if (present(status)) status = istat
      return
    end if
    sc = weighted_scale(x, wn, center=.true., scale=.true.)
    if (sc%status /= corpcor_success) then
      if (present(status)) status = sc%status
      return
    end if
    w2 = sum(wn * wn)
    h1w2 = w2 / (1.0_dp - w2)
    allocate(xsw(n, p), xs2w(n, p), column_norm2(p))
    do i = 1, n
      xsw(i, :) = sqrt(wn(i)) * sc%x(i, :)
      xs2w(i, :) = sqrt(wn(i)) * sc%x(i, :) ** 2
    end do

    ! ||X'X||_F = ||XX'||_F.  Use the smaller Gram matrix and remove
    ! diagonal terms, retaining the upstream O(min(n,p)^2 + np) memory profile.
    if (n <= p) then
      gram = matmul(xsw, transpose(xsw))
    else
      gram = matmul(transpose(xsw), xsw)
    end if
    column_norm2 = sum(xsw * xsw, dim=1)
    s_e2r = sum(gram * gram) - sum(column_norm2 * column_norm2)
    s_er2 = sum(sum(xs2w, dim=2) ** 2 - sum(xs2w * xs2w, dim=2))

    denominator = s_e2r
    numerator = s_er2 - s_e2r
    if (abs(denominator) <= epsilon(1.0_dp)) then
      lambda = 1.0_dp
    else
      lambda = max(0.0_dp, min(1.0_dp, numerator / denominator * h1w2))
    end if
    if (present(status)) status = corpcor_success
  end function estimate_lambda

  function variance_shrinkage(x, lambda_var, w) result(res)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(in), optional :: lambda_var
    real(dp), intent(in), optional :: w(:)
    type(vector_shrinkage_result) :: res
    type_mom: block
      use corpcor_types, only : moments_result
      type(moments_result) :: mom
      real(dp) :: target
      integer :: istat
      if (present(lambda_var)) then
        res%lambda_var = max(0.0_dp, min(1.0_dp, lambda_var))
        res%lambda_var_estimated = .false.
      else
        res%lambda_var = estimate_lambda_var(x, w, istat)
        res%lambda_var_estimated = .true.
        if (istat /= corpcor_success) then
          allocate(res%value(0))
          res%status = istat
          exit type_mom
        end if
      end if
      mom = weighted_moments(x, w)
      if (mom%status /= corpcor_success) then
        allocate(res%value(0))
        res%status = mom%status
        exit type_mom
      end if
      target = median_value(mom%variance)
      res%value = res%lambda_var * target + (1.0_dp - res%lambda_var) * mom%variance
      res%status = corpcor_success
    end block type_mom
  end function variance_shrinkage

  function empirical_correlation(x, w, status) result(r)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(in), optional :: w(:)
    integer, intent(out) :: status
    real(dp), allocatable :: r(:, :)
    real(dp), allocatable :: wn(:), xsw(:, :)
    real(dp) :: h1
    integer :: n, p, i, j, istat
    type_scale: block
      use corpcor_types, only : scale_result
      type(scale_result) :: sc
      n = size(x, 1)
      p = size(x, 2)
      allocate(r(p, p))
      r = 0.0_dp
      wn = normalized_weights(n, w, istat)
      if (istat /= corpcor_success) then
        status = istat
        exit type_scale
      end if
      sc = weighted_scale(x, wn, center=.true., scale=.true.)
      if (sc%status /= corpcor_success) then
        status = sc%status
        exit type_scale
      end if
      h1 = 1.0_dp / (1.0_dp - sum(wn * wn))
      allocate(xsw(n, p))
      do i = 1, n
        xsw(i, :) = sqrt(wn(i)) * sc%x(i, :)
      end do
      r = h1 * matmul(transpose(xsw), xsw)
      do j = 1, p
        r(j, j) = 1.0_dp
        if (sc%zero_scale(j)) then
          r(j, :) = 0.0_dp
          r(:, j) = 0.0_dp
          r(j, j) = 1.0_dp
        end if
      end do
      r = 0.5_dp * (r + transpose(r))
      status = corpcor_success
    end block type_scale
  end function empirical_correlation

  function correlation_power_shrinkage(x, alpha, lambda, w) result(res)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(in) :: alpha
    real(dp), intent(in), optional :: lambda
    real(dp), intent(in), optional :: w(:)
    type(matrix_shrinkage_result) :: res
    real(dp), allocatable :: r0(:, :), wn(:), utwu(:, :), cmat(:, :)
    real(dp), allocatable :: fmat(:, :), tmp(:, :)
    real(dp) :: h1
    integer :: n, p, m, i, j, istat
    type(scale_result) :: sc
    type(svd_result) :: s

    n = size(x, 1)
    p = size(x, 2)
    if (present(lambda)) then
      res%lambda = max(0.0_dp, min(1.0_dp, lambda))
      res%lambda_estimated = .false.
    else
      res%lambda = estimate_lambda(x, w, istat)
      res%lambda_estimated = .true.
      if (istat /= corpcor_success) then
        allocate(res%value(0, 0))
        res%status = istat
        return
      end if
    end if
    if (abs(res%lambda - 1.0_dp) <= epsilon(1.0_dp) .or. abs(alpha) <= epsilon(1.0_dp)) then
      res%value = identity_matrix(p)
      res%status = corpcor_success
      return
    end if
    if (abs(alpha - 1.0_dp) <= epsilon(1.0_dp)) then
      r0 = empirical_correlation(x, w, istat)
      if (istat /= corpcor_success) then
        allocate(res%value(0, 0))
        res%status = istat
        return
      end if
      res%value = (1.0_dp - res%lambda) * r0
      do j = 1, p
        res%value(j, j) = 1.0_dp
      end do
      res%status = corpcor_success
      return
    end if

    wn = normalized_weights(n, w, istat)
    sc = weighted_scale(x, wn, center=.true., scale=.true.)
    if (istat /= corpcor_success .or. sc%status /= corpcor_success) then
      allocate(res%value(0, 0))
      res%status = max(istat, sc%status)
      return
    end if
    h1 = 1.0_dp / (1.0_dp - sum(wn * wn))
    s = fast_svd(sc%x)
    m = s%rank
    if (m == 0) then
      res%value = identity_matrix(p)
      res%status = corpcor_success
      return
    end if
    allocate(utwu(m, m), cmat(m, m))
    utwu = 0.0_dp
    do i = 1, n
      utwu = utwu + wn(i) * spread(s%u(i, :), 2, m) * spread(s%u(i, :), 1, m)
    end do
    cmat = (1.0_dp - res%lambda) * h1 * &
      (spread(s%d, 2, m) * utwu * spread(s%d, 1, m))
    cmat = 0.5_dp * (cmat + transpose(cmat))

    if (res%lambda <= epsilon(1.0_dp)) then
      tmp = matrix_power(cmat, alpha, pseudo=.true., status=istat)
      res%value = matmul(s%v, matmul(tmp, transpose(s%v)))
    else
      tmp = matrix_power(cmat / res%lambda + identity_matrix(m), alpha, status=istat)
      fmat = identity_matrix(m) - tmp
      res%value = (identity_matrix(p) - matmul(s%v, matmul(fmat, transpose(s%v)))) * &
        res%lambda ** alpha
    end if
    if (istat /= corpcor_success) then
      res%status = istat
      return
    end if
    do j = 1, p
      if (sc%zero_scale(j)) then
        res%value(j, :) = 0.0_dp
        res%value(:, j) = 0.0_dp
        res%value(j, j) = 1.0_dp
      end if
    end do
    res%value = 0.5_dp * (res%value + transpose(res%value))
    res%status = corpcor_success
  end function correlation_power_shrinkage

  function correlation_shrinkage(x, lambda, w) result(res)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(in), optional :: lambda
    real(dp), intent(in), optional :: w(:)
    type(matrix_shrinkage_result) :: res
    if (present(lambda)) then
      res = correlation_power_shrinkage(x, 1.0_dp, lambda=lambda, w=w)
    else
      res = correlation_power_shrinkage(x, 1.0_dp, w=w)
    end if
  end function correlation_shrinkage

  function inverse_correlation_shrinkage(x, lambda, w) result(res)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(in), optional :: lambda
    real(dp), intent(in), optional :: w(:)
    type(matrix_shrinkage_result) :: res
    if (present(lambda)) then
      res = correlation_power_shrinkage(x, -1.0_dp, lambda=lambda, w=w)
    else
      res = correlation_power_shrinkage(x, -1.0_dp, w=w)
    end if
  end function inverse_correlation_shrinkage

  function covariance_shrinkage(x, lambda, lambda_var, w) result(res)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(in), optional :: lambda, lambda_var
    real(dp), intent(in), optional :: w(:)
    type(matrix_shrinkage_result) :: res
    type(vector_shrinkage_result) :: vr
    type(matrix_shrinkage_result) :: cr

    if (present(lambda_var)) then
      vr = variance_shrinkage(x, lambda_var=lambda_var, w=w)
    else
      vr = variance_shrinkage(x, w=w)
    end if
    if (present(lambda)) then
      cr = correlation_shrinkage(x, lambda=lambda, w=w)
    else
      cr = correlation_shrinkage(x, w=w)
    end if
    if (vr%status /= corpcor_success .or. cr%status /= corpcor_success) then
      allocate(res%value(0, 0))
      res%status = max(vr%status, cr%status)
      return
    end if
    res%value = rebuild_covariance(cr%value, vr%value, res%status)
    res%lambda = cr%lambda
    res%lambda_estimated = cr%lambda_estimated
    res%lambda_var = vr%lambda_var
    res%lambda_var_estimated = vr%lambda_var_estimated
  end function covariance_shrinkage

  function inverse_covariance_shrinkage(x, lambda, lambda_var, w) result(res)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(in), optional :: lambda, lambda_var
    real(dp), intent(in), optional :: w(:)
    type(matrix_shrinkage_result) :: res
    type(vector_shrinkage_result) :: vr
    type(matrix_shrinkage_result) :: ir
    real(dp), allocatable :: invsd(:)
    integer :: i, j, p

    if (present(lambda_var)) then
      vr = variance_shrinkage(x, lambda_var=lambda_var, w=w)
    else
      vr = variance_shrinkage(x, w=w)
    end if
    if (present(lambda)) then
      ir = inverse_correlation_shrinkage(x, lambda=lambda, w=w)
    else
      ir = inverse_correlation_shrinkage(x, w=w)
    end if
    if (vr%status /= corpcor_success .or. ir%status /= corpcor_success .or. any(vr%value <= 0.0_dp)) then
      allocate(res%value(0, 0))
      res%status = corpcor_invalid_argument
      return
    end if
    p = size(vr%value)
    allocate(res%value(p, p), invsd(p))
    invsd = 1.0_dp / sqrt(vr%value)
    do j = 1, p
      do i = 1, p
        res%value(i, j) = ir%value(i, j) * invsd(i) * invsd(j)
      end do
    end do
    res%lambda = ir%lambda
    res%lambda_estimated = ir%lambda_estimated
    res%lambda_var = vr%lambda_var
    res%lambda_var_estimated = vr%lambda_var_estimated
    res%status = corpcor_success
  end function inverse_covariance_shrinkage

  function correlation_to_partial(correlation, tolerance, status) result(partial)
    real(dp), intent(in) :: correlation(:, :)
    real(dp), intent(in), optional :: tolerance
    integer, intent(out), optional :: status
    real(dp), allocatable :: partial(:, :)
    real(dp), allocatable :: work(:, :)
    integer :: i, n, istat

    if (present(tolerance)) then
      work = -pseudoinverse(correlation, tolerance, istat)
    else
      work = -pseudoinverse(correlation, status=istat)
    end if
    n = size(work, 1)
    do i = 1, n
      work(i, i) = -work(i, i)
    end do
    partial = covariance_to_correlation(work, istat)
    if (present(status)) status = istat
  end function correlation_to_partial

  function partial_to_correlation(partial, tolerance, status) result(correlation)
    real(dp), intent(in) :: partial(:, :)
    real(dp), intent(in), optional :: tolerance
    integer, intent(out), optional :: status
    real(dp), allocatable :: correlation(:, :)
    real(dp), allocatable :: work(:, :)
    integer :: i, n, istat

    work = -partial
    n = size(work, 1)
    do i = 1, n
      work(i, i) = -work(i, i)
    end do
    if (present(tolerance)) then
      work = pseudoinverse(work, tolerance, istat)
    else
      work = pseudoinverse(work, status=istat)
    end if
    correlation = covariance_to_correlation(work, istat)
    if (present(status)) status = istat
  end function partial_to_correlation

  function partial_correlation_shrinkage(x, lambda, w) result(res)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(in), optional :: lambda
    real(dp), intent(in), optional :: w(:)
    type(matrix_shrinkage_result) :: res
    type(matrix_shrinkage_result) :: ir
    real(dp), allocatable :: work(:, :)
    integer :: i, n, istat

    if (present(lambda)) then
      ir = inverse_correlation_shrinkage(x, lambda=lambda, w=w)
    else
      ir = inverse_correlation_shrinkage(x, w=w)
    end if
    if (ir%status /= corpcor_success) then
      allocate(res%value(0, 0))
      res%status = ir%status
      return
    end if
    work = -ir%value
    n = size(work, 1)
    allocate(res%standardized_partial_variance(n))
    do i = 1, n
      work(i, i) = -work(i, i)
      res%standardized_partial_variance(i) = 1.0_dp / work(i, i)
    end do
    res%value = covariance_to_correlation(work, istat)
    res%lambda = ir%lambda
    res%lambda_estimated = ir%lambda_estimated
    res%status = istat
  end function partial_correlation_shrinkage

  function partial_variance_shrinkage(x, lambda, lambda_var, w) result(res)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(in), optional :: lambda, lambda_var
    real(dp), intent(in), optional :: w(:)
    type(vector_shrinkage_result) :: res
    type(matrix_shrinkage_result) :: prec
    integer :: i, p

    if (present(lambda) .and. present(lambda_var)) then
      prec = inverse_covariance_shrinkage(x, lambda=lambda, lambda_var=lambda_var, w=w)
    else if (present(lambda)) then
      prec = inverse_covariance_shrinkage(x, lambda=lambda, w=w)
    else if (present(lambda_var)) then
      prec = inverse_covariance_shrinkage(x, lambda_var=lambda_var, w=w)
    else
      prec = inverse_covariance_shrinkage(x, w=w)
    end if
    if (prec%status /= corpcor_success) then
      allocate(res%value(0))
      res%status = prec%status
      return
    end if
    p = size(prec%value, 1)
    allocate(res%value(p))
    do i = 1, p
      res%value(i) = 1.0_dp / prec%value(i, i)
    end do
    res%lambda_var = prec%lambda_var
    res%lambda_var_estimated = prec%lambda_var_estimated
    res%status = corpcor_success
  end function partial_variance_shrinkage

  function crossprod_correlation_power_shrinkage(x, y, alpha, lambda, w) result(res)
    real(dp), intent(in) :: x(:, :), y(:, :)
    real(dp), intent(in) :: alpha
    real(dp), intent(in), optional :: lambda
    real(dp), intent(in), optional :: w(:)
    type(matrix_shrinkage_result) :: res
    real(dp), allocatable :: wn(:), utwu(:, :), cmat(:, :), fmat(:, :), tmp(:, :), vty(:, :)
    real(dp) :: h1
    integer :: n, p, m, i, j, istat
    type(scale_result) :: sc
    type(svd_result) :: s
    type(matrix_shrinkage_result) :: pr

    n = size(x, 1)
    p = size(x, 2)
    if (size(y, 1) /= p) then
      allocate(res%value(0, 0))
      res%status = corpcor_dimension_error
      return
    end if
    if (present(lambda)) then
      res%lambda = max(0.0_dp, min(1.0_dp, lambda))
      res%lambda_estimated = .false.
    else
      res%lambda = estimate_lambda(x, w, istat)
      res%lambda_estimated = .true.
      if (istat /= corpcor_success) then
        allocate(res%value(0, 0))
        res%status = istat
        return
      end if
    end if
    if (abs(res%lambda - 1.0_dp) <= epsilon(1.0_dp) .or. abs(alpha) <= epsilon(1.0_dp)) then
      res%value = y
      res%status = corpcor_success
      return
    end if
    if (abs(alpha - 1.0_dp) <= epsilon(1.0_dp)) then
      pr = correlation_power_shrinkage(x, alpha, lambda=res%lambda, w=w)
      res%value = matmul(pr%value, y)
      res%status = pr%status
      return
    end if

    wn = normalized_weights(n, w, istat)
    sc = weighted_scale(x, wn, center=.true., scale=.true.)
    if (istat /= corpcor_success .or. sc%status /= corpcor_success) then
      allocate(res%value(0, 0))
      res%status = max(istat, sc%status)
      return
    end if
    h1 = 1.0_dp / (1.0_dp - sum(wn * wn))
    s = fast_svd(sc%x)
    m = s%rank
    if (m == 0) then
      res%value = y
      res%status = corpcor_success
      return
    end if
    allocate(utwu(m, m), cmat(m, m))
    utwu = 0.0_dp
    do i = 1, n
      utwu = utwu + wn(i) * spread(s%u(i, :), 2, m) * spread(s%u(i, :), 1, m)
    end do
    cmat = (1.0_dp - res%lambda) * h1 * &
      (spread(s%d, 2, m) * utwu * spread(s%d, 1, m))
    cmat = 0.5_dp * (cmat + transpose(cmat))
    vty = matmul(transpose(s%v), y)
    if (res%lambda <= epsilon(1.0_dp)) then
      tmp = matrix_power(cmat, alpha, pseudo=.true., status=istat)
      res%value = matmul(s%v, matmul(tmp, vty))
    else
      tmp = matrix_power(cmat / res%lambda + identity_matrix(m), alpha, status=istat)
      fmat = identity_matrix(m) - tmp
      res%value = (y - matmul(s%v, matmul(fmat, vty))) * res%lambda ** alpha
    end if
    if (istat /= corpcor_success) then
      res%status = istat
      return
    end if
    do j = 1, p
      if (sc%zero_scale(j)) res%value(j, :) = y(j, :)
    end do
    res%status = corpcor_success
  end function crossprod_correlation_power_shrinkage

end module corpcor_shrinkage
