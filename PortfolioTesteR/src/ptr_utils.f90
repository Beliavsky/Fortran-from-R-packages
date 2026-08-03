! SPDX-License-Identifier: MIT
! PortfolioTesteR modern Fortran translation
module ptr_utils
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
  use ptr_kinds, only : dp
  implicit none
  private
  public :: nan_dp, is_finite, safe_divide, finite_mean, finite_sd, finite_sum
  public :: covariance_matrix, correlation_matrix, solve_linear, matrix_inverse
  public :: symmetric_eigen_jacobi, rank_vector, percentile, normalize_nonnegative
  public :: sort_indices, median_value, argsort_descending, clamp

contains

  pure real(dp) function nan_dp()
    nan_dp = ieee_value(0.0_dp, ieee_quiet_nan)
  end function nan_dp

  elemental logical function is_finite(x)
    real(dp), intent(in) :: x
    is_finite = ieee_is_finite(x)
  end function is_finite

  elemental real(dp) function safe_divide(a, b, fallback)
    real(dp), intent(in) :: a, b
    real(dp), intent(in), optional :: fallback
    real(dp) :: fb
    fb = nan_dp()
    if (present(fallback)) fb = fallback
    if (is_finite(a) .and. is_finite(b) .and. abs(b) > tiny(1.0_dp)) then
      safe_divide = a / b
    else
      safe_divide = fb
    end if
  end function safe_divide

  pure real(dp) function clamp(x, lo, hi)
    real(dp), intent(in) :: x, lo, hi
    clamp = min(hi, max(lo, x))
  end function clamp

  real(dp) function finite_sum(x)
    real(dp), intent(in) :: x(:)
    integer :: i
    finite_sum = 0.0_dp
    do i = 1, size(x)
      if (is_finite(x(i))) finite_sum = finite_sum + x(i)
    end do
  end function finite_sum

  real(dp) function finite_mean(x)
    real(dp), intent(in) :: x(:)
    integer :: i, n
    finite_mean = 0.0_dp
    n = 0
    do i = 1, size(x)
      if (is_finite(x(i))) then
        finite_mean = finite_mean + x(i)
        n = n + 1
      end if
    end do
    if (n > 0) then
      finite_mean = finite_mean / real(n, dp)
    else
      finite_mean = nan_dp()
    end if
  end function finite_mean

  real(dp) function finite_sd(x, sample)
    real(dp), intent(in) :: x(:)
    logical, intent(in), optional :: sample
    logical :: use_sample
    real(dp) :: mu, ss
    integer :: i, n, denom
    use_sample = .true.
    if (present(sample)) use_sample = sample
    mu = finite_mean(x)
    if (.not. is_finite(mu)) then
      finite_sd = nan_dp()
      return
    end if
    ss = 0.0_dp
    n = 0
    do i = 1, size(x)
      if (is_finite(x(i))) then
        ss = ss + (x(i) - mu)**2
        n = n + 1
      end if
    end do
    denom = n
    if (use_sample) denom = n - 1
    if (denom > 0) then
      finite_sd = sqrt(max(0.0_dp, ss / real(denom, dp)))
    else
      finite_sd = nan_dp()
    end if
  end function finite_sd

  subroutine covariance_matrix(x, cov, means, status)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable, intent(out) :: cov(:,:)
    real(dp), allocatable, intent(out), optional :: means(:)
    integer, intent(out), optional :: status
    integer :: n, p, i, j, k, cnt
    real(dp), allocatable :: mu(:)
    real(dp) :: s
    n = size(x, 1); p = size(x, 2)
    allocate(cov(p,p), mu(p))
    do j = 1, p
      mu(j) = finite_mean(x(:,j))
    end do
    cov = 0.0_dp
    do j = 1, p
      do k = j, p
        s = 0.0_dp; cnt = 0
        do i = 1, n
          if (is_finite(x(i,j)) .and. is_finite(x(i,k))) then
            s = s + (x(i,j) - mu(j)) * (x(i,k) - mu(k))
            cnt = cnt + 1
          end if
        end do
        if (cnt > 1) then
          cov(j,k) = s / real(cnt - 1, dp)
        else
          cov(j,k) = 0.0_dp
        end if
        cov(k,j) = cov(j,k)
      end do
    end do
    if (present(means)) then
      allocate(means(p)); means = mu
    end if
    if (present(status)) status = 0
  end subroutine covariance_matrix

  subroutine correlation_matrix(x, corr)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable, intent(out) :: corr(:,:)
    real(dp), allocatable :: cov(:,:)
    real(dp) :: den
    integer :: i, j, p
    call covariance_matrix(x, cov)
    p = size(cov,1)
    allocate(corr(p,p)); corr = 0.0_dp
    do i = 1, p
      corr(i,i) = 1.0_dp
      do j = i + 1, p
        den = sqrt(max(0.0_dp, cov(i,i) * cov(j,j)))
        if (den > tiny(1.0_dp)) corr(i,j) = cov(i,j) / den
        corr(j,i) = corr(i,j)
      end do
    end do
  end subroutine correlation_matrix

  subroutine solve_linear(a, b, x, status)
    real(dp), intent(in) :: a(:,:), b(:)
    real(dp), allocatable, intent(out) :: x(:)
    integer, intent(out) :: status
    real(dp), allocatable :: aug(:,:)
    real(dp) :: pivot, factor, tmp
    integer :: n, i, j, k, piv
    n = size(b)
    status = 0
    allocate(x(n), aug(n,n+1))
    if (size(a,1) /= n .or. size(a,2) /= n) then
      status = 1; x = 0.0_dp; return
    end if
    aug(:,1:n) = a; aug(:,n+1) = b
    do k = 1, n
      piv = k
      do i = k + 1, n
        if (abs(aug(i,k)) > abs(aug(piv,k))) piv = i
      end do
      if (abs(aug(piv,k)) <= 100.0_dp * epsilon(1.0_dp)) then
        status = 2; x = 0.0_dp; return
      end if
      if (piv /= k) then
        do j = k, n + 1
          tmp = aug(k,j); aug(k,j) = aug(piv,j); aug(piv,j) = tmp
        end do
      end if
      pivot = aug(k,k)
      aug(k,k:n+1) = aug(k,k:n+1) / pivot
      do i = 1, n
        if (i == k) cycle
        factor = aug(i,k)
        aug(i,k:n+1) = aug(i,k:n+1) - factor * aug(k,k:n+1)
      end do
    end do
    x = aug(:,n+1)
  end subroutine solve_linear

  subroutine matrix_inverse(a, ainv, status)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: ainv(:,:)
    integer, intent(out) :: status
    real(dp), allocatable :: e(:), col(:)
    integer :: n, j, st
    n = size(a,1)
    allocate(ainv(n,n), e(n)); ainv = 0.0_dp
    status = 0
    do j = 1, n
      e = 0.0_dp; e(j) = 1.0_dp
      call solve_linear(a, e, col, st)
      if (st /= 0) then
        status = st; ainv = 0.0_dp; return
      end if
      ainv(:,j) = col
    end do
  end subroutine matrix_inverse

  subroutine symmetric_eigen_jacobi(a, values, vectors, status, tol, max_iter)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: values(:), vectors(:,:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: tol
    integer, intent(in), optional :: max_iter
    real(dp), allocatable :: d(:,:)
    real(dp) :: threshold, app, aqq, apq, tau, t, c, s, tmp, maxoff
    integer :: n, p, q, i, j, iter, itmax
    n = size(a,1)
    allocate(d(n,n), values(n), vectors(n,n))
    d = 0.5_dp * (a + transpose(a)); vectors = 0.0_dp
    do i = 1, n; vectors(i,i) = 1.0_dp; end do
    threshold = 1.0e-12_dp; if (present(tol)) threshold = tol
    itmax = max(100, 50*n*n); if (present(max_iter)) itmax = max_iter
    status = 1
    do iter = 1, itmax
      maxoff = 0.0_dp; p = 1; q = min(2,n)
      do i = 1, n - 1
        do j = i + 1, n
          if (abs(d(i,j)) > maxoff) then
            maxoff = abs(d(i,j)); p = i; q = j
          end if
        end do
      end do
      if (maxoff <= threshold) then
        status = 0; exit
      end if
      app = d(p,p); aqq = d(q,q); apq = d(p,q)
      tau = (aqq - app) / (2.0_dp * apq)
      if (tau >= 0.0_dp) then
        t = 1.0_dp / (tau + sqrt(1.0_dp + tau*tau))
      else
        t = -1.0_dp / (-tau + sqrt(1.0_dp + tau*tau))
      end if
      c = 1.0_dp / sqrt(1.0_dp + t*t); s = t*c
      do i = 1, n
        if (i /= p .and. i /= q) then
          tmp = d(i,p)
          d(i,p) = c*tmp - s*d(i,q); d(p,i) = d(i,p)
          d(i,q) = s*tmp + c*d(i,q); d(q,i) = d(i,q)
        end if
      end do
      d(p,p) = c*c*app - 2.0_dp*s*c*apq + s*s*aqq
      d(q,q) = s*s*app + 2.0_dp*s*c*apq + c*c*aqq
      d(p,q) = 0.0_dp; d(q,p) = 0.0_dp
      do i = 1, n
        tmp = vectors(i,p)
        vectors(i,p) = c*tmp - s*vectors(i,q)
        vectors(i,q) = s*tmp + c*vectors(i,q)
      end do
    end do
    do i = 1, n; values(i) = d(i,i); end do
    call sort_eigenpairs(values, vectors)
  contains
    subroutine sort_eigenpairs(vals, vecs)
      real(dp), intent(inout) :: vals(:), vecs(:,:)
      real(dp) :: tv
      real(dp), allocatable :: tmpv(:)
      integer :: ii, jj, m
      allocate(tmpv(size(vals)))
      do ii = 1, size(vals)-1
        m = ii
        do jj = ii+1, size(vals)
          if (vals(jj) < vals(m)) m = jj
        end do
        if (m /= ii) then
          tv = vals(ii); vals(ii) = vals(m); vals(m) = tv
          tmpv = vecs(:,ii); vecs(:,ii) = vecs(:,m); vecs(:,m) = tmpv
        end if
      end do
    end subroutine sort_eigenpairs
  end subroutine symmetric_eigen_jacobi

  subroutine sort_indices(x, idx, ascending)
    real(dp), intent(in) :: x(:)
    integer, allocatable, intent(out) :: idx(:)
    logical, intent(in), optional :: ascending
    logical :: asc
    integer :: i, j, key
    asc = .true.; if (present(ascending)) asc = ascending
    allocate(idx(size(x))); idx = [(i, i=1,size(x))]
    do i = 2, size(idx)
      key = idx(i); j = i - 1
      do while (j >= 1)
        if (asc) then
          if (x(idx(j)) <= x(key)) exit
        else
          if (x(idx(j)) >= x(key)) exit
        end if
        idx(j+1) = idx(j); j = j - 1
      end do
      idx(j+1) = key
    end do
  end subroutine sort_indices

  subroutine argsort_descending(x, idx)
    real(dp), intent(in) :: x(:)
    integer, allocatable, intent(out) :: idx(:)
    call sort_indices(x, idx, .false.)
  end subroutine argsort_descending

  subroutine rank_vector(x, ranks, descending, normalize)
    real(dp), intent(in) :: x(:)
    real(dp), allocatable, intent(out) :: ranks(:)
    logical, intent(in), optional :: descending, normalize
    logical :: desc, norm
    integer, allocatable :: idx(:)
    integer :: i, j, k, nvalid
    real(dp) :: avg
    desc = .false.; if (present(descending)) desc = descending
    norm = .false.; if (present(normalize)) norm = normalize
    allocate(ranks(size(x))); ranks = nan_dp()
    nvalid = count(is_finite(x))
    if (nvalid == 0) return
    allocate(idx(nvalid)); k = 0
    do i = 1, size(x)
      if (is_finite(x(i))) then; k = k + 1; idx(k) = i; end if
    end do
    call insertion_sort_valid(x, idx, desc)
    i = 1
    do while (i <= nvalid)
      j = i
      do while (j < nvalid)
        if (abs(x(idx(j+1)) - x(idx(i))) > &
            10.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(x(idx(i))))) exit
        j = j + 1
      end do
      avg = 0.5_dp * real(i + j, dp)
      do k = i, j; ranks(idx(k)) = avg; end do
      i = j + 1
    end do
    if (norm .and. nvalid > 1) then
      do i = 1, size(ranks)
        if (is_finite(ranks(i))) ranks(i) = (ranks(i) - 1.0_dp) / real(nvalid - 1, dp)
      end do
    else if (norm .and. nvalid == 1) then
      where (is_finite(ranks)) ranks = 1.0_dp
    end if
  contains
    subroutine insertion_sort_valid(v, id, dsc)
      real(dp), intent(in) :: v(:)
      integer, intent(inout) :: id(:)
      logical, intent(in) :: dsc
      integer :: aa, bb, kk
      do aa = 2, size(id)
        kk = id(aa); bb = aa - 1
        do while (bb >= 1)
          if (.not. dsc) then
            if (v(id(bb)) <= v(kk)) exit
          else
            if (v(id(bb)) >= v(kk)) exit
          end if
          id(bb+1) = id(bb); bb = bb - 1
        end do
        id(bb+1) = kk
      end do
    end subroutine insertion_sort_valid
  end subroutine rank_vector

  real(dp) function percentile(x, prob)
    real(dp), intent(in) :: x(:), prob
    real(dp), allocatable :: v(:)
    integer, allocatable :: idx(:)
    real(dp) :: pos, frac
    integer :: i, n, lo, hi
    n = count(is_finite(x))
    if (n == 0) then; percentile = nan_dp(); return; end if
    allocate(v(n)); n = 0
    do i = 1, size(x)
      if (is_finite(x(i))) then; n = n + 1; v(n) = x(i); end if
    end do
    call sort_indices(v, idx, .true.)
    pos = 1.0_dp + clamp(prob,0.0_dp,1.0_dp) * real(n-1,dp)
    lo = floor(pos); hi = ceiling(pos); frac = pos - real(lo,dp)
    percentile = (1.0_dp-frac)*v(idx(lo)) + frac*v(idx(hi))
  end function percentile

  real(dp) function median_value(x)
    real(dp), intent(in) :: x(:)
    median_value = percentile(x, 0.5_dp)
  end function median_value

  subroutine normalize_nonnegative(x)
    real(dp), intent(inout) :: x(:)
    real(dp) :: s
    where (.not. is_finite(x) .or. x < 0.0_dp) x = 0.0_dp
    s = sum(x)
    if (s > tiny(1.0_dp)) x = x / s
  end subroutine normalize_nonnegative

end module ptr_utils
