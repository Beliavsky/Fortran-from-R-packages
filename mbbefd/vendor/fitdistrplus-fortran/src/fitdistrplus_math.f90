! SPDX-License-Identifier: GPL-2.0-or-later
module fitdistrplus_math
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fitdistrplus_kinds, only : dp
  use fitdistrplus_types, only : fit_success, fit_invalid_argument, fit_singular
  implicit none
  private

  real(dp), parameter, public :: pi_dp = acos(-1.0_dp)
  public :: normal_pdf, normal_cdf, normal_quantile
  public :: regularized_gamma_p, regularized_gamma_q, regularized_beta
  public :: chi_square_survival, type7_quantile, sort_real, sample_mean
  public :: sample_variance, weighted_mean, weighted_variance, weighted_quantile
  public :: invert_matrix, numerical_hessian, covariance_from_hessian
  public :: random_normal, random_gamma, random_poisson, seed_rng
  public :: log1mexp, clamp_probability

contains

  pure elemental function normal_pdf(x) result(value)
    real(dp), intent(in) :: x
    real(dp) :: value
    value = exp(-0.5_dp * x * x) / sqrt(2.0_dp * pi_dp)
  end function normal_pdf

  pure elemental function normal_cdf(x) result(value)
    real(dp), intent(in) :: x
    real(dp) :: value
    value = 0.5_dp * erfc(-x / sqrt(2.0_dp))
  end function normal_cdf

  function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp) :: x, q, r
    real(dp), parameter :: a(6) = [ &
      -3.969683028665376e1_dp, 2.209460984245205e2_dp, &
      -2.759285104469687e2_dp, 1.383577518672690e2_dp, &
      -3.066479806614716e1_dp, 2.506628277459239_dp ]
    real(dp), parameter :: b(5) = [ &
      -5.447609879822406e1_dp, 1.615858368580409e2_dp, &
      -1.556989798598866e2_dp, 6.680131188771972e1_dp, &
      -1.328068155288572e1_dp ]
    real(dp), parameter :: c(6) = [ &
      -7.784894002430293e-3_dp, -3.223964580411365e-1_dp, &
      -2.400758277161838_dp, -2.549732539343734_dp, &
      4.374664141464968_dp, 2.938163982698783_dp ]
    real(dp), parameter :: d(4) = [ &
      7.784695709041462e-3_dp, 3.224671290700398e-1_dp, &
      2.445134137142996_dp, 3.754408661907416_dp ]

    if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
    else if (p < 0.02425_dp) then
      q = sqrt(-2.0_dp * log(p))
      x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
        ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else if (p > 0.97575_dp) then
      q = sqrt(-2.0_dp * log(1.0_dp-p))
      x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
        ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else
      q = p - 0.5_dp
      r = q * q
      x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q / &
        (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
    end if
  end function normal_quantile

  pure elemental function clamp_probability(p) result(value)
    real(dp), intent(in) :: p
    real(dp) :: value
    value = min(max(p, tiny(1.0_dp)), 1.0_dp - epsilon(1.0_dp))
  end function clamp_probability

  pure elemental function log1mexp(logp) result(value)
    real(dp), intent(in) :: logp
    real(dp) :: value
    if (logp >= 0.0_dp) then
      value = -huge(1.0_dp)
    else if (logp < log(0.5_dp)) then
      value = log(1.0_dp - exp(logp))
    else
      value = log(-(exp(logp) - 1.0_dp))
    end if
  end function log1mexp

  pure recursive function regularized_gamma_p(a, x) result(value)
    real(dp), intent(in) :: a, x
    real(dp) :: value, sumv, del, ap
    integer :: n

    if (a <= 0.0_dp .or. x < 0.0_dp) then
      value = 0.0_dp
      return
    end if
    if (x <= tiny(1.0_dp)) then
      value = 0.0_dp
      return
    end if
    if (x >= a + 1.0_dp) then
      value = 1.0_dp - regularized_gamma_q(a, x)
      return
    end if

    ap = a
    sumv = 1.0_dp / a
    del = sumv
    do n = 1, 10000
      ap = ap + 1.0_dp
      del = del * x / ap
      sumv = sumv + del
      if (abs(del) <= abs(sumv) * 2.0e-15_dp) exit
    end do
    value = sumv * exp(-x + a * log(x) - log_gamma(a))
    value = min(max(value, 0.0_dp), 1.0_dp)
  end function regularized_gamma_p

  pure recursive function regularized_gamma_q(a, x) result(value)
    real(dp), intent(in) :: a, x
    real(dp) :: value, b, c, d, h, an, del
    integer :: i

    if (a <= 0.0_dp .or. x < 0.0_dp) then
      value = 1.0_dp
      return
    end if
    if (x <= tiny(1.0_dp)) then
      value = 1.0_dp
      return
    end if
    if (x < a + 1.0_dp) then
      value = 1.0_dp - regularized_gamma_p(a, x)
      return
    end if

    b = x + 1.0_dp - a
    c = 1.0_dp / tiny(1.0_dp)
    d = 1.0_dp / b
    h = d
    do i = 1, 10000
      an = -real(i, dp) * (real(i, dp) - a)
      b = b + 2.0_dp
      d = an * d + b
      if (abs(d) < tiny(1.0_dp)) d = tiny(1.0_dp)
      c = b + an / c
      if (abs(c) < tiny(1.0_dp)) c = tiny(1.0_dp)
      d = 1.0_dp / d
      del = d * c
      h = h * del
      if (abs(del - 1.0_dp) <= 2.0e-15_dp) exit
    end do
    value = exp(-x + a * log(x) - log_gamma(a)) * h
    value = min(max(value, 0.0_dp), 1.0_dp)
  end function regularized_gamma_q

  pure function regularized_beta(x, a, b) result(value)
    real(dp), intent(in) :: x, a, b
    real(dp) :: value, bt

    if (x <= 0.0_dp) then
      value = 0.0_dp
      return
    else if (x >= 1.0_dp) then
      value = 1.0_dp
      return
    else if (a <= 0.0_dp .or. b <= 0.0_dp) then
      value = 0.0_dp
      return
    end if

    bt = exp(log_gamma(a+b) - log_gamma(a) - log_gamma(b) + &
      a*log(x) + b*log(1.0_dp-x))
    if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
      value = bt * beta_continued_fraction(a, b, x) / a
    else
      value = 1.0_dp - bt * beta_continued_fraction(b, a, 1.0_dp-x) / b
    end if
    value = min(max(value, 0.0_dp), 1.0_dp)
  end function regularized_beta

  pure function beta_continued_fraction(a, b, x) result(value)
    real(dp), intent(in) :: a, b, x
    real(dp) :: value, qab, qap, qam, c, d, h, aa, del
    integer :: m, m2

    qab = a + b
    qap = a + 1.0_dp
    qam = a - 1.0_dp
    c = 1.0_dp
    d = 1.0_dp - qab*x/qap
    if (abs(d) < tiny(1.0_dp)) d = tiny(1.0_dp)
    d = 1.0_dp/d
    h = d
    do m = 1, 10000
      m2 = 2*m
      aa = real(m,dp)*(b-real(m,dp))*x / &
        ((qam+real(m2,dp))*(a+real(m2,dp)))
      d = 1.0_dp + aa*d
      if (abs(d) < tiny(1.0_dp)) d = tiny(1.0_dp)
      c = 1.0_dp + aa/c
      if (abs(c) < tiny(1.0_dp)) c = tiny(1.0_dp)
      d = 1.0_dp/d
      h = h*d*c
      aa = -(a+real(m,dp))*(qab+real(m,dp))*x / &
        ((a+real(m2,dp))*(qap+real(m2,dp)))
      d = 1.0_dp + aa*d
      if (abs(d) < tiny(1.0_dp)) d = tiny(1.0_dp)
      c = 1.0_dp + aa/c
      if (abs(c) < tiny(1.0_dp)) c = tiny(1.0_dp)
      d = 1.0_dp/d
      del = d*c
      h = h*del
      if (abs(del-1.0_dp) <= 2.0e-15_dp) exit
    end do
    value = h
  end function beta_continued_fraction

  function chi_square_survival(x, df) result(value)
    real(dp), intent(in) :: x
    integer, intent(in) :: df
    real(dp) :: value
    if (x <= 0.0_dp) then
      value = 1.0_dp
    else if (df <= 0) then
      value = -1.0_dp
    else
      value = regularized_gamma_q(0.5_dp*real(df,dp), 0.5_dp*x)
    end if
  end function chi_square_survival

  subroutine sort_real(x)
    real(dp), intent(inout) :: x(:)
    integer :: i, j
    real(dp) :: key
    do i = 2, size(x)
      key = x(i)
      j = i - 1
      do while (j >= 1)
        if (x(j) <= key) exit
        x(j+1) = x(j)
        j = j - 1
      end do
      x(j+1) = key
    end do
  end subroutine sort_real

  function type7_quantile(x, p) result(value)
    real(dp), intent(in) :: x(:), p
    real(dp) :: value, h, frac
    real(dp), allocatable :: work(:)
    integer :: lo

    if (size(x) == 0 .or. p < 0.0_dp .or. p > 1.0_dp) then
      value = huge(1.0_dp)
      return
    end if
    work = x
    call sort_real(work)
    if (size(work) == 1) then
      value = work(1)
      return
    end if
    h = 1.0_dp + real(size(work)-1,dp)*p
    lo = floor(h)
    frac = h - real(lo,dp)
    if (lo >= size(work)) then
      value = work(size(work))
    else
      value = (1.0_dp-frac)*work(lo) + frac*work(lo+1)
    end if
  end function type7_quantile

  pure function sample_mean(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value
    if (size(x) == 0) then
      value = 0.0_dp
    else
      value = sum(x)/real(size(x),dp)
    end if
  end function sample_mean

  pure function sample_variance(x, unbiased) result(value)
    real(dp), intent(in) :: x(:)
    logical, intent(in), optional :: unbiased
    real(dp) :: value, meanx, denom
    logical :: ub
    ub = .true.
    if (present(unbiased)) ub = unbiased
    if (size(x) < merge(2,1,ub)) then
      value = 0.0_dp
      return
    end if
    meanx = sum(x)/real(size(x),dp)
    denom = real(size(x)-merge(1,0,ub),dp)
    value = sum((x-meanx)**2)/denom
  end function sample_variance

  pure function weighted_mean(x, w) result(value)
    real(dp), intent(in) :: x(:), w(:)
    real(dp) :: value
    if (size(x) /= size(w) .or. sum(w) <= 0.0_dp) then
      value = huge(1.0_dp)
    else
      value = sum(w*x)/sum(w)
    end if
  end function weighted_mean

  pure function weighted_variance(x, w, unbiased) result(value)
    real(dp), intent(in) :: x(:), w(:)
    logical, intent(in), optional :: unbiased
    real(dp) :: value, meanx, sw, denom
    logical :: ub
    ub = .true.
    if (present(unbiased)) ub = unbiased
    if (size(x) /= size(w) .or. sum(w) <= 0.0_dp) then
      value = huge(1.0_dp)
      return
    end if
    sw = sum(w)
    meanx = sum(w*x)/sw
    denom = sw
    if (ub) denom = sw - 1.0_dp
    if (denom <= 0.0_dp) then
      value = 0.0_dp
    else
      value = sum(w*(x-meanx)**2)/denom
    end if
  end function weighted_variance

  function weighted_quantile(x, w, p) result(value)
    real(dp), intent(in) :: x(:), w(:), p
    real(dp) :: value, target, sw, cumulative, frac
    real(dp), allocatable :: xs(:), ws(:)
    integer :: i, j, n
    real(dp) :: tx, tw

    n = size(x)
    if (n == 0 .or. size(w) /= n .or. any(w < 0.0_dp) .or. &
        sum(w) <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
      value = huge(1.0_dp)
      return
    end if
    xs = x
    ws = w
    do i = 2, n
      tx = xs(i); tw = ws(i); j = i-1
      do while (j >= 1)
        if (xs(j) <= tx) exit
        xs(j+1)=xs(j); ws(j+1)=ws(j); j=j-1
      end do
      xs(j+1)=tx; ws(j+1)=tw
    end do
    sw = sum(ws)
    target = 1.0_dp + (sw-1.0_dp)*p
    cumulative = 0.0_dp
    do i = 1, n
      cumulative = cumulative + ws(i)
      if (cumulative >= target) then
        if (i == 1 .or. abs(cumulative-target) > 1.0e-12_dp) then
          value = xs(i)
        else if (i < n) then
          frac = target - floor(target)
          value = (1.0_dp-frac)*xs(i) + frac*xs(i+1)
        else
          value = xs(i)
        end if
        return
      end if
    end do
    value = xs(n)
  end function weighted_quantile

  subroutine invert_matrix(a, inverse, status)
    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable, intent(out) :: inverse(:, :)
    integer, intent(out) :: status
    real(dp), allocatable :: aug(:, :), row(:)
    real(dp) :: pivot, factor
    integer :: n, i, k, pivot_row

    n = size(a,1)
    if (n == 0 .or. size(a,2) /= n) then
      allocate(inverse(0,0)); status = fit_invalid_argument; return
    end if
    allocate(aug(n,2*n), row(2*n), inverse(n,n))
    aug = 0.0_dp
    aug(:,1:n) = a
    do i = 1, n
      aug(i,n+i) = 1.0_dp
    end do
    do k = 1, n
      pivot_row = k
      do i = k+1, n
        if (abs(aug(i,k)) > abs(aug(pivot_row,k))) pivot_row = i
      end do
      if (abs(aug(pivot_row,k)) <= 100.0_dp*epsilon(1.0_dp)) then
        inverse = 0.0_dp; status = fit_singular; return
      end if
      if (pivot_row /= k) then
        row = aug(k,:); aug(k,:) = aug(pivot_row,:); aug(pivot_row,:) = row
      end if
      pivot = aug(k,k)
      aug(k,:) = aug(k,:)/pivot
      do i = 1, n
        if (i == k) cycle
        factor = aug(i,k)
        aug(i,:) = aug(i,:) - factor*aug(k,:)
      end do
    end do
    inverse = aug(:,n+1:2*n)
    status = fit_success
  end subroutine invert_matrix

  subroutine numerical_hessian(fun, context, x, hessian)
    interface
      function fun(v, context) result(value)
        import dp
        real(dp), intent(in) :: v(:)
        class(*), intent(inout) :: context
        real(dp) :: value
      end function fun
    end interface
    class(*), intent(inout) :: context
    real(dp), intent(in) :: x(:)
    real(dp), allocatable, intent(out) :: hessian(:, :)
    real(dp), allocatable :: xp(:), xm(:), xpp(:), xpm(:), xmp(:), xmm(:)
    real(dp) :: hi, hj, f0
    integer :: i, j, n

    n = size(x)
    allocate(hessian(n,n), xp(n), xm(n), xpp(n), xpm(n), xmp(n), xmm(n))
    hessian = 0.0_dp
    f0 = fun(x,context)
    do i = 1, n
      hi = epsilon(1.0_dp)**0.25_dp * max(abs(x(i)), 1.0_dp)
      xp=x; xm=x; xp(i)=xp(i)+hi; xm(i)=xm(i)-hi
      hessian(i,i) = (fun(xp,context)-2.0_dp*f0+fun(xm,context))/(hi*hi)
      do j = i+1, n
        hj = epsilon(1.0_dp)**0.25_dp * max(abs(x(j)), 1.0_dp)
        xpp=x; xpm=x; xmp=x; xmm=x
        xpp(i)=xpp(i)+hi; xpp(j)=xpp(j)+hj
        xpm(i)=xpm(i)+hi; xpm(j)=xpm(j)-hj
        xmp(i)=xmp(i)-hi; xmp(j)=xmp(j)+hj
        xmm(i)=xmm(i)-hi; xmm(j)=xmm(j)-hj
        hessian(i,j)=(fun(xpp,context)-fun(xpm,context)-fun(xmp,context)+fun(xmm,context))/(4.0_dp*hi*hj)
        hessian(j,i)=hessian(i,j)
      end do
    end do
  end subroutine numerical_hessian

  subroutine covariance_from_hessian(hessian, covariance, status)
    real(dp), intent(in) :: hessian(:, :)
    real(dp), allocatable, intent(out) :: covariance(:, :)
    integer, intent(out) :: status
    call invert_matrix(hessian, covariance, status)
  end subroutine covariance_from_hessian

  subroutine seed_rng(seed)
    integer, intent(in), optional :: seed
    integer, allocatable :: values(:)
    integer :: n, i, s
    call random_seed(size=n)
    allocate(values(n))
    s = 104729
    if (present(seed)) s = seed
    do i = 1, n
      values(i) = modulo(s + 104729*i + 37*i*i, huge(1)-1)
      if (values(i) <= 0) values(i) = i + 1
    end do
    call random_seed(put=values)
  end subroutine seed_rng

  function random_normal() result(value)
    real(dp) :: value, u1, u2
    call random_number(u1)
    call random_number(u2)
    u1 = max(u1, tiny(1.0_dp))
    value = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi_dp*u2)
  end function random_normal

  recursive function random_gamma(shape, scale) result(value)
    real(dp), intent(in) :: shape, scale
    real(dp) :: value, d, c, x, v, u
    if (shape <= 0.0_dp .or. scale < 0.0_dp) then
      value = 0.0_dp
      return
    end if
    if (shape < 1.0_dp) then
      call random_number(u)
      value = random_gamma(shape+1.0_dp,scale)*u**(1.0_dp/shape)
      return
    end if
    d = shape - 1.0_dp/3.0_dp
    c = 1.0_dp/sqrt(9.0_dp*d)
    do
      x = random_normal()
      v = (1.0_dp+c*x)**3
      if (v <= 0.0_dp) cycle
      call random_number(u)
      if (u < 1.0_dp-0.0331_dp*x**4) exit
      if (log(u) < 0.5_dp*x*x+d*(1.0_dp-v+log(v))) exit
    end do
    value = scale*d*v
  end function random_gamma

  function random_poisson(lambda) result(value)
    real(dp), intent(in) :: lambda
    integer :: value, k
    real(dp) :: threshold, product, u, proposal
    if (lambda <= 0.0_dp) then
      value = 0
    else if (lambda < 30.0_dp) then
      threshold = exp(-lambda)
      product = 1.0_dp
      k = 0
      do
        k = k + 1
        call random_number(u)
        product = product*u
        if (product <= threshold) exit
      end do
      value = k-1
    else
      do
        proposal = lambda + sqrt(lambda)*random_normal()
        value = nint(proposal)
        if (value >= 0) exit
      end do
    end if
  end function random_poisson

end module fitdistrplus_math
