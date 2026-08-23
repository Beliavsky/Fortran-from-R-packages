module joker_asymptotics
  use joker_special, only: dp, trigamma_j, digamma_j
  implicit none
  private

  public :: finf_bern, finf_beta, finf_binom, finf_cat, finf_cauchy
  public :: finf_chisq, finf_dir, finf_exp, finf_gamma, finf_geom
  public :: finf_laplace, finf_lnorm, finf_multigam, finf_multinom
  public :: finf_nbinom, finf_norm, finf_pois

  public :: avar_bern_mle, avar_bern_me
  public :: avar_beta_mle, avar_beta_me, avar_beta_same
  public :: avar_binom_mle, avar_binom_me
  public :: avar_cat_mle, avar_cat_me
  public :: avar_cauchy_mle
  public :: avar_chisq_mle, avar_chisq_me
  public :: avar_dir_me, avar_dir_same
  public :: avar_exp_mle, avar_exp_me
  public :: avar_gamma_mle, avar_gamma_me, avar_gamma_same
  public :: avar_geom_mle, avar_geom_me
  public :: avar_laplace_mle, avar_laplace_me
  public :: avar_lnorm_mle, avar_lnorm_me
  public :: avar_multigam_me, avar_multigam_same
  public :: avar_multinom_mle, avar_multinom_me
  public :: avar_nbinom_mle, avar_nbinom_me
  public :: avar_norm_me
  public :: avar_pois_mle, avar_pois_me

contains

  pure real(dp) function finf_bern(prob) result(v)
    real(dp), intent(in) :: prob
    v = 1.0_dp / (prob * (1.0_dp - prob))
  end function finf_bern

  pure function finf_beta(shape1, shape2) result(d)
    real(dp), intent(in) :: shape1, shape2
    real(dp) :: d(2,2), t
    t = trigamma_j(shape1 + shape2)
    d(1,1) = trigamma_j(shape1) - t
    d(1,2) = -t
    d(2,1) = -t
    d(2,2) = trigamma_j(shape2) - t
  end function finf_beta

  pure real(dp) function finf_binom(sizep, prob) result(v)
    integer, intent(in) :: sizep
    real(dp), intent(in) :: prob
    v = real(sizep, dp) / (prob * (1.0_dp - prob))
  end function finf_binom

  pure function finf_cat(prob) result(d)
    real(dp), intent(in) :: prob(:)
    real(dp), allocatable :: d(:,:)
    integer :: k, i, j
    k = size(prob)
    allocate(d(k-1,k-1))
    do i = 1, k-1
      do j = 1, k-1
        d(i,j) = 1.0_dp / prob(k)
      end do
      d(i,i) = d(i,i) + 1.0_dp / prob(i)
    end do
  end function finf_cat

  pure function finf_cauchy(scale) result(d)
    real(dp), intent(in) :: scale
    real(dp) :: d(2,2)
    d = 0.0_dp
    d(1,1) = 1.0_dp / (2.0_dp * scale**2)
    d(2,2) = d(1,1)
  end function finf_cauchy

  pure real(dp) function finf_chisq(df) result(v)
    real(dp), intent(in) :: df
    v = 0.25_dp * trigamma_j(df / 2.0_dp)
  end function finf_chisq

  pure function finf_dir(alpha) result(d)
    real(dp), intent(in) :: alpha(:)
    real(dp), allocatable :: d(:,:)
    real(dp) :: t0
    integer :: k, i, j
    k = size(alpha)
    allocate(d(k,k))
    t0 = trigamma_j(sum(alpha))
    do i = 1, k
      do j = 1, k
        d(i,j) = -t0
      end do
      d(i,i) = d(i,i) + trigamma_j(alpha(i))
    end do
  end function finf_dir

  pure real(dp) function finf_exp(rate) result(v)
    real(dp), intent(in) :: rate
    v = 1.0_dp / rate**2
  end function finf_exp

  pure function finf_gamma(shape, scale) result(d)
    real(dp), intent(in) :: shape, scale
    real(dp) :: d(2,2)
    d(1,1) = trigamma_j(shape)
    d(1,2) = 1.0_dp / scale
    d(2,1) = d(1,2)
    d(2,2) = shape / scale**2
  end function finf_gamma

  pure real(dp) function finf_geom(prob) result(v)
    real(dp), intent(in) :: prob
    v = 1.0_dp / (prob**2 * (1.0_dp - prob))
  end function finf_geom

  pure function finf_laplace(scale) result(d)
    real(dp), intent(in) :: scale
    real(dp) :: d(2,2)
    ! Correct Fisher information for Laplace(mu, scale).
    d = 0.0_dp
    d(1,1) = 1.0_dp / scale**2
    d(2,2) = 1.0_dp / scale**2
  end function finf_laplace

  pure function finf_lnorm(sdlog) result(d)
    real(dp), intent(in) :: sdlog
    real(dp) :: d(2,2)
    ! Same information as Normal(meanlog, sdlog) after log transformation.
    d = 0.0_dp
    d(1,1) = 1.0_dp / sdlog**2
    d(2,2) = 2.0_dp / sdlog**2
  end function finf_lnorm

  pure function finf_multigam(shape, scale) result(d)
    real(dp), intent(in) :: shape(:), scale
    real(dp), allocatable :: d(:,:)
    integer :: k, i
    k = size(shape)
    allocate(d(k+1,k+1))
    d = 0.0_dp
    do i = 1, k
      d(i,i) = trigamma_j(shape(i))
      d(i,k+1) = 1.0_dp / scale
      d(k+1,i) = d(i,k+1)
    end do
    d(k+1,k+1) = sum(shape) / scale**2
  end function finf_multigam

  pure function finf_multinom(sizep, prob) result(d)
    integer, intent(in) :: sizep
    real(dp), intent(in) :: prob(:)
    real(dp), allocatable :: d(:,:)
    integer :: k, i, j
    ! Correct constrained multinomial Fisher information.  The upstream
    ! R source has a minus sign before the 1/p_k term; that is a typo.
    k = size(prob)
    allocate(d(k-1,k-1))
    do i = 1, k-1
      do j = 1, k-1
        d(i,j) = real(sizep,dp) / prob(k)
      end do
      d(i,i) = d(i,i) + real(sizep,dp) / prob(i)
    end do
  end function finf_multinom

  pure real(dp) function finf_nbinom(sizep, prob) result(v)
    real(dp), intent(in) :: sizep, prob
    v = sizep / (prob**2 * (1.0_dp - prob))
  end function finf_nbinom

  pure function finf_norm(sd) result(d)
    real(dp), intent(in) :: sd
    real(dp) :: d(2,2)
    d = 0.0_dp
    d(1,1) = 1.0_dp / sd**2
    d(2,2) = 2.0_dp / sd**2
  end function finf_norm

  pure real(dp) function finf_pois(lambda) result(v)
    real(dp), intent(in) :: lambda
    v = 1.0_dp / lambda
  end function finf_pois

  pure real(dp) function avar_bern_mle(prob) result(v)
    real(dp), intent(in) :: prob
    v = prob * (1.0_dp - prob)
  end function avar_bern_mle

  pure real(dp) function avar_bern_me(prob) result(v)
    real(dp), intent(in) :: prob
    v = avar_bern_mle(prob)
  end function avar_bern_me

  pure function avar_beta_mle(shape1, shape2) result(d)
    real(dp), intent(in) :: shape1, shape2
    real(dp) :: d(2,2), f(2,2), det
    f = finf_beta(shape1, shape2)
    det = f(1,1) * f(2,2) - f(1,2) * f(2,1)
    d(1,1) = f(2,2) / det
    d(1,2) = -f(1,2) / det
    d(2,1) = -f(2,1) / det
    d(2,2) = f(1,1) / det
  end function avar_beta_mle

  pure function avar_beta_me(shape1, shape2) result(d)
    real(dp), intent(in) :: shape1, shape2
    real(dp) :: d(2,2), prd, th, th2, s2, s4, m3, m4, den, e
    real(dp) :: s11, s22, s12
    prd = shape1 * shape2
    th = shape1 + shape2
    th2 = th**2
    s2 = prd / (th2 * (th + 1.0_dp))
    s4 = s2**2
    m3 = 2.0_dp * (shape2 - shape1) * s2 / (th * (th + 2.0_dp))
    m4 = 3.0_dp * prd * (prd * (th + 2.0_dp) + 2.0_dp * (shape2-shape1)**2) / &
      (th**4 * (th + 1.0_dp) * (th + 2.0_dp) * (th + 3.0_dp))
    den = (th + 1.0_dp)**2 * (th + 2.0_dp)**2 * s2
    e = (th + 1.0_dp)**3 * (m4 - s4 - m3**2 / s2) / s2
    s11 = (shape1 * (shape1 + 1.0_dp))**2 / den + shape1 * e / shape2
    s22 = (shape2 * (shape2 + 1.0_dp))**2 / den + shape2 * e / shape1
    s12 = -shape1 * (shape1 + 1.0_dp) * shape2 * (shape2 + 1.0_dp) / den + e
    d = reshape([s11, s12, s12, s22], [2,2])
  end function avar_beta_me

  pure function avar_beta_same(shape1, shape2) result(d)
    real(dp), intent(in) :: shape1, shape2
    real(dp) :: d(2,2), prd, th, th2, s2, cfac
    real(dp) :: m1(2,2), m2(2,2)
    prd = shape1 * shape2
    th = shape1 + shape2
    th2 = th**2
    s2 = prd / (th2 * (th + 1.0_dp))
    m1 = reshape([shape1**2, prd, prd, shape2**2], [2,2])
    m2 = reshape([prd, th2-prd, th2-prd, prd], [2,2])
    cfac = s2 * th2 * (trigamma_j(shape1) + trigamma_j(shape2)) + 1.0_dp
    d = cfac * m1 - m2 / (th + 1.0_dp)
  end function avar_beta_same

  pure real(dp) function avar_binom_mle(sizep, prob) result(v)
    integer, intent(in) :: sizep
    real(dp), intent(in) :: prob
    v = prob * (1.0_dp - prob) / real(sizep,dp)
  end function avar_binom_mle

  pure real(dp) function avar_binom_me(sizep, prob) result(v)
    integer, intent(in) :: sizep
    real(dp), intent(in) :: prob
    v = avar_binom_mle(sizep, prob)
  end function avar_binom_me

  pure function avar_cat_mle(prob) result(d)
    real(dp), intent(in) :: prob(:)
    real(dp), allocatable :: d(:,:)
    integer :: k, i, j
    k = size(prob)
    allocate(d(k-1,k-1))
    do i = 1, k-1
      do j = 1, k-1
        d(i,j) = -prob(i) * prob(j)
      end do
      d(i,i) = d(i,i) + prob(i)
    end do
  end function avar_cat_mle

  pure function avar_cat_me(prob) result(d)
    real(dp), intent(in) :: prob(:)
    real(dp), allocatable :: d(:,:)
    d = avar_cat_mle(prob)
  end function avar_cat_me

  pure function avar_cauchy_mle(scale) result(d)
    real(dp), intent(in) :: scale
    real(dp) :: d(2,2)
    d = 0.0_dp
    d(1,1) = 2.0_dp * scale**2
    d(2,2) = d(1,1)
  end function avar_cauchy_mle

  pure real(dp) function avar_chisq_mle(df) result(v)
    real(dp), intent(in) :: df
    v = 4.0_dp / trigamma_j(df / 2.0_dp)
  end function avar_chisq_mle

  pure real(dp) function avar_chisq_me(df) result(v)
    real(dp), intent(in) :: df
    v = 2.0_dp * df
  end function avar_chisq_me

  pure function avar_dir_me(alpha) result(d)
    real(dp), intent(in) :: alpha(:)
    real(dp), allocatable :: d(:,:)
    real(dp) :: a0, dn, a1, a2, a3, s2, s3, c0
    integer :: k, i, j
    k = size(alpha)
    allocate(d(k,k))
    a0 = sum(alpha)
    dn = a0**2 - sum(alpha**2)
    a1 = a0 + 1.0_dp
    a2 = a0 + 2.0_dp
    a3 = a0 + 3.0_dp
    s2 = sum(alpha**2)
    s3 = sum(alpha**3)
    c0 = (-4.0_dp*a0*(a0-1.0_dp)*a1**2*s3 + &
      (2.0_dp*a0**3+a0**2+a0)*s2**2 + &
      (2.0_dp*a0**5+2.0_dp*a0**4-6.0_dp*a0**3-4.0_dp*a0**2-2.0_dp*a0)*s2 + &
      a0**6+a0**5+2.0_dp*a0**3) / (dn**2*a1*a2*a3)
    do i = 1, k
      do j = 1, k
        d(i,j) = 2.0_dp*a0/(dn*a2) * &
          (alpha(i)*alpha(j)**2 + alpha(i)**2*alpha(j)) + &
          c0*alpha(i)*alpha(j)
      end do
      d(i,i) = d(i,i) + (a0/a1)*alpha(i)
    end do
  end function avar_dir_me

  pure function avar_dir_same(alpha) result(d)
    real(dp), intent(in) :: alpha(:)
    real(dp), allocatable :: d(:,:)
    real(dp) :: a0, par1, par2, c, fac
    integer :: k, i, j
    k = size(alpha)
    allocate(d(k,k))
    a0 = sum(alpha)
    par1 = 1.0_dp / (real(k-1,dp) * (a0 + 1.0_dp))
    par2 = 1.0_dp / (real(k-1,dp)**2 * (a0 + 1.0_dp))
    c = -par2*sum(alpha**2*[(trigamma_j(alpha(i)), i=1,k)]) + &
      a0*par2*sum(alpha*[(trigamma_j(alpha(i)), i=1,k)]) + &
      (a0+2.0_dp)*par1
    do i = 1, k
      do j = 1, k
        fac = c - (1.0_dp/alpha(i) + 1.0_dp/alpha(j))*a0*par1
        if (i == j) fac = fac + a0/(alpha(i)*(a0+1.0_dp))
        d(i,j) = alpha(i)*alpha(j)*fac
      end do
    end do
  end function avar_dir_same

  pure real(dp) function avar_exp_mle(rate) result(v)
    real(dp), intent(in) :: rate
    v = rate**2
  end function avar_exp_mle

  pure real(dp) function avar_exp_me(rate) result(v)
    real(dp), intent(in) :: rate
    v = avar_exp_mle(rate)
  end function avar_exp_me

  pure function avar_gamma_mle(shape, scale) result(d)
    real(dp), intent(in) :: shape, scale
    real(dp) :: d(2,2), c
    c = shape * trigamma_j(shape) - 1.0_dp
    d(1,1) = shape / c
    d(1,2) = -scale / c
    d(2,1) = d(1,2)
    d(2,2) = scale**2 * trigamma_j(shape) / c
  end function avar_gamma_mle

  pure function avar_gamma_me(shape, scale) result(d)
    real(dp), intent(in) :: shape, scale
    real(dp) :: d(2,2), s11, s22, s12
    s11 = 2.0_dp * shape * (shape + 1.0_dp)
    s22 = scale**2 * (2.0_dp*shape + 3.0_dp) / shape
    s12 = -2.0_dp * scale * (shape + 1.0_dp)
    d = reshape([s11, s12, s12, s22], [2,2])
  end function avar_gamma_me

  pure function avar_gamma_same(shape, scale) result(d)
    real(dp), intent(in) :: shape, scale
    real(dp) :: d(2,2), c1, c2
    c1 = 1.0_dp + shape * trigamma_j(shape + 1.0_dp)
    c2 = 1.0_dp + shape * trigamma_j(shape)
    d(1,1) = shape**2 * c1
    d(1,2) = -shape * scale * c1
    d(2,1) = d(1,2)
    d(2,2) = scale**2 * c2
  end function avar_gamma_same

  pure real(dp) function avar_geom_mle(prob) result(v)
    real(dp), intent(in) :: prob
    v = prob**2 * (1.0_dp - prob)
  end function avar_geom_mle

  pure real(dp) function avar_geom_me(prob) result(v)
    real(dp), intent(in) :: prob
    v = avar_geom_mle(prob)
  end function avar_geom_me

  pure function avar_laplace_mle(scale) result(d)
    real(dp), intent(in) :: scale
    real(dp) :: d(2,2)
    d = 0.0_dp
    d(1,1) = scale**2
    d(2,2) = scale**2
  end function avar_laplace_mle

  pure function avar_laplace_me(scale) result(d)
    real(dp), intent(in) :: scale
    real(dp) :: d(2,2)
    d = avar_laplace_mle(scale)
  end function avar_laplace_me

  pure function avar_lnorm_mle(sdlog) result(d)
    real(dp), intent(in) :: sdlog
    real(dp) :: d(2,2)
    d = 0.0_dp
    d(1,1) = sdlog**2
    d(2,2) = sdlog**2 / 2.0_dp
  end function avar_lnorm_mle

  pure function avar_lnorm_me(sdlog) result(d)
    real(dp), intent(in) :: sdlog
    real(dp) :: d(2,2)
    d = avar_lnorm_mle(sdlog)
  end function avar_lnorm_me

  pure function avar_multigam_me(shape, scale) result(d)
    real(dp), intent(in) :: shape(:), scale
    real(dp), allocatable :: d(:,:), a(:,:), b(:,:)
    integer :: k, i, j
    k = size(shape)
    allocate(a(k+1,2*k), b(2*k,2*k), d(k+1,k+1))
    a = 0.0_dp
    b = 0.0_dp
    do i = 1, k
      do j = 1, k
        a(i,j) = shape(i) * (2.0_dp + 1.0_dp/shape(j)) / &
          (real(k,dp)*scale)
        if (i == j) a(i,j) = a(i,j) + 1.0_dp/scale
        a(i,k+j) = -shape(i) / &
          (shape(j)*real(k,dp)*scale**2)
      end do
      a(k+1,i) = -(2.0_dp + 1.0_dp/shape(i)) / real(k,dp)
      a(k+1,k+i) = 1.0_dp / (shape(i)*real(k,dp)*scale)

      b(i,i) = shape(i)*scale**2
      b(k+i,k+i) = 2.0_dp*shape(i)*(shape(i)+1.0_dp)*scale**4 * &
        (2.0_dp*shape(i)+3.0_dp)
      b(i,k+i) = 2.0_dp*shape(i)*(shape(i)+1.0_dp)*scale**3
      b(k+i,i) = b(i,k+i)
    end do
    d = matmul(a, matmul(b, transpose(a)))
    d = 0.5_dp * (d + transpose(d))
  end function avar_multigam_me

  pure function avar_multigam_same(shape, scale) result(d)
    real(dp), intent(in) :: shape(:), scale
    real(dp), allocatable :: d(:,:), a(:,:), b(:,:)
    real(dp) :: a21, b11, b22, b33, b12, b13, b23
    integer :: k, i, j
    k = size(shape)
    allocate(a(k+1,3*k), b(3*k,3*k), d(k+1,k+1))
    a = 0.0_dp
    b = 0.0_dp
    do i = 1, k
      do j = 1, k
        a21 = -(digamma_j(shape(j)) + log(scale)) / real(k,dp)
        a(i,j) = -shape(i)*a21/scale
        if (i == j) a(i,j) = a(i,j) + 1.0_dp/scale
        a(i,k+j) = shape(i)*shape(j)/real(k,dp)
        a(i,2*k+j) = -shape(i)/(real(k,dp)*scale)
      end do
      a(k+1,i) = -(digamma_j(shape(i)) + log(scale))/real(k,dp)
      a(k+1,k+i) = -shape(i)*scale/real(k,dp)
      a(k+1,2*k+i) = 1.0_dp/real(k,dp)

      b11 = shape(i)*scale**2
      b22 = trigamma_j(shape(i))
      b33 = shape(i)*(shape(i)+1.0_dp)*scale**2 * &
        (trigamma_j(shape(i)+2.0_dp) + &
        (digamma_j(shape(i)+2.0_dp)+log(scale))**2) - &
        (shape(i)*scale)**2 * &
        (digamma_j(shape(i)+1.0_dp)+log(scale))**2
      b12 = scale
      b13 = shape(i)*(shape(i)+1.0_dp)*scale**2 * &
        (digamma_j(shape(i)+2.0_dp)+log(scale)) - &
        (shape(i)*scale)**2 * &
        (digamma_j(shape(i)+1.0_dp)+log(scale))
      b23 = shape(i)*scale * &
        (trigamma_j(shape(i)+1.0_dp) + &
        (digamma_j(shape(i)+1.0_dp)+log(scale))**2) - &
        shape(i)*scale*(digamma_j(shape(i))+log(scale)) * &
        (digamma_j(shape(i)+1.0_dp)+log(scale))

      b(i,i) = b11
      b(k+i,k+i) = b22
      b(2*k+i,2*k+i) = b33
      b(i,k+i) = b12
      b(k+i,i) = b12
      b(i,2*k+i) = b13
      b(2*k+i,i) = b13
      b(k+i,2*k+i) = b23
      b(2*k+i,k+i) = b23
    end do
    d = matmul(a, matmul(b, transpose(a)))
    d = 0.5_dp * (d + transpose(d))
  end function avar_multigam_same

  pure function avar_multinom_mle(sizep, prob) result(d)
    integer, intent(in) :: sizep
    real(dp), intent(in) :: prob(:)
    real(dp), allocatable :: d(:,:)
    d = avar_cat_mle(prob) / real(sizep,dp)
  end function avar_multinom_mle

  pure function avar_multinom_me(sizep, prob) result(d)
    integer, intent(in) :: sizep
    real(dp), intent(in) :: prob(:)
    real(dp), allocatable :: d(:,:)
    d = avar_multinom_mle(sizep, prob)
  end function avar_multinom_me

  pure real(dp) function avar_nbinom_mle(sizep, prob) result(v)
    real(dp), intent(in) :: sizep, prob
    v = prob**2 * (1.0_dp - prob) / sizep
  end function avar_nbinom_mle

  pure real(dp) function avar_nbinom_me(sizep, prob) result(v)
    real(dp), intent(in) :: sizep, prob
    v = avar_nbinom_mle(sizep, prob)
  end function avar_nbinom_me

  pure function avar_norm_me(sd) result(d)
    real(dp), intent(in) :: sd
    real(dp) :: d(2,2)
    d = 0.0_dp
    d(1,1) = sd**2
    d(2,2) = sd**2 / 2.0_dp
  end function avar_norm_me

  pure real(dp) function avar_pois_mle(lambda) result(v)
    real(dp), intent(in) :: lambda
    v = lambda
  end function avar_pois_mle

  pure real(dp) function avar_pois_me(lambda) result(v)
    real(dp), intent(in) :: lambda
    v = lambda
  end function avar_pois_me

end module joker_asymptotics
