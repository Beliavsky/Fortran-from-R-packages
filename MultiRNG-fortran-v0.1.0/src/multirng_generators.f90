! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module multirng_generators
  use multirng_kinds, only : dp
  use multirng_rng, only : rng_uniform, rng_normal, rng_gamma, rng_chisq, rng_poisson
  use multirng_rng, only : rng_binomial, rng_hypergeometric
  use multirng_linalg, only : cholesky_lower, invert_spd, is_symmetric_pd
  use multirng_math, only : normal_cdf
  implicit none
  private

  public :: draw_d_variate_normal
  public :: draw_d_variate_t
  public :: draw_d_variate_uniform
  public :: draw_dirichlet
  public :: draw_multinomial
  public :: draw_dirichlet_multinomial
  public :: draw_multivariate_hypergeometric
  public :: generate_point_in_sphere
  public :: draw_multivariate_laplace
  public :: draw_wishart, draw_wishart_flat
  public :: draw_inv_wishart, draw_inv_wishart_flat
  public :: draw_inv_wishart_legacy
  public :: draw_correlated_binary
  public :: loc_min

contains

  function draw_d_variate_normal(no_row, d, mean_vec, cov_mat) result(x)
    integer, intent(in) :: no_row, d
    real(dp), intent(in) :: mean_vec(:), cov_mat(:, :)
    real(dp), allocatable :: x(:, :)
    real(dp), allocatable :: l(:, :), z(:)
    integer :: i, j, info

    call validate_common(no_row, d)
    if (size(mean_vec) /= d) error stop "draw_d_variate_normal: wrong mean dimension"
    if (size(cov_mat, 1) /= d .or. size(cov_mat, 2) /= d) error stop "draw_d_variate_normal: wrong covariance dimension"
    if (.not. is_symmetric_pd(cov_mat)) error stop "draw_d_variate_normal: covariance must be symmetric positive definite"

    allocate(x(no_row, d), l(d, d), z(d))
    call cholesky_lower(cov_mat, l, info)
    if (info /= 0) error stop "draw_d_variate_normal: Cholesky failure"
    do i = 1, no_row
      do j = 1, d
        z(j) = rng_normal()
      end do
      x(i, :) = mean_vec + matmul(l, z)
    end do
  end function draw_d_variate_normal

  function draw_d_variate_t(dof, no_row, d, mean_vec, cov_mat) result(x)
    real(dp), intent(in) :: dof
    integer, intent(in) :: no_row, d
    real(dp), intent(in) :: mean_vec(:), cov_mat(:, :)
    real(dp), allocatable :: x(:, :)
    real(dp), allocatable :: l(:, :), z(:)
    real(dp) :: scale
    integer :: i, j, info

    call validate_common(no_row, d)
    if (dof <= 1.0_dp) error stop "draw_d_variate_t: degrees of freedom must exceed 1"
    if (size(mean_vec) /= d) error stop "draw_d_variate_t: wrong mean dimension"
    if (size(cov_mat, 1) /= d .or. size(cov_mat, 2) /= d) error stop "draw_d_variate_t: wrong covariance dimension"
    if (.not. is_symmetric_pd(cov_mat)) error stop "draw_d_variate_t: covariance must be symmetric positive definite"

    allocate(x(no_row, d), l(d, d), z(d))
    call cholesky_lower(cov_mat, l, info)
    if (info /= 0) error stop "draw_d_variate_t: Cholesky failure"

    ! Preserve MultiRNG 1.2.4 behavior: one chi-square draw is shared by
    ! every row in a call, rather than drawing one independent scale per row.
    scale = sqrt(dof / rng_chisq(dof))
    do i = 1, no_row
      do j = 1, d
        z(j) = rng_normal()
      end do
      x(i, :) = mean_vec + scale * matmul(l, z)
    end do
  end function draw_d_variate_t

  function draw_d_variate_uniform(no_row, d, cov_mat) result(x)
    integer, intent(in) :: no_row, d
    real(dp), intent(in) :: cov_mat(:, :)
    real(dp), allocatable :: x(:, :), z(:, :), zero(:)

    call validate_common(no_row, d)
    if (size(cov_mat, 1) /= d .or. size(cov_mat, 2) /= d) error stop "draw_d_variate_uniform: wrong covariance dimension"
    if (.not. is_symmetric_pd(cov_mat)) error stop "draw_d_variate_uniform: covariance must be symmetric positive definite"
    allocate(zero(d))
    zero = 0.0_dp
    z = draw_d_variate_normal(no_row, d, zero, cov_mat)
    allocate(x(no_row, d))
    x = normal_cdf(z)
  end function draw_d_variate_uniform

  function draw_dirichlet(no_row, d, alpha, beta) result(x)
    integer, intent(in) :: no_row, d
    real(dp), intent(in) :: alpha(:), beta
    real(dp), allocatable :: x(:, :)
    real(dp) :: s
    integer :: i, j

    call validate_common(no_row, d)
    if (size(alpha) /= d) error stop "draw_dirichlet: wrong alpha dimension"
    if (minval(alpha) <= 0.0_dp) error stop "draw_dirichlet: alpha must be positive"
    if (beta <= 0.0_dp) error stop "draw_dirichlet: beta must be positive"

    allocate(x(no_row, d))
    do i = 1, no_row
      do j = 1, d
        x(i, j) = rng_gamma(alpha(j), beta)
      end do
      s = sum(x(i, :))
      x(i, :) = x(i, :) / s
    end do
  end function draw_dirichlet

  function draw_multinomial(no_row, d, theta_in, n) result(x)
    integer, intent(in) :: no_row, d, n
    real(dp), intent(in) :: theta_in(:)
    integer, allocatable :: x(:, :)
    real(dp), allocatable :: theta(:)
    real(dp) :: denom, pj
    integer :: i, j, left

    if (no_row < 1) error stop "draw_multinomial: no_row must be positive"
    if (d < 1) error stop "draw_multinomial: d must be positive"
    if (size(theta_in) /= d) error stop "draw_multinomial: wrong theta dimension"
    if (minval(theta_in) < 0.0_dp) error stop "draw_multinomial: negative probability"
    if (abs(sum(theta_in) - 1.0_dp) > 100.0_dp * epsilon(1.0_dp)) error stop "draw_multinomial: probabilities must sum to one"
    if (n < 2) error stop "draw_multinomial: N must be at least 2"

    allocate(x(no_row, d), theta(d))
    x = 0
    do i = 1, no_row
      theta = theta_in
      left = n
      if (d == 1) then
        x(i, 1) = n
        cycle
      end if
      x(i, 1) = rng_binomial(left, theta(1))
      do j = 2, d - 1
        left = n - sum(x(i, 1:j - 1))
        denom = sum(theta(j:d))
        if (denom <= 0.0_dp) then
          pj = 0.0_dp
        else
          pj = theta(j) / denom
        end if
        x(i, j) = rng_binomial(left, pj)
      end do
      x(i, d) = n - sum(x(i, 1:d - 1))
    end do
  end function draw_multinomial

  function draw_dirichlet_multinomial(no_row, d, alpha, beta, n) result(x)
    integer, intent(in) :: no_row, d, n
    real(dp), intent(in) :: alpha(:), beta
    integer, allocatable :: x(:, :)
    real(dp), allocatable :: dir_draws(:, :), theta(:)
    integer :: j

    call validate_common(no_row, d)
    if (size(alpha) /= d) error stop "draw_dirichlet_multinomial: wrong alpha dimension"
    if (minval(alpha) <= 0.0_dp) error stop "draw_dirichlet_multinomial: alpha must be positive"
    if (beta <= 0.0_dp) error stop "draw_dirichlet_multinomial: beta must be positive"
    if (n < 2) error stop "draw_dirichlet_multinomial: N must be at least 2"

    ! Preserve the upstream implementation: draw no_row Dirichlet vectors,
    ! average them columnwise, then use that one averaged probability vector
    ! for all multinomial rows in this call.
    dir_draws = draw_dirichlet(no_row, d, alpha, beta)
    allocate(theta(d))
    do j = 1, d
      theta(j) = sum(dir_draws(:, j)) / real(no_row, dp)
    end do
    theta(d) = max(0.0_dp, 1.0_dp - sum(theta(1:d - 1)))
    theta = theta / sum(theta)
    x = draw_multinomial(no_row, d, theta, n)
  end function draw_dirichlet_multinomial

  function draw_multivariate_hypergeometric(no_row, d, mean_vec, k) result(x)
    integer, intent(in) :: no_row, d, mean_vec(:), k
    integer, allocatable :: x(:, :)
    integer :: i, j, remaining_draws, remaining_items

    call validate_common(no_row, d)
    if (size(mean_vec) /= d) error stop "draw_multivariate_hypergeometric: wrong item-count dimension"
    if (minval(mean_vec) <= 0) error stop "draw_multivariate_hypergeometric: counts must be positive"
    if (k <= 0) error stop "draw_multivariate_hypergeometric: k must be positive"
    if (k > sum(mean_vec)) error stop "draw_multivariate_hypergeometric: k exceeds population"

    allocate(x(no_row, d))
    x = 0
    do i = 1, no_row
      remaining_draws = k
      remaining_items = sum(mean_vec)
      do j = 1, d - 1
        x(i, j) = rng_hypergeometric(mean_vec(j), remaining_items - mean_vec(j), remaining_draws)
        remaining_draws = remaining_draws - x(i, j)
        remaining_items = remaining_items - mean_vec(j)
      end do
      x(i, d) = remaining_draws
    end do
  end function draw_multivariate_hypergeometric

  function generate_point_in_sphere(no_row, d) result(x)
    integer, intent(in) :: no_row, d
    real(dp), allocatable :: x(:, :)
    real(dp), allocatable :: u(:)
    real(dp) :: s, s1, s2, u1, u2, u3, u4
    integer :: i, j

    call validate_common(no_row, d)
    allocate(x(no_row, d), u(d))
    x = 0.0_dp

    if (d == 3) then
      do i = 1, no_row
        do
          u1 = 2.0_dp * rng_uniform() - 1.0_dp
          u2 = 2.0_dp * rng_uniform() - 1.0_dp
          s1 = u1 * u1 + u2 * u2
          if (s1 <= 1.0_dp) exit
        end do
        x(i, 1) = 2.0_dp * u1 * sqrt(max(0.0_dp, 1.0_dp - s1))
        x(i, 2) = 2.0_dp * u2 * sqrt(max(0.0_dp, 1.0_dp - s1))
        x(i, 3) = 1.0_dp - 2.0_dp * s1
      end do
    else if (d == 4) then
      do i = 1, no_row
        do
          u1 = 2.0_dp * rng_uniform() - 1.0_dp
          u2 = 2.0_dp * rng_uniform() - 1.0_dp
          u3 = 2.0_dp * rng_uniform() - 1.0_dp
          u4 = 2.0_dp * rng_uniform() - 1.0_dp
          s1 = u1 * u1 + u2 * u2
          s2 = u3 * u3 + u4 * u4
          if (s1 <= 1.0_dp .and. s2 <= 1.0_dp .and. s2 > tiny(1.0_dp)) exit
        end do
        x(i, 1) = u1
        x(i, 2) = u2
        x(i, 3) = u3 * sqrt((1.0_dp - s1) / s2)
        x(i, 4) = u4 * sqrt((1.0_dp - s1) / s2)
      end do
    else
      do i = 1, no_row
        do
          do j = 1, d
            u(j) = (2.0_dp * rng_uniform() - 1.0_dp)
          end do
          s = sum(u * u)
          if (s <= 1.0_dp .and. s > tiny(1.0_dp)) exit
        end do
        x(i, :) = u / sqrt(s)
      end do
    end if
  end function generate_point_in_sphere

  function draw_multivariate_laplace(no_row, d, gamma_shape, mu, sigma) result(x)
    integer, intent(in) :: no_row, d
    real(dp), intent(in) :: gamma_shape, mu(:), sigma(:, :)
    real(dp), allocatable :: x(:, :), s(:, :), l(:, :)
    real(dp) :: radius
    integer :: i, info

    if (no_row < 2) error stop "draw_multivariate_laplace: no_row must be at least 2"
    if (d < 2) error stop "draw_multivariate_laplace: d must be at least 2"
    if (gamma_shape <= 0.0_dp) error stop "draw_multivariate_laplace: gamma must be positive"
    if (size(mu) /= d) error stop "draw_multivariate_laplace: wrong mean dimension"
    if (size(sigma, 1) /= d .or. size(sigma, 2) /= d) error stop "draw_multivariate_laplace: wrong covariance dimension"
    if (.not. is_symmetric_pd(sigma)) error stop "draw_multivariate_laplace: Sigma must be symmetric positive definite"

    allocate(x(no_row, d), l(d, d))
    call cholesky_lower(sigma, l, info)
    if (info /= 0) error stop "draw_multivariate_laplace: Cholesky failure"
    s = generate_point_in_sphere(no_row, d)
    do i = 1, no_row
      radius = rng_gamma(real(d, dp), 1.0_dp) ** (1.0_dp / gamma_shape)
      x(i, :) = mu + radius * matmul(l, s(i, :))
    end do
  end function draw_multivariate_laplace

  function draw_wishart(no_row, d, nu, sigma) result(w)
    integer, intent(in) :: no_row, d, nu
    real(dp), intent(in) :: sigma(:, :)
    real(dp), allocatable :: w(:, :, :)
    real(dp), allocatable :: z(:, :), mean0(:)
    integer :: i

    call validate_common(no_row, d)
    if (nu < d) error stop "draw_wishart: nu must be at least d"
    if (size(sigma, 1) /= d .or. size(sigma, 2) /= d) error stop "draw_wishart: wrong scale dimension"
    if (.not. is_symmetric_pd(sigma)) error stop "draw_wishart: scale must be symmetric positive definite"

    allocate(w(no_row, d, d), mean0(d))
    mean0 = 0.0_dp
    do i = 1, no_row
      z = draw_d_variate_normal(nu, d, mean0, sigma)
      w(i, :, :) = matmul(transpose(z), z)
    end do
  end function draw_wishart

  function draw_wishart_flat(no_row, d, nu, sigma) result(wflat)
    integer, intent(in) :: no_row, d, nu
    real(dp), intent(in) :: sigma(:, :)
    real(dp), allocatable :: wflat(:, :), w(:, :, :)
    integer :: i, j, k, pos

    w = draw_wishart(no_row, d, nu, sigma)
    allocate(wflat(no_row, d * d))
    do i = 1, no_row
      pos = 0
      do j = 1, d
        do k = 1, d
          pos = pos + 1
          wflat(i, pos) = w(i, j, k)
        end do
      end do
    end do
  end function draw_wishart_flat

  function draw_inv_wishart(no_row, d, nu, inv_sigma) result(iw)
    integer, intent(in) :: no_row, d, nu
    real(dp), intent(in) :: inv_sigma(:, :)
    real(dp), allocatable :: iw(:, :, :), sigma(:, :), w(:, :, :)
    integer :: i, info

    call validate_common(no_row, d)
    if (nu < d) error stop "draw_inv_wishart: nu must be at least d"
    if (size(inv_sigma, 1) /= d .or. size(inv_sigma, 2) /= d) error stop "draw_inv_wishart: wrong inverse-scale dimension"
    if (.not. is_symmetric_pd(inv_sigma)) error stop "draw_inv_wishart: inverse scale must be symmetric positive definite"

    allocate(sigma(d, d), iw(no_row, d, d))
    sigma = inv_sigma
    w = draw_wishart(no_row, d, nu, sigma)
    do i = 1, no_row
      call invert_spd(w(i, :, :), iw(i, :, :), info)
      if (info /= 0) error stop "draw_inv_wishart: Wishart draw inversion failed"
    end do
  end function draw_inv_wishart

  function draw_inv_wishart_flat(no_row, d, nu, inv_sigma) result(wflat)
    integer, intent(in) :: no_row, d, nu
    real(dp), intent(in) :: inv_sigma(:, :)
    real(dp), allocatable :: wflat(:, :), w(:, :, :)
    integer :: i, j, k, pos

    w = draw_inv_wishart(no_row, d, nu, inv_sigma)
    allocate(wflat(no_row, d * d))
    do i = 1, no_row
      pos = 0
      do j = 1, d
        do k = 1, d
          pos = pos + 1
          wflat(i, pos) = w(i, j, k)
        end do
      end do
    end do
  end function draw_inv_wishart_flat

  function draw_inv_wishart_legacy(no_row, d, nu, inv_sigma) result(w)
    integer, intent(in) :: no_row, d, nu
    real(dp), intent(in) :: inv_sigma(:, :)
    real(dp), allocatable :: w(:, :, :), sigma(:, :)
    integer :: info

    ! Exact behavior of MultiRNG 1.2.4's R implementation. It inverts the
    ! supplied inverse-scale matrix, draws Wishart matrices, but does not invert
    ! the resulting draws. This is kept only for reproducibility of the upstream
    ! implementation; draw_inv_wishart() above implements the documented law.
    allocate(sigma(d, d))
    call invert_spd(inv_sigma, sigma, info)
    if (info /= 0) error stop "draw_inv_wishart_legacy: inverse scale inversion failed"
    w = draw_wishart(no_row, d, nu, sigma)
  end function draw_inv_wishart_legacy

  function loc_min(my_mat, d) result(where)
    real(dp), intent(in) :: my_mat(:, :)
    integer, intent(in) :: d
    integer :: where(2)
    real(dp) :: best
    integer :: i, j

    if (size(my_mat, 1) /= d .or. size(my_mat, 2) /= d) error stop "loc_min: matrix dimension mismatch"
    best = huge(1.0_dp)
    where = 0
    do i = 1, d
      do j = 1, d
        if (my_mat(i, j) > 0.0_dp .and. my_mat(i, j) < best) then
          best = my_mat(i, j)
          where = [i, j]
        end if
      end do
    end do
  end function loc_min

  function draw_correlated_binary(no_row, d, prop_vec, corr_mat) result(y)
    integer, intent(in) :: no_row, d
    real(dp), intent(in) :: prop_vec(:), corr_mat(:, :)
    integer, allocatable :: y(:, :)
    real(dp), allocatable :: alpha(:, :), lambda(:)
    logical, allocatable :: incidence(:, :), active(:)
    real(dp) :: lim, rate, tol
    integer :: i, j, k, ncomp, idx(2), pdraw

    call validate_common(no_row, d)
    if (size(prop_vec) /= d) error stop "draw_correlated_binary: wrong mean dimension"
    if (size(corr_mat, 1) /= d .or. size(corr_mat, 2) /= d) error stop "draw_correlated_binary: wrong correlation dimension"
    if (minval(prop_vec) <= 0.0_dp .or. maxval(prop_vec) >= 1.0_dp) error stop "draw_correlated_binary: means must lie in (0,1)"
    if (maxval(abs(corr_mat - transpose(corr_mat))) > 100.0_dp * epsilon(1.0_dp)) then
      error stop "draw_correlated_binary: correlation matrix must be symmetric"
    end if
    if (maxval(abs([(corr_mat(i, i) - 1.0_dp, i=1,d)])) > 100.0_dp * epsilon(1.0_dp)) then
      error stop "draw_correlated_binary: diagonal correlations must equal one"
    end if
    if (minval(corr_mat) < 0.0_dp .or. maxval(corr_mat) > 1.0_dp) then
      error stop "draw_correlated_binary: correlations must lie in [0,1]"
    end if

    do i = 1, d
      do j = 1, d
        lim = min(sqrt(prop_vec(j) * (1.0_dp - prop_vec(i)) / &
          (prop_vec(i) * (1.0_dp - prop_vec(j)))), &
          sqrt(prop_vec(i) * (1.0_dp - prop_vec(j)) / &
          (prop_vec(j) * (1.0_dp - prop_vec(i)))))
        if (corr_mat(i, j) > lim + 1000.0_dp * epsilon(1.0_dp)) then
          error stop "draw_correlated_binary: correlation exceeds Bernoulli upper limit"
        end if
      end do
    end do

    allocate(alpha(d, d), lambda(d * d), incidence(d, d * d), active(d), y(no_row, d))
    do i = 1, d
      do j = 1, d
        alpha(i, j) = log(1.0_dp + corr_mat(i, j) * &
          sqrt((1.0_dp - prop_vec(i)) * (1.0_dp - prop_vec(j)) / &
          (prop_vec(i) * prop_vec(j))))
      end do
    end do

    tol = 1000.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(alpha))
    incidence = .false.
    lambda = 0.0_dp
    ncomp = 0

    do while (sum(max(alpha, 0.0_dp)) > tol)
      idx = loc_min_positive(alpha, tol)
      if (idx(1) == 0) exit
      rate = alpha(idx(1), idx(2))
      active = .false.
      active(idx(1)) = .true.
      active(idx(2)) = .true.
      do j = 1, d
        if (minval(alpha(:, j)) > tol) active(j) = .true.
      end do

      ncomp = ncomp + 1
      if (ncomp > d * d) error stop "draw_correlated_binary: decomposition failed to terminate"
      lambda(ncomp) = rate
      incidence(:, ncomp) = active

      if (alpha(idx(1), idx(1)) <= tol .or. alpha(idx(2), idx(2)) <= tol) then
        error stop "draw_correlated_binary: unsupported parameter configuration"
      end if

      do i = 1, d
        if (.not. active(i)) cycle
        do j = 1, d
          if (active(j)) alpha(i, j) = alpha(i, j) - rate
        end do
      end do
      where (abs(alpha) <= tol) alpha = 0.0_dp
      if (minval(alpha) < -100.0_dp * tol) error stop "draw_correlated_binary: negative decomposition remainder"
      where (alpha < 0.0_dp) alpha = 0.0_dp
    end do

    y = 1
    do k = 1, no_row
      do j = 1, ncomp
        pdraw = rng_poisson(lambda(j))
        if (pdraw <= 0) cycle
        do i = 1, d
          if (incidence(i, j)) y(k, i) = 0
        end do
      end do
    end do
  end function draw_correlated_binary

  function loc_min_positive(a, tol) result(where)
    real(dp), intent(in) :: a(:, :), tol
    integer :: where(2)
    real(dp) :: best
    integer :: i, j

    best = huge(1.0_dp)
    where = 0
    do i = 1, size(a, 1)
      do j = 1, size(a, 2)
        if (a(i, j) > tol .and. a(i, j) < best) then
          best = a(i, j)
          where = [i, j]
        end if
      end do
    end do
  end function loc_min_positive

  subroutine validate_common(no_row, d)
    integer, intent(in) :: no_row, d
    if (no_row < 1) error stop "MultiRNG: no_row must be a positive integer"
    if (d < 2) error stop "MultiRNG: d must be at least 2"
  end subroutine validate_common

end module multirng_generators
