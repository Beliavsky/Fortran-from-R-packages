module rlrsim_api
  use rlrsim_kinds, only : dp
  use rlrsim_numerics, only : chisq_random, residualize_against, seed_random_number, &
    singular_values_squared, sort_descending
  implicit none
  private

  type, public :: rlrsim_result
    real(dp), allocatable :: statistic(:)
    real(dp), allocatable :: lambda(:)
  end type rlrsim_result

  type, public :: exact_test_result
    real(dp) :: statistic = 0.0_dp
    real(dp) :: p_value = 1.0_dp
    real(dp), allocatable :: sample(:)
    real(dp), allocatable :: lambda(:)
  end type exact_test_result

  public :: dp
  public :: lrt_sim
  public :: rlrt_sim
  public :: exact_lrt_from_design
  public :: exact_rlrt_from_design

contains

  subroutine lrt_sim(x, z, q, sqrt_sigma, result, seed, nsim, log_grid_hi, log_grid_lo, gridlength)
    real(dp), intent(in) :: x(:, :) !! Fixed-effect design matrix under the alternative, shape (n,p).
    real(dp), intent(in) :: z(:, :) !! Random-effect design matrix under the alternative, shape (n,k).
    integer, intent(in) :: q !! Number of fixed-effect restrictions tested jointly with the variance component; must be nonnegative.
    real(dp), intent(in) :: sqrt_sigma(:, :) !! Square root of random-effect correlation, shape (k,k).
    type(rlrsim_result), intent(out) :: result !! Simulated LRT values and maximizing lambda values for each simulation.
    integer, optional, intent(in) :: seed !! Optional deterministic seed for Fortran's intrinsic RNG; stream differs from R's RNG.
    integer, optional, intent(in) :: nsim !! Number of null simulations; defaults to 10000 and must be positive.
    real(dp), optional, intent(in) :: log_grid_hi !! Log-scale upper lambda-grid endpoint; default 8.
    real(dp), optional, intent(in) :: log_grid_lo !! Log-scale lower lambda-grid endpoint; default -10.
    integer, optional, intent(in) :: gridlength !! Total lambda-grid length including zero; defaults to 200 and must be at least 2.
    real(dp), allocatable :: mu(:)
    real(dp), allocatable :: xi(:)
    real(dp), allocatable :: zr(:, :)
    real(dp), allocatable :: a(:, :)
    real(dp), allocatable :: lambda_grid(:)
    real(dp) :: hi
    real(dp) :: lo
    real(dp) :: denom
    integer :: ns
    integer :: gl
    integer :: n
    integer :: p
    integer :: k

    n = size(x, 1)
    p = size(x, 2)
    k = size(z, 2)
    if (size(z, 1) /= n) error stop "lrt_sim: x and z must have the same row count"
    if (size(sqrt_sigma, 1) /= k .or. size(sqrt_sigma, 2) /= k) &
      error stop "lrt_sim: sqrt_sigma must be k by k"
    if (q < 0) error stop "lrt_sim: q must be nonnegative"
    if (n - p - k < 0) error stop "lrt_sim: n - p - k must be nonnegative"

    ns = 10000
    if (present(nsim)) ns = nsim
    gl = 200
    if (present(gridlength)) gl = gridlength
    lo = -10.0_dp
    if (present(log_grid_lo)) lo = log_grid_lo
    hi = 8.0_dp
    if (present(log_grid_hi)) hi = log_grid_hi
    if (ns <= 0) error stop "lrt_sim: nsim must be positive"
    if (gl < 2) error stop "lrt_sim: gridlength must be at least 2"
    if (present(seed)) call seed_random_number(seed)

    call residualize_against(x, z, zr)
    a = matmul(sqrt_sigma, transpose(zr))
    call singular_values_squared(a, mu)
    a = matmul(sqrt_sigma, transpose(z))
    call singular_values_squared(a, xi)
    if (size(mu) /= k .or. size(xi) /= k) &
      error stop "lrt_sim: design requires k <= n, matching the upstream exactLRT requirements"

    denom = max(maxval(mu), maxval(xi))
    if (denom <= 0.0_dp) error stop "lrt_sim: all singular values are zero"
    mu = mu / denom
    denom = max(maxval(mu), maxval(xi))
    if (denom <= 0.0_dp) error stop "lrt_sim: invalid singular-value normalization"
    xi = xi / denom

    call regular_lambda_grid(0.0_dp, gl, lo, hi, lambda_grid)
    call rlrsim_kernel(p, k, n, ns, q, mu, lambda_grid, 0.0_dp, xi, .false., result)
  end subroutine lrt_sim

  subroutine rlrt_sim(x, z, sqrt_sigma, result, lambda0, seed, nsim, use_approx, log_grid_hi, log_grid_lo, gridlength)
    real(dp), intent(in) :: x(:, :) !! Fixed-effect design matrix, shape (n,p), used to residualize z.
    real(dp), intent(in) :: z(:, :) !! Random-effect design matrix, shape (n,k0), for the tested variance component.
    real(dp), intent(in) :: sqrt_sigma(:, :) !! Correlation square root, shape (k0,k0); identity is allowed.
    type(rlrsim_result), intent(out) :: result !! Simulated RLRT values and maximizing lambda values for each simulation.
    real(dp), optional, intent(in) :: lambda0 !! Null variance ratio; defaults to zero and must be nonnegative.
    integer, optional, intent(in) :: seed !! Optional deterministic seed for Fortran's intrinsic RNG; stream differs from R's RNG.
    integer, optional, intent(in) :: nsim !! Number of null simulations; defaults to 10000 and must be positive.
    real(dp), optional, intent(in) :: use_approx !! R-style eigenvalue threshold; zero disables approximation.
    real(dp), optional, intent(in) :: log_grid_hi !! Upper lambda-grid endpoint on natural-log scale; defaults to 8.
    real(dp), optional, intent(in) :: log_grid_lo !! Lower lambda-grid endpoint on natural-log scale; defaults to -10.
    integer, optional, intent(in) :: gridlength !! Lambda-grid length; default 200, minimum 4 if lambda0 > 0.
    real(dp), allocatable :: mu(:)
    real(dp), allocatable :: zr(:, :)
    real(dp), allocatable :: a(:, :)
    real(dp), allocatable :: lambda_grid(:)
    real(dp) :: hi
    real(dp) :: lo
    real(dp) :: lam0
    real(dp) :: approx_ratio
    real(dp) :: total
    integer :: ns
    integer :: gl
    integer :: n
    integer :: p
    integer :: k0
    integer :: k
    integer :: new_k

    n = size(x, 1)
    p = size(x, 2)
    k0 = size(z, 2)
    if (size(z, 1) /= n) error stop "rlrt_sim: x and z must have the same row count"
    if (size(sqrt_sigma, 1) /= k0 .or. size(sqrt_sigma, 2) /= k0) &
      error stop "rlrt_sim: sqrt_sigma must match the number of columns of z"

    lam0 = 0.0_dp
    if (present(lambda0)) lam0 = lambda0
    if (lam0 < 0.0_dp) error stop "rlrt_sim: lambda0 must be nonnegative"
    ns = 10000
    if (present(nsim)) ns = nsim
    gl = 200
    if (present(gridlength)) gl = gridlength
    lo = -10.0_dp
    if (present(log_grid_lo)) lo = log_grid_lo
    hi = 8.0_dp
    if (present(log_grid_hi)) hi = log_grid_hi
    approx_ratio = 0.0_dp
    if (present(use_approx)) approx_ratio = use_approx
    if (ns <= 0) error stop "rlrt_sim: nsim must be positive"
    if (gl < 2) error stop "rlrt_sim: gridlength must be at least 2"
    if (lam0 > exp(hi)) hi = log(10.0_dp * lam0)
    if (lam0 > 0.0_dp .and. lam0 < exp(lo)) lo = log(lam0 / 10.0_dp)
    if (present(seed)) call seed_random_number(seed)

    call residualize_against(x, z, zr)
    a = matmul(sqrt_sigma, transpose(zr))
    call singular_values_squared(a, mu)
    k = min(n, k0)
    if (size(mu) /= k) error stop "rlrt_sim: unexpected singular-value count"
    if (maxval(mu) <= 0.0_dp) error stop "rlrt_sim: all residualized singular values are zero"
    mu = mu / maxval(mu)

    if (approx_ratio > 0.0_dp) then
      if (balanced_anova_pattern(mu) .and. n - p - k + 1 > 0) then
        call rlrt_balanced_anova_approx(n, p, k, mu(1), ns, result)
        return
      end if
      total = sum(mu)
      if (mu(1) / total > approx_ratio .and. n - p - 1 > 0) then
        call rlrt_dominating_eigenvalue_approx(n, p, mu(1), ns, result)
        return
      end if
      if (total > 0.0_dp) then
        new_k = count(cumulative_fraction(mu) < approx_ratio)
        new_k = max(new_k, 1)
        if (new_k < k) then
          mu = mu(1:new_k)
          k = new_k
        end if
      end if
    end if

    call regular_lambda_grid(lam0, gl, lo, hi, lambda_grid)
    call rlrsim_kernel(p, k, n, ns, 0, mu, lambda_grid, lam0, mu, .true., result)
  end subroutine rlrt_sim

  subroutine exact_lrt_from_design(x, z, q, sqrt_sigma, observed_lrt, result, seed, nsim, log_grid_hi, log_grid_lo, gridlength)
    real(dp), intent(in) :: x(:, :) !! Fixed-effect design matrix under the alternative, shape (n,p).
    real(dp), intent(in) :: z(:, :) !! Random-effect design matrix under the alternative, shape (n,k).
    integer, intent(in) :: q !! Number of fixed-effect restrictions imposed by the null model.
    real(dp), intent(in) :: sqrt_sigma(:, :) !! Square root of the random-effect correlation matrix, shape (k,k).
    real(dp), intent(in) :: observed_lrt !! Observed twice-log-likelihood difference, truncated to zero internally when negative.
    type(exact_test_result), intent(out) :: result !! Observed statistic, simulation p-value, null sample, and maximizing lambdas.
    integer, optional, intent(in) :: seed !! Optional deterministic seed for Fortran's intrinsic RNG.
    integer, optional, intent(in) :: nsim !! Number of null simulations; defaults to 10000.
    real(dp), optional, intent(in) :: log_grid_hi !! Upper lambda-grid endpoint on natural-log scale; defaults to 8.
    real(dp), optional, intent(in) :: log_grid_lo !! Lower lambda-grid endpoint on natural-log scale; defaults to -10.
    integer, optional, intent(in) :: gridlength !! Total grid length including zero; defaults to 200.
    type(rlrsim_result) :: sim

    call lrt_sim(x, z, q, sqrt_sigma, sim, seed, nsim, log_grid_hi, log_grid_lo, gridlength)
    result%statistic = max(0.0_dp, observed_lrt)
    result%sample = sim%statistic
    result%lambda = sim%lambda
    result%p_value = real(count(result%statistic < result%sample), dp) / real(size(result%sample), dp)
  end subroutine exact_lrt_from_design

  subroutine exact_rlrt_from_design(x, z, sqrt_sigma, observed_rlrt, result, lambda0, seed, &
    nsim, use_approx, log_grid_hi, log_grid_lo, gridlength)
    real(dp), intent(in) :: x(:, :) !! Fixed-effect design matrix, shape (n,p), for the tested mixed model.
    real(dp), intent(in) :: z(:, :) !! Random-effect design matrix, shape (n,k), for the tested variance component.
    real(dp), intent(in) :: sqrt_sigma(:, :) !! Square root of the random-effect correlation matrix, shape (k,k).
    real(dp), intent(in) :: observed_rlrt !! Observed restricted likelihood-ratio statistic; negative values are truncated to zero.
    type(exact_test_result), intent(out) :: result !! Observed statistic, simulation p-value, null sample, and maximizing lambdas.
    real(dp), optional, intent(in) :: lambda0 !! Null variance ratio; defaults to zero.
    integer, optional, intent(in) :: seed !! Optional deterministic seed for Fortran's intrinsic RNG.
    integer, optional, intent(in) :: nsim !! Number of null simulations; defaults to 10000 when the observed statistic is positive.
    real(dp), optional, intent(in) :: use_approx !! Eigenvalue truncation threshold corresponding to R's use.approx argument.
    real(dp), optional, intent(in) :: log_grid_hi !! Upper lambda-grid endpoint on natural-log scale; defaults to 8.
    real(dp), optional, intent(in) :: log_grid_lo !! Lower lambda-grid endpoint on natural-log scale; defaults to -10.
    integer, optional, intent(in) :: gridlength !! Requested lambda-grid length; defaults to 200.
    type(rlrsim_result) :: sim

    result%statistic = max(0.0_dp, observed_rlrt)
    if (result%statistic == 0.0_dp) then
      result%p_value = 1.0_dp
      allocate(result%sample(0), result%lambda(0))
      return
    end if
    call rlrt_sim(x, z, sqrt_sigma, sim, lambda0, seed, nsim, use_approx, log_grid_hi, log_grid_lo, gridlength)
    result%sample = sim%statistic
    result%lambda = sim%lambda
    result%p_value = real(count(result%statistic < result%sample), dp) / real(size(result%sample), dp)
  end subroutine exact_rlrt_from_design

  subroutine rlrt_balanced_anova_approx(n, p, k, mu1, nsim, result)
    integer, intent(in) :: n !! Number of observations in the mixed model.
    integer, intent(in) :: p !! Number of fixed-effect parameters removed by REML.
    integer, intent(in) :: k !! Number of normalized eigenvalues, including the negligible final eigenvalue.
    real(dp), intent(in) :: mu1 !! Common leading normalized eigenvalue for the balanced-ANOVA approximation.
    integer, intent(in) :: nsim !! Number of null simulations to generate.
    type(rlrsim_result), intent(out) :: result !! Simulated RLRT statistics and analytic maximizing variance ratios.
    real(dp) :: w_k
    real(dp) :: w_n
    real(dp) :: ratio
    integer :: df_k
    integer :: df_n
    integer :: i

    df_k = k - 1
    df_n = n - p - k + 1
    if (df_k <= 0 .or. df_n <= 0) error stop "rlrt_balanced_anova_approx: invalid degrees of freedom"
    allocate(result%statistic(nsim), result%lambda(nsim))
    do i = 1, nsim
      w_k = chisq_random(real(df_k, dp))
      w_n = chisq_random(real(df_n, dp))
      ratio = real(df_n, dp) * w_k / (real(df_k, dp) * w_n)
      result%lambda(i) = max(0.0_dp, (ratio - 1.0_dp) / mu1)
      if (result%lambda(i) == 0.0_dp) then
        result%statistic(i) = 0.0_dp
      else
        result%statistic(i) = real(n - p, dp) * log((w_k + w_n) / real(n - p, dp)) &
          - real(df_n, dp) * log(w_n / real(df_n, dp)) &
          - real(df_k, dp) * log(w_k / real(df_k, dp))
      end if
    end do
  end subroutine rlrt_balanced_anova_approx

  subroutine rlrt_dominating_eigenvalue_approx(n, p, mu1, nsim, result)
    integer, intent(in) :: n !! Number of observations in the mixed model.
    integer, intent(in) :: p !! Number of fixed-effect parameters removed by REML.
    real(dp), intent(in) :: mu1 !! Dominating normalized squared singular value.
    integer, intent(in) :: nsim !! Number of null simulations to generate.
    type(rlrsim_result), intent(out) :: result !! Simulated RLRT statistics and analytic maximizing variance ratios.
    real(dp) :: w1
    real(dp) :: w_n
    real(dp) :: ratio
    integer :: df_n
    integer :: i

    df_n = n - p - 1
    if (df_n <= 0) error stop "rlrt_dominating_eigenvalue_approx: invalid degrees of freedom"
    allocate(result%statistic(nsim), result%lambda(nsim))
    do i = 1, nsim
      w1 = chisq_random(1.0_dp)
      w_n = chisq_random(real(df_n, dp))
      ratio = real(df_n, dp) * w1 / w_n
      result%lambda(i) = max(0.0_dp, (ratio - 1.0_dp) / mu1)
      if (result%lambda(i) == 0.0_dp) then
        result%statistic(i) = 0.0_dp
      else
        result%statistic(i) = real(n - p, dp) * log((w1 + w_n) / real(n - p, dp)) &
          - log(w1) - real(df_n, dp) * log(w_n / real(df_n, dp))
      end if
    end do
  end subroutine rlrt_dominating_eigenvalue_approx

  pure logical function balanced_anova_pattern(mu) result(is_balanced)
    real(dp), intent(in) :: mu(:) !! Descending normalized squared singular values checked after six-decimal rounding.
    real(dp) :: first_value
    real(dp) :: second_value
    real(dp) :: rounded
    logical :: have_second
    integer :: i

    is_balanced = .false.
    if (size(mu) < 2) return
    first_value = anint(mu(1) * 1.0e6_dp) / 1.0e6_dp
    second_value = 0.0_dp
    have_second = .false.
    do i = 2, size(mu)
      rounded = anint(mu(i) * 1.0e6_dp) / 1.0e6_dp
      if (rounded == first_value) cycle
      if (.not. have_second) then
        second_value = rounded
        have_second = .true.
      else if (rounded /= second_value) then
        return
      end if
    end do
    if (.not. have_second) return
    is_balanced = 1000.0_dp * mu(size(mu)) < mu(1)
  end function balanced_anova_pattern

  subroutine rlrsim_kernel(p, k, n, nsim, q, mu, lambda, lambda0, xi, reml, result)
    integer, intent(in) :: p !! Number of fixed-effect parameters in the alternative model.
    integer, intent(in) :: k !! Number of retained random-effect singular values.
    integer, intent(in) :: n !! Number of observations.
    integer, intent(in) :: nsim !! Number of null simulations to generate.
    integer, intent(in) :: q !! Number of tested fixed-effect restrictions; used only for ordinary LRT simulation.
    real(dp), intent(in) :: mu(:) !! Normalized squared singular values governing the numerator terms, length k.
    real(dp), intent(in) :: lambda(:) !! Candidate variance-ratio grid, ordered increasingly and including zero.
    real(dp), intent(in) :: lambda0 !! Null variance ratio used in numerator and denominator factors.
    real(dp), intent(in) :: xi(:) !! Normalized squared singular values used in the log-determinant term, length k.
    logical, intent(in) :: reml !! True for restricted likelihood simulation and false for ordinary likelihood simulation.
    type(rlrsim_result), intent(out) :: result !! Simulated statistics and selected lambda grid values.
    real(dp), allocatable :: f_n(:, :)
    real(dp), allocatable :: f_d(:, :)
    real(dp), allocatable :: sumlog(:)
    real(dp), allocatable :: chi1(:)
    real(dp) :: chi_k
    real(dp) :: chi_sum
    real(dp) :: lr
    real(dp) :: num
    real(dp) :: den
    integer :: df_chi_k
    integer :: n0
    integer :: is
    integer :: ig
    integer :: ik
    integer :: best

    if (size(mu) /= k .or. size(xi) /= k) error stop "rlrsim_kernel: eigenvalue length mismatch"
    allocate(f_n(size(lambda), k), f_d(size(lambda), k), sumlog(size(lambda)))
    do ig = 1, size(lambda)
      sumlog(ig) = 0.0_dp
      do ik = 1, k
        f_n(ig, ik) = ((lambda(ig) - lambda0) * mu(ik)) / (1.0_dp + lambda(ig) * mu(ik))
        f_d(ig, ik) = (1.0_dp + lambda0 * mu(ik)) / (1.0_dp + lambda(ig) * mu(ik))
        sumlog(ig) = sumlog(ig) + log1p_safe(lambda(ig) * xi(ik))
      end do
    end do

    df_chi_k = max(n - p - k, 0)
    if (reml) then
      n0 = n - p
    else
      n0 = n
    end if
    allocate(result%statistic(nsim), result%lambda(nsim), chi1(k))
    result%statistic = 0.0_dp
    result%lambda = lambda(1)

    do is = 1, nsim
      chi_k = chisq_random(real(df_chi_k, dp))
      do ik = 1, k
        chi1(ik) = chisq_random(1.0_dp)
      end do
      if (reml) then
        chi_sum = 0.0_dp
      else
        chi_sum = sum(chi1)
      end if
      best = 1
      do ig = 1, size(lambda)
        num = dot_product(f_n(ig, :), chi1)
        den = dot_product(f_d(ig, :), chi1) + chi_k
        if (den <= 0.0_dp .or. 1.0_dp + num / den <= 0.0_dp) then
          lr = -huge(1.0_dp)
        else
          lr = real(n0, dp) * log1p_safe(num / den) - sumlog(ig)
        end if
        if (lr >= result%statistic(is)) then
          result%statistic(is) = lr
          best = ig
        else
          exit
        end if
      end do
      result%lambda(is) = lambda(best)
      if (.not. reml) then
        den = chi_sum + chi_k
        if (q > 0 .and. den > 0.0_dp) then
          result%statistic(is) = result%statistic(is) + &
            real(n, dp) * log1p_safe(chisq_random(real(q, dp)) / den)
        end if
      end if
    end do
  end subroutine rlrsim_kernel

  subroutine regular_lambda_grid(lambda0, gridlength, log_grid_lo, log_grid_hi, grid)
    real(dp), intent(in) :: lambda0 !! Null variance ratio; zero requests the standard grid.
    integer, intent(in) :: gridlength !! Requested total grid length; must be at least 2.
    real(dp), intent(in) :: log_grid_lo !! Natural-log lower endpoint for the positive grid.
    real(dp), intent(in) :: log_grid_hi !! Natural-log upper endpoint for the positive grid.
    real(dp), allocatable, intent(out) :: grid(:) !! Increasing lambda grid beginning at zero and containing lambda0 when positive.
    real(dp), allocatable :: leftgrid(:)
    real(dp), allocatable :: rightgrid(:)
    real(dp) :: base
    real(dp) :: leftdistance
    real(dp) :: rightdistance
    real(dp) :: leftratio
    real(dp) :: eps10
    integer :: i
    integer :: leftlength
    integer :: rightlength

    if (lambda0 == 0.0_dp) then
      allocate(grid(gridlength))
      grid(1) = 0.0_dp
      do i = 2, gridlength
        grid(i) = exp(log_grid_lo + real(i - 2, dp) * (log_grid_hi - log_grid_lo) / real(gridlength - 2, dp))
      end do
      return
    end if
    if (gridlength < 4) error stop "regular_lambda_grid: gridlength must be at least 4 when lambda0 is positive"

    leftratio = min(max(log(lambda0) / (log_grid_hi - log_grid_lo), 0.2_dp), 0.8_dp)
    leftlength = max(nint(leftratio * real(gridlength, dp)) - 1, 2)
    leftdistance = lambda0 - exp(log_grid_lo)
    eps10 = 10.0_dp * epsilon(1.0_dp)
    if (leftdistance < real(leftlength, dp) * eps10) then
      leftlength = max(nint(leftdistance / eps10), 2)
    end if
    rightlength = gridlength - leftlength
    allocate(leftgrid(leftlength - 1), rightgrid(rightlength - 1))

    if (abs(leftdistance - 1.0_dp) < 0.3_dp) then
      do i = 1, leftlength - 1
        leftgrid(i) = exp(log_grid_lo) + real(i - 1, dp) * leftdistance / real(leftlength, dp)
      end do
    else if (leftdistance > 1.0_dp) then
      base = leftdistance ** (1.0_dp / real(leftlength, dp))
      do i = 1, leftlength - 1
        leftgrid(i) = lambda0 - (leftdistance ** (real(leftlength - i + 1, dp) / real(leftlength, dp)) - base)
      end do
    else
      base = leftdistance ** real(leftlength, dp)
      do i = 1, leftlength - 1
        leftgrid(i) = lambda0 - (leftdistance ** real(i, dp) - base)
      end do
    end if

    rightdistance = exp(log_grid_hi) - lambda0
    base = rightdistance ** (1.0_dp / real(rightlength, dp))
    do i = 1, rightlength - 1
      rightgrid(i) = lambda0 + rightdistance ** (real(i + 1, dp) / real(rightlength, dp)) - base
    end do

    allocate(grid(gridlength))
    grid(1) = 0.0_dp
    grid(2:leftlength) = leftgrid
    grid(leftlength + 1) = lambda0
    grid(leftlength + 2:) = rightgrid
  end subroutine regular_lambda_grid

  pure function cumulative_fraction(x) result(frac)
    real(dp), intent(in) :: x(:) !! Nonnegative vector whose cumulative fraction of the total is requested.
    real(dp) :: frac(size(x))
    real(dp) :: total
    integer :: i

    total = sum(x)
    if (total <= 0.0_dp) then
      frac = 0.0_dp
      return
    end if
    frac(1) = x(1) / total
    do i = 2, size(x)
      frac(i) = frac(i - 1) + x(i) / total
    end do
  end function cumulative_fraction

  pure elemental real(dp) function log1p_safe(x) result(y)
    real(dp), intent(in) :: x !! Increment in log(1+x); must exceed -1.

    if (abs(x) > 1.0e-4_dp) then
      y = log(1.0_dp + x)
    else
      y = x * (1.0_dp + x * (-0.5_dp + x * (1.0_dp / 3.0_dp + x * (-0.25_dp + 0.2_dp * x))))
    end if
  end function log1p_safe

end module rlrsim_api
