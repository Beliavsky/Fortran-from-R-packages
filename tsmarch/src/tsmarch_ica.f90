! SPDX-License-Identifier: GPL-2.0-only
module tsmarch_ica
  use ghyp_kinds, only : dp, i8
  use tsmarch_types, only : ica_result, tsm_success, tsm_invalid_argument, tsm_no_convergence, tsm_numerical_failure
  use tsmarch_linalg
  implicit none
  private
  public :: fastica, radical, whiten_data

contains

  subroutine whiten_data(x, centered, whitening, mean, ok)
    real(dp), intent(in) :: x(:, :)
    real(dp), allocatable, intent(out) :: centered(:, :), whitening(:, :), mean(:)
    logical, intent(out) :: ok
    real(dp), allocatable :: cov(:, :)
    integer :: n
    n = size(x, 1)
    if (n < 2 .or. size(x, 2) < 1) then
      ok = .false.
      allocate(centered(0, 0), whitening(0, 0), mean(0))
      return
    end if
    mean = column_mean(x)
    allocate(centered(n, size(x, 2)))
    centered = x - spread(mean, 1, n)
    cov = sample_covariance(centered)
    whitening = symmetric_inverse_sqrt(cov, ok, 1.0e-10_dp)
    if (ok) centered = matmul(centered, whitening)
  end subroutine whiten_data

  function fastica(x, max_iterations, tolerance, seed, symmetric) result(out)
    real(dp), intent(in) :: x(:, :)
    integer, intent(in), optional :: max_iterations
    real(dp), intent(in), optional :: tolerance
    integer(i8), intent(in), optional :: seed
    logical, intent(in), optional :: symmetric
    type(ica_result) :: out
    real(dp), allocatable :: xwhite(:, :), whitening(:, :), mean(:), w(:, :), wnew(:, :)
    real(dp), allocatable :: y(:, :), g(:, :), gp(:), cross(:, :), invmix(:, :)
    real(dp) :: tol, convergence
    integer :: maxit, iter, n, m, i
    logical :: ok, symm

    n = size(x, 1)
    m = size(x, 2)
    if (n < max(10, m + 1) .or. m < 1) then
      out%status = tsm_invalid_argument
      out%message = 'FastICA requires more observations than variables'
      return
    end if
    maxit = 500
    if (present(max_iterations)) maxit = max_iterations
    tol = 1.0e-7_dp
    if (present(tolerance)) tol = tolerance
    symm = .true.
    if (present(symmetric)) symm = symmetric
    call whiten_data(x, xwhite, whitening, mean, ok)
    if (.not. ok) then
      out%status = tsm_numerical_failure
      out%message = 'whitening failed'
      return
    end if
    if (present(seed)) then
      w = random_normal_matrix(m, m, seed)
    else
      w = random_normal_matrix(m, m)
    end if
    if (symm) then
      call symmetric_decorrelation(w, ok)
      if (.not. ok) then
        out%status = tsm_numerical_failure
        out%message = 'initial orthogonalization failed'
        return
      end if
      allocate(wnew(m, m), y(n, m), g(n, m), gp(m), cross(m, m))
      do iter = 1, maxit
        y = matmul(xwhite, transpose(w))
        g = tanh(y)
        do i = 1, m
          gp(i) = sum(1.0_dp - g(:, i) ** 2) / real(n, dp)
        end do
        wnew = matmul(transpose(g), xwhite) / real(n, dp) - spread(gp, 2, m) * w
        call symmetric_decorrelation(wnew, ok)
        if (.not. ok) exit
        cross = matmul(wnew, transpose(w))
        convergence = maxval(abs(abs([(cross(i, i), i = 1, m)]) - 1.0_dp))
        w = wnew
        if (convergence < tol) exit
      end do
    else
      call fastica_deflation(xwhite, w, maxit, tol, iter, ok)
    end if
    invmix = matmul(w, whitening)
    call matrix_inverse_general(invmix, out%mixing, ok)
    if (.not. ok) then
      out%status = tsm_numerical_failure
      out%message = 'estimated unmixing matrix is singular'
      return
    end if
    out%unmixing = invmix
    out%whitening = whitening
    out%components = matmul(x - spread(mean, 1, n), transpose(invmix))
    out%mean = mean
    out%iterations = iter
    if (iter <= maxit .and. ok) then
      out%status = tsm_success
      out%message = 'ok'
    else
      out%status = tsm_no_convergence
      out%message = 'FastICA reached the iteration limit'
    end if
  end function fastica

  function radical(x, sweeps, angle_grid, seed) result(out)
    real(dp), intent(in) :: x(:, :)
    integer, intent(in), optional :: sweeps, angle_grid
    integer(i8), intent(in), optional :: seed
    type(ica_result) :: out
    real(dp), allocatable :: xwhite(:, :), whitening(:, :), mean(:), rotation(:, :), pair(:, :)
    real(dp), allocatable :: invmix(:, :)
    real(dp) :: angle, best_angle, score, best_score, c, s
    integer :: nsweep, ngrid, sweep, i, j, k, n, m
    logical :: ok

    n = size(x, 1)
    m = size(x, 2)
    if (n < max(20, 2 * m) .or. m < 1) then
      out%status = tsm_invalid_argument
      out%message = 'RADICAL requires at least 20 observations and more observations than variables'
      return
    end if
    nsweep = 5
    if (present(sweeps)) nsweep = max(1, sweeps)
    ngrid = 61
    if (present(angle_grid)) ngrid = max(11, angle_grid)
    if (present(seed)) call set_random_seed(seed)
    call whiten_data(x, xwhite, whitening, mean, ok)
    if (.not. ok) then
      out%status = tsm_numerical_failure
      out%message = 'whitening failed'
      return
    end if
    allocate(rotation(m, m), pair(n, 2))
    rotation = 0.0_dp
    do i = 1, m
      rotation(i, i) = 1.0_dp
    end do
    do sweep = 1, nsweep
      do i = 1, m - 1
        do j = i + 1, m
          pair(:, 1) = xwhite(:, i)
          pair(:, 2) = xwhite(:, j)
          best_score = huge(1.0_dp)
          best_angle = 0.0_dp
          do k = 0, ngrid - 1
            angle = -0.25_dp * acos(-1.0_dp) + 0.5_dp * acos(-1.0_dp) * real(k, dp) / real(ngrid - 1, dp)
            c = cos(angle)
            s = sin(angle)
            score = spacing_entropy(c * pair(:, 1) + s * pair(:, 2)) + &
              spacing_entropy(-s * pair(:, 1) + c * pair(:, 2))
            if (score < best_score) then
              best_score = score
              best_angle = angle
            end if
          end do
          c = cos(best_angle)
          s = sin(best_angle)
          xwhite(:, i) = c * pair(:, 1) + s * pair(:, 2)
          xwhite(:, j) = -s * pair(:, 1) + c * pair(:, 2)
          call rotate_rows(rotation, i, j, c, s)
        end do
      end do
    end do
    invmix = matmul(rotation, whitening)
    call matrix_inverse_general(invmix, out%mixing, ok)
    if (.not. ok) then
      out%status = tsm_numerical_failure
      out%message = 'RADICAL unmixing matrix is singular'
      return
    end if
    out%unmixing = invmix
    out%whitening = whitening
    out%components = matmul(x - spread(mean, 1, n), transpose(invmix))
    out%mean = mean
    out%iterations = nsweep
    out%status = tsm_success
    out%message = 'ok'
  end function radical

  subroutine symmetric_decorrelation(w, ok)
    real(dp), intent(inout) :: w(:, :)
    logical, intent(out) :: ok
    real(dp), allocatable :: invsqrt(:, :)
    invsqrt = symmetric_inverse_sqrt(matmul(w, transpose(w)), ok, 1.0e-12_dp)
    if (ok) w = matmul(invsqrt, w)
  end subroutine symmetric_decorrelation

  subroutine fastica_deflation(x, w, maxit, tol, iterations, ok)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(inout) :: w(:, :)
    integer, intent(in) :: maxit
    real(dp), intent(in) :: tol
    integer, intent(out) :: iterations
    logical, intent(out) :: ok
    real(dp), allocatable :: wi(:), old(:), y(:), g(:)
    real(dp) :: gp, normw
    integer :: component, iter, j, n, m
    n = size(x, 1)
    m = size(x, 2)
    ok = .true.
    iterations = 0
    do component = 1, m
      allocate(wi(m), old(m), y(n), g(n))
      wi = w(component, :)
      normw = sqrt(sum(wi ** 2))
      wi = wi / max(normw, tiny(1.0_dp))
      do iter = 1, maxit
        old = wi
        y = matmul(x, wi)
        g = tanh(y)
        gp = sum(1.0_dp - g ** 2) / real(n, dp)
        wi = matmul(transpose(x), g) / real(n, dp) - gp * wi
        do j = 1, component - 1
          wi = wi - dot_product(wi, w(j, :)) * w(j, :)
        end do
        normw = sqrt(sum(wi ** 2))
        if (normw <= tiny(1.0_dp)) then
          ok = .false.
          return
        end if
        wi = wi / normw
        if (abs(abs(dot_product(wi, old)) - 1.0_dp) < tol) exit
      end do
      w(component, :) = wi
      iterations = max(iterations, iter)
      deallocate(wi, old, y, g)
    end do
  end subroutine fastica_deflation

  subroutine rotate_rows(a, i, j, c, s)
    real(dp), intent(inout) :: a(:, :)
    integer, intent(in) :: i, j
    real(dp), intent(in) :: c, s
    real(dp), allocatable :: ri(:), rj(:)
    ri = a(i, :)
    rj = a(j, :)
    a(i, :) = c * ri + s * rj
    a(j, :) = -s * ri + c * rj
  end subroutine rotate_rows

  function spacing_entropy(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value, key, spacing
    real(dp), allocatable :: sorted(:)
    integer :: n, i, j, m
    n = size(x)
    allocate(sorted(n))
    sorted = x
    do i = 2, n
      key = sorted(i)
      j = i - 1
      do while (j >= 1)
        if (sorted(j) <= key) exit
        sorted(j + 1) = sorted(j)
        j = j - 1
      end do
      sorted(j + 1) = key
    end do
    m = max(1, int(sqrt(real(n, dp))))
    value = 0.0_dp
    do i = m + 1, n - m
      spacing = max(sorted(i + m) - sorted(i - m), tiny(1.0_dp))
      value = value + log(real(n, dp) * spacing / real(2 * m, dp))
    end do
    value = value / real(max(1, n - 2 * m), dp)
  end function spacing_entropy

end module tsmarch_ica
