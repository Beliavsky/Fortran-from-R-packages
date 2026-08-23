module mev_spatial
  use mev_kinds, only: dp, pi
  implicit none
  private
  public :: power_vario, powerexp_cor, schlather_vario
  public :: distg, dgeoaniso, lambda2cov
contains
  pure real(dp) function power_vario(h, alpha, scale) result(v)
    real(dp), intent(in) :: h, alpha, scale
    if (scale <= 0.0_dp .or. alpha < 0.0_dp .or. alpha >= 2.0_dp) then
      v = huge(1.0_dp)
    else
      v = (h / scale)**alpha
    end if
  end function power_vario

  pure real(dp) function powerexp_cor(h, alpha, scale) result(v)
    real(dp), intent(in) :: h, alpha, scale
    if (scale <= 0.0_dp .or. alpha < 0.0_dp .or. alpha >= 2.0_dp) then
      v = 0.0_dp
    else
      v = exp(-(h / scale)**alpha)
    end if
  end function powerexp_cor

  pure real(dp) function schlather_vario(h, alpha, beta, scale) result(v)
    real(dp), intent(in) :: h, alpha, beta, scale
    real(dp) :: powr
    if (alpha <= 0.0_dp .or. alpha >= 2.0_dp .or. beta >= 2.0_dp .or. scale <= 0.0_dp) then
      v = huge(1.0_dp)
      return
    end if
    if (abs(beta) < 1.0e-5_dp) then
      v = log(1.0_dp + (h / scale)**alpha) / log(2.0_dp)
    else
      powr = beta / alpha
      v = ((1.0_dp + (h / scale)**alpha)**powr - 1.0_dp) / (exp(powr * log(2.0_dp)) - 1.0_dp)
    end if
  end function schlather_vario

  subroutine distg(loc, scale, rho, dmat, info)
    real(dp), intent(in) :: loc(:, :), scale, rho
    real(dp), intent(out) :: dmat(:, :)
    integer, intent(out), optional :: info
    integer :: n, i, j
    real(dp) :: a11, a12, a21, a22, dx, dy, tx, ty
    if (present(info)) info = 0
    n = size(loc, 1)
    if (size(loc, 2) /= 2 .or. size(dmat, 1) /= n .or. size(dmat, 2) /= n .or. &
        scale < 1.0_dp .or. abs(rho) > 0.5_dp*pi) then
      if (present(info)) info = 1
      dmat = 0.0_dp
      return
    end if
    a11 = cos(rho)
    a12 = sin(rho)
    a21 = -scale * sin(rho)
    a22 = scale * cos(rho)
    do i = 1, n
      dmat(i, i) = 0.0_dp
      do j = i + 1, n
        dx = loc(i, 1) - loc(j, 1)
        dy = loc(i, 2) - loc(j, 2)
        tx = a11*dx + a12*dy
        ty = a21*dx + a22*dy
        dmat(i, j) = sqrt(tx*tx + ty*ty)
        dmat(j, i) = dmat(i, j)
      end do
    end do
  end subroutine distg

  subroutine dgeoaniso(loc, theta, dmat, info)
    real(dp), intent(in) :: loc(:, :), theta(:)
    real(dp), intent(out) :: dmat(:, :)
    integer, intent(out), optional :: info
    integer :: n, i, j
    real(dp) :: rho, sqr, a11, a12, a21, a22, dx, dy, tx, ty
    if (present(info)) info = 0
    n = size(loc, 1)
    if (size(loc, 2) /= 2 .or. size(theta) /= 2 .or. size(dmat, 1) /= n .or. size(dmat, 2) /= n) then
      if (present(info)) info = 1
      dmat = 0.0_dp
      return
    end if
    rho = 0.5_dp * atan2(theta(1), theta(2))
    sqr = sqrt(1.0_dp + theta(1)*theta(1) + theta(2)*theta(2))
    a11 = sqr*cos(rho)
    a12 = -sqr*sin(rho)
    a21 = sin(rho)/sqr
    a22 = cos(rho)/sqr
    do i = 1, n
      dmat(i, i) = 0.0_dp
      do j = i + 1, n
        dx = loc(i, 1) - loc(j, 1)
        dy = loc(i, 2) - loc(j, 2)
        tx = a11*dx + a12*dy
        ty = a21*dx + a22*dy
        dmat(i, j) = sqrt(tx*tx + ty*ty)
        dmat(j, i) = dmat(i, j)
      end do
    end do
  end subroutine dgeoaniso

  subroutine lambda2cov(lambda, co, suba, subb, cov, info)
    real(dp), intent(in) :: lambda(:, :)
    integer, intent(in) :: co, suba(:), subb(:)
    real(dp), intent(out) :: cov(:, :)
    integer, intent(out), optional :: info
    integer :: i, j, n
    if (present(info)) info = 0
    n = size(lambda, 1)
    if (size(lambda, 2) /= n .or. co < 1 .or. co > n .or. &
        size(cov, 1) /= size(suba) .or. size(cov, 2) /= size(subb) .or. &
        any(suba < 1) .or. any(suba > n) .or. any(subb < 1) .or. any(subb > n) .or. &
        any(suba == co) .or. any(subb == co)) then
      if (present(info)) info = 1
      cov = 0.0_dp
      return
    end if
    do i = 1, size(suba)
      do j = 1, size(subb)
        cov(i, j) = 2.0_dp * (lambda(co, suba(i)) + lambda(co, subb(j)) - lambda(suba(i), subb(j)))
      end do
    end do
  end subroutine lambda2cov
end module mev_spatial
