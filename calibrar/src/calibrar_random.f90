! SPDX-License-Identifier: GPL-2.0-only
module calibrar_random
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
  use calibrar_kinds, only : dp, pi_dp
  implicit none
  private
  public :: set_seed, rand_uniform, rand_normal, rand_exponential
  public :: rtnorm_sample, rtnorm_matrix, dmvnorm_pdf, gaussian_kernel_2d

contains

  subroutine set_seed(seed)
    integer, intent(in) :: seed
    integer :: n, i
    integer, allocatable :: put(:)
    call random_seed(size=n)
    allocate(put(n))
    do i = 1, n
      put(i) = modulo(abs(seed) + 104729*i + 7919*i*i, huge(1)-1)
      if (put(i) == 0) put(i) = i
    end do
    call random_seed(put=put)
  end subroutine set_seed

  function rand_uniform() result(u)
    real(dp) :: u
    call random_number(u)
  end function rand_uniform

  function rand_normal() result(z)
    real(dp) :: z, u1, u2
    call random_number(u1)
    call random_number(u2)
    u1 = max(u1, tiny(1.0_dp))
    z = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi_dp*u2)
  end function rand_normal

  function rand_exponential(rate) result(x)
    real(dp), intent(in) :: rate
    real(dp) :: x, u
    call random_number(u)
    u = max(u, tiny(1.0_dp))
    x = -log(u)/rate
  end function rand_exponential

  function rtnorm_sample(mu, sd, lower, upper) result(x)
    real(dp), intent(in) :: mu, sd, lower, upper
    real(dp) :: x, z, a, b, u, rho
    integer :: it
    if (sd <= 0.0_dp) then
      x = max(lower, min(upper, mu))
      return
    end if
    a = (lower-mu)/sd
    b = (upper-mu)/sd
    if (a >= b) then
      x = ieee_value(x, ieee_quiet_nan)
      return
    end if
    ! Direct rejection is efficient for central and moderately truncated cases.
    if ((a < 0.0_dp .and. b > 0.0_dp) .or. (.not. ieee_is_finite(a) .and. b > 0.0_dp) .or. &
        (a < 0.0_dp .and. .not. ieee_is_finite(b))) then
      do it = 1, 100000
        z = rand_normal()
        if (z >= a .and. z <= b) then
          x = mu + sd*z
          return
        end if
      end do
    end if
    ! Robert-style exponential proposal for one-sided tails.
    if (a >= 0.0_dp .and. .not. ieee_is_finite(b)) then
      rho = 0.5_dp*(a + sqrt(a*a + 4.0_dp))
      do it = 1, 100000
        z = a + rand_exponential(rho)
        u = rand_uniform()
        if (u <= exp(-0.5_dp*(z-rho)**2)) then
          x = mu + sd*z
          return
        end if
      end do
    else if (.not. ieee_is_finite(a) .and. b <= 0.0_dp) then
      rho = 0.5_dp*(-b + sqrt(b*b + 4.0_dp))
      do it = 1, 100000
        z = -b + rand_exponential(rho)
        u = rand_uniform()
        if (u <= exp(-0.5_dp*(z-rho)**2)) then
          x = mu - sd*z
          return
        end if
      end do
    end if
    ! Uniform rejection, matching the fallback structure of calibrar::.rtnormx.
    do it = 1, 200000
      if (ieee_is_finite(a) .and. ieee_is_finite(b)) then
        z = a + (b-a)*rand_uniform()
        if (a > 0.0_dp) then
          rho = exp(0.5_dp*(a*a-z*z))
        else if (b < 0.0_dp) then
          rho = exp(0.5_dp*(b*b-z*z))
        else
          rho = exp(-0.5_dp*z*z)
        end if
        if (rand_uniform() <= rho) then
          x = mu + sd*z
          return
        end if
      else
        z = rand_normal()
        if (z >= a .and. z <= b) then
          x = mu + sd*z
          return
        end if
      end if
    end do
    x = max(lower, min(upper, mu))
  end function rtnorm_sample

  subroutine rtnorm_matrix(n, mean, sd, lower, upper, out)
    integer, intent(in) :: n
    real(dp), intent(in) :: mean(:), sd(:), lower(:), upper(:)
    real(dp), intent(out) :: out(:,:)
    integer :: i, j
    if (size(out,1) /= size(mean) .or. size(out,2) /= n) error stop "rtnorm_matrix: size mismatch"
    do j = 1, n
      do i = 1, size(mean)
        if (sd(i) <= 0.0_dp) then
          out(i,j) = max(lower(i), min(upper(i), mean(i)))
        else if (truncation_mass(mean(i),sd(i),lower(i),upper(i)) > 0.025_dp) then
          out(i,j) = rtnorm_sample(mean(i), sd(i), lower(i), upper(i))
        else
          out(i,j) = mean(i) + sd(i)*rand_normal()
          out(i,j) = max(lower(i), min(upper(i), out(i,j)))
        end if
      end do
    end do
  end subroutine rtnorm_matrix

  function truncation_mass(mu,sd,lower,upper) result(pout)
    real(dp), intent(in) :: mu,sd,lower,upper
    real(dp) :: pout,a,b,pa,pb
    if (sd <= 0.0_dp) then
      pout=0.0_dp
      return
    end if
    a=(lower-mu)/sd
    b=(upper-mu)/sd
    if (ieee_is_finite(a)) then
      pa=0.5_dp*(1.0_dp+erf(a/sqrt(2.0_dp)))
    else if (a < 0.0_dp) then
      pa=0.0_dp
    else
      pa=1.0_dp
    end if
    if (ieee_is_finite(b)) then
      pb=0.5_dp*(1.0_dp+erf(b/sqrt(2.0_dp)))
    else if (b > 0.0_dp) then
      pb=1.0_dp
    else
      pb=0.0_dp
    end if
    pout=max(0.0_dp,min(1.0_dp,pa+1.0_dp-pb))
  end function truncation_mass

  function dmvnorm_pdf(x, mean, sigma, log_pdf) result(v)
    real(dp), intent(in) :: x(:), mean(:), sigma(:,:)
    logical, intent(in), optional :: log_pdf
    real(dp) :: v
    real(dp), allocatable :: l(:,:), y(:), d(:)
    real(dp) :: logv
    integer :: n, i, j, k
    logical :: want_log
    n = size(x)
    if (size(mean) /= n .or. size(sigma,1) /= n .or. size(sigma,2) /= n) error stop "dmvnorm_pdf: size mismatch"
    allocate(l(n,n), y(n), d(n))
    l = 0.0_dp
    do i = 1, n
      do j = 1, i
        v = sigma(i,j)
        do k = 1, j-1
          v = v - l(i,k)*l(j,k)
        end do
        if (i == j) then
          if (v <= 0.0_dp) then
            v = -huge(1.0_dp)
            if (.not. present(log_pdf)) v = 0.0_dp
            return
          end if
          l(i,j) = sqrt(v)
        else
          l(i,j) = v/l(j,j)
        end if
      end do
    end do
    d = x - mean
    do i = 1, n
      y(i) = d(i)
      do k = 1, i-1
        y(i) = y(i) - l(i,k)*y(k)
      end do
      y(i) = y(i)/l(i,i)
    end do
    logv = -sum(log([(l(i,i), i=1,n)])) - 0.5_dp*real(n,dp)*log(2.0_dp*pi_dp) - 0.5_dp*sum(y*y)
    want_log = .false.
    if (present(log_pdf)) want_log = log_pdf
    if (want_log) then
      v = logv
    else
      v = exp(logv)
    end if
  end function dmvnorm_pdf

  subroutine gaussian_kernel_2d(mean, sigma, lower, upper, nx, ny, xgrid, ygrid, z)
    real(dp), intent(in) :: mean(2), sigma(2,2), lower(2), upper(2)
    integer, intent(in) :: nx, ny
    real(dp), intent(out) :: xgrid(nx), ygrid(ny), z(nx,ny)
    integer :: i, j
    real(dp) :: pt(2)
    if (nx < 2 .or. ny < 2) error stop "gaussian_kernel_2d: grid dimensions must be >= 2"
    do i = 1, nx
      xgrid(i) = lower(1) + real(i-1,dp)*(upper(1)-lower(1))/real(nx-1,dp)
    end do
    do j = 1, ny
      ygrid(j) = lower(2) + real(j-1,dp)*(upper(2)-lower(2))/real(ny-1,dp)
    end do
    do j = 1, ny
      do i = 1, nx
        pt = [xgrid(i), ygrid(j)]
        z(i,j) = dmvnorm_pdf(pt, mean, sigma)
      end do
    end do
  end subroutine gaussian_kernel_2d
end module calibrar_random
