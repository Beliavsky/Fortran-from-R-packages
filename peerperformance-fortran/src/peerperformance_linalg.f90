! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from PeerPerformance 2.4.0, copyright 2012-2023 David Ardia and Kris Boudt.
module peerperformance_linalg
  use peerperformance_kinds, only: dp
  implicit none
  private
  public :: solve_linear, invert_matrix, sample_covariance, long_run_covariance
  public :: ols_fit, column_means, outer_product, symmetrize

contains

  subroutine solve_linear(a, b, x, ok)
    real(dp), intent(in) :: a(:,:), b(:)
    real(dp), intent(out) :: x(size(b))
    logical, intent(out) :: ok
    real(dp), allocatable :: aa(:,:), bb(:)
    real(dp) :: pivot, factor, tmp, scale
    integer :: n, i, j, k, p
    n = size(b)
    x = 0.0_dp
    ok = size(a,1) == n .and. size(a,2) == n
    if (.not. ok) return
    allocate(aa(n,n),bb(n))
    aa = a
    bb = b
    scale = max(1.0_dp,maxval(abs(aa)))
    do k = 1, n
      p = k
      do i = k+1, n
        if (abs(aa(i,k)) > abs(aa(p,k))) p = i
      end do
      if (abs(aa(p,k)) <= 100.0_dp*epsilon(1.0_dp)*scale) then
        ok = .false.
        return
      end if
      if (p /= k) then
        do j = k, n
          tmp = aa(k,j)
          aa(k,j) = aa(p,j)
          aa(p,j) = tmp
        end do
        tmp = bb(k)
        bb(k) = bb(p)
        bb(p) = tmp
      end if
      pivot = aa(k,k)
      do i = k+1, n
        factor = aa(i,k)/pivot
        aa(i,k) = 0.0_dp
        if (k < n) aa(i,k+1:n) = aa(i,k+1:n)-factor*aa(k,k+1:n)
        bb(i) = bb(i)-factor*bb(k)
      end do
    end do
    do i = n, 1, -1
      if (i < n) then
        x(i) = (bb(i)-dot_product(aa(i,i+1:n),x(i+1:n)))/aa(i,i)
      else
        x(i) = bb(i)/aa(i,i)
      end if
    end do
  end subroutine solve_linear

  subroutine invert_matrix(a, ainv, ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: ainv(size(a,1),size(a,2))
    logical, intent(out) :: ok
    real(dp), allocatable :: e(:), x(:)
    integer :: n, j
    n = size(a,1)
    ainv = 0.0_dp
    ok = size(a,2) == n
    if (.not. ok) return
    allocate(e(n),x(n))
    do j = 1, n
      e = 0.0_dp
      e(j) = 1.0_dp
      call solve_linear(a,e,x,ok)
      if (.not. ok) return
      ainv(:,j) = x
    end do
    call symmetrize(ainv)
  end subroutine invert_matrix

  pure function column_means(x) result(mu)
    real(dp), intent(in) :: x(:,:)
    real(dp) :: mu(size(x,2))
    if (size(x,1) > 0) then
      mu = sum(x,dim=1)/real(size(x,1),dp)
    else
      mu = 0.0_dp
    end if
  end function column_means

  pure function sample_covariance(x) result(cov)
    real(dp), intent(in) :: x(:,:)
    real(dp) :: cov(size(x,2),size(x,2))
    real(dp) :: mu(size(x,2)), d(size(x,2))
    integer :: i, n
    n = size(x,1)
    cov = 0.0_dp
    if (n <= 1) return
    mu = column_means(x)
    do i = 1, n
      d = x(i,:)-mu
      cov = cov+outer_product(d,d)
    end do
    cov = cov/real(n-1,dp)
  end function sample_covariance

  pure function outer_product(x, y) result(a)
    real(dp), intent(in) :: x(:), y(:)
    real(dp) :: a(size(x),size(y))
    a = spread(x,2,size(y))*spread(y,1,size(x))
  end function outer_product

  pure real(dp) function parzen_kernel(x) result(value)
    real(dp), intent(in) :: x
    if (abs(x) <= 0.5_dp) then
      value = 1.0_dp-6.0_dp*x*x+6.0_dp*abs(x)**3
    else if (abs(x) <= 1.0_dp) then
      value = 2.0_dp*(1.0_dp-abs(x))**3
    else
      value = 0.0_dp
    end if
  end function parzen_kernel

  pure subroutine ar1_parameters(x, rho, sigma)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: rho, sigma
    real(dp) :: mx, my, den, resid
    integer :: n, i
    n = size(x)
    rho = 0.0_dp
    sigma = 0.0_dp
    if (n <= 3) return
    mx = sum(x(1:n-1))/real(n-1,dp)
    my = sum(x(2:n))/real(n-1,dp)
    den = sum((x(1:n-1)-mx)**2)
    if (den <= tiny(1.0_dp)) return
    rho = sum((x(1:n-1)-mx)*(x(2:n)-my))/den
    rho = min(0.999_dp,max(-0.999_dp,rho))
    resid = 0.0_dp
    do i = 2, n
      resid = resid+(x(i)-my-rho*(x(i-1)-mx))**2
    end do
    sigma = sqrt(max(0.0_dp,resid/real(max(1,n-2),dp)))
  end subroutine ar1_parameters

  function long_run_covariance(v, hac) result(psi)
    real(dp), intent(in) :: v(:,:)
    logical, intent(in) :: hac
    real(dp) :: psi(size(v,2),size(v,2))
    real(dp) :: gamma(size(v,2),size(v,2))
    real(dp) :: alpha_hat, numerator, denominator, rho, sigma, bandwidth
    integer :: n, p, i, j, t
    n = size(v,1)
    p = size(v,2)
    psi = 0.0_dp
    if (n <= 1) return
    if (.not. hac) then
      psi = sample_covariance(v)
      return
    end if
    numerator = 0.0_dp
    denominator = 0.0_dp
    do i = 1, p
      call ar1_parameters(v(:,i),rho,sigma)
      numerator = numerator+4.0_dp*rho*rho*sigma**4/max((1.0_dp-rho)**8,tiny(1.0_dp))
      denominator = denominator+sigma**4/max((1.0_dp-rho)**4,tiny(1.0_dp))
    end do
    if (denominator > tiny(1.0_dp)) then
      alpha_hat = max(0.0_dp,numerator/denominator)
    else
      alpha_hat = 1.0_dp
    end if
    bandwidth = min(real(n-1,dp),2.6614_dp*(max(alpha_hat*real(n,dp),tiny(1.0_dp)))**0.2_dp)
    do t = 1, n
      psi = psi+outer_product(v(t,:),v(t,:))
    end do
    psi = psi/real(n,dp)
    j = 1
    do while (real(j,dp) < bandwidth .and. j < n)
      gamma = 0.0_dp
      do t = j+1, n
        gamma = gamma+outer_product(v(t,:),v(t-j,:))
      end do
      gamma = gamma/real(n,dp)
      psi = psi+parzen_kernel(real(j,dp)/bandwidth)*(gamma+transpose(gamma))
      j = j+1
    end do
    if (n > 4) psi = real(n,dp)/real(n-4,dp)*psi
    call symmetrize(psi)
  end function long_run_covariance

  subroutine ols_fit(design, y, beta, residuals, covariance, standard_error, &
                     tstat, pvalue, hac, ok)
    use peerperformance_math, only: two_sided_t_pvalue
    real(dp), intent(in) :: design(:,:), y(:)
    real(dp), intent(out) :: beta(size(design,2))
    real(dp), intent(out) :: residuals(size(y))
    real(dp), intent(out) :: covariance(size(design,2),size(design,2))
    real(dp), intent(out) :: standard_error(size(design,2))
    real(dp), intent(out) :: tstat(size(design,2)), pvalue(size(design,2))
    logical, intent(in) :: hac
    logical, intent(out) :: ok
    real(dp), allocatable :: xtx(:,:), xtx_inv(:,:), xty(:), scores(:,:), meat(:,:)
    real(dp) :: sigma2
    integer :: n, p, i
    n = size(design,1)
    p = size(design,2)
    beta = 0.0_dp
    residuals = 0.0_dp
    covariance = 0.0_dp
    standard_error = 0.0_dp
    tstat = 0.0_dp
    pvalue = 1.0_dp
    ok = size(y) == n .and. n > p
    if (.not. ok) return
    allocate(xtx(p,p),xtx_inv(p,p),xty(p),scores(n,p),meat(p,p))
    xtx = matmul(transpose(design),design)
    xty = matmul(transpose(design),y)
    call solve_linear(xtx,xty,beta,ok)
    if (.not. ok) return
    call invert_matrix(xtx,xtx_inv,ok)
    if (.not. ok) return
    residuals = y-matmul(design,beta)
    if (hac) then
      do i = 1, n
        scores(i,:) = design(i,:)*residuals(i)
      end do
      meat = long_run_covariance(scores,.true.)
      covariance = matmul(xtx_inv,matmul(real(n,dp)*meat,xtx_inv))
    else
      sigma2 = sum(residuals*residuals)/real(n-p,dp)
      covariance = sigma2*xtx_inv
    end if
    call symmetrize(covariance)
    do i = 1, p
      standard_error(i) = sqrt(max(0.0_dp,covariance(i,i)))
      if (standard_error(i) > tiny(1.0_dp)) then
        tstat(i) = beta(i)/standard_error(i)
        pvalue(i) = two_sided_t_pvalue(tstat(i),real(n-p,dp))
      else
        ok = .false.
        return
      end if
    end do
  end subroutine ols_fit

  pure subroutine symmetrize(a)
    real(dp), intent(inout) :: a(:,:)
    real(dp) :: value
    integer :: i, j, n
    n = min(size(a,1),size(a,2))
    do i = 1, n
      do j = i+1, n
        value = 0.5_dp*(a(i,j)+a(j,i))
        a(i,j) = value
        a(j,i) = value
      end do
    end do
  end subroutine symmetrize

end module peerperformance_linalg
