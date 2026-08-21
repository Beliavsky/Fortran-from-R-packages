! SPDX-License-Identifier: GPL-2.0-or-later
module tsa_garch_diagnostics
  use tsa_kinds, only : dp
  use tsa_types, only : tsa_test_result
  use tsa_statistics, only : autocorrelation, kurtosis
  use tsa_utils, only : mean_value
  use tseries_types, only : garch_result
  use tseries_linalg, only : invert_matrix
  use tseries_special, only : chi_square_cdf
  implicit none
  private
  public :: gbox_test
contains

  function gbox_test(x, model, lag, absolute) result(res)
    real(dp), intent(in) :: x(:)
    type(garch_result), intent(in) :: model
    integer, intent(in) :: lag
    logical, intent(in), optional :: absolute
    type(tsa_test_result) :: res
    real(dp), allocatable :: m(:,:), e(:,:), jmat(:,:), lambda(:,:), hmat(:,:)
    real(dp), allocatable :: mtm(:,:), invm(:,:), ac(:), acv(:), invh(:,:)
    real(dp), allocatable :: eps(:), hv(:), beta(:)
    real(dp) :: sigma2, kval, tau, nu, factor
    integer :: p, q, n, i, j, k, status
    logical :: abs_method

    abs_method = .false.
    if (present(absolute)) abs_method = absolute
    p = model%p
    q = model%q
    n = min(size(x), size(model%residuals))
    if (lag < 1 .or. n <= lag + max(p,q) .or. &
        .not. allocated(model%conditional_variance)) then
      res%status = 1
      return
    end if

    allocate(eps(n), hv(n))
    eps = model%residuals(:n)
    hv = model%conditional_variance(:n)
    allocate(beta(p))
    if (p > 0) beta = model%coefficients(2+q:1+q+p)

    allocate(m(n,1+q+p))
    m = 0.0_dp
    sigma2 = mean_value(x(:n)**2)
    if (p > 0) then
      m(:,1) = 1.0_dp / max(1.0_dp-sum(beta), tiny(1.0_dp))
    else
      m(:,1) = 1.0_dp
    end if
    do j = 1, q
      m(:,1+j) = sigma2
      if (j < n) m(j+1:,1+j) = x(:n-j)**2
    end do
    do j = 1, p
      m(:,1+q+j) = sigma2
      if (j < n) m(j+1:,1+q+j) = hv(:n-j)
    end do

    if (p > 0) then
      do j = 2, size(m,2)
        do i = 2, n
          do k = 1, min(p,i-1)
            m(i,j) = m(i,j) + beta(k)*m(i-k,j)
          end do
        end do
      end do
    end if
    do j = 1, size(m,2)
      m(:,j) = m(:,j) / max(hv, tiny(1.0_dp))
    end do

    allocate(e(n,lag))
    e = 0.0_dp
    if (abs_method) then
      tau = mean_value(abs(eps))
      do j = 1, lag
        if (j < n) e(j+1:,j) = abs(eps(:n-j)) - tau
      end do
    else
      do j = 1, lag
        if (j < n) e(j+1:,j) = eps(:n-j)**2 - mean_value(eps**2)
      end do
    end if

    allocate(jmat(lag,size(m,2)))
    jmat = matmul(transpose(e),m) / real(n,dp)
    allocate(mtm(size(m,2),size(m,2)), invm(size(m,2),size(m,2)))
    mtm = matmul(transpose(m),m)
    call invert_matrix(mtm,invm,status)
    if (status /= 0) then
      res%status = 2
      return
    end if
    kval = kurtosis(eps) + 2.0_dp
    allocate(lambda(size(m,2),size(m,2)))
    lambda = 2.0_dp*invm

    allocate(hmat(lag,lag))
    if (.not. abs_method) then
      hmat = -matmul(matmul(jmat,lambda),transpose(jmat)) * &
        real(n,dp) / max(2.0_dp*kval,tiny(1.0_dp))
    else
      nu = mean_value(abs(eps)**3)
      lambda = lambda*kval/2.0_dp
      factor = (kval/2.0_dp*tau*tau/4.0_dp - &
        tau*(nu-tau)/2.0_dp) / max((1.0_dp-tau*tau)**2,tiny(1.0_dp))
      hmat = matmul(matmul(jmat,lambda),transpose(jmat))*real(n,dp)*factor
    end if
    do i = 1, lag
      hmat(i,i) = 1.0_dp + hmat(i,i)
    end do

    if (abs_method) then
      call autocorrelation(abs(eps),lag,ac)
    else
      call autocorrelation(eps**2,lag,ac)
    end if
    allocate(invh(lag,lag), acv(lag))
    call invert_matrix(hmat,invh,status)
    if (status /= 0) then
      res%status = 3
      return
    end if
    acv = ac
    res%statistic = real(n,dp)*dot_product(acv,matmul(invh,acv))
    res%df = lag
    res%lag = lag
    res%p_value = 1.0_dp-chi_square_cdf(res%statistic,real(lag,dp))
    res%status = 0
  end function gbox_test
end module tsa_garch_diagnostics
