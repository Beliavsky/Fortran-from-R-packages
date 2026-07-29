! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module dowd_portfolio
  use dowd_kinds, only: dp
  use dowd_math, only: normal_pdf, normal_quantile
  use dowd_risk, only: cornish_fisher_var
  implicit none
  private

  public :: portfolio_mean, portfolio_variance
  public :: variance_covariance_var, variance_covariance_es
  public :: adjusted_variance_covariance_var, adjusted_variance_covariance_es
  public :: normal_var_hotspots, normal_es_hotspots
  public :: adjusted_var_hotspots, adjusted_es_hotspots
  public :: symmetric_eigen_jacobi, pca_prelim, pca_var, pca_es

contains

  pure subroutine validate_dimensions(covariance, mu, positions)
    real(dp), intent(in) :: covariance(:,:), mu(:), positions(:)
    integer :: n
    n = size(positions)
    if (size(mu) /= n .or. size(covariance,1) /= n .or. size(covariance,2) /= n) &
      error stop "portfolio inputs have incompatible dimensions"
  end subroutine validate_dimensions

  pure real(dp) function portfolio_mean(mu, positions) result(value)
    real(dp), intent(in) :: mu(:), positions(:)
    if (size(mu) /= size(positions)) error stop "portfolio_mean: incompatible dimensions"
    value = dot_product(mu,positions)
  end function portfolio_mean

  pure real(dp) function portfolio_variance(covariance, positions) result(value)
    real(dp), intent(in) :: covariance(:,:), positions(:)
    real(dp), allocatable :: temp(:)
    if (size(covariance,1) /= size(positions) .or. size(covariance,2) /= size(positions)) &
      error stop "portfolio_variance: incompatible dimensions"
    temp = matmul(covariance,positions)
    value = max(0.0_dp,dot_product(positions,temp))
  end function portfolio_variance

  pure real(dp) function variance_covariance_var(covariance, mu, positions, cl, hp) result(value)
    real(dp), intent(in) :: covariance(:,:), mu(:), positions(:), cl, hp
    real(dp) :: mean_p, sd_p
    call validate_dimensions(covariance,mu,positions)
    if (cl <= 0.0_dp .or. cl >= 1.0_dp .or. hp <= 0.0_dp) error stop "variance_covariance_var: invalid input"
    mean_p = portfolio_mean(mu,positions)
    sd_p = sqrt(portfolio_variance(covariance,positions))
    value = -mean_p*hp-normal_quantile(1.0_dp-cl)*sd_p*sqrt(hp)
  end function variance_covariance_var

  pure real(dp) function variance_covariance_es(covariance, mu, positions, cl, hp) result(value)
    real(dp), intent(in) :: covariance(:,:), mu(:), positions(:), cl, hp
    real(dp) :: mean_p, sd_p, z
    call validate_dimensions(covariance,mu,positions)
    if (cl <= 0.0_dp .or. cl >= 1.0_dp .or. hp <= 0.0_dp) error stop "variance_covariance_es: invalid input"
    mean_p = portfolio_mean(mu,positions)
    sd_p = sqrt(portfolio_variance(covariance,positions))
    z = normal_quantile(cl)
    value = sd_p*sqrt(hp)*normal_pdf(z)/(1.0_dp-cl)-mean_p*hp
  end function variance_covariance_es

  pure real(dp) function adjusted_variance_covariance_var(covariance, mu, positions, &
      skewness, kurtosis, cl, hp) result(value)
    real(dp), intent(in) :: covariance(:,:), mu(:), positions(:)
    real(dp), intent(in) :: skewness, kurtosis, cl, hp
    real(dp) :: mean_p, sd_p
    call validate_dimensions(covariance,mu,positions)
    mean_p = portfolio_mean(mu,positions)*hp
    sd_p = sqrt(portfolio_variance(covariance,positions)*hp)
    value = cornish_fisher_var(mean_p,sd_p,skewness,kurtosis,cl)
  end function adjusted_variance_covariance_var

  real(dp) function adjusted_variance_covariance_es(covariance, mu, positions, &
      skewness, kurtosis, cl, hp, n_slices) result(value)
    real(dp), intent(in) :: covariance(:,:), mu(:), positions(:)
    real(dp), intent(in) :: skewness, kurtosis, cl, hp
    integer, intent(in), optional :: n_slices
    integer :: i, n
    real(dp) :: p
    n = 1000
    if (present(n_slices)) n = max(100,n_slices)
    value = 0.0_dp
    do i = 1, n
      p = cl+(1.0_dp-cl)*(real(i,dp)-0.5_dp)/real(n,dp)
      value = value+adjusted_variance_covariance_var(covariance,mu,positions,skewness,kurtosis,p,hp)
    end do
    value = value/real(n,dp)
  end function adjusted_variance_covariance_es

  subroutine normal_var_hotspots(covariance, mu, positions, cl, hp, incremental)
    real(dp), intent(in) :: covariance(:,:), mu(:), positions(:), cl, hp
    real(dp), intent(out) :: incremental(:)
    real(dp), allocatable :: reduced(:)
    real(dp) :: full_value
    integer :: i, n
    call validate_dimensions(covariance,mu,positions)
    n = size(positions)
    if (size(incremental) /= n) error stop "normal_var_hotspots: wrong output size"
    allocate(reduced(n))
    full_value = variance_covariance_var(covariance,mu,positions,cl,hp)
    do i = 1, n
      reduced = positions
      reduced(i) = 0.0_dp
      incremental(i) = full_value-variance_covariance_var(covariance,mu,reduced,cl,hp)
    end do
  end subroutine normal_var_hotspots

  subroutine normal_es_hotspots(covariance, mu, positions, cl, hp, incremental)
    real(dp), intent(in) :: covariance(:,:), mu(:), positions(:), cl, hp
    real(dp), intent(out) :: incremental(:)
    real(dp), allocatable :: reduced(:)
    real(dp) :: full_value
    integer :: i, n
    call validate_dimensions(covariance,mu,positions)
    n = size(positions)
    if (size(incremental) /= n) error stop "normal_es_hotspots: wrong output size"
    allocate(reduced(n))
    full_value = variance_covariance_es(covariance,mu,positions,cl,hp)
    do i = 1, n
      reduced = positions
      reduced(i) = 0.0_dp
      incremental(i) = full_value-variance_covariance_es(covariance,mu,reduced,cl,hp)
    end do
  end subroutine normal_es_hotspots

  subroutine adjusted_var_hotspots(covariance, mu, positions, skewness, kurtosis, cl, hp, incremental)
    real(dp), intent(in) :: covariance(:,:), mu(:), positions(:), skewness, kurtosis, cl, hp
    real(dp), intent(out) :: incremental(:)
    real(dp), allocatable :: reduced(:)
    real(dp) :: full_value
    integer :: i, n
    n = size(positions)
    if (size(incremental) /= n) error stop "adjusted_var_hotspots: wrong output size"
    allocate(reduced(n))
    full_value = adjusted_variance_covariance_var(covariance,mu,positions,skewness,kurtosis,cl,hp)
    do i = 1, n
      reduced = positions
      reduced(i) = 0.0_dp
      incremental(i) = full_value-adjusted_variance_covariance_var(covariance,mu,reduced,skewness,kurtosis,cl,hp)
    end do
  end subroutine adjusted_var_hotspots

  subroutine adjusted_es_hotspots(covariance, mu, positions, skewness, kurtosis, cl, hp, incremental)
    real(dp), intent(in) :: covariance(:,:), mu(:), positions(:), skewness, kurtosis, cl, hp
    real(dp), intent(out) :: incremental(:)
    real(dp), allocatable :: reduced(:)
    real(dp) :: full_value
    integer :: i, n
    n = size(positions)
    if (size(incremental) /= n) error stop "adjusted_es_hotspots: wrong output size"
    allocate(reduced(n))
    full_value = adjusted_variance_covariance_es(covariance,mu,positions,skewness,kurtosis,cl,hp)
    do i = 1, n
      reduced = positions
      reduced(i) = 0.0_dp
      incremental(i) = full_value-adjusted_variance_covariance_es(covariance,mu,reduced,skewness,kurtosis,cl,hp)
    end do
  end subroutine adjusted_es_hotspots

  subroutine symmetric_eigen_jacobi(matrix, eigenvalues, eigenvectors)
    real(dp), intent(in) :: matrix(:,:)
    real(dp), intent(out) :: eigenvalues(:), eigenvectors(:,:)
    real(dp), allocatable :: a(:,:)
    real(dp) :: app, aqq, apq, tau, t, c, s, aip, aiq, vip, viq
    real(dp) :: max_off, temp
    integer :: n, p, q, i, iter, max_iter, j
    n = size(matrix,1)
    if (size(matrix,2) /= n .or. size(eigenvalues) /= n .or. &
        size(eigenvectors,1) /= n .or. size(eigenvectors,2) /= n) &
      error stop "symmetric_eigen_jacobi: incompatible dimensions"
    allocate(a(n,n))
    a = 0.5_dp*(matrix+transpose(matrix))
    eigenvectors = 0.0_dp
    do i = 1, n
      eigenvectors(i,i) = 1.0_dp
    end do
    max_iter = max(100,50*n*n)
    do iter = 1, max_iter
      max_off = 0.0_dp
      p = 1
      q = min(2,n)
      do i = 1, n-1
        do j = i+1, n
          if (abs(a(i,j)) > max_off) then
            max_off = abs(a(i,j))
            p = i
            q = j
          end if
        end do
      end do
      if (max_off <= 100.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(a)))) exit
      app = a(p,p)
      aqq = a(q,q)
      apq = a(p,q)
      tau = (aqq-app)/(2.0_dp*apq)
      if (tau >= 0.0_dp) then
        t = 1.0_dp/(tau+sqrt(1.0_dp+tau*tau))
      else
        t = -1.0_dp/(-tau+sqrt(1.0_dp+tau*tau))
      end if
      c = 1.0_dp/sqrt(1.0_dp+t*t)
      s = t*c
      do i = 1, n
        if (i /= p .and. i /= q) then
          aip = a(i,p)
          aiq = a(i,q)
          a(i,p) = c*aip-s*aiq
          a(p,i) = a(i,p)
          a(i,q) = s*aip+c*aiq
          a(q,i) = a(i,q)
        end if
      end do
      a(p,p) = c*c*app-2.0_dp*s*c*apq+s*s*aqq
      a(q,q) = s*s*app+2.0_dp*s*c*apq+c*c*aqq
      a(p,q) = 0.0_dp
      a(q,p) = 0.0_dp
      do i = 1, n
        vip = eigenvectors(i,p)
        viq = eigenvectors(i,q)
        eigenvectors(i,p) = c*vip-s*viq
        eigenvectors(i,q) = s*vip+c*viq
      end do
    end do
    do i = 1, n
      eigenvalues(i) = a(i,i)
    end do
    do i = 1, n-1
      p = i
      do j = i+1, n
        if (eigenvalues(j) > eigenvalues(p)) p = j
      end do
      if (p /= i) then
        temp = eigenvalues(i)
        eigenvalues(i) = eigenvalues(p)
        eigenvalues(p) = temp
        do j = 1, n
          temp = eigenvectors(j,i)
          eigenvectors(j,i) = eigenvectors(j,p)
          eigenvectors(j,p) = temp
        end do
      end if
    end do
  end subroutine symmetric_eigen_jacobi

  subroutine pca_prelim(covariance, eigenvalues, eigenvectors, explained_fraction)
    real(dp), intent(in) :: covariance(:,:)
    real(dp), intent(out) :: eigenvalues(:), eigenvectors(:,:), explained_fraction(:)
    real(dp) :: total
    call symmetric_eigen_jacobi(covariance,eigenvalues,eigenvectors)
    total = sum(max(eigenvalues,0.0_dp))
    if (total <= 0.0_dp) then
      explained_fraction = 0.0_dp
    else
      explained_fraction = max(eigenvalues,0.0_dp)/total
    end if
  end subroutine pca_prelim

  real(dp) function pca_variance(covariance, positions, number_factors) result(value)
    real(dp), intent(in) :: covariance(:,:), positions(:)
    integer, intent(in) :: number_factors
    real(dp), allocatable :: eigenvalues(:), eigenvectors(:,:), exposures(:)
    integer :: n, k
    n = size(positions)
    allocate(eigenvalues(n),eigenvectors(n,n),exposures(n))
    call symmetric_eigen_jacobi(covariance,eigenvalues,eigenvectors)
    exposures = matmul(transpose(eigenvectors),positions)
    k = min(max(1,number_factors),n)
    value = sum(max(eigenvalues(1:k),0.0_dp)*exposures(1:k)**2)
  end function pca_variance

  real(dp) function pca_var(covariance, mu, positions, number_factors, cl, hp) result(value)
    real(dp), intent(in) :: covariance(:,:), mu(:), positions(:), cl, hp
    integer, intent(in) :: number_factors
    real(dp) :: mean_p, sd_p
    call validate_dimensions(covariance,mu,positions)
    mean_p = portfolio_mean(mu,positions)
    sd_p = sqrt(max(0.0_dp,pca_variance(covariance,positions,number_factors)))
    value = -mean_p*hp-normal_quantile(1.0_dp-cl)*sd_p*sqrt(hp)
  end function pca_var

  real(dp) function pca_es(covariance, mu, positions, number_factors, cl, hp) result(value)
    real(dp), intent(in) :: covariance(:,:), mu(:), positions(:), cl, hp
    integer, intent(in) :: number_factors
    real(dp) :: mean_p, sd_p, z
    call validate_dimensions(covariance,mu,positions)
    mean_p = portfolio_mean(mu,positions)
    sd_p = sqrt(max(0.0_dp,pca_variance(covariance,positions,number_factors)))
    z = normal_quantile(cl)
    value = sd_p*sqrt(hp)*normal_pdf(z)/(1.0_dp-cl)-mean_p*hp
  end function pca_es

end module dowd_portfolio
