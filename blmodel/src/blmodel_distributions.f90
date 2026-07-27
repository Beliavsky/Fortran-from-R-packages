! SPDX-License-Identifier: GPL-3.0-only
! Derived from BLModel 1.0.2, Copyright (C) 2017 Andrzej Palczewski and Jan Palczewski.
module blmodel_distributions
  use blmodel_kinds, only : dp, pi
  use blmodel_linalg, only : cholesky_lower, quadratic_form_cholesky, logdet_cholesky
  implicit none
  private

  abstract interface
    subroutine view_density_interface(points, q, covmat, params, density, info)
      import dp
      real(dp), intent(in) :: points(:,:), q(:), covmat(:,:)
      real(dp), intent(in), optional :: params(:)
      real(dp), allocatable, intent(out) :: density(:)
      integer, intent(out) :: info
    end subroutine view_density_interface
  end interface

  public :: view_density_interface
  public :: observ_normal, observ_powerexp, observ_student_t, observ_ts

contains

  subroutine observ_normal(points, q, covmat, params, density, info)
    real(dp), intent(in) :: points(:,:), q(:), covmat(:,:)
    real(dp), intent(in), optional :: params(:)
    real(dp), allocatable, intent(out) :: density(:)
    integer, intent(out) :: info
    real(dp), allocatable :: lower(:,:)
    real(dp) :: log_constant, quad
    integer :: j, k

    if (present(params)) continue
    call validate_density_inputs(points, q, covmat, info)
    if (info /= 0) then
      allocate(density(0))
      return
    end if

    call cholesky_lower(covmat, lower, info)
    if (info /= 0) then
      allocate(density(0))
      info = 10 + info
      return
    end if

    k = size(q)
    log_constant = -0.5_dp * (real(k, dp) * log(2.0_dp * pi) + logdet_cholesky(lower))
    allocate(density(size(points, 2)))
    do j = 1, size(points, 2)
      quad = quadratic_form_cholesky(lower, points(:, j) - q)
      density(j) = exp(log_constant - 0.5_dp * quad)
    end do
  end subroutine observ_normal

  subroutine observ_student_t(points, q, covmat, params, density, info)
    real(dp), intent(in) :: points(:,:), q(:), covmat(:,:)
    real(dp), intent(in), optional :: params(:)
    real(dp), allocatable, intent(out) :: density(:)
    integer, intent(out) :: info
    real(dp), allocatable :: omega(:,:), lower(:,:)
    real(dp) :: df, log_constant, quad
    integer :: j, k

    df = 5.0_dp
    if (present(params)) then
      if (size(params) >= 1) df = params(1)
    end if
    call validate_density_inputs(points, q, covmat, info)
    if (info /= 0 .or. df <= 2.0_dp) then
      allocate(density(0))
      if (info == 0) info = 4
      return
    end if

    k = size(q)
    allocate(omega(k, k))
    omega = ((df - 2.0_dp) / df) * covmat
    call cholesky_lower(omega, lower, info)
    if (info /= 0) then
      allocate(density(0))
      info = 10 + info
      return
    end if

    log_constant = log_gamma(0.5_dp * (df + real(k, dp))) - log_gamma(0.5_dp * df) &
      - 0.5_dp * (real(k, dp) * log(df * pi) + logdet_cholesky(lower))
    allocate(density(size(points, 2)))
    do j = 1, size(points, 2)
      quad = quadratic_form_cholesky(lower, points(:, j) - q)
      density(j) = exp(log_constant - 0.5_dp * (df + real(k, dp)) * log(1.0_dp + quad / df))
    end do
  end subroutine observ_student_t

  subroutine observ_ts(points, q, covmat, params, density, info)
    real(dp), intent(in) :: points(:,:), q(:), covmat(:,:)
    real(dp), intent(in), optional :: params(:)
    real(dp), allocatable, intent(out) :: density(:)
    integer, intent(out) :: info

    if (present(params)) then
      call observ_student_t(points, q, covmat, params, density, info)
    else
      call observ_student_t(points, q, covmat, density=density, info=info)
    end if
  end subroutine observ_ts

  subroutine observ_powerexp(points, q, covmat, params, density, info)
    real(dp), intent(in) :: points(:,:), q(:), covmat(:,:)
    real(dp), intent(in), optional :: params(:)
    real(dp), allocatable, intent(out) :: density(:)
    integer, intent(out) :: info
    real(dp), allocatable :: omega(:,:), lower(:,:)
    real(dp) :: beta, scale, log_constant, quad, rk
    integer :: j, k

    beta = 0.6_dp
    if (present(params)) then
      if (size(params) >= 1) beta = params(1)
    end if
    call validate_density_inputs(points, q, covmat, info)
    if (info /= 0 .or. beta <= 0.0_dp) then
      allocate(density(0))
      if (info == 0) info = 4
      return
    end if

    k = size(q)
    rk = real(k, dp)
    scale = exp(log(rk) + log_gamma(rk / (2.0_dp * beta)) - log_gamma((rk + 2.0_dp) / (2.0_dp * beta)) &
      - log(2.0_dp) / beta)
    allocate(omega(k, k))
    omega = scale * covmat
    call cholesky_lower(omega, lower, info)
    if (info /= 0) then
      allocate(density(0))
      info = 10 + info
      return
    end if

    log_constant = log(rk) + log_gamma(0.5_dp * rk) - 0.5_dp * rk * log(pi) &
      - log_gamma(1.0_dp + rk / (2.0_dp * beta)) &
      - (1.0_dp + rk / (2.0_dp * beta)) * log(2.0_dp) - 0.5_dp * logdet_cholesky(lower)
    allocate(density(size(points, 2)))
    do j = 1, size(points, 2)
      quad = quadratic_form_cholesky(lower, points(:, j) - q)
      density(j) = exp(log_constant - 0.5_dp * quad**beta)
    end do
  end subroutine observ_powerexp

  subroutine validate_density_inputs(points, q, covmat, info)
    real(dp), intent(in) :: points(:,:), q(:), covmat(:,:)
    integer, intent(out) :: info

    info = 0
    if (size(points, 1) /= size(q)) info = 1
    if (size(covmat, 1) /= size(q) .or. size(covmat, 2) /= size(q)) info = 2
    if (size(points, 2) == 0 .or. size(q) == 0) info = 3
  end subroutine validate_density_inputs

end module blmodel_distributions
