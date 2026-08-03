! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module multiatsm_jll
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use multiatsm_kinds, only : dp
  use multiatsm_linalg, only : least_squares, inverse_matrix, cholesky_lower, eye
  use multiatsm_var, only : fit_var
  use multiatsm_types, only : var_model, jll_model
  implicit none
  private

  public :: fit_jll, jll_feedback_restrictions, jll_cholesky_mask

contains

  subroutine regression_no_intercept(y, x, beta, residuals, info)
    real(dp), intent(in) :: y(:, :), x(:, :)
    real(dp), allocatable, intent(out) :: beta(:, :), residuals(:, :)
    integer, intent(out) :: info
    real(dp), allocatable :: sol(:, :)
    integer :: ny, nx, t

    ny = size(y, 1)
    nx = size(x, 1)
    t = size(y, 2)
    if (size(x, 2) /= t) then
      info = -1
      allocate(beta(0, 0), residuals(0, 0))
      return
    end if
    if (nx == 0) then
      allocate(beta(ny, 0), residuals(ny, t))
      residuals = y
      info = 0
      return
    end if
    call least_squares(transpose(x), transpose(y), sol, info)
    if (info /= 0) then
      allocate(beta(0, 0), residuals(0, 0))
      return
    end if
    allocate(beta(ny, nx), residuals(ny, t))
    beta = transpose(sol)
    residuals = y - matmul(beta, x)
  end subroutine regression_no_intercept

  subroutine jll_feedback_restrictions(g, m, n, c, dominant_index, restrictions, free_mask, info)
    integer, intent(in) :: g, m, n, c, dominant_index
    real(dp), allocatable, intent(out) :: restrictions(:, :)
    logical, allocatable, intent(out) :: free_mask(:, :)
    integer, intent(out) :: info
    integer :: k, i, r0, r1, q0, q1, block
    real(dp) :: nan_value

    k = g + c * (m + n)
    if (g < 0 .or. m < 0 .or. n < 1 .or. c < 1 .or. dominant_index < 0 .or. dominant_index > c) then
      info = -1
      allocate(restrictions(0, 0), free_mask(0, 0))
      return
    end if
    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    allocate(restrictions(k, k + 1), free_mask(k, k + 1))
    restrictions = 0.0_dp
    free_mask = .false.
    free_mask(:, 1) = .true.
    if (g > 0) free_mask(:, 2:g+1) = .true.

    block = m + n
    if (dominant_index == 0) then
      do i = 1, c
        r0 = g + (i - 1) * block + 1
        r1 = g + i * block
        q0 = 1 + r0
        q1 = 1 + r1
        free_mask(r0:r1, q0:q1) = .true.
      end do
    else
      ! fit_jll reorders the dominant country first
      q0 = g + 2
      q1 = g + block + 1
      free_mask(:, q0:q1) = .true.
      do i = 2, c
        r0 = g + (i - 1) * block + 1
        r1 = g + i * block
        q0 = 1 + r0
        q1 = 1 + r1
        free_mask(r0:r1, q0:q1) = .true.
      end do
    end if
    where (free_mask) restrictions = nan_value
    info = 0
  end subroutine jll_feedback_restrictions

  subroutine jll_cholesky_mask(g, m, n, c, dominant_index, free_mask, info)
    integer, intent(in) :: g, m, n, c, dominant_index
    logical, allocatable, intent(out) :: free_mask(:, :)
    integer, intent(out) :: info
    integer :: k, i, j, block, ri0, rmi1, rpi0, rpi1
    integer :: cj0, cmj1, cpj0, cpj1

    k = g + c * (m + n)
    if (g < 0 .or. m < 0 .or. n < 1 .or. c < 1 .or. dominant_index < 0 .or. dominant_index > c) then
      info = -1
      allocate(free_mask(0, 0))
      return
    end if
    allocate(free_mask(k, k))
    free_mask = .false.
    do i = 1, k
      free_mask(i, 1:i) = .true.
    end do
    block = m + n

    ! Global shocks do not load directly on pricing-factor rows.
    if (g > 0) then
      do i = 1, c
        rpi0 = g + (i - 1) * block + m + 1
        rpi1 = rpi0 + n - 1
        free_mask(rpi0:rpi1, 1:g) = .false.
      end do
    end if

    ! Macro/pricing cross-block restrictions.
    do i = 1, c
      ri0 = g + (i - 1) * block + 1
      rmi1 = ri0 + m - 1
      rpi0 = rmi1 + 1
      rpi1 = rpi0 + n - 1
      do j = 1, c
        cj0 = g + (j - 1) * block + 1
        cmj1 = cj0 + m - 1
        cpj0 = cmj1 + 1
        cpj1 = cpj0 + n - 1
        if (m > 0) then
          free_mask(rpi0:rpi1, cj0:cmj1) = .false.
          free_mask(ri0:rmi1, cpj0:cpj1) = .false.
        end if
        if (dominant_index > 0 .and. i /= j .and. i > 1 .and. j > 1) then
          if (m > 0) free_mask(ri0:rmi1, cj0:cmj1) = .false.
          free_mask(rpi0:rpi1, cpj0:cpj1) = .false.
        end if
      end do
    end do
    info = 0
  end subroutine jll_cholesky_mask

  subroutine fit_jll(nonorthogonal_factors, g, m, n, c, dominant_index, model, info)
    real(dp), intent(in) :: nonorthogonal_factors(:, :)
    integer, intent(in) :: g, m, n, c, dominant_index
    type(jll_model), intent(out) :: model
    integer, intent(out) :: info
    integer :: k, t, block, i, orig, r0, r1, macro1, price0, price1
    integer :: dom_macro0, dom_macro1, dom_price0, dom_price1
    integer, allocatable :: order(:)
    real(dp), allocatable :: y(:, :), global(:, :), macro(:, :, :), pricing(:, :, :)
    real(dp), allocatable :: b(:, :, :), aw(:, :, :), adu(:, :, :), cc(:, :, :)
    real(dp), allocatable :: me(:, :, :), pe(:, :, :), me_final(:, :, :), pe_final(:, :, :)
    real(dp), allocatable :: beta(:, :), resid(:, :), xreg(:, :)
    real(dp), allocatable :: restrictions(:, :), invpi(:, :), l(:, :)
    type(var_model) :: vm
    logical, allocatable :: fmask(:, :), cmask(:, :)

    k = g + c * (m + n)
    t = size(nonorthogonal_factors, 2)
    block = m + n
    if (size(nonorthogonal_factors, 1) /= k .or. t < 2 .or. n < 1 .or. c < 1 .or. &
        dominant_index < 0 .or. dominant_index > c) then
      info = -1
      return
    end if

    allocate(order(c))
    if (dominant_index > 0) then
      order(1) = dominant_index
      r0 = 1
      do i = 1, c
        if (i /= dominant_index) then
          r0 = r0 + 1
          order(r0) = i
        end if
      end do
    else
      order = [(i, i = 1, c)]
    end if

    allocate(y(k, t), global(g, t), macro(c, m, t), pricing(c, n, t))
    if (g > 0) then
      global = nonorthogonal_factors(1:g, :)
      y(1:g, :) = global
    end if
    do i = 1, c
      orig = order(i)
      r0 = g + (orig - 1) * block + 1
      if (m > 0) macro(i, :, :) = nonorthogonal_factors(r0:r0+m-1, :)
      pricing(i, :, :) = nonorthogonal_factors(r0+m:r0+block-1, :)
      r1 = g + (i - 1) * block + 1
      if (m > 0) y(r1:r1+m-1, :) = macro(i, :, :)
      y(r1+m:r1+block-1, :) = pricing(i, :, :)
    end do

    allocate(b(c, n, m), aw(c, m, g), adu(c, m, m), cc(c, n, n))
    allocate(me(c, m, t), pe(c, n, t), me_final(c, m, t), pe_final(c, n, t))
    b = 0.0_dp
    aw = 0.0_dp
    adu = 0.0_dp
    cc = 0.0_dp
    do i = 1, c
      call regression_no_intercept(pricing(i, :, :), macro(i, :, :), beta, resid, info)
      if (info /= 0) return
      if (m > 0) b(i, :, :) = beta
      pe(i, :, :) = resid
      if (m > 0) then
        call regression_no_intercept(macro(i, :, :), global, beta, resid, info)
        if (info /= 0) return
        if (g > 0) aw(i, :, :) = beta
        me(i, :, :) = resid
      end if
    end do
    me_final = me
    pe_final = pe

    if (dominant_index > 0) then
      do i = 2, c
        if (m > 0) then
          allocate(xreg(g + m, t))
          if (g > 0) xreg(1:g, :) = global
          xreg(g+1:g+m, :) = me(1, :, :)
          call regression_no_intercept(macro(i, :, :), xreg, beta, resid, info)
          if (info /= 0) return
          if (g > 0) aw(i, :, :) = beta(:, 1:g)
          adu(i, :, :) = beta(:, g+1:g+m)
          me_final(i, :, :) = resid
          deallocate(xreg)
        end if
        call regression_no_intercept(pe(i, :, :), pe(1, :, :), beta, resid, info)
        if (info /= 0) return
        cc(i, :, :) = beta
        pe_final(i, :, :) = resid
      end do
    end if

    allocate(model%pi_matrix(k, k), model%orthogonal_factors(k, t))
    model%pi_matrix = eye(k)
    model%orthogonal_factors = 0.0_dp
    if (g > 0) model%orthogonal_factors(1:g, :) = global

    ! PIb contributions: pricing = b * macro + pricing residual.
    do i = 1, c
      r0 = g + (i - 1) * block + 1
      macro1 = r0 + m - 1
      price0 = macro1 + 1
      price1 = price0 + n - 1
      if (m > 0) model%pi_matrix(price0:price1, r0:macro1) = b(i, :, :)
      if (m > 0) model%orthogonal_factors(r0:macro1, :) = me_final(i, :, :)
      model%orthogonal_factors(price0:price1, :) = pe_final(i, :, :)
    end do

    ! PIac contributions are applied through a separate matrix.
    allocate(xreg(k, k))
    xreg = eye(k)
    do i = 1, c
      r0 = g + (i - 1) * block + 1
      macro1 = r0 + m - 1
      price0 = macro1 + 1
      price1 = price0 + n - 1
      if (m > 0 .and. g > 0) xreg(r0:macro1, 1:g) = aw(i, :, :)
      if (dominant_index > 0 .and. i > 1) then
        dom_macro0 = g + 1
        dom_macro1 = g + m
        dom_price0 = dom_macro1 + 1
        dom_price1 = dom_price0 + n - 1
        if (m > 0) xreg(r0:macro1, dom_macro0:dom_macro1) = adu(i, :, :)
        xreg(price0:price1, dom_price0:dom_price1) = cc(i, :, :)
      end if
    end do
    model%pi_matrix = matmul(model%pi_matrix, xreg)

    call jll_feedback_restrictions(g, m, n, c, dominant_index, restrictions, fmask, info)
    if (info /= 0) return
    call fit_var(model%orthogonal_factors, vm, info, restrictions)
    if (info /= 0) return
    allocate(model%k0_e(k), model%k1_e(k, k), model%k0(k), model%k1(k, k))
    model%k0_e = vm%intercept
    model%k1_e = vm%phi
    call inverse_matrix(model%pi_matrix, invpi, info, 1.0e-12_dp)
    if (info /= 0) return
    model%k0 = matmul(model%pi_matrix, model%k0_e)
    model%k1 = matmul(model%pi_matrix, matmul(model%k1_e, invpi))
    allocate(model%feedback_free(k, k + 1))
    model%feedback_free = fmask

    call jll_cholesky_mask(g, m, n, c, dominant_index, cmask, info)
    if (info /= 0) return
    call cholesky_lower(vm%sigma, l, info, 1.0e-12_dp)
    if (info /= 0) return
    where (.not. cmask) l = 0.0_dp
    allocate(model%chol_ortho(k, k), model%chol_nonortho(k, k))
    allocate(model%sigma_ortho(k, k), model%sigma_nonortho(k, k), model%chol_free(k, k))
    model%chol_ortho = l
    model%sigma_ortho = matmul(l, transpose(l))
    model%chol_nonortho = matmul(model%pi_matrix, l)
    model%sigma_nonortho = matmul(model%chol_nonortho, transpose(model%chol_nonortho))
    model%chol_free = cmask
    info = 0
  end subroutine fit_jll

end module multiatsm_jll
