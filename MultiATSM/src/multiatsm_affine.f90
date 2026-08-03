! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module multiatsm_affine
  use multiatsm_kinds, only : dp
  use multiatsm_linalg, only : inverse_matrix, pseudo_inverse, eye
  use multiatsm_types, only : affine_loadings
  implicit none
  private

  public :: affine_yield_loadings, multicountry_affine_loadings
  public :: pricing_factor_loadings, rotate_latent_to_observed
  public :: estimate_long_run_short_rate, build_yield_intercepts

contains

  subroutine affine_yield_loadings(maturities, k1q, r0, sigma, short_rate_loadings, loadings, info)
    integer, intent(in) :: maturities(:)
    real(dp), intent(in) :: k1q(:, :), r0, sigma(:, :), short_rate_loadings(:)
    type(affine_loadings), intent(out) :: loadings
    integer, intent(out) :: info
    real(dp), allocatable :: b_all(:, :), a_all(:)
    integer :: n, maxmat, h, j

    n = size(k1q, 1)
    if (size(k1q, 2) /= n .or. size(sigma, 1) /= n .or. size(sigma, 2) /= n .or. &
        size(short_rate_loadings) /= n .or. size(maturities) < 1 .or. any(maturities < 1)) then
      info = -1
      return
    end if
    maxmat = maxval(maturities)
    allocate(b_all(0:maxmat, n), a_all(0:maxmat))
    b_all = 0.0_dp
    a_all = 0.0_dp
    do h = 1, maxmat
      b_all(h, :) = short_rate_loadings + matmul(b_all(h - 1, :), k1q)
      a_all(h) = r0 + a_all(h - 1) - 0.5_dp * dot_product(b_all(h - 1, :), &
        matmul(sigma, b_all(h - 1, :)))
    end do
    allocate(loadings%a(size(maturities)), loadings%b(size(maturities), n), &
      loadings%b_adjustment(size(maturities)))
    do j = 1, size(maturities)
      loadings%a(j) = a_all(maturities(j)) / real(maturities(j), dp)
      loadings%b(j, :) = b_all(maturities(j), :) / real(maturities(j), dp)
      loadings%b_adjustment(j) = loadings%a(j) - r0
    end do
    info = 0
  end subroutine affine_yield_loadings

  subroutine multicountry_affine_loadings(maturities, k1q, r0, sigma, n_countries, loadings, info)
    integer, intent(in) :: maturities(:)
    real(dp), intent(in) :: k1q(:, :), r0(:), sigma(:, :)
    integer, intent(in) :: n_countries
    type(affine_loadings), intent(out) :: loadings
    integer, intent(out) :: info
    real(dp), allocatable :: d(:), b_all(:, :), a_all(:)
    integer :: nc, n_total, n, j, c, maxmat, h, row

    nc = n_countries
    n_total = size(k1q, 1)
    if (nc < 1 .or. mod(n_total, nc) /= 0 .or. size(k1q, 2) /= n_total .or. &
        size(sigma, 1) /= n_total .or. size(sigma, 2) /= n_total .or. size(r0) /= nc .or. &
        size(maturities) < 1 .or. any(maturities < 1)) then
      info = -1
      return
    end if
    n = n_total / nc
    maxmat = maxval(maturities)
    allocate(loadings%a(nc * size(maturities)), loadings%b(nc * size(maturities), n_total), &
      loadings%b_adjustment(nc * size(maturities)))
    allocate(d(n_total), b_all(0:maxmat, n_total), a_all(0:maxmat))

    do c = 1, nc
      d = 0.0_dp
      d((c - 1) * n + 1:c * n) = 1.0_dp
      b_all = 0.0_dp
      a_all = 0.0_dp
      do h = 1, maxmat
        b_all(h, :) = d + matmul(b_all(h - 1, :), k1q)
        a_all(h) = r0(c) + a_all(h - 1) - 0.5_dp * dot_product(b_all(h - 1, :), &
          matmul(sigma, b_all(h - 1, :)))
      end do
      do j = 1, size(maturities)
        row = (c - 1) * size(maturities) + j
        loadings%a(row) = a_all(maturities(j)) / real(maturities(j), dp)
        loadings%b(row, :) = b_all(maturities(j), :) / real(maturities(j), dp)
        loadings%b_adjustment(row) = loadings%a(row) - r0(c)
      end do
    end do
    info = 0
  end subroutine multicountry_affine_loadings

  subroutine pricing_factor_loadings(latent_b, wpca, observed_b, rotation, info)
    real(dp), intent(in) :: latent_b(:, :), wpca(:, :)
    real(dp), allocatable, intent(out) :: observed_b(:, :), rotation(:, :)
    integer, intent(out) :: info
    real(dp), allocatable :: invrot(:, :)

    if (size(wpca, 2) /= size(latent_b, 1)) then
      info = -1
      allocate(observed_b(0, 0), rotation(0, 0))
      return
    end if
    allocate(rotation(size(wpca, 1), size(latent_b, 2)))
    rotation = matmul(wpca, latent_b)
    if (size(rotation, 1) == size(rotation, 2)) then
      call inverse_matrix(rotation, invrot, info, 1.0e-12_dp)
    else
      call pseudo_inverse(rotation, invrot, info)
    end if
    if (info /= 0) then
      allocate(observed_b(0, 0))
      return
    end if
    allocate(observed_b(size(latent_b, 1), size(invrot, 2)))
    observed_b = matmul(latent_b, invrot)
  end subroutine pricing_factor_loadings

  subroutine build_yield_intercepts(latent_a, observed_b, wpca, observed_a, info)
    real(dp), intent(in) :: latent_a(:), observed_b(:, :), wpca(:, :)
    real(dp), allocatable, intent(out) :: observed_a(:)
    integer, intent(out) :: info

    if (size(observed_b, 1) /= size(latent_a) .or. size(wpca, 2) /= size(latent_a) .or. &
        size(observed_b, 2) /= size(wpca, 1)) then
      info = -1
      allocate(observed_a(0))
      return
    end if
    allocate(observed_a(size(latent_a)))
    observed_a = latent_a - matmul(observed_b, matmul(wpca, latent_a))
    info = 0
  end subroutine build_yield_intercepts

  subroutine rotate_latent_to_observed(latent, u1, u0, observed, observed_sigma, latent_sigma, info)
    type(affine_loadings), intent(in) :: latent
    real(dp), intent(in) :: u1(:, :), u0(:), latent_sigma(:, :)
    type(affine_loadings), intent(out) :: observed
    real(dp), allocatable, intent(out) :: observed_sigma(:, :)
    integer, intent(out) :: info
    real(dp), allocatable :: invu1(:, :)

    if (size(u1, 1) /= size(u1, 2) .or. size(u0) /= size(u1, 1) .or. &
        size(latent%b, 2) /= size(u1, 1) .or. size(latent_sigma, 1) /= size(u1, 1) .or. &
        size(latent_sigma, 2) /= size(u1, 1)) then
      info = -1
      allocate(observed_sigma(0, 0))
      return
    end if
    call inverse_matrix(u1, invu1, info, 1.0e-12_dp)
    if (info /= 0) return
    allocate(observed%b(size(latent%b, 1), size(latent%b, 2)), observed%a(size(latent%a)), &
      observed%b_adjustment(size(latent%b_adjustment)), observed_sigma(size(u1, 1), size(u1, 1)))
    observed%b = matmul(latent%b, invu1)
    observed%a = latent%a - matmul(observed%b, u0)
    observed%b_adjustment = observed%a
    observed_sigma = matmul(u1, matmul(latent_sigma, transpose(u1)))
    info = 0
  end subroutine rotate_latent_to_observed

  subroutine estimate_long_run_short_rate(yields, factors, observed_b, b_adjustment, wpca, we, &
      n_countries, r0, info)
    real(dp), intent(in) :: yields(:, :), factors(:, :), observed_b(:, :), b_adjustment(:)
    real(dp), intent(in) :: wpca(:, :), we(:, :)
    integer, intent(in) :: n_countries
    real(dp), allocatable, intent(out) :: r0(:)
    integer, intent(out) :: info
    real(dp), allocatable :: a0(:), a1(:), aper(:), mean_resid(:), ones(:)
    integer :: nc, j, nf, ne, t, c, y0, y1, f0, f1, e0, e1
    real(dp) :: num, den

    nc = n_countries
    if (nc < 1 .or. mod(size(yields, 1), nc) /= 0 .or. mod(size(factors, 1), nc) /= 0 .or. &
        size(yields, 2) /= size(factors, 2) .or. size(observed_b, 1) /= size(yields, 1) .or. &
        size(observed_b, 2) /= size(factors, 1) .or. size(b_adjustment) /= size(yields, 1) .or. &
        size(wpca, 2) /= size(yields, 1) .or. size(we, 2) /= size(yields, 1)) then
      info = -1
      allocate(r0(0))
      return
    end if
    j = size(yields, 1) / nc
    nf = size(factors, 1) / nc
    ne = size(we, 1) / nc
    t = size(yields, 2)
    allocate(r0(nc))
    do c = 1, nc
      y0 = (c - 1) * j + 1
      y1 = c * j
      f0 = (c - 1) * nf + 1
      f1 = c * nf
      e0 = (c - 1) * ne + 1
      e1 = c * ne
      allocate(a0(j), a1(j), aper(j), mean_resid(j), ones(j))
      ones = 1.0_dp
      a0 = b_adjustment(y0:y1) - matmul(observed_b(y0:y1, f0:f1), &
        matmul(wpca(f0:f1, y0:y1), b_adjustment(y0:y1)))
      a1 = 1.0_dp - matmul(observed_b(y0:y1, f0:f1), &
        matmul(wpca(f0:f1, y0:y1), ones))
      aper = matmul(transpose(we(e0:e1, y0:y1)), matmul(we(e0:e1, y0:y1), a1))
      mean_resid = sum(yields(y0:y1, 2:t) - matmul(observed_b(y0:y1, f0:f1), &
        factors(f0:f1, 2:t)), dim=2) / real(t - 1, dp) - a0
      num = dot_product(aper, a1)
      den = dot_product(aper, mean_resid)
      if (abs(num) <= tiny(1.0_dp)) then
        info = -2
        return
      end if
      r0(c) = den / num
      deallocate(a0, a1, aper, mean_resid, ones)
    end do
    info = 0
  end subroutine estimate_long_run_short_rate

end module multiatsm_affine
