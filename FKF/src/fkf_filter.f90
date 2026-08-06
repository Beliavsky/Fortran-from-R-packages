! SPDX-License-Identifier: GPL-2.0-or-later
module fkf_filter
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
  use fkf_kinds, only : dp
  use fkf_linalg, only : spd_inverse_logdet, symmetrize
  use fkf_types, only : fkf_model, fkf_result, fkf_success, fkf_invalid_input, fkf_non_pos_def
  implicit none
  private
  public :: fkf, kalman_filter, validate_model

  real(dp), parameter :: log_2pi = log(2.0_dp * acos(-1.0_dp))

contains

  subroutine fkf(a0, p0, dt, ct, tt, zt, hht, ggt, yt, result, corrected_missing_likelihood)
    real(dp), intent(in) :: a0(:), p0(:, :), dt(:, :), ct(:, :)
    real(dp), intent(in) :: tt(:, :, :), zt(:, :, :), hht(:, :, :), ggt(:, :, :), yt(:, :)
    type(fkf_result), intent(out) :: result
    logical, intent(in), optional :: corrected_missing_likelihood
    type(fkf_model) :: model

    model%a0 = a0
    model%p0 = p0
    model%dt = dt
    model%ct = ct
    model%tt = tt
    model%zt = zt
    model%hht = hht
    model%ggt = ggt
    call kalman_filter(model, yt, result, corrected_missing_likelihood)
  end subroutine fkf

  subroutine kalman_filter(model, yt, result, corrected_missing_likelihood)
    type(fkf_model), intent(in) :: model
    real(dp), intent(in) :: yt(:, :)
    type(fkf_result), intent(out) :: result
    logical, intent(in), optional :: corrected_missing_likelihood

    real(dp), allocatable :: c(:), dvec(:), tmat(:, :), z(:, :), q(:, :), rmat(:, :)
    real(dp), allocatable :: yobs(:), cobs(:), zobs(:, :), robs(:, :)
    real(dp), allocatable :: v(:), f(:, :), finv(:, :), k(:, :), work(:, :)
    real(dp) :: logdet, mahal, nan_value
    integer, allocatable :: pos(:)
    integer :: d, info, m, n, nobs, t
    logical :: corrected

    corrected = .false.
    if (present(corrected_missing_likelihood)) corrected = corrected_missing_likelihood

    result%status = fkf_success
    result%failure_time = 0
    result%message = 'success'
    call validate_model(model, yt, result%status, result%message)
    if (result%status /= fkf_success) return

    m = size(model%a0)
    d = size(yt, 1)
    n = size(yt, 2)
    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)

    allocate(result%att(m, n), result%at(m, n + 1))
    allocate(result%ptt(m, m, n), result%pt(m, m, n + 1))
    allocate(result%vt(d, n), result%ft(d, d, n), result%ftinv(d, d, n), result%kt(m, d, n))
    result%att = nan_value
    result%at = nan_value
    result%ptt = nan_value
    result%pt = nan_value
    result%vt = nan_value
    result%ft = nan_value
    result%ftinv = nan_value
    result%kt = nan_value

    result%at(:, 1) = model%a0
    result%pt(:, :, 1) = model%p0
    result%log_likelihood = 0.0_dp
    if (.not. corrected) result%log_likelihood = -0.5_dp * real(n * d, dp) * log_2pi

    allocate(c(d), dvec(m), tmat(m, m), z(d, m), q(m, m), rmat(d, d), pos(d))

    do t = 1, n
      call get_vector_slice(model%ct, t, c)
      call get_vector_slice(model%dt, t, dvec)
      call get_matrix_slice(model%tt, t, tmat)
      call get_matrix_slice(model%zt, t, z)
      call get_matrix_slice(model%hht, t, q)
      call get_matrix_slice(model%ggt, t, rmat)

      nobs = 0
      do info = 1, d
        if (ieee_is_finite(yt(info, t))) then
          nobs = nobs + 1
          pos(nobs) = info
        end if
      end do

      if (nobs == 0) then
        result%att(:, t) = result%at(:, t)
        result%ptt(:, :, t) = result%pt(:, :, t)
      else
        allocate(yobs(nobs), cobs(nobs), zobs(nobs, m), robs(nobs, nobs))
        allocate(v(nobs), f(nobs, nobs), finv(nobs, nobs), k(m, nobs), work(m, m))
        yobs = yt(pos(1:nobs), t)
        cobs = c(pos(1:nobs))
        zobs = z(pos(1:nobs), :)
        robs = rmat(pos(1:nobs), pos(1:nobs))

        v = yobs - cobs - matmul(zobs, result%at(:, t))
        f = matmul(matmul(zobs, result%pt(:, :, t)), transpose(zobs)) + robs
        call symmetrize(f)
        call spd_inverse_logdet(f, finv, logdet, info)
        if (info /= 0) then
          result%status = fkf_non_pos_def
          result%failure_time = t
          result%message = 'prediction-error covariance is not positive definite'
          return
        end if

        k = matmul(matmul(result%pt(:, :, t), transpose(zobs)), finv)
        result%att(:, t) = result%at(:, t) + matmul(k, v)
        work = matmul(matmul(result%pt(:, :, t), transpose(zobs)), transpose(k))
        result%ptt(:, :, t) = result%pt(:, :, t) - work
        call symmetrize(result%ptt(:, :, t))

        result%vt(pos(1:nobs), t) = v
        result%ft(pos(1:nobs), pos(1:nobs), t) = f
        result%ftinv(pos(1:nobs), pos(1:nobs), t) = finv
        result%kt(:, pos(1:nobs), t) = k

        mahal = dot_product(v, matmul(finv, v))
        if (corrected) then
          result%log_likelihood = result%log_likelihood - &
            0.5_dp * (real(nobs, dp) * log_2pi + logdet + mahal)
        else
          result%log_likelihood = result%log_likelihood - 0.5_dp * (logdet + mahal)
        end if
        deallocate(yobs, cobs, zobs, robs, v, f, finv, k, work)
      end if

      result%at(:, t + 1) = dvec + matmul(tmat, result%att(:, t))
      result%pt(:, :, t + 1) = matmul(matmul(tmat, result%ptt(:, :, t)), transpose(tmat)) + q
      call symmetrize(result%pt(:, :, t + 1))
    end do
  end subroutine kalman_filter

  subroutine validate_model(model, yt, status, message)
    type(fkf_model), intent(in) :: model
    real(dp), intent(in) :: yt(:, :)
    integer, intent(out) :: status
    character(len=:), allocatable, intent(out) :: message
    integer :: d, m, n

    status = fkf_invalid_input
    message = 'invalid model'
    if (.not. allocated(model%a0) .or. .not. allocated(model%p0) .or. &
        .not. allocated(model%dt) .or. .not. allocated(model%ct) .or. &
        .not. allocated(model%tt) .or. .not. allocated(model%zt) .or. &
        .not. allocated(model%hht) .or. .not. allocated(model%ggt)) return

    m = size(model%a0)
    d = size(yt, 1)
    n = size(yt, 2)
    if (m < 1 .or. d < 1 .or. n < 1) then
      message = 'all dimensions must be positive'
      return
    end if
    if (size(model%p0, 1) /= m .or. size(model%p0, 2) /= m) then
      message = 'p0 has incompatible dimensions'
      return
    end if
    if (size(model%dt, 1) /= m .or. .not. valid_time_extent(size(model%dt, 2), n)) then
      message = 'dt must have shape (m,1) or (m,n)'
      return
    end if
    if (size(model%ct, 1) /= d .or. .not. valid_time_extent(size(model%ct, 2), n)) then
      message = 'ct must have shape (d,1) or (d,n)'
      return
    end if
    if (.not. valid_cube(model%tt, m, m, n)) then
      message = 'tt must have shape (m,m,1) or (m,m,n)'
      return
    end if
    if (.not. valid_cube(model%zt, d, m, n)) then
      message = 'zt must have shape (d,m,1) or (d,m,n)'
      return
    end if
    if (.not. valid_cube(model%hht, m, m, n)) then
      message = 'hht must have shape (m,m,1) or (m,m,n)'
      return
    end if
    if (.not. valid_cube(model%ggt, d, d, n)) then
      message = 'ggt must have shape (d,d,1) or (d,d,n)'
      return
    end if
    status = fkf_success
    message = 'success'
  end subroutine validate_model

  pure logical function valid_time_extent(k, n)
    integer, intent(in) :: k, n
    valid_time_extent = k == 1 .or. k == n
  end function valid_time_extent

  pure logical function valid_cube(a, n1, n2, nt)
    real(dp), intent(in) :: a(:, :, :)
    integer, intent(in) :: n1, n2, nt
    valid_cube = size(a, 1) == n1 .and. size(a, 2) == n2 .and. &
      (size(a, 3) == 1 .or. size(a, 3) == nt)
  end function valid_cube

  subroutine get_vector_slice(a, t, out)
    real(dp), intent(in) :: a(:, :)
    integer, intent(in) :: t
    real(dp), intent(out) :: out(:)
    integer :: it
    it = merge(t, 1, size(a, 2) > 1)
    out = a(:, it)
  end subroutine get_vector_slice

  subroutine get_matrix_slice(a, t, out)
    real(dp), intent(in) :: a(:, :, :)
    integer, intent(in) :: t
    real(dp), intent(out) :: out(:, :)
    integer :: it
    it = merge(t, 1, size(a, 3) > 1)
    out = a(:, :, it)
  end subroutine get_matrix_slice

end module fkf_filter
