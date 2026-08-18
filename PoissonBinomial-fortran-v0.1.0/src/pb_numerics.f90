! SPDX-License-Identifier: GPL-3.0-only
module pb_numerics
  use pb_kinds, only : dp, pi
  implicit none
  private
  public :: normalize_pmf, convolve_naive, convolve_fft, gcd_vector
  public :: sample_from_pmf, expand_probs, expand_gpb

contains

  subroutine normalize_pmf(p)
    real(dp), intent(inout) :: p(:)
    real(dp) :: s
    where (p < 0.0_dp .and. p > -1.0e-13_dp) p = 0.0_dp
    where (p > 1.0_dp .and. p < 1.0_dp + 1.0e-13_dp) p = 1.0_dp
    s = sum(p)
    if (s > 0.0_dp) p = p/s
  end subroutine normalize_pmf

  function convolve_naive(a, b) result(c)
    real(dp), intent(in) :: a(:), b(:)
    real(dp), allocatable :: c(:)
    integer :: i, j, la, lb
    la = size(a)
    lb = size(b)
    allocate(c(0:la+lb-2))
    c = 0.0_dp
    do i = 1, la
      if (abs(a(i)) <= tiny(1.0_dp)) cycle
      do j = 1, lb
        if (abs(b(j)) <= tiny(1.0_dp)) cycle
        c(i+j-2) = c(i+j-2) + a(i)*b(j)
      end do
    end do
  end function convolve_naive

  subroutine fft_inplace(z, inverse)
    complex(dp), intent(inout) :: z(:)
    logical, intent(in) :: inverse
    integer :: n, i, j, m, len, half, k
    complex(dp) :: tmp, w, wlen
    real(dp) :: ang

    n = size(z)
    j = 1
    do i = 2, n
      m = n/2
      do while (j > m .and. m >= 1)
        j = j - m
        m = m/2
      end do
      j = j + m
      if (i < j) then
        tmp = z(i)
        z(i) = z(j)
        z(j) = tmp
      end if
    end do

    len = 2
    do while (len <= n)
      ang = 2.0_dp*pi/real(len,dp)
      if (.not. inverse) ang = -ang
      wlen = cmplx(cos(ang), sin(ang), kind=dp)
      half = len/2
      do i = 1, n, len
        w = cmplx(1.0_dp, 0.0_dp, kind=dp)
        do k = 0, half-1
          tmp = w*z(i+k+half)
          z(i+k+half) = z(i+k) - tmp
          z(i+k) = z(i+k) + tmp
          w = w*wlen
        end do
      end do
      len = 2*len
    end do
    if (inverse) z = z/real(n,dp)
  end subroutine fft_inplace

  function convolve_fft(a, b) result(c)
    real(dp), intent(in) :: a(:), b(:)
    real(dp), allocatable :: c(:)
    complex(dp), allocatable :: za(:), zb(:)
    integer :: need, nfft, i

    need = size(a) + size(b) - 1
    if (need <= 64) then
      c = convolve_naive(a,b)
      return
    end if
    nfft = 1
    do while (nfft < need)
      nfft = 2*nfft
    end do
    allocate(za(nfft), zb(nfft))
    za = cmplx(0.0_dp, 0.0_dp, kind=dp)
    zb = cmplx(0.0_dp, 0.0_dp, kind=dp)
    do i = 1, size(a)
      za(i) = cmplx(a(i),0.0_dp,kind=dp)
    end do
    do i = 1, size(b)
      zb(i) = cmplx(b(i),0.0_dp,kind=dp)
    end do
    call fft_inplace(za,.false.)
    call fft_inplace(zb,.false.)
    za = za*zb
    call fft_inplace(za,.true.)
    allocate(c(0:need-1))
    do i = 0, need-1
      c(i) = real(za(i+1),dp)
    end do
    where (abs(c) < 5.0e-15_dp) c = 0.0_dp
  end function convolve_fft

  pure integer function gcd_pair(a, b) result(g)
    integer, intent(in) :: a, b
    integer :: x, y, r
    x = abs(a)
    y = abs(b)
    do while (y /= 0)
      r = mod(x,y)
      x = y
      y = r
    end do
    g = x
  end function gcd_pair

  pure integer function gcd_vector(x) result(g)
    integer, intent(in) :: x(:)
    integer :: i
    g = 0
    do i = 1, size(x)
      if (x(i) == 0) cycle
      if (g == 0) then
        g = abs(x(i))
      else
        g = gcd_pair(g,x(i))
      end if
      if (g == 1) exit
    end do
  end function gcd_vector

  integer function sample_from_pmf(pmf, lower) result(x)
    real(dp), intent(in) :: pmf(:)
    integer, intent(in) :: lower
    real(dp) :: u, s
    integer :: i
    call random_number(u)
    s = 0.0_dp
    do i = 1, size(pmf)
      s = s + pmf(i)
      if (u <= s) then
        x = lower + i - 1
        return
      end if
    end do
    x = lower + size(pmf) - 1
  end function sample_from_pmf

  subroutine expand_probs(probs, wts, out)
    real(dp), intent(in) :: probs(:)
    integer, intent(in), optional :: wts(:)
    real(dp), allocatable, intent(out) :: out(:)
    integer :: i, j, k, n
    if (.not. present(wts)) then
      allocate(out(size(probs)))
      out = probs
      return
    end if
    if (size(wts) /= size(probs)) error stop "weights and probabilities differ in length"
    if (any(wts < 0)) error stop "weights must be nonnegative"
    n = sum(wts)
    allocate(out(n))
    k = 0
    do i = 1, size(probs)
      do j = 1, wts(i)
        k = k + 1
        out(k) = probs(i)
      end do
    end do
  end subroutine expand_probs

  subroutine expand_gpb(probs, val_p, val_q, wts, p_out, vp_out, vq_out)
    real(dp), intent(in) :: probs(:)
    integer, intent(in) :: val_p(:), val_q(:)
    integer, intent(in), optional :: wts(:)
    real(dp), allocatable, intent(out) :: p_out(:)
    integer, allocatable, intent(out) :: vp_out(:), vq_out(:)
    integer :: i, j, k, n
    if (size(val_p) /= size(probs) .or. size(val_q) /= size(probs)) then
      error stop "probs, val_p and val_q must have equal length"
    end if
    if (.not. present(wts)) then
      allocate(p_out(size(probs)), vp_out(size(probs)), vq_out(size(probs)))
      p_out = probs
      vp_out = val_p
      vq_out = val_q
      return
    end if
    if (size(wts) /= size(probs)) error stop "weights and probabilities differ in length"
    if (any(wts < 0)) error stop "weights must be nonnegative"
    n = sum(wts)
    allocate(p_out(n), vp_out(n), vq_out(n))
    k = 0
    do i = 1, size(probs)
      do j = 1, wts(i)
        k = k + 1
        p_out(k) = probs(i)
        vp_out(k) = val_p(i)
        vq_out(k) = val_q(i)
      end do
    end do
  end subroutine expand_gpb

end module pb_numerics
