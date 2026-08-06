! SPDX-License-Identifier: BSD-3-Clause
! Modern Fortran computational translation of waveslim.
module waveslim_dualtree
  use waveslim_kinds, only : dp, sqrt2
  use waveslim_status, only : clear_status, set_status, waveslim_invalid_level
  use waveslim_types, only : complex_wavelet_transform
  implicit none
  private

  type, public :: dual_filter_pair
    real(dp), allocatable :: af1(:,:), af2(:,:)
    real(dp), allocatable :: sf1(:,:), sf2(:,:)
  end type dual_filter_pair

  type, public :: complex_wavelet_level_2d
    complex(dp), allocatable :: lh(:,:), hl(:,:), hh(:,:)
  end type complex_wavelet_level_2d

  type, public :: complex_wavelet_transform_2d
    type(complex_wavelet_level_2d), allocatable :: level(:)
    complex(dp), allocatable :: smooth(:,:)
    integer :: original_shape(2) = 0
    character(len=16) :: method = 'dualtree2d'
  end type complex_wavelet_transform_2d

  public :: farras_filters, fsfarras_filters, dualfilt1_filters
  public :: afb, sfb, dualtree, idualtree
  public :: afb2d, sfb2d, dualtree_2d, idualtree_2d
  public :: circular_shift_dt, circular_shift_2d_dt, plus_minus

  interface plus_minus
    module procedure plus_minus_1d
    module procedure plus_minus_2d
  end interface plus_minus

contains

  function farras_filters() result(f)
    type(dual_filter_pair) :: f
    allocate(f%af1(10,2), f%af2(10,2), f%sf1(10,2), f%sf2(10,2))
    f%af1 = reshape([ &
      0.0_dp, -0.01122679215254_dp, &
      0.0_dp,  0.01122679215254_dp, &
     -0.08838834764832_dp, 0.08838834764832_dp, &
      0.08838834764832_dp, 0.08838834764832_dp, &
      0.69587998903400_dp,-0.69587998903400_dp, &
      0.69587998903400_dp, 0.69587998903400_dp, &
      0.08838834764832_dp,-0.08838834764832_dp, &
     -0.08838834764832_dp,-0.08838834764832_dp, &
      0.01122679215254_dp, 0.0_dp, &
      0.01122679215254_dp, 0.0_dp ], [10,2], order=[2,1])
    ! The single-tree Farras bank is duplicated here for a convenient pair API.
    f%af2 = f%af1
    f%sf1 = f%af1(10:1:-1,:)
    f%sf2 = f%sf1
  end function farras_filters

  function fsfarras_filters() result(f)
    type(dual_filter_pair) :: f
    allocate(f%af1(10,2), f%af2(10,2), f%sf1(10,2), f%sf2(10,2))
    f%af1 = reshape([ &
      0.0_dp,0.0_dp, &
     -0.08838834764832_dp,-0.01122679215254_dp, &
      0.08838834764832_dp, 0.01122679215254_dp, &
      0.69587998903400_dp, 0.08838834764832_dp, &
      0.69587998903400_dp, 0.08838834764832_dp, &
      0.08838834764832_dp,-0.69587998903400_dp, &
     -0.08838834764832_dp, 0.69587998903400_dp, &
      0.01122679215254_dp,-0.08838834764832_dp, &
      0.01122679215254_dp,-0.08838834764832_dp, &
      0.0_dp,0.0_dp ], [10,2], order=[2,1])
    f%af2 = reshape([ &
      0.01122679215254_dp,0.0_dp, &
      0.01122679215254_dp,0.0_dp, &
     -0.08838834764832_dp,-0.08838834764832_dp, &
      0.08838834764832_dp,-0.08838834764832_dp, &
      0.69587998903400_dp, 0.69587998903400_dp, &
      0.69587998903400_dp,-0.69587998903400_dp, &
      0.08838834764832_dp, 0.08838834764832_dp, &
     -0.08838834764832_dp, 0.08838834764832_dp, &
      0.0_dp, 0.01122679215254_dp, &
      0.0_dp,-0.01122679215254_dp ], [10,2], order=[2,1])
    f%sf1 = f%af1(10:1:-1,:)
    f%sf2 = f%af2(10:1:-1,:)
  end function fsfarras_filters

  function dualfilt1_filters() result(f)
    type(dual_filter_pair) :: f
    allocate(f%af1(10,2), f%af2(10,2), f%sf1(10,2), f%sf2(10,2))
    f%af1 = reshape([ &
      0.03516384_dp,0.0_dp, &
      0.0_dp,0.0_dp, &
     -0.08832942_dp,-0.11430184_dp, &
      0.23389032_dp,0.0_dp, &
      0.76027237_dp,0.58751830_dp, &
      0.58751830_dp,-0.76027237_dp, &
      0.0_dp,0.23389032_dp, &
     -0.11430184_dp,0.08832942_dp, &
      0.0_dp,0.0_dp, &
      0.0_dp,-0.03516384_dp ], [10,2], order=[2,1])
    f%af2 = reshape([ &
      0.0_dp,-0.03516384_dp, &
      0.0_dp,0.0_dp, &
     -0.11430184_dp,0.08832942_dp, &
      0.0_dp,0.23389032_dp, &
      0.58751830_dp,-0.76027237_dp, &
      0.76027237_dp,0.58751830_dp, &
      0.23389032_dp,0.0_dp, &
     -0.08832942_dp,-0.11430184_dp, &
      0.0_dp,0.0_dp, &
      0.03516384_dp,0.0_dp ], [10,2], order=[2,1])
    f%sf1 = f%af1(10:1:-1,:)
    f%sf2 = f%af2(10:1:-1,:)
  end function dualfilt1_filters

  function circular_shift_dt(x, m) result(y)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: m
    real(dp), allocatable :: y(:)
    integer :: i, n
    n = size(x)
    allocate(y(n))
    do i = 1, n
      y(i) = x(modulo(i - 1 - m, n) + 1)
    end do
  end function circular_shift_dt

  function circular_shift_2d_dt(x, m, dim) result(y)
    real(dp), intent(in) :: x(:,:)
    integer, intent(in) :: m
    integer, intent(in), optional :: dim
    real(dp), allocatable :: y(:,:)
    integer :: i, j, d
    d = 1
    if (present(dim)) d = dim
    allocate(y(size(x,1), size(x,2)))
    if (d == 1) then
      do j = 1, size(x,2)
        do i = 1, size(x,1)
          y(i,j) = x(modulo(i - 1 - m, size(x,1)) + 1,j)
        end do
      end do
    else
      do j = 1, size(x,2)
        do i = 1, size(x,1)
          y(i,j) = x(i,modulo(j - 1 - m, size(x,2)) + 1)
        end do
      end do
    end if
  end function circular_shift_2d_dt

  subroutine plus_minus_1d(a, b, u, v)
    real(dp), intent(in) :: a(:), b(:)
    real(dp), allocatable, intent(out) :: u(:), v(:)
    allocate(u(size(a)), v(size(a)))
    u = (a + b) / sqrt2
    v = (a - b) / sqrt2
  end subroutine plus_minus_1d

  subroutine plus_minus_2d(a, b, u, v)
    real(dp), intent(in) :: a(:,:), b(:,:)
    real(dp), allocatable, intent(out) :: u(:,:), v(:,:)
    allocate(u(size(a,1),size(a,2)), v(size(a,1),size(a,2)))
    u = (a + b) / sqrt2
    v = (a - b) / sqrt2
  end subroutine plus_minus_2d

  function open_convolve_r(x, h) result(y)
    real(dp), intent(in) :: x(:), h(:)
    real(dp), allocatable :: y(:), z(:)
    integer :: i, j, n, m
    n = size(x)
    m = size(h)
    allocate(z(n+m-1))
    z = 0.0_dp
    do i = 1, n
      do j = 1, m
        z(i+j-1) = z(i+j-1) + x(i)*h(j)
      end do
    end do
    ! R's open convolve ordering is the ordinary convolution shifted by m-1.
    y = circular_shift_dt(z, m-1)
  end function open_convolve_r

  subroutine afb(x, af, lo, hi)
    real(dp), intent(in) :: x(:), af(:,:)
    real(dp), allocatable, intent(out) :: lo(:), hi(:)
    real(dp), allocatable :: work(:), z(:)
    integer :: n, l
    n = size(x)
    l = size(af,1)/2
    work = circular_shift_dt(x, -l)
    z = open_convolve_r(work, af(:,1))
    z = circular_shift_dt(z, -(2*l-1))
    lo = z(1:size(z):2)
    lo(1:l) = lo(1:l) + lo(n/2+1:n/2+l)
    lo = lo(1:n/2)
    z = open_convolve_r(work, af(:,2))
    z = circular_shift_dt(z, -(2*l-1))
    hi = z(1:size(z):2)
    hi(1:l) = hi(1:l) + hi(n/2+1:n/2+l)
    hi = hi(1:n/2)
  end subroutine afb

  function sfb(lo_in, hi_in, sf) result(y)
    real(dp), intent(in) :: lo_in(:), hi_in(:), sf(:,:)
    real(dp), allocatable :: y(:)
    real(dp), allocatable :: ulo(:), uhi(:), lo(:), hi(:)
    integer :: n, l, i
    n = 2*size(lo_in)
    l = size(sf,1)
    allocate(ulo(n), uhi(n))
    ulo = 0.0_dp
    uhi = 0.0_dp
    do i = 1, size(lo_in)
      ulo(2*i) = lo_in(i)
      uhi(2*i) = hi_in(i)
    end do
    lo = open_convolve_r(ulo, sf(:,1))
    hi = open_convolve_r(uhi, sf(:,2))
    lo = circular_shift_dt(lo, -l)
    hi = circular_shift_dt(hi, -l)
    y = lo + hi
    if (l > 2) y(1:l-2) = y(1:l-2) + y(n+1:n+l-2)
    y = y(1:n)
    y = circular_shift_dt(y, 1-l/2)
  end function sfb

  function dualtree(x, n_levels) result(w)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: n_levels
    type(complex_wavelet_transform) :: w
    type(dual_filter_pair) :: first, later
    real(dp), allocatable :: x1(:), x2(:), lo(:), hi(:)
    integer :: j
    w%method = 'dualtree'
    w%wavelet = 'kingsbury'
    w%boundary = 'periodic'
    w%original_length = size(x)
    if (n_levels < 1 .or. mod(size(x),2**n_levels) /= 0) then
      call set_status(w%status, waveslim_invalid_level, &
        'input length must be divisible by 2**levels')
      return
    end if
    first = fsfarras_filters()
    later = dualfilt1_filters()
    allocate(w%detail(n_levels))
    x1 = x/sqrt2
    x2 = x/sqrt2
    call afb(x1, first%af1, lo, hi)
    x1 = lo
    w%detail(1)%values = cmplx(hi, 0.0_dp, dp)
    call afb(x2, first%af2, lo, hi)
    x2 = lo
    w%detail(1)%values = cmplx(real(w%detail(1)%values,dp), hi, dp)
    do j = 2, n_levels
      call afb(x1, later%af1, lo, hi)
      x1 = lo
      w%detail(j)%values = cmplx(hi, 0.0_dp, dp)
      call afb(x2, later%af2, lo, hi)
      x2 = lo
      w%detail(j)%values = cmplx(real(w%detail(j)%values,dp), hi, dp)
    end do
    w%smooth = cmplx(x1, x2, dp)
    call clear_status(w%status)
  end function dualtree

  function idualtree(w) result(x)
    type(complex_wavelet_transform), intent(in) :: w
    real(dp), allocatable :: x(:)
    type(dual_filter_pair) :: first, later
    real(dp), allocatable :: y1(:), y2(:)
    integer :: j
    if (.not. allocated(w%detail) .or. .not. allocated(w%smooth)) then
      allocate(x(0))
      return
    end if
    first = fsfarras_filters()
    later = dualfilt1_filters()
    y1 = real(w%smooth,dp)
    y2 = aimag(w%smooth)
    do j = size(w%detail), 2, -1
      y1 = sfb(y1, real(w%detail(j)%values,dp), later%sf1)
      y2 = sfb(y2, aimag(w%detail(j)%values), later%sf2)
    end do
    y1 = sfb(y1, real(w%detail(1)%values,dp), first%sf1)
    y2 = sfb(y2, aimag(w%detail(1)%values), first%sf2)
    x = (y1+y2)/sqrt2
    if (w%original_length > 0 .and. size(x) > w%original_length) then
      x = x(1:w%original_length)
    end if
  end function idualtree

  subroutine afb2d(x, af1, af2, lo, lh, hl, hh)
    real(dp), intent(in) :: x(:,:), af1(:,:)
    real(dp), intent(in), optional :: af2(:,:)
    real(dp), allocatable, intent(out) :: lo(:,:), lh(:,:), hl(:,:), hh(:,:)
    real(dp), allocatable :: low_rows(:,:), high_rows(:,:)
    real(dp), allocatable :: a(:), b(:)
    integer :: i, j, nr, nc
    nr = size(x,1)
    nc = size(x,2)
    allocate(low_rows(nr/2,nc), high_rows(nr/2,nc))
    do j = 1, nc
      call afb(x(:,j), af1, a, b)
      low_rows(:,j) = a
      high_rows(:,j) = b
    end do
    allocate(lo(nr/2,nc/2), lh(nr/2,nc/2), hl(nr/2,nc/2), hh(nr/2,nc/2))
    do i = 1, nr/2
      if (present(af2)) then
        call afb(low_rows(i,:), af2, a, b)
      else
        call afb(low_rows(i,:), af1, a, b)
      end if
      lo(i,:) = a
      lh(i,:) = b
      if (present(af2)) then
        call afb(high_rows(i,:), af2, a, b)
      else
        call afb(high_rows(i,:), af1, a, b)
      end if
      hl(i,:) = a
      hh(i,:) = b
    end do
  end subroutine afb2d

  function sfb2d(lo, lh, hl, hh, sf1, sf2) result(x)
    real(dp), intent(in) :: lo(:,:), lh(:,:), hl(:,:), hh(:,:), sf1(:,:)
    real(dp), intent(in), optional :: sf2(:,:)
    real(dp), allocatable :: x(:,:)
    real(dp), allocatable :: low_rows(:,:), high_rows(:,:), a(:)
    integer :: i, j, nr, nc
    nr = 2*size(lo,1)
    nc = 2*size(lo,2)
    allocate(low_rows(nr/2,nc), high_rows(nr/2,nc))
    do i = 1, nr/2
      if (present(sf2)) then
        low_rows(i,:) = sfb(lo(i,:), lh(i,:), sf2)
        high_rows(i,:) = sfb(hl(i,:), hh(i,:), sf2)
      else
        low_rows(i,:) = sfb(lo(i,:), lh(i,:), sf1)
        high_rows(i,:) = sfb(hl(i,:), hh(i,:), sf1)
      end if
    end do
    allocate(x(nr,nc))
    do j = 1, nc
      a = sfb(low_rows(:,j), high_rows(:,j), sf1)
      x(:,j) = a
    end do
  end function sfb2d

  function dualtree_2d(x, n_levels) result(w)
    real(dp), intent(in) :: x(:,:)
    integer, intent(in) :: n_levels
    type(complex_wavelet_transform_2d) :: w
    type(dual_filter_pair) :: first, later
    real(dp), allocatable :: x1(:,:), x2(:,:), lo(:,:), lh(:,:), hl(:,:), hh(:,:)
    integer :: j
    if (n_levels < 1 .or. mod(size(x,1),2**n_levels) /= 0 .or. &
        mod(size(x,2),2**n_levels) /= 0) return
    first = fsfarras_filters()
    later = dualfilt1_filters()
    allocate(w%level(n_levels))
    w%original_shape = shape(x)
    x1 = x/sqrt2
    x2 = x/sqrt2
    call afb2d(x1, first%af1, first%af1, lo, lh, hl, hh)
    x1 = lo
    w%level(1)%lh = cmplx(lh,0.0_dp,dp)
    w%level(1)%hl = cmplx(hl,0.0_dp,dp)
    w%level(1)%hh = cmplx(hh,0.0_dp,dp)
    call afb2d(x2, first%af2, first%af2, lo, lh, hl, hh)
    x2 = lo
    w%level(1)%lh = cmplx(real(w%level(1)%lh,dp),lh,dp)
    w%level(1)%hl = cmplx(real(w%level(1)%hl,dp),hl,dp)
    w%level(1)%hh = cmplx(real(w%level(1)%hh,dp),hh,dp)
    do j = 2, n_levels
      call afb2d(x1, later%af1, later%af1, lo, lh, hl, hh)
      x1 = lo
      w%level(j)%lh = cmplx(lh,0.0_dp,dp)
      w%level(j)%hl = cmplx(hl,0.0_dp,dp)
      w%level(j)%hh = cmplx(hh,0.0_dp,dp)
      call afb2d(x2, later%af2, later%af2, lo, lh, hl, hh)
      x2 = lo
      w%level(j)%lh = cmplx(real(w%level(j)%lh,dp),lh,dp)
      w%level(j)%hl = cmplx(real(w%level(j)%hl,dp),hl,dp)
      w%level(j)%hh = cmplx(real(w%level(j)%hh,dp),hh,dp)
    end do
    w%smooth = cmplx(x1,x2,dp)
  end function dualtree_2d

  function idualtree_2d(w) result(x)
    type(complex_wavelet_transform_2d), intent(in) :: w
    real(dp), allocatable :: x(:,:)
    type(dual_filter_pair) :: first, later
    real(dp), allocatable :: x1(:,:), x2(:,:)
    integer :: j
    if (.not. allocated(w%level) .or. .not. allocated(w%smooth)) then
      allocate(x(0,0))
      return
    end if
    first = fsfarras_filters()
    later = dualfilt1_filters()
    x1 = real(w%smooth,dp)
    x2 = aimag(w%smooth)
    do j = size(w%level), 2, -1
      x1 = sfb2d(x1, real(w%level(j)%lh,dp), real(w%level(j)%hl,dp), &
        real(w%level(j)%hh,dp), later%sf1, later%sf1)
      x2 = sfb2d(x2, aimag(w%level(j)%lh), aimag(w%level(j)%hl), &
        aimag(w%level(j)%hh), later%sf2, later%sf2)
    end do
    x1 = sfb2d(x1, real(w%level(1)%lh,dp), real(w%level(1)%hl,dp), &
      real(w%level(1)%hh,dp), first%sf1, first%sf1)
    x2 = sfb2d(x2, aimag(w%level(1)%lh), aimag(w%level(1)%hl), &
      aimag(w%level(1)%hh), first%sf2, first%sf2)
    x = (x1+x2)/sqrt2
  end function idualtree_2d

end module waveslim_dualtree
