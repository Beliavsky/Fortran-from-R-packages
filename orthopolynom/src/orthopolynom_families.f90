module orthopolynom_families
  use polynom, only : dp, polynomial_t, polylist_t, polynomial
  use orthopolynom_types, only : recurrence_t
  use orthopolynom_core, only : pochhammer, orthogonal_polynomials_from_recurrence, &
    orthonormal_polynomials_from_recurrence
  implicit none
  private

  real(dp), parameter :: pi = acos(-1.0_dp)
  real(dp), parameter :: parameter_tol = 1.0e-6_dp

  public :: chebyshev_c_inner_products, chebyshev_c_recurrences, chebyshev_c_polynomials, chebyshev_c_weight
  public :: chebyshev_s_inner_products, chebyshev_s_recurrences, chebyshev_s_polynomials, chebyshev_s_weight
  public :: chebyshev_t_inner_products, chebyshev_t_recurrences, chebyshev_t_polynomials, chebyshev_t_weight
  public :: chebyshev_u_inner_products, chebyshev_u_recurrences, chebyshev_u_polynomials, chebyshev_u_weight
  public :: schebyshev_t_inner_products, schebyshev_t_recurrences, schebyshev_t_polynomials, schebyshev_t_weight
  public :: schebyshev_u_inner_products, schebyshev_u_recurrences, schebyshev_u_polynomials, schebyshev_u_weight
  public :: gegenbauer_inner_products, gegenbauer_recurrences, gegenbauer_polynomials, gegenbauer_weight
  public :: ultraspherical_inner_products, ultraspherical_recurrences, ultraspherical_polynomials, ultraspherical_weight
  public :: hermite_h_inner_products, hermite_h_recurrences, hermite_h_polynomials, hermite_h_weight
  public :: hermite_he_inner_products, hermite_he_recurrences, hermite_he_polynomials, hermite_he_weight
  public :: ghermite_h_inner_products, ghermite_h_recurrences, ghermite_h_polynomials, ghermite_h_weight
  public :: glaguerre_inner_products, glaguerre_recurrences, glaguerre_polynomials, glaguerre_weight
  public :: laguerre_inner_products, laguerre_recurrences, laguerre_polynomials, laguerre_weight
  public :: legendre_inner_products, legendre_recurrences, legendre_polynomials, legendre_weight
  public :: slegendre_inner_products, slegendre_recurrences, slegendre_polynomials, slegendre_weight
  public :: spherical_inner_products, spherical_recurrences, spherical_polynomials, spherical_weight
  public :: jacobi_p_inner_products, jacobi_p_recurrences, jacobi_p_polynomials, jacobi_p_weight
  public :: jacobi_g_inner_products, jacobi_g_recurrences, jacobi_g_polynomials, jacobi_g_weight

contains

  function allocate_recurrence(n) result(r)
    integer, intent(in) :: n
    type(recurrence_t) :: r
    call check_order(n)
    allocate(r%c(n + 1), r%d(n + 1), r%e(n + 1), r%f(n + 1))
    r%c = 0.0_dp
    r%d = 0.0_dp
    r%e = 0.0_dp
    r%f = 0.0_dp
  end function allocate_recurrence

  subroutine check_order(n)
    integer, intent(in) :: n
    if (n < 0) error stop 'orthopolynom: negative highest polynomial order'
  end subroutine check_order

  function make_family_polynomials(r, normalized, h0) result(polys)
    type(recurrence_t), intent(in) :: r
    logical, intent(in) :: normalized
    real(dp), intent(in) :: h0
    type(polylist_t) :: polys
    type(polynomial_t) :: p0

    if (normalized) then
      if (h0 <= 0.0_dp) error stop 'orthopolynom: nonpositive zeroth norm'
      p0 = polynomial(1.0_dp / sqrt(h0))
      polys = orthonormal_polynomials_from_recurrence(r, p0)
    else
      polys = orthogonal_polynomials_from_recurrence(r)
    end if
  end function make_family_polynomials

  pure real(dp) function factorial_real(k) result(v)
    integer, intent(in) :: k
    if (k < 0) error stop 'factorial_real: negative argument'
    v = exp(log_gamma(real(k + 1, dp)))
  end function factorial_real

  ! ---- Chebyshev C on (-2, 2) ----

  function chebyshev_c_inner_products(n) result(h)
    integer, intent(in) :: n
    real(dp), allocatable :: h(:)
    call check_order(n)
    allocate(h(n + 1), source=4.0_dp * pi)
    h(1) = 8.0_dp * pi
  end function chebyshev_c_inner_products

  function chebyshev_c_recurrences(n, normalized) result(r)
    integer, intent(in) :: n
    logical, intent(in), optional :: normalized
    type(recurrence_t) :: r
    logical :: norm
    integer :: j

    norm = .false.; if (present(normalized)) norm = normalized
    r = allocate_recurrence(n)
    r%c = 1.0_dp
    if (norm) then
      do j = 0, n
        if (j == 0) then
          r%e(j + 1) = 0.5_dp * sqrt(2.0_dp)
          r%f(j + 1) = 0.0_dp
        else
          r%e(j + 1) = 1.0_dp
          r%f(j + 1) = merge(sqrt(2.0_dp), 1.0_dp, j == 1)
        end if
      end do
    else
      r%e = 1.0_dp
      r%f = 1.0_dp
      r%e(1) = 0.5_dp
    end if
  end function chebyshev_c_recurrences

  function chebyshev_c_polynomials(n, normalized) result(polys)
    integer, intent(in) :: n
    logical, intent(in), optional :: normalized
    type(polylist_t) :: polys
    logical :: norm
    norm = .false.; if (present(normalized)) norm = normalized
    polys = make_family_polynomials(chebyshev_c_recurrences(n, norm), norm, 8.0_dp * pi)
  end function chebyshev_c_polynomials

  elemental real(dp) function chebyshev_c_weight(x) result(w)
    real(dp), intent(in) :: x
    w = 0.0_dp
    if (x > -2.0_dp .and. x < 2.0_dp) w = 1.0_dp / sqrt(1.0_dp - 0.25_dp * x * x)
  end function chebyshev_c_weight

  ! ---- Chebyshev S on (-2, 2) ----

  function chebyshev_s_inner_products(n) result(h)
    integer, intent(in) :: n
    real(dp), allocatable :: h(:)
    call check_order(n)
    allocate(h(n + 1), source=pi)
  end function chebyshev_s_inner_products

  function chebyshev_s_recurrences(n, normalized) result(r)
    integer, intent(in) :: n
    logical, intent(in), optional :: normalized
    type(recurrence_t) :: r
    logical :: norm
    norm = .false.; if (present(normalized)) norm = normalized
    r = allocate_recurrence(n)
    r%c = 1.0_dp
    r%e = 1.0_dp
    r%f = 1.0_dp
    if (norm) r%f(1) = 0.0_dp
  end function chebyshev_s_recurrences

  function chebyshev_s_polynomials(n, normalized) result(polys)
    integer, intent(in) :: n
    logical, intent(in), optional :: normalized
    type(polylist_t) :: polys
    logical :: norm
    norm = .false.; if (present(normalized)) norm = normalized
    polys = make_family_polynomials(chebyshev_s_recurrences(n, norm), norm, pi)
  end function chebyshev_s_polynomials

  elemental real(dp) function chebyshev_s_weight(x) result(w)
    real(dp), intent(in) :: x
    w = 0.0_dp
    if (x > -2.0_dp .and. x < 2.0_dp) w = sqrt(1.0_dp - 0.25_dp * x * x)
  end function chebyshev_s_weight

  ! ---- Chebyshev T on (-1, 1) ----

  function chebyshev_t_inner_products(n) result(h)
    integer, intent(in) :: n
    real(dp), allocatable :: h(:)
    call check_order(n)
    allocate(h(n + 1), source=0.5_dp * pi)
    h(1) = pi
  end function chebyshev_t_inner_products

  function chebyshev_t_recurrences(n, normalized) result(r)
    integer, intent(in) :: n
    logical, intent(in), optional :: normalized
    type(recurrence_t) :: r
    logical :: norm
    integer :: j

    norm = .false.; if (present(normalized)) norm = normalized
    r = allocate_recurrence(n)
    r%c = 1.0_dp
    if (norm) then
      do j = 0, n
        r%e(j + 1) = merge(sqrt(2.0_dp), 2.0_dp, j == 0)
        if (j == 0) then
          r%f(j + 1) = 0.0_dp
        else
          r%f(j + 1) = merge(sqrt(2.0_dp), 1.0_dp, j == 1)
        end if
      end do
    else
      r%e = 2.0_dp
      r%e(1) = 1.0_dp
      r%f = 1.0_dp
    end if
  end function chebyshev_t_recurrences

  function chebyshev_t_polynomials(n, normalized) result(polys)
    integer, intent(in) :: n
    logical, intent(in), optional :: normalized
    type(polylist_t) :: polys
    logical :: norm
    norm = .false.; if (present(normalized)) norm = normalized
    polys = make_family_polynomials(chebyshev_t_recurrences(n, norm), norm, pi)
  end function chebyshev_t_polynomials

  elemental real(dp) function chebyshev_t_weight(x) result(w)
    real(dp), intent(in) :: x
    w = 0.0_dp
    if (x > -1.0_dp .and. x < 1.0_dp) w = 1.0_dp / sqrt(1.0_dp - x * x)
  end function chebyshev_t_weight

  ! ---- Chebyshev U on (-1, 1) ----

  function chebyshev_u_inner_products(n) result(h)
    integer, intent(in) :: n
    real(dp), allocatable :: h(:)
    call check_order(n)
    allocate(h(n + 1), source=0.5_dp * pi)
  end function chebyshev_u_inner_products

  function chebyshev_u_recurrences(n, normalized) result(r)
    integer, intent(in) :: n
    logical, intent(in), optional :: normalized
    type(recurrence_t) :: r
    logical :: norm
    norm = .false.; if (present(normalized)) norm = normalized
    r = allocate_recurrence(n)
    r%c = 1.0_dp
    r%e = 2.0_dp
    r%f = 1.0_dp
    if (norm) r%f(1) = 0.0_dp
  end function chebyshev_u_recurrences

  function chebyshev_u_polynomials(n, normalized) result(polys)
    integer, intent(in) :: n
    logical, intent(in), optional :: normalized
    type(polylist_t) :: polys
    logical :: norm
    norm = .false.; if (present(normalized)) norm = normalized
    polys = make_family_polynomials(chebyshev_u_recurrences(n, norm), norm, 0.5_dp * pi)
  end function chebyshev_u_polynomials

  elemental real(dp) function chebyshev_u_weight(x) result(w)
    real(dp), intent(in) :: x
    w = 0.0_dp
    if (x > -1.0_dp .and. x < 1.0_dp) w = sqrt(1.0_dp - x * x)
  end function chebyshev_u_weight

  ! ---- Shifted Chebyshev on (0, 1) ----

  function schebyshev_t_inner_products(n) result(h)
    integer, intent(in) :: n
    real(dp), allocatable :: h(:)
    h = chebyshev_t_inner_products(n)
  end function schebyshev_t_inner_products

  function schebyshev_t_recurrences(n, normalized) result(r)
    integer, intent(in) :: n
    logical, intent(in), optional :: normalized
    type(recurrence_t) :: r
    logical :: norm
    integer :: j

    norm = .false.; if (present(normalized)) norm = normalized
    r = allocate_recurrence(n)
    r%c = 1.0_dp
    if (norm) then
      do j = 0, n
        if (j == 0) then
          r%d(1) = -sqrt(2.0_dp)
          r%e(1) = 2.0_dp * sqrt(2.0_dp)
          r%f(1) = 0.0_dp
        else
          r%d(j + 1) = -2.0_dp
          r%e(j + 1) = 4.0_dp
          r%f(j + 1) = merge(sqrt(2.0_dp), 1.0_dp, j == 1)
        end if
      end do
    else
      r%d = -2.0_dp
      r%d(1) = -1.0_dp
      r%e = 4.0_dp
      r%e(1) = 2.0_dp
      r%f = 1.0_dp
    end if
  end function schebyshev_t_recurrences

  function schebyshev_t_polynomials(n, normalized) result(polys)
    integer, intent(in) :: n
    logical, intent(in), optional :: normalized
    type(polylist_t) :: polys
    logical :: norm
    norm = .false.; if (present(normalized)) norm = normalized
    polys = make_family_polynomials(schebyshev_t_recurrences(n, norm), norm, pi)
  end function schebyshev_t_polynomials

  elemental real(dp) function schebyshev_t_weight(x) result(w)
    real(dp), intent(in) :: x
    w = 0.0_dp
    if (x > 0.0_dp .and. x < 1.0_dp) w = 1.0_dp / sqrt(x - x * x)
  end function schebyshev_t_weight

  function schebyshev_u_inner_products(n) result(h)
    integer, intent(in) :: n
    real(dp), allocatable :: h(:)
    call check_order(n)
    allocate(h(n + 1), source=pi / 8.0_dp)
  end function schebyshev_u_inner_products

  function schebyshev_u_recurrences(n, normalized) result(r)
    integer, intent(in) :: n
    logical, intent(in), optional :: normalized
    type(recurrence_t) :: r
    logical :: norm
    norm = .false.; if (present(normalized)) norm = normalized
    r = allocate_recurrence(n)
    r%c = 1.0_dp
    r%d = -2.0_dp
    r%e = 4.0_dp
    r%f = 1.0_dp
    if (norm) r%f(1) = 0.0_dp
  end function schebyshev_u_recurrences

  function schebyshev_u_polynomials(n, normalized) result(polys)
    integer, intent(in) :: n
    logical, intent(in), optional :: normalized
    type(polylist_t) :: polys
    logical :: norm
    norm = .false.; if (present(normalized)) norm = normalized
    polys = make_family_polynomials(schebyshev_u_recurrences(n, norm), norm, pi / 8.0_dp)
  end function schebyshev_u_polynomials

  elemental real(dp) function schebyshev_u_weight(x) result(w)
    real(dp), intent(in) :: x
    w = 0.0_dp
    if (x > 0.0_dp .and. x < 1.0_dp) w = sqrt(x - x * x)
  end function schebyshev_u_weight

  ! ---- Legendre and shifted Legendre ----

  function legendre_inner_products(n) result(h)
    integer, intent(in) :: n
    real(dp), allocatable :: h(:)
    integer :: k
    call check_order(n)
    allocate(h(n + 1))
    do k = 0, n
      h(k + 1) = 2.0_dp / real(2 * k + 1, dp)
    end do
  end function legendre_inner_products

  function legendre_recurrences(n, normalized) result(r)
    integer, intent(in) :: n
    logical, intent(in), optional :: normalized
    type(recurrence_t) :: r
    logical :: norm
    integer :: j

    norm = .false.; if (present(normalized)) norm = normalized
    r = allocate_recurrence(n)
    do j = 0, n
      r%c(j + 1) = real(j + 1, dp)
      if (norm) then
        r%e(j + 1) = sqrt(real((2 * j + 1) * (2 * j + 3), dp))
        if (j > 0) r%f(j + 1) = real(j, dp) * sqrt(real(2 * j + 3, dp) / real(2 * j - 1, dp))
      else
        r%e(j + 1) = real(2 * j + 1, dp)
        r%f(j + 1) = real(j, dp)
      end if
    end do
  end function legendre_recurrences

  function legendre_polynomials(n, normalized) result(polys)
    integer, intent(in) :: n
    logical, intent(in), optional :: normalized
    type(polylist_t) :: polys
    logical :: norm
    norm = .false.; if (present(normalized)) norm = normalized
    polys = make_family_polynomials(legendre_recurrences(n, norm), norm, 2.0_dp)
  end function legendre_polynomials

  elemental real(dp) function legendre_weight(x) result(w)
    real(dp), intent(in) :: x
    w = merge(1.0_dp, 0.0_dp, x > -1.0_dp .and. x < 1.0_dp)
  end function legendre_weight

  function slegendre_inner_products(n) result(h)
    integer, intent(in) :: n
    real(dp), allocatable :: h(:)
    integer :: k
    call check_order(n)
    allocate(h(n + 1))
    do k = 0, n
      h(k + 1) = 1.0_dp / real(2 * k + 1, dp)
    end do
  end function slegendre_inner_products

  function slegendre_recurrences(n, normalized) result(r)
    integer, intent(in) :: n
    logical, intent(in), optional :: normalized
    type(recurrence_t) :: r
    logical :: norm
    integer :: j

    norm = .false.; if (present(normalized)) norm = normalized
    r = allocate_recurrence(n)
    do j = 0, n
      r%c(j + 1) = real(j + 1, dp)
      if (norm) then
        r%d(j + 1) = -sqrt(real((2 * j + 3) * (2 * j + 1), dp))
        r%e(j + 1) = 2.0_dp * sqrt(real((2 * j + 3) * (2 * j + 1), dp))
        if (j > 0) r%f(j + 1) = real(j, dp) * sqrt(real(2 * j + 3, dp) / real(2 * j - 1, dp))
      else
        r%d(j + 1) = -real(2 * j + 1, dp)
        r%e(j + 1) = real(4 * j + 2, dp)
        r%f(j + 1) = real(j, dp)
      end if
    end do
  end function slegendre_recurrences

  function slegendre_polynomials(n, normalized) result(polys)
    integer, intent(in) :: n
    logical, intent(in), optional :: normalized
    type(polylist_t) :: polys
    logical :: norm
    norm = .false.; if (present(normalized)) norm = normalized
    polys = make_family_polynomials(slegendre_recurrences(n, norm), norm, 1.0_dp)
  end function slegendre_polynomials

  elemental real(dp) function slegendre_weight(x) result(w)
    real(dp), intent(in) :: x
    w = merge(1.0_dp, 0.0_dp, x > 0.0_dp .and. x < 1.0_dp)
  end function slegendre_weight

  function spherical_inner_products(n) result(h)
    integer, intent(in) :: n
    real(dp), allocatable :: h(:)
    h = legendre_inner_products(n)
  end function spherical_inner_products

  function spherical_recurrences(n, normalized) result(r)
    integer, intent(in) :: n
    logical, intent(in), optional :: normalized
    type(recurrence_t) :: r
    logical :: norm
    norm = .false.; if (present(normalized)) norm = normalized
    r = legendre_recurrences(n, norm)
  end function spherical_recurrences

  function spherical_polynomials(n, normalized) result(polys)
    integer, intent(in) :: n
    logical, intent(in), optional :: normalized
    type(polylist_t) :: polys
    logical :: norm
    norm = .false.; if (present(normalized)) norm = normalized
    polys = legendre_polynomials(n, norm)
  end function spherical_polynomials

  elemental real(dp) function spherical_weight(x) result(w)
    real(dp), intent(in) :: x
    w = legendre_weight(x)
  end function spherical_weight

  ! ---- Hermite ----

  function hermite_h_inner_products(n) result(h)
    integer, intent(in) :: n
    real(dp), allocatable :: h(:)
    integer :: k
    call check_order(n)
    allocate(h(n + 1))
    do k = 0, n
      h(k + 1) = sqrt(pi) * 2.0_dp**k * factorial_real(k)
    end do
  end function hermite_h_inner_products

  function hermite_h_recurrences(n, normalized) result(r)
    integer, intent(in) :: n
    logical, intent(in), optional :: normalized
    type(recurrence_t) :: r
    logical :: norm
    integer :: j

    norm = .false.; if (present(normalized)) norm = normalized
    r = allocate_recurrence(n)
    r%c = 1.0_dp
    do j = 0, n
      if (norm) then
        r%e(j + 1) = sqrt(2.0_dp / real(j + 1, dp))
        if (j > 0) r%f(j + 1) = sqrt(real(j, dp) / real(j + 1, dp))
      else
        r%e(j + 1) = 2.0_dp
        r%f(j + 1) = 2.0_dp * real(j, dp)
      end if
    end do
  end function hermite_h_recurrences

  function hermite_h_polynomials(n, normalized) result(polys)
    integer, intent(in) :: n
    logical, intent(in), optional :: normalized
    type(polylist_t) :: polys
    logical :: norm
    norm = .false.; if (present(normalized)) norm = normalized
    polys = make_family_polynomials(hermite_h_recurrences(n, norm), norm, sqrt(pi))
  end function hermite_h_polynomials

  elemental real(dp) function hermite_h_weight(x) result(w)
    real(dp), intent(in) :: x
    w = exp(-x * x)
  end function hermite_h_weight

  function hermite_he_inner_products(n) result(h)
    integer, intent(in) :: n
    real(dp), allocatable :: h(:)
    integer :: k
    call check_order(n)
    allocate(h(n + 1))
    do k = 0, n
      h(k + 1) = sqrt(2.0_dp * pi) * factorial_real(k)
    end do
  end function hermite_he_inner_products

  function hermite_he_recurrences(n, normalized) result(r)
    integer, intent(in) :: n
    logical, intent(in), optional :: normalized
    type(recurrence_t) :: r
    logical :: norm
    integer :: j

    norm = .false.; if (present(normalized)) norm = normalized
    r = allocate_recurrence(n)
    r%c = 1.0_dp
    do j = 0, n
      if (norm) then
        r%e(j + 1) = 1.0_dp / sqrt(real(j + 1, dp))
        if (j > 0) r%f(j + 1) = sqrt(real(j, dp) / real(j + 1, dp))
      else
        r%e(j + 1) = 1.0_dp
        r%f(j + 1) = real(j, dp)
      end if
    end do
  end function hermite_he_recurrences

  function hermite_he_polynomials(n, normalized) result(polys)
    integer, intent(in) :: n
    logical, intent(in), optional :: normalized
    type(polylist_t) :: polys
    logical :: norm
    norm = .false.; if (present(normalized)) norm = normalized
    polys = make_family_polynomials(hermite_he_recurrences(n, norm), norm, sqrt(2.0_dp * pi))
  end function hermite_he_polynomials

  elemental real(dp) function hermite_he_weight(x) result(w)
    real(dp), intent(in) :: x
    w = exp(-0.5_dp * x * x)
  end function hermite_he_weight

  ! ---- Generalized Hermite ----

  function ghermite_h_inner_products(n, mu) result(h)
    integer, intent(in) :: n
    real(dp), intent(in) :: mu
    real(dp), allocatable :: h(:)
    integer :: k, floor_k, floor_kp1
    real(dp) :: log_h

    call check_order(n)
    if (mu <= -0.5_dp) error stop 'ghermite_h_inner_products: mu must be greater than -0.5'
    if (abs(mu) < parameter_tol) then
      h = hermite_h_inner_products(n)
      return
    end if
    allocate(h(n + 1))
    do k = 0, n
      floor_k = k / 2
      floor_kp1 = (k + 1) / 2
      log_h = 2.0_dp * real(k, dp) * log(2.0_dp) + log_gamma(real(floor_k + 1, dp)) + &
        log_gamma(real(floor_kp1, dp) + mu + 0.5_dp)
      h(k + 1) = exp(log_h)
    end do
  end function ghermite_h_inner_products

  function ghermite_h_recurrences(n, mu, normalized) result(r)
    integer, intent(in) :: n
    real(dp), intent(in) :: mu
    logical, intent(in), optional :: normalized
    type(recurrence_t) :: r
    logical :: norm
    integer :: j, theta
    real(dp) :: two_mu

    call check_order(n)
    if (mu <= -0.5_dp) error stop 'ghermite_h_recurrences: mu must be greater than -0.5'
    norm = .false.; if (present(normalized)) norm = normalized
    if (abs(mu) < parameter_tol) then
      r = hermite_h_recurrences(n, norm)
      return
    end if
    r = allocate_recurrence(n)
    r%c = 1.0_dp
    two_mu = 2.0_dp * mu
    do j = 0, n
      theta = mod(j, 2)
      if (norm) then
        if (theta == 0) then
          r%e(j + 1) = sqrt(2.0_dp / (real(j + 1, dp) + two_mu))
          if (j > 0) r%f(j + 1) = sqrt(real(j, dp) / (real(j + 1, dp) + two_mu))
        else
          r%e(j + 1) = sqrt(2.0_dp / real(j + 1, dp))
          if (j > 0) r%f(j + 1) = sqrt((real(j, dp) + two_mu) / real(j + 1, dp))
        end if
      else
        r%e(j + 1) = 2.0_dp
        r%f(j + 1) = 2.0_dp * (real(j, dp) + two_mu * real(theta, dp))
      end if
    end do
  end function ghermite_h_recurrences

  function ghermite_h_polynomials(n, mu, normalized) result(polys)
    integer, intent(in) :: n
    real(dp), intent(in) :: mu
    logical, intent(in), optional :: normalized
    type(polylist_t) :: polys
    logical :: norm
    norm = .false.; if (present(normalized)) norm = normalized
    if (mu <= -0.5_dp) error stop 'ghermite_h_polynomials: mu must be greater than -0.5'
    polys = make_family_polynomials(ghermite_h_recurrences(n, mu, norm), norm, gamma(mu + 0.5_dp))
  end function ghermite_h_polynomials

  elemental real(dp) function ghermite_h_weight(x, mu) result(w)
    real(dp), intent(in) :: x, mu
    if (mu <= -0.5_dp) error stop 'ghermite_h_weight: mu must be greater than -0.5'
    if (abs(mu) < parameter_tol) then
      w = exp(-x * x)
    else
      w = abs(x)**(2.0_dp * mu) * exp(-x * x)
    end if
  end function ghermite_h_weight

  ! ---- Generalized Laguerre and Laguerre ----

  function glaguerre_inner_products(n, alpha) result(h)
    integer, intent(in) :: n
    real(dp), intent(in) :: alpha
    real(dp), allocatable :: h(:)
    integer :: k

    call check_order(n)
    if (alpha <= -1.0_dp) error stop 'glaguerre_inner_products: alpha must be greater than -1'
    allocate(h(n + 1))
    do k = 0, n
      h(k + 1) = exp(log_gamma(alpha + real(k + 1, dp)) - log_gamma(real(k + 1, dp)))
    end do
  end function glaguerre_inner_products

  function glaguerre_recurrences(n, alpha, normalized) result(r)
    integer, intent(in) :: n
    real(dp), intent(in) :: alpha
    logical, intent(in), optional :: normalized
    type(recurrence_t) :: r
    logical :: norm
    integer :: j
    real(dp) :: rho

    call check_order(n)
    if (alpha <= -1.0_dp) error stop 'glaguerre_recurrences: alpha must be greater than -1'
    norm = .false.; if (present(normalized)) norm = normalized
    r = allocate_recurrence(n)
    do j = 0, n
      r%c(j + 1) = real(j + 1, dp)
      if (norm) then
        rho = sqrt(real(j + 1, dp) / (alpha + real(j + 1, dp)))
        r%d(j + 1) = (2.0_dp * real(j, dp) + alpha + 1.0_dp) * rho
        r%e(j + 1) = -rho
        if (j > 0) then
          r%f(j + 1) = sqrt(real(j * (j + 1), dp) * (real(j, dp) + alpha) / &
            (alpha + real(j + 1, dp)))
        end if
      else
        r%d(j + 1) = 2.0_dp * real(j, dp) + alpha + 1.0_dp
        r%e(j + 1) = -1.0_dp
        r%f(j + 1) = real(j, dp) + alpha
      end if
    end do
  end function glaguerre_recurrences

  function glaguerre_polynomials(n, alpha, normalized) result(polys)
    integer, intent(in) :: n
    real(dp), intent(in) :: alpha
    logical, intent(in), optional :: normalized
    type(polylist_t) :: polys
    logical :: norm
    norm = .false.; if (present(normalized)) norm = normalized
    if (alpha <= -1.0_dp) error stop 'glaguerre_polynomials: alpha must be greater than -1'
    polys = make_family_polynomials(glaguerre_recurrences(n, alpha, norm), norm, gamma(alpha + 1.0_dp))
  end function glaguerre_polynomials

  elemental real(dp) function glaguerre_weight(x, alpha) result(w)
    real(dp), intent(in) :: x, alpha
    w = 0.0_dp
    if (x > 0.0_dp) w = exp(-x) * x**alpha
  end function glaguerre_weight

  function laguerre_inner_products(n) result(h)
    integer, intent(in) :: n
    real(dp), allocatable :: h(:)
    h = glaguerre_inner_products(n, 0.0_dp)
  end function laguerre_inner_products

  function laguerre_recurrences(n, normalized) result(r)
    integer, intent(in) :: n
    logical, intent(in), optional :: normalized
    type(recurrence_t) :: r
    logical :: norm
    norm = .false.; if (present(normalized)) norm = normalized
    r = glaguerre_recurrences(n, 0.0_dp, norm)
  end function laguerre_recurrences

  function laguerre_polynomials(n, normalized) result(polys)
    integer, intent(in) :: n
    logical, intent(in), optional :: normalized
    type(polylist_t) :: polys
    logical :: norm
    norm = .false.; if (present(normalized)) norm = normalized
    polys = glaguerre_polynomials(n, 0.0_dp, norm)
  end function laguerre_polynomials

  elemental real(dp) function laguerre_weight(x) result(w)
    real(dp), intent(in) :: x
    w = 0.0_dp
    if (x > 0.0_dp) w = exp(-x)
  end function laguerre_weight

  ! ---- Gegenbauer / ultraspherical ----

  function gegenbauer_inner_products(n, alpha) result(h)
    integer, intent(in) :: n
    real(dp), intent(in) :: alpha
    real(dp), allocatable :: h(:)
    integer :: k
    real(dp) :: coef

    call check_order(n)
    if (alpha <= -0.5_dp) error stop 'gegenbauer_inner_products: alpha must be greater than -0.5'
    if (abs(alpha - 0.5_dp) < parameter_tol) then
      h = legendre_inner_products(n)
      return
    end if
    if (abs(alpha - 1.0_dp) < parameter_tol) then
      h = chebyshev_u_inner_products(n)
      return
    end if
    allocate(h(n + 1))
    if (abs(alpha) < parameter_tol) then
      h(1) = pi
      do k = 1, n
        h(k + 1) = 2.0_dp * pi / real(k * k, dp)
      end do
      return
    end if
    coef = pi * 2.0_dp**(1.0_dp - 2.0_dp * alpha)
    do k = 0, n
      h(k + 1) = coef * gamma(real(k, dp) + 2.0_dp * alpha) / &
        (factorial_real(k) * (real(k, dp) + alpha) * gamma(alpha)**2)
    end do
  end function gegenbauer_inner_products

  function gegenbauer_recurrences(n, alpha, normalized) result(r)
    integer, intent(in) :: n
    real(dp), intent(in) :: alpha
    logical, intent(in), optional :: normalized
    type(recurrence_t) :: r
    logical :: norm
    integer :: j
    real(dp) :: rho_j, rho_jm1, two_alpha

    call check_order(n)
    if (alpha <= -0.5_dp) error stop 'gegenbauer_recurrences: alpha must be greater than -0.5'
    norm = .false.; if (present(normalized)) norm = normalized
    if (abs(alpha - 1.0_dp) < parameter_tol) then
      r = chebyshev_u_recurrences(n, norm)
      return
    end if
    if (abs(alpha - 0.5_dp) < parameter_tol) then
      r = legendre_recurrences(n, norm)
      return
    end if
    r = allocate_recurrence(n)
    if (abs(alpha) < parameter_tol) then
      do j = 0, n
        r%c(j + 1) = real(j + 1, dp)
        if (norm) then
          r%e(j + 1) = merge(sqrt(2.0_dp), 2.0_dp * real(j + 1, dp), j == 0)
          if (j == 0) then
            r%f(j + 1) = 0.0_dp
          else if (j == 1) then
            r%f(j + 1) = 2.0_dp * sqrt(2.0_dp)
          else
            r%f(j + 1) = real(j + 1, dp)
          end if
        else
          r%e(j + 1) = merge(2.0_dp, 2.0_dp * real(j, dp), j == 0)
          if (j == 0) then
            r%f(j + 1) = 0.0_dp
          else if (j == 1) then
            r%f(j + 1) = 2.0_dp
          else
            r%f(j + 1) = real(j - 1, dp)
          end if
        end if
      end do
      return
    end if
    two_alpha = 2.0_dp * alpha
    do j = 0, n
      r%c(j + 1) = real(j + 1, dp)
      if (norm) then
        rho_j = sqrt((alpha + real(j + 1, dp)) * real(j + 1, dp) * &
          gamma(two_alpha + real(j, dp)) / ((alpha + real(j, dp)) * &
          gamma(two_alpha + real(j + 1, dp))))
        r%e(j + 1) = 2.0_dp * (alpha + real(j, dp)) * rho_j
        if (j > 0) then
          rho_jm1 = sqrt(real(j * (j + 1), dp) * (alpha + real(j + 1, dp)) * &
            gamma(two_alpha + real(j - 1, dp)) / ((alpha + real(j - 1, dp)) * &
            gamma(two_alpha + real(j + 1, dp))))
          r%f(j + 1) = (real(j, dp) + two_alpha - 1.0_dp) * rho_jm1
        end if
      else
        r%e(j + 1) = 2.0_dp * (real(j, dp) + alpha)
        r%f(j + 1) = real(j, dp) + two_alpha - 1.0_dp
      end if
    end do
  end function gegenbauer_recurrences

  function gegenbauer_polynomials(n, alpha, normalized) result(polys)
    integer, intent(in) :: n
    real(dp), intent(in) :: alpha
    logical, intent(in), optional :: normalized
    type(polylist_t) :: polys
    logical :: norm
    real(dp) :: h0

    norm = .false.; if (present(normalized)) norm = normalized
    if (alpha <= -0.5_dp) error stop 'gegenbauer_polynomials: alpha must be greater than -0.5'
    if (abs(alpha - 0.5_dp) < parameter_tol) then
      polys = legendre_polynomials(n, norm)
      return
    end if
    if (abs(alpha - 1.0_dp) < parameter_tol) then
      polys = chebyshev_u_polynomials(n, norm)
      return
    end if
    if (abs(alpha) < parameter_tol) then
      h0 = pi
    else
      h0 = sqrt(pi) * gamma(alpha + 0.5_dp) / gamma(alpha + 1.0_dp)
    end if
    polys = make_family_polynomials(gegenbauer_recurrences(n, alpha, norm), norm, h0)
  end function gegenbauer_polynomials

  elemental real(dp) function gegenbauer_weight(x, alpha) result(w)
    real(dp), intent(in) :: x, alpha
    w = 0.0_dp
    if (x > -1.0_dp .and. x < 1.0_dp) w = (1.0_dp - x * x)**(alpha - 0.5_dp)
  end function gegenbauer_weight

  function ultraspherical_inner_products(n, alpha) result(h)
    integer, intent(in) :: n
    real(dp), intent(in) :: alpha
    real(dp), allocatable :: h(:)
    h = gegenbauer_inner_products(n, alpha)
  end function ultraspherical_inner_products

  function ultraspherical_recurrences(n, alpha, normalized) result(r)
    integer, intent(in) :: n
    real(dp), intent(in) :: alpha
    logical, intent(in), optional :: normalized
    type(recurrence_t) :: r
    logical :: norm
    norm = .false.; if (present(normalized)) norm = normalized
    r = gegenbauer_recurrences(n, alpha, norm)
  end function ultraspherical_recurrences

  function ultraspherical_polynomials(n, alpha, normalized) result(polys)
    integer, intent(in) :: n
    real(dp), intent(in) :: alpha
    logical, intent(in), optional :: normalized
    type(polylist_t) :: polys
    logical :: norm
    norm = .false.; if (present(normalized)) norm = normalized
    polys = gegenbauer_polynomials(n, alpha, norm)
  end function ultraspherical_polynomials

  elemental real(dp) function ultraspherical_weight(x, alpha) result(w)
    real(dp), intent(in) :: x, alpha
    if (abs(x) < 1.0_dp) then
      w = exp((alpha - 0.5_dp) * log(1.0_dp - x * x))
    else
      w = 0.0_dp
    end if
  end function ultraspherical_weight

  ! ---- Jacobi P ----

  function jacobi_p_inner_products(n, alpha, beta) result(h)
    integer, intent(in) :: n
    real(dp), intent(in) :: alpha, beta
    real(dp), allocatable :: h(:)
    integer :: k
    real(dp) :: ab, abp1, coef

    call check_order(n)
    if (alpha <= -1.0_dp .or. beta <= -1.0_dp) error stop 'jacobi_p_inner_products: alpha and beta must exceed -1'
    if (abs(alpha) < parameter_tol .and. abs(beta) < parameter_tol) then
      h = legendre_inner_products(n)
      return
    end if
    if (abs(alpha - beta) < parameter_tol) then
      h = gegenbauer_inner_products(n, alpha + 0.5_dp)
      return
    end if
    ab = alpha + beta
    abp1 = ab + 1.0_dp
    coef = 2.0_dp**abp1
    allocate(h(n + 1))
    do k = 0, n
      h(k + 1) = coef * gamma(real(k, dp) + alpha + 1.0_dp) * gamma(real(k, dp) + beta + 1.0_dp) / &
        ((2.0_dp * real(k, dp) + abp1) * factorial_real(k) * gamma(real(k, dp) + abp1))
    end do
  end function jacobi_p_inner_products

  function jacobi_p_recurrences(n, alpha, beta, normalized) result(r)
    integer, intent(in) :: n
    real(dp), intent(in) :: alpha, beta
    logical, intent(in), optional :: normalized
    type(recurrence_t) :: r
    logical :: norm, ab_zero
    integer :: j
    real(dp) :: ab, aabb, c0, d0, e0, f0, rho_j, rho_jm1, num, den

    norm = .false.; if (present(normalized)) norm = normalized
    if (abs(alpha) < parameter_tol .and. abs(beta) < parameter_tol) then
      r = legendre_recurrences(n, norm)
      return
    end if
    if (abs(alpha - beta) < parameter_tol) then
      r = gegenbauer_recurrences(n, alpha + 0.5_dp, norm)
      return
    end if
    call check_order(n)
    if (alpha <= -1.0_dp .or. beta <= -1.0_dp) error stop 'jacobi_p_recurrences: alpha and beta must exceed -1'
    r = allocate_recurrence(n)
    ab = alpha + beta
    ab_zero = abs(ab) < parameter_tol
    aabb = alpha * alpha - beta * beta
    do j = 0, n
      c0 = 2.0_dp * real(j + 1, dp) * (real(j + 1, dp) + ab) * (2.0_dp * real(j, dp) + ab)
      d0 = (2.0_dp * real(j, dp) + ab + 1.0_dp) * aabb
      e0 = pochhammer(2.0_dp * real(j, dp) + ab, 3)
      f0 = 2.0_dp * (real(j, dp) + alpha) * (real(j, dp) + beta) * &
        (2.0_dp * real(j, dp) + ab + 2.0_dp)
      if (ab_zero .and. j == 0) then
        c0 = 1.0_dp; d0 = alpha; e0 = 1.0_dp; f0 = 0.0_dp
      end if
      if (norm) then
        if (j == 0) then
          rho_j = sqrt((ab + 3.0_dp) / ((alpha + 1.0_dp) * (beta + 1.0_dp)))
        else
          num = real(j + 1, dp) * (2.0_dp * real(j, dp) + ab + 3.0_dp) * (real(j + 1, dp) + ab)
          den = (real(j + 1, dp) + alpha) * (real(j + 1, dp) + beta) * &
            (2.0_dp * real(j, dp) + ab + 1.0_dp)
          rho_j = sqrt(num / den)
        end if
        r%c(j + 1) = c0
        r%d(j + 1) = d0 * rho_j
        r%e(j + 1) = e0 * rho_j
        if (j == 0) then
          r%f(j + 1) = 0.0_dp
        else
          if (j == 1) then
            num = 2.0_dp * (ab + 5.0_dp) * (ab + 2.0_dp)
            den = (alpha + 1.0_dp) * (alpha + 2.0_dp) * (beta + 1.0_dp) * (beta + 2.0_dp)
          else
            num = real(j * (j + 1), dp) * (2.0_dp * real(j, dp) + ab + 3.0_dp) * &
              gamma(real(j, dp) + ab + 2.0_dp)
            den = (2.0_dp * real(j, dp) + ab - 1.0_dp) * gamma(real(j, dp) + ab) * &
              (real(j, dp) + alpha + 1.0_dp) * (real(j, dp) + beta + 1.0_dp) * &
              (real(j, dp) + alpha) * (real(j, dp) + beta)
          end if
          rho_jm1 = sqrt(num / den)
          r%f(j + 1) = f0 * rho_jm1
        end if
      else
        r%c(j + 1) = c0
        r%d(j + 1) = d0
        r%e(j + 1) = e0
        r%f(j + 1) = f0
      end if
    end do
  end function jacobi_p_recurrences

  function jacobi_p_polynomials(n, alpha, beta, normalized) result(polys)
    integer, intent(in) :: n
    real(dp), intent(in) :: alpha, beta
    logical, intent(in), optional :: normalized
    type(polylist_t) :: polys
    logical :: norm
    real(dp) :: h0

    norm = .false.; if (present(normalized)) norm = normalized
    if (abs(alpha) < parameter_tol .and. abs(beta) < parameter_tol) then
      polys = legendre_polynomials(n, norm)
      return
    end if
    if (abs(alpha - beta) < parameter_tol) then
      polys = gegenbauer_polynomials(n, alpha + 0.5_dp, norm)
      return
    end if
    if (alpha <= -1.0_dp .or. beta <= -1.0_dp) error stop 'jacobi_p_polynomials: alpha and beta must exceed -1'
    h0 = 2.0_dp**(alpha + beta + 1.0_dp) * gamma(alpha + 1.0_dp) * gamma(beta + 1.0_dp) / &
      gamma(alpha + beta + 2.0_dp)
    polys = make_family_polynomials(jacobi_p_recurrences(n, alpha, beta, norm), norm, h0)
  end function jacobi_p_polynomials

  elemental real(dp) function jacobi_p_weight(x, alpha, beta) result(w)
    real(dp), intent(in) :: x, alpha, beta
    if (alpha < -1.0_dp) error stop 'jacobi_p_weight: alpha is less than -1'
    if (beta < -1.0_dp) error stop 'jacobi_p_weight: beta is less than -1'
    w = 0.0_dp
    if (x > -1.0_dp .and. x < 1.0_dp) w = (1.0_dp - x)**alpha * (1.0_dp + x)**beta
  end function jacobi_p_weight

  ! ---- Jacobi G on (0, 1) ----

  function jacobi_g_inner_products(n, p, q) result(h)
    integer, intent(in) :: n
    real(dp), intent(in) :: p, q
    real(dp), allocatable :: h(:)
    integer :: k
    real(dp) :: pmq

    call check_order(n)
    if (p - q <= -1.0_dp) error stop 'jacobi_g_inner_products: p-q must exceed -1'
    if (q <= 0.0_dp) error stop 'jacobi_g_inner_products: q must be positive'
    pmq = p - q
    allocate(h(n + 1))
    h(1) = gamma(q) * gamma(pmq + 1.0_dp) / gamma(p + 1.0_dp)
    do k = 1, n
      h(k + 1) = factorial_real(k) * gamma(real(k, dp) + q) * gamma(real(k, dp) + p) * &
        gamma(real(k, dp) + pmq + 1.0_dp) / ((2.0_dp * real(k, dp) + p) * &
        gamma(2.0_dp * real(k, dp) + p)**2)
    end do
  end function jacobi_g_inner_products

  function jacobi_g_recurrences(n, p, q, normalized) result(r)
    integer, intent(in) :: n
    real(dp), intent(in) :: p, q
    logical, intent(in), optional :: normalized
    type(recurrence_t) :: r
    logical :: norm
    real(dp), allocatable :: norms(:)
    integer :: j, k
    real(dp) :: c0, d0, e0, f0, rho_j, rho_jm1

    call check_order(n)
    if (p - q <= -1.0_dp) error stop 'jacobi_g_recurrences: p-q must exceed -1'
    if (q <= 0.0_dp) error stop 'jacobi_g_recurrences: q must be positive'
    norm = .false.; if (present(normalized)) norm = normalized
    r = allocate_recurrence(n)
    if (norm) norms = sqrt(jacobi_g_inner_products(n + 1, p, q))

    if (abs(p) < parameter_tol) then
      c0 = 1.0_dp; d0 = -q; e0 = 1.0_dp; f0 = 0.0_dp
      if (norm) then
        rho_j = sqrt(2.0_dp / (q * (1.0_dp - q)))
        r%c(1) = c0; r%d(1) = d0 * rho_j; r%e(1) = e0 * rho_j; r%f(1) = 0.0_dp
      else
        r%c(1) = c0; r%d(1) = d0; r%e(1) = e0; r%f(1) = 0.0_dp
      end if
      if (n >= 1) then
        c0 = 1.0_dp; d0 = (q - 2.0_dp) / 3.0_dp; e0 = 1.0_dp; f0 = -q * (q - 1.0_dp) / 2.0_dp
        if (norm) then
          rho_j = 6.0_dp / sqrt((1.0_dp + q) * (2.0_dp - q))
          rho_jm1 = 6.0_dp * sqrt(2.0_dp / (q * (2.0_dp - q) * (1.0_dp - q * q)))
          r%c(2) = c0; r%d(2) = d0 * rho_j; r%e(2) = e0 * rho_j; r%f(2) = f0 * rho_jm1
        else
          r%c(2) = c0; r%d(2) = d0; r%e(2) = e0; r%f(2) = f0
        end if
      end if
      do j = 2, n
        k = j + 1
        c0 = pochhammer(2.0_dp * real(j, dp) - 2.0_dp, 4) * (2.0_dp * real(j, dp) - 1.0_dp)
        d0 = -(2.0_dp * real(j * j, dp) - q) * pochhammer(2.0_dp * real(j, dp) - 2.0_dp, 3)
        e0 = c0
        f0 = real(j, dp) * (real(j, dp) + q - 1.0_dp) * real(j - 1, dp) * &
          (real(j, dp) - q) * (2.0_dp * real(j, dp) + 1.0_dp)
        if (norm) then
          rho_j = norms(k) / norms(k + 1)
          rho_jm1 = norms(k - 1) / norms(k + 1)
          r%c(k) = c0; r%d(k) = d0 * rho_j; r%e(k) = e0 * rho_j; r%f(k) = f0 * rho_jm1
        else
          r%c(k) = c0; r%d(k) = d0; r%e(k) = e0; r%f(k) = f0
        end if
      end do
      return
    end if

    do j = 0, n
      k = j + 1
      if (j == 0) then
        c0 = 1.0_dp
        d0 = -q / (p + 1.0_dp)
        e0 = 1.0_dp
        f0 = 0.0_dp
      else
        c0 = pochhammer(2.0_dp * real(j, dp) + p - 2.0_dp, 4) * (2.0_dp * real(j, dp) + p - 1.0_dp)
        d0 = -(2.0_dp * real(j, dp) * (real(j, dp) + p) + q * (p - 1.0_dp)) * &
          pochhammer(2.0_dp * real(j, dp) + p - 2.0_dp, 3)
        e0 = c0
        f0 = real(j, dp) * (real(j, dp) + q - 1.0_dp) * (real(j, dp) + p - 1.0_dp) * &
          (real(j, dp) + p - q) * (2.0_dp * real(j, dp) + p + 1.0_dp)
      end if
      if (norm) then
        rho_j = norms(k) / norms(k + 1)
        r%c(k) = c0; r%d(k) = d0 * rho_j; r%e(k) = e0 * rho_j
        if (j == 0) then
          r%f(k) = 0.0_dp
        else
          rho_jm1 = norms(k - 1) / norms(k + 1)
          r%f(k) = f0 * rho_jm1
        end if
      else
        r%c(k) = c0; r%d(k) = d0; r%e(k) = e0; r%f(k) = f0
      end if
    end do
  end function jacobi_g_recurrences

  function jacobi_g_polynomials(n, p, q, normalized) result(polys)
    integer, intent(in) :: n
    real(dp), intent(in) :: p, q
    logical, intent(in), optional :: normalized
    type(polylist_t) :: polys
    logical :: norm
    real(dp) :: h0

    norm = .false.; if (present(normalized)) norm = normalized
    if (p - q <= -1.0_dp) error stop 'jacobi_g_polynomials: p-q must exceed -1'
    if (q <= 0.0_dp) error stop 'jacobi_g_polynomials: q must be positive'
    h0 = gamma(q) * gamma(p - q + 1.0_dp) / gamma(p + 1.0_dp)
    polys = make_family_polynomials(jacobi_g_recurrences(n, p, q, norm), norm, h0)
  end function jacobi_g_polynomials

  elemental real(dp) function jacobi_g_weight(x, p, q) result(w)
    real(dp), intent(in) :: x, p, q
    w = 0.0_dp
    if (x > 0.0_dp .and. x < 1.0_dp) w = (1.0_dp - x)**(p - q) * x**(q - 1.0_dp)
  end function jacobi_g_weight

end module orthopolynom_families
