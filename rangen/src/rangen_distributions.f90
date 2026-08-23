module rangen_distributions
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use rangen_kinds, only : dp, i8
    use rangen_pcg32, only : pcg32_state
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: euler_gamma = 0.577215664901532860606512090082402431_dp

    type(pcg32_state), save :: dist_rng
    type(pcg32_state), save :: runif_rng
    type(pcg32_state), save :: normal_rng
    logical, save :: initialized = .false.
    logical, save :: has_spare = .false.
    real(dp), save :: spare_normal = 0.0_dp

    public :: set_seed
    public :: runif, rbeta, rexp, rchisq, rgamma, rgeom, rcauchy, rt
    public :: rpareto, rfrechet, rlaplace, rgumbel, rgumble, rarcsine, rnorm
    public :: runif_mat, rbeta_mat, rexp_mat, rchisq_mat, rgamma_mat, rgeom_mat
    public :: rcauchy_mat, rt_mat, rpareto_mat, rfrechet_mat, rlaplace_mat
    public :: rgumbel_mat, rgumble_mat, rarcsine_mat, rnorm_mat
    public :: col_runif, col_rbeta, col_rexp, col_rchisq, col_rgamma, col_rgeom
    public :: col_rcauchy, col_rt, col_rpareto, col_rfrechet, col_rlaplace
    public :: col_rgumbel, col_rgumble, col_rarcsine, col_rnorm
    public :: euler_gamma

contains

    subroutine ensure_initialized()
        integer(i8) :: seed_value
        integer :: count
        if (initialized) return
        call system_clock(count)
        seed_value = int(count, i8)
        call set_seed(seed_value)
    end subroutine ensure_initialized

    subroutine set_seed(seed_value)
        integer(i8), intent(in) :: seed_value
        call dist_rng%seed(seed_value, 1442695040888963407_i8)
        call runif_rng%seed(seed_value, 6364136223846793005_i8)
        call normal_rng%seed(seed_value, 3202034522624059733_i8)
        has_spare = .false.
        initialized = .true.
    end subroutine set_seed

    function nan_value() result(x)
        real(dp) :: x
        x = ieee_value(0.0_dp, ieee_quiet_nan)
    end function nan_value

    function uniform_open() result(u)
        real(dp) :: u
        call ensure_initialized()
        u = dist_rng%uniform_open()
    end function uniform_open

    function normal_standard() result(z)
        real(dp) :: z
        real(dp) :: u1, u2, r
        call ensure_initialized()
        if (has_spare) then
            z = spare_normal
            has_spare = .false.
            return
        end if
        u1 = normal_rng%uniform_open()
        u2 = normal_rng%uniform_open()
        r = sqrt(-2.0_dp * log(u1))
        z = r * cos(2.0_dp * pi * u2)
        spare_normal = r * sin(2.0_dp * pi * u2)
        has_spare = .true.
    end function normal_standard

    recursive function gamma_scalar(shape, rate) result(x)
        real(dp), intent(in) :: shape, rate
        real(dp) :: x
        real(dp) :: d, c, z, v, u

        if (shape <= 0.0_dp .or. rate <= 0.0_dp) then
            x = nan_value()
            return
        end if
        if (shape < 1.0_dp) then
            x = gamma_scalar(shape + 1.0_dp, rate) * uniform_open() ** (1.0_dp / shape)
            return
        end if
        d = shape - 1.0_dp / 3.0_dp
        c = 1.0_dp / sqrt(9.0_dp * d)
        do
            z = normal_standard()
            v = 1.0_dp + c * z
            if (v <= 0.0_dp) cycle
            v = v * v * v
            u = uniform_open()
            if (u < 1.0_dp - 0.0331_dp * z ** 4) exit
            if (log(u) < 0.5_dp * z * z + d * (1.0_dp - v + log(v))) exit
        end do
        x = d * v / rate
    end function gamma_scalar

    function beta_scalar(alpha, beta) result(x)
        real(dp), intent(in) :: alpha, beta
        real(dp) :: x, a, b
        if (alpha <= 0.0_dp .or. beta <= 0.0_dp) then
            x = nan_value()
            return
        end if
        a = gamma_scalar(alpha, 1.0_dp)
        b = gamma_scalar(beta, 1.0_dp)
        x = a / (a + b)
    end function beta_scalar

    function runif_scalar(minval, maxval) result(x)
        real(dp), intent(in) :: minval, maxval
        real(dp) :: x
        call ensure_initialized()
        if (maxval < minval) then
            x = nan_value()
        else
            x = minval + (maxval - minval) * runif_rng%uniform_closed()
        end if
    end function runif_scalar

    function rexp_scalar(rate) result(x)
        real(dp), intent(in) :: rate
        real(dp) :: x
        if (rate <= 0.0_dp) then
            x = nan_value()
        else
            x = -log(uniform_open()) / rate
        end if
    end function rexp_scalar

    function rgeom_scalar(prob) result(x)
        real(dp), intent(in) :: prob
        real(dp) :: x
        if (prob <= 0.0_dp .or. prob > 1.0_dp) then
            x = nan_value()
        else if (prob >= 1.0_dp) then
            x = 0.0_dp
        else
            x = floor(log(uniform_open()) / log(1.0_dp - prob))
        end if
    end function rgeom_scalar

    function rcauchy_scalar(location, scale) result(x)
        real(dp), intent(in) :: location, scale
        real(dp) :: x
        if (scale < 0.0_dp) then
            x = nan_value()
        else
            x = location + scale * tan(pi * (uniform_open() - 0.5_dp))
        end if
    end function rcauchy_scalar

    function rt_scalar(df, ncp) result(x)
        real(dp), intent(in) :: df, ncp
        real(dp) :: x, y
        if (df <= 0.0_dp) then
            x = nan_value()
            return
        end if
        y = gamma_scalar(0.5_dp * df, 0.5_dp)
        x = (normal_standard() + ncp) / sqrt(y / df)
    end function rt_scalar

    function rpareto_scalar(shape, scale) result(x)
        real(dp), intent(in) :: shape, scale
        real(dp) :: x
        if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
            x = nan_value()
        else
            x = scale * (1.0_dp - uniform_open()) ** (-1.0_dp / shape)
        end if
    end function rpareto_scalar

    function rfrechet_scalar(lambda, mu, sigma) result(x)
        real(dp), intent(in) :: lambda, mu, sigma
        real(dp) :: x
        if (lambda <= 0.0_dp .or. sigma <= 0.0_dp) then
            x = nan_value()
        else
            x = mu + sigma * (-log(uniform_open())) ** (-1.0_dp / lambda)
        end if
    end function rfrechet_scalar

    function rlaplace_scalar(mu, sigma) result(x)
        real(dp), intent(in) :: mu, sigma
        real(dp) :: x, u
        if (sigma < 0.0_dp) then
            x = nan_value()
            return
        end if
        u = uniform_open() - 0.5_dp
        if (u >= 0.0_dp) then
            x = mu - sigma * log(1.0_dp - 2.0_dp * u)
        else
            x = mu + sigma * log(1.0_dp + 2.0_dp * u)
        end if
    end function rlaplace_scalar

    function rgumbel_scalar(mu, sigma) result(x)
        real(dp), intent(in) :: mu, sigma
        real(dp) :: x
        if (sigma <= 0.0_dp) then
            x = nan_value()
        else
            x = mu - sigma * log(-log(uniform_open()))
        end if
    end function rgumbel_scalar

    function rarcsine_scalar(minval, maxval) result(x)
        real(dp), intent(in) :: minval, maxval
        real(dp) :: x, s
        if (maxval < minval) then
            x = nan_value()
        else
            s = sin(0.5_dp * pi * uniform_open())
            x = minval + (maxval - minval) * s * s
        end if
    end function rarcsine_scalar

    function rnorm_scalar(mean, sd) result(x)
        real(dp), intent(in) :: mean, sd
        real(dp) :: x
        if (sd < 0.0_dp) then
            x = nan_value()
        else
            x = mean + sd * normal_standard()
        end if
    end function rnorm_scalar

    function runif(n, minval, maxval) result(x)
        integer, intent(in) :: n
        real(dp), intent(in), optional :: minval, maxval
        real(dp), allocatable :: x(:)
        real(dp) :: a, b
        integer :: i
        a = 0.0_dp
        b = 1.0_dp
        if (present(minval)) a = minval
        if (present(maxval)) b = maxval
        allocate(x(max(n, 0)))
        do i = 1, size(x)
            x(i) = runif_scalar(a, b)
        end do
    end function runif

    function rbeta(n, alpha, beta) result(x)
        integer, intent(in) :: n
        real(dp), intent(in) :: alpha, beta
        real(dp), allocatable :: x(:)
        integer :: i
        allocate(x(max(n, 0)))
        do i = 1, size(x)
            x(i) = beta_scalar(alpha, beta)
        end do
    end function rbeta

    function rexp(n, rate) result(x)
        integer, intent(in) :: n
        real(dp), intent(in), optional :: rate
        real(dp), allocatable :: x(:)
        real(dp) :: r
        integer :: i
        r = 1.0_dp
        if (present(rate)) r = rate
        allocate(x(max(n, 0)))
        do i = 1, size(x)
            x(i) = rexp_scalar(r)
        end do
    end function rexp

    function rchisq(n, df) result(x)
        integer, intent(in) :: n
        real(dp), intent(in) :: df
        real(dp), allocatable :: x(:)
        integer :: i
        allocate(x(max(n, 0)))
        do i = 1, size(x)
            x(i) = gamma_scalar(0.5_dp * df, 0.5_dp)
        end do
    end function rchisq

    function rgamma(n, shape, rate) result(x)
        integer, intent(in) :: n
        real(dp), intent(in) :: shape
        real(dp), intent(in), optional :: rate
        real(dp), allocatable :: x(:)
        real(dp) :: r
        integer :: i
        r = 1.0_dp
        if (present(rate)) r = rate
        allocate(x(max(n, 0)))
        do i = 1, size(x)
            x(i) = gamma_scalar(shape, r)
        end do
    end function rgamma

    function rgeom(n, prob) result(x)
        integer, intent(in) :: n
        real(dp), intent(in) :: prob
        real(dp), allocatable :: x(:)
        integer :: i
        allocate(x(max(n, 0)))
        do i = 1, size(x)
            x(i) = rgeom_scalar(prob)
        end do
    end function rgeom

    function rcauchy(n, location, scale) result(x)
        integer, intent(in) :: n
        real(dp), intent(in), optional :: location, scale
        real(dp), allocatable :: x(:)
        real(dp) :: loc, sc
        integer :: i
        loc = 0.0_dp
        sc = 1.0_dp
        if (present(location)) loc = location
        if (present(scale)) sc = scale
        allocate(x(max(n, 0)))
        do i = 1, size(x)
            x(i) = rcauchy_scalar(loc, sc)
        end do
    end function rcauchy

    function rt(n, df, ncp) result(x)
        integer, intent(in) :: n
        real(dp), intent(in) :: df
        real(dp), intent(in), optional :: ncp
        real(dp), allocatable :: x(:)
        real(dp) :: nc
        integer :: i
        nc = 0.0_dp
        if (present(ncp)) nc = ncp
        allocate(x(max(n, 0)))
        do i = 1, size(x)
            x(i) = rt_scalar(df, nc)
        end do
    end function rt

    function rpareto(n, shape, scale) result(x)
        integer, intent(in) :: n
        real(dp), intent(in), optional :: shape, scale
        real(dp), allocatable :: x(:)
        real(dp) :: sh, sc
        integer :: i
        sh = 1.0_dp
        sc = 1.0_dp
        if (present(shape)) sh = shape
        if (present(scale)) sc = scale
        allocate(x(max(n, 0)))
        do i = 1, size(x)
            x(i) = rpareto_scalar(sh, sc)
        end do
    end function rpareto

    function rfrechet(n, lambda, mu, sigma) result(x)
        integer, intent(in) :: n
        real(dp), intent(in), optional :: lambda, mu, sigma
        real(dp), allocatable :: x(:)
        real(dp) :: la, m, s
        integer :: i
        la = 1.0_dp
        m = 0.0_dp
        s = 1.0_dp
        if (present(lambda)) la = lambda
        if (present(mu)) m = mu
        if (present(sigma)) s = sigma
        allocate(x(max(n, 0)))
        do i = 1, size(x)
            x(i) = rfrechet_scalar(la, m, s)
        end do
    end function rfrechet

    function rlaplace(n, mu, sigma) result(x)
        integer, intent(in) :: n
        real(dp), intent(in), optional :: mu, sigma
        real(dp), allocatable :: x(:)
        real(dp) :: m, s
        integer :: i
        m = 0.0_dp
        s = 1.0_dp
        if (present(mu)) m = mu
        if (present(sigma)) s = sigma
        allocate(x(max(n, 0)))
        do i = 1, size(x)
            x(i) = rlaplace_scalar(m, s)
        end do
    end function rlaplace

    function rgumbel(n, mu, sigma) result(x)
        integer, intent(in) :: n
        real(dp), intent(in), optional :: mu, sigma
        real(dp), allocatable :: x(:)
        real(dp) :: m, s
        integer :: i
        m = 0.0_dp
        s = 1.0_dp
        if (present(mu)) m = mu
        if (present(sigma)) s = sigma
        allocate(x(max(n, 0)))
        do i = 1, size(x)
            x(i) = rgumbel_scalar(m, s)
        end do
    end function rgumbel

    function rgumble(n, mu, sigma) result(x)
        integer, intent(in) :: n
        real(dp), intent(in), optional :: mu, sigma
        real(dp), allocatable :: x(:)
        if (present(mu)) then
            if (present(sigma)) then
                x = rgumbel(n, mu, sigma)
            else
                x = rgumbel(n, mu)
            end if
        else if (present(sigma)) then
            x = rgumbel(n, sigma=sigma)
        else
            x = rgumbel(n)
        end if
    end function rgumble

    function rarcsine(n, minval, maxval) result(x)
        integer, intent(in) :: n
        real(dp), intent(in), optional :: minval, maxval
        real(dp), allocatable :: x(:)
        real(dp) :: a, b
        integer :: i
        a = 0.0_dp
        b = 1.0_dp
        if (present(minval)) a = minval
        if (present(maxval)) b = maxval
        allocate(x(max(n, 0)))
        do i = 1, size(x)
            x(i) = rarcsine_scalar(a, b)
        end do
    end function rarcsine

    function rnorm(n, mean, sd) result(x)
        integer, intent(in) :: n
        real(dp), intent(in), optional :: mean, sd
        real(dp), allocatable :: x(:)
        real(dp) :: m, s
        integer :: i
        m = 0.0_dp
        s = 1.0_dp
        if (present(mean)) m = mean
        if (present(sd)) s = sd
        allocate(x(max(n, 0)))
        do i = 1, size(x)
            x(i) = rnorm_scalar(m, s)
        end do
    end function rnorm

    function runif_mat(nrow, ncol, minval, maxval) result(x)
        integer, intent(in) :: nrow, ncol
        real(dp), intent(in), optional :: minval, maxval
        real(dp), allocatable :: x(:,:)
        real(dp), allocatable :: v(:)
        v = runif(max(nrow,0) * max(ncol,0), minval, maxval)
        allocate(x(max(nrow,0), max(ncol,0)))
        if (size(v) > 0) x = reshape(v, shape(x))
    end function runif_mat

    function rbeta_mat(nrow, ncol, alpha, beta) result(x)
        integer, intent(in) :: nrow, ncol
        real(dp), intent(in) :: alpha, beta
        real(dp), allocatable :: x(:,:)
        real(dp), allocatable :: v(:)
        v = rbeta(max(nrow,0) * max(ncol,0), alpha, beta)
        allocate(x(max(nrow,0), max(ncol,0)))
        if (size(v) > 0) x = reshape(v, shape(x))
    end function rbeta_mat

    function rexp_mat(nrow, ncol, rate) result(x)
        integer, intent(in) :: nrow, ncol
        real(dp), intent(in), optional :: rate
        real(dp), allocatable :: x(:,:)
        real(dp), allocatable :: v(:)
        v = rexp(max(nrow,0) * max(ncol,0), rate)
        allocate(x(max(nrow,0), max(ncol,0)))
        if (size(v) > 0) x = reshape(v, shape(x))
    end function rexp_mat

    function rchisq_mat(nrow, ncol, df) result(x)
        integer, intent(in) :: nrow, ncol
        real(dp), intent(in) :: df
        real(dp), allocatable :: x(:,:)
        real(dp), allocatable :: v(:)
        v = rchisq(max(nrow,0) * max(ncol,0), df)
        allocate(x(max(nrow,0), max(ncol,0)))
        if (size(v) > 0) x = reshape(v, shape(x))
    end function rchisq_mat

    function rgamma_mat(nrow, ncol, shape_par, rate) result(x)
        integer, intent(in) :: nrow, ncol
        real(dp), intent(in) :: shape_par
        real(dp), intent(in), optional :: rate
        real(dp), allocatable :: x(:,:)
        real(dp), allocatable :: v(:)
        v = rgamma(max(nrow,0) * max(ncol,0), shape_par, rate)
        allocate(x(max(nrow,0), max(ncol,0)))
        if (size(v) > 0) x = reshape(v, shape(x))
    end function rgamma_mat

    function rgeom_mat(nrow, ncol, prob) result(x)
        integer, intent(in) :: nrow, ncol
        real(dp), intent(in) :: prob
        real(dp), allocatable :: x(:,:)
        real(dp), allocatable :: v(:)
        v = rgeom(max(nrow,0) * max(ncol,0), prob)
        allocate(x(max(nrow,0), max(ncol,0)))
        if (size(v) > 0) x = reshape(v, shape(x))
    end function rgeom_mat

    function rcauchy_mat(nrow, ncol, location, scale) result(x)
        integer, intent(in) :: nrow, ncol
        real(dp), intent(in), optional :: location, scale
        real(dp), allocatable :: x(:,:)
        real(dp), allocatable :: v(:)
        v = rcauchy(max(nrow,0) * max(ncol,0), location, scale)
        allocate(x(max(nrow,0), max(ncol,0)))
        if (size(v) > 0) x = reshape(v, shape(x))
    end function rcauchy_mat

    function rt_mat(nrow, ncol, df, ncp) result(x)
        integer, intent(in) :: nrow, ncol
        real(dp), intent(in) :: df
        real(dp), intent(in), optional :: ncp
        real(dp), allocatable :: x(:,:)
        real(dp), allocatable :: v(:)
        v = rt(max(nrow,0) * max(ncol,0), df, ncp)
        allocate(x(max(nrow,0), max(ncol,0)))
        if (size(v) > 0) x = reshape(v, shape(x))
    end function rt_mat

    function rpareto_mat(nrow, ncol, shape_par, scale) result(x)
        integer, intent(in) :: nrow, ncol
        real(dp), intent(in), optional :: shape_par, scale
        real(dp), allocatable :: x(:,:)
        real(dp), allocatable :: v(:)
        v = rpareto(max(nrow,0) * max(ncol,0), shape_par, scale)
        allocate(x(max(nrow,0), max(ncol,0)))
        if (size(v) > 0) x = reshape(v, shape(x))
    end function rpareto_mat

    function rfrechet_mat(nrow, ncol, lambda, mu, sigma) result(x)
        integer, intent(in) :: nrow, ncol
        real(dp), intent(in), optional :: lambda, mu, sigma
        real(dp), allocatable :: x(:,:)
        real(dp), allocatable :: v(:)
        v = rfrechet(max(nrow,0) * max(ncol,0), lambda, mu, sigma)
        allocate(x(max(nrow,0), max(ncol,0)))
        if (size(v) > 0) x = reshape(v, shape(x))
    end function rfrechet_mat

    function rlaplace_mat(nrow, ncol, mu, sigma) result(x)
        integer, intent(in) :: nrow, ncol
        real(dp), intent(in), optional :: mu, sigma
        real(dp), allocatable :: x(:,:)
        real(dp), allocatable :: v(:)
        v = rlaplace(max(nrow,0) * max(ncol,0), mu, sigma)
        allocate(x(max(nrow,0), max(ncol,0)))
        if (size(v) > 0) x = reshape(v, shape(x))
    end function rlaplace_mat

    function rgumbel_mat(nrow, ncol, mu, sigma) result(x)
        integer, intent(in) :: nrow, ncol
        real(dp), intent(in), optional :: mu, sigma
        real(dp), allocatable :: x(:,:)
        real(dp), allocatable :: v(:)
        v = rgumbel(max(nrow,0) * max(ncol,0), mu, sigma)
        allocate(x(max(nrow,0), max(ncol,0)))
        if (size(v) > 0) x = reshape(v, shape(x))
    end function rgumbel_mat

    function rgumble_mat(nrow, ncol, mu, sigma) result(x)
        integer, intent(in) :: nrow, ncol
        real(dp), intent(in), optional :: mu, sigma
        real(dp), allocatable :: x(:,:)
        if (present(mu)) then
            if (present(sigma)) then
                x = rgumbel_mat(nrow, ncol, mu, sigma)
            else
                x = rgumbel_mat(nrow, ncol, mu)
            end if
        else if (present(sigma)) then
            x = rgumbel_mat(nrow, ncol, sigma=sigma)
        else
            x = rgumbel_mat(nrow, ncol)
        end if
    end function rgumble_mat

    function rarcsine_mat(nrow, ncol, minval, maxval) result(x)
        integer, intent(in) :: nrow, ncol
        real(dp), intent(in), optional :: minval, maxval
        real(dp), allocatable :: x(:,:)
        real(dp), allocatable :: v(:)
        v = rarcsine(max(nrow,0) * max(ncol,0), minval, maxval)
        allocate(x(max(nrow,0), max(ncol,0)))
        if (size(v) > 0) x = reshape(v, shape(x))
    end function rarcsine_mat

    function rnorm_mat(nrow, ncol, mean, sd) result(x)
        integer, intent(in) :: nrow, ncol
        real(dp), intent(in), optional :: mean, sd
        real(dp), allocatable :: x(:,:)
        real(dp), allocatable :: v(:)
        v = rnorm(max(nrow,0) * max(ncol,0), mean, sd)
        allocate(x(max(nrow,0), max(ncol,0)))
        if (size(v) > 0) x = reshape(v, shape(x))
    end function rnorm_mat

    pure function param_at(x, j) result(v)
        real(dp), intent(in) :: x(:)
        integer, intent(in) :: j
        real(dp) :: v
        v = x(1 + modulo(j - 1, size(x)))
    end function param_at

    function col_runif(nrow, ncol, minval, maxval) result(x)
        integer, intent(in) :: nrow, ncol
        real(dp), intent(in) :: minval(:), maxval(:)
        real(dp), allocatable :: x(:,:)
        integer :: j
        allocate(x(max(nrow,0), max(ncol,0)))
        do j = 1, size(x,2)
            x(:,j) = runif(size(x,1), param_at(minval,j), param_at(maxval,j))
        end do
    end function col_runif

    function col_rbeta(nrow, ncol, alpha, beta) result(x)
        integer, intent(in) :: nrow, ncol
        real(dp), intent(in) :: alpha(:), beta(:)
        real(dp), allocatable :: x(:,:)
        integer :: j
        allocate(x(max(nrow,0), max(ncol,0)))
        do j = 1, size(x,2)
            x(:,j) = rbeta(size(x,1), param_at(alpha,j), param_at(beta,j))
        end do
    end function col_rbeta

    function col_rexp(nrow, ncol, rate) result(x)
        integer, intent(in) :: nrow, ncol
        real(dp), intent(in) :: rate(:)
        real(dp), allocatable :: x(:,:)
        integer :: j
        allocate(x(max(nrow,0), max(ncol,0)))
        do j = 1, size(x,2)
            x(:,j) = rexp(size(x,1), param_at(rate,j))
        end do
    end function col_rexp

    function col_rchisq(nrow, ncol, df) result(x)
        integer, intent(in) :: nrow, ncol
        real(dp), intent(in) :: df(:)
        real(dp), allocatable :: x(:,:)
        integer :: j
        allocate(x(max(nrow,0), max(ncol,0)))
        do j = 1, size(x,2)
            x(:,j) = rchisq(size(x,1), param_at(df,j))
        end do
    end function col_rchisq

    function col_rgamma(nrow, ncol, shape, rate) result(x)
        integer, intent(in) :: nrow, ncol
        real(dp), intent(in) :: shape(:), rate(:)
        real(dp), allocatable :: x(:,:)
        integer :: j
        allocate(x(max(nrow,0), max(ncol,0)))
        do j = 1, size(x,2)
            x(:,j) = rgamma(size(x,1), param_at(shape,j), param_at(rate,j))
        end do
    end function col_rgamma

    function col_rgeom(nrow, ncol, prob) result(x)
        integer, intent(in) :: nrow, ncol
        real(dp), intent(in) :: prob(:)
        real(dp), allocatable :: x(:,:)
        integer :: j
        allocate(x(max(nrow,0), max(ncol,0)))
        do j = 1, size(x,2)
            x(:,j) = rgeom(size(x,1), param_at(prob,j))
        end do
    end function col_rgeom

    function col_rcauchy(nrow, ncol, location, scale) result(x)
        integer, intent(in) :: nrow, ncol
        real(dp), intent(in) :: location(:), scale(:)
        real(dp), allocatable :: x(:,:)
        integer :: j
        allocate(x(max(nrow,0), max(ncol,0)))
        do j = 1, size(x,2)
            x(:,j) = rcauchy(size(x,1), param_at(location,j), param_at(scale,j))
        end do
    end function col_rcauchy

    function col_rt(nrow, ncol, df, ncp) result(x)
        integer, intent(in) :: nrow, ncol
        real(dp), intent(in) :: df(:), ncp(:)
        real(dp), allocatable :: x(:,:)
        integer :: j
        allocate(x(max(nrow,0), max(ncol,0)))
        do j = 1, size(x,2)
            x(:,j) = rt(size(x,1), param_at(df,j), param_at(ncp,j))
        end do
    end function col_rt

    function col_rpareto(nrow, ncol, shape, scale) result(x)
        integer, intent(in) :: nrow, ncol
        real(dp), intent(in) :: shape(:), scale(:)
        real(dp), allocatable :: x(:,:)
        integer :: j
        allocate(x(max(nrow,0), max(ncol,0)))
        do j = 1, size(x,2)
            x(:,j) = rpareto(size(x,1), param_at(shape,j), param_at(scale,j))
        end do
    end function col_rpareto

    function col_rfrechet(nrow, ncol, lambda, mu, sigma) result(x)
        integer, intent(in) :: nrow, ncol
        real(dp), intent(in) :: lambda(:), mu(:), sigma(:)
        real(dp), allocatable :: x(:,:)
        integer :: j
        allocate(x(max(nrow,0), max(ncol,0)))
        do j = 1, size(x,2)
            x(:,j) = rfrechet(size(x,1), param_at(lambda,j), param_at(mu,j), param_at(sigma,j))
        end do
    end function col_rfrechet

    function col_rlaplace(nrow, ncol, mu, sigma) result(x)
        integer, intent(in) :: nrow, ncol
        real(dp), intent(in) :: mu(:), sigma(:)
        real(dp), allocatable :: x(:,:)
        integer :: j
        allocate(x(max(nrow,0), max(ncol,0)))
        do j = 1, size(x,2)
            x(:,j) = rlaplace(size(x,1), param_at(mu,j), param_at(sigma,j))
        end do
    end function col_rlaplace

    function col_rgumbel(nrow, ncol, mu, sigma) result(x)
        integer, intent(in) :: nrow, ncol
        real(dp), intent(in) :: mu(:), sigma(:)
        real(dp), allocatable :: x(:,:)
        integer :: j
        allocate(x(max(nrow,0), max(ncol,0)))
        do j = 1, size(x,2)
            x(:,j) = rgumbel(size(x,1), param_at(mu,j), param_at(sigma,j))
        end do
    end function col_rgumbel

    function col_rgumble(nrow, ncol, mu, sigma) result(x)
        integer, intent(in) :: nrow, ncol
        real(dp), intent(in) :: mu(:), sigma(:)
        real(dp), allocatable :: x(:,:)
        x = col_rgumbel(nrow, ncol, mu, sigma)
    end function col_rgumble

    function col_rarcsine(nrow, ncol, minval, maxval) result(x)
        integer, intent(in) :: nrow, ncol
        real(dp), intent(in) :: minval(:), maxval(:)
        real(dp), allocatable :: x(:,:)
        integer :: j
        allocate(x(max(nrow,0), max(ncol,0)))
        do j = 1, size(x,2)
            x(:,j) = rarcsine(size(x,1), param_at(minval,j), param_at(maxval,j))
        end do
    end function col_rarcsine

    function col_rnorm(nrow, ncol, mean, sd) result(x)
        integer, intent(in) :: nrow, ncol
        real(dp), intent(in) :: mean(:), sd(:)
        real(dp), allocatable :: x(:,:)
        integer :: j
        allocate(x(max(nrow,0), max(ncol,0)))
        do j = 1, size(x,2)
            x(:,j) = rnorm(size(x,1), param_at(mean,j), param_at(sd,j))
        end do
    end function col_rnorm

end module rangen_distributions
