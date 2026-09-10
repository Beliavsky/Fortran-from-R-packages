! SPDX-License-Identifier: GPL-2.0-or-later
module dlm_models
    use dlm_types, only : dp, dlm_model, dlm_success, dlm_invalid_shape, dlm_invalid_argument
    implicit none
    private

    public :: ar_trans_pars
    public :: block_diag2
    public :: dlm_add
    public :: dlm_mod_arma
    public :: dlm_mod_poly
    public :: dlm_mod_reg
    public :: dlm_mod_seas
    public :: dlm_mod_trig
    public :: dlm_sum

contains

    pure subroutine block_diag2(a, b, c)
        real(dp), intent(in) :: a(:, :) !! First matrix placed in the upper-left block.
        real(dp), intent(in) :: b(:, :) !! Second matrix placed in the lower-right block.
        real(dp), allocatable, intent(out) :: c(:, :) !! Allocated block-diagonal matrix containing `a` and `b`.

        integer :: ra, ca, rb, cb

        ra = size(a, 1)
        ca = size(a, 2)
        rb = size(b, 1)
        cb = size(b, 2)
        allocate(c(ra + rb, ca + cb))
        c = 0.0_dp
        if (ra > 0 .and. ca > 0) c(1:ra, 1:ca) = a
        if (rb > 0 .and. cb > 0) c(ra + 1:ra + rb, ca + 1:ca + cb) = b
    end subroutine block_diag2

    pure subroutine ar_trans_pars(raw, coefficients)
        real(dp), intent(in) :: raw(:) !! Unconstrained real parameters mapped to a stationary AR coefficient vector.
        real(dp), intent(out) :: coefficients(size(raw)) !! Stationary autoregression coefficients in lag order.

        real(dp) :: work(size(raw)), reflection
        integer :: j, k

        coefficients = tanh(raw)
        work = coefficients
        do j = 2, size(raw)
            reflection = coefficients(j)
            do k = 1, j - 1
                work(k) = work(k) - reflection * coefficients(j - k)
            end do
            coefficients(1:j - 1) = work(1:j - 1)
        end do
    end subroutine ar_trans_pars

    pure subroutine dlm_mod_poly(order, model, info, d_v, d_w, m0, c0)
        integer, intent(in) :: order !! Polynomial trend order and resulting state dimension; must be positive.
        type(dlm_model), intent(out) :: model !! Constructed polynomial-trend DLM.
        integer, intent(out) :: info !! Zero on success or `dlm_invalid_argument`/`dlm_invalid_shape`.
        real(dp), intent(in), optional :: d_v !! Scalar observation variance; defaults to one.
        real(dp), intent(in), optional :: d_w(:) !! Innovation variances; length `order`, default `[0,...,0,1]`.
        real(dp), intent(in), optional :: m0(:) !! Initial state mean; length `order`, default all zero.
        real(dp), intent(in), optional :: c0(:, :) !! Initial covariance; shape `(order,order)`, default `1e7 I`.

        integer :: i

        info = dlm_success
        if (order < 1) then
            info = dlm_invalid_argument
            return
        end if
        if (present(d_w)) then
            if (size(d_w) /= order) then
                info = dlm_invalid_shape
                return
            end if
        end if
        if (present(m0)) then
            if (size(m0) /= order) then
                info = dlm_invalid_shape
                return
            end if
        end if
        if (present(c0)) then
            if (size(c0, 1) /= order .or. size(c0, 2) /= order) then
                info = dlm_invalid_shape
                return
            end if
        end if

        allocate(model%m0(order), model%c0(order, order), model%ff(1, order), model%v(1, 1))
        allocate(model%gg(order, order), model%w(order, order))
        model%m0 = 0.0_dp
        if (present(m0)) model%m0 = m0
        model%c0 = 0.0_dp
        do i = 1, order
            model%c0(i, i) = 1.0e7_dp
        end do
        if (present(c0)) model%c0 = c0
        model%ff = 0.0_dp
        model%ff(1, 1) = 1.0_dp
        model%v(1, 1) = 1.0_dp
        if (present(d_v)) model%v(1, 1) = d_v
        model%gg = 0.0_dp
        do i = 1, order
            model%gg(i, i) = 1.0_dp
            if (i < order) model%gg(i, i + 1) = 1.0_dp
        end do
        model%w = 0.0_dp
        if (present(d_w)) then
            do i = 1, order
                model%w(i, i) = d_w(i)
            end do
        else
            model%w(order, order) = 1.0_dp
        end if
    end subroutine dlm_mod_poly

    pure subroutine dlm_mod_seas(frequency, model, info, d_v, d_w, m0, c0)
        integer, intent(in) :: frequency !! Seasonal frequency; must be at least two, producing `frequency-1` states.
        type(dlm_model), intent(out) :: model !! Constructed seasonal-factor DLM.
        integer, intent(out) :: info !! Zero on success or a `dlm_*` status code.
        real(dp), intent(in), optional :: d_v !! Scalar observation variance; defaults to one.
        real(dp), intent(in), optional :: d_w(:) !! Innovation variances of length `frequency-1`; default `(1,0,...)`.
        real(dp), intent(in), optional :: m0(:) !! Initial state mean of length `frequency-1`; defaults to zero.
        real(dp), intent(in), optional :: c0(:, :) !! Initial covariance of shape `(p,p)`; defaults to `1e7 I`.

        integer :: i, p

        info = dlm_success
        if (frequency < 2) then
            info = dlm_invalid_argument
            return
        end if
        p = frequency - 1
        if (present(d_w)) then
            if (size(d_w) /= p) then
                info = dlm_invalid_shape
                return
            end if
        end if
        if (present(m0)) then
            if (size(m0) /= p) then
                info = dlm_invalid_shape
                return
            end if
        end if
        if (present(c0)) then
            if (size(c0, 1) /= p .or. size(c0, 2) /= p) then
                info = dlm_invalid_shape
                return
            end if
        end if

        allocate(model%m0(p), model%c0(p, p), model%ff(1, p), model%v(1, 1))
        allocate(model%gg(p, p), model%w(p, p))
        model%m0 = 0.0_dp
        if (present(m0)) model%m0 = m0
        model%c0 = 0.0_dp
        do i = 1, p
            model%c0(i, i) = 1.0e7_dp
        end do
        if (present(c0)) model%c0 = c0
        model%ff = 0.0_dp
        model%ff(1, 1) = 1.0_dp
        model%v(1, 1) = 1.0_dp
        if (present(d_v)) model%v(1, 1) = d_v
        model%gg = 0.0_dp
        model%gg(1, :) = -1.0_dp
        do i = 2, p
            model%gg(i, i - 1) = 1.0_dp
        end do
        model%w = 0.0_dp
        if (present(d_w)) then
            do i = 1, p
                model%w(i, i) = d_w(i)
            end do
        else
            model%w(1, 1) = 1.0_dp
        end if
    end subroutine dlm_mod_seas

    pure subroutine dlm_mod_reg(x, model, info, add_intercept, d_v, d_w, m0, c0)
        real(dp), intent(in) :: x(:, :) !! Regression design matrix with observations in rows and predictors in columns.
        type(dlm_model), intent(out) :: model !! Constructed time-varying regression DLM using `JFF` and `X`.
        integer, intent(out) :: info !! Zero on success or a `dlm_*` status code.
        logical, intent(in), optional :: add_intercept !! Include a constant state when true; defaults to true.
        real(dp), intent(in), optional :: d_v !! Scalar observation variance; defaults to one.
        real(dp), intent(in), optional :: d_w(:) !! State innovation variances of length number of states; defaults to zero.
        real(dp), intent(in), optional :: m0(:) !! Initial coefficient mean; length number of states, defaults to zero.
        real(dp), intent(in), optional :: c0(:, :) !! Initial coefficient covariance; defaults to `1e7 I`.

        logical :: intercept
        integer :: i, j, p

        intercept = .true.
        if (present(add_intercept)) intercept = add_intercept
        p = size(x, 2)
        if (intercept) p = p + 1
        info = dlm_success
        if (p < 1) then
            info = dlm_invalid_argument
            return
        end if
        if (present(d_w)) then
            if (size(d_w) /= p) then
                info = dlm_invalid_shape
                return
            end if
        end if
        if (present(m0)) then
            if (size(m0) /= p) then
                info = dlm_invalid_shape
                return
            end if
        end if
        if (present(c0)) then
            if (size(c0, 1) /= p .or. size(c0, 2) /= p) then
                info = dlm_invalid_shape
                return
            end if
        end if

        allocate(model%m0(p), model%c0(p, p), model%ff(1, p), model%v(1, 1))
        allocate(model%gg(p, p), model%w(p, p), model%jff(1, p))
        model%x = x
        model%m0 = 0.0_dp
        if (present(m0)) model%m0 = m0
        model%c0 = 0.0_dp
        do i = 1, p
            model%c0(i, i) = 1.0e7_dp
        end do
        if (present(c0)) model%c0 = c0
        model%ff = 1.0_dp
        model%v(1, 1) = 1.0_dp
        if (present(d_v)) model%v(1, 1) = d_v
        model%gg = 0.0_dp
        do i = 1, p
            model%gg(i, i) = 1.0_dp
        end do
        model%w = 0.0_dp
        if (present(d_w)) then
            do i = 1, p
                model%w(i, i) = d_w(i)
            end do
        end if
        model%jff = 0
        if (intercept) then
            do j = 2, p
                model%jff(1, j) = j - 1
            end do
        else
            do j = 1, p
                model%jff(1, j) = j
            end do
        end if
    end subroutine dlm_mod_reg

    pure subroutine dlm_mod_trig(model, info, period, q, omega, tau, d_v, d_w, m0, c0)
        type(dlm_model), intent(out) :: model !! Constructed Fourier/trigonometric DLM.
        integer, intent(out) :: info !! Zero on success or a `dlm_*` status code.
        integer, intent(in), optional :: period !! Integer seasonal period corresponding to R argument `s`.
        integer, intent(in), optional :: q !! Number of harmonics; required with `omega`/`tau`, optional with `period`.
        real(dp), intent(in), optional :: omega !! Base angular frequency in radians per step; mutually exclusive with `period`.
        real(dp), intent(in), optional :: tau !! Period in steps converted to `2*pi/tau`; mutually exclusive with `period`.
        real(dp), intent(in), optional :: d_v !! Scalar observation variance; defaults to one.
        real(dp), intent(in), optional :: d_w(:) !! Innovation variance vector of length one or state dimension; defaults to zero.
        real(dp), intent(in), optional :: m0(:) !! Initial state mean with one entry per trigonometric state; defaults to zero.
        real(dp), intent(in), optional :: c0(:, :) !! Initial covariance with state-dimension square shape; defaults to `1e7 I`.

        real(dp), parameter :: pi = acos(-1.0_dp)
        real(dp) :: angle, h(2, 2), hp(2, 2)
        logical :: even_all
        integer :: block, harmonic_count, i, p, pos, s_half

        info = dlm_success
        angle = 0.0_dp
        if (present(period)) then
            if (present(omega) .or. present(tau) .or. period < 2) then
                info = dlm_invalid_argument
                return
            end if
            s_half = period / 2
            harmonic_count = s_half
            if (present(q)) harmonic_count = q
            if (harmonic_count < 1 .or. harmonic_count > s_half) then
                info = dlm_invalid_argument
                return
            end if
            angle = 2.0_dp * pi / real(period, dp)
            even_all = mod(period, 2) == 0 .and. harmonic_count == s_half
        else
            if (present(omega) .eqv. present(tau)) then
                info = dlm_invalid_argument
                return
            end if
            if (.not. present(q)) then
                info = dlm_invalid_argument
                return
            end if
            harmonic_count = q
            if (harmonic_count < 1) then
                info = dlm_invalid_argument
                return
            end if
            if (present(omega)) angle = omega
            if (present(tau)) then
                if (tau <= 0.0_dp) then
                    info = dlm_invalid_argument
                    return
                end if
                angle = 2.0_dp * pi / tau
            end if
            even_all = .false.
        end if

        p = 2 * harmonic_count
        if (even_all) p = p - 1
        if (present(d_w)) then
            if (size(d_w) /= 1 .and. size(d_w) /= p) then
                info = dlm_invalid_shape
                return
            end if
        end if
        if (present(m0)) then
            if (size(m0) /= p) then
                info = dlm_invalid_shape
                return
            end if
        end if
        if (present(c0)) then
            if (size(c0, 1) /= p .or. size(c0, 2) /= p) then
                info = dlm_invalid_shape
                return
            end if
        end if

        allocate(model%m0(p), model%c0(p, p), model%ff(1, p), model%v(1, 1))
        allocate(model%gg(p, p), model%w(p, p))
        model%m0 = 0.0_dp
        if (present(m0)) model%m0 = m0
        model%c0 = 0.0_dp
        do i = 1, p
            model%c0(i, i) = 1.0e7_dp
        end do
        if (present(c0)) model%c0 = c0
        model%ff = 0.0_dp
        model%v(1, 1) = 1.0_dp
        if (present(d_v)) model%v(1, 1) = d_v
        model%gg = 0.0_dp
        model%w = 0.0_dp

        h = 0.0_dp
        h(1, 1) = cos(angle)
        h(2, 2) = h(1, 1)
        h(1, 2) = sin(angle)
        h(2, 1) = -h(1, 2)
        hp = h
        pos = 1
        do block = 1, harmonic_count
            if (even_all .and. block == harmonic_count) then
                model%gg(pos, pos) = -1.0_dp
                model%ff(1, pos) = 1.0_dp
                pos = pos + 1
            else
                model%gg(pos:pos + 1, pos:pos + 1) = hp
                model%ff(1, pos) = 1.0_dp
                pos = pos + 2
                hp = matmul(hp, h)
            end if
        end do
        if (present(d_w)) then
            if (size(d_w) == 1) then
                do i = 1, p
                    model%w(i, i) = d_w(1)
                end do
            else
                do i = 1, p
                    model%w(i, i) = d_w(i)
                end do
            end if
        end if
    end subroutine dlm_mod_trig

    pure subroutine dlm_mod_arma(model, info, ar, ma, sigma2, d_v, m0, c0)
        type(dlm_model), intent(out) :: model !! Constructed univariate ARMA state-space DLM.
        integer, intent(out) :: info !! Zero on success or a `dlm_*` status code.
        real(dp), intent(in), optional :: ar(:) !! Autoregressive coefficients in lag order; omitted means AR order zero.
        real(dp), intent(in), optional :: ma(:) !! Moving-average coefficients in lag order; omitted means MA order zero.
        real(dp), intent(in), optional :: sigma2 !! Innovation variance; defaults to one.
        real(dp), intent(in), optional :: d_v !! Additional observation variance; defaults to zero.
        real(dp), intent(in), optional :: m0(:) !! Initial state mean of length `max(p,q+1)`; defaults to zero.
        real(dp), intent(in), optional :: c0(:, :) !! Initial covariance of matching state dimension; defaults to `1e7 I`.

        real(dp), allocatable :: innovation_loading(:)
        real(dp) :: variance
        integer :: i, p_ar, p_ma, r

        p_ar = 0
        p_ma = 0
        if (present(ar)) p_ar = size(ar)
        if (present(ma)) p_ma = size(ma)
        r = max(p_ar, p_ma + 1)
        if (r < 1) r = 1
        info = dlm_success
        if (present(m0)) then
            if (size(m0) /= r) then
                info = dlm_invalid_shape
                return
            end if
        end if
        if (present(c0)) then
            if (size(c0, 1) /= r .or. size(c0, 2) /= r) then
                info = dlm_invalid_shape
                return
            end if
        end if

        allocate(model%m0(r), model%c0(r, r), model%ff(1, r), model%v(1, 1))
        allocate(model%gg(r, r), model%w(r, r), innovation_loading(r))
        model%m0 = 0.0_dp
        if (present(m0)) model%m0 = m0
        model%c0 = 0.0_dp
        do i = 1, r
            model%c0(i, i) = 1.0e7_dp
        end do
        if (present(c0)) model%c0 = c0
        model%ff = 0.0_dp
        model%ff(1, 1) = 1.0_dp
        model%v(1, 1) = 0.0_dp
        if (present(d_v)) model%v(1, 1) = d_v
        model%gg = 0.0_dp
        if (present(ar)) model%gg(1:p_ar, 1) = ar
        do i = 1, r - 1
            model%gg(i, i + 1) = 1.0_dp
        end do
        innovation_loading = 0.0_dp
        innovation_loading(1) = 1.0_dp
        if (present(ma)) innovation_loading(2:p_ma + 1) = ma
        variance = 1.0_dp
        if (present(sigma2)) variance = sigma2
        model%w = variance * spread(innovation_loading, 2, r) * spread(innovation_loading, 1, r)
    end subroutine dlm_mod_arma

    pure subroutine dlm_sum(left, right, model, info)
        type(dlm_model), intent(in) :: left !! First constant DLM whose observation vector is stacked with `right`.
        type(dlm_model), intent(in) :: right !! Second constant DLM whose observation vector is stacked with `left`.
        type(dlm_model), intent(out) :: model !! Outer-sum DLM with block-diagonal observation and state matrices.
        integer, intent(out) :: info !! Zero on success; time-varying inputs return `dlm_invalid_argument`.

        integer :: p1, p2

        info = dlm_success
        if (allocated(left%x) .or. allocated(right%x)) then
            info = dlm_invalid_argument
            return
        end if
        p1 = size(left%ff, 2)
        p2 = size(right%ff, 2)
        allocate(model%m0(p1 + p2))
        model%m0 = [left%m0, right%m0]
        call block_diag2(left%c0, right%c0, model%c0)
        call block_diag2(left%ff, right%ff, model%ff)
        call block_diag2(left%v, right%v, model%v)
        call block_diag2(left%gg, right%gg, model%gg)
        call block_diag2(left%w, right%w, model%w)
    end subroutine dlm_sum

    pure subroutine dlm_add(left, right, model, info)
        type(dlm_model), intent(in) :: left !! First constant DLM contributing additively to the same observation vector.
        type(dlm_model), intent(in) :: right !! Second constant DLM with the same observation dimension as `left`.
        type(dlm_model), intent(out) :: model !! Additive latent-state DLM corresponding to R's `+.dlm` constant-model case.
        integer, intent(out) :: info !! Zero on success or `dlm_invalid_shape`/`dlm_invalid_argument`.

        integer :: m_obs, p1, p2

        info = dlm_success
        if (allocated(left%x) .or. allocated(right%x)) then
            info = dlm_invalid_argument
            return
        end if
        m_obs = size(left%ff, 1)
        if (size(right%ff, 1) /= m_obs) then
            info = dlm_invalid_shape
            return
        end if
        p1 = size(left%ff, 2)
        p2 = size(right%ff, 2)
        allocate(model%m0(p1 + p2), model%ff(m_obs, p1 + p2), model%v(m_obs, m_obs))
        model%m0 = [left%m0, right%m0]
        model%ff(:, 1:p1) = left%ff
        model%ff(:, p1 + 1:p1 + p2) = right%ff
        model%v = left%v + right%v
        call block_diag2(left%c0, right%c0, model%c0)
        call block_diag2(left%gg, right%gg, model%gg)
        call block_diag2(left%w, right%w, model%w)
    end subroutine dlm_add

end module dlm_models
