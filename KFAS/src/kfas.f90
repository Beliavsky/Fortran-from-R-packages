! SPDX-License-Identifier: GPL-2.0-or-later
!
! Modern public Fortran interface for the computational kernels derived from
! KFAS 1.6.0 by Jouni Helske. See NOTICE.md and UPSTREAM.md.

module kfas
    use kfas_kinds, only: dp
    implicit none
    private

    public :: dp
    integer, parameter, public :: kfas_gaussian = 1
    integer, parameter, public :: kfas_poisson = 2
    integer, parameter, public :: kfas_binomial = 3
    integer, parameter, public :: kfas_gamma = 4
    integer, parameter, public :: kfas_negative_binomial = 5

    type, public :: kfas_model
        real(dp), allocatable :: y(:, :)
        integer, allocatable :: missing(:, :)
        real(dp), allocatable :: z(:, :, :)
        real(dp), allocatable :: h(:, :, :)
        real(dp), allocatable :: tmat(:, :, :)
        real(dp), allocatable :: rmat(:, :, :)
        real(dp), allocatable :: q(:, :, :)
        real(dp), allocatable :: a1(:)
        real(dp), allocatable :: p1(:, :)
        real(dp), allocatable :: p1inf(:, :)
        integer :: time_varying(5) = 0
        integer :: diffuse_rank = 0
    contains
        procedure :: validate => validate_model
    end type kfas_model

    type, public :: kfas_filter_result
        real(dp), allocatable :: a(:, :)
        real(dp), allocatable :: p(:, :, :)
        real(dp), allocatable :: pinf(:, :, :)
        real(dp), allocatable :: att(:, :)
        real(dp), allocatable :: ptt(:, :, :)
        real(dp), allocatable :: v(:, :)
        real(dp), allocatable :: f(:, :)
        real(dp), allocatable :: finf(:, :)
        real(dp), allocatable :: k(:, :, :)
        real(dp), allocatable :: kinf(:, :, :)
        real(dp), allocatable :: signal(:, :)
        real(dp), allocatable :: signal_var(:, :, :)
        real(dp) :: loglik = 0.0_dp
        integer :: diffuse_end = 0
        integer :: diffuse_component = 0
        integer :: remaining_diffuse_rank = 0
        logical :: transformed_h = .false.
    end type kfas_filter_result

    type, public :: kfas_approx_result
        real(dp), allocatable :: theta(:, :)
        real(dp), allocatable :: ytilde(:, :)
        real(dp), allocatable :: htilde(:, :, :)
        real(dp) :: gaussian_loglik = 0.0_dp
        real(dp) :: difference = 0.0_dp
        real(dp) :: h_tol = 1.0e15_dp
        integer :: iterations = 0
        integer :: info = 0
    end type kfas_approx_result

    type, public :: kfas_smooth_result
        real(dp), allocatable :: state(:, :)
        real(dp), allocatable :: state_var(:, :, :)
        real(dp), allocatable :: obs_disturbance(:, :)
        real(dp), allocatable :: obs_disturbance_var(:, :)
        real(dp), allocatable :: state_disturbance(:, :)
        real(dp), allocatable :: state_disturbance_var(:, :, :)
        real(dp), allocatable :: signal(:, :)
        real(dp), allocatable :: signal_var(:, :, :)
        real(dp), allocatable :: rvec(:, :)
        real(dp), allocatable :: nmat(:, :, :)
        real(dp) :: loglik = 0.0_dp
        integer :: diffuse_end = 0
        logical :: transformed_h = .false.
        logical :: observation_disturbance_transformed = .false.
    end type kfas_smooth_result

    public :: kfas_approximate_nongaussian
    public :: kfas_ar_transform
    public :: kfas_gaussian_filter
    public :: kfas_gaussian_loglik
    public :: kfas_gaussian_smooth
    public :: kfas_init_theta
    public :: kfas_ldl_factor
    public :: kfas_ldl_transform
    public :: kfas_nongaussian_loglik
    public :: kfas_weighted_mean_cov

    interface
        subroutine artransform(phi, p)
            import :: dp
            integer, intent(in) :: p
            real(dp), intent(inout) :: phi(p)
        end subroutine artransform

        subroutine covmeanwprotect(x, w, m, n, k, meanx, covx)
            import :: dp
            integer, intent(in) :: m, n, k
            real(dp), intent(in) :: x(m, n, k)
            real(dp), intent(in) :: w(k)
            real(dp), intent(inout) :: meanx(m, n)
            real(dp), intent(inout) :: covx(m, m, n)
        end subroutine covmeanwprotect

        subroutine kfilter2(yt, ymiss, timevar, zt, ht, tt, rt, qt, a1, p1, p1inf, &
                p, n, m, r, d, j, at, pt, vt, ft, kt, pinf, finf, kinf, lik, tol, &
                rankp, theta, thetavar, filtersignal, att, ptt)
            import :: dp
            integer, intent(in) :: p, n, m, r, filtersignal
            integer, intent(inout) :: d, j, rankp
            integer, intent(in) :: ymiss(n, p)
            integer, intent(in) :: timevar(5)
            real(dp), intent(in) :: yt(n, p)
            real(dp), intent(in) :: zt(p, m, (n - 1) * timevar(1) + 1)
            real(dp), intent(in) :: ht(p, p, (n - 1) * timevar(2) + 1)
            real(dp), intent(in) :: tt(m, m, (n - 1) * timevar(3) + 1)
            real(dp), intent(in) :: rt(m, r, (n - 1) * timevar(4) + 1)
            real(dp), intent(in) :: qt(r, r, (n - 1) * timevar(5) + 1)
            real(dp), intent(in) :: a1(m)
            real(dp), intent(in) :: p1(m, m), p1inf(m, m)
            real(dp), intent(inout) :: at(m, n + 1)
            real(dp), intent(inout) :: pt(m, m, n + 1), pinf(m, m, n + 1)
            real(dp), intent(inout) :: vt(p, n), ft(p, n), finf(p, n)
            real(dp), intent(inout) :: kt(m, p, n), kinf(m, p, n)
            real(dp), intent(inout) :: lik
            real(dp), intent(in) :: tol
            real(dp), intent(inout) :: theta(n, p)
            real(dp), intent(inout) :: thetavar(p, p, n)
            real(dp), intent(inout) :: att(m, n), ptt(m, m, n)
        end subroutine kfilter2

        subroutine ldl(a, n, tol, info)
            import :: dp
            integer, intent(in) :: n
            integer, intent(inout) :: info
            real(dp), intent(inout) :: a(n, n)
            real(dp), intent(in) :: tol
        end subroutine ldl

        subroutine ldlssm(yt, ydimt, yobs, timevar, zt, p, m, n, ichols, nh, &
                hchol, dim, info, hobs, tol)
            import :: dp
            integer, intent(in) :: p, m, n, nh
            integer, intent(in) :: ydimt(n), yobs(p, n), timevar(5), dim(nh), hobs(p, nh)
            integer, intent(inout) :: hchol(n), info
            real(dp), intent(inout) :: yt(p, n)
            real(dp), intent(inout) :: zt(p, m, (n - 1) * timevar(1) + 1)
            real(dp), intent(inout) :: ichols(p, p, nh)
            real(dp), intent(in) :: tol
        end subroutine ldlssm

        subroutine approx(yt, ymiss, timevar, zt, tt, rtv, ht, qt, a1, p1, p1inf, &
                p, n, m, r, theta, u, ytilde, dist, maxiter, tol, rankp, convtol, &
                diff, lik, info, expected, htol)
            import :: dp
            integer, intent(in) :: p, m, r, n, rankp, expected
            integer, intent(in) :: ymiss(n, p), timevar(5), dist(p)
            integer, intent(inout) :: maxiter, info
            real(dp), intent(in) :: tol, convtol, u(n, p), yt(n, p)
            real(dp), intent(in) :: zt(p, m, (n - 1) * timevar(1) + 1)
            real(dp), intent(in) :: tt(m, m, (n - 1) * timevar(3) + 1)
            real(dp), intent(in) :: rtv(m, r, (n - 1) * timevar(4) + 1)
            real(dp), intent(in) :: qt(r, r, (n - 1) * timevar(5) + 1)
            real(dp), intent(in) :: a1(m), p1(m, m), p1inf(m, m)
            real(dp), intent(inout) :: theta(n, p), ytilde(n, p), ht(p, p, n)
            real(dp), intent(inout) :: diff, lik, htol
        end subroutine approx

        subroutine ngloglik(yt, ymiss, timevar, zt, tt, rtv, qt, a1, p1, p1inf, &
                p, m, r, n, lik, theta, u, dist, maxiter, rankp, convtol, nnd, nsim, &
                epsplus, etaplus, aplus1, c, tol, info, antit, sim, nsim2, diff, &
                marginal, expected, htol)
            import :: dp
            integer, intent(in) :: p, m, r, n, nnd, antit, nsim, sim, nsim2, rankp, expected
            integer, intent(in) :: ymiss(n, p), dist(p), timevar(5)
            integer, intent(inout) :: maxiter, marginal, info
            real(dp), intent(in) :: convtol, tol, u(n, p), yt(n, p)
            real(dp), intent(in) :: zt(p, m, (n - 1) * timevar(1) + 1)
            real(dp), intent(in) :: tt(m, m, (n - 1) * timevar(3) + 1)
            real(dp), intent(in) :: rtv(m, r, (n - 1) * timevar(4) + 1)
            real(dp), intent(in) :: qt(r, r, (n - 1) * timevar(5) + 1)
            real(dp), intent(in) :: a1(m), p1(m, m), p1inf(m, m), c(nsim)
            real(dp), intent(inout) :: aplus1(m, nsim), epsplus(p, n, nsim)
            real(dp), intent(inout) :: etaplus(r, n, nsim), theta(n, p)
            real(dp), intent(inout) :: lik, diff, htol
        end subroutine ngloglik

        subroutine gsmoothall(ymiss, timevar, zt, ht, tt, rtv, qt, p, n, m, r, &
                d, j, at, pt, vt, ft, kt, rt, rt0, rt1, nt, nt0, nt1, nt2, pinf, &
                kinf, finf, ahat, vvt, epshat, epshatvar, etahat, etahatvar, &
                thetahat, thetahatvar, ldlsignal, zorig, zorigtv, aug, state, dist, signal)
            import :: dp
            integer, intent(in) :: d, j, p, r, m, n, aug, state, dist, signal
            integer, intent(in) :: ldlsignal, zorigtv
            integer, intent(in) :: ymiss(n, p), timevar(5)
            real(dp), intent(in) :: zt(p, m, (n - 1) * timevar(1) + 1)
            real(dp), intent(in) :: ht(p, p, (n - 1) * timevar(2) + 1)
            real(dp), intent(in) :: tt(m, m, (n - 1) * timevar(3) + 1)
            real(dp), intent(in) :: rtv(m, r, (n - 1) * timevar(4) + 1)
            real(dp), intent(in) :: qt(r, r, (n - 1) * timevar(5) + 1)
            real(dp), intent(in) :: at(m, n + 1), pt(m, m, n + 1)
            real(dp), intent(in) :: vt(p, n), ft(p, n), kt(m, p, n)
            real(dp), intent(in) :: pinf(m, m, d + 1), kinf(m, p, d), finf(p, d)
            real(dp), intent(inout) :: rt(m, n + 1), rt0(m, d + 1), rt1(m, d + 1)
            real(dp), intent(inout) :: nt(m, m, n + 1), nt0(m, m, d + 1)
            real(dp), intent(inout) :: nt1(m, m, d + 1), nt2(m, m, d + 1)
            real(dp), intent(inout) :: ahat(m * state, n * state)
            real(dp), intent(inout) :: vvt(m, m, n)
            real(dp), intent(inout) :: epshat(p * dist * aug, n * dist * aug)
            real(dp), intent(inout) :: epshatvar(p * dist * aug, n * dist * aug)
            real(dp), intent(inout) :: etahat(r * dist, n * dist)
            real(dp), intent(inout) :: etahatvar(r * dist, r * dist, n * dist)
            real(dp), intent(inout) :: thetahat(p * signal, n * signal)
            real(dp), intent(inout) :: thetahatvar(p * signal, p * signal, n * signal)
            real(dp), intent(in) :: zorig(ldlsignal * p, ldlsignal * m, &
                ldlsignal * ((n - 1) * zorigtv + 1))
        end subroutine gsmoothall
    end interface

contains

    subroutine validate_model(self, ok, message)
        class(kfas_model), intent(in) :: self
        logical, intent(out) :: ok
        character(len = :), allocatable, intent(out) :: message
        integer :: m, n, p, r

        ok = .false.
        message = ""

        if (.not. allocated(self%y)) then
            message = "y is not allocated"
            return
        end if
        if (.not. allocated(self%missing)) then
            message = "missing is not allocated"
            return
        end if
        if (.not. allocated(self%z) .or. .not. allocated(self%h)) then
            message = "z and h must be allocated"
            return
        end if
        if (.not. allocated(self%tmat) .or. .not. allocated(self%rmat) .or. &
                .not. allocated(self%q)) then
            message = "tmat, rmat, and q must be allocated"
            return
        end if
        if (.not. allocated(self%a1) .or. .not. allocated(self%p1) .or. &
                .not. allocated(self%p1inf)) then
            message = "a1, p1, and p1inf must be allocated"
            return
        end if
        if (any(self%time_varying < 0) .or. any(self%time_varying > 1)) then
            message = "time_varying entries must be 0 or 1"
            return
        end if

        n = size(self%y, 1)
        p = size(self%y, 2)
        m = size(self%a1)
        r = size(self%q, 1)

        if (n < 1 .or. p < 1 .or. m < 1 .or. r < 1) then
            message = "model dimensions must be positive"
            return
        end if
        if (size(self%missing, 1) /= n .or. size(self%missing, 2) /= p) then
            message = "missing must have shape (n, p)"
            return
        end if
        if (any(self%missing /= 0 .and. self%missing /= 1)) then
            message = "missing entries must be 0 or 1"
            return
        end if
        if (size(self%z, 1) /= p .or. size(self%z, 2) /= m .or. &
                size(self%z, 3) /= (n - 1) * self%time_varying(1) + 1) then
            message = "z has incompatible shape"
            return
        end if
        if (size(self%h, 1) /= p .or. size(self%h, 2) /= p .or. &
                size(self%h, 3) /= (n - 1) * self%time_varying(2) + 1) then
            message = "h has incompatible shape"
            return
        end if
        if (size(self%tmat, 1) /= m .or. size(self%tmat, 2) /= m .or. &
                size(self%tmat, 3) /= (n - 1) * self%time_varying(3) + 1) then
            message = "tmat has incompatible shape"
            return
        end if
        if (size(self%rmat, 1) /= m .or. size(self%rmat, 2) /= r .or. &
                size(self%rmat, 3) /= (n - 1) * self%time_varying(4) + 1) then
            message = "rmat has incompatible shape"
            return
        end if
        if (size(self%q, 1) /= r .or. size(self%q, 2) /= r .or. &
                size(self%q, 3) /= (n - 1) * self%time_varying(5) + 1) then
            message = "q has incompatible shape"
            return
        end if
        if (size(self%p1, 1) /= m .or. size(self%p1, 2) /= m .or. &
                size(self%p1inf, 1) /= m .or. size(self%p1inf, 2) /= m) then
            message = "p1 and p1inf must have shape (m, m)"
            return
        end if
        if (self%diffuse_rank < 0 .or. self%diffuse_rank > m) then
            message = "diffuse_rank must be between 0 and m"
            return
        end if

        ok = .true.
    end subroutine validate_model

    subroutine kfas_approximate_nongaussian(model, u, distribution, result, theta, &
            maxiter, filter_tol, convergence_tol, expected, h_tol, info)
        type(kfas_model), intent(in) :: model
        real(dp), intent(in) :: u(:, :)
        integer, intent(in) :: distribution(:)
        type(kfas_approx_result), intent(out) :: result
        real(dp), intent(in), optional :: theta(:, :)
        integer, intent(in), optional :: maxiter
        real(dp), intent(in), optional :: filter_tol, convergence_tol, h_tol
        logical, intent(in), optional :: expected
        integer, intent(out), optional :: info
        logical :: expected_local, ok
        character(len = :), allocatable :: message
        integer :: expected_flag, local_info, max_iterations, m, n, p, r
        real(dp) :: conv_tolerance, kalman_tolerance

        call model%validate(ok, message)
        if (.not. ok) then
            result%info = 1
            if (present(info)) info = result%info
            return
        end if
        n = size(model%y, 1)
        p = size(model%y, 2)
        m = size(model%a1)
        r = size(model%q, 1)
        if (size(u, 1) /= n .or. size(u, 2) /= p .or. size(distribution) /= p) then
            result%info = 2
            if (present(info)) info = result%info
            return
        end if
        if (any(distribution < kfas_gaussian) .or. &
                any(distribution > kfas_negative_binomial)) then
            result%info = 3
            if (present(info)) info = result%info
            return
        end if

        allocate(result%theta(n, p), result%ytilde(n, p), result%htilde(p, p, n))
        if (present(theta)) then
            if (size(theta, 1) /= n .or. size(theta, 2) /= p) then
                result%info = 4
                if (present(info)) info = result%info
                return
            end if
            result%theta = theta
        else
            call kfas_init_theta(model%y, u, distribution, model%missing, result%theta)
        end if
        result%ytilde = 0.0_dp
        result%htilde = 0.0_dp
        result%gaussian_loglik = 0.0_dp
        result%difference = 0.0_dp
        result%h_tol = 1.0e15_dp
        if (present(h_tol)) result%h_tol = h_tol
        max_iterations = 50
        if (present(maxiter)) max_iterations = maxiter
        if (max_iterations < 1) then
            result%info = 5
            if (present(info)) info = result%info
            return
        end if
        kalman_tolerance = 1.0e-8_dp
        if (present(filter_tol)) kalman_tolerance = filter_tol
        conv_tolerance = 1.0e-8_dp
        if (present(convergence_tol)) conv_tolerance = convergence_tol
        expected_local = .false.
        if (present(expected)) expected_local = expected
        expected_flag = merge(1, 0, expected_local)
        local_info = 0

        call approx(model%y, model%missing, model%time_varying, model%z, model%tmat, &
            model%rmat, result%htilde, model%q, model%a1, model%p1, model%p1inf, &
            p, n, m, r, result%theta, u, result%ytilde, distribution, max_iterations, &
            kalman_tolerance, model%diffuse_rank, conv_tolerance, result%difference, &
            result%gaussian_loglik, local_info, expected_flag, result%h_tol)
        result%iterations = max_iterations
        result%info = local_info
        if (present(info)) info = local_info
    end subroutine kfas_approximate_nongaussian

    subroutine kfas_nongaussian_loglik(model, u, distribution, loglik, theta, maxiter, &
            filter_tol, convergence_tol, marginal, expected, h_tol, difference, info)
        type(kfas_model), intent(in) :: model
        real(dp), intent(in) :: u(:, :)
        integer, intent(in) :: distribution(:)
        real(dp), intent(out) :: loglik
        real(dp), intent(inout), optional :: theta(:, :)
        integer, intent(in), optional :: maxiter
        real(dp), intent(in), optional :: filter_tol, convergence_tol, h_tol
        logical, intent(in), optional :: marginal, expected
        real(dp), intent(out), optional :: difference
        integer, intent(out), optional :: info
        logical :: expected_local, marginal_local, ok
        character(len = :), allocatable :: message
        integer :: antit, expected_flag, i, local_info, marginal_flag, max_iterations
        integer :: m, n, nnd, nsim, nsim2, p, r, sim
        real(dp) :: conv_tolerance, diff, h_tolerance, kalman_tolerance
        real(dp), allocatable :: aplus1(:, :), c(:), epsplus(:, :, :), etaplus(:, :, :)
        real(dp), allocatable :: theta_work(:, :)

        loglik = -huge(1.0_dp)
        call model%validate(ok, message)
        if (.not. ok) then
            if (present(info)) info = 1
            return
        end if
        n = size(model%y, 1)
        p = size(model%y, 2)
        m = size(model%a1)
        r = size(model%q, 1)
        if (size(u, 1) /= n .or. size(u, 2) /= p .or. size(distribution) /= p) then
            if (present(info)) info = 2
            return
        end if
        if (any(distribution < kfas_gaussian) .or. &
                any(distribution > kfas_negative_binomial)) then
            if (present(info)) info = 3
            return
        end if

        allocate(theta_work(n, p))
        if (present(theta)) then
            if (size(theta, 1) /= n .or. size(theta, 2) /= p) then
                if (present(info)) info = 4
                return
            end if
            theta_work = theta
        else
            call kfas_init_theta(model%y, u, distribution, model%missing, theta_work)
        end if

        max_iterations = 50
        if (present(maxiter)) max_iterations = maxiter
        if (max_iterations < 1) then
            if (present(info)) info = 5
            return
        end if
        kalman_tolerance = 1.0e-8_dp
        if (present(filter_tol)) kalman_tolerance = filter_tol
        conv_tolerance = 1.0e-8_dp
        if (present(convergence_tol)) conv_tolerance = convergence_tol
        h_tolerance = 1.0e15_dp
        if (present(h_tol)) h_tolerance = h_tol
        marginal_local = .false.
        if (present(marginal)) marginal_local = marginal
        expected_local = .false.
        if (present(expected)) expected_local = expected
        marginal_flag = merge(1, 0, marginal_local)
        expected_flag = merge(1, 0, expected_local)

        nsim = 1
        nsim2 = 1
        sim = 0
        antit = 0
        nnd = 0
        do i = 1, m
            if (model%p1(i, i) > kalman_tolerance) nnd = nnd + 1
        end do
        allocate(epsplus(p, n, nsim), etaplus(r, n, nsim), aplus1(size(model%a1), nsim), c(nsim))
        epsplus = 0.0_dp
        etaplus = 0.0_dp
        aplus1 = 0.0_dp
        c = 0.0_dp
        diff = 0.0_dp
        local_info = 0
        loglik = 0.0_dp

        call ngloglik(model%y, model%missing, model%time_varying, model%z, model%tmat, &
            model%rmat, model%q, model%a1, model%p1, model%p1inf, p, size(model%a1), &
            r, n, loglik, theta_work, u, distribution, max_iterations, model%diffuse_rank, &
            conv_tolerance, nnd, nsim, epsplus, etaplus, aplus1, c, kalman_tolerance, &
            local_info, antit, sim, nsim2, diff, marginal_flag, expected_flag, h_tolerance)

        if (present(theta)) theta = theta_work
        if (present(difference)) difference = diff
        if (local_info /= 0 .and. local_info /= 3) loglik = -huge(1.0_dp)
        if (present(info)) info = local_info
    end subroutine kfas_nongaussian_loglik

    subroutine kfas_gaussian_filter(model, result, tol, filter_signal, info)
        type(kfas_model), intent(in) :: model
        type(kfas_filter_result), intent(out) :: result
        real(dp), intent(in), optional :: tol
        logical, intent(in), optional :: filter_signal
        integer, intent(out), optional :: info
        type(kfas_model) :: work_model
        logical :: transformed, want_signal
        integer :: local_info
        real(dp) :: tolerance

        tolerance = 1.0e-8_dp
        if (present(tol)) tolerance = tol
        want_signal = .false.
        if (present(filter_signal)) want_signal = filter_signal

        call prepare_gaussian_model(model, work_model, transformed, local_info)
        if (local_info /= 0) then
            if (present(info)) info = local_info
            return
        end if

        call gaussian_filter_core(work_model, result, tolerance)
        result%transformed_h = transformed
        if (want_signal) call compute_original_signal(model, result)
        if (present(info)) info = 0
    end subroutine kfas_gaussian_filter

    function kfas_gaussian_loglik(model, tol, info) result(loglik)
        type(kfas_model), intent(in) :: model
        real(dp), intent(in), optional :: tol
        integer, intent(out), optional :: info
        real(dp) :: loglik
        type(kfas_filter_result) :: result
        integer :: local_info

        if (present(tol)) then
            call kfas_gaussian_filter(model, result, tol = tol, info = local_info)
        else
            call kfas_gaussian_filter(model, result, info = local_info)
        end if
        if (local_info == 0) then
            loglik = result%loglik
        else
            loglik = -huge(1.0_dp)
        end if
        if (present(info)) info = local_info
    end function kfas_gaussian_loglik

    subroutine kfas_gaussian_smooth(model, result, tol, info)
        type(kfas_model), intent(in) :: model
        type(kfas_smooth_result), intent(out) :: result
        real(dp), intent(in), optional :: tol
        integer, intent(out), optional :: info
        type(kfas_model) :: work_model
        type(kfas_filter_result) :: filter_result
        logical :: transformed
        integer :: d, j, m, n, p, r, local_info
        real(dp) :: tolerance
        real(dp), allocatable :: rt0(:, :), rt1(:, :)
        real(dp), allocatable :: nt0(:, :, :), nt1(:, :, :), nt2(:, :, :)

        tolerance = 1.0e-8_dp
        if (present(tol)) tolerance = tol
        call prepare_gaussian_model(model, work_model, transformed, local_info)
        if (local_info /= 0) then
            if (present(info)) info = local_info
            return
        end if

        call gaussian_filter_core(work_model, filter_result, tolerance)
        n = size(work_model%y, 1)
        p = size(work_model%y, 2)
        m = size(work_model%a1)
        r = size(work_model%q, 1)
        d = filter_result%diffuse_end
        j = filter_result%diffuse_component

        allocate(result%state(m, n), result%state_var(m, m, n))
        allocate(result%obs_disturbance(p, n), result%obs_disturbance_var(p, n))
        allocate(result%state_disturbance(r, n), result%state_disturbance_var(r, r, n))
        allocate(result%signal(p, n), result%signal_var(p, p, n))
        allocate(result%rvec(m, n + 1), result%nmat(m, m, n + 1))
        allocate(rt0(m, d + 1), rt1(m, d + 1))
        allocate(nt0(m, m, d + 1), nt1(m, m, d + 1), nt2(m, m, d + 1))
        result%state = 0.0_dp
        result%state_var = 0.0_dp
        result%obs_disturbance = 0.0_dp
        result%obs_disturbance_var = 0.0_dp
        result%state_disturbance = 0.0_dp
        result%state_disturbance_var = 0.0_dp
        result%signal = 0.0_dp
        result%signal_var = 0.0_dp
        result%rvec = 0.0_dp
        result%nmat = 0.0_dp
        rt0 = 0.0_dp
        rt1 = 0.0_dp
        nt0 = 0.0_dp
        nt1 = 0.0_dp
        nt2 = 0.0_dp

        call gsmoothall(work_model%missing, work_model%time_varying, work_model%z, &
            work_model%h, work_model%tmat, work_model%rmat, work_model%q, p, n, m, r, &
            d, j, filter_result%a, filter_result%p, filter_result%v, filter_result%f, &
            filter_result%k, result%rvec, rt0, rt1, result%nmat, nt0, nt1, nt2, &
            filter_result%pinf(:, :, 1:d + 1), filter_result%kinf(:, :, 1:d), &
            filter_result%finf(:, 1:d), result%state, result%state_var, &
            result%obs_disturbance, result%obs_disturbance_var, result%state_disturbance, &
            result%state_disturbance_var, result%signal, result%signal_var, 1, model%z, &
            model%time_varying(1), 1, 1, 1, 1)

        result%loglik = filter_result%loglik
        result%diffuse_end = d
        result%transformed_h = transformed
        result%observation_disturbance_transformed = transformed
        if (present(info)) info = 0
    end subroutine kfas_gaussian_smooth

    subroutine kfas_ldl_transform(model, transformed, tol, info)
        type(kfas_model), intent(in) :: model
        type(kfas_model), intent(out) :: transformed
        real(dp), intent(in), optional :: tol
        integer, intent(out), optional :: info
        integer :: i, j, m, n, p, t, local_info
        integer :: timevar(5)
        integer, allocatable :: dim(:), hchol(:), hobs(:, :), ydim(:), yobs(:, :)
        real(dp) :: max_diag, tolerance
        real(dp), allocatable :: ichols(:, :, :), yt(:, :), zt(:, :, :)
        logical :: ok
        character(len = :), allocatable :: message

        call model%validate(ok, message)
        if (.not. ok) then
            if (present(info)) info = 1
            return
        end if

        n = size(model%y, 1)
        p = size(model%y, 2)
        m = size(model%a1)
        transformed = model
        if (p == 1) then
            if (present(info)) info = 0
            return
        end if

        max_diag = 0.0_dp
        do t = 1, size(model%h, 3)
            do i = 1, p
                max_diag = max(max_diag, abs(model%h(i, i, t)))
            end do
        end do
        tolerance = max(100.0_dp, max_diag) * epsilon(1.0_dp)
        if (present(tol)) tolerance = tol

        timevar = model%time_varying
        timevar(1) = 1
        timevar(2) = 1
        allocate(yt(p, n), zt(p, m, n), ichols(p, p, n))
        allocate(ydim(n), yobs(p, n), hchol(n), dim(n), hobs(p, n))
        yt = transpose(model%y)

        do t = 1, n
            zt(:, :, t) = model%z(:, :, (t - 1) * model%time_varying(1) + 1)
            ichols(:, :, t) = model%h(:, :, (t - 1) * model%time_varying(2) + 1)
            ydim(t) = count(model%missing(t, :) == 0)
            j = 0
            do i = 1, p
                if (model%missing(t, i) == 0) then
                    j = j + 1
                    yobs(j, t) = i
                end if
            end do
            do i = 1, p
                if (model%missing(t, i) /= 0) then
                    j = j + 1
                    yobs(j, t) = i
                end if
            end do
            hobs(:, t) = yobs(:, t)
            dim(t) = ydim(t)
            hchol(t) = t
        end do

        local_info = 0
        call ldlssm(yt, ydim, yobs, timevar, zt, p, m, n, ichols, n, hchol, &
            dim, local_info, hobs, tolerance)
        if (local_info /= 0) then
            if (present(info)) info = 10 + local_info
            return
        end if

        transformed%y = transpose(yt)
        transformed%z = zt
        deallocate(transformed%h)
        allocate(transformed%h(p, p, n))
        transformed%h = 0.0_dp
        do t = 1, n
            do i = 1, p
                transformed%h(i, i, t) = ichols(i, i, t)
            end do
        end do
        transformed%time_varying = timevar
        if (present(info)) info = 0
    end subroutine kfas_ldl_transform

    subroutine kfas_init_theta(y, u, distribution, missing, theta)
        real(dp), intent(in) :: y(:, :), u(:, :)
        integer, intent(in) :: distribution(:)
        integer, intent(in) :: missing(:, :)
        real(dp), intent(out) :: theta(:, :)
        integer :: i, j, n, p
        real(dp) :: prob, x

        n = size(y, 1)
        p = size(y, 2)
        if (size(u, 1) /= n .or. size(u, 2) /= p) error stop "u has incompatible shape"
        if (size(missing, 1) /= n .or. size(missing, 2) /= p) then
            error stop "missing has incompatible shape"
        end if
        if (size(distribution) /= p) error stop "distribution has incompatible shape"
        if (size(theta, 1) /= n .or. size(theta, 2) /= p) then
            error stop "theta has incompatible shape"
        end if

        theta = y
        do j = 1, p
            select case (distribution(j))
            case (kfas_gaussian)
                do i = 1, n
                    if (missing(i, j) /= 0) theta(i, j) = 0.0_dp
                end do
            case (kfas_poisson)
                do i = 1, n
                    if (missing(i, j) /= 0) then
                        x = 0.1_dp
                    else
                        x = y(i, j) / u(i, j)
                        if (x < 0.1_dp) x = 0.1_dp
                    end if
                    theta(i, j) = log(x)
                end do
            case (kfas_binomial)
                do i = 1, n
                    if (missing(i, j) /= 0) then
                        prob = 1.0_dp / (u(i, j) + 1.0_dp)
                    else
                        prob = (y(i, j) + 0.5_dp) / (u(i, j) + 1.0_dp)
                    end if
                    prob = min(max(prob, tiny(1.0_dp)), 1.0_dp - epsilon(1.0_dp))
                    theta(i, j) = log(prob / (1.0_dp - prob))
                end do
            case (kfas_gamma)
                do i = 1, n
                    if (missing(i, j) /= 0) then
                        x = 1.0_dp
                    else
                        x = max(y(i, j), 1.0_dp)
                    end if
                    theta(i, j) = log(x)
                end do
            case (kfas_negative_binomial)
                do i = 1, n
                    if (missing(i, j) /= 0) then
                        x = 1.0_dp / 6.0_dp
                    else
                        x = max(y(i, j), 1.0_dp / 6.0_dp)
                    end if
                    theta(i, j) = log(x)
                end do
            case default
                error stop "unknown KFAS distribution code"
            end select
        end do
    end subroutine kfas_init_theta

    subroutine kfas_ar_transform(phi)
        real(dp), intent(inout) :: phi(:)

        call artransform(phi, size(phi))
    end subroutine kfas_ar_transform

    subroutine kfas_ldl_factor(a, tol, info)
        real(dp), intent(inout) :: a(:, :)
        real(dp), intent(in), optional :: tol
        integer, intent(out), optional :: info
        integer :: local_info
        real(dp) :: tolerance

        if (size(a, 1) /= size(a, 2)) then
            if (present(info)) info = -2
            return
        end if
        tolerance = 1.0e-12_dp
        if (present(tol)) tolerance = tol
        local_info = 0
        call ldl(a, size(a, 1), tolerance, local_info)
        if (present(info)) info = local_info
    end subroutine kfas_ldl_factor

    subroutine kfas_weighted_mean_cov(x, w, meanx, covx)
        real(dp), intent(in) :: x(:, :, :)
        real(dp), intent(in) :: w(:)
        real(dp), intent(out) :: meanx(:, :)
        real(dp), intent(out) :: covx(:, :, :)
        integer :: k, m, n

        m = size(x, 1)
        n = size(x, 2)
        k = size(x, 3)
        meanx = 0.0_dp
        covx = 0.0_dp
        call covmeanwprotect(x, w, m, n, k, meanx, covx)
    end subroutine kfas_weighted_mean_cov

    subroutine prepare_gaussian_model(model, work_model, transformed, info)
        type(kfas_model), intent(in) :: model
        type(kfas_model), intent(out) :: work_model
        logical, intent(out) :: transformed
        integer, intent(out) :: info
        logical :: ok
        character(len = :), allocatable :: message
        real(dp) :: transform_tolerance

        call model%validate(ok, message)
        if (.not. ok) then
            info = 1
            transformed = .false.
            return
        end if

        transform_tolerance = default_transform_tolerance(model)
        transformed = has_nondiagonal_h(model, transform_tolerance)
        if (transformed) then
            call kfas_ldl_transform(model, work_model, tol = transform_tolerance, info = info)
        else
            work_model = model
            info = 0
        end if
    end subroutine prepare_gaussian_model

    real(dp) function default_transform_tolerance(model) result(tol)
        type(kfas_model), intent(in) :: model
        integer :: i, t
        real(dp) :: max_diag

        max_diag = 0.0_dp
        do t = 1, size(model%h, 3)
            do i = 1, size(model%h, 1)
                max_diag = max(max_diag, abs(model%h(i, i, t)))
            end do
        end do
        tol = max(100.0_dp, max_diag) * epsilon(1.0_dp)
    end function default_transform_tolerance

    logical function has_nondiagonal_h(model, tol) result(found)
        type(kfas_model), intent(in) :: model
        real(dp), intent(in) :: tol
        integer :: i, j, t

        found = .false.
        do t = 1, size(model%h, 3)
            do j = 1, size(model%h, 2)
                do i = 1, size(model%h, 1)
                    if (i /= j .and. abs(model%h(i, j, t)) > tol) then
                        found = .true.
                        return
                    end if
                end do
            end do
        end do
    end function has_nondiagonal_h

    subroutine gaussian_filter_core(model, result, tol)
        type(kfas_model), intent(in) :: model
        type(kfas_filter_result), intent(out) :: result
        real(dp), intent(in) :: tol
        integer :: d, j, m, n, p, r, rankp

        n = size(model%y, 1)
        p = size(model%y, 2)
        m = size(model%a1)
        r = size(model%q, 1)

        allocate(result%a(m, n + 1), result%p(m, m, n + 1), result%pinf(m, m, n + 1))
        allocate(result%att(m, n), result%ptt(m, m, n))
        allocate(result%v(p, n), result%f(p, n), result%finf(p, n))
        allocate(result%k(m, p, n), result%kinf(m, p, n))
        allocate(result%signal(n, p), result%signal_var(p, p, n))
        result%a = 0.0_dp
        result%p = 0.0_dp
        result%pinf = 0.0_dp
        result%att = 0.0_dp
        result%ptt = 0.0_dp
        result%v = 0.0_dp
        result%f = 0.0_dp
        result%finf = 0.0_dp
        result%k = 0.0_dp
        result%kinf = 0.0_dp
        result%signal = 0.0_dp
        result%signal_var = 0.0_dp
        result%loglik = 0.0_dp

        d = 0
        j = 0
        rankp = model%diffuse_rank
        call kfilter2(model%y, model%missing, model%time_varying, model%z, model%h, &
            model%tmat, model%rmat, model%q, model%a1, model%p1, model%p1inf, &
            p, n, m, r, d, j, result%a, result%p, result%v, result%f, result%k, &
            result%pinf, result%finf, result%kinf, result%loglik, tol, rankp, &
            result%signal, result%signal_var, 0, result%att, result%ptt)

        result%diffuse_end = d
        result%diffuse_component = j
        result%remaining_diffuse_rank = rankp
    end subroutine gaussian_filter_core

    subroutine compute_original_signal(model, result)
        type(kfas_model), intent(in) :: model
        type(kfas_filter_result), intent(inout) :: result
        integer :: iz, n, t
        real(dp), allocatable :: temp(:, :)

        n = size(model%y, 1)
        allocate(temp(size(model%z, 1), size(model%a1)))
        do t = 1, n
            iz = (t - 1) * model%time_varying(1) + 1
            result%signal(t, :) = matmul(model%z(:, :, iz), result%a(:, t))
            temp = matmul(model%z(:, :, iz), result%p(:, :, t))
            result%signal_var(:, :, t) = matmul(temp, transpose(model%z(:, :, iz)))
        end do
    end subroutine compute_original_signal

end module kfas
