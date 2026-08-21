! SPDX-License-Identifier: GPL-2.0-or-later
! Translation of dirmult 0.1.3-5 by Torben Tvedebrink.
! See LICENSE and provenance/upstream/DESCRIPTION.

module dirmult_profiles
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use dirmult_types, only : dp, dirmult_fit_type, profile_fit_type, profile_grid_type, &
        count_table_type, real_vector_type, equal_theta_fit_type
    use dirmult_core, only : clean_counts, score_function, observed_fim, dirmult_loglik, fit_dirmult
    use dirmult_linalg, only : solve_linear
    implicit none
    private

    public :: estimate_profile_loglik, grid_profile, adaptive_grid_profile, fit_equal_theta

contains

    subroutine estimate_profile_loglik(data, theta, fit, epsilon, trace, init_pi, maxit)
        integer, intent(in) :: data(:,:)
        real(dp), intent(in) :: theta
        type(profile_fit_type), intent(out) :: fit
        real(dp), intent(in), optional :: epsilon
        logical, intent(in), optional :: trace
        real(dp), intent(in), optional :: init_pi(:)
        integer, intent(in), optional :: maxit
        integer, allocatable :: x(:,:)
        real(dp), allocatable :: gamma(:), state(:), jac(:,:), rhs(:), delta(:)
        real(dp) :: gamma_plus, eps, lik1, lik2, gsum
        integer :: k, i, ite, lim, stat
        logical :: tr, stop_after

        fit = profile_fit_type()
        if (theta <= 0.0_dp .or. theta >= 1.0_dp) then
            fit%info = 1
            return
        end if
        eps = 1.0e-4_dp
        if (present(epsilon)) eps = epsilon
        tr = .true.
        if (present(trace)) tr = trace
        lim = 1000
        if (present(maxit)) lim = maxit

        gamma_plus = (1.0_dp-theta) / theta
        call clean_counts(data, .false., .true., x)
        k = size(x,2)
        if (k == 0 .or. sum(x) <= 0) then
            fit%info = 2
            return
        end if
        allocate(gamma(k), state(k+1), jac(k+1,k+1), rhs(k+1), delta(k+1))
        if (present(init_pi)) then
            if (size(init_pi) /= k) then
                fit%info = 3
                return
            end if
            gamma = init_pi * gamma_plus
        else
            do i = 1, k
                gamma(i) = real(sum(x(:,i)),dp) / real(sum(x),dp) * gamma_plus
            end do
        end if
        state(1:k) = gamma
        state(k+1) = 1.0_dp

        lik1 = 0.0_dp
        lik2 = eps * 10.0_dp
        fit%converged = .false.
        do ite = 1, lim
            stop_after = abs(lik2-lik1) < eps
            jac = -1.0_dp
            jac(1:k,1:k) = observed_fim(x, gamma)
            jac(k+1,k+1) = 0.0_dp
            lik1 = dirmult_loglik(x, gamma) + state(k+1) * (gamma_plus-sum(gamma))

            ! This deliberately follows dirmult 0.1.3-5: profU subtracts
            ! gamma_plus from every gamma score. Under the sum constraint,
            ! a common score shift changes only the multiplier, not gamma.
            rhs(1:k) = score_function(x, gamma) - gamma_plus
            rhs(k+1) = gamma_plus - sum(gamma)
            call solve_linear(jac, rhs, delta, stat)
            if (stat /= 0) then
                fit%info = 10 + stat
                exit
            end if
            state = state - delta
            gamma = state(1:k)
            where (gamma < 0.0_dp) gamma = 0.001_dp
            if (tr) write(*,'(a,i0,a,es24.15)') 'Iteration ', ite, &
                ': Log-likelihood value: ', lik1
            lik2 = dirmult_loglik(x, gamma) + state(k+1) * (gamma_plus-sum(gamma))
            if (stop_after) then
                fit%converged = .true.
                exit
            end if
        end do

        fit%iterations = min(ite,lim)
        fit%loglik = lik1
        allocate(fit%gamma(k), fit%pi(k))
        ! Upstream returns gamlambda[1:K], rather than the temporary clamped
        ! gamma vector. Preserve that behavior for parity.
        fit%gamma = state(1:k)
        gsum = sum(fit%gamma)
        if (gsum > 0.0_dp .and. all(fit%gamma > 0.0_dp)) then
            fit%pi = fit%gamma / gsum
            fit%theta = 1.0_dp / (gsum+1.0_dp)
        else
            fit%pi = ieee_value(0.0_dp, ieee_quiet_nan)
            fit%theta = ieee_value(0.0_dp, ieee_quiet_nan)
        end if
        fit%lambda = state(k+1)
        if (.not. fit%converged .and. fit%info == 0) fit%info = 4
    end subroutine estimate_profile_loglik

    subroutine grid_profile(data, theta, from, to, len, grid, epsilon, maxit)
        integer, intent(in) :: data(:,:)
        real(dp), intent(in) :: theta, from, to
        integer, intent(in) :: len
        type(profile_grid_type), intent(out) :: grid
        real(dp), intent(in), optional :: epsilon
        integer, intent(in), optional :: maxit
        type(profile_fit_type) :: pf
        real(dp) :: z, eps
        integer :: i, lim

        grid = profile_grid_type()
        if (len <= 0) then
            grid%info = 1
            return
        end if
        eps = 1.0e-4_dp
        if (present(epsilon)) eps = epsilon
        lim = 1000
        if (present(maxit)) lim = maxit
        allocate(grid%theta(len), grid%loglik(len), grid%success(len))
        do i = 1, len
            if (len == 1) then
                z = theta + from
            else
                z = theta + from + real(i-1,dp) * (to-from) / real(len-1,dp)
            end if
            grid%theta(i) = z
            call estimate_profile_loglik(data, z, pf, epsilon=eps, trace=.false., maxit=lim)
            grid%success(i) = pf%converged .and. pf%info == 0
            if (allocated(pf%gamma)) then
                grid%loglik(i) = pf%loglik
            else
                grid%loglik(i) = ieee_value(0.0_dp, ieee_quiet_nan)
            end if
        end do
    end subroutine grid_profile

    subroutine adaptive_grid_profile(data, delta_loglik, grid, stepsize, epsilon, max_points)
        integer, intent(in) :: data(:,:)
        real(dp), intent(in) :: delta_loglik
        type(profile_grid_type), intent(out) :: grid
        integer, intent(in), optional :: stepsize, max_points
        real(dp), intent(in), optional :: epsilon
        type(dirmult_fit_type) :: mle
        type(profile_fit_type) :: pf
        real(dp), allocatable :: th(:), ll(:), pi_plus(:), pi_minus(:)
        logical, allocatable :: ok(:)
        real(dp) :: step, candidate, eps, tv, lv
        logical :: ov
        integer :: ss, cap, n, k, iter, i, j

        grid = profile_grid_type()
        ss = 50
        if (present(stepsize)) ss = stepsize
        cap = 10000
        if (present(max_points)) cap = max_points
        eps = 1.0e-8_dp
        if (present(epsilon)) eps = epsilon
        if (ss <= 0 .or. cap < 3 .or. delta_loglik < 0.0_dp) then
            grid%info = 1
            return
        end if
        call fit_dirmult(data, mle, epsilon=eps, trace=.false.)
        if (.not. allocated(mle%pi) .or. mle%theta <= 0.0_dp) then
            grid%info = 2
            return
        end if
        step = mle%theta / real(ss,dp)
        if (step <= 0.0_dp) then
            grid%info = 3
            return
        end if
        k = size(mle%pi)
        allocate(th(cap), ll(cap), ok(cap), pi_plus(k), pi_minus(k))
        n = 1
        th(1) = mle%theta
        ll(1) = mle%loglik
        ok(1) = .true.
        pi_plus = mle%pi
        pi_minus = mle%pi

        iter = 1
        do while (n < cap)
            candidate = mle%theta + real(iter,dp) * step
            if (candidate >= 1.0_dp) exit
            call estimate_profile_loglik(data, candidate, pf, epsilon=eps, trace=.false., &
                init_pi=pi_plus, maxit=1000)
            n = n + 1
            th(n) = candidate
            ll(n) = pf%loglik
            ok(n) = pf%converged .and. pf%info == 0
            if (.not. allocated(pf%pi) .or. pf%iterations > 300) then
                grid%info = 4
                exit
            end if
            pi_plus = pf%pi
            if (.not. (abs(pf%loglik-mle%loglik) <= delta_loglik)) exit
            iter = iter + 1
        end do

        iter = 1
        do while (n < cap)
            candidate = mle%theta - real(iter,dp) * step
            if (candidate <= 0.0_dp) exit
            call estimate_profile_loglik(data, candidate, pf, epsilon=eps, trace=.false., &
                init_pi=pi_minus, maxit=1000)
            n = n + 1
            th(n) = candidate
            ll(n) = pf%loglik
            ok(n) = pf%converged .and. pf%info == 0
            if (.not. allocated(pf%pi) .or. pf%iterations > 300) then
                grid%info = 5
                exit
            end if
            pi_minus = pf%pi
            if (.not. (abs(pf%loglik-mle%loglik) <= delta_loglik)) exit
            iter = iter + 1
        end do

        ! R returns rows sorted by theta.
        do i = 2, n
            tv = th(i)
            lv = ll(i)
            ov = ok(i)
            j = i - 1
            do while (j >= 1)
                if (th(j) <= tv) exit
                th(j+1) = th(j)
                ll(j+1) = ll(j)
                ok(j+1) = ok(j)
                j = j - 1
            end do
            th(j+1) = tv
            ll(j+1) = lv
            ok(j+1) = ov
        end do
        allocate(grid%theta(n), grid%loglik(n), grid%success(n))
        grid%theta = th(1:n)
        grid%loglik = ll(1:n)
        grid%success = ok(1:n)
        if (n == cap .and. grid%info == 0) grid%info = 6
    end subroutine adaptive_grid_profile

    subroutine fit_equal_theta(data, theta, fit, epsilon, trace, init_pi, maxit)
        type(count_table_type), intent(in) :: data(:)
        real(dp), intent(in) :: theta
        type(equal_theta_fit_type), intent(out) :: fit
        real(dp), intent(in), optional :: epsilon
        logical, intent(in), optional :: trace
        type(real_vector_type), intent(in), optional :: init_pi(:)
        integer, intent(in), optional :: maxit
        type(count_table_type), allocatable :: x(:)
        integer, allocatable :: kval(:), off(:)
        real(dp), allocatable :: state(:), gamma(:), jac(:,:), rhs(:), delta(:), h(:,:), sc(:)
        real(dp) :: gamma_plus, eps, lik1, lik2, lambda, gsum
        logical :: tr, stop_after
        integer :: l, i, k, nt, total_local, npar, ite, lim, stat, a, b, lam_idx, gp_idx

        fit = equal_theta_fit_type()
        if (theta <= 0.0_dp .or. theta >= 1.0_dp .or. size(data) == 0) then
            fit%info = 1
            return
        end if
        eps = 1.0e-4_dp
        if (present(epsilon)) eps = epsilon
        tr = .true.
        if (present(trace)) tr = trace
        lim = 1000
        if (present(maxit)) lim = maxit
        nt = size(data)
        if (present(init_pi)) then
            if (size(init_pi) /= nt) then
                fit%info = 2
                return
            end if
        end if

        gamma_plus = (1.0_dp-theta) / theta
        allocate(x(nt), kval(nt), off(nt+1))
        off(1) = 0
        total_local = 0
        do l = 1, nt
            if (.not. allocated(data(l)%x)) then
                fit%info = 3
                return
            end if
            call clean_counts(data(l)%x, .false., .true., x(l)%x)
            kval(l) = size(x(l)%x,2)
            if (kval(l) == 0 .or. sum(x(l)%x) <= 0) then
                fit%info = 4
                return
            end if
            total_local = total_local + kval(l) + 1
            off(l+1) = total_local
        end do
        npar = total_local + 1
        gp_idx = npar
        allocate(state(npar), jac(npar,npar), rhs(npar), delta(npar))
        state = 0.0_dp
        do l = 1, nt
            a = off(l) + 1
            k = kval(l)
            lam_idx = a + k
            if (present(init_pi)) then
                if (.not. allocated(init_pi(l)%value)) then
                    fit%info = 5
                    return
                end if
                if (size(init_pi(l)%value) /= k) then
                    fit%info = 6
                    return
                end if
                state(a:a+k-1) = init_pi(l)%value * gamma_plus
            else
                do i = 1, k
                    state(a+i-1) = real(sum(x(l)%x(:,i)),dp) / &
                        real(sum(x(l)%x),dp) * gamma_plus
                end do
            end if
            state(lam_idx) = 1.0_dp
        end do
        state(gp_idx) = gamma_plus

        lik1 = 0.0_dp
        lik2 = eps * 10.0_dp
        fit%converged = .false.
        do ite = 1, lim
            stop_after = abs(lik2-lik1) < eps
            jac = 0.0_dp
            rhs = 0.0_dp
            lik1 = 0.0_dp
            do l = 1, nt
                a = off(l) + 1
                k = kval(l)
                b = a + k - 1
                lam_idx = b + 1
                allocate(gamma(k), h(k,k), sc(k))
                gamma = state(a:b)
                where (gamma < 0.0_dp) gamma = 0.01_dp
                lambda = state(lam_idx)
                h = observed_fim(x(l)%x, gamma)
                sc = score_function(x(l)%x, gamma)
                jac(a:b,a:b) = h
                jac(a:b,lam_idx) = -1.0_dp
                jac(lam_idx,a:b) = -1.0_dp
                jac(lam_idx,gp_idx) = 1.0_dp
                jac(gp_idx,lam_idx) = 1.0_dp
                rhs(a:b) = sc - lambda
                rhs(lam_idx) = state(gp_idx) - sum(gamma)
                rhs(gp_idx) = rhs(gp_idx) + lambda
                lik1 = lik1 + dirmult_loglik(x(l)%x, gamma) + &
                    lambda * (state(gp_idx)-sum(gamma))
                deallocate(gamma, h, sc)
            end do
            call solve_linear(jac, rhs, delta, stat)
            if (stat /= 0) then
                fit%info = 10 + stat
                exit
            end if
            state = state - delta
            if (tr) write(*,'(a,i0,a,es24.15)') 'Iteration ', ite, &
                ': Log-likelihood value: ', lik1

            lik2 = 0.0_dp
            do l = 1, nt
                a = off(l) + 1
                k = kval(l)
                b = a + k - 1
                lam_idx = b + 1
                allocate(gamma(k))
                gamma = state(a:b)
                where (gamma < 0.0_dp) gamma = 0.01_dp
                lik2 = lik2 + dirmult_loglik(x(l)%x, gamma) + &
                    state(lam_idx) * (state(gp_idx)-sum(gamma))
                deallocate(gamma)
            end do
            if (stop_after) then
                fit%converged = .true.
                exit
            end if
        end do

        fit%iterations = min(ite,lim)
        fit%loglik = lik1
        fit%common_gamma_sum = state(gp_idx)
        allocate(fit%table(nt))
        do l = 1, nt
            a = off(l) + 1
            k = kval(l)
            b = a + k - 1
            lam_idx = b + 1
            allocate(fit%table(l)%gamma(k), fit%table(l)%pi(k))
            fit%table(l)%gamma = state(a:b)
            where (fit%table(l)%gamma < 0.0_dp) fit%table(l)%gamma = 0.01_dp
            gsum = sum(fit%table(l)%gamma)
            fit%table(l)%pi = fit%table(l)%gamma / gsum
            fit%table(l)%theta = 1.0_dp / (1.0_dp+gsum)
            fit%table(l)%lambda = state(lam_idx)
        end do
        if (.not. fit%converged .and. fit%info == 0) fit%info = 7
    end subroutine fit_equal_theta

end module dirmult_profiles
