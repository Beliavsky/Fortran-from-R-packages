! SPDX-License-Identifier: GPL-2.0-or-later
! Translation of dirmult 0.1.3-5 by Torben Tvedebrink.
! See LICENSE and provenance/upstream/DESCRIPTION.

module dirmult_simulation
    use dirmult_types, only : dp, dirmult_fit_type, mom_result_type, sim_pop_result_type, &
        null_test_result_type
    use dirmult_rng, only : seed_rng, normal_random, gamma_random, multinomial_random
    use dirmult_core, only : fit_dirmult, weir_mom, multinomial_loglik
    implicit none
    private

    public :: random_dirichlet, sim_pop_sizes, sim_pop_equal_n, null_test

contains

    subroutine random_dirichlet(n, alpha, samples, info)
        integer, intent(in) :: n
        real(dp), intent(in) :: alpha(:)
        real(dp), allocatable, intent(out) :: samples(:,:)
        integer, intent(out), optional :: info
        real(dp) :: s
        integer :: i, j, stat

        stat = 0
        if (n < 0 .or. size(alpha) == 0 .or. any(alpha <= 0.0_dp)) then
            stat = 1
            allocate(samples(0,0))
            if (present(info)) info = stat
            return
        end if
        allocate(samples(n,size(alpha)))
        do i = 1, n
            do j = 1, size(alpha)
                samples(i,j) = gamma_random(alpha(j))
            end do
            s = sum(samples(i,:))
            if (s <= 0.0_dp) then
                stat = 2
                samples(i,:) = 1.0_dp / real(size(alpha),dp)
            else
                samples(i,:) = samples(i,:) / s
            end if
        end do
        if (present(info)) info = stat
    end subroutine random_dirichlet

    subroutine sim_pop_sizes(n, theta, result, pi, k, seed)
        integer, intent(in) :: n(:)
        real(dp), intent(in) :: theta
        type(sim_pop_result_type), intent(out) :: result
        real(dp), intent(in), optional :: pi(:)
        integer, intent(in), optional :: k, seed
        real(dp), allocatable :: probs(:), alpha(:), pdraw(:,:)
        integer :: j, kk, stat

        result = sim_pop_result_type()
        if (present(seed)) call seed_rng(seed)
        if (theta <= 0.0_dp .or. theta >= 1.0_dp .or. size(n) == 0 .or. any(n < 0)) then
            result%info = 1
            return
        end if
        if (present(pi)) then
            kk = size(pi)
        else
            if (.not. present(k)) then
                result%info = 2
                return
            end if
            kk = k
        end if
        if (kk <= 0) then
            result%info = 3
            return
        end if
        allocate(probs(kk), alpha(kk))
        if (present(pi)) then
            probs = pi
        else
            do j = 1, kk
                probs(j) = normal_random(14.0_dp, 4.0_dp)
            end do
        end if
        if (any(probs <= 0.0_dp) .or. sum(probs) <= 0.0_dp) then
            result%info = 4
            return
        end if
        probs = probs / sum(probs)
        alpha = probs * (1.0_dp-theta) / theta
        call random_dirichlet(size(n), alpha, pdraw, stat)
        if (stat /= 0) then
            result%info = 5
            return
        end if
        allocate(result%pi(kk), result%data(size(n),kk))
        result%theta = theta
        result%pi = probs
        do j = 1, size(n)
            call multinomial_random(n(j), pdraw(j,:), result%data(j,:))
        end do
    end subroutine sim_pop_sizes

    subroutine sim_pop_equal_n(j, k, n, theta, result, pi, seed)
        integer, intent(in) :: j, k, n
        real(dp), intent(in) :: theta
        type(sim_pop_result_type), intent(out) :: result
        real(dp), intent(in), optional :: pi(:)
        integer, intent(in), optional :: seed
        integer, allocatable :: sizes(:)

        if (j <= 0) then
            result = sim_pop_result_type()
            result%info = 1
            return
        end if
        allocate(sizes(j))
        sizes = n
        if (present(pi) .and. present(seed)) then
            call sim_pop_sizes(sizes, theta, result, pi=pi, seed=seed)
        else if (present(pi)) then
            call sim_pop_sizes(sizes, theta, result, pi=pi)
        else if (present(seed)) then
            call sim_pop_sizes(sizes, theta, result, k=k, seed=seed)
        else
            call sim_pop_sizes(sizes, theta, result, k=k)
        end if
    end subroutine sim_pop_equal_n

    subroutine null_test(data, result, m, prec, store_data, seed, maxit)
        integer, intent(in) :: data(:,:)
        type(null_test_result_type), intent(out) :: result
        integer, intent(in), optional :: m, prec, seed, maxit
        logical, intent(in), optional :: store_data
        real(dp), allocatable :: pi_null(:)
        integer, allocatable :: sim(:,:)
        type(dirmult_fit_type) :: fit
        type(mom_result_type) :: mr
        integer :: nsim, digits, lim, i, j, k, nrow, ncol
        real(dp) :: eps
        logical :: keep

        result = null_test_result_type()
        nsim = 1000
        if (present(m)) nsim = m
        digits = 6
        if (present(prec)) digits = prec
        lim = 1000
        if (present(maxit)) lim = maxit
        keep = .true.
        if (present(store_data)) keep = store_data
        if (present(seed)) call seed_rng(seed)
        if (nsim < 0 .or. sum(data) <= 0) then
            result%info = 1
            return
        end if
        eps = 10.0_dp ** (-digits)
        nrow = size(data,1)
        ncol = size(data,2)
        allocate(pi_null(ncol), sim(nrow,ncol))
        do k = 1, ncol
            pi_null(k) = real(sum(data(:,k)),dp) / real(sum(data),dp)
        end do
        allocate(result%mle_theta(nsim+1), result%dm_loglik(nsim+1))
        allocate(result%mom(nsim+1), result%mn_loglik(nsim+1), result%converged(nsim+1))
        if (keep) allocate(result%simulated(nrow,ncol,nsim))

        do i = 1, nsim
            do j = 1, nrow
                call multinomial_random(sum(data(j,:)), pi_null, sim(j,:))
            end do
            if (keep) result%simulated(:,:,i) = sim
            call fit_dirmult(sim, fit, epsilon=eps, trace=.false., maxit=lim)
            if (allocated(fit%gamma)) then
                result%mle_theta(i) = fit%theta
                result%dm_loglik(i) = fit%loglik
            else
                result%mle_theta(i) = 0.0_dp
                result%dm_loglik(i) = 0.0_dp
            end if
            result%converged(i) = fit%converged
            call weir_mom(sim, mr)
            result%mom(i) = mr%theta
            result%mn_loglik(i) = multinomial_loglik(sim)
        end do

        call fit_dirmult(data, fit, epsilon=eps, trace=.false., maxit=lim)
        if (allocated(fit%gamma)) then
            result%mle_theta(nsim+1) = fit%theta
            result%dm_loglik(nsim+1) = fit%loglik
        else
            result%mle_theta(nsim+1) = 0.0_dp
            result%dm_loglik(nsim+1) = 0.0_dp
        end if
        result%converged(nsim+1) = fit%converged
        call weir_mom(data, mr)
        result%mom(nsim+1) = mr%theta
        result%mn_loglik(nsim+1) = multinomial_loglik(data)
    end subroutine null_test

end module dirmult_simulation
