! SPDX-License-Identifier: GPL-2.0-or-later
!
! Modern Fortran translation of the computational parts of RCEIM 0.3.
module rceim
    use rceim_kinds, only : dp
    use rceim_random, only : rceim_set_seed, fill_uniform, fill_normal
    use rceim_utils, only : enforce_domain, sort_population_by_score, sample_mean, sample_sd
    implicit none
    private

    public :: dp
    public :: rceim_options, rceim_result, rceim_objective
    public :: ceim_optimize
    public :: enforce_domain, sort_population_by_score

    abstract interface
        function rceim_objective(x) result(value)
            import :: dp
            real(dp), intent(in) :: x(:)
            real(dp) :: value
        end function rceim_objective
    end interface

    type :: rceim_options
        integer :: n_total = 1000
        integer :: n_elite = 0
        integer :: n_super = 1
        real(dp) :: alpha = 1.0_dp
        real(dp) :: epsilon = 0.1_dp
        real(dp) :: q = 2.0_dp
        integer :: max_iter = 50
        integer :: wait_gen = 0
        integer :: chaos_gen = 0
        logical :: minimize = .true.
        logical :: verbose = .false.
        integer :: seed = -1
    end type rceim_options

    type :: rceim_result
        real(dp), allocatable :: x(:)
        real(dp) :: value = huge(1.0_dp)
        real(dp) :: score = huge(1.0_dp)
        logical :: converged = .false.
        integer :: iterations = 0
        character(len=:), allocatable :: criterion
        real(dp), allocatable :: elite_x(:,:)
        real(dp), allocatable :: elite_score(:)
        real(dp), allocatable :: history_best(:)
        real(dp), allocatable :: history_mean(:)
        real(dp), allocatable :: history_sd(:)
        integer :: history_n = 0
    end type rceim_result
contains
    subroutine evaluate_population(fn, population, mfactor, score)
        procedure(rceim_objective) :: fn
        real(dp), intent(in) :: population(:,:)
        real(dp), intent(in) :: mfactor
        real(dp), intent(out) :: score(:)
        integer :: i
        if (size(score) /= size(population,1)) error stop "evaluate_population: incompatible dimensions"
        do i = 1, size(score)
            score(i) = mfactor*fn(population(i,:))
        end do
    end subroutine evaluate_population

    subroutine validate_inputs(lower, upper, opt, ne, ns, waitg, chaosg)
        real(dp), intent(in) :: lower(:), upper(:)
        type(rceim_options), intent(in) :: opt
        integer, intent(out) :: ne, ns, waitg, chaosg

        if (size(lower) == 0 .or. size(lower) /= size(upper)) then
            error stop "ceim_optimize: invalid boundary dimensions"
        end if
        if (any(lower >= upper)) error stop "ceim_optimize: every lower bound must be less than upper bound"
        if (opt%n_total < 2) error stop "ceim_optimize: n_total must be at least 2"
        ne = opt%n_elite
        if (ne <= 0) ne = opt%n_total/4
        if (ne < 2 .or. ne > opt%n_total) error stop "ceim_optimize: n_elite must be in [2,n_total]"
        ns = opt%n_super
        if (ns < 0 .or. ns > ne .or. ns >= opt%n_total) error stop "ceim_optimize: invalid n_super"
        if (opt%max_iter < 2) error stop "ceim_optimize: max_iter must be at least 2"
        if (opt%epsilon <= 0.0_dp) error stop "ceim_optimize: epsilon must be positive"
        if (opt%alpha < 0.0_dp .or. opt%alpha > 1.0_dp) error stop "ceim_optimize: alpha must lie in [0,1]"
        waitg = opt%wait_gen
        if (waitg <= 0) waitg = opt%max_iter
        chaosg = opt%chaos_gen
        if (chaosg <= 0) chaosg = opt%max_iter
    end subroutine validate_inputs

    subroutine ceim_optimize(fn, lower, upper, result, options)
        procedure(rceim_objective) :: fn
        real(dp), intent(in) :: lower(:), upper(:)
        type(rceim_result), intent(out) :: result
        type(rceim_options), intent(in), optional :: options

        type(rceim_options) :: opt
        real(dp), allocatable :: pop(:,:), score(:), mu(:), sig(:)
        real(dp), allocatable :: elite_x(:,:), elite_score(:)
        real(dp) :: mfactor, fracsig, alpha_d, old_elite_s
        real(dp) :: best_s, mean_s, sd_s, denom
        integer :: npar, ntot, nelite, nsuper, nnew
        integer :: waitg, chaosg, wait_counter, chaos_counter
        integer :: iter, i, imax, hist_n
        logical :: same_criterion

        if (present(options)) then
            opt = options
        else
            opt = rceim_options()
        end if
        call validate_inputs(lower, upper, opt, nelite, nsuper, waitg, chaosg)
        if (opt%seed >= 0) call rceim_set_seed(opt%seed)

        npar = size(lower)
        ntot = opt%n_total
        nnew = ntot - nsuper
        mfactor = merge(1.0_dp, -1.0_dp, opt%minimize)

        allocate(pop(ntot,npar), score(ntot), mu(npar), sig(npar))
        allocate(elite_x(nelite,npar), elite_score(nelite))
        allocate(result%history_best(opt%max_iter), result%history_mean(opt%max_iter), &
                 result%history_sd(opt%max_iter))
        result%history_best = 0.0_dp
        result%history_mean = 0.0_dp
        result%history_sd = 0.0_dp

        do i = 1, npar
            call fill_uniform(pop(:,i), lower(i), upper(i))
            mu(i) = sample_mean(pop(:,i))
            sig(i) = sample_sd(pop(:,i))
        end do

        fracsig = 10.0_dp*opt%epsilon
        iter = 1
        hist_n = 0
        wait_counter = 0
        chaos_counter = 0
        old_elite_s = 0.0_dp
        same_criterion = .false.

        do while (iter < opt%max_iter .and. opt%epsilon < fracsig .and. &
                  wait_counter <= waitg .and. .not. same_criterion)
            call evaluate_population(fn, pop, mfactor, score)
            call sort_population_by_score(pop, score)
            elite_x = pop(1:nelite,:)
            elite_score = score(1:nelite)

            best_s = elite_score(1)
            mean_s = sample_mean(elite_score)
            sd_s = sample_sd(elite_score)
            hist_n = hist_n + 1
            result%history_best(hist_n) = best_s
            result%history_mean(hist_n) = mean_s
            result%history_sd(hist_n) = sd_s

            if (abs(old_elite_s-best_s) > 0.0_dp) then
                wait_counter = 0
                chaos_counter = 0
            else
                wait_counter = wait_counter + 1
                chaos_counter = chaos_counter + 1
            end if
            old_elite_s = best_s
            same_criterion = abs(best_s-mean_s) <= 0.0_dp

            if (opt%verbose) then
                write(*,'(a,i0,a,es14.6,a,es14.6,a,es14.6)') &
                    'iteration ', iter, ': best=', best_s, ' mean=', mean_s, ' sd=', sd_s
            end if

            if (opt%epsilon < minval(sig)) then
                alpha_d = opt%alpha - opt%alpha*(1.0_dp - 1.0_dp/real(iter,dp))**opt%q
                do i = 1, npar
                    mu(i) = opt%alpha*sample_mean(elite_x(:,i)) + (1.0_dp-opt%alpha)*mu(i)
                    sig(i) = alpha_d*sample_sd(elite_x(:,i)) + (1.0_dp-alpha_d)*sig(i)
                end do
                pop = 0.0_dp
                do i = 1, npar
                    call fill_normal(pop(1:nnew,i), mu(i), sig(i))
                end do
                if (nsuper > 0) pop(nnew+1:ntot,:) = elite_x(1:nsuper,:)
                imax = maxloc(sig, dim=1)
                denom = abs(mu(imax))
                if (denom > tiny(1.0_dp)) then
                    fracsig = sig(imax)/denom
                else if (sig(imax) > 0.0_dp) then
                    fracsig = huge(1.0_dp)
                else
                    fracsig = 0.0_dp
                end if
            end if

            if (chaos_counter >= chaosg) then
                do i = 1, npar
                    call fill_normal(pop(1:nnew,i), mu(i), (upper(i)-lower(i))/2.0_dp)
                end do
                if (nsuper > 0) pop(nnew+1:ntot,:) = elite_x(1:nsuper,:)
                chaos_counter = 0
            end if
            call enforce_domain(pop, lower, upper)
            iter = iter + 1
        end do

        allocate(result%x(npar), result%elite_x(nelite,npar), result%elite_score(nelite))
        result%x = elite_x(1,:)
        result%score = elite_score(1)
        result%value = mfactor*result%score
        result%elite_x = elite_x
        result%elite_score = elite_score
        result%iterations = iter
        result%history_n = hist_n

        if (iter >= opt%max_iter) then
            result%converged = .false.
            result%criterion = "Maximum iterations reached."
        else if (wait_counter >= waitg) then
            result%converged = .false.
            result%criterion = "Maximum wait generations reached without changes in the best individual."
        else if (same_criterion) then
            result%converged = .true.
            result%criterion = "Convergence: elite mean and best are equal."
        else
            result%converged = .true.
            result%criterion = "Convergence: fractional sigma criterion reached."
        end if
    end subroutine ceim_optimize
end module rceim
