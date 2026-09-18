module tclust_core
    use tclust_kinds, only : dp
    use tclust_types, only : tclust_result
    use tclust_rng, only : set_tclust_seed, sample_index, rand_uniform
    use tclust_linalg, only : covariance_matrix, weighted_mean_cov, log_mvnorm, &
                              mahalanobis_sq, sort_real_with_index
    use tclust_restrictions, only : restrict_covariances, unrestricted_factor
    implicit none
    private

    type :: iteration_state
        real(dp), allocatable :: centers(:, :)
        real(dp), allocatable :: cov(:, :, :)
        integer, allocatable :: cluster(:)
        real(dp), allocatable :: size(:)
        real(dp), allocatable :: weights(:)
        real(dp), allocatable :: posterior(:, :)
        real(dp) :: obj = -huge(1.0_dp)
        real(dp) :: nlogl = huge(1.0_dp)
        integer :: code = 0
    end type iteration_state

    public :: tclust
    public :: tclust_refine
    public :: tclust_information_criteria

contains

    subroutine tclust(x, k, result, alpha, nstart, niter1, niter2, nkeep, equal_weights, &
                      restriction, restr_fact, cshape, opt, zero_tol, seed)
        real(dp), intent(in) :: x(:, :) !! Numeric data matrix with observations in rows and variables in columns.
        integer, intent(in) :: k !! Initial number of clusters; must be at least one.
        type(tclust_result), intent(out) :: result !! Fitted trimmed-clustering result, including parameters and assignments.
        real(dp), intent(in), optional :: alpha !! Trimming proportion in [0,1); default 0.05.
        integer, intent(in), optional :: nstart !! Number of random initializations; default 500.
        integer, intent(in), optional :: niter1 !! Concentration steps for each initial solution; default 3.
        integer, intent(in), optional :: niter2 !! Refinement steps for retained solutions; default 20.
        integer, intent(in), optional :: nkeep !! Number of phase-one solutions retained for refinement; default 5.
        logical, intent(in), optional :: equal_weights !! Whether cluster mixing weights are fixed equal; default false.
        character(len=*), intent(in), optional :: restriction !! Scatter restriction, either 'eigen' or 'deter'; default 'eigen'.
        real(dp), intent(in), optional :: restr_fact !! Across-cluster scatter restriction factor, at least one; default 12.
        real(dp), intent(in), optional :: cshape !! Within-cluster shape restriction for determinant mode; default 1e10.
        character(len=*), intent(in), optional :: opt !! Objective mode 'HARD' or 'MIXT'; default 'HARD'.
        real(dp), intent(in), optional :: zero_tol !! Nonnegative numerical tolerance for singularity decisions; default 1e-16.
        integer, intent(in), optional :: seed !! Optional deterministic RNG seed; exact R RNG stream matching is not attempted.
        type(iteration_state) :: state
        type(iteration_state) :: candidate
        type(iteration_state) :: best
        integer, allocatable :: starts(:, :)
        integer, allocatable :: order(:)
        real(dp), allocatable :: objectives(:)
        real(dp) :: a
        real(dp) :: cs
        real(dp) :: rf
        real(dp) :: ztol
        integer :: ns
        integer :: ni1
        integer :: ni2
        integer :: nk
        integer :: j
        integer :: n
        integer :: no_trim
        logical :: eq
        character(len=8) :: rest
        character(len=4) :: objective_mode

        n = size(x, 1)
        if (n < 1 .or. size(x, 2) < 1) error stop 'tclust: x must be a nonempty matrix'
        if (k < 1) error stop 'tclust: k must be >= 1'

        a = 0.05_dp
        if (present(alpha)) a = alpha
        if (a < 0.0_dp .or. a >= 1.0_dp) error stop 'tclust: alpha must be in [0,1)'
        no_trim = floor(real(n, dp) * (1.0_dp - a))
        if (no_trim < k) error stop 'tclust: too few untrimmed observations for k clusters'

        ns = 500
        if (present(nstart)) ns = nstart
        ni1 = 3
        if (present(niter1)) ni1 = niter1
        ni2 = 20
        if (present(niter2)) ni2 = niter2
        nk = 5
        if (present(nkeep)) nk = nkeep
        if (ns < 1 .or. ni1 < 0 .or. ni2 < 0 .or. nk < 1 .or. nk > ns) error stop 'tclust: invalid iteration counts'

        eq = .false.
        if (present(equal_weights)) eq = equal_weights
        rest = 'eigen'
        if (present(restriction)) rest = trim(adjustl(restriction))
        if (trim(rest) /= 'eigen' .and. trim(rest) /= 'deter') then
            error stop 'tclust: this Fortran translation supports restriction=eigen or deter; GPCM is not yet translated'
        end if
        rf = 12.0_dp
        if (present(restr_fact)) rf = restr_fact
        if (rf < 1.0_dp) error stop 'tclust: restr_fact must be >= 1'
        cs = 1.0e10_dp
        if (present(cshape)) cs = cshape
        if (cs < 1.0_dp) error stop 'tclust: cshape must be >= 1'
        objective_mode = 'HARD'
        if (present(opt)) objective_mode = trim(adjustl(opt))
        if (objective_mode /= 'HARD' .and. objective_mode /= 'MIXT') error stop 'tclust: opt must be HARD or MIXT'
        if (eq .and. objective_mode == 'MIXT') eq = .false.
        ztol = 1.0e-16_dp
        if (present(zero_tol)) ztol = zero_tol
        if (ztol < 0.0_dp) error stop 'tclust: zero_tol must be nonnegative'
        if (present(seed)) call set_tclust_seed(seed)

        allocate(starts(n, ns), objectives(ns), order(ns))
        do j = 1, ns
            call initialize_state(x, k, no_trim, eq, state)
            call concentration_steps(x, state, a, rest, rf, cs, objective_mode, eq, ztol, ni1)
            starts(:, j) = state%cluster
            objectives(j) = state%obj
        end do
        call sort_real_with_index(objectives, order, ascending=.false.)

        best%obj = -huge(1.0_dp)
        do j = 1, nk
            call state_from_cluster(x, k, no_trim, eq, starts(:, order(j)), candidate)
            call concentration_steps(x, candidate, a, rest, rf, cs, objective_mode, eq, ztol, ni2)
            if (candidate%obj > best%obj) best = candidate
        end do

        call finalize_result(x, best, a, rest, rf, cs, objective_mode, eq, ztol, result)
    end subroutine tclust

    subroutine tclust_refine(x, cluster, result, alpha, restriction, restr_fact, cshape, opt, &
                             equal_weights, niter, zero_tol)
        real(dp), intent(in) :: x(:, :) !! Numeric data matrix with observations in rows and variables in columns.
        integer, intent(in) :: cluster(:) !! Initial labels of length n; zero denotes trimmed observations.
        type(tclust_result), intent(out) :: result !! Refined trimmed-clustering result.
        real(dp), intent(in), optional :: alpha !! Trimming proportion; inferred from labels when absent.
        character(len=*), intent(in), optional :: restriction !! Scatter restriction 'eigen' or 'deter'.
        real(dp), intent(in), optional :: restr_fact !! Across-cluster restriction factor; default 12.
        real(dp), intent(in), optional :: cshape !! Shape restriction for determinant mode; default 1e10.
        character(len=*), intent(in), optional :: opt !! Objective mode 'HARD' or 'MIXT'; default HARD.
        logical, intent(in), optional :: equal_weights !! Fix equal cluster weights when true; default false.
        integer, intent(in), optional :: niter !! Maximum concentration steps; default 20.
        real(dp), intent(in), optional :: zero_tol !! Numerical tolerance; default 1e-16.
        type(iteration_state) :: state
        real(dp) :: a
        real(dp) :: cs
        real(dp) :: rf
        real(dp) :: ztol
        integer :: nit
        integer :: k
        integer :: no_trim
        logical :: eq
        character(len=8) :: rest
        character(len=4) :: objective_mode

        if (size(cluster) /= size(x, 1)) error stop 'tclust_refine: cluster length mismatch'
        k = maxval(cluster)
        if (k < 1) error stop 'tclust_refine: at least one positive cluster is required'
        a = real(count(cluster == 0), dp) / real(size(cluster), dp)
        if (present(alpha)) a = alpha
        no_trim = floor(real(size(cluster), dp) * (1.0_dp - a))
        eq = .false.
        if (present(equal_weights)) eq = equal_weights
        rest = 'eigen'
        if (present(restriction)) rest = trim(adjustl(restriction))
        rf = 12.0_dp
        if (present(restr_fact)) rf = restr_fact
        cs = 1.0e10_dp
        if (present(cshape)) cs = cshape
        objective_mode = 'HARD'
        if (present(opt)) objective_mode = trim(adjustl(opt))
        ztol = 1.0e-16_dp
        if (present(zero_tol)) ztol = zero_tol
        nit = 20
        if (present(niter)) nit = niter
        call state_from_cluster(x, k, no_trim, eq, cluster, state)
        call concentration_steps(x, state, a, rest, rf, cs, objective_mode, eq, ztol, nit)
        call finalize_result(x, state, a, rest, rf, cs, objective_mode, eq, ztol, result)
    end subroutine tclust_refine

    subroutine initialize_state(x, k, no_trim, equal_weights, state)
        real(dp), intent(in) :: x(:, :) !! Data matrix used for random parameter initialization.
        integer, intent(in) :: k !! Number of clusters to initialize.
        integer, intent(in) :: no_trim !! Number of observations retained after trimming.
        logical, intent(in) :: equal_weights !! Whether initial cluster weights are fixed equal.
        type(iteration_state), intent(out) :: state !! Initialized working clustering state.
        real(dp), allocatable :: sub(:, :)
        real(dp) :: sw
        integer :: i
        integer :: j
        integer :: n
        integer :: p

        n = size(x, 1)
        p = size(x, 2)
        allocate(state%centers(k, p), state%cov(p, p, k), state%cluster(n), state%size(k), &
                 state%weights(k), state%posterior(n, k), sub(p + 1, p))
        state%cluster = 0
        state%posterior = 0.0_dp
        do j = 1, k
            do i = 1, p + 1
                sub(i, :) = x(sample_index(n), :)
            end do
            state%centers(j, :) = sub(1, :)
            call covariance_matrix(sub, state%cov(:, :, j), unbiased=.false.)
        end do
        if (equal_weights) then
            state%weights = 1.0_dp / real(k, dp)
            state%size = real(no_trim / k, dp)
        else
            do j = 1, k
                state%weights(j) = rand_uniform(0.0_dp, 1.0_dp)
            end do
            sw = sum(state%weights)
            if (sw <= tiny(1.0_dp)) then
                state%weights = 1.0_dp / real(k, dp)
            else
                state%weights = state%weights / sw
            end if
            state%size = anint(real(n, dp) * state%weights)
        end if
        state%code = 0
    end subroutine initialize_state

    subroutine state_from_cluster(x, k, no_trim, equal_weights, cluster, state)
        real(dp), intent(in) :: x(:, :) !! Data matrix used to reconstruct parameters from a hard assignment.
        integer, intent(in) :: k !! Number of possible positive cluster labels.
        integer, intent(in) :: no_trim !! Number of observations expected to remain untrimmed.
        logical, intent(in) :: equal_weights !! Whether to use equal cluster weights rather than empirical proportions.
        integer, intent(in) :: cluster(:) !! Hard labels of length n; zero denotes trimmed observations.
        type(iteration_state), intent(out) :: state !! Reconstructed working state ready for concentration steps.
        integer :: i
        integer :: n
        integer :: p

        n = size(x, 1)
        p = size(x, 2)
        allocate(state%centers(k, p), state%cov(p, p, k), state%cluster(n), state%size(k), &
                 state%weights(k), state%posterior(n, k))
        state%cluster = cluster
        state%posterior = 0.0_dp
        do i = 1, n
            if (cluster(i) >= 1 .and. cluster(i) <= k) state%posterior(i, cluster(i)) = 1.0_dp
        end do
        state%size = sum(state%posterior, dim=1)
        if (equal_weights) then
            state%weights = 1.0_dp / real(k, dp)
        else
            state%weights = state%size / real(max(no_trim, 1), dp)
            where (state%weights <= 0.0_dp) state%weights = tiny(1.0_dp)
            state%weights = state%weights / sum(state%weights)
        end if
        call estimate_parameters(x, state, 0.0_dp)
        state%code = 0
    end subroutine state_from_cluster

    subroutine concentration_steps(x, state, alpha, restriction, restr_fact, cshape, opt, &
                                   equal_weights, zero_tol, niter)
        real(dp), intent(in) :: x(:, :) !! Data matrix for concentration updates.
        type(iteration_state), intent(inout) :: state !! Working clustering state modified by each concentration step.
        real(dp), intent(in) :: alpha !! Trimming proportion used to determine the retained sample size.
        character(len=*), intent(in) :: restriction !! Scatter restriction mode, 'eigen' or 'deter'.
        real(dp), intent(in) :: restr_fact !! Across-cluster restriction factor.
        real(dp), intent(in) :: cshape !! Shape restriction used for determinant mode.
        character(len=*), intent(in) :: opt !! Objective mode, HARD or MIXT.
        logical, intent(in) :: equal_weights !! Whether mixture proportions remain equal.
        real(dp), intent(in) :: zero_tol !! Numerical singularity tolerance.
        integer, intent(in) :: niter !! Maximum number of concentration steps.
        integer :: code
        integer :: i
        integer :: j
        integer :: p

        p = size(x, 2)
        do i = 1, niter
            call restrict_covariances(state%cov, state%size, restriction, restr_fact, cshape, zero_tol, code)
            if (code == 0) then
                if (i == 1) then
                    state%cov = 0.0_dp
                    do j = 1, size(state%cov, 3)
                        state%cov(:, :, j) = identity_matrix(p)
                    end do
                else
                    exit
                end if
            end if
            call find_assignment(x, state, alpha, opt, equal_weights, zero_tol)
            if (state%code == 2 .and. trim(opt) == 'HARD') exit
            call estimate_parameters(x, state, zero_tol)
        end do
        call restrict_covariances(state%cov, state%size, restriction, restr_fact, cshape, zero_tol, code)
        call calculate_objective(x, state, opt, zero_tol)
    end subroutine concentration_steps

    subroutine find_assignment(x, state, alpha, opt, equal_weights, zero_tol)
        real(dp), intent(in) :: x(:, :) !! Data matrix whose observations are assigned and trimmed.
        type(iteration_state), intent(inout) :: state !! Working state receiving labels, posteriors, sizes, and weights.
        real(dp), intent(in) :: alpha !! Trimming proportion controlling the number of excluded observations.
        character(len=*), intent(in) :: opt !! HARD uses one-hot membership; MIXT uses soft posterior membership.
        logical, intent(in) :: equal_weights !! Whether weights are fixed or updated from effective cluster sizes.
        real(dp), intent(in) :: zero_tol !! Numerical covariance floor passed to Gaussian density evaluation.
        real(dp), allocatable :: ll(:, :)
        real(dp), allocatable :: logdens(:)
        real(dp), allocatable :: density(:)
        real(dp), allocatable :: pre_z(:)
        real(dp), allocatable :: post(:, :)
        integer, allocatable :: order(:)
        integer, allocatable :: old_cluster(:)
        real(dp) :: mx
        integer :: i
        integer :: j
        integer :: k
        integer :: n
        integer :: no_trim
        integer :: bestj

        n = size(x, 1)
        k = size(state%weights)
        no_trim = floor(real(n, dp) * (1.0_dp - alpha))
        allocate(ll(n, k), logdens(n), density(n), pre_z(n), post(n, k), order(n), old_cluster(n))
        old_cluster = state%cluster
        do j = 1, k
            call log_mvnorm(x, state%centers(j, :), state%cov(:, :, j), logdens, zero_tol)
            ll(:, j) = log(max(state%weights(j), tiny(1.0_dp))) + logdens
        end do
        do i = 1, n
            mx = maxval(ll(i, :))
            density(i) = sum(exp(ll(i, :) - mx))
            post(i, :) = exp(ll(i, :) - mx) / density(i)
            if (trim(opt) == 'HARD') then
                pre_z(i) = mx
            else
                pre_z(i) = log(density(i)) + mx
            end if
            bestj = maxloc(ll(i, :), dim=1)
            state%cluster(i) = bestj
        end do
        call sort_real_with_index(pre_z, order, ascending=.false.)
        do i = no_trim + 1, n
            state%cluster(order(i)) = 0
        end do

        state%posterior = 0.0_dp
        if (trim(opt) == 'MIXT') then
            state%posterior = post
            do i = 1, n
                if (state%cluster(i) == 0) state%posterior(i, :) = 0.0_dp
            end do
        else
            do i = 1, n
                if (state%cluster(i) > 0) state%posterior(i, state%cluster(i)) = 1.0_dp
            end do
        end if
        state%size = sum(state%posterior, dim=1)
        if (.not. equal_weights) then
            state%weights = state%size / real(max(no_trim, 1), dp)
            where (state%weights <= 0.0_dp) state%weights = tiny(1.0_dp)
            state%weights = state%weights / sum(state%weights)
        end if
        state%code = 1
        if (trim(opt) == 'HARD' .and. all(old_cluster == state%cluster)) state%code = 2
    end subroutine find_assignment

    subroutine estimate_parameters(x, state, zero_tol)
        real(dp), intent(in) :: x(:, :) !! Data matrix used for weighted parameter estimation.
        type(iteration_state), intent(inout) :: state !! Working state receiving updated means and covariance matrices.
        real(dp), intent(in) :: zero_tol !! Effective-size threshold below which a cluster is treated as empty.
        integer :: j
        integer :: p

        p = size(x, 2)
        do j = 1, size(state%size)
            if (state%size(j) > zero_tol) then
                call weighted_mean_cov(x, state%posterior(:, j), state%centers(j, :), state%cov(:, :, j))
            else
                state%centers(j, :) = 0.0_dp
                state%cov(:, :, j) = identity_matrix(p)
            end if
        end do
    end subroutine estimate_parameters

    subroutine calculate_objective(x, state, opt, zero_tol)
        real(dp), intent(in) :: x(:, :) !! Data matrix at which the fitted objective is evaluated.
        type(iteration_state), intent(inout) :: state !! Working state receiving objective and classification log likelihood.
        character(len=*), intent(in) :: opt !! HARD or MIXT objective convention.
        real(dp), intent(in) :: zero_tol !! Numerical covariance floor for density evaluation.
        real(dp), allocatable :: ll(:, :)
        real(dp), allocatable :: logdens(:)
        real(dp) :: mx
        real(dp) :: mixlog
        real(dp) :: classlog
        integer :: i
        integer :: j
        integer :: k
        integer :: n

        n = size(x, 1)
        k = size(state%weights)
        allocate(ll(n, k), logdens(n))
        do j = 1, k
            call log_mvnorm(x, state%centers(j, :), state%cov(:, :, j), logdens, zero_tol)
            ll(:, j) = log(max(state%weights(j), tiny(1.0_dp))) + logdens
        end do
        classlog = 0.0_dp
        mixlog = 0.0_dp
        do i = 1, n
            if (state%cluster(i) > 0) classlog = classlog + ll(i, state%cluster(i))
            mx = maxval(ll(i, :))
            mixlog = mixlog + mx + log(sum(exp(ll(i, :) - mx)))
        end do
        state%nlogl = -2.0_dp * classlog
        if (trim(opt) == 'HARD') then
            state%obj = classlog
        else
            state%obj = mixlog
        end if
    end subroutine calculate_objective

    subroutine finalize_result(x, state, alpha, restriction, restr_fact, cshape, opt, &
                               equal_weights, zero_tol, result)
        real(dp), intent(in) :: x(:, :) !! Original unscaled data matrix copied into the fitted result.
        type(iteration_state), intent(in) :: state !! Best working state to convert into the public result type.
        real(dp), intent(in) :: alpha !! Trimming proportion used by the fit.
        character(len=*), intent(in) :: restriction !! Scatter restriction mode used by the fit.
        real(dp), intent(in) :: restr_fact !! Across-cluster restriction factor used by the fit.
        real(dp), intent(in) :: cshape !! Shape restriction used by determinant mode.
        character(len=*), intent(in) :: opt !! Objective mode used by the fit.
        logical, intent(in) :: equal_weights !! Whether equal mixing proportions were imposed.
        real(dp), intent(in) :: zero_tol !! Numerical tolerance used for Mahalanobis calculations.
        type(tclust_result), intent(out) :: result !! Public result with active clusters sorted by decreasing size.
        integer, allocatable :: active(:)
        integer, allocatable :: map(:)
        integer :: i
        integer :: j
        integer :: k0
        integer :: kreal
        integer :: m
        integer :: n
        integer :: p

        n = size(x, 1)
        p = size(x, 2)
        k0 = size(state%size)
        allocate(active(k0), map(0:k0))
        kreal = 0
        do j = 1, k0
            if (state%size(j) > zero_tol) then
                kreal = kreal + 1
                active(kreal) = j
            end if
        end do
        if (kreal < 1) error stop 'finalize_result: no nonempty cluster remained'
        do i = 1, kreal - 1
            m = i
            do j = i + 1, kreal
                if (state%size(active(j)) > state%size(active(m))) m = j
            end do
            if (m /= i) call swap_integer(active(i), active(m))
        end do
        map = 0
        do j = 1, kreal
            map(active(j)) = j
        end do

        result%n = n
        result%p = p
        result%k = kreal
        result%code = state%code
        result%alpha = alpha
        result%obj = state%obj
        result%nlogl = state%nlogl
        result%restr_fact = restr_fact
        result%cshape = cshape
        result%restriction = trim(adjustl(restriction))
        result%opt = trim(adjustl(opt))
        result%equal_weights = equal_weights
        allocate(result%cluster(n), result%size(kreal), result%weights(kreal), result%centers(p, kreal), &
                 result%cov(p, p, kreal), result%posterior(n, kreal), result%mah(n), result%x(n, p))
        result%x = x
        do i = 1, n
            if (state%cluster(i) > 0) then
                result%cluster(i) = map(state%cluster(i))
            else
                result%cluster(i) = 0
            end if
        end do
        do j = 1, kreal
            result%size(j) = state%size(active(j))
            result%weights(j) = state%weights(active(j))
            result%centers(:, j) = state%centers(active(j), :)
            result%cov(:, :, j) = state%cov(:, :, active(j))
            result%posterior(:, j) = state%posterior(:, active(j))
        end do
        result%mah = huge(1.0_dp)
        do i = 1, n
            if (result%cluster(i) == 0) then
                result%mah(i) = huge(1.0_dp)
            else
                j = result%cluster(i)
                result%mah(i) = mahalanobis_sq(x(i, :), result%centers(:, j), result%cov(:, :, j), zero_tol)
            end if
        end do
        result%unrestr_fact = unrestricted_factor(x, result%cluster, kreal, trim(restriction) == 'deter')
        call tclust_information_criteria(result)
    end subroutine finalize_result

    subroutine tclust_information_criteria(result)
        type(tclust_result), intent(inout) :: result !! Fitted result whose CLACLA/MIXMIX/MIXCLA fields are updated in place.
        real(dp) :: h
        real(dp) :: logh
        real(dp) :: npar
        integer :: k
        integer :: p

        p = result%p
        k = result%k
        h = floor((1.0_dp - result%alpha) * real(result%n, dp))
        npar = real(p * k, dp)
        if (.not. result%equal_weights) npar = npar + real(k - 1, dp)
        select case (trim(result%restriction))
        case ('eigen')
            npar = npar + 0.5_dp * real(p * (p - 1) * k, dp) + &
                   real(p * k - 1, dp) * (1.0_dp - 1.0_dp / result%restr_fact) + 1.0_dp
        case ('deter')
            npar = npar + 0.5_dp * real(p * (p - 1) * k, dp) + &
                   real(k - 1, dp) * (1.0_dp - 1.0_dp / result%restr_fact ** (1.0_dp / real(p, dp))) + &
                   1.0_dp + real(k * (p - 1), dp) * (1.0_dp - 1.0_dp / result%cshape)
        case default
            error stop 'tclust_information_criteria: unsupported restriction'
        end select
        logh = log(max(h, 1.0_dp))
        if (trim(result%opt) == 'HARD') then
            result%clacla = result%nlogl + npar * logh
            result%mixmix = huge(1.0_dp)
            result%mixcla = huge(1.0_dp)
        else
            result%mixcla = result%nlogl + npar * logh
            result%mixmix = -2.0_dp * result%obj + npar * logh
            result%clacla = huge(1.0_dp)
        end if
    end subroutine tclust_information_criteria

    pure function identity_matrix(n) result(a)
        integer, intent(in) :: n !! Order of the identity matrix to create.
        real(dp) :: a(n, n)
        integer :: i

        a = 0.0_dp
        do i = 1, n
            a(i, i) = 1.0_dp
        end do
    end function identity_matrix

    subroutine swap_integer(a, b)
        integer, intent(inout) :: a !! First integer to swap.
        integer, intent(inout) :: b !! Second integer to swap.
        integer :: t

        t = a
        a = b
        b = t
    end subroutine swap_integer

end module tclust_core
