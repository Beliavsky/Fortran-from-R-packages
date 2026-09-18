module tclust_diagnostics
    use tclust_kinds, only : dp
    use tclust_types, only : tclust_result, discr_fact_result, ctlcurves_result
    use tclust_linalg, only : log_mvnorm, sort_real_with_index, qchisq
    use tclust_core, only : tclust
    implicit none
    private

    public :: discr_fact
    public :: ctlcurves

contains

    subroutine discr_fact(fit, result, threshold)
        type(tclust_result), intent(in) :: fit !! Fitted tclust model supplying data, sizes, centers, and covariance matrices.
        type(discr_fact_result), intent(out) :: result !! Per-observation best/second likelihoods and discriminant factors.
        real(dp), intent(in), optional :: threshold !! Relative-likelihood cutoff in (0,1); default 0.1 and stored on log scale.
        real(dp), allocatable :: logll(:, :)
        real(dp), allocatable :: bestlog(:)
        real(dp), allocatable :: work(:)
        real(dp), allocatable :: one_x(:, :)
        real(dp) :: one_log(1)
        integer, allocatable :: order(:)
        real(dp) :: thr
        real(dp) :: logcut
        real(dp) :: minfinite
        real(dp) :: s
        integer :: n
        integer :: k
        integer :: no_trim
        integer :: trimn
        integer :: i
        integer :: j
        integer :: best
        integer :: second

        n = fit%n
        k = fit%k
        if (.not. allocated(fit%x)) error stop 'discr_fact: fit must retain its data matrix'
        if (n < 1 .or. k < 1) error stop 'discr_fact: invalid fitted object'
        thr = 0.1_dp
        if (present(threshold)) thr = threshold
        if (thr <= 0.0_dp .or. thr >= 1.0_dp) error stop 'discr_fact: threshold must be in (0,1)'
        no_trim = floor(real(n, dp) * (1.0_dp - fit%alpha))
        trimn = n - no_trim
        allocate(logll(n, k), bestlog(n), work(n), order(n), one_x(1, fit%p))
        allocate(result%assignment(n), result%second_assignment(n), result%likelihood(n), &
                 result%second_likelihood(n), result%factor(n), result%mean_factor(0:k))
        do j = 1, k
            if (fit%size(j) > 0.0_dp) then
                do i = 1, n
                    one_x(1, :) = fit%x(i, :)
                    call log_mvnorm(one_x, fit%centers(:, j), fit%cov(:, :, j), one_log)
                    logll(i, j) = log(max(fit%size(j) / real(no_trim, dp), tiny(1.0_dp))) + one_log(1)
                end do
            else
                logll(:, j) = -huge(1.0_dp)
            end if
        end do
        do i = 1, n
            best = maxloc(logll(i, :), dim=1)
            bestlog(i) = logll(i, best)
            result%assignment(i) = best
            if (k > 1) then
                work(1:k) = logll(i, :)
                work(best) = -huge(1.0_dp)
                second = maxloc(work(1:k), dim=1)
                result%second_assignment(i) = second
                result%second_likelihood(i) = exp(max(logll(i, second), log(tiny(1.0_dp))))
                result%factor(i) = logll(i, second) - bestlog(i)
            else
                result%second_assignment(i) = 0
                result%second_likelihood(i) = exp(max(bestlog(i) + log(thr), log(tiny(1.0_dp))))
                result%factor(i) = log(thr)
            end if
            result%likelihood(i) = exp(max(bestlog(i), log(tiny(1.0_dp))))
        end do
        call sort_real_with_index(bestlog, order, ascending=.true.)
        logcut = bestlog(order(min(n, trimn + 1)))
        do i = 1, n
            if (bestlog(i) < logcut) then
                result%second_assignment(i) = result%assignment(i)
                result%second_likelihood(i) = result%likelihood(i)
                result%likelihood(i) = exp(max(logcut, log(tiny(1.0_dp))))
                result%assignment(i) = 0
                result%factor(i) = bestlog(i) - logcut
            end if
        end do
        minfinite = 0.0_dp
        if (n > 0) minfinite = minval(result%factor)
        result%ylimmin = 1.5_dp * minfinite
        result%threshold = log(thr)
        do j = 0, k
            s = 0.0_dp
            second = 0
            do i = 1, n
                if (result%assignment(i) == j) then
                    s = s + result%factor(i)
                    second = second + 1
                end if
            end do
            if (second > 0) then
                result%mean_factor(j) = s / real(second, dp)
            else
                result%mean_factor(j) = 0.0_dp
            end if
        end do
    end subroutine discr_fact

    subroutine ctlcurves(x, k_values, alpha_values, result, restr_fact, nstart, niter1, niter2, nkeep, seed)
        real(dp), intent(in) :: x(:, :) !! Numeric data matrix evaluated over the k-by-alpha grid.
        integer, intent(in) :: k_values(:) !! Positive cluster counts to evaluate.
        real(dp), intent(in) :: alpha_values(:) !! Trimming proportions to evaluate, each in [0,1).
        type(ctlcurves_result), intent(out) :: result !! Objective, weight, restriction, doubtful, and outlier diagnostic grids.
        real(dp), intent(in), optional :: restr_fact !! Eigenvalue restriction factor passed to tclust; default 50.
        integer, intent(in), optional :: nstart !! Random starts per tclust fit; default 50 in the Fortran grid helper.
        integer, intent(in), optional :: niter1 !! Initial concentration steps per fit; default 3.
        integer, intent(in), optional :: niter2 !! Refinement steps per retained fit; default 20.
        integer, intent(in), optional :: nkeep !! Number of retained starts per fit; default 5.
        integer, intent(in), optional :: seed !! Optional base RNG seed; grid cells use deterministic offsets.
        type(tclust_result), allocatable :: fit
        type(discr_fact_result) :: disc
        real(dp) :: rf
        real(dp) :: q95
        integer :: ns
        integer :: ni1
        integer :: ni2
        integer :: nk
        integer :: base_seed
        integer :: i
        integer :: j

        if (any(k_values < 1)) error stop 'ctlcurves: k_values must be positive'
        if (any(alpha_values < 0.0_dp) .or. any(alpha_values >= 1.0_dp)) error stop 'ctlcurves: alpha_values must lie in [0,1)'
        rf = 50.0_dp
        if (present(restr_fact)) rf = restr_fact
        ns = 50
        if (present(nstart)) ns = nstart
        ni1 = 3
        if (present(niter1)) ni1 = niter1
        ni2 = 20
        if (present(niter2)) ni2 = niter2
        nk = min(5, ns)
        if (present(nkeep)) nk = nkeep
        allocate(fit)
        base_seed = 104729
        if (present(seed)) base_seed = seed
        allocate(result%k_values(size(k_values)), result%alpha_values(size(alpha_values)), &
                 result%obj(size(k_values), size(alpha_values)), result%min_weights(size(k_values), size(alpha_values)), &
                 result%unrestr_fact(size(k_values), size(alpha_values)), result%doubtful(size(k_values), size(alpha_values)), &
                 result%outlying(size(k_values), size(alpha_values)))
        result%k_values = k_values
        result%alpha_values = alpha_values
        q95 = qchisq(0.95_dp, real(size(x, 2), dp))
        do i = 1, size(k_values)
            do j = 1, size(alpha_values)
                call tclust(x, k_values(i), fit, alpha=alpha_values(j), restr_fact=rf, &
                            nstart=ns, niter1=ni1, niter2=ni2, nkeep=nk, &
                            seed=base_seed + 1009 * i + 9176 * j)
                call discr_fact(fit, disc)
                result%obj(i, j) = fit%obj
                result%min_weights(i, j) = minval(fit%weights)
                result%unrestr_fact(i, j) = fit%unrestr_fact
                result%doubtful(i, j) = count(disc%factor > disc%threshold)
                result%outlying(i, j) = count(fit%cluster > 0 .and. fit%mah > q95)
            end do
        end do
    end subroutine ctlcurves

end module tclust_diagnostics
