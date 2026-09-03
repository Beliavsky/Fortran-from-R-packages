! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 computational code; see NOTICE.md and upstream/.
module mcmcglmm_distributions
    use r_kinds, only : dp
    use r_linalg, only : cholesky_factor, inverse_matrix, solve_spd
    use mcmcglmm_rng, only : rng_state, rng_chisq, rng_exponential, rng_normal, rng_uniform
    use mcmcglmm_math, only : normal_cdf
    use mcmcglmm_matrix, only : kronecker_product, mvn_log_density, sample_mvn_precision
    implicit none
    private

    public :: truncated_normal_sample
    public :: conditional_mvn_parameters
    public :: conditional_mvn_log_density
    public :: truncated_conditional_mvn_sample
    public :: truncated_conditional_mvn_log_probability
    public :: inverse_wishart_sample
    public :: riw_mcmcglmm
    public :: riw_mcmcglmm_conditioned
    public :: pkk_probability

contains

    pure subroutine truncated_normal_sample(state, mean_value, sd_value, lower, upper, sample, info)
        type(rng_state), intent(inout) :: state !! Generator state consumed by the rejection sampler.
        real(dp), intent(in) :: mean_value !! Untruncated normal mean.
        real(dp), intent(in) :: sd_value !! Positive untruncated normal standard deviation.
        real(dp), intent(in) :: lower !! Lower truncation bound; use a very negative finite value for an open lower tail.
        real(dp), intent(in) :: upper !! Upper truncation bound; use a very positive finite value for an open upper tail.
        real(dp), intent(out) :: sample !! Draw from the truncated normal, or the interval midpoint when lower is not below upper.
        integer, intent(out) :: info !! Zero on success; one for nonpositive standard deviation.
        real(dp) :: alpha
        real(dp) :: pz
        real(dp) :: slower
        real(dp) :: supper
        real(dp) :: tail_cut
        real(dp) :: tr
        real(dp) :: u
        real(dp) :: z
        logical :: accepted
        logical :: lower_open
        logical :: reflected_tail
        logical :: upper_open

        info = 0
        if (sd_value <= 0.0_dp) then
            sample = mean_value
            info = 1
            return
        end if
        if (lower >= upper) then
            sample = 0.5_dp * (lower + upper)
            return
        end if

        tail_cut = 1.0e32_dp
        lower_open = lower < -tail_cut
        upper_open = upper > tail_cut
        accepted = .false.

        if (lower_open .or. upper_open) then
            if (lower_open .and. upper_open) then
                call rng_normal(state, z)
                sample = mean_value + sd_value * z
                return
            end if
            if (upper_open) then
                tr = (lower - mean_value) / sd_value
            else
                tr = (mean_value - upper) / sd_value
            end if
            if (tr < 0.0_dp) then
                do while (.not. accepted)
                    call rng_normal(state, z)
                    accepted = z > tr
                end do
            else
                alpha = 0.5_dp * (tr + sqrt(tr * tr + 4.0_dp))
                do while (.not. accepted)
                    call rng_exponential(state, 1.0_dp / alpha, z)
                    z = z + tr
                    pz = -0.5_dp * (alpha - z) ** 2
                    call rng_exponential(state, 1.0_dp, u)
                    accepted = -u <= pz
                end do
            end if
            if (lower_open) then
                sample = mean_value - z * sd_value
            else
                sample = mean_value + z * sd_value
            end if
            return
        end if

        slower = (lower - mean_value) / sd_value
        supper = (upper - mean_value) / sd_value
        if (slower > 0.0_dp .or. supper < 0.0_dp) then
            reflected_tail = supper < 0.0_dp
            if (reflected_tail) then
                tr = -supper
                tail_cut = -slower
            else
                tr = slower
                tail_cut = supper
            end if
            alpha = 0.5_dp * (tr + sqrt(tr * tr + 4.0_dp))
            do while (.not. accepted)
                call rng_exponential(state, 1.0_dp / alpha, z)
                z = z + tr
                if (z >= tail_cut) cycle
                pz = -0.5_dp * (alpha - z) ** 2
                call rng_exponential(state, 1.0_dp, u)
                accepted = -u <= pz
            end do
            if (reflected_tail) z = -z
        else
            tr = normal_cdf(supper) - normal_cdf(slower)
            if (tr > 0.5_dp) then
                do while (.not. accepted)
                    call rng_normal(state, z)
                    accepted = z > slower .and. z < supper
                end do
            else
                do while (.not. accepted)
                    call rng_uniform(state, u)
                    z = slower + (supper - slower) * u
                    pz = -0.5_dp * z * z
                    call rng_exponential(state, 1.0_dp, u)
                    accepted = -u < pz
                end do
            end if
        end if
        sample = mean_value + sd_value * z
    end subroutine truncated_normal_sample

    pure subroutine conditional_mvn_parameters(mean_value, covariance, x, keep, cond, conditional_mean, &
                                               conditional_covariance, info)
        real(dp), intent(in) :: mean_value(:) !! Full multivariate-normal mean vector.
        real(dp), intent(in) :: covariance(:, :) !! Full symmetric positive-definite covariance matrix.
        real(dp), intent(in) :: x(:) !! Full conditioning/evaluation vector; only indices in cond are used to condition.
        integer, intent(in) :: keep(:) !! One-based indices of coordinates whose conditional distribution is requested.
        integer, intent(in) :: cond(:) !! One-based indices of observed conditioning coordinates; may be empty.
        real(dp), allocatable, intent(out) :: conditional_mean(:) !! Conditional mean for coordinates in keep.
        real(dp), allocatable, intent(out) :: conditional_covariance(:, :) !! Conditional covariance for coordinates in keep.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid indices, dimensions, or SPD solve failure.
        real(dp), allocatable :: cdev(:)
        real(dp), allocatable :: s11(:, :)
        real(dp), allocatable :: s12(:, :)
        real(dp), allocatable :: s22(:, :)
        real(dp), allocatable :: solved(:, :)
        integer :: i
        integer :: j
        integer :: n
        integer :: nk
        integer :: nc

        info = 0
        n = size(mean_value)
        nk = size(keep)
        nc = size(cond)
        if (size(x) /= n .or. size(covariance, 1) /= n .or. size(covariance, 2) /= n) then
            allocate(conditional_mean(0), conditional_covariance(0, 0))
            info = 1
            return
        end if
        if (any(keep < 1) .or. any(keep > n) .or. any(cond < 1) .or. any(cond > n)) then
            allocate(conditional_mean(0), conditional_covariance(0, 0))
            info = 2
            return
        end if
        do i = 1, nk
            if (any(cond == keep(i))) then
                allocate(conditional_mean(0), conditional_covariance(0, 0))
                info = 3
                return
            end if
        end do

        allocate(s11(nk, nk), conditional_mean(nk), conditional_covariance(nk, nk))
        do j = 1, nk
            do i = 1, nk
                s11(i, j) = covariance(keep(i), keep(j))
            end do
        end do
        conditional_mean = mean_value(keep)
        conditional_covariance = s11
        if (nc == 0) return

        allocate(s12(nk, nc), s22(nc, nc), cdev(nc))
        do j = 1, nc
            cdev(j) = x(cond(j)) - mean_value(cond(j))
            do i = 1, nk
                s12(i, j) = covariance(keep(i), cond(j))
            end do
        end do
        do j = 1, nc
            do i = 1, nc
                s22(i, j) = covariance(cond(i), cond(j))
            end do
        end do
        allocate(solved(nc, nk))
        call solve_spd(s22, transpose(s12), solved, info)
        if (info /= 0) then
            conditional_mean = 0.0_dp
            conditional_covariance = 0.0_dp
            return
        end if
        conditional_mean = mean_value(keep) + matmul(transpose(solved), cdev)
        conditional_covariance = s11 - matmul(s12, solved)
    end subroutine conditional_mvn_parameters

    pure subroutine conditional_mvn_log_density(x, mean_value, covariance, keep, cond, log_density, info)
        real(dp), intent(in) :: x(:) !! Full vector containing both evaluated and conditioning coordinates.
        real(dp), intent(in) :: mean_value(:) !! Full multivariate-normal mean vector.
        real(dp), intent(in) :: covariance(:, :) !! Full symmetric positive-definite covariance matrix.
        integer, intent(in) :: keep(:) !! One-based indices whose conditional density is evaluated.
        integer, intent(in) :: cond(:) !! One-based indices conditioned upon; may be empty.
        real(dp), intent(out) :: log_density !! Log conditional multivariate-normal density of x(keep) given x(cond).
        integer, intent(out) :: info !! Zero on success; nonzero for invalid inputs or linear-algebra failure.
        real(dp), allocatable :: cmean(:)
        real(dp), allocatable :: ccov(:, :)

        call conditional_mvn_parameters(mean_value, covariance, x, keep, cond, cmean, ccov, info)
        if (info /= 0) then
            log_density = -huge(1.0_dp)
            return
        end if
        call mvn_log_density(x(keep), cmean, ccov, log_density, info)
    end subroutine conditional_mvn_log_density

    pure subroutine truncated_conditional_mvn_sample(state, mean_value, covariance, x, keep, lower, upper, sample, info)
        type(rng_state), intent(inout) :: state !! Generator state consumed by the conditional truncated-normal sampler.
        real(dp), intent(in) :: mean_value(:) !! Full multivariate-normal mean vector.
        real(dp), intent(in) :: covariance(:, :) !! Full symmetric positive-definite covariance matrix.
        real(dp), intent(in) :: x(:) !! Current full vector; all coordinates except keep condition the requested coordinate.
        integer, intent(in) :: keep !! One-based coordinate to sample conditionally.
        real(dp), intent(in) :: lower !! Lower truncation bound for the selected coordinate.
        real(dp), intent(in) :: upper !! Upper truncation bound for the selected coordinate.
        real(dp), intent(out) :: sample !! Conditional truncated-normal draw.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid inputs or covariance failure.
        real(dp), allocatable :: cmean(:)
        real(dp), allocatable :: ccov(:, :)
        integer, allocatable :: cond(:)
        integer :: i
        integer :: index

        if (keep < 1 .or. keep > size(mean_value)) then
            sample = 0.0_dp
            info = 1
            return
        end if
        allocate(cond(size(mean_value) - 1))
        index = 0
        do i = 1, size(mean_value)
            if (i == keep) cycle
            index = index + 1
            cond(index) = i
        end do
        call conditional_mvn_parameters(mean_value, covariance, x, [keep], cond, cmean, ccov, info)
        if (info /= 0) then
            sample = 0.0_dp
            return
        end if
        call truncated_normal_sample(state, cmean(1), sqrt(ccov(1, 1)), lower, upper, sample, info)
    end subroutine truncated_conditional_mvn_sample

    pure subroutine truncated_conditional_mvn_log_probability(mean_value, covariance, x, keep, lower, upper, &
                                                              log_probability, info)
        real(dp), intent(in) :: mean_value(:) !! Full multivariate-normal mean vector.
        real(dp), intent(in) :: covariance(:, :) !! Full symmetric positive-definite covariance matrix.
        real(dp), intent(in) :: x(:) !! Current full vector; all coordinates except keep condition the requested coordinate.
        integer, intent(in) :: keep !! One-based coordinate whose interval probability is requested.
        real(dp), intent(in) :: lower !! Lower interval bound for the selected coordinate.
        real(dp), intent(in) :: upper !! Upper interval bound for the selected coordinate.
        real(dp), intent(out) :: log_probability !! Log probability that the conditional normal lies between lower and upper.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid inputs or covariance failure.
        real(dp), allocatable :: cmean(:)
        real(dp), allocatable :: ccov(:, :)
        real(dp) :: probability
        real(dp) :: sd_value
        integer, allocatable :: cond(:)
        integer :: i
        integer :: index

        if (keep < 1 .or. keep > size(mean_value)) then
            log_probability = -huge(1.0_dp)
            info = 1
            return
        end if
        allocate(cond(size(mean_value) - 1))
        index = 0
        do i = 1, size(mean_value)
            if (i == keep) cycle
            index = index + 1
            cond(index) = i
        end do
        call conditional_mvn_parameters(mean_value, covariance, x, [keep], cond, cmean, ccov, info)
        if (info /= 0) then
            log_probability = -huge(1.0_dp)
            return
        end if
        sd_value = sqrt(ccov(1, 1))
        probability = normal_cdf((upper - cmean(1)) / sd_value) - normal_cdf((lower - cmean(1)) / sd_value)
        log_probability = log(max(probability, tiny(1.0_dp)))
    end subroutine truncated_conditional_mvn_log_probability

    pure subroutine inverse_wishart_sample(state, scale_matrix, degrees_freedom, sample, info)
        type(rng_state), intent(inout) :: state !! Generator state consumed by the Bartlett decomposition.
        real(dp), intent(in) :: scale_matrix(:, :) !! Positive-definite inverse-Wishart scale matrix Psi.
        real(dp), intent(in) :: degrees_freedom !! Degrees of freedom, required to exceed dimension minus one.
        real(dp), allocatable, intent(out) :: sample(:, :) !! Allocated inverse-Wishart draw IW(Psi, degrees_freedom).
        integer, intent(out) :: info !! Zero on success; nonzero for invalid shape, degrees of freedom, or factorization failure.
        real(dp), allocatable :: inverse_scale(:, :)
        real(dp), allocatable :: factor(:, :)
        real(dp), allocatable :: bartlett(:, :)
        real(dp), allocatable :: wishart(:, :)
        real(dp), allocatable :: product(:, :)
        real(dp) :: chi
        real(dp) :: z
        integer :: i
        integer :: j
        integer :: n

        n = size(scale_matrix, 1)
        if (size(scale_matrix, 2) /= n .or. degrees_freedom <= real(n - 1, dp)) then
            allocate(sample(0, 0))
            info = 1
            return
        end if
        call inverse_matrix(scale_matrix, inverse_scale, info)
        if (info /= 0) then
            allocate(sample(0, 0))
            return
        end if
        call cholesky_factor(inverse_scale, factor, info)
        if (info /= 0) then
            allocate(sample(0, 0))
            return
        end if
        allocate(bartlett(n, n))
        bartlett = 0.0_dp
        do i = 1, n
            call rng_chisq(state, degrees_freedom - real(i - 1, dp), chi)
            bartlett(i, i) = sqrt(chi)
            do j = 1, i - 1
                call rng_normal(state, z)
                bartlett(i, j) = z
            end do
        end do
        product = matmul(factor, bartlett)
        wishart = matmul(product, transpose(product))
        call inverse_matrix(wishart, sample, info)
    end subroutine inverse_wishart_sample

    pure subroutine riw_mcmcglmm(state, v_matrix, nu, sample, info)
        type(rng_state), intent(inout) :: state !! Generator state consumed by the MCMCglmm-compatible inverse-Wishart sampler.
        real(dp), intent(in) :: v_matrix(:, :) !! MCMCglmm rIW V argument; the effective IW scale is nu times V.
        real(dp), intent(in) :: nu !! Positive MCMCglmm rIW degrees-of-freedom parameter.
        real(dp), allocatable, intent(out) :: sample(:, :) !! Draw matching unconditioned rIW(V, nu) semantics.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid input or factorization failure.

        call inverse_wishart_sample(state, nu * v_matrix, nu, sample, info)
    end subroutine riw_mcmcglmm


    pure subroutine riw_mcmcglmm_conditioned(state, v_matrix, nu, fix_index, conditioned_matrix, sample, info)
        type(rng_state), intent(inout) :: state !! Generator state consumed by the conditional inverse-Wishart sampler.
        real(dp), intent(in) :: v_matrix(:, :) !! MCMCglmm rIW V argument; nu times V is the corresponding unconditional IW scale.
        real(dp), intent(in) :: nu !! Positive MCMCglmm rIW degrees-of-freedom parameter.
        integer, intent(in) :: fix_index !! One-based first coordinate of the fixed lower-right rIW covariance block.
        real(dp), intent(in) :: conditioned_matrix(:, :) !! Fixed lower-right covariance block CM beginning at fix_index.
        real(dp), allocatable, intent(out) :: sample(:, :) !! Conditional inverse-Wishart draw with lower-right block fixed.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid shape, degrees of freedom, or factorization.
        real(dp), allocatable :: a_matrix(:, :)
        real(dp), allocatable :: a11(:, :)
        real(dp), allocatable :: a11_inverse(:, :)
        real(dp), allocatable :: a12(:, :)
        real(dp), allocatable :: a22(:, :)
        real(dp), allocatable :: half(:, :)
        real(dp), allocatable :: half_c(:, :)
        real(dp), allocatable :: iw11(:, :)
        real(dp), allocatable :: noise(:)
        real(dp), allocatable :: noise_mean(:)
        real(dp), allocatable :: sampled_mean(:)
        real(dp), allocatable :: precision(:, :)
        real(dp), allocatable :: schur(:, :)
        real(dp), allocatable :: schur_inverse(:, :)
        real(dp), allocatable :: t1(:, :)
        real(dp), allocatable :: t1_inverse(:, :)
        integer :: n
        integer :: nc
        integer :: split

        info = 0
        n = size(v_matrix, 1)
        if (size(v_matrix, 2) /= n .or. nu <= 0.0_dp .or. fix_index < 1 .or. fix_index > n) then
            allocate(sample(0, 0))
            info = 1
            return
        end if
        if (fix_index == 1) then
            allocate(sample(n, n))
            sample = v_matrix
            return
        end if
        split = fix_index - 1
        nc = n - split
        if (size(conditioned_matrix, 1) /= nc .or. size(conditioned_matrix, 2) /= nc) then
            allocate(sample(0, 0))
            info = 2
            return
        end if

        call inverse_matrix(nu * v_matrix, a_matrix, info)
        if (info /= 0) then
            allocate(sample(0, 0))
            return
        end if
        a11 = a_matrix(1:split, 1:split)
        a12 = a_matrix(1:split, split + 1:n)
        a22 = a_matrix(split + 1:n, split + 1:n)
        call inverse_matrix(a11, a11_inverse, info)
        if (info /= 0) then
            allocate(sample(0, 0))
            return
        end if
        schur = a22 - matmul(transpose(a12), matmul(a11_inverse, a12))
        call inverse_matrix(schur, schur_inverse, info)
        if (info /= 0) then
            allocate(sample(0, 0))
            return
        end if
        call inverse_wishart_sample(state, a11_inverse, nu, t1_inverse, info)
        if (info /= 0) then
            allocate(sample(0, 0))
            return
        end if
        call inverse_matrix(t1_inverse, t1, info)
        if (info /= 0) then
            allocate(sample(0, 0))
            return
        end if
        half = matmul(a11_inverse, a12)
        call kronecker_product(t1, schur_inverse, precision)
        allocate(noise_mean(split * nc))
        noise_mean = 0.0_dp
        call sample_mvn_precision(state, noise_mean, precision, noise, sampled_mean, info)
        if (info /= 0) then
            allocate(sample(0, 0))
            return
        end if
        half = -(half + reshape(noise, [split, nc]))
        half_c = matmul(conditioned_matrix, transpose(half))
        iw11 = t1_inverse + matmul(half, half_c)
        allocate(sample(n, n))
        sample = 0.0_dp
        sample(1:split, 1:split) = iw11
        sample(1:split, split + 1:n) = transpose(half_c)
        sample(split + 1:n, 1:split) = half_c
        sample(split + 1:n, split + 1:n) = conditioned_matrix
    end subroutine riw_mcmcglmm_conditioned

    pure recursive subroutine pkk_recurse(probability, exponent, start_index, depth, cumulative, accumulator)
        real(dp), intent(in) :: probability(:) !! Normalized positive category probabilities.
        real(dp), intent(in) :: exponent !! Positive occupancy sample size, permitted to be noninteger as in upstream pkk.
        integer, intent(in) :: start_index !! First category index still available to the inclusion-exclusion recursion.
        integer, intent(in) :: depth !! Current zero-based subset depth.
        real(dp), intent(in) :: cumulative !! Sum of probabilities selected on the current recursion path.
        real(dp), intent(inout) :: accumulator !! Inclusion-exclusion sum accumulated across all nonempty category subsets.
        real(dp) :: next_sum
        integer :: i
        integer :: k

        k = size(probability)
        do i = start_index, k
            next_sum = cumulative + probability(i)
            accumulator = accumulator + (-1.0_dp) ** (k - depth + 1) * next_sum ** exponent
            call pkk_recurse(probability, exponent, i + 1, depth + 1, next_sum, accumulator)
        end do
    end subroutine pkk_recurse

    pure real(dp) function pkk_probability(probability, size_value) result(value)
        real(dp), intent(in) :: probability(:) !! Nonnegative category weights; normalization is performed internally.
        real(dp), intent(in) :: size_value !! Occupancy sample size; pkk is zero below the category count.
        real(dp), allocatable :: normalized(:)
        real(dp) :: total
        integer :: k

        k = size(probability)
        if (k == 0 .or. any(probability < 0.0_dp)) then
            value = 0.0_dp
            return
        end if
        if (size_value < real(k, dp)) then
            value = 0.0_dp
            return
        end if
        total = sum(probability)
        if (total <= 0.0_dp) then
            value = 0.0_dp
            return
        end if
        normalized = probability / total
        value = 0.0_dp
        call pkk_recurse(normalized, size_value, 1, 0, 0.0_dp, value)
    end function pkk_probability

end module mcmcglmm_distributions
