module otrimle_simulation
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, ieee_value, ieee_quiet_nan, ieee_positive_inf
    use otrimle_kinds, only : dp
    use otrimle_types, only : otrimle_fit, otrimle_grid_point, otrimle_grid_result, &
        otrimle_sim_result, otrimle_sim_summary
    use otrimle_core, only : rimle, otrimle_fit_grid
    use otrimle_density, only : kerndenscluster
    use otrimle_linalg, only : mvn_sample, scale_tau2_location_scale
    use otrimle_rng, only : sample_weighted_index
    implicit none
    private

    public :: generator_otrimle
    public :: otrimleg
    public :: otrimlesimg
    public :: summarize_otrimlesimgdens

contains

    subroutine generator_otrimle(data, fit, generated, clustering)
        real(dp), intent(in) :: data(:, :) !! Original observations used as the empirical source for simulated noise points.
        type(otrimle_fit), intent(in) :: fit !! Successful fitted mixture defining class probabilities and Gaussian components.
        real(dp), allocatable, intent(out) :: generated(:, :) !! Simulated data matrix with the same n by p shape as data.
        integer, allocatable, intent(out) :: clustering(:) !! Generated component labels, with zero denoting the noise component.
        integer, allocatable :: count_by_component(:)
        real(dp), allocatable :: noise_weights(:)
        real(dp), allocatable :: z(:)
        integer :: comp
        integer :: i
        integer :: j
        integer :: n
        integer :: p
        integer :: pos
        integer :: src

        if (.not. allocated(fit%exproportion) .or. .not. allocated(fit%tau)) then
            error stop 'generator_otrimle: fit is incomplete'
        end if
        n = size(data, 1)
        p = size(data, 2)
        allocate(generated(n, p), clustering(n), count_by_component(fit%g + 1), noise_weights(n), z(p))
        count_by_component = 0
        do i = 1, n
            comp = sample_weighted_index(fit%exproportion)
            count_by_component(comp) = count_by_component(comp) + 1
        end do
        noise_weights = fit%tau(:, 1)
        pos = 0
        do i = 1, count_by_component(1)
            pos = pos + 1
            src = sample_weighted_index(noise_weights)
            generated(pos, :) = data(src, :)
            clustering(pos) = 0
        end do
        do j = 1, fit%g
            do i = 1, count_by_component(j + 1)
                pos = pos + 1
                call mvn_sample(fit%mean(:, j), fit%cov(:, :, j), z)
                generated(pos, :) = z
                clustering(pos) = j
            end do
        end do
    end subroutine generator_otrimle

    subroutine otrimleg(dataset, g_values, result, erc, beta0, fixlogicd, dmaxq, kernn)
        real(dp), intent(in) :: dataset(:, :) !! Observations in rows used for every candidate cluster count.
        integer, intent(in) :: g_values(:) !! Positive candidate numbers of Gaussian clusters.
        type(otrimle_grid_result), intent(out) :: result !! Fits and OTRIMLE/density criteria for all requested cluster counts.
        real(dp), intent(in), optional :: erc !! Eigenratio constraint forwarded to fitting; default 20.
        real(dp), intent(in), optional :: beta0 !! OTRIMLE criterion noise penalty used when logicd is estimated; default zero.
        real(dp), intent(in), optional :: fixlogicd(:) !! Optional fixed log improper densities, one per candidate or one shared
                                                    !! value.
        real(dp), intent(in), optional :: dmaxq !! Kernel-density grid endpoint; default qnorm(0.9995).
        integer, intent(in), optional :: kernn !! Even kernel-density grid size; default 100.
        type(otrimle_grid_point), allocatable :: opt(:)
        real(dp), allocatable :: dd(:)
        real(dp) :: ratio
        real(dp) :: bpen
        real(dp) :: dens
        real(dp) :: ld
        integer :: i
        integer :: k
        integer :: maxg
        integer :: n
        integer :: p

        ratio = 20.0_dp
        if (present(erc)) ratio = erc
        bpen = 0.0_dp
        if (present(beta0)) bpen = beta0
        n = size(dataset, 1)
        p = size(dataset, 2)
        maxg = maxval(g_values)
        allocate(result%g_values(size(g_values)), result%solution(size(g_values)))
        allocate(result%npar(size(g_values)), result%ibic(size(g_values)), result%criterion(size(g_values)))
        allocate(result%iloglik(size(g_values)), result%noiseprob(size(g_values)), result%logicd(size(g_values)))
        allocate(result%denscrit(size(g_values)), result%ddpm(size(g_values), maxg))
        result%g_values = g_values
        result%npar = 0.0_dp
        result%ibic = huge(1.0_dp)
        result%criterion = huge(1.0_dp)
        result%iloglik = -huge(1.0_dp)
        result%noiseprob = 0.0_dp
        result%logicd = 0.0_dp
        result%denscrit = huge(1.0_dp)
        result%ddpm = 0.0_dp
        do i = 1, size(g_values)
            k = g_values(i)
            if (present(fixlogicd)) then
                if (size(fixlogicd) == 1) then
                    ld = fixlogicd(1)
                else
                    if (size(fixlogicd) /= size(g_values)) error stop 'otrimleg: fixlogicd size mismatch'
                    ld = fixlogicd(i)
                end if
                call rimle(dataset, k, result%solution(i), logicd=ld, erc=ratio)
                if (result%solution(i)%code == 0) then
                    call otrimle_fit_grid(dataset, k, result%solution(i), opt, erc=ratio, beta=bpen)
                end if
            else
                call otrimle_fit_grid(dataset, k, result%solution(i), opt, erc=ratio, beta=bpen)
            end if
            result%npar(i) = real(k + k * p + k * (p + 1) * p / 2, dp)
            if (result%solution(i)%code > 0) then
                result%iloglik(i) = result%solution(i)%iloglik
                result%ibic(i) = -2.0_dp * result%iloglik(i) + result%npar(i) * log(real(n, dp))
                result%criterion(i) = result%solution(i)%criterion
                result%logicd(i) = result%solution(i)%logicd
                result%noiseprob(i) = result%solution(i)%exproportion(1)
                if (present(dmaxq) .and. present(kernn)) then
                    call kerndenscluster(dataset, result%solution(i), dens, dd, maxq=dmaxq, kernn=kernn)
                else if (present(dmaxq)) then
                    call kerndenscluster(dataset, result%solution(i), dens, dd, maxq=dmaxq)
                else if (present(kernn)) then
                    call kerndenscluster(dataset, result%solution(i), dens, dd, kernn=kernn)
                else
                    call kerndenscluster(dataset, result%solution(i), dens, dd)
                end if
                result%denscrit(i) = dens
                result%ddpm(i, 1:k) = dd
            end if
        end do
    end subroutine otrimleg

    subroutine otrimlesimg(dataset, g_values, simruns, result, erc, beta0, sim_est_logicd)
        real(dp), intent(in) :: dataset(:, :) !! Original observations whose candidate fits seed the simulation study.
        integer, intent(in) :: g_values(:) !! Candidate cluster counts fitted to the original and simulated data.
        integer, intent(in) :: simruns !! Positive number of simulated data sets generated per candidate model.
        type(otrimle_sim_result), intent(out) :: result !! Original fits and simulated density criteria.
        real(dp), intent(in), optional :: erc !! Eigenratio constraint used in all fits; default 20.
        real(dp), intent(in), optional :: beta0 !! OTRIMLE noise penalty used in all fits; default zero.
        logical, intent(in), optional :: sim_est_logicd !! If true re-estimate logicd in simulations; otherwise reuse each
                                                        !! original fit.
        type(otrimle_grid_result) :: one
        real(dp), allocatable :: generated(:, :)
        real(dp), allocatable :: fixed(:)
        integer, allocatable :: labels(:)
        integer :: i
        integer :: r
        logical :: estimate_ld
        real(dp) :: ratio
        real(dp) :: bpen

        if (simruns < 1) error stop 'otrimlesimg: simruns must be positive'
        ratio = 20.0_dp
        if (present(erc)) ratio = erc
        bpen = 0.0_dp
        if (present(beta0)) bpen = beta0
        estimate_ld = .false.
        if (present(sim_est_logicd)) estimate_ld = sim_est_logicd
        call otrimleg(dataset, g_values, result%result, erc=ratio, beta0=bpen)
        result%simruns = simruns
        allocate(result%sim_denscrit(size(g_values), simruns), result%sim_ok(size(g_values), simruns))
        result%sim_denscrit = huge(1.0_dp)
        result%sim_ok = .false.
        do i = 1, size(g_values)
            if (result%result%solution(i)%code == 0) cycle
            do r = 1, simruns
                call generator_otrimle(dataset, result%result%solution(i), generated, labels)
                if (estimate_ld) then
                    call otrimleg(generated, [g_values(i)], one, erc=ratio, beta0=bpen)
                else
                    allocate(fixed(1))
                    fixed(1) = result%result%solution(i)%logicd
                    call otrimleg(generated, [g_values(i)], one, erc=ratio, beta0=bpen, fixlogicd=fixed)
                    deallocate(fixed)
                end if
                if (one%solution(1)%code > 0 .and. ieee_is_finite(one%denscrit(1))) then
                    result%sim_denscrit(i, r) = one%denscrit(1)
                    result%sim_ok(i, r) = .true.
                end if
            end do
        end do
    end subroutine otrimlesimg

    subroutine summarize_otrimlesimgdens(object, summary, noisepenalty, sdcutoff)
        type(otrimle_sim_result), intent(in) :: object !! Simulation-study result returned by otrimlesimg().
        type(otrimle_sim_summary), intent(out) :: summary !! Density-standardization and penalized model-selection summary.
        real(dp), intent(in), optional :: noisepenalty !! Noise-proportion penalty denominator; default 0.05.
        real(dp), intent(in), optional :: sdcutoff !! Maximum acceptable standardized density discrepancy; default 2.
        real(dp), allocatable :: vals(:)
        real(dp) :: penalty
        real(dp) :: cutoff
        real(dp) :: location
        real(dp) :: scale
        real(dp), allocatable :: fcriterion(:)
        real(dp) :: diff
        integer :: i
        integer :: j
        integer :: m
        integer :: bestpos
        integer :: nvalid

        penalty = 0.05_dp
        if (present(noisepenalty)) penalty = noisepenalty
        cutoff = 2.0_dp
        if (present(sdcutoff)) cutoff = sdcutoff
        m = size(object%result%g_values)
        allocate(summary%g_values(m), summary%npr(m), summary%nprdiff(m), summary%logicd(m), summary%denscrit(m))
        allocate(summary%mean_dens(m), summary%sd_dens(m), summary%standardized_dens(m), summary%penalized_g(m))
        allocate(summary%penorder(m))
        summary%g_values = object%result%g_values
        summary%sdcutoff = cutoff
        summary%npr = 0.0_dp
        summary%nprdiff = 0.0_dp
        summary%logicd = 0.0_dp
        summary%denscrit = object%result%denscrit
        summary%mean_dens = 0.0_dp
        summary%sd_dens = 0.0_dp
        summary%standardized_dens = 0.0_dp
        do i = 1, m
            if (object%result%solution(i)%code == 0) cycle
            summary%npr(i) = object%result%solution(i)%exproportion(1)
            summary%nprdiff(i) = minval(object%result%solution(i)%exproportion(2:)) - summary%npr(i)
            summary%logicd(i) = object%result%solution(i)%logicd
            nvalid = count(object%sim_ok(i, :))
            if (nvalid > 0) then
                allocate(vals(nvalid))
                nvalid = 0
                do j = 1, object%simruns
                    if (object%sim_ok(i, j)) then
                        nvalid = nvalid + 1
                        vals(nvalid) = object%sim_denscrit(i, j)
                    end if
                end do
                call scale_tau2_location_scale(vals, location, scale)
                summary%mean_dens(i) = location
                summary%sd_dens(i) = scale
                diff = summary%denscrit(i) - summary%mean_dens(i)
                if (scale > 0.0_dp) then
                    summary%standardized_dens(i) = diff / scale
                else if (diff > 0.0_dp .or. diff < 0.0_dp) then
                    summary%standardized_dens(i) = sign(ieee_value(0.0_dp, ieee_positive_inf), diff)
                else
                    summary%standardized_dens(i) = ieee_value(0.0_dp, ieee_quiet_nan)
                end if
                deallocate(vals)
            end if
            summary%penalized_g(i) = real(summary%g_values(i), dp) + summary%npr(i) / penalty
        end do
        call order_penalty(summary%penalized_g, summary%penorder)
        allocate(fcriterion(m))
        fcriterion = summary%standardized_dens
        where (ieee_is_nan(fcriterion)) fcriterion = 0.0_dp
        bestpos = minloc(fcriterion, dim=1)
        do j = 1, m
            i = summary%penorder(j)
            if (fcriterion(i) <= cutoff) then
                bestpos = i
                exit
            end if
        end do
        summary%best_g = summary%g_values(bestpos)
        if (allocated(object%result%solution(bestpos)%cluster)) then
            allocate(summary%cluster(size(object%result%solution(bestpos)%cluster)))
            summary%cluster = object%result%solution(bestpos)%cluster
        end if
    end subroutine summarize_otrimlesimgdens

    subroutine order_penalty(x, order)
        real(dp), intent(in) :: x(:) !! Penalized model-size values to rank from smallest to largest.
        integer, intent(out) :: order(:) !! Permutation of candidate positions sorted by x.
        integer :: i
        integer :: j
        integer :: k
        integer :: t

        do i = 1, size(x)
            order(i) = i
        end do
        do i = 1, size(x) - 1
            k = i
            do j = i + 1, size(x)
                if (x(order(j)) < x(order(k))) k = j
            end do
            if (k /= i) then
                t = order(i)
                order(i) = order(k)
                order(k) = t
            end if
        end do
    end subroutine order_penalty

end module otrimle_simulation
