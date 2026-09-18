program test_tclust
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use tclust_mod
    use tclust_linalg, only : symmetric_eigen, determinant_sym
    implicit none

    call test_indices()
    call test_clusterers()
    call test_diagnostics_and_grids()
    call test_simulators()
    print '(a)', 'All tclust tests passed.'

contains

    subroutine test_indices()
        integer :: tab(3, 3)
        integer :: c1(10)
        integer :: c2(10)
        type(rand_index_result) :: ri
        type(fm_index_result) :: fm

        tab = reshape([1, 1, 0, 1, 2, 0, 0, 1, 4], [3, 3])
        call rand_index_table(tab, ri)
        call assert_close(ri%adjusted_rand, 0.31257344300822565_dp, 2.0e-13_dp, 'adjusted Rand table')
        call assert_close(ri%rand, 0.7111111111111111_dp, 2.0e-13_dp, 'Rand table')
        call assert_close(ri%mirkin, 0.2888888888888889_dp, 2.0e-13_dp, 'Mirkin table')
        call assert_close(ri%hubert, 0.4222222222222222_dp, 2.0e-13_dp, 'Hubert table')
        call fowlkes_mallows_table(tab, fm)
        call assert_close(fm%adjusted, 0.3128799327301514_dp, 2.0e-13_dp, 'adjusted FM table')
        call assert_close(fm%raw, 0.5188745216627708_dp, 2.0e-13_dp, 'FM table')
        call assert_close(fm%expected, 0.29979416807182313_dp, 2.0e-13_dp, 'FM expectation')
        call assert_close(fm%variance, 0.011152586390681604_dp, 2.0e-13_dp, 'FM variance')

        c1 = [1, 1, 2, 2, 2, 3, 3, 3, 3, 0]
        c2 = [1, 2, 1, 2, 2, 3, 3, 3, 3, 0]
        call rand_index_labels(c1, c2, ri, noisecluster=0)
        call assert_true(ri%rand > 0.6_dp .and. ri%rand <= 1.0_dp, 'label Rand index range')
        call fowlkes_mallows_labels(c1, c2, fm, noisecluster=0)
        call assert_true(fm%raw > 0.0_dp .and. fm%raw <= 1.0_dp, 'label FM index range')
    end subroutine test_indices

    subroutine test_clusterers()
        real(dp) :: x(30, 2)
        real(dp) :: y(36, 2)
        real(dp) :: eig(2)
        real(dp) :: vec(2, 2)
        real(dp) :: evmin
        real(dp) :: evmax
        real(dp) :: detmin
        real(dp) :: detmax
        integer :: truec(30)
        integer :: truey(36)
        integer :: d(2)
        integer :: i
        integer :: j
        type(tclust_result) :: fit
        type(tclust_result) :: detfit
        type(tclust_result) :: mixfit
        type(tkmeans_result) :: km
        type(rlg_result) :: rg
        type(rand_index_result) :: ri

        call make_cluster_data(x, truec)
        call tclust(x, 2, fit, alpha=0.2_dp, nstart=30, niter1=3, niter2=20, nkeep=5, &
                    restr_fact=20.0_dp, seed=1234)
        call rand_index_labels(truec, fit%cluster, ri, noisecluster=0)
        call assert_close(ri%adjusted_rand, 1.0_dp, 1.0e-12_dp, 'tclust separated-data partition')
        call assert_true(count(fit%cluster == 0) == 6, 'tclust trims requested count')
        call assert_close(sum(fit%weights), 1.0_dp, 1.0e-12_dp, 'tclust weights sum')
        call assert_true(ieee_is_finite(fit%obj), 'tclust finite objective')
        evmin = huge(1.0_dp)
        evmax = 0.0_dp
        do j = 1, fit%k
            call symmetric_eigen(fit%cov(:, :, j), eig, vec)
            evmin = min(evmin, minval(eig))
            evmax = max(evmax, maxval(eig))
        end do
        call assert_true(evmax / max(evmin, tiny(1.0_dp)) <= 20.0_dp * (1.0_dp + 1.0e-10_dp), &
                         'eigenvalue restriction is enforced')

        call tclust_refine(x, fit%cluster, detfit, alpha=0.2_dp, restriction='deter', restr_fact=10.0_dp, &
                           cshape=5.0_dp, niter=8)
        detmin = huge(1.0_dp)
        detmax = 0.0_dp
        do j = 1, detfit%k
            detmin = min(detmin, determinant_sym(detfit%cov(:, :, j)))
            detmax = max(detmax, determinant_sym(detfit%cov(:, :, j)))
        end do
        call assert_true(detmax / max(detmin, tiny(1.0_dp)) <= 10.0_dp * (1.0_dp + 1.0e-9_dp), &
                         'determinant restriction is enforced')

        call tclust_refine(x, fit%cluster, mixfit, alpha=0.2_dp, restriction='eigen', restr_fact=20.0_dp, &
                           opt='MIXT', niter=5)
        call assert_true(ieee_is_finite(mixfit%obj), 'MIXT objective is finite')
        do i = 1, mixfit%n
            if (mixfit%cluster(i) > 0) then
                call assert_close(sum(mixfit%posterior(i, :)), 1.0_dp, 2.0e-12_dp, 'MIXT posterior row sum')
            else
                call assert_close(sum(mixfit%posterior(i, :)), 0.0_dp, 2.0e-12_dp, 'trimmed posterior row sum')
            end if
        end do

        call tkmeans(x, 2, km, alpha=0.2_dp, nstart=30, niter1=3, niter2=20, nkeep=5, seed=4321)
        call rand_index_labels(truec, km%cluster, ri, noisecluster=0)
        call assert_close(ri%adjusted_rand, 1.0_dp, 1.0e-12_dp, 'tkmeans separated-data partition')
        call assert_true(count(km%cluster == 0) == 6, 'tkmeans trims requested count')

        call make_line_data(y, truey)
        d = [1, 1]
        call rlg(y, d, rg, alpha=real(6, dp) / 36.0_dp, nstart=30, niter1=3, niter2=20, nkeep=5, seed=987)
        call rand_index_labels(truey, rg%cluster, ri, noisecluster=0)
        call assert_close(ri%adjusted_rand, 1.0_dp, 1.0e-12_dp, 'rlg separated-lines partition')
        call assert_true(rg%obj < 1.0e-9_dp, 'rlg exact line residual objective')
    end subroutine test_clusterers

    subroutine test_diagnostics_and_grids()
        real(dp) :: x(30, 2)
        integer :: truec(30)
        integer :: kvals(2)
        real(dp) :: avals(2)
        real(dp) :: cvals(2)
        integer :: bestloc2(2)
        type(tclust_result) :: fit
        type(discr_fact_result) :: disc
        type(ctlcurves_result) :: ctl
        type(tclust_ic_result) :: ic
        type(tclust_ic_solution_result) :: sol

        call make_cluster_data(x, truec)
        call tclust(x, 2, fit, alpha=0.2_dp, nstart=20, nkeep=4, restr_fact=20.0_dp, seed=765)
        call discr_fact(fit, disc)
        call assert_close(disc%threshold, log(0.1_dp), 1.0e-14_dp, 'DiscrFact threshold')
        call assert_true(all(disc%factor <= 1.0e-12_dp), 'DiscrFact second likelihood is not better than first')
        call assert_true(size(disc%mean_factor) == 3, 'DiscrFact includes outlier plus cluster means')

        kvals = [1, 2]
        avals = [0.0_dp, 0.2_dp]
        call ctlcurves(x, kvals, avals, ctl, restr_fact=20.0_dp, nstart=8, nkeep=3, seed=100)
        call assert_true(all(ieee_is_finite(ctl%obj)), 'ctlcurves objectives finite')
        call assert_true(all(ctl%min_weights >= 0.0_dp), 'ctlcurves minimum weights nonnegative')

        cvals = [1.0_dp, 20.0_dp]
        call tclust_ic(x, kvals, cvals, ic, alpha=0.2_dp, nstart=8, nkeep=3, seed=200)
        call assert_true(all(ieee_is_finite(ic%clacla)), 'tclustIC CLA finite')
        call assert_true(all(ieee_is_finite(ic%mixmix)), 'tclustIC BIC finite')
        call assert_true(all(ieee_is_finite(ic%mixcla)), 'tclustIC ICL finite')
        call tclust_ic_solutions(ic, sol, criterion='MIXMIX', nsol=2)
        bestloc2 = minloc(ic%mixmix)
        call assert_true(sol%best_k(1) == kvals(bestloc2(1)), 'tclustICsol best k matches global criterion minimum')
        call assert_close(sol%best_c(1), cvals(bestloc2(2)), 0.0_dp, 'tclustICsol best c matches global criterion minimum')
    end subroutine test_diagnostics_and_grids

    subroutine test_simulators()
        type(simulation_result) :: st
        type(simulation_result) :: sr

        call simulate_tclust(100, st, p=4, k=3, type_id=2, balanced=1, seed=11)
        call assert_true(size(st%x, 1) == 100 .and. size(st%x, 2) == 4, 'simula.tclust dimensions')
        call assert_true(count(st%true_cluster == 0) == 10, 'simula.tclust ten-percent contamination')
        call assert_true(all(st%true_cluster >= 0 .and. st%true_cluster <= 3), 'simula.tclust labels')

        call simulate_rlg(sr, q=2, p=5, n=100, variance=0.01_dp, alpha=0.1_dp, seed=12)
        call assert_true(size(sr%x, 1) == 100 .and. size(sr%x, 2) == 5, 'simula.rlg dimensions')
        call assert_true(count(sr%true_cluster == 0) == 10, 'simula.rlg requested contamination')
        call assert_true(all(sr%true_cluster >= 0 .and. sr%true_cluster <= 3), 'simula.rlg labels')
    end subroutine test_simulators

    subroutine make_cluster_data(x, labels)
        real(dp), intent(out) :: x(:, :) !! Deterministic two-cluster data with six extreme contaminants.
        integer, intent(out) :: labels(:) !! True labels; zero marks the six contaminants.
        integer :: i

        if (size(x, 1) /= 30 .or. size(x, 2) /= 2 .or. size(labels) /= 30) error stop 'make_cluster_data: bad shape'
        do i = 1, 12
            x(i, 1) = -5.0_dp + 0.15_dp * real(i - 6, dp)
            x(i, 2) = -4.0_dp + 0.10_dp * real(mod(i, 4) - 2, dp)
            labels(i) = 1
        end do
        do i = 13, 24
            x(i, 1) = 5.0_dp + 0.12_dp * real(i - 18, dp)
            x(i, 2) = 4.0_dp + 0.08_dp * real(mod(i, 5) - 2, dp)
            labels(i) = 2
        end do
        do i = 25, 30
            x(i, 1) = 30.0_dp + real(i, dp)
            x(i, 2) = -30.0_dp - 2.0_dp * real(i, dp)
            labels(i) = 0
        end do
    end subroutine make_cluster_data

    subroutine make_line_data(x, labels)
        real(dp), intent(out) :: x(:, :) !! Deterministic two-line data with six distant contaminants.
        integer, intent(out) :: labels(:) !! True labels; zero marks contaminants.
        integer :: i

        if (size(x, 1) /= 36 .or. size(x, 2) /= 2 .or. size(labels) /= 36) error stop 'make_line_data: bad shape'
        do i = 1, 15
            x(i, :) = [-4.0_dp + 0.4_dp * real(i, dp), -2.0_dp + 0.8_dp * real(i, dp)]
            labels(i) = 1
        end do
        do i = 16, 30
            x(i, :) = [-5.0_dp + 0.45_dp * real(i - 15, dp), 8.0_dp - 0.65_dp * real(i - 15, dp)]
            labels(i) = 2
        end do
        do i = 31, 36
            x(i, :) = [30.0_dp + real(i, dp), -20.0_dp - real(i, dp)]
            labels(i) = 0
        end do
    end subroutine make_line_data

    subroutine assert_close(actual, expected, tol, label)
        real(dp), intent(in) :: actual !! Computed scalar value.
        real(dp), intent(in) :: expected !! Reference scalar value.
        real(dp), intent(in) :: tol !! Maximum permitted absolute error.
        character(len=*), intent(in) :: label !! Human-readable assertion name printed on failure.

        if (.not. ieee_is_finite(actual) .or. abs(actual - expected) > tol) then
            write (*, '(a,2(1x,es24.16))') 'FAIL '//trim(label)//':', actual, expected
            error stop 1
        end if
    end subroutine assert_close

    subroutine assert_true(condition, label)
        logical, intent(in) :: condition !! Boolean condition that must evaluate true.
        character(len=*), intent(in) :: label !! Human-readable assertion name printed on failure.

        if (.not. condition) then
            write (*, '(a)') 'FAIL '//trim(label)
            error stop 1
        end if
    end subroutine assert_true

end program test_tclust
