module tclust_grid
    use tclust_kinds, only : dp
    use tclust_types, only : tclust_result, tclust_ic_result, tclust_ic_solution_result
    use tclust_core, only : tclust
    use tclust_metrics, only : rand_index_labels, fowlkes_mallows_labels
    use tclust_types, only : rand_index_result, fm_index_result
    implicit none
    private

    public :: tclust_ic
    public :: tclust_ic_solutions

contains

    subroutine tclust_ic(x, k_values, c_values, result, alpha, nstart, niter1, niter2, nkeep, restriction, seed)
        real(dp), intent(in) :: x(:, :) !! Numeric data matrix used for all information-criterion fits.
        integer, intent(in) :: k_values(:) !! Positive component counts forming the rows of the criterion grids.
        real(dp), intent(in) :: c_values(:) !! Restriction factors at least one forming the grid columns.
        type(tclust_ic_result), intent(out) :: result !! CLA, BIC, ICL values and associated HARD/MIXT classifications.
        real(dp), intent(in), optional :: alpha !! Common trimming fraction for all fits; default 0.05.
        integer, intent(in), optional :: nstart !! Number of random starts per fit; default 50 in this grid helper.
        integer, intent(in), optional :: niter1 !! Initial concentration steps per fit; default 3.
        integer, intent(in), optional :: niter2 !! Refinement steps per fit; default 20.
        integer, intent(in), optional :: nkeep !! Number of retained starts per fit; default 5.
        character(len=*), intent(in), optional :: restriction !! Scatter restriction passed to tclust; default 'eigen'.
        integer, intent(in), optional :: seed !! Optional base seed; individual grid fits use deterministic offsets.
        type(tclust_result), allocatable :: hard
        type(tclust_result), allocatable :: mixt
        real(dp) :: a
        integer :: ns
        integer :: ni1
        integer :: ni2
        integer :: nk
        integer :: base_seed
        integer :: i
        integer :: j
        character(len=8) :: rest

        if (any(k_values < 1)) error stop 'tclust_ic: k_values must be positive'
        if (any(c_values < 1.0_dp)) error stop 'tclust_ic: c_values must be at least one'
        a = 0.05_dp
        if (present(alpha)) a = alpha
        ns = 50
        if (present(nstart)) ns = nstart
        ni1 = 3
        if (present(niter1)) ni1 = niter1
        ni2 = 20
        if (present(niter2)) ni2 = niter2
        nk = min(5, ns)
        if (present(nkeep)) nk = nkeep
        rest = 'eigen'
        if (present(restriction)) rest = trim(adjustl(restriction))
        allocate(hard, mixt)
        base_seed = 834927
        if (present(seed)) base_seed = seed

        result%n = size(x, 1)
        allocate(result%k_values(size(k_values)), result%c_values(size(c_values)), &
                 result%clacla(size(k_values), size(c_values)), result%mixmix(size(k_values), size(c_values)), &
                 result%mixcla(size(k_values), size(c_values)), &
                 result%idxcla(size(x, 1), size(k_values), size(c_values)), &
                 result%idxmix(size(x, 1), size(k_values), size(c_values)))
        result%k_values = k_values
        result%c_values = c_values
        do i = 1, size(k_values)
            do j = 1, size(c_values)
                call tclust(x, k_values(i), hard, alpha=a, nstart=ns, niter1=ni1, niter2=ni2, nkeep=nk, &
                            restriction=rest, restr_fact=c_values(j), opt='HARD', &
                            seed=base_seed + 10007 * i + 101 * j)
                call tclust(x, k_values(i), mixt, alpha=a, nstart=ns, niter1=ni1, niter2=ni2, nkeep=nk, &
                            restriction=rest, restr_fact=c_values(j), opt='MIXT', &
                            seed=base_seed + 20011 * i + 211 * j)
                result%clacla(i, j) = hard%clacla
                result%mixmix(i, j) = mixt%mixmix
                result%mixcla(i, j) = mixt%mixcla
                result%idxcla(:, i, j) = hard%cluster
                result%idxmix(:, i, j) = mixt%cluster
            end do
        end do
    end subroutine tclust_ic

    subroutine tclust_ic_solutions(grid, result, criterion, nsol, index_name, threshold)
        type(tclust_ic_result), intent(in) :: grid !! Previously computed tclust information-criterion grid.
        type(tclust_ic_solution_result), intent(out) :: result !! Ranked grid solutions plus adjacent-partition stability matrices.
        character(len=*), intent(in), optional :: criterion !! Criterion to rank: MIXMIX, MIXCLA, or CLACLA; default MIXMIX.
        integer, intent(in), optional :: nsol !! Maximum number of ranked solutions returned; default five.
        character(len=*), intent(in), optional :: index_name !! Stability index 'Rand' or 'FM'; default Rand.
        real(dp), intent(in), optional :: threshold !! Adjacent-partition stability cutoff; default 0.7.
        real(dp), allocatable :: values(:)
        integer, allocatable :: ki(:)
        integer, allocatable :: ci(:)
        logical, allocatable :: used(:)
        type(rand_index_result) :: ri
        type(fm_index_result) :: fm
        character(len=8) :: crit
        character(len=4) :: idx
        real(dp) :: cut
        real(dp) :: score
        real(dp) :: best
        integer :: want
        integer :: nk
        integer :: nc
        integer :: total
        integer :: i
        integer :: j
        integer :: z
        integer :: m

        nk = size(grid%k_values)
        nc = size(grid%c_values)
        crit = 'MIXMIX'
        if (present(criterion)) crit = trim(adjustl(criterion))
        idx = 'Rand'
        if (present(index_name)) idx = trim(adjustl(index_name))
        cut = 0.7_dp
        if (present(threshold)) cut = threshold
        want = 5
        if (present(nsol)) want = nsol
        want = min(want, nk * nc)
        if (want < 1) error stop 'tclust_ic_solutions: nsol must be positive'
        if (crit /= 'MIXMIX' .and. crit /= 'MIXCLA' .and. crit /= 'CLACLA') then
            error stop 'tclust_ic_solutions: criterion must be MIXMIX, MIXCLA, or CLACLA'
        end if
        if (idx /= 'Rand' .and. idx /= 'FM') error stop 'tclust_ic_solutions: index_name must be Rand or FM'

        allocate(result%ari_mix(nk, max(0, nc - 1)), result%ari_cla(nk, max(0, nc - 1)))
        result%ari_mix = 1.0_dp
        result%ari_cla = 1.0_dp
        if (nc > 1) then
            do i = 1, nk
                do j = 2, nc
                    if (idx == 'Rand') then
                        call rand_index_labels(grid%idxmix(:, i, j - 1), grid%idxmix(:, i, j), ri)
                        result%ari_mix(i, j - 1) = ri%adjusted_rand
                        call rand_index_labels(grid%idxcla(:, i, j - 1), grid%idxcla(:, i, j), ri)
                        result%ari_cla(i, j - 1) = ri%adjusted_rand
                    else
                        call fowlkes_mallows_labels(grid%idxmix(:, i, j - 1), grid%idxmix(:, i, j), fm)
                        result%ari_mix(i, j - 1) = fm%adjusted
                        call fowlkes_mallows_labels(grid%idxcla(:, i, j - 1), grid%idxcla(:, i, j), fm)
                        result%ari_cla(i, j - 1) = fm%adjusted
                    end if
                end do
            end do
        end if

        total = nk * nc
        allocate(values(total), ki(total), ci(total), used(total))
        m = 0
        do i = 1, nk
            do j = 1, nc
                m = m + 1
                ki(m) = i
                ci(m) = j
                select case (crit)
                case ('MIXMIX')
                    values(m) = grid%mixmix(i, j)
                case ('MIXCLA')
                    values(m) = grid%mixcla(i, j)
                case default
                    values(m) = grid%clacla(i, j)
                end select
            end do
        end do
        used = .false.
        result%nsolutions = want
        allocate(result%best_k(want), result%best_c(want), result%spurious(want), result%ari_best(want, 1))
        do z = 1, want
            best = huge(1.0_dp)
            m = 0
            do i = 1, total
                if (.not. used(i) .and. values(i) < best) then
                    best = values(i)
                    m = i
                end if
            end do
            if (m == 0) error stop 'tclust_ic_solutions: internal ranking failure'
            used(m) = .true.
            result%best_k(z) = grid%k_values(ki(m))
            result%best_c(z) = grid%c_values(ci(m))
            score = 1.0_dp
            if (ci(m) > 1) then
                if (crit == 'CLACLA') then
                    score = result%ari_cla(ki(m), ci(m) - 1)
                else
                    score = result%ari_mix(ki(m), ci(m) - 1)
                end if
            end if
            result%ari_best(z, 1) = score
            result%spurious(z) = score < cut
        end do
    end subroutine tclust_ic_solutions

end module tclust_grid
