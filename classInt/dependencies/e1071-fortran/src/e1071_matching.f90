module e1071_matching
    use e1071_kinds, only: dp
    use e1071_rng, only: rng_state, rng_integer
    use e1071_utils, only: permutations
    use proxy, only: proxy_dist_cross, gower_cross_similarity, proxy_gower_logical, proxy_gower_factor, &
                     proxy_gower_metric, proxy_gower_ordinal
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
    implicit none
    private

    integer, parameter, public :: match_rowmax = 1
    integer, parameter, public :: match_greedy = 2
    integer, parameter, public :: match_exact = 3

    type, public :: agreement_matrix_result
        real(dp), allocatable :: diagonal(:, :)
        real(dp), allocatable :: kappa(:, :)
        real(dp), allocatable :: rand(:, :)
        real(dp), allocatable :: corrected_rand(:, :)
    end type agreement_matrix_result

    type, public :: class_agreement_result
        real(dp) :: diagonal = 0.0_dp
        real(dp) :: kappa = 0.0_dp
        real(dp) :: rand = 0.0_dp
        real(dp) :: corrected_rand = 0.0_dp
    end type class_agreement_result

    public :: class_agreement, match_classes, compare_matched_classes, match_controls_numeric, match_controls_gower
    public :: proxy_gower_logical, proxy_gower_factor, proxy_gower_metric, proxy_gower_ordinal

contains

    function class_agreement(table) result(result)
        integer, intent(in) :: table(:, :) !! Nonnegative contingency table comparing two class assignments.
        type(class_agreement_result) :: result
        real(dp), allocatable :: row_sum(:)
        real(dp), allocatable :: col_sum(:)
        real(dp) :: n
        real(dp) :: p0
        real(dp) :: pc
        real(dp) :: n2
        real(dp) :: nis2
        real(dp) :: njs2
        real(dp) :: numerator
        real(dp) :: denominator
        integer :: m
        integer :: i
        integer :: j

        if (any(table < 0)) error stop "class_agreement: contingency counts must be nonnegative"
        n = real(sum(table), dp)
        if (n <= 0.0_dp) error stop "class_agreement: empty contingency table"
        allocate(row_sum(size(table, 1)), col_sum(size(table, 2)))
        row_sum = real(sum(table, dim=2), dp)
        col_sum = real(sum(table, dim=1), dp)
        m = min(size(table, 1), size(table, 2))
        p0 = 0.0_dp
        pc = 0.0_dp
        do i = 1, m
            p0 = p0 + real(table(i, i), dp) / n
            pc = pc + (row_sum(i) / n) * (col_sum(i) / n)
        end do
        n2 = n * (n - 1.0_dp) / 2.0_dp
        if (n2 > 0.0_dp) then
            result%rand = 1.0_dp
            result%rand = result%rand + (sum(real(table, dp)**2) &
                - 0.5_dp * (sum(row_sum**2) + sum(col_sum**2))) / n2
            nis2 = 0.0_dp
            do i = 1, size(row_sum)
                if (row_sum(i) > 1.0_dp) nis2 = nis2 + row_sum(i) * (row_sum(i) - 1.0_dp) / 2.0_dp
            end do
            njs2 = 0.0_dp
            do j = 1, size(col_sum)
                if (col_sum(j) > 1.0_dp) njs2 = njs2 + col_sum(j) * (col_sum(j) - 1.0_dp) / 2.0_dp
            end do
            numerator = 0.0_dp
            do j = 1, size(table, 2)
                do i = 1, size(table, 1)
                    if (table(i, j) > 1) numerator = numerator + real(table(i, j) * (table(i, j) - 1) / 2, dp)
                end do
            end do
            numerator = numerator - nis2 * njs2 / n2
            denominator = 0.5_dp * (nis2 + njs2) - nis2 * njs2 / n2
            if (abs(denominator) > tiny(1.0_dp)) result%corrected_rand = numerator / denominator
        end if
        result%diagonal = p0
        if (abs(1.0_dp - pc) > tiny(1.0_dp)) result%kappa = (p0 - pc) / (1.0_dp - pc)
    end function class_agreement

    subroutine match_classes(table, mapping, method, iterations, max_exact, rng)
        integer, intent(in) :: table(:, :) !! Contingency table whose rows are reference classes and columns are classes to relabel.
        integer, allocatable, intent(out) :: mapping(:) !! One-based matched column index for every row class.
        integer, intent(in), optional :: method !! match_rowmax, match_greedy, or match_exact; defaults to rowmax.
        integer, intent(in), optional :: iterations !! Number of randomized greedy attempts; defaults to one.
        integer, intent(in), optional :: max_exact !! Largest unmatched square problem allowed for exhaustive permutation;
        !! defaults to nine.
        type(rng_state), intent(inout), optional :: rng !! Optional RNG used by randomized greedy matching.
        integer :: use_method
        integer :: use_iter
        integer :: use_max_exact
        integer :: n
        integer :: i
        integer :: trial
        integer :: row
        integer :: col
        integer :: best_col
        integer, allocatable :: candidate(:)
        integer, allocatable :: best(:)
        integer, allocatable :: available(:)
        integer, allocatable :: perm(:, :)
        real(dp) :: score
        real(dp) :: best_score

        use_method = match_rowmax
        if (present(method)) use_method = method
        use_iter = 1
        if (present(iterations)) use_iter = iterations
        use_max_exact = 9
        if (present(max_exact)) use_max_exact = max_exact
        if (use_method == match_rowmax) then
            allocate(mapping(size(table, 1)))
            do row = 1, size(table, 1)
                mapping(row) = maxloc(table(row, :), dim=1)
            end do
            return
        end if
        if (size(table, 1) /= size(table, 2)) error stop "match_classes: unique matching requires a square table"
        n = size(table, 1)
        if (use_method == match_exact .and. n <= use_max_exact) then
            perm = permutations(n)
            best_score = -1.0_dp
            allocate(best(n))
            do trial = 1, size(perm, 1)
                score = 0.0_dp
                do row = 1, n
                    score = score + real(table(row, perm(trial, row)), dp)
                end do
                if (score > best_score) then
                    best_score = score
                    best = perm(trial, :)
                end if
            end do
            mapping = best
            return
        end if

        allocate(candidate(n), best(n), available(n))
        best_score = -1.0_dp
        do trial = 1, max(1, use_iter)
            available = [(i, i = 1, n)]
            candidate = 0
            do i = 1, n
                if (present(rng)) then
                    row = nth_unmatched(candidate, rng_integer(rng, n - i + 1))
                else
                    row = nth_unmatched(candidate, 1)
                end if
                best_col = 0
                do col = 1, n
                    if (.not. any(available(:n - i + 1) == col)) cycle
                    if (best_col == 0) then
                        best_col = col
                    else if (table(row, col) > table(row, best_col)) then
                        best_col = col
                    end if
                end do
                candidate(row) = best_col
                call remove_available(available, n - i + 1, best_col)
            end do
            score = 0.0_dp
            do row = 1, n
                score = score + real(table(row, candidate(row)), dp)
            end do
            if (score > best_score) then
                best_score = score
                best = candidate
            end if
        end do
        mapping = best
    end subroutine match_classes


    subroutine compare_matched_classes(x, result, y, method, iterations, max_exact, rng)
        integer, intent(in) :: x(:, :) !! First observation-by-clustering label matrix; columns are cluster solutions to compare.
        type(agreement_matrix_result), intent(out) :: result !! Pairwise diagonal, kappa, Rand, and adjusted-Rand agreement
        !! matrices.
        integer, intent(in), optional :: y(:, :) !! Optional second label matrix; absence compares distinct column pairs within x.
        integer, intent(in), optional :: method !! Matching rule passed to match_classes; defaults to match_rowmax.
        integer, intent(in), optional :: iterations !! Greedy matching repetitions forwarded to match_classes; defaults to one.
        integer, intent(in), optional :: max_exact !! Largest exact square matching problem; defaults to nine.
        type(rng_state), intent(inout), optional :: rng !! Optional RNG used by randomized greedy class matching.
        integer, allocatable :: table(:, :)
        integer, allocatable :: mapping(:)
        integer, allocatable :: aligned(:, :)
        type(class_agreement_result) :: one
        real(dp) :: nan
        integer :: nx
        integer :: ny
        integer :: i
        integer :: j

        nan = ieee_value(0.0_dp, ieee_quiet_nan)
        nx = size(x, 2)
        if (present(y)) then
            if (size(y, 1) /= size(x, 1)) error stop "compare_matched_classes: x/y row mismatch"
            ny = size(y, 2)
        else
            ny = nx
        end if
        allocate(result%diagonal(nx, ny), result%kappa(nx, ny), result%rand(nx, ny), result%corrected_rand(nx, ny))
        result%diagonal = nan
        result%kappa = nan
        result%rand = nan
        result%corrected_rand = nan
        if (present(y)) then
            do i = 1, nx
                do j = 1, ny
                    call contingency_table(x(:, i), y(:, j), table)
                    call match_classes(table, mapping, method, iterations, max_exact, rng)
                    call aligned_table(table, mapping, aligned)
                    one = class_agreement(aligned)
                    call store_agreement(result, i, j, one)
                end do
            end do
        else
            do i = 1, nx - 1
                do j = i + 1, nx
                    call contingency_table(x(:, i), x(:, j), table)
                    call match_classes(table, mapping, method, iterations, max_exact, rng)
                    call aligned_table(table, mapping, aligned)
                    one = class_agreement(aligned)
                    call store_agreement(result, i, j, one)
                end do
            end do
        end if
    end subroutine compare_matched_classes

    subroutine contingency_table(a, b, table)
        integer, intent(in) :: a(:) !! First integer class assignment vector.
        integer, intent(in) :: b(:) !! Second integer class assignment vector with the same observation count as a.
        integer, allocatable, intent(out) :: table(:, :) !! Contingency counts over sorted distinct labels of a and b.
        integer, allocatable :: la(:)
        integer, allocatable :: lb(:)
        integer :: i
        integer :: ia
        integer :: ib

        if (size(a) /= size(b)) error stop "contingency_table: label vectors have different lengths"
        call unique_labels_local(a, la)
        call unique_labels_local(b, lb)
        allocate(table(size(la), size(lb)))
        table = 0
        do i = 1, size(a)
            ia = find_label_local(la, a(i))
            ib = find_label_local(lb, b(i))
            table(ia, ib) = table(ia, ib) + 1
        end do
    end subroutine contingency_table

    subroutine aligned_table(table, mapping, aligned)
        integer, intent(in) :: table(:, :) !! Original contingency table before class-column matching.
        integer, intent(in) :: mapping(:) !! One-based source column selected for each aligned output column.
        integer, allocatable, intent(out) :: aligned(:, :) !! Table with columns reordered or repeated according to mapping.
        integer :: j

        allocate(aligned(size(table, 1), size(mapping)))
        do j = 1, size(mapping)
            if (mapping(j) < 1 .or. mapping(j) > size(table, 2)) error stop "aligned_table: invalid mapping"
            aligned(:, j) = table(:, mapping(j))
        end do
    end subroutine aligned_table

    subroutine store_agreement(result, i, j, value)
        type(agreement_matrix_result), intent(inout) :: result !! Agreement matrices receiving one pairwise comparison.
        integer, intent(in) :: i !! Row index of the comparison result to update.
        integer, intent(in) :: j !! Column index of the comparison result to update.
        type(class_agreement_result), intent(in) :: value !! Scalar agreement measures copied into the four result matrices.

        result%diagonal(i, j) = value%diagonal
        result%kappa(i, j) = value%kappa
        result%rand(i, j) = value%rand
        result%corrected_rand(i, j) = value%corrected_rand
    end subroutine store_agreement

    subroutine unique_labels_local(x, labels)
        integer, intent(in) :: x(:) !! Integer labels from which sorted unique values are extracted.
        integer, allocatable, intent(out) :: labels(:) !! Sorted unique label values present in x.
        integer, allocatable :: work(:)
        integer :: i
        integer :: j
        integer :: n
        integer :: key

        if (size(x) == 0) then
            allocate(labels(0))
            return
        end if
        work = x
        do i = 2, size(work)
            key = work(i)
            j = i - 1
            do while (j >= 1)
                if (work(j) <= key) exit
                work(j + 1) = work(j)
                j = j - 1
            end do
            work(j + 1) = key
        end do
        n = 1
        do i = 2, size(work)
            if (work(i) /= work(n)) then
                n = n + 1
                work(n) = work(i)
            end if
        end do
        allocate(labels(n))
        labels = work(:n)
    end subroutine unique_labels_local

    pure function find_label_local(labels, value) result(index_value)
        integer, intent(in) :: labels(:) !! Sorted distinct labels searched for an exact integer match.
        integer, intent(in) :: value !! Label value whose one-based position is returned.
        integer :: index_value
        integer :: i

        index_value = 0
        do i = 1, size(labels)
            if (labels(i) == value) then
                index_value = i
                return
            end if
        end do
    end function find_label_local

    subroutine match_controls_numeric(x, is_case, is_control, controls, replace, method)
        real(dp), intent(in) :: x(:, :) !! Observation-by-variable numeric matching covariates.
        logical, intent(in) :: is_case(:) !! Row mask identifying observations that require matched controls.
        logical, intent(in) :: is_control(:) !! Row mask identifying candidate control observations.
        integer, allocatable, intent(out) :: controls(:) !! Original one-based row index of the selected control for each case
        !! in case order.
        logical, intent(in), optional :: replace !! If true, controls may be reused; default false.
        character(len=*), intent(in), optional :: method !! proxy distance name used for matching; defaults to Euclidean.
        real(dp), allocatable :: cases_x(:, :)
        real(dp), allocatable :: controls_x(:, :)
        real(dp), allocatable :: distance(:, :)
        integer, allocatable :: case_index(:)
        integer, allocatable :: control_index(:)
        logical, allocatable :: available(:)
        logical :: do_replace
        character(len=:), allocatable :: use_method
        integer :: ncases
        integer :: ncontrols
        integer :: i
        integer :: j
        integer :: best
        integer :: status

        if (size(is_case) /= size(x, 1) .or. size(is_control) /= size(x, 1)) then
            error stop "match_controls_numeric: mask length mismatch"
        end if
        ncases = count(is_case)
        ncontrols = count(is_control)
        if (ncases < 1 .or. ncontrols < 1) error stop "match_controls_numeric: cases and controls must be nonempty"
        do_replace = .false.
        if (present(replace)) do_replace = replace
        if (.not. do_replace .and. ncontrols < ncases) then
            error stop "match_controls_numeric: insufficient controls without replacement"
        end if
        use_method = "euclidean"
        if (present(method)) use_method = trim(adjustl(method))
        allocate(cases_x(ncases, size(x, 2)), controls_x(ncontrols, size(x, 2)))
        allocate(case_index(ncases), control_index(ncontrols), available(ncontrols), controls(ncases))
        i = 0
        j = 0
        do best = 1, size(x, 1)
            if (is_case(best)) then
                i = i + 1
                case_index(i) = best
                cases_x(i, :) = x(best, :)
            end if
            if (is_control(best)) then
                j = j + 1
                control_index(j) = best
                controls_x(j, :) = x(best, :)
            end if
        end do
        call proxy_dist_cross(cases_x, controls_x, use_method, distance, status=status)
        if (status /= 0) error stop "match_controls_numeric: proxy distance evaluation failed"
        available = .true.
        do i = 1, ncases
            best = 0
            do j = 1, ncontrols
                if (.not. available(j)) cycle
                if (best == 0 .or. distance(i, j) < distance(i, best)) best = j
            end do
            controls(i) = control_index(best)
            if (.not. do_replace) available(best) = .false.
        end do
    end subroutine match_controls_numeric


    subroutine match_controls_gower(x, types, is_case, is_control, controls, replace, rng)
        real(dp), intent(in) :: x(:, :) !! Mixed matching covariates numerically encoded for proxy Gower semantics.
        integer, intent(in) :: types(:) !! proxy Gower type code for every column of x, including logical/factor/metric/ordinal.
        logical, intent(in) :: is_case(:) !! Row mask identifying observations that require matched controls.
        logical, intent(in) :: is_control(:) !! Row mask identifying candidate control observations.
        integer, allocatable, intent(out) :: controls(:) !! Original one-based row index of the selected control for each case.
        logical, intent(in), optional :: replace !! If true, controls may be reused; default false, matching e1071.
        type(rng_state), intent(inout), optional :: rng !! Optional RNG for exact distance ties; absence chooses the first minimum.
        real(dp), allocatable :: cases_x(:, :)
        real(dp), allocatable :: controls_x(:, :)
        real(dp), allocatable :: similarity(:, :)
        real(dp), allocatable :: distance(:)
        integer, allocatable :: control_index(:)
        logical, allocatable :: available(:)
        integer, allocatable :: tied(:)
        logical :: do_replace
        integer :: ncases
        integer :: ncontrols
        integer :: i
        integer :: j
        integer :: row
        integer :: icase
        integer :: icontrol
        integer :: ntied
        integer :: pick
        real(dp) :: best

        if (size(types) /= size(x, 2)) error stop "match_controls_gower: type count does not match columns"
        if (size(is_case) /= size(x, 1) .or. size(is_control) /= size(x, 1)) then
            error stop "match_controls_gower: mask length mismatch"
        end if
        ncases = count(is_case)
        ncontrols = count(is_control)
        if (ncases < 1 .or. ncontrols < 1) error stop "match_controls_gower: cases and controls must be nonempty"
        do_replace = .false.
        if (present(replace)) do_replace = replace
        if (.not. do_replace .and. ncontrols < ncases) then
            error stop "match_controls_gower: insufficient controls without replacement"
        end if
        allocate(cases_x(ncases, size(x, 2)), controls_x(ncontrols, size(x, 2)))
        allocate(control_index(ncontrols), available(ncontrols), controls(ncases), tied(ncontrols), distance(ncontrols))
        icase = 0
        icontrol = 0
        do row = 1, size(x, 1)
            if (is_case(row)) then
                icase = icase + 1
                cases_x(icase, :) = x(row, :)
            end if
            if (is_control(row)) then
                icontrol = icontrol + 1
                controls_x(icontrol, :) = x(row, :)
                control_index(icontrol) = row
            end if
        end do
        call gower_cross_similarity(cases_x, controls_x, types, similarity)
        available = .true.
        do i = 1, ncases
            distance = 1.0_dp - similarity(i, :)
            best = huge(1.0_dp)
            do j = 1, ncontrols
                if (.not. available(j) .or. .not. ieee_is_finite(distance(j))) cycle
                best = min(best, distance(j))
            end do
            if (best >= huge(1.0_dp)) error stop "match_controls_gower: no finite eligible control distance"
            ntied = 0
            do j = 1, ncontrols
                if (.not. available(j) .or. .not. ieee_is_finite(distance(j))) cycle
                if (abs(distance(j) - best) <= 0.0_dp) then
                    ntied = ntied + 1
                    tied(ntied) = j
                end if
            end do
            pick = 1
            if (ntied > 1 .and. present(rng)) pick = rng_integer(rng, ntied)
            j = tied(pick)
            controls(i) = control_index(j)
            if (.not. do_replace) available(j) = .false.
        end do
    end subroutine match_controls_gower

    pure function nth_unmatched(candidate, rank) result(row)
        integer, intent(in) :: candidate(:) !! Partial row-to-column mapping where zero marks an unmatched row.
        integer, intent(in) :: rank !! One-based rank among currently unmatched rows.
        integer :: row
        integer :: seen
        integer :: i

        seen = 0
        row = 0
        do i = 1, size(candidate)
            if (candidate(i) /= 0) cycle
            seen = seen + 1
            if (seen == rank) then
                row = i
                return
            end if
        end do
    end function nth_unmatched

    subroutine remove_available(available, nactive, value)
        integer, intent(inout) :: available(:) !! Prefix list of currently available column indices, compacted in place.
        integer, intent(in) :: nactive !! Number of valid entries in the available prefix before removal.
        integer, intent(in) :: value !! Column index to remove from the active prefix.
        integer :: i
        integer :: pos

        pos = 0
        do i = 1, nactive
            if (available(i) == value) then
                pos = i
                exit
            end if
        end do
        if (pos == 0) return
        do i = pos, nactive - 1
            available(i) = available(i + 1)
        end do
    end subroutine remove_available

end module e1071_matching
