program test_dbscan
    use dbscan_api
    implicit none

    real(dp), parameter :: tol = 2.0e-10_dp
    real(dp) :: xline(10, 1)
    real(dp) :: x(12, 2)
    real(dp) :: x10(10, 2)
    type(knn_result) :: knr
    type(frnn_result) :: fnr
    type(clustering_result) :: cl
    type(clustering_result) :: snc
    type(clustering_result) :: jpc
    type(optics_result) :: op
    type(hdbscan_result) :: hd
    type(dbcv_result) :: dv
    type(dendrogram_result) :: tree
    type(frnn_result) :: adj
    integer, allocatable :: comps(:)
    integer, allocatable :: fosc_cl(:)
    real(dp), allocatable :: fosc_stability(:)
    real(dp), allocatable :: fosc_constraint_score(:)
    real(dp), allocatable :: fosc_total_score(:)
    real(dp), allocatable :: score(:)
    real(dp), allocatable :: dens(:)
    real(dp), allocatable :: kd(:)
    real(dp), allocatable :: cd(:)
    real(dp), allocatable :: mr(:, :)
    real(dp), allocatable :: edges(:, :)
    logical, allocatable :: core(:)
    integer :: i
    integer :: status
    integer :: failures
    integer :: cluster10(10)
    integer :: constraints(12, 12)
    integer, parameter :: expected_db(12) = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 0, 0]
    real(dp), parameter :: expected_lof(12) = [ &
        0.966444658838_dp, 1.015469530053_dp, 1.105473972125_dp, 0.981537555829_dp, &
        0.981537555829_dp, 0.966444658838_dp, 1.015469530053_dp, 1.105473972125_dp, &
        0.981537555829_dp, 0.981537555829_dp, 18.304671027975_dp, 38.166614016032_dp]
    real(dp), parameter :: expected_hprob(12) = [ &
        0.938956389958643_dp, 0.948122906869829_dp, 0.938114369354383_dp, &
        0.948122906869829_dp, 0.943125974314807_dp, 0.013606076167856_dp, &
        0.161726355715091_dp, 0.0_dp, 0.161726355715091_dp, 0.080981722382741_dp, &
        0.0_dp, 0.0_dp]
    real(dp), parameter :: expected_glosh(12) = [ &
        0.1501634144012_dp, 0.0_dp, 0.1617263557151_dp, 0.0_dp, 0.08785965992069_dp, &
        0.1501634144012_dp, 0.0_dp, 0.1617263557151_dp, 0.0_dp, 0.08785965992069_dp, &
        0.9481229068698_dp, 0.9755590949365_dp]

    failures = 0
    do i = 1, 10
        xline(i, 1) = real(i, dp)
    end do
    x = reshape([ &
        0.00_dp, 0.00_dp, 0.10_dp, 0.00_dp, 0.00_dp, 0.12_dp, 0.12_dp, 0.10_dp, 0.20_dp, 0.05_dp, &
        3.00_dp, 3.00_dp, 3.10_dp, 3.00_dp, 3.00_dp, 3.12_dp, 3.12_dp, 3.10_dp, 3.20_dp, 3.05_dp, &
        1.50_dp, 1.50_dp, 6.00_dp, 0.00_dp], shape(x), order = [2, 1])
    x10 = x(1:10, :)

    call knn(xline, 5, knr, status)
    call check(status == 0, 'kNN status', failures)
    call check(all(knr%id(1, :) == [2, 3, 4, 5, 6]), 'kNN row 1', failures)
    call check(all(knr%id(5, :) == [4, 6, 3, 7, 2]), 'kNN row 5 tie order', failures)
    call check(all(knr%id(10, :) == [9, 8, 7, 6, 5]), 'kNN row 10', failures)
    kd = knn_dist(xline, 2, status)
    call check(status == 0 .and. abs(kd(1) - 2.0_dp) <= tol, 'kNNdist', failures)

    call frnn(xline, 2.0_dp, fnr, status)
    call check(status == 0, 'frNN status', failures)
    call check(all(fnr%id(fnr%offset(1):fnr%offset(2) - 1) == [2, 3]), 'frNN row 1', failures)
    call check(all(fnr%id(fnr%offset(5):fnr%offset(6) - 1) == [4, 6, 3, 7]), 'frNN row 5', failures)
    adj = adjacencylist_knn(knr)
    call check(adj%n == 10 .and. size(adj%id) == 50, 'adjacencylist kNN', failures)
    call comps_knn(knr, comps)
    call check(size(comps) == 10 .and. all(comps == comps(1)), 'connected components', failures)

    call dbscan(x, 0.25_dp, 3, cl, status = status)
    call check(status == 0, 'DBSCAN status', failures)
    call check(same_partition(cl%cluster, expected_db), 'DBSCAN partition', failures)
    call check(ncluster(cl) == 2 .and. nnoise(cl) == 2, 'DBSCAN cluster/noise counts', failures)
    core = is_corepoint(x, 0.25_dp, 3, status = status)
    call check(status == 0 .and. count(core) == 10, 'core points', failures)

    score = lof(x, 3, status)
    call check(status == 0, 'LOF status', failures)
    call check(maxval(abs(score - expected_lof)) < 2.0e-9_dp, 'LOF reference values', failures)
    dens = pointdensity(x, 0.25_dp, 'frequency', status)
    call check(status == 0 .and. all(dens(1:10) >= 3.0_dp), 'point density frequency', failures)

    call snn(x, 4, knr, status = status)
    call check(status == 0 .and. allocated(knr%shared), 'sNN shared counts', failures)
    call snnclust(x, 4, 2, 2, snc, status = status)
    call check(status == 0 .and. size(snc%cluster) == 12, 'sNN clustering', failures)
    call jpclust(x, 4, 2, jpc, status)
    call check(status == 0 .and. size(jpc%cluster) == 12, 'Jarvis-Patrick clustering', failures)

    call optics(x, 3, op, eps = 0.5_dp, status = status)
    call check(status == 0, 'OPTICS status', failures)
    call check(all(sort_copy(op%order) == [(i, i = 1, 12)]), 'OPTICS order permutation', failures)
    call check(abs(op%reachdist(2) - 0.12_dp) < tol, 'OPTICS reachability', failures)
    call extract_dbscan(op, 0.25_dp)
    call check(same_partition(op%cluster, expected_db), 'extractDBSCAN partition', failures)
    call optics(x, 3, op, eps = 5.0_dp, status = status)
    call extract_xi(op, 0.1_dp, .true., .true., status)
    call check(status == 0 .and. size(op%cluster) == 12, 'extractXi', failures)

    cd = coredist(x, 3, status)
    call check(status == 0 .and. abs(cd(1) - 0.12_dp) < tol, 'core distance', failures)
    mr = mrdist(x(1:3, :), 2, status)
    call check(status == 0 .and. size(mr, 1) == 3, 'mutual reachability', failures)
    edges = mst_dense(reshape([0.0_dp, 1.0_dp, sqrt(2.0_dp), 1.0_dp, 0.0_dp, 1.0_dp, &
        sqrt(2.0_dp), 1.0_dp, 0.0_dp], [3, 3]), status)
    call check(status == 0 .and. size(edges, 1) == 2, 'MST', failures)

    call hdbscan(x, 3, hd, status)
    call check(status == 0, 'HDBSCAN status', failures)
    call check(ncluster(hd) == 2 .and. nnoise(hd) == 1, 'HDBSCAN counts', failures)
    call check(same_partition(hd%cluster, [2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 2, 0]), &
        'HDBSCAN partition', failures)
    call check(maxval(abs(hd%membership_prob - expected_hprob)) < 2.0e-9_dp, &
        'HDBSCAN membership probabilities', failures)
    call extract_fosc(hd%hierarchy, 3, fosc_cl, fosc_stability, status)
    call check(status == 0 .and. same_partition(fosc_cl, hd%cluster), 'FOSC EOM extraction', failures)

    constraints = 0
    do i = 1, 5
        constraints(i, 6:10) = -1
        constraints(6:10, i) = -1
    end do
    constraints(1, 2) = 1
    constraints(2, 1) = 1
    constraints(6, 7) = 1
    constraints(7, 6) = 1
    call extract_fosc(hd%hierarchy, 3, fosc_cl, fosc_stability, status, constraints = constraints, &
        alpha = 0.5_dp, prune_unstable = .true., constraint_score = fosc_constraint_score, &
        total_score = fosc_total_score)
    call check(status == 0 .and. same_partition(fosc_cl, hd%cluster), 'FOSC constrained partition', failures)
    call check(size(fosc_constraint_score) == 2 .and. &
        maxval(abs(fosc_constraint_score - 0.5_dp)) < tol, 'FOSC constraint scores', failures)
    call check(size(fosc_total_score) == 2 .and. all(abs(fosc_total_score) < huge(1.0_dp)), &
        'FOSC mixed objective scores', failures)

    score = glosh(x, 3, status)
    call check(status == 0 .and. size(score) == 12, 'GLOSH scores', failures)
    call check(maxval(abs(score - expected_glosh)) < 2.0e-9_dp, 'GLOSH reference values', failures)
    call check_hdbscan_epsilon(failures)

    cluster10 = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2]
    call dbcv(x10, cluster10, dv, status = status)
    call check(status == 0, 'DBCV status', failures)
    call check(abs(dv%score - 0.96866265790062955_dp) < 2.0e-12_dp, 'DBCV score', failures)

    call optics(x, 3, op, eps = 0.5_dp, status = status)
    call optics_to_dendrogram(op, tree)
    call check(size(tree%height) == 11 .and. size(tree%order) == 12, 'OPTICS to dendrogram', failures)
    call dendrogram_to_reachability(tree, op)
    call check(size(op%order) == 12 .and. size(op%reachdist) == 12, 'dendrogram to reachability', failures)

    if (failures /= 0) then
        error stop 'dbscan tests failed'
    end if
    print '(a)', 'All dbscan tests passed.'

contains


    subroutine check_hdbscan_epsilon(failures)
        integer, intent(inout) :: failures !! Running number of failed assertions.
        real(dp), parameter :: xcoord(118) = [ &
            0.08_dp, 0.46_dp, 0.46_dp, 2.95_dp, 3.5_dp, 1.49_dp, &
            6.89_dp, 6.87_dp, 0.21_dp, 0.15_dp, 0.15_dp, 0.39_dp, &
            0.8_dp, 0.8_dp, 0.37_dp, 3.63_dp, 0.35_dp, 0.3_dp, &
            0.64_dp, 0.59_dp, 1.2_dp, 1.22_dp, 1.42_dp, 0.95_dp, &
            2.7_dp, 6.36_dp, 6.36_dp, 6.36_dp, 6.6_dp, 0.04_dp, &
            0.71_dp, 0.57_dp, 0.24_dp, 0.24_dp, 0.04_dp, 0.04_dp, &
            1.35_dp, 0.82_dp, 1.04_dp, 0.62_dp, 0.26_dp, 5.98_dp, &
            1.67_dp, 1.67_dp, 0.48_dp, 0.15_dp, 6.67_dp, 6.67_dp, &
            1.2_dp, 0.21_dp, 3.99_dp, 0.12_dp, 0.19_dp, 0.15_dp, &
            6.96_dp, 0.26_dp, 0.08_dp, 0.3_dp, 1.04_dp, 1.04_dp, &
            1.04_dp, 0.62_dp, 0.04_dp, 0.04_dp, 0.04_dp, 0.82_dp, &
            0.82_dp, 1.29_dp, 1.35_dp, 0.46_dp, 0.46_dp, 0.04_dp, &
            0.04_dp, 5.98_dp, 5.98_dp, 6.87_dp, 0.37_dp, 6.47_dp, &
            6.47_dp, 6.47_dp, 6.67_dp, 0.3_dp, 1.49_dp, 3.21_dp, &
            3.21_dp, 0.75_dp, 0.75_dp, 0.46_dp, 0.46_dp, 0.46_dp, &
            0.46_dp, 3.63_dp, 0.39_dp, 3.65_dp, 4.09_dp, 4.01_dp, &
            3.36_dp, 1.43_dp, 3.28_dp, 5.94_dp, 6.35_dp, 6.87_dp, &
            5.6_dp, 5.99_dp, 0.12_dp, 0.0_dp, 0.32_dp, 0.39_dp, &
            0.0_dp, 1.63_dp, 1.36_dp, 5.67_dp, 5.6_dp, 5.79_dp, &
            1.1_dp, 2.99_dp, 0.39_dp, 0.18_dp &
        ]
        real(dp), parameter :: ycoord(118) = [ &
            7.41_dp, 8.01_dp, 8.01_dp, 5.44_dp, 7.11_dp, 7.13_dp, &
            1.83_dp, 1.83_dp, 8.22_dp, 8.08_dp, 8.08_dp, 7.2_dp, &
            7.83_dp, 7.83_dp, 8.29_dp, 5.99_dp, 8.32_dp, 8.22_dp, &
            7.38_dp, 7.69_dp, 8.22_dp, 7.31_dp, 8.25_dp, 8.39_dp, &
            6.34_dp, 0.16_dp, 0.16_dp, 0.16_dp, 1.66_dp, 7.55_dp, &
            7.9_dp, 8.18_dp, 8.32_dp, 8.32_dp, 7.97_dp, 7.97_dp, &
            8.15_dp, 8.43_dp, 7.83_dp, 8.32_dp, 8.29_dp, 1.03_dp, &
            7.27_dp, 7.27_dp, 8.08_dp, 7.27_dp, 0.79_dp, 0.79_dp, &
            8.22_dp, 7.73_dp, 6.62_dp, 7.62_dp, 8.39_dp, 8.36_dp, &
            1.73_dp, 8.29_dp, 8.04_dp, 8.22_dp, 7.83_dp, 7.83_dp, &
            7.83_dp, 8.32_dp, 8.11_dp, 7.69_dp, 7.55_dp, 7.2_dp, &
            7.2_dp, 8.01_dp, 8.15_dp, 7.55_dp, 7.55_dp, 7.97_dp, &
            7.97_dp, 1.03_dp, 1.03_dp, 1.24_dp, 7.2_dp, 0.47_dp, &
            0.47_dp, 0.47_dp, 0.79_dp, 8.22_dp, 7.13_dp, 6.48_dp, &
            6.48_dp, 7.1_dp, 7.1_dp, 8.01_dp, 8.01_dp, 8.01_dp, &
            8.01_dp, 5.99_dp, 8.04_dp, 5.22_dp, 5.82_dp, 5.14_dp, &
            4.81_dp, 7.62_dp, 5.73_dp, 0.55_dp, 1.31_dp, 0.05_dp, &
            0.95_dp, 1.59_dp, 7.99_dp, 7.48_dp, 8.38_dp, 7.12_dp, &
            2.01_dp, 1.4_dp, 0.0_dp, 9.69_dp, 9.47_dp, 9.25_dp, &
            2.63_dp, 6.89_dp, 0.56_dp, 3.11_dp &
        ]
        real(dp) :: xeps(118, 2)
        type(hdbscan_result) :: hdeps
        integer :: status_eps

        xeps(:, 1) = xcoord
        xeps(:, 2) = ycoord
        call hdbscan(xeps, 3, hdeps, status_eps, cluster_selection_epsilon = 1.0_dp)
        call check(status_eps == 0, 'HDBSCAN epsilon status', failures)
        call check(ncluster(hdeps) == 5 .and. nnoise(hdeps) == 0, &
            'HDBSCAN epsilon upstream regression', failures)
    end subroutine check_hdbscan_epsilon

    subroutine check(condition, name, failures)
        logical, intent(in) :: condition !! Whether the tested condition is satisfied.
        character(len=*), intent(in) :: name !! Human-readable name of the tested behavior.
        integer, intent(inout) :: failures !! Running number of failed assertions.
        if (.not. condition) then
            failures = failures + 1
            print '(a)', 'FAIL: ' // trim(name)
        end if
    end subroutine check

    pure logical function same_partition(a, b) result(equal)
        integer, intent(in) :: a(:) !! First zero-noise cluster labeling.
        integer, intent(in) :: b(:) !! Second zero-noise cluster labeling to compare modulo positive labels.
        integer :: i
        integer :: j

        equal = size(a) == size(b)
        if (.not. equal) return
        do i = 1, size(a)
            if ((a(i) == 0) .neqv. (b(i) == 0)) then
                equal = .false.
                return
            end if
            do j = 1, size(a)
                if ((a(i) == a(j)) .neqv. (b(i) == b(j))) then
                    equal = .false.
                    return
                end if
            end do
        end do
    end function same_partition

    function sort_copy(a) result(b)
        integer, intent(in) :: a(:) !! Integer vector to copy and sort in ascending order.
        integer, allocatable :: b(:)
        integer :: i
        integer :: j
        integer :: key

        b = a
        do i = 2, size(b)
            key = b(i)
            j = i - 1
            do while (j >= 1)
                if (b(j) <= key) exit
                b(j + 1) = b(j)
                j = j - 1
            end do
            b(j + 1) = key
        end do
    end function sort_copy

end program test_dbscan
