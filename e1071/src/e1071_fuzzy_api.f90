module e1071_fuzzy_api
    use e1071_kinds, only: dp
    use e1071_rng, only: rng_state, rng_integer
    use e1071_fuzzy, only: fuzzy_cluster_model, cmeans_fit, cshell_fit
    implicit none
    private

    public :: cmeans_fit_k, cshell_fit_k

contains

    subroutine cmeans_fit_k(x, k, rng, model, weights, m, distance, iter_max, reltol, ufcl, rate_par)
        real(dp), intent(in) :: x(:, :) !! Observation-by-variable matrix from which k distinct initial centers are sampled.
        integer, intent(in) :: k !! Number of fuzzy clusters; must not exceed the number of distinct observation rows.
        type(rng_state), intent(inout) :: rng !! Mutable RNG used for center sampling and the e1071-style row permutation.
        type(fuzzy_cluster_model), intent(out) :: model !! Fitted fuzzy c-means or UFCL model with memberships restored to
        !! input order.
        real(dp), intent(in), optional :: weights(:) !! Optional nonnegative observation weights in original row order.
        real(dp), intent(in), optional :: m !! Fuzzification exponent greater than one; defaults to two in cmeans_fit.
        integer, intent(in), optional :: distance !! Dissimilarity code, fuzzy_euclidean or fuzzy_manhattan.
        integer, intent(in), optional :: iter_max !! Positive maximum update count; default is inherited from cmeans_fit.
        real(dp), intent(in), optional :: reltol !! Positive relative convergence tolerance; default is inherited from cmeans_fit.
        logical, intent(in), optional :: ufcl !! True requests the on-line UFCL update path; default is batch c-means.
        real(dp), intent(in), optional :: rate_par !! UFCL learning-rate multiplier; ignored for batch c-means.
        real(dp), allocatable :: initial(:, :)
        real(dp), allocatable :: shuffled_x(:, :)
        real(dp), allocatable :: shuffled_w(:)
        real(dp), allocatable :: membership(:, :)
        integer, allocatable :: cluster(:)
        integer, allocatable :: permutation(:)
        integer :: i

        if (k < 1 .or. k > size(x, 1)) error stop "cmeans_fit_k: invalid cluster count"
        call sample_distinct_rows(x, k, rng, initial)
        call random_permutation(size(x, 1), rng, permutation)
        allocate(shuffled_x(size(x, 1), size(x, 2)))
        do i = 1, size(x, 1)
            shuffled_x(i, :) = x(permutation(i), :)
        end do
        if (present(weights)) then
            if (size(weights) /= size(x, 1)) error stop "cmeans_fit_k: weights has wrong length"
            allocate(shuffled_w(size(weights)))
            do i = 1, size(weights)
                shuffled_w(i) = weights(permutation(i))
            end do
            call cmeans_fit(shuffled_x, initial, model, shuffled_w, m, distance, iter_max, reltol, ufcl, rate_par)
        else
            call cmeans_fit(shuffled_x, initial, model, m=m, distance=distance, iter_max=iter_max, reltol=reltol, &
                            ufcl=ufcl, rate_par=rate_par)
        end if
        membership = model%membership
        cluster = model%cluster
        do i = 1, size(permutation)
            model%membership(permutation(i), :) = membership(i, :)
            model%cluster(permutation(i)) = cluster(i)
        end do
        model%size = 0
        do i = 1, size(model%cluster)
            model%size(model%cluster(i)) = model%size(model%cluster(i)) + 1
        end do
    end subroutine cmeans_fit_k

    subroutine cshell_fit_k(x, k, rng, model, radius, m, distance, iter_max)
        real(dp), intent(in) :: x(:, :) !! Observation-by-variable matrix from which distinct initial shell centers are sampled.
        integer, intent(in) :: k !! Number of fuzzy shells; must not exceed the number of distinct observation rows.
        type(rng_state), intent(inout) :: rng !! Mutable RNG used to sample distinct starting centers.
        type(fuzzy_cluster_model), intent(out) :: model !! Fitted fuzzy c-shell model.
        real(dp), intent(in), optional :: radius(:) !! Optional initial shell radii; length must equal k when supplied.
        real(dp), intent(in), optional :: m !! Fuzzification exponent greater than one; defaults to two.
        integer, intent(in), optional :: distance !! Dissimilarity code, fuzzy_euclidean or fuzzy_manhattan.
        integer, intent(in), optional :: iter_max !! Positive maximum shell-update count; default is inherited from cshell_fit.
        real(dp), allocatable :: initial(:, :)

        if (k < 1 .or. k > size(x, 1)) error stop "cshell_fit_k: invalid shell count"
        call sample_distinct_rows(x, k, rng, initial)
        call cshell_fit(x, initial, model, radius, m, distance, iter_max)
    end subroutine cshell_fit_k

    subroutine sample_distinct_rows(x, k, rng, centers)
        real(dp), intent(in) :: x(:, :) !! Candidate data rows used for random distinct-center selection.
        integer, intent(in) :: k !! Number of distinct rows required.
        type(rng_state), intent(inout) :: rng !! Mutable RNG used to randomize candidate-row order.
        real(dp), allocatable, intent(out) :: centers(:, :) !! Selected k-by-nvar initial-center matrix with no duplicate rows.
        integer, allocatable :: order(:)
        integer :: i
        integer :: j
        integer :: nselected
        logical :: duplicate

        call random_permutation(size(x, 1), rng, order)
        allocate(centers(k, size(x, 2)))
        nselected = 0
        do i = 1, size(order)
            duplicate = .false.
            do j = 1, nselected
                if (all(abs(x(order(i), :) - centers(j, :)) <= 0.0_dp)) then
                    duplicate = .true.
                    exit
                end if
            end do
            if (duplicate) cycle
            nselected = nselected + 1
            centers(nselected, :) = x(order(i), :)
            if (nselected == k) return
        end do
        error stop "sample_distinct_rows: fewer distinct data rows than requested centers"
    end subroutine sample_distinct_rows

    subroutine random_permutation(n, rng, permutation)
        integer, intent(in) :: n !! Number of one-based row indices placed into random order.
        type(rng_state), intent(inout) :: rng !! Mutable RNG used for the Fisher-Yates permutation.
        integer, allocatable, intent(out) :: permutation(:) !! Random permutation of integers one through n.
        integer :: i
        integer :: j
        integer :: tmp

        allocate(permutation(n))
        permutation = [(i, i = 1, n)]
        do i = 1, n
            j = i + rng_integer(rng, n - i + 1) - 1
            tmp = permutation(i)
            permutation(i) = permutation(j)
            permutation(j) = tmp
        end do
    end subroutine random_permutation

end module e1071_fuzzy_api
