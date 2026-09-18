module tclust_rlg
    use tclust_kinds, only : dp
    use tclust_types, only : rlg_result
    use tclust_rng, only : set_tclust_seed, sample_index
    use tclust_linalg, only : covariance_matrix, symmetric_eigen, sort_real_with_index
    implicit none
    private

    public :: rlg
    public :: rlg_refine

contains

    subroutine rlg(x, dimensions, result, alpha, nstart, niter1, niter2, nkeep, seed)
        real(dp), intent(in) :: x(:, :) !! Numeric data matrix with observations in rows and variables in columns.
        integer, intent(in) :: dimensions(:) !! Intrinsic affine dimension for each cluster; each value must be in [0,p-1].
        type(rlg_result), intent(out) :: result !! Best robust linear-grouping fit with centers and affine bases.
        real(dp), intent(in), optional :: alpha !! Fraction of observations trimmed; default 0.05.
        integer, intent(in), optional :: nstart !! Number of random initializations; default 500.
        integer, intent(in), optional :: niter1 !! Concentration steps for each random start; default 3.
        integer, intent(in), optional :: niter2 !! Refinement steps for retained starts; default 20.
        integer, intent(in), optional :: nkeep !! Number of phase-one starts retained; default 5.
        integer, intent(in), optional :: seed !! Optional deterministic seed for the intrinsic Fortran RNG.
        integer, allocatable :: labels(:, :)
        integer, allocatable :: order(:)
        integer, allocatable :: lab(:)
        real(dp), allocatable :: obj(:)
        real(dp), allocatable :: centers(:, :)
        real(dp), allocatable :: basis(:, :, :)
        real(dp) :: a
        real(dp) :: current_obj
        real(dp) :: best_obj
        integer :: n
        integer :: p
        integer :: k
        integer :: ns
        integer :: ni1
        integer :: ni2
        integer :: nk
        integer :: no_trim
        integer :: maxd
        integer :: j
        logical :: found

        n = size(x, 1)
        p = size(x, 2)
        k = size(dimensions)
        if (n < 1 .or. p < 1 .or. k < 1) error stop 'rlg: x and dimensions must be nonempty'
        if (any(dimensions < 0) .or. any(dimensions >= p)) error stop 'rlg: each intrinsic dimension must be in [0,p-1]'
        a = 0.05_dp
        if (present(alpha)) a = alpha
        if (a < 0.0_dp .or. a >= 1.0_dp) error stop 'rlg: alpha must be in [0,1)'
        no_trim = n - floor(a * real(n, dp))
        if (no_trim < k) error stop 'rlg: fewer retained observations than clusters'
        ns = 500
        if (present(nstart)) ns = nstart
        ni1 = 3
        if (present(niter1)) ni1 = niter1
        ni2 = 20
        if (present(niter2)) ni2 = niter2
        nk = 5
        if (present(nkeep)) nk = nkeep
        if (ns < 1 .or. nk < 1 .or. nk > ns .or. ni1 < 0 .or. ni2 < 0) error stop 'rlg: invalid iteration counts'
        if (present(seed)) call set_tclust_seed(seed)
        maxd = max(1, maxval(dimensions))
        allocate(labels(n, ns), order(ns), obj(ns), lab(n), centers(k, p), basis(p, maxd, k))

        do j = 1, ns
            call initialize_affine(x, dimensions, centers, basis)
            lab = 0
            call rlg_steps(x, dimensions, centers, basis, lab, no_trim, ni1, current_obj)
            labels(:, j) = lab
            obj(j) = current_obj
        end do
        call sort_real_with_index(obj, order, ascending=.true.)

        best_obj = huge(1.0_dp)
        found = .false.
        do j = 1, nk
            lab = labels(:, order(j))
            if (.not. all_clusters_present(lab, k)) cycle
            call fit_affine_from_labels(x, dimensions, lab, centers, basis)
            call rlg_steps(x, dimensions, centers, basis, lab, no_trim, ni2, current_obj)
            if (current_obj < best_obj) then
                best_obj = current_obj
                found = .true.
                call fill_rlg_result(x, dimensions, centers, basis, lab, a, current_obj, result)
            end if
        end do
        if (.not. found) error stop 'rlg: no retained initialization contained every requested cluster'
    end subroutine rlg

    subroutine rlg_refine(x, dimensions, cluster, result, alpha, niter)
        real(dp), intent(in) :: x(:, :) !! Numeric data matrix with observations in rows and variables in columns.
        integer, intent(in) :: dimensions(:) !! Intrinsic affine dimension for each cluster.
        integer, intent(in) :: cluster(:) !! Starting labels of length n; zero denotes trimmed observations.
        type(rlg_result), intent(out) :: result !! Refined robust linear-grouping fit.
        real(dp), intent(in), optional :: alpha !! Trimming fraction; inferred from zero labels when absent.
        integer, intent(in), optional :: niter !! Maximum refinement concentration steps; default 20.
        real(dp), allocatable :: centers(:, :)
        real(dp), allocatable :: basis(:, :, :)
        integer, allocatable :: lab(:)
        real(dp) :: a
        real(dp) :: obj
        integer :: nit
        integer :: no_trim
        integer :: maxd

        if (size(cluster) /= size(x, 1)) error stop 'rlg_refine: cluster length mismatch'
        if (maxval(cluster) > size(dimensions)) error stop 'rlg_refine: cluster label exceeds dimensions length'
        a = real(count(cluster == 0), dp) / real(size(cluster), dp)
        if (present(alpha)) a = alpha
        no_trim = size(cluster) - floor(a * real(size(cluster), dp))
        nit = 20
        if (present(niter)) nit = niter
        maxd = max(1, maxval(dimensions))
        allocate(centers(size(dimensions), size(x, 2)), basis(size(x, 2), maxd, size(dimensions)), lab(size(cluster)))
        lab = cluster
        call fit_affine_from_labels(x, dimensions, lab, centers, basis)
        call rlg_steps(x, dimensions, centers, basis, lab, no_trim, nit, obj)
        call fill_rlg_result(x, dimensions, centers, basis, lab, a, obj, result)
    end subroutine rlg_refine

    subroutine initialize_affine(x, dimensions, centers, basis)
        real(dp), intent(in) :: x(:, :) !! Data matrix used to sample initial points and affine subspaces.
        integer, intent(in) :: dimensions(:) !! Intrinsic dimension requested for each cluster.
        real(dp), intent(out) :: centers(:, :) !! Initial k-by-p location matrix.
        real(dp), intent(out) :: basis(:, :, :) !! Initial p-by-max(d,1)-by-k orthonormal direction arrays.
        real(dp), allocatable :: sub(:, :)
        real(dp), allocatable :: cov(:, :)
        real(dp), allocatable :: vals(:)
        real(dp), allocatable :: vecs(:, :)
        integer :: j
        integer :: i
        integer :: p
        integer :: d

        p = size(x, 2)
        basis = 0.0_dp
        allocate(cov(p, p), vals(p), vecs(p, p))
        do j = 1, size(dimensions)
            d = dimensions(j)
            if (d == 0) then
                centers(j, :) = x(sample_index(size(x, 1)), :)
            else
                allocate(sub(d + 1, p))
                do i = 1, d + 1
                    sub(i, :) = x(sample_index(size(x, 1)), :)
                end do
                centers(j, :) = sum(sub, dim=1) / real(d + 1, dp)
                call covariance_matrix(sub, cov, unbiased=.false.)
                call symmetric_eigen(cov, vals, vecs)
                do i = 1, d
                    basis(:, i, j) = vecs(:, p - i + 1)
                end do
                deallocate(sub)
            end if
        end do
    end subroutine initialize_affine

    subroutine fit_affine_from_labels(x, dimensions, labels, centers, basis)
        real(dp), intent(in) :: x(:, :) !! Data matrix used to fit each currently assigned affine subspace.
        integer, intent(in) :: dimensions(:) !! Intrinsic dimension for each cluster.
        integer, intent(in) :: labels(:) !! Current labels; zero denotes trimmed observations.
        real(dp), intent(out) :: centers(:, :) !! Fitted k-by-p locations.
        real(dp), intent(out) :: basis(:, :, :) !! Fitted p-by-max(d,1)-by-k direction arrays.
        real(dp), allocatable :: sub(:, :)
        real(dp), allocatable :: cov(:, :)
        real(dp), allocatable :: vals(:)
        real(dp), allocatable :: vecs(:, :)
        integer :: j
        integer :: i
        integer :: m
        integer :: p
        integer :: pos

        p = size(x, 2)
        allocate(cov(p, p), vals(p), vecs(p, p))
        basis = 0.0_dp
        centers = 0.0_dp
        do j = 1, size(dimensions)
            m = count(labels == j)
            if (m == 0) cycle
            allocate(sub(m, p))
            pos = 0
            do i = 1, size(x, 1)
                if (labels(i) == j) then
                    pos = pos + 1
                    sub(pos, :) = x(i, :)
                end if
            end do
            centers(j, :) = sum(sub, dim=1) / real(m, dp)
            if (dimensions(j) > 0 .and. m > 1) then
                call covariance_matrix(sub, cov, unbiased=.false.)
                call symmetric_eigen(cov, vals, vecs)
                do i = 1, min(dimensions(j), p)
                    basis(:, i, j) = vecs(:, p - i + 1)
                end do
            end if
            deallocate(sub)
        end do
    end subroutine fit_affine_from_labels

    subroutine rlg_steps(x, dimensions, centers, basis, labels, no_trim, niter, obj)
        real(dp), intent(in) :: x(:, :) !! Data matrix assigned by squared orthogonal affine-subspace distance.
        integer, intent(in) :: dimensions(:) !! Intrinsic dimension for each candidate cluster.
        real(dp), intent(inout) :: centers(:, :) !! Current k-by-p locations, updated by concentration steps.
        real(dp), intent(inout) :: basis(:, :, :) !! Current affine direction arrays, updated by concentration steps.
        integer, intent(inout) :: labels(:) !! Current labels; zero marks trimmed observations.
        integer, intent(in) :: no_trim !! Number of observations retained after trimming.
        integer, intent(in) :: niter !! Maximum number of concentration steps.
        real(dp), intent(out) :: obj !! Sum of retained squared orthogonal distances.
        real(dp), allocatable :: dist(:, :)
        real(dp), allocatable :: mind(:)
        integer, allocatable :: order(:)
        integer, allocatable :: old(:)
        integer :: it
        integer :: i

        allocate(dist(size(x, 1), size(dimensions)), mind(size(x, 1)), order(size(x, 1)), old(size(labels)))
        do it = 1, max(1, niter)
            old = labels
            call affine_distances(x, dimensions, centers, basis, dist)
            do i = 1, size(x, 1)
                labels(i) = minloc(dist(i, :), dim=1)
                mind(i) = minval(dist(i, :))
            end do
            call sort_real_with_index(mind, order, ascending=.true.)
            if (no_trim < size(x, 1)) labels(order(no_trim + 1:)) = 0
            if (niter == 0 .or. all(labels == old)) exit
            call fit_affine_from_labels(x, dimensions, labels, centers, basis)
        end do
        call affine_distances(x, dimensions, centers, basis, dist)
        do i = 1, size(x, 1)
            if (labels(i) > 0) then
                mind(i) = dist(i, labels(i))
            else
                mind(i) = 0.0_dp
            end if
        end do
        obj = sum(mind)
    end subroutine rlg_steps

    subroutine affine_distances(x, dimensions, centers, basis, dist)
        real(dp), intent(in) :: x(:, :) !! Data matrix whose affine residual distances are requested.
        integer, intent(in) :: dimensions(:) !! Intrinsic dimension of each affine cluster.
        real(dp), intent(in) :: centers(:, :) !! k-by-p location matrix.
        real(dp), intent(in) :: basis(:, :, :) !! p-by-max(d,1)-by-k orthonormal direction arrays.
        real(dp), intent(out) :: dist(:, :) !! Output n-by-k squared orthogonal residual distances.
        real(dp), allocatable :: z(:)
        real(dp) :: proj
        integer :: i
        integer :: j
        integer :: q

        allocate(z(size(x, 2)))
        do j = 1, size(dimensions)
            do i = 1, size(x, 1)
                z = x(i, :) - centers(j, :)
                dist(i, j) = dot_product(z, z)
                do q = 1, dimensions(j)
                    proj = dot_product(z, basis(:, q, j))
                    dist(i, j) = dist(i, j) - proj * proj
                end do
                dist(i, j) = max(0.0_dp, dist(i, j))
            end do
        end do
    end subroutine affine_distances

    pure logical function all_clusters_present(labels, k) result(ok)
        integer, intent(in) :: labels(:) !! Cluster labels checked for complete positive-cluster representation.
        integer, intent(in) :: k !! Number of required positive cluster labels.
        integer :: j

        ok = .true.
        do j = 1, k
            if (count(labels == j) == 0) then
                ok = .false.
                return
            end if
        end do
    end function all_clusters_present

    subroutine fill_rlg_result(x, dimensions, centers, basis, labels, alpha, obj, result)
        real(dp), intent(in) :: x(:, :) !! Original data matrix stored in the returned result.
        integer, intent(in) :: dimensions(:) !! Intrinsic dimension for each fitted cluster.
        real(dp), intent(in) :: centers(:, :) !! Final k-by-p location matrix.
        real(dp), intent(in) :: basis(:, :, :) !! Final p-by-max(d,1)-by-k direction arrays.
        integer, intent(in) :: labels(:) !! Final assignments; zero denotes trimming.
        real(dp), intent(in) :: alpha !! Requested trimming fraction.
        real(dp), intent(in) :: obj !! Final retained squared-residual objective.
        type(rlg_result), intent(out) :: result !! Populated public robust linear-grouping result.

        result%n = size(x, 1)
        result%p = size(x, 2)
        result%k = size(dimensions)
        result%alpha = alpha
        result%obj = obj
        allocate(result%dimensions(result%k), result%cluster(result%n), result%centers(result%p, result%k), &
                 result%basis(result%p, size(basis, 2), result%k), result%x(result%n, result%p))
        result%dimensions = dimensions
        result%cluster = labels
        result%centers = transpose(centers)
        result%basis = basis
        result%x = x
    end subroutine fill_rlg_result

end module tclust_rlg
