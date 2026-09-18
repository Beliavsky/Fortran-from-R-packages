module tclust_tkmeans
    use tclust_kinds, only : dp
    use tclust_types, only : tkmeans_result
    use tclust_rng, only : set_tclust_seed, sample_index
    use tclust_linalg, only : sort_real_with_index
    implicit none
    private

    public :: tkmeans
    public :: tkmeans_refine

contains

    subroutine tkmeans(x, k, result, alpha, nstart, niter1, niter2, nkeep, points, zero_tol, seed)
        real(dp), intent(in) :: x(:, :) !! Numeric data matrix with observations in rows and variables in columns.
        integer, intent(in) :: k !! Number of spherical clusters initially searched for; must be positive.
        type(tkmeans_result), intent(out) :: result !! Best trimmed k-means fit, including assignments and centers.
        real(dp), intent(in), optional :: alpha !! Fraction of observations trimmed; default 0.05 and valid in [0,1).
        integer, intent(in), optional :: nstart !! Number of random initializations; default 500, ignored when points are supplied.
        integer, intent(in), optional :: niter1 !! Concentration steps in the initial phase; default 3.
        integer, intent(in), optional :: niter2 !! Maximum refinement steps for retained starts; default 20.
        integer, intent(in), optional :: nkeep !! Number of phase-one starts retained; default 5.
        real(dp), intent(in), optional :: points(:, :) !! Optional k-by-p initial center matrix; supplying it forces one start.
        real(dp), intent(in), optional :: zero_tol !! Nonnegative empty-cluster threshold. Default 1e-16.
        integer, intent(in), optional :: seed !! Optional deterministic seed for the intrinsic Fortran RNG.
        real(dp), allocatable :: obj(:)
        real(dp), allocatable :: centers(:, :, :)
        integer, allocatable :: labels(:, :)
        integer, allocatable :: order(:)
        real(dp), allocatable :: c(:, :)
        real(dp) :: a
        real(dp) :: ztol
        real(dp) :: best_obj
        integer :: n
        integer :: p
        integer :: ns
        integer :: ni1
        integer :: ni2
        integer :: nk
        integer :: no_trim
        integer :: j
        integer :: best_code
        integer, allocatable :: lab(:)

        n = size(x, 1)
        p = size(x, 2)
        if (n < 1 .or. p < 1) error stop 'tkmeans: x must be nonempty'
        if (k < 1) error stop 'tkmeans: k must be positive'
        a = 0.05_dp
        if (present(alpha)) a = alpha
        if (a < 0.0_dp .or. a >= 1.0_dp) error stop 'tkmeans: alpha must be in [0,1)'
        no_trim = floor(real(n, dp) * (1.0_dp - a))
        if (no_trim < k) error stop 'tkmeans: fewer retained observations than clusters'
        ns = 500
        if (present(nstart)) ns = nstart
        ni1 = 3
        if (present(niter1)) ni1 = niter1
        ni2 = 20
        if (present(niter2)) ni2 = niter2
        nk = 5
        if (present(nkeep)) nk = nkeep
        if (present(points)) then
            if (size(points, 1) /= k .or. size(points, 2) /= p) error stop 'tkmeans: points must have shape k by p'
            ns = 1
            nk = 1
        end if
        if (ns < 1 .or. nk < 1 .or. nk > ns .or. ni1 < 0 .or. ni2 < 0) error stop 'tkmeans: invalid iteration counts'
        ztol = 1.0e-16_dp
        if (present(zero_tol)) ztol = zero_tol
        if (ztol < 0.0_dp) error stop 'tkmeans: zero_tol must be nonnegative'
        if (present(seed)) call set_tclust_seed(seed)

        allocate(obj(ns), centers(k, p, ns), labels(n, ns), order(ns), c(k, p), lab(n))
        do j = 1, ns
            if (present(points)) then
                c = points
            else
                call initialize_centers(x, c)
            end if
            lab = 0
            call tkmeans_steps(x, c, lab, no_trim, ni1, ztol, obj(j), best_code)
            centers(:, :, j) = c
            labels(:, j) = lab
        end do
        call sort_real_with_index(obj, order, ascending=.true.)

        best_obj = huge(1.0_dp)
        best_code = 0
        do j = 1, nk
            lab = labels(:, order(j))
            call centers_from_labels(x, k, lab, c, ztol)
            call tkmeans_steps(x, c, lab, no_trim, ni2, ztol, obj(order(j)), best_code)
            if (obj(order(j)) < best_obj) then
                best_obj = obj(order(j))
                call fill_result(x, c, lab, a, best_obj, best_code, result)
            end if
        end do
    end subroutine tkmeans

    subroutine tkmeans_refine(x, cluster, result, alpha, niter, zero_tol)
        real(dp), intent(in) :: x(:, :) !! Numeric data matrix with observations in rows and variables in columns.
        integer, intent(in) :: cluster(:) !! Starting labels of length n; zero means trimmed.
        type(tkmeans_result), intent(out) :: result !! Refined trimmed k-means fit.
        real(dp), intent(in), optional :: alpha !! Trimming fraction; inferred from zero labels when absent.
        integer, intent(in), optional :: niter !! Maximum concentration steps; default 20.
        real(dp), intent(in), optional :: zero_tol !! Empty-cluster tolerance; default 1e-16.
        real(dp), allocatable :: centers(:, :)
        integer, allocatable :: lab(:)
        real(dp) :: a
        real(dp) :: obj
        real(dp) :: ztol
        integer :: k
        integer :: nit
        integer :: no_trim
        integer :: code

        if (size(cluster) /= size(x, 1)) error stop 'tkmeans_refine: cluster length mismatch'
        k = maxval(cluster)
        if (k < 1) error stop 'tkmeans_refine: at least one positive cluster is required'
        a = real(count(cluster == 0), dp) / real(size(cluster), dp)
        if (present(alpha)) a = alpha
        no_trim = floor(real(size(cluster), dp) * (1.0_dp - a))
        nit = 20
        if (present(niter)) nit = niter
        ztol = 1.0e-16_dp
        if (present(zero_tol)) ztol = zero_tol
        allocate(centers(k, size(x, 2)), lab(size(cluster)))
        lab = cluster
        call centers_from_labels(x, k, lab, centers, ztol)
        call tkmeans_steps(x, centers, lab, no_trim, nit, ztol, obj, code)
        call fill_result(x, centers, lab, a, obj, code, result)
    end subroutine tkmeans_refine

    subroutine initialize_centers(x, centers)
        real(dp), intent(in) :: x(:, :) !! Data matrix from which random observation centers are sampled.
        real(dp), intent(out) :: centers(:, :) !! Output k-by-p initial center matrix.
        integer :: j

        do j = 1, size(centers, 1)
            centers(j, :) = x(sample_index(size(x, 1)), :)
        end do
    end subroutine initialize_centers

    subroutine centers_from_labels(x, k, labels, centers, zero_tol)
        real(dp), intent(in) :: x(:, :) !! Data matrix used to recompute centers.
        integer, intent(in) :: k !! Number of candidate clusters.
        integer, intent(in) :: labels(:) !! Current labels; zero denotes trimmed observations.
        real(dp), intent(out) :: centers(:, :) !! Recomputed k-by-p center matrix.
        real(dp), intent(in) :: zero_tol !! Threshold for recognizing empty clusters.
        real(dp) :: nj
        integer :: j
        integer :: i

        centers = 0.0_dp
        do j = 1, k
            nj = real(count(labels == j), dp)
            if (nj > zero_tol) then
                do i = 1, size(x, 1)
                    if (labels(i) == j) centers(j, :) = centers(j, :) + x(i, :)
                end do
                centers(j, :) = centers(j, :) / nj
            end if
        end do
    end subroutine centers_from_labels

    subroutine tkmeans_steps(x, centers, labels, no_trim, niter, zero_tol, obj, code)
        real(dp), intent(in) :: x(:, :) !! Data matrix clustered by squared Euclidean distance.
        real(dp), intent(inout) :: centers(:, :) !! Current k-by-p centers, updated after each concentration step.
        integer, intent(inout) :: labels(:) !! Current assignments; zero denotes trimmed observations.
        integer, intent(in) :: no_trim !! Number of observations retained after trimming.
        integer, intent(in) :: niter !! Maximum number of concentration steps.
        real(dp), intent(in) :: zero_tol !! Threshold used when deciding whether a cluster is empty.
        real(dp), intent(out) :: obj !! Mean retained squared distance objective.
        integer, intent(out) :: code !! Convergence code; two indicates unchanged consecutive assignments.
        real(dp), allocatable :: dist(:, :)
        real(dp), allocatable :: mind(:)
        integer, allocatable :: idx(:)
        integer, allocatable :: old(:)
        integer :: i
        integer :: j
        integer :: it
        integer :: bestj

        allocate(dist(size(x, 1), size(centers, 1)), mind(size(x, 1)), idx(size(x, 1)), old(size(labels)))
        code = 0
        do it = 1, niter
            old = labels
            do j = 1, size(centers, 1)
                do i = 1, size(x, 1)
                    dist(i, j) = sum((x(i, :) - centers(j, :))**2)
                end do
            end do
            do i = 1, size(x, 1)
                bestj = minloc(dist(i, :), dim=1)
                labels(i) = bestj
                mind(i) = dist(i, bestj)
            end do
            call sort_real_with_index(mind, idx, ascending=.true.)
            if (no_trim < size(x, 1)) labels(idx(no_trim + 1:)) = 0
            if (all(labels == old)) then
                code = 2
                exit
            end if
            call centers_from_labels(x, size(centers, 1), labels, centers, zero_tol)
        end do
        if (niter == 0) then
            call tkmeans_assign(x, centers, labels, no_trim, mind)
        else
            call tkmeans_assign(x, centers, labels, no_trim, mind)
        end if
        obj = sum(mind, mask=labels > 0) / real(no_trim, dp)
    end subroutine tkmeans_steps

    subroutine tkmeans_assign(x, centers, labels, no_trim, mind)
        real(dp), intent(in) :: x(:, :) !! Data matrix assigned to the nearest centers.
        real(dp), intent(in) :: centers(:, :) !! Current k-by-p center matrix.
        integer, intent(out) :: labels(:) !! Nearest-center labels with trimmed observations set to zero.
        integer, intent(in) :: no_trim !! Number of observations retained after trimming.
        real(dp), intent(out) :: mind(:) !! Squared distance to the nearest center for each observation.
        real(dp), allocatable :: d(:)
        integer, allocatable :: idx(:)
        integer :: i
        integer :: j

        allocate(d(size(centers, 1)), idx(size(x, 1)))
        do i = 1, size(x, 1)
            do j = 1, size(centers, 1)
                d(j) = sum((x(i, :) - centers(j, :))**2)
            end do
            labels(i) = minloc(d, dim=1)
            mind(i) = minval(d)
        end do
        call sort_real_with_index(mind, idx, ascending=.true.)
        if (no_trim < size(x, 1)) labels(idx(no_trim + 1:)) = 0
    end subroutine tkmeans_assign

    subroutine fill_result(x, centers, labels, alpha, obj, code, result)
        real(dp), intent(in) :: x(:, :) !! Original data matrix stored in the result.
        real(dp), intent(in) :: centers(:, :) !! Final k-by-p center matrix.
        integer, intent(in) :: labels(:) !! Final labels with zero denoting trimming.
        real(dp), intent(in) :: alpha !! Requested trimming fraction.
        real(dp), intent(in) :: obj !! Final trimmed k-means objective.
        integer, intent(in) :: code !! Convergence code returned by the concentration steps.
        type(tkmeans_result), intent(out) :: result !! Populated public result object.
        integer :: j
        integer :: no_trim

        result%n = size(x, 1)
        result%p = size(x, 2)
        result%k = size(centers, 1)
        result%alpha = alpha
        result%obj = obj
        result%code = code
        allocate(result%cluster(result%n), result%size(result%k), result%weights(result%k), &
                 result%centers(result%p, result%k), result%x(result%n, result%p))
        result%cluster = labels
        result%x = x
        result%centers = transpose(centers)
        no_trim = count(labels > 0)
        do j = 1, result%k
            result%size(j) = real(count(labels == j), dp)
        end do
        if (no_trim > 0) then
            result%weights = result%size / real(no_trim, dp)
        else
            result%weights = 0.0_dp
        end if
    end subroutine fill_result

end module tclust_tkmeans
