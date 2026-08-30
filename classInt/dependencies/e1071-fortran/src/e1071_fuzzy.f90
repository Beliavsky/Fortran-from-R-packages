module e1071_fuzzy
    use e1071_kinds, only: dp
    implicit none
    private

    integer, parameter, public :: fuzzy_euclidean = 0
    integer, parameter, public :: fuzzy_manhattan = 1

    type, public :: fuzzy_cluster_model
        real(dp), allocatable :: centers(:, :)
        real(dp), allocatable :: radius(:)
        real(dp), allocatable :: membership(:, :)
        integer, allocatable :: cluster(:)
        integer, allocatable :: size(:)
        integer :: iterations = 0
        real(dp) :: within_error = 0.0_dp
        real(dp) :: m = 2.0_dp
        integer :: distance = fuzzy_euclidean
        logical :: shell = .false.
    end type fuzzy_cluster_model

    public :: cmeans_fit, cmeans_assign
    public :: cshell_fit, cshell_assign

contains

    subroutine cmeans_fit(x, initial_centers, model, weights, m, distance, iter_max, reltol, ufcl, rate_par)
        real(dp), intent(in) :: x(:, :) !! Observation-by-variable matrix used for fuzzy clustering.
        real(dp), intent(in) :: initial_centers(:, :) !! Initial cluster-center matrix; rows are clusters and columns match x.
        type(fuzzy_cluster_model), intent(out) :: model !! Fitted centers, memberships, hard assignments, sizes, and objective
        !! value.
        real(dp), intent(in), optional :: weights(:) !! Nonnegative observation weights; normalized internally and default to
        !! equal weights.
        real(dp), intent(in), optional :: m !! Fuzzification exponent greater than one; defaults to 2.
        integer, intent(in), optional :: distance !! Dissimilarity code: fuzzy_euclidean=0 or fuzzy_manhattan=1; default Euclidean.
        integer, intent(in), optional :: iter_max !! Maximum update sweeps; defaults to 100 and must be positive.
        real(dp), intent(in), optional :: reltol !! Relative objective convergence tolerance; defaults to sqrt(machine epsilon).
        logical, intent(in), optional :: ufcl !! If true, use on-line UFCL prototype updates rather than batch c-means; default
        !! false.
        real(dp), intent(in), optional :: rate_par !! UFCL learning-rate multiplier; defaults to 0.3 and is ignored for batch
        !! c-means.
        real(dp), allocatable :: w(:)
        real(dp), allocatable :: d(:, :)
        real(dp), allocatable :: u(:, :)
        real(dp), allocatable :: centers(:, :)
        real(dp) :: fuzz
        real(dp) :: exponent
        real(dp) :: tol
        real(dp) :: rate
        real(dp) :: old_value
        real(dp) :: new_value
        real(dp) :: lrate
        integer :: dist
        integer :: maxit
        integer :: iter
        integer :: i
        logical :: online

        if (size(initial_centers, 2) /= size(x, 2)) error stop "cmeans_fit: center dimension mismatch"
        if (size(initial_centers, 1) < 1 .or. size(initial_centers, 1) > size(x, 1)) then
            error stop "cmeans_fit: invalid number of centers"
        end if
        fuzz = 2.0_dp
        if (present(m)) fuzz = m
        if (fuzz <= 1.0_dp) error stop "cmeans_fit: m must exceed one"
        dist = fuzzy_euclidean
        if (present(distance)) dist = distance
        if (dist /= fuzzy_euclidean .and. dist /= fuzzy_manhattan) error stop "cmeans_fit: invalid distance code"
        maxit = 100
        if (present(iter_max)) maxit = iter_max
        if (maxit < 1) error stop "cmeans_fit: iter_max must be positive"
        tol = sqrt(epsilon(1.0_dp))
        if (present(reltol)) tol = reltol
        if (tol <= 0.0_dp) error stop "cmeans_fit: reltol must be positive"
        online = .false.
        if (present(ufcl)) online = ufcl
        rate = 0.3_dp
        if (present(rate_par)) rate = rate_par

        allocate(w(size(x, 1)))
        if (present(weights)) then
            if (size(weights) /= size(x, 1)) error stop "cmeans_fit: weights has wrong length"
            if (any(weights < 0.0_dp) .or. all(weights <= 0.0_dp)) error stop "cmeans_fit: invalid weights"
            w = weights / sum(weights)
        else
            w = 1.0_dp / real(size(x, 1), dp)
        end if
        centers = initial_centers
        allocate(d(size(x, 1), size(centers, 1)), u(size(x, 1), size(centers, 1)))
        exponent = 1.0_dp / (fuzz - 1.0_dp)
        call fuzzy_dissimilarities(x, centers, dist, d)
        call fuzzy_memberships(d, exponent, u)
        old_value = fuzzy_error(u, d, w, fuzz)
        new_value = old_value

        do iter = 1, maxit
            if (online) then
                lrate = rate * (1.0_dp - real(iter, dp) / real(maxit, dp))
                do i = 1, size(x, 1)
                    call fuzzy_dissimilarities_row(x, centers, dist, i, d)
                    call fuzzy_memberships_row(d, exponent, i, u)
                    call ufcl_update_prototypes(x, u, w, fuzz, dist, lrate, i, centers)
                end do
            else
                call update_prototypes(x, u, w, fuzz, dist, centers)
                call fuzzy_dissimilarities(x, centers, dist, d)
                call fuzzy_memberships(d, exponent, u)
            end if
            if (online) call fuzzy_dissimilarities(x, centers, dist, d)
            new_value = fuzzy_error(u, d, w, fuzz)
            if (abs(old_value - new_value) < tol * (old_value + tol)) exit
            old_value = new_value
        end do

        model%centers = centers
        model%membership = u
        model%m = fuzz
        model%distance = dist
        model%shell = .false.
        model%iterations = min(iter, maxit)
        model%within_error = new_value
        call derive_hard_clusters(u, model%cluster, model%size)
    end subroutine cmeans_fit

    subroutine cmeans_assign(model, x, membership, cluster)
        type(fuzzy_cluster_model), intent(in) :: model !! Previously fitted non-shell fuzzy-cluster model supplying centers and
        !! fuzzifier.
        real(dp), intent(in) :: x(:, :) !! New observation-by-variable matrix to assign to the fitted centers.
        real(dp), allocatable, intent(out) :: membership(:, :) !! Fuzzy memberships; rows correspond to x and columns to clusters.
        integer, allocatable, intent(out), optional :: cluster(:) !! Optional one-based hard cluster obtained by maximum membership.
        real(dp), allocatable :: d(:, :)

        if (model%shell) error stop "cmeans_assign: shell model supplied"
        if (size(x, 2) /= size(model%centers, 2)) error stop "cmeans_assign: variable count mismatch"
        allocate(d(size(x, 1), size(model%centers, 1)), membership(size(x, 1), size(model%centers, 1)))
        call fuzzy_dissimilarities(x, model%centers, model%distance, d)
        call fuzzy_memberships(d, 1.0_dp / (model%m - 1.0_dp), membership)
        if (present(cluster)) call hard_labels(membership, cluster)
    end subroutine cmeans_assign

    subroutine cshell_fit(x, initial_centers, model, radius, m, distance, iter_max)
        real(dp), intent(in) :: x(:, :) !! Observation-by-variable matrix used for fuzzy c-shell clustering.
        real(dp), intent(in) :: initial_centers(:, :) !! Initial shell centers; rows are shells and columns match x.
        type(fuzzy_cluster_model), intent(out) :: model !! Fitted shell centers, radii, memberships, hard assignments, and
        !! objective value.
        real(dp), intent(in), optional :: radius(:) !! Initial shell radii; defaults to 0.2 for every shell.
        real(dp), intent(in), optional :: m !! Fuzzification exponent greater than one; defaults to 2.
        integer, intent(in), optional :: distance !! Distance code: fuzzy_euclidean=0 or fuzzy_manhattan=1; default Euclidean.
        integer, intent(in), optional :: iter_max !! Maximum shell-update sweeps including optimization refinements; defaults to
        !! 100.
        real(dp), allocatable :: centers(:, :)
        real(dp), allocatable :: radii(:)
        real(dp), allocatable :: u(:, :)
        real(dp), allocatable :: old_u(:, :)
        real(dp) :: fuzz
        real(dp) :: error_value
        real(dp) :: conv
        integer :: dist
        integer :: maxit
        integer :: iter
        integer :: flag
        integer :: c

        if (size(initial_centers, 2) /= size(x, 2)) error stop "cshell_fit: center dimension mismatch"
        if (size(initial_centers, 1) < 1) error stop "cshell_fit: at least one center is required"
        fuzz = 2.0_dp
        if (present(m)) fuzz = m
        if (fuzz <= 1.0_dp) error stop "cshell_fit: m must exceed one"
        dist = fuzzy_euclidean
        if (present(distance)) dist = distance
        if (dist /= fuzzy_euclidean .and. dist /= fuzzy_manhattan) error stop "cshell_fit: invalid distance"
        maxit = 100
        if (present(iter_max)) maxit = iter_max
        if (maxit < 1) error stop "cshell_fit: iter_max must be positive"

        centers = initial_centers
        allocate(radii(size(centers, 1)), u(size(x, 1), size(centers, 1)), old_u(size(x, 1), size(centers, 1)))
        if (present(radius)) then
            if (size(radius) /= size(radii)) error stop "cshell_fit: radius has wrong length"
            radii = radius
        else
            radii = 0.2_dp
        end if
        call shell_memberships(x, centers, radii, fuzz, dist, u)
        old_u = u
        flag = 0
        error_value = 0.0_dp

        do iter = 1, maxit
            if (flag == 0 .or. flag == 5) call shell_update_geometry(x, u, fuzz, dist, centers, radii)
            old_u = u
            call shell_memberships(x, centers, radii, fuzz, dist, u)
            call shell_error(x, centers, radii, u, fuzz, dist, error_value)
            conv = sum(abs(u - old_u))
            if (conv <= real(size(x, 1) * size(x, 2), dp) * 0.002_dp) then
                flag = 2
                exit
            else if (conv <= real(size(x, 1) * size(x, 2), dp) * 0.2_dp) then
                flag = 1
                do c = 1, size(centers, 1)
                    call refine_shell_center(x, u(:, c), fuzz, centers(c, :), radii(c))
                end do
                flag = 5
            else
                flag = 5
            end if
        end do

        model%centers = centers
        model%radius = radii
        model%membership = u
        model%m = fuzz
        model%distance = dist
        model%shell = .true.
        model%iterations = min(iter, maxit)
        model%within_error = error_value
        call derive_hard_clusters(u, model%cluster, model%size)
    end subroutine cshell_fit

    subroutine cshell_assign(model, x, membership, cluster)
        type(fuzzy_cluster_model), intent(in) :: model !! Previously fitted shell model supplying centers, radii, distance, and
        !! fuzzifier.
        real(dp), intent(in) :: x(:, :) !! New observation-by-variable matrix to assign to shells.
        real(dp), allocatable, intent(out) :: membership(:, :) !! Fuzzy shell memberships for each new observation.
        integer, allocatable, intent(out), optional :: cluster(:) !! Optional one-based hard shell assignment from maximum
        !! membership.

        if (.not. model%shell) error stop "cshell_assign: non-shell model supplied"
        if (size(x, 2) /= size(model%centers, 2)) error stop "cshell_assign: variable count mismatch"
        allocate(membership(size(x, 1), size(model%centers, 1)))
        call shell_memberships(x, model%centers, model%radius, model%m, model%distance, membership)
        if (present(cluster)) call hard_labels(membership, cluster)
    end subroutine cshell_assign

    subroutine fuzzy_dissimilarities(x, centers, distance, d)
        real(dp), intent(in) :: x(:, :) !! Observation-by-variable data matrix.
        real(dp), intent(in) :: centers(:, :) !! Cluster-center matrix.
        integer, intent(in) :: distance !! Distance code selecting squared Euclidean or Manhattan dissimilarity.
        real(dp), intent(out) :: d(:, :) !! Observation-by-center dissimilarity matrix.
        integer :: i

        do i = 1, size(x, 1)
            call fuzzy_dissimilarities_row(x, centers, distance, i, d)
        end do
    end subroutine fuzzy_dissimilarities

    subroutine fuzzy_dissimilarities_row(x, centers, distance, row, d)
        real(dp), intent(in) :: x(:, :) !! Observation-by-variable data matrix.
        real(dp), intent(in) :: centers(:, :) !! Cluster-center matrix.
        integer, intent(in) :: distance !! Distance code selecting squared Euclidean or Manhattan dissimilarity.
        integer, intent(in) :: row !! One-based observation row whose dissimilarities are updated.
        real(dp), intent(inout) :: d(:, :) !! Dissimilarity matrix receiving the selected row.
        integer :: c

        do c = 1, size(centers, 1)
            if (distance == fuzzy_euclidean) then
                d(row, c) = sum((x(row, :) - centers(c, :))**2)
            else
                d(row, c) = sum(abs(x(row, :) - centers(c, :)))
            end if
        end do
    end subroutine fuzzy_dissimilarities_row

    subroutine fuzzy_memberships(d, exponent, u)
        real(dp), intent(in) :: d(:, :) !! Observation-by-center dissimilarity matrix.
        real(dp), intent(in) :: exponent !! Membership exponent 1/(m-1), positive when m exceeds one.
        real(dp), intent(out) :: u(:, :) !! Observation-by-center membership matrix whose rows sum to one.
        integer :: i

        do i = 1, size(d, 1)
            call fuzzy_memberships_row(d, exponent, i, u)
        end do
    end subroutine fuzzy_memberships

    subroutine fuzzy_memberships_row(d, exponent, row, u)
        real(dp), intent(in) :: d(:, :) !! Observation-by-center dissimilarity matrix.
        real(dp), intent(in) :: exponent !! Membership exponent 1/(m-1).
        integer, intent(in) :: row !! One-based observation row whose memberships are updated.
        real(dp), intent(inout) :: u(:, :) !! Membership matrix receiving the selected row.
        integer :: c
        integer :: nzero
        real(dp) :: total

        nzero = count(d(row, :) <= 0.0_dp)
        if (nzero > 0) then
            u(row, :) = 0.0_dp
            where (d(row, :) <= 0.0_dp) u(row, :) = 1.0_dp / real(nzero, dp)
        else
            do c = 1, size(d, 2)
                u(row, c) = 1.0_dp / d(row, c)**exponent
            end do
            total = sum(u(row, :))
            u(row, :) = u(row, :) / total
        end if
    end subroutine fuzzy_memberships_row

    subroutine update_prototypes(x, u, weights, m, distance, centers)
        real(dp), intent(in) :: x(:, :) !! Observation-by-variable data matrix.
        real(dp), intent(in) :: u(:, :) !! Observation-by-center membership matrix.
        real(dp), intent(in) :: weights(:) !! Normalized nonnegative observation weights.
        real(dp), intent(in) :: m !! Fuzzification exponent applied to memberships.
        integer, intent(in) :: distance !! Distance code selecting mean or weighted-median prototype updates.
        real(dp), intent(inout) :: centers(:, :) !! Cluster centers replaced by their updated prototypes.
        real(dp), allocatable :: wlocal(:)
        real(dp), allocatable :: xlocal(:)
        real(dp) :: total
        real(dp) :: weight
        integer :: c
        integer :: j
        integer :: i

        allocate(wlocal(size(x, 1)), xlocal(size(x, 1)))
        do c = 1, size(centers, 1)
            if (distance == fuzzy_euclidean) then
                centers(c, :) = 0.0_dp
                total = 0.0_dp
                do i = 1, size(x, 1)
                    weight = weights(i) * u(i, c)**m
                    total = total + weight
                    centers(c, :) = centers(c, :) + weight * x(i, :)
                end do
                if (total > 0.0_dp) centers(c, :) = centers(c, :) / total
            else
                do j = 1, size(x, 2)
                    xlocal = x(:, j)
                    wlocal = weights * u(:, c)**m
                    centers(c, j) = weighted_median(xlocal, wlocal)
                end do
            end if
        end do
    end subroutine update_prototypes

    subroutine ufcl_update_prototypes(x, u, weights, m, distance, learning_rate, row, centers)
        real(dp), intent(in) :: x(:, :) !! Observation-by-variable data matrix.
        real(dp), intent(in) :: u(:, :) !! Current observation-by-center memberships.
        real(dp), intent(in) :: weights(:) !! Normalized observation weights.
        real(dp), intent(in) :: m !! Fuzzification exponent applied to memberships.
        integer, intent(in) :: distance !! Distance code selecting raw or sign-gradient prototype steps.
        real(dp), intent(in) :: learning_rate !! Current UFCL learning rate after iteration decay.
        integer, intent(in) :: row !! One-based observation whose stochastic prototype update is applied.
        real(dp), intent(inout) :: centers(:, :) !! Cluster centers updated in place.
        real(dp) :: grad
        integer :: c
        integer :: j

        do c = 1, size(centers, 1)
            do j = 1, size(centers, 2)
                grad = x(row, j) - centers(c, j)
                if (distance == fuzzy_manhattan) grad = sign(1.0_dp, grad)
                if (abs(x(row, j) - centers(c, j)) <= 0.0_dp) grad = 0.0_dp
                centers(c, j) = centers(c, j) + learning_rate * weights(row) * u(row, c)**m * grad
            end do
        end do
    end subroutine ufcl_update_prototypes

    function fuzzy_error(u, d, weights, m) result(value)
        real(dp), intent(in) :: u(:, :) !! Observation-by-center membership matrix.
        real(dp), intent(in) :: d(:, :) !! Observation-by-center dissimilarity matrix.
        real(dp), intent(in) :: weights(:) !! Normalized nonnegative observation weights.
        real(dp), intent(in) :: m !! Fuzzification exponent applied to memberships.
        real(dp) :: value
        integer :: i

        value = 0.0_dp
        do i = 1, size(u, 1)
            value = value + weights(i) * sum(u(i, :)**m * d(i, :))
        end do
    end function fuzzy_error

    function weighted_median(x, weights) result(value)
        real(dp), intent(in) :: x(:) !! Values whose weighted median is requested; sorted internally.
        real(dp), intent(in) :: weights(:) !! Nonnegative weights paired with x; at least one weight must be positive.
        real(dp) :: value
        real(dp), allocatable :: sx(:)
        real(dp), allocatable :: sw(:)
        integer, allocatable :: order(:)
        real(dp) :: total
        real(dp) :: cumulative
        real(dp) :: cumulative_wx
        real(dp) :: score
        real(dp) :: best
        real(dp) :: key
        integer :: i
        integer :: j
        integer :: ikey

        if (size(x) /= size(weights)) error stop "weighted_median: size mismatch"
        if (any(weights < 0.0_dp) .or. sum(weights) <= 0.0_dp) error stop "weighted_median: invalid weights"
        allocate(sx(size(x)), sw(size(x)), order(size(x)))
        sx = x
        order = [(i, i = 1, size(x))]
        do i = 2, size(x)
            key = sx(i)
            ikey = order(i)
            j = i - 1
            do while (j >= 1)
                if (sx(j) <= key) exit
                sx(j + 1) = sx(j)
                order(j + 1) = order(j)
                j = j - 1
            end do
            sx(j + 1) = key
            order(j + 1) = ikey
        end do
        total = sum(weights)
        do i = 1, size(x)
            sw(i) = weights(order(i)) / total
        end do
        cumulative = 0.0_dp
        cumulative_wx = 0.0_dp
        best = huge(1.0_dp)
        value = sx(1)
        do i = 1, size(x)
            cumulative = cumulative + sw(i)
            cumulative_wx = cumulative_wx + sw(i) * sx(i)
            score = sx(i) * (cumulative - 0.5_dp) - cumulative_wx
            if (score < best) then
                best = score
                value = sx(i)
            end if
        end do
    end function weighted_median

    subroutine shell_update_geometry(x, u, m, distance, centers, radius)
        real(dp), intent(in) :: x(:, :) !! Observation-by-variable data matrix.
        real(dp), intent(in) :: u(:, :) !! Observation-by-shell fuzzy memberships.
        real(dp), intent(in) :: m !! Fuzzification exponent applied to membership weights.
        integer, intent(in) :: distance !! Distance code controlling Euclidean-radius or Manhattan-radius updates.
        real(dp), intent(inout) :: centers(:, :) !! Shell centers replaced by membership-weighted means.
        real(dp), intent(inout) :: radius(:) !! Shell radii replaced by membership-weighted mean radial distances.
        real(dp) :: weight
        real(dp) :: total
        real(dp) :: dist_value
        integer :: c
        integer :: i

        do c = 1, size(centers, 1)
            centers(c, :) = 0.0_dp
            total = 0.0_dp
            do i = 1, size(x, 1)
                weight = u(i, c)**m
                total = total + weight
                centers(c, :) = centers(c, :) + weight * x(i, :)
            end do
            if (total > 0.0_dp) centers(c, :) = centers(c, :) / total
            radius(c) = 0.0_dp
            do i = 1, size(x, 1)
                weight = u(i, c)**m
                if (distance == fuzzy_euclidean) then
                    dist_value = sqrt(sum((x(i, :) - centers(c, :))**2))
                else
                    dist_value = sum(abs(x(i, :) - centers(c, :)))
                end if
                radius(c) = radius(c) + weight * dist_value
            end do
            if (total > 0.0_dp) radius(c) = radius(c) / total
        end do
    end subroutine shell_update_geometry

    subroutine shell_memberships(x, centers, radius, m, distance, u)
        real(dp), intent(in) :: x(:, :) !! Observation-by-variable data matrix.
        real(dp), intent(in) :: centers(:, :) !! Current shell centers.
        real(dp), intent(in) :: radius(:) !! Current shell radii, one per center.
        real(dp), intent(in) :: m !! Fuzzification exponent greater than one.
        integer, intent(in) :: distance !! Distance code selecting Euclidean or Manhattan radial distance.
        real(dp), intent(out) :: u(:, :) !! Observation-by-shell memberships whose rows are normalized by radial ratios.
        real(dp), allocatable :: shell_distance(:)
        real(dp) :: exponent
        real(dp) :: total
        integer :: i
        integer :: c
        integer :: nzero

        exponent = 2.0_dp / (m - 1.0_dp)
        allocate(shell_distance(size(centers, 1)))
        do i = 1, size(x, 1)
            do c = 1, size(centers, 1)
                if (distance == fuzzy_euclidean) then
                    shell_distance(c) = abs(sqrt(sum((x(i, :) - centers(c, :))**2)) - radius(c))
                else
                    shell_distance(c) = abs(sum(abs(x(i, :) - centers(c, :))) - radius(c))
                end if
            end do
            nzero = count(shell_distance <= sqrt(tiny(1.0_dp)))
            if (nzero > 0) then
                u(i, :) = 0.0_dp
                where (shell_distance <= sqrt(tiny(1.0_dp))) u(i, :) = 1.0_dp / real(nzero, dp)
            else
                do c = 1, size(centers, 1)
                    u(i, c) = 1.0_dp / shell_distance(c)**exponent
                end do
                total = sum(u(i, :))
                u(i, :) = u(i, :) / total
            end if
        end do
    end subroutine shell_memberships

    subroutine shell_error(x, centers, radius, u, m, distance, value)
        real(dp), intent(in) :: x(:, :) !! Observation-by-variable data matrix.
        real(dp), intent(in) :: centers(:, :) !! Current shell centers.
        real(dp), intent(in) :: radius(:) !! Current shell radii.
        real(dp), intent(in) :: u(:, :) !! Current observation-by-shell memberships.
        real(dp), intent(in) :: m !! Fuzzification exponent used in the shell objective.
        integer, intent(in) :: distance !! Distance code selecting Euclidean or Manhattan radial errors.
        real(dp), intent(out) :: value !! Sum of membership-weighted squared radial deviations.
        real(dp) :: radial
        integer :: i
        integer :: c

        value = 0.0_dp
        do c = 1, size(centers, 1)
            do i = 1, size(x, 1)
                if (distance == fuzzy_euclidean) then
                    radial = abs(sqrt(sum((x(i, :) - centers(c, :))**2)) - radius(c))
                else
                    radial = abs(sum(abs(x(i, :) - centers(c, :))) - radius(c))
                end if
                value = value + u(i, c)**m * radial**2
            end do
        end do
    end subroutine shell_error

    subroutine refine_shell_center(x, membership, m, center, radius)
        real(dp), intent(in) :: x(:, :) !! Observation-by-variable data used in the R-side c-shell optimization objective.
        real(dp), intent(in) :: membership(:) !! Memberships for one shell; length equals the number of observations.
        real(dp), intent(in) :: m !! Fuzzification exponent applied to memberships in the objective.
        real(dp), intent(inout) :: center(:) !! Shell center refined by a finite-difference nonlinear conjugate-gradient search.
        real(dp), intent(inout) :: radius !! Shell radius refined together with center; no positivity constraint is imposed
        !! upstream.
        real(dp), allocatable :: param(:)
        real(dp), allocatable :: grad(:)
        real(dp), allocatable :: old_grad(:)
        real(dp), allocatable :: direction(:)
        real(dp), allocatable :: trial(:)
        real(dp) :: f0
        real(dp) :: ftrial
        real(dp) :: step
        real(dp) :: beta
        integer :: iter

        allocate(param(size(center) + 1), grad(size(center) + 1), old_grad(size(center) + 1), &
                 direction(size(center) + 1), trial(size(center) + 1))
        param(:size(center)) = center
        param(size(param)) = radius
        call shell_system_gradient(x, membership, m, param, grad)
        direction = -grad
        do iter = 1, 40
            if (sqrt(sum(grad * grad)) < 1.0e-6_dp) exit
            f0 = shell_system_objective(x, membership, m, param)
            step = 1.0_dp
            do
                trial = param + step * direction
                ftrial = shell_system_objective(x, membership, m, trial)
                if (ftrial < f0 .or. step < 1.0e-8_dp) exit
                step = 0.5_dp * step
            end do
            if (step < 1.0e-8_dp) exit
            param = trial
            old_grad = grad
            call shell_system_gradient(x, membership, m, param, grad)
            beta = max(0.0_dp, dot_product(grad, grad - old_grad) / max(dot_product(old_grad, old_grad), tiny(1.0_dp)))
            direction = -grad + beta * direction
            if (dot_product(direction, grad) >= 0.0_dp) direction = -grad
        end do
        center = param(:size(center))
        radius = param(size(param))
    end subroutine refine_shell_center

    function shell_system_objective(x, membership, m, param) result(value)
        real(dp), intent(in) :: x(:, :) !! Observation-by-variable data matrix in the Euclidean R-side c-shell system objective.
        real(dp), intent(in) :: membership(:) !! Membership values for the shell being optimized.
        real(dp), intent(in) :: m !! Fuzzification exponent applied to memberships.
        real(dp), intent(in) :: param(:) !! Candidate center coordinates followed by the candidate radius.
        real(dp) :: value
        real(dp), allocatable :: equation_sum(:)
        real(dp) :: radial
        real(dp) :: factor
        real(dp) :: weight
        integer :: i
        integer :: j

        allocate(equation_sum(size(x, 2)))
        equation_sum = 0.0_dp
        value = 0.0_dp
        do i = 1, size(x, 1)
            radial = sqrt(sum((x(i, :) - param(:size(x, 2)))**2))
            weight = membership(i)**m
            value = value + (weight * (radial - param(size(param))))**2
            if (radial > sqrt(tiny(1.0_dp))) then
                factor = weight * (1.0_dp - param(size(param)) / radial)
                do j = 1, size(x, 2)
                    equation_sum(j) = equation_sum(j) + factor * (x(i, j) - param(j))
                end do
            end if
        end do
        value = value + sum(equation_sum**2)
    end function shell_system_objective

    subroutine shell_system_gradient(x, membership, m, param, gradient)
        real(dp), intent(in) :: x(:, :) !! Observation-by-variable data matrix used by the c-shell refinement objective.
        real(dp), intent(in) :: membership(:) !! Membership vector for the shell being refined.
        real(dp), intent(in) :: m !! Fuzzification exponent applied in the refinement objective.
        real(dp), intent(in) :: param(:) !! Candidate center coordinates followed by radius.
        real(dp), intent(out) :: gradient(:) !! Central finite-difference gradient with the same length as param.
        real(dp), allocatable :: plus(:)
        real(dp), allocatable :: minus(:)
        real(dp) :: h
        integer :: j

        allocate(plus(size(param)), minus(size(param)))
        do j = 1, size(param)
            h = 1.0e-5_dp * max(1.0_dp, abs(param(j)))
            plus = param
            minus = param
            plus(j) = plus(j) + h
            minus(j) = minus(j) - h
            gradient(j) = (shell_system_objective(x, membership, m, plus) &
                           - shell_system_objective(x, membership, m, minus)) / (2.0_dp * h)
        end do
    end subroutine shell_system_gradient

    subroutine derive_hard_clusters(u, labels, sizes)
        real(dp), intent(in) :: u(:, :) !! Observation-by-cluster membership matrix.
        integer, allocatable, intent(out) :: labels(:) !! One-based maximum-membership cluster label for every observation.
        integer, allocatable, intent(out) :: sizes(:) !! Number of hard-assigned observations in every cluster.
        integer :: i
        integer :: c

        allocate(labels(size(u, 1)), sizes(size(u, 2)))
        sizes = 0
        do i = 1, size(u, 1)
            labels(i) = maxloc(u(i, :), dim=1)
            c = labels(i)
            sizes(c) = sizes(c) + 1
        end do
    end subroutine derive_hard_clusters

    subroutine hard_labels(u, labels)
        real(dp), intent(in) :: u(:, :) !! Observation-by-cluster membership matrix.
        integer, allocatable, intent(out) :: labels(:) !! One-based maximum-membership cluster label for each observation.
        integer :: i

        allocate(labels(size(u, 1)))
        do i = 1, size(u, 1)
            labels(i) = maxloc(u(i, :), dim=1)
        end do
    end subroutine hard_labels

end module e1071_fuzzy
