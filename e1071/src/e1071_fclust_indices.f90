module e1071_fclust_indices
    use e1071_kinds, only: dp
    use e1071_fuzzy, only: fuzzy_cluster_model
    implicit none
    private

    type, public :: fclust_indices_result
        real(dp) :: fuzzy_hypervolume = 0.0_dp
        real(dp) :: average_partition_density = 0.0_dp
        real(dp) :: partition_density = 0.0_dp
        real(dp) :: xie_beni = 0.0_dp
        real(dp) :: fukuyama_sugeno = 0.0_dp
        real(dp) :: partition_coefficient = 0.0_dp
        real(dp) :: partition_entropy = 0.0_dp
        real(dp) :: proportion_exponent = 0.0_dp
        real(dp) :: separation_index = 0.0_dp
    end type fclust_indices_result

    public :: fclust_indices

contains

    function fclust_indices(model, x) result(result)
        type(fuzzy_cluster_model), intent(in) :: model !! Fitted fuzzy clustering model containing centers, memberships,
        !! assignments, and objective.
        real(dp), intent(in) :: x(:, :) !! Original observation-by-variable matrix corresponding rowwise to model memberships.
        type(fclust_indices_result) :: result

        if (size(x, 1) /= size(model%membership, 1) .or. size(x, 2) /= size(model%centers, 2)) then
            error stop "fclust_indices: data/model shape mismatch"
        end if
        call gath_geva(model, x, result%fuzzy_hypervolume, result%average_partition_density, result%partition_density)
        result%xie_beni = xie_beni(model)
        result%fukuyama_sugeno = fukuyama_sugeno(model)
        result%partition_coefficient = sum(model%membership**2) / real(size(model%membership, 1), dp)
        result%partition_entropy = partition_entropy(model%membership)
        result%proportion_exponent = proportion_exponent(model%membership, size(model%centers, 2))
        result%separation_index = separation_index(model, x)
    end function fclust_indices

    subroutine gath_geva(model, x, fhv, apd, pd)
        type(fuzzy_cluster_model), intent(in) :: model !! Fuzzy clustering model whose membership exponent is treated as m=2
        !! upstream.
        real(dp), intent(in) :: x(:, :) !! Original data used to form per-cluster scatter matrices.
        real(dp), intent(out) :: fhv !! Sum of square roots of cluster scatter determinants.
        real(dp), intent(out) :: apd !! Average partition density from points inside each scatter ellipsoid.
        real(dp), intent(out) :: pd !! Overall partition density normalized by fuzzy hypervolume.
        real(dp), allocatable :: scatter(:, :)
        real(dp), allocatable :: inverse(:, :)
        real(dp), allocatable :: diff(:)
        real(dp) :: denom
        real(dp) :: inside_weight
        real(dp) :: det
        real(dp) :: control
        integer :: c
        integer :: i
        integer :: info

        fhv = 0.0_dp
        apd = 0.0_dp
        pd = 0.0_dp
        allocate(scatter(size(x, 2), size(x, 2)), inverse(size(x, 2), size(x, 2)), diff(size(x, 2)))
        do c = 1, size(model%centers, 1)
            scatter = 0.0_dp
            denom = sum(model%membership(:, c))
            if (denom <= 0.0_dp) cycle
            do i = 1, size(x, 1)
                diff = x(i, :) - model%centers(c, :)
                scatter = scatter + model%membership(i, c) * outer_product(diff, diff)
            end do
            scatter = scatter / denom
            call invert_and_det(scatter, inverse, det, info)
            if (info /= 0 .or. det <= 0.0_dp) cycle
            inside_weight = 0.0_dp
            do i = 1, size(x, 1)
                diff = x(i, :) - model%centers(c, :)
                control = dot_product(diff, matmul(inverse, diff))
                if (control < 1.0_dp) inside_weight = inside_weight + model%membership(i, c)
            end do
            fhv = fhv + sqrt(det)
            apd = apd + inside_weight / sqrt(det)
            pd = pd + inside_weight
        end do
        if (size(model%centers, 1) > 0) apd = apd / real(size(model%centers, 1), dp)
        if (fhv > 0.0_dp) pd = pd / fhv
    end subroutine gath_geva

    function xie_beni(model) result(value)
        type(fuzzy_cluster_model), intent(in) :: model !! Fitted fuzzy clustering model supplying within-error and center locations.
        real(dp) :: value
        real(dp) :: minimum
        real(dp) :: d2
        integer :: i
        integer :: j

        minimum = huge(1.0_dp)
        do i = 1, size(model%centers, 1) - 1
            do j = i + 1, size(model%centers, 1)
                d2 = sum((model%centers(i, :) - model%centers(j, :))**2)
                minimum = min(minimum, d2)
            end do
        end do
        if (minimum <= 0.0_dp .or. minimum >= huge(1.0_dp)) then
            value = huge(1.0_dp)
        else
            value = model%within_error / (real(size(model%membership, 1), dp) * minimum)
        end if
    end function xie_beni

    function fukuyama_sugeno(model) result(value)
        type(fuzzy_cluster_model), intent(in) :: model !! Fitted fuzzy model supplying centers, memberships, and within-cluster
        !! error.
        real(dp) :: value
        real(dp), allocatable :: mean_center(:)
        real(dp) :: membership_weight
        integer :: c

        allocate(mean_center(size(model%centers, 2)))
        mean_center = sum(model%centers, dim=1) / real(size(model%centers, 1), dp)
        value = model%within_error
        do c = 1, size(model%centers, 1)
            membership_weight = sum(model%membership(:, c)**2)
            value = value - membership_weight * sum((model%centers(c, :) - mean_center)**2)
        end do
    end function fukuyama_sugeno

    function partition_entropy(membership) result(value)
        real(dp), intent(in) :: membership(:, :) !! Observation-by-cluster fuzzy membership matrix.
        real(dp) :: value
        integer :: i
        integer :: c

        value = 0.0_dp
        do i = 1, size(membership, 1)
            do c = 1, size(membership, 2)
                if (membership(i, c) > 0.0_dp) value = value + membership(i, c) * log(membership(i, c))
            end do
        end do
        value = -value / real(size(membership, 1), dp)
    end function partition_entropy

    function proportion_exponent(membership, k) result(value)
        real(dp), intent(in) :: membership(:, :) !! Observation-by-cluster membership matrix.
        integer, intent(in) :: k !! e1071 exponent parameter, equal upstream to the number of center columns.
        real(dp) :: value
        real(dp) :: product_value
        real(dp) :: aexp
        real(dp) :: umax
        integer :: i
        integer :: l
        integer :: upper

        product_value = 1.0_dp
        do i = 1, size(membership, 1)
            umax = maxval(membership(i, :))
            if (umax <= 0.0_dp) cycle
            upper = int(1.0_dp / umax)
            aexp = 0.0_dp
            do l = 1, upper
                if (l > k) cycle
                aexp = aexp + (-1.0_dp)**(l + 1) * binomial_real(k, l) * max(0.0_dp, 1.0_dp - real(l, dp) * umax)**(k - 1)
            end do
            product_value = product_value * aexp
        end do
        if (product_value > 0.0_dp) then
            value = -log(product_value)
        else
            value = huge(1.0_dp)
        end if
    end function proportion_exponent

    function separation_index(model, x) result(value)
        type(fuzzy_cluster_model), intent(in) :: model !! Fitted fuzzy model supplying hard cluster labels.
        real(dp), intent(in) :: x(:, :) !! Original observation-by-variable matrix used for within/between Euclidean distances.
        real(dp) :: value
        real(dp), allocatable :: diameter(:)
        real(dp) :: max_diameter
        real(dp) :: d
        integer :: c
        integer :: c2
        integer :: i
        integer :: j

        allocate(diameter(size(model%centers, 1)))
        diameter = 0.0_dp
        do c = 1, size(diameter)
            do i = 1, size(x, 1) - 1
                if (model%cluster(i) /= c) cycle
                do j = i + 1, size(x, 1)
                    if (model%cluster(j) /= c) cycle
                    d = sqrt(sum((x(i, :) - x(j, :))**2))
                    diameter(c) = max(diameter(c), d)
                end do
            end do
        end do
        max_diameter = maxval(diameter)
        if (max_diameter <= 0.0_dp) then
            value = huge(1.0_dp)
            return
        end if
        value = huge(1.0_dp)
        do c = 1, size(diameter) - 1
            do c2 = c + 1, size(diameter)
                do i = 1, size(x, 1)
                    if (model%cluster(i) /= c) cycle
                    do j = 1, size(x, 1)
                        if (model%cluster(j) /= c2) cycle
                        d = sqrt(sum((x(i, :) - x(j, :))**2)) / max_diameter
                        value = min(value, d)
                    end do
                end do
            end do
        end do
    end function separation_index

    pure function outer_product(x, y) result(matrix)
        real(dp), intent(in) :: x(:) !! Column vector forming the left factor of an outer product.
        real(dp), intent(in) :: y(:) !! Row vector forming the right factor of an outer product.
        real(dp) :: matrix(size(x), size(y))
        integer :: i

        do i = 1, size(x)
            matrix(i, :) = x(i) * y
        end do
    end function outer_product

    subroutine invert_and_det(a, inverse, determinant, info)
        real(dp), intent(in) :: a(:, :) !! Square matrix to invert and take the determinant of.
        real(dp), intent(out) :: inverse(:, :) !! Matrix inverse when info is zero.
        real(dp), intent(out) :: determinant !! Determinant including row-swap sign; zero when singular.
        integer, intent(out) :: info !! Zero on success; nonzero if the matrix is singular or not square.
        real(dp), allocatable :: work(:, :)
        real(dp) :: pivot
        real(dp) :: factor
        real(dp) :: max_value
        real(dp), allocatable :: tmp(:)
        integer :: n
        integer :: i
        integer :: j
        integer :: p
        integer :: sign_det

        n = size(a, 1)
        info = 0
        determinant = 0.0_dp
        if (size(a, 2) /= n .or. size(inverse, 1) /= n .or. size(inverse, 2) /= n) then
            info = 1
            return
        end if
        allocate(work(n, n), tmp(n))
        work = a
        inverse = 0.0_dp
        do i = 1, n
            inverse(i, i) = 1.0_dp
        end do
        determinant = 1.0_dp
        sign_det = 1
        do i = 1, n
            p = i
            max_value = abs(work(i, i))
            do j = i + 1, n
                if (abs(work(j, i)) > max_value) then
                    p = j
                    max_value = abs(work(j, i))
                end if
            end do
            if (max_value <= sqrt(tiny(1.0_dp))) then
                info = 2
                determinant = 0.0_dp
                return
            end if
            if (p /= i) then
                tmp = work(i, :)
                work(i, :) = work(p, :)
                work(p, :) = tmp
                tmp = inverse(i, :)
                inverse(i, :) = inverse(p, :)
                inverse(p, :) = tmp
                sign_det = -sign_det
            end if
            pivot = work(i, i)
            determinant = determinant * pivot
            work(i, :) = work(i, :) / pivot
            inverse(i, :) = inverse(i, :) / pivot
            do j = 1, n
                if (j == i) cycle
                factor = work(j, i)
                work(j, :) = work(j, :) - factor * work(i, :)
                inverse(j, :) = inverse(j, :) - factor * inverse(i, :)
            end do
        end do
        determinant = real(sign_det, dp) * determinant
    end subroutine invert_and_det

    pure function binomial_real(n, k) result(value)
        integer, intent(in) :: n !! Nonnegative upper argument of the binomial coefficient.
        integer, intent(in) :: k !! Lower argument; values outside 0..n return zero.
        real(dp) :: value
        integer :: i
        integer :: kk

        if (k < 0 .or. k > n) then
            value = 0.0_dp
            return
        end if
        kk = min(k, n - k)
        value = 1.0_dp
        do i = 1, kk
            value = value * real(n - kk + i, dp) / real(i, dp)
        end do
    end function binomial_real

end module e1071_fclust_indices
