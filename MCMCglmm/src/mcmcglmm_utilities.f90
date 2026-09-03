! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 computational code; see NOTICE.md and upstream/.
module mcmcglmm_utilities
    use r_kinds, only : dp
    use r_linalg, only : symmetric_eigen
    use mcmcglmm_matrix, only : covariance_to_correlation, kronecker_product
    implicit none
    private

    real(dp), parameter :: radians_to_degrees = 57.295779513082320876798154814105_dp

    public :: triangle_to_matrix
    public :: matrix_to_triangle
    public :: uniform_central_moment
    public :: central_moment_tensor
    public :: krzanowski_compare
    public :: symmetrizer_matrix
    public :: normal_moment_matrix

contains

    pure subroutine triangle_to_matrix(values, lower_triangle, include_diagonal, mirror, matrix_value, info)
        real(dp), intent(in) :: values(:) !! Packed triangular values in Fortran/R column-major triangular traversal order.
        logical, intent(in) :: lower_triangle !! True to fill the lower triangle; false to fill the upper triangle.
        logical, intent(in) :: include_diagonal !! True when packed values include diagonal entries.
        logical, intent(in) :: mirror !! True to copy the populated triangle into the opposite off-diagonal triangle.
        real(dp), allocatable, intent(out) :: matrix_value(:, :) !! Allocated square matrix reconstructed from packed values.
        integer, intent(out) :: info !! Zero on success; nonzero when packed length is not triangular.
        integer :: count
        integer :: i
        integer :: j
        integer :: n
        integer :: needed

        info = 0
        n = 0
        do
            n = n + 1
            if (include_diagonal) then
                needed = n * (n + 1) / 2
            else
                needed = n * (n - 1) / 2
            end if
            if (needed >= size(values)) exit
        end do
        if (needed /= size(values)) then
            allocate(matrix_value(0, 0))
            info = 1
            return
        end if
        allocate(matrix_value(n, n))
        matrix_value = 0.0_dp
        count = 0
        if (lower_triangle) then
            do j = 1, n
                do i = 1, n
                    if (i < j) cycle
                    if (.not. include_diagonal .and. i == j) cycle
                    count = count + 1
                    matrix_value(i, j) = values(count)
                end do
            end do
            if (mirror) then
                do j = 1, n
                    do i = 1, j - 1
                        matrix_value(i, j) = matrix_value(j, i)
                    end do
                end do
            end if
        else
            do j = 1, n
                do i = 1, n
                    if (i > j) cycle
                    if (.not. include_diagonal .and. i == j) cycle
                    count = count + 1
                    matrix_value(i, j) = values(count)
                end do
            end do
            if (mirror) then
                do j = 1, n
                    do i = j + 1, n
                        matrix_value(i, j) = matrix_value(j, i)
                    end do
                end do
            end if
        end if
    end subroutine triangle_to_matrix

    pure subroutine matrix_to_triangle(matrix_value, lower_triangle, include_diagonal, values, info)
        real(dp), intent(in) :: matrix_value(:, :) !! Square matrix to pack by triangular traversal.
        logical, intent(in) :: lower_triangle !! True to pack lower-triangle entries; false for upper-triangle entries.
        logical, intent(in) :: include_diagonal !! True to include diagonal entries in the packed vector.
        real(dp), allocatable, intent(out) :: values(:) !! Allocated packed triangular values in column-major traversal order.
        integer, intent(out) :: info !! Zero on success; nonzero when the input is not square.
        integer :: count
        integer :: i
        integer :: j
        integer :: n
        integer :: nvalue

        info = 0
        n = size(matrix_value, 1)
        if (size(matrix_value, 2) /= n) then
            allocate(values(0))
            info = 1
            return
        end if
        if (include_diagonal) then
            nvalue = n * (n + 1) / 2
        else
            nvalue = n * (n - 1) / 2
        end if
        allocate(values(nvalue))
        count = 0
        do j = 1, n
            do i = 1, n
                if (lower_triangle .and. i < j) cycle
                if (.not. lower_triangle .and. i > j) cycle
                if (.not. include_diagonal .and. i == j) cycle
                count = count + 1
                values(count) = matrix_value(i, j)
            end do
        end do
    end subroutine matrix_to_triangle

    pure elemental real(dp) function uniform_central_moment(minimum, maximum, order) result(value)
        real(dp), intent(in) :: minimum !! Lower endpoint of a continuous uniform distribution.
        real(dp), intent(in) :: maximum !! Upper endpoint of a continuous uniform distribution.
        integer, intent(in) :: order !! Nonnegative central-moment order.

        if (order < 0 .or. maximum < minimum) then
            value = 0.0_dp
        else
            value = ((minimum - maximum) ** order + (maximum - minimum) ** order) / &
                (2.0_dp ** (order + 1) * real(order + 1, dp))
        end if
    end function uniform_central_moment

    pure subroutine central_moment_tensor(x, order, moments, info)
        real(dp), intent(in) :: x(:, :) !! Observation-by-variable data matrix; each variable is centered internally.
        integer, intent(in) :: order !! Positive tensor order k; output contains p^k central moments.
        real(dp), allocatable, intent(out) :: moments(:) !! Flattened column-major p^k Ptensor central moments.
        integer, intent(out) :: info !! Zero on success; nonzero for empty data or nonpositive order.
        real(dp), allocatable :: centered(:, :)
        real(dp), allocatable :: means(:)
        real(dp) :: product_value
        integer :: combination
        integer :: exponent
        integer :: i
        integer :: index_value
        integer :: n
        integer :: p
        integer :: variable

        info = 0
        n = size(x, 1)
        p = size(x, 2)
        if (n < 1 .or. p < 1 .or. order < 1) then
            allocate(moments(0))
            info = 1
            return
        end if
        means = sum(x, dim=1) / real(n, dp)
        allocate(centered(n, p))
        do variable = 1, p
            centered(:, variable) = x(:, variable) - means(variable)
        end do
        allocate(moments(p ** order))
        do combination = 0, p ** order - 1
            moments(combination + 1) = 0.0_dp
            do i = 1, n
                product_value = 1.0_dp
                index_value = combination
                do exponent = 1, order
                    variable = modulo(index_value, p) + 1
                    index_value = index_value / p
                    product_value = product_value * centered(i, variable)
                end do
                moments(combination + 1) = moments(combination + 1) + product_value
            end do
            moments(combination + 1) = moments(combination + 1) / real(n, dp)
        end do
    end subroutine central_moment_tensor

    subroutine krzanowski_compare(covariance_a, covariance_b, vectors_a, vectors_b, use_correlation, &
                                       sum_of_s, angles_degrees, bisectors, info)
        real(dp), intent(in) :: covariance_a(:, :) !! First symmetric covariance or correlation matrix.
        real(dp), intent(in) :: covariance_b(:, :) !! Second symmetric covariance or correlation matrix of the same dimension.
        integer, intent(in) :: vectors_a(:) !! One-based eigenvector positions selected from the first matrix.
        integer, intent(in) :: vectors_b(:) !! One-based eigenvector positions selected from the second matrix.
        logical, intent(in) :: use_correlation !! If true, convert both inputs to correlation matrices before eigendecomposition.
        real(dp), intent(out) :: sum_of_s !! Sum of eigenvalues of the subspace-overlap matrix.
        real(dp), allocatable, intent(out) :: angles_degrees(:) !! Principal angles between selected subspaces in degrees.
        real(dp), allocatable, intent(out) :: bisectors(:, :) !! Normalized Krzanowski bisector vectors, one per selected component.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid dimensions/indices or a failed eigendecomposition.
        real(dp), allocatable :: a_matrix(:, :)
        real(dp), allocatable :: a_values(:)
        real(dp), allocatable :: a_vectors(:, :)
        real(dp), allocatable :: b_matrix(:, :)
        real(dp), allocatable :: b_values(:)
        real(dp), allocatable :: b_vectors(:, :)
        real(dp), allocatable :: ca(:, :)
        real(dp), allocatable :: cb(:, :)
        real(dp), allocatable :: eig_values(:)
        real(dp), allocatable :: eig_vectors(:, :)
        real(dp), allocatable :: overlap(:, :)
        real(dp), allocatable :: projection_b(:, :)
        real(dp), allocatable :: work(:)
        real(dp) :: norm_value
        real(dp) :: root_value
        integer :: i
        integer :: n
        integer :: nv

        info = 0
        n = size(covariance_a, 1)
        nv = size(vectors_a)
        if (size(covariance_a, 2) /= n .or. size(covariance_b, 1) /= n .or. size(covariance_b, 2) /= n .or. &
            size(vectors_b) /= nv .or. nv < 1) then
            allocate(angles_degrees(0), bisectors(0, 0))
            sum_of_s = 0.0_dp
            info = 1
            return
        end if
        if (any(vectors_a < 1) .or. any(vectors_a > n) .or. any(vectors_b < 1) .or. any(vectors_b > n)) then
            allocate(angles_degrees(0), bisectors(0, 0))
            sum_of_s = 0.0_dp
            info = 2
            return
        end if
        a_matrix = covariance_a
        b_matrix = covariance_b
        if (use_correlation) then
            call covariance_to_correlation(covariance_a, a_matrix, info)
            if (info /= 0) return
            call covariance_to_correlation(covariance_b, b_matrix, info)
            if (info /= 0) return
        end if
        call symmetric_eigen(a_matrix, a_values, a_vectors, info, descending=.true.)
        if (info /= 0) return
        call symmetric_eigen(b_matrix, b_values, b_vectors, info, descending=.true.)
        if (info /= 0) return
        ca = a_vectors(:, vectors_a)
        cb = b_vectors(:, vectors_b)
        projection_b = matmul(cb, transpose(cb))
        overlap = matmul(transpose(ca), matmul(projection_b, ca))
        call symmetric_eigen(overlap, eig_values, eig_vectors, info, descending=.true.)
        if (info /= 0) return
        allocate(angles_degrees(nv), bisectors(n, nv), work(n))
        sum_of_s = sum(eig_values)
        do i = 1, nv
            root_value = sqrt(max(eig_values(i), 0.0_dp))
            angles_degrees(i) = acos(min(max(root_value, 0.0_dp), 1.0_dp)) * radians_to_degrees
            if (root_value <= tiny(1.0_dp)) then
                work = matmul(ca, eig_vectors(:, i))
            else
                work = matmul(ca, eig_vectors(:, i)) + matmul(projection_b, matmul(ca, eig_vectors(:, i))) / root_value
            end if
            norm_value = sqrt(dot_product(work, work))
            if (norm_value > 0.0_dp) then
                bisectors(:, i) = work / norm_value
            else
                bisectors(:, i) = 0.0_dp
            end if
        end do
    end subroutine krzanowski_compare

    pure subroutine symmetrizer_matrix(dimension, order, matrix_value, info)
        integer, intent(in) :: dimension !! Tensor mode dimension m; must be positive.
        integer, intent(in) :: order !! Tensor order k to symmetrize; must be nonnegative and computationally tractable.
        real(dp), allocatable, intent(out) :: matrix_value(:, :) !! m^k by m^k KPPM symmetrizer matrix.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid or impractically large dimensions.
        integer, allocatable :: current(:)
        integer, allocatable :: permutations(:, :)
        logical, allocatable :: used(:)
        integer :: column
        integer :: count
        integer :: factorial_order
        integer :: index_value
        integer :: n_state
        integer :: permutation
        integer :: row
        integer, allocatable :: tuple(:)
        integer, allocatable :: tuple_permuted(:)

        info = 0
        if (dimension < 1 .or. order < 0) then
            allocate(matrix_value(0, 0))
            info = 1
            return
        end if
        call safe_integer_power(dimension, order, n_state, info)
        if (info /= 0) then
            allocate(matrix_value(0, 0))
            return
        end if
        call safe_factorial(order, factorial_order, info)
        if (info /= 0) then
            allocate(matrix_value(0, 0))
            return
        end if
        if (n_state > 4096 .or. factorial_order > 40320) then
            allocate(matrix_value(0, 0))
            info = 2
            return
        end if
        allocate(matrix_value(n_state, n_state))
        matrix_value = 0.0_dp
        if (order == 0) then
            matrix_value(1, 1) = 1.0_dp
            return
        end if
        allocate(permutations(order, factorial_order), current(order), used(order))
        count = 0
        used = .false.
        current = 0
        call enumerate_permutations(order, 1, used, current, permutations, count)
        allocate(tuple(order), tuple_permuted(order))
        do column = 1, n_state
            index_value = column - 1
            call decode_tensor_index(index_value, dimension, tuple)
            do permutation = 1, factorial_order
                tuple_permuted = tuple(permutations(:, permutation))
                row = encode_tensor_index(tuple_permuted, dimension) + 1
                matrix_value(row, column) = matrix_value(row, column) + 1.0_dp / real(factorial_order, dp)
            end do
        end do
    end subroutine symmetrizer_matrix

    pure subroutine normal_moment_matrix(covariance, order, moments, info)
        real(dp), intent(in) :: covariance(:, :) !! Symmetric covariance matrix V of a zero-mean multivariate normal vector.
        integer, intent(in) :: order !! Even tensor order requested by MCMCglmm knorm; must be positive.
        real(dp), allocatable, intent(out) :: moments(:, :) !! Flattened m^(order/2) square normal central-moment tensor.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid dimensions/order or excessive tensor size.
        real(dp), allocatable :: kron_component(:, :)
        real(dp), allocatable :: mm(:, :)
        real(dp), allocatable :: outer_covariance(:, :)
        real(dp), allocatable :: symmetrizer(:, :)
        real(dp), allocatable :: vector_covariance(:)
        real(dp) :: coefficient
        integer :: half_order
        integer :: h
        integer :: ir
        integer :: is
        integer :: m
        integer :: n_state
        integer :: status

        info = 0
        m = size(covariance, 1)
        if (m < 1 .or. size(covariance, 2) /= m .or. order < 2 .or. modulo(order, 2) /= 0) then
            allocate(moments(0, 0))
            info = 1
            return
        end if
        half_order = order / 2
        call safe_integer_power(m, half_order, n_state, info)
        if (info /= 0 .or. n_state > 4096) then
            allocate(moments(0, 0))
            if (info == 0) info = 2
            return
        end if
        call symmetrizer_matrix(m, half_order, symmetrizer, status)
        if (status /= 0) then
            allocate(moments(0, 0))
            info = status
            return
        end if
        allocate(mm(n_state, n_state))
        mm = 0.0_dp
        vector_covariance = reshape(covariance, [m * m])
        allocate(outer_covariance(m * m, m * m))
        outer_covariance = spread(vector_covariance, 2, m * m) * spread(vector_covariance, 1, m * m)
        do h = 1, half_order / 2 + 1
            is = h - 1
            ir = half_order - 2 * is
            coefficient = factorial_real(half_order) ** 2 / &
                (factorial_real(ir) * factorial_real(is) ** 2 * 2.0_dp ** (2 * is))
            call matrix_kronecker_power(covariance, ir, kron_component, status)
            if (status /= 0) then
                allocate(moments(0, 0))
                info = status
                return
            end if
            if (is > 0) then
                call append_kronecker_power(kron_component, outer_covariance, is, status)
                if (status /= 0) then
                    allocate(moments(0, 0))
                    info = status
                    return
                end if
            end if
            mm = mm + coefficient * kron_component
        end do
        moments = matmul(symmetrizer, matmul(mm, symmetrizer))
    end subroutine normal_moment_matrix

    pure recursive subroutine enumerate_permutations(order, position, used, current, permutations, count)
        integer, intent(in) :: order !! Number of positions in each permutation.
        integer, intent(in) :: position !! One-based position currently being filled.
        logical, intent(inout) :: used(:) !! Flags indicating which source positions are already present.
        integer, intent(inout) :: current(:) !! Partially constructed permutation.
        integer, intent(out) :: permutations(:, :) !! Completed permutations stored columnwise.
        integer, intent(inout) :: count !! Number of completed permutations written so far.
        integer :: candidate

        if (position > order) then
            count = count + 1
            permutations(:, count) = current
            return
        end if
        do candidate = 1, order
            if (used(candidate)) cycle
            used(candidate) = .true.
            current(position) = candidate
            call enumerate_permutations(order, position + 1, used, current, permutations, count)
            used(candidate) = .false.
        end do
    end subroutine enumerate_permutations

    pure subroutine decode_tensor_index(index_value, dimension, tuple)
        integer, intent(in) :: index_value !! Zero-based flattened tensor index.
        integer, intent(in) :: dimension !! Common tensor mode dimension.
        integer, intent(out) :: tuple(:) !! One-based tensor indices, fastest varying in the first mode.
        integer :: i
        integer :: remainder

        remainder = index_value
        do i = 1, size(tuple)
            tuple(i) = modulo(remainder, dimension) + 1
            remainder = remainder / dimension
        end do
    end subroutine decode_tensor_index

    pure integer function encode_tensor_index(tuple, dimension) result(index_value)
        integer, intent(in) :: tuple(:) !! One-based tensor indices, fastest varying in the first mode.
        integer, intent(in) :: dimension !! Common tensor mode dimension.
        integer :: i
        integer :: multiplier

        index_value = 0
        multiplier = 1
        do i = 1, size(tuple)
            index_value = index_value + (tuple(i) - 1) * multiplier
            multiplier = multiplier * dimension
        end do
    end function encode_tensor_index

    pure subroutine safe_integer_power(base, exponent, value, info)
        integer, intent(in) :: base !! Positive integer base.
        integer, intent(in) :: exponent !! Nonnegative integer exponent.
        integer, intent(out) :: value !! base raised to exponent when representable.
        integer, intent(out) :: info !! Zero on success; nonzero on invalid input or integer overflow.
        integer :: i

        info = 0
        if (base < 1 .or. exponent < 0) then
            value = 0
            info = 1
            return
        end if
        value = 1
        do i = 1, exponent
            if (value > huge(value) / base) then
                value = 0
                info = 2
                return
            end if
            value = value * base
        end do
    end subroutine safe_integer_power

    pure subroutine safe_factorial(order, value, info)
        integer, intent(in) :: order !! Nonnegative factorial argument.
        integer, intent(out) :: value !! order factorial when representable.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid input or integer overflow.
        integer :: i

        info = 0
        if (order < 0) then
            value = 0
            info = 1
            return
        end if
        value = 1
        do i = 2, order
            if (value > huge(value) / i) then
                value = 0
                info = 2
                return
            end if
            value = value * i
        end do
    end subroutine safe_factorial

    pure real(dp) function factorial_real(order) result(value)
        integer, intent(in) :: order !! Nonnegative factorial argument represented in real arithmetic.
        integer :: i

        value = 1.0_dp
        do i = 2, order
            value = value * real(i, dp)
        end do
    end function factorial_real

    pure subroutine matrix_kronecker_power(matrix_value, exponent, result_matrix, info)
        real(dp), intent(in) :: matrix_value(:, :) !! Square matrix forming the Kronecker power.
        integer, intent(in) :: exponent !! Nonnegative Kronecker exponent.
        real(dp), allocatable, intent(out) :: result_matrix(:, :) !! Kronecker power, with exponent zero represented by [1].
        integer, intent(out) :: info !! Zero on success; nonzero for a nonsquare matrix or excessive allocation size.
        real(dp), allocatable :: next_matrix(:, :)
        integer :: i
        integer :: n

        info = 0
        n = size(matrix_value, 1)
        if (n < 1 .or. size(matrix_value, 2) /= n .or. exponent < 0) then
            allocate(result_matrix(0, 0))
            info = 1
            return
        end if
        allocate(result_matrix(1, 1))
        result_matrix(1, 1) = 1.0_dp
        do i = 1, exponent
            call kronecker_product(result_matrix, matrix_value, next_matrix)
            call move_alloc(next_matrix, result_matrix)
        end do
    end subroutine matrix_kronecker_power

    pure subroutine append_kronecker_power(left_matrix, factor, exponent, info)
        real(dp), allocatable, intent(inout) :: left_matrix(:, :) !! Left Kronecker factor, replaced by the expanded product.
        real(dp), intent(in) :: factor(:, :) !! Square matrix appended repeatedly on the right.
        integer, intent(in) :: exponent !! Number of copies of factor to append.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid dimensions or exponent.
        real(dp), allocatable :: next_matrix(:, :)
        integer :: i

        info = 0
        if (.not. allocated(left_matrix) .or. size(factor, 1) < 1 .or. size(factor, 1) /= size(factor, 2) .or. exponent < 0) then
            info = 1
            return
        end if
        do i = 1, exponent
            call kronecker_product(left_matrix, factor, next_matrix)
            call move_alloc(next_matrix, left_matrix)
        end do
    end subroutine append_kronecker_power

end module mcmcglmm_utilities
