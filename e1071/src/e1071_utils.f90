module e1071_utils
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
    use e1071_kinds, only: dp
    use e1071_constants, only: e1071_pi
    use e1071_rng, only: rng_state, rng_normal
    implicit none
    private

    type, public :: interpolation_axis
        real(dp), allocatable :: values(:)
    end type interpolation_axis

    public :: e1071_moment, e1071_skewness, e1071_kurtosis
    public :: sigmoid, dsigmoid, d2sigmoid
    public :: hanning_window, hamming_window, rectangle_window
    public :: hamming_distance_vector, hamming_distance_matrix
    public :: bincombinations, permutations
    public :: impute_columns, scale_matrix
    public :: interpolate_nd
    public :: rwiener, rbridge

contains

    function e1071_moment(x, order, center, absolute, na_rm) result(value)
        real(dp), intent(in) :: x(:) !! Sample values whose raw or centered moment is requested; NaN represents missing data.
        integer, intent(in), optional :: order !! Moment order; defaults to one and may be any nonnegative integer.
        logical, intent(in), optional :: center !! If true, subtract the complete-case mean before exponentiation; default false.
        logical, intent(in), optional :: absolute !! If true, take absolute values before exponentiation; default false.
        logical, intent(in), optional :: na_rm !! If true, omit NaN values; otherwise any NaN makes the result NaN.
        real(dp) :: value
        integer :: p
        integer :: i
        integer :: n
        logical :: do_center
        logical :: do_absolute
        logical :: remove_na
        real(dp) :: mu
        real(dp) :: z

        p = 1
        if (present(order)) p = order
        do_center = .false.
        if (present(center)) do_center = center
        do_absolute = .false.
        if (present(absolute)) do_absolute = absolute
        remove_na = .false.
        if (present(na_rm)) remove_na = na_rm
        if (p < 0) error stop "e1071_moment: order must be nonnegative"

        if (.not. remove_na .and. any(ieee_is_nan(x))) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if

        mu = 0.0_dp
        n = 0
        if (do_center) then
            do i = 1, size(x)
                if (ieee_is_nan(x(i))) cycle
                mu = mu + x(i)
                n = n + 1
            end do
            if (n == 0) then
                value = ieee_value(0.0_dp, ieee_quiet_nan)
                return
            end if
            mu = mu / real(n, dp)
        end if

        value = 0.0_dp
        n = 0
        do i = 1, size(x)
            if (ieee_is_nan(x(i))) cycle
            z = x(i) - mu
            if (do_absolute) z = abs(z)
            value = value + z**p
            n = n + 1
        end do
        if (n == 0) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
        else
            value = value / real(n, dp)
        end if
    end function e1071_moment

    function e1071_skewness(x, type, na_rm) result(value)
        real(dp), intent(in) :: x(:) !! Sample values used to compute skewness; NaN values are missing observations.
        integer, intent(in), optional :: type !! e1071 skewness convention 1, 2, or 3; defaults to type 3.
        logical, intent(in), optional :: na_rm !! If true, omit NaNs; otherwise a NaN input yields NaN.
        real(dp) :: value
        integer :: kind
        integer :: i
        integer :: n
        logical :: remove_na
        real(dp) :: mu
        real(dp) :: s2
        real(dp) :: s3

        kind = 3
        if (present(type)) kind = type
        if (kind < 1 .or. kind > 3) error stop "e1071_skewness: type must be 1, 2, or 3"
        remove_na = .false.
        if (present(na_rm)) remove_na = na_rm
        if (.not. remove_na .and. any(ieee_is_nan(x))) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        call complete_mean(x, mu, n)
        if (n == 0) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        s2 = 0.0_dp
        s3 = 0.0_dp
        do i = 1, size(x)
            if (ieee_is_nan(x(i))) cycle
            s2 = s2 + (x(i) - mu)**2
            s3 = s3 + (x(i) - mu)**3
        end do
        if (s2 <= 0.0_dp) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        value = sqrt(real(n, dp)) * s3 / s2**1.5_dp
        if (kind == 2) then
            if (n < 3) error stop "e1071_skewness: type 2 needs at least three observations"
            value = value * sqrt(real(n * (n - 1), dp)) / real(n - 2, dp)
        else if (kind == 3) then
            value = value * (1.0_dp - 1.0_dp / real(n, dp))**1.5_dp
        end if
    end function e1071_skewness

    function e1071_kurtosis(x, type, na_rm) result(value)
        real(dp), intent(in) :: x(:) !! Sample values used to compute excess kurtosis; NaN values are missing observations.
        integer, intent(in), optional :: type !! e1071 kurtosis convention 1, 2, or 3; defaults to type 3.
        logical, intent(in), optional :: na_rm !! If true, omit NaNs; otherwise a NaN input yields NaN.
        real(dp) :: value
        integer :: kind
        integer :: i
        integer :: n
        logical :: remove_na
        real(dp) :: mu
        real(dp) :: s2
        real(dp) :: s4
        real(dp) :: r

        kind = 3
        if (present(type)) kind = type
        if (kind < 1 .or. kind > 3) error stop "e1071_kurtosis: type must be 1, 2, or 3"
        remove_na = .false.
        if (present(na_rm)) remove_na = na_rm
        if (.not. remove_na .and. any(ieee_is_nan(x))) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        call complete_mean(x, mu, n)
        if (n == 0) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        s2 = 0.0_dp
        s4 = 0.0_dp
        do i = 1, size(x)
            if (ieee_is_nan(x(i))) cycle
            s2 = s2 + (x(i) - mu)**2
            s4 = s4 + (x(i) - mu)**4
        end do
        if (s2 <= 0.0_dp) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        r = real(n, dp) * s4 / (s2 * s2)
        if (kind == 1) then
            value = r - 3.0_dp
        else if (kind == 2) then
            if (n < 4) error stop "e1071_kurtosis: type 2 needs at least four observations"
            value = ((real(n + 1, dp) * (r - 3.0_dp) + 6.0_dp) * real(n - 1, dp)) &
                / real((n - 2) * (n - 3), dp)
        else
            value = r * (1.0_dp - 1.0_dp / real(n, dp))**2 - 3.0_dp
        end if
    end function e1071_kurtosis

    elemental function sigmoid(x) result(value)
        real(dp), intent(in) :: x !! Real argument of the logistic sigmoid function.
        real(dp) :: value

        if (x >= 0.0_dp) then
            value = 1.0_dp / (1.0_dp + exp(-x))
        else
            value = exp(x) / (1.0_dp + exp(x))
        end if
    end function sigmoid

    elemental function dsigmoid(x) result(value)
        real(dp), intent(in) :: x !! Real argument at which the first logistic derivative is evaluated.
        real(dp) :: value
        real(dp) :: s

        s = sigmoid(x)
        value = s * (1.0_dp - s)
    end function dsigmoid

    elemental function d2sigmoid(x) result(value)
        real(dp), intent(in) :: x !! Real argument at which the second logistic derivative is evaluated.
        real(dp) :: value
        real(dp) :: s

        s = sigmoid(x)
        value = s * (1.0_dp - s) * (1.0_dp - 2.0_dp * s)
    end function d2sigmoid

    function hanning_window(n) result(values)
        integer, intent(in) :: n !! Number of window coefficients; must be positive.
        real(dp), allocatable :: values(:)
        integer :: i

        if (n <= 0) error stop "hanning_window: n must be positive"
        allocate(values(n))
        if (n == 1) then
            values = 1.0_dp
        else
            do i = 1, n
                values(i) = 0.5_dp - 0.5_dp * cos(2.0_dp * e1071_pi * real(i - 1, dp) / real(n - 1, dp))
            end do
        end if
    end function hanning_window

    function hamming_window(n) result(values)
        integer, intent(in) :: n !! Number of window coefficients; must be positive.
        real(dp), allocatable :: values(:)
        integer :: i

        if (n <= 0) error stop "hamming_window: n must be positive"
        allocate(values(n))
        if (n == 1) then
            values = 1.0_dp
        else
            do i = 1, n
                values(i) = 0.54_dp - 0.46_dp * cos(2.0_dp * e1071_pi * real(i - 1, dp) / real(n - 1, dp))
            end do
        end if
    end function hamming_window

    function rectangle_window(n) result(values)
        integer, intent(in) :: n !! Number of unit window coefficients; must be nonnegative.
        real(dp), allocatable :: values(:)

        if (n < 0) error stop "rectangle_window: n must be nonnegative"
        allocate(values(n))
        values = 1.0_dp
    end function rectangle_window

    function hamming_distance_vector(x, y) result(distance)
        integer, intent(in) :: x(:) !! First integer or encoded-categorical vector.
        integer, intent(in) :: y(:) !! Second vector; must have the same length as x.
        integer :: distance

        if (size(x) /= size(y)) error stop "hamming_distance_vector: size mismatch"
        distance = count(x /= y)
    end function hamming_distance_vector

    function hamming_distance_matrix(x) result(distance)
        integer, intent(in) :: x(:, :) !! Observations by variables matrix of integer or encoded-categorical values.
        integer, allocatable :: distance(:, :)
        integer :: i
        integer :: j

        allocate(distance(size(x, 1), size(x, 1)))
        distance = 0
        do i = 1, size(x, 1) - 1
            do j = i + 1, size(x, 1)
                distance(i, j) = count(x(i, :) /= x(j, :))
                distance(j, i) = distance(i, j)
            end do
        end do
    end function hamming_distance_matrix

    function bincombinations(p) result(values)
        integer, intent(in) :: p !! Number of binary columns; must be nonnegative and small enough that 2**p fits an integer.
        integer, allocatable :: values(:, :)
        integer :: nrow
        integer :: i
        integer :: j
        integer :: block

        if (p < 0) error stop "bincombinations: p must be nonnegative"
        nrow = 2**p
        allocate(values(nrow, p))
        do j = 1, p
            block = 2**(p - j)
            do i = 1, nrow
                values(i, j) = modulo((i - 1) / block, 2)
            end do
        end do
    end function bincombinations

    function permutations(n) result(values)
        integer, intent(in) :: n !! Number of labels to permute; must be a positive integer and is practical only for small n.
        integer, allocatable :: values(:, :)
        integer, allocatable :: work(:)
        integer :: nperm
        integer :: row

        if (n < 1) error stop "permutations: n must be positive"
        nperm = factorial_int(n)
        allocate(values(nperm, n))
        allocate(work(n))
        work = [(row, row = 1, n)]
        row = 0
        call enumerate_permutations(work, 1, values, row)
    end function permutations

    subroutine impute_columns(x, method)
        real(dp), intent(inout) :: x(:, :) !! Matrix modified in place; NaN entries are replaced column by column.
        character(len=*), intent(in), optional :: method !! Replacement statistic, either "median" or "mean"; defaults to median.
        character(len=:), allocatable :: use_method
        integer :: j
        real(dp) :: replacement

        use_method = "median"
        if (present(method)) use_method = trim(adjustl(method))
        do j = 1, size(x, 2)
            if (use_method == "mean") then
                replacement = complete_column_mean(x(:, j))
            else if (use_method == "median") then
                replacement = complete_column_median(x(:, j))
            else
                error stop "impute_columns: method must be mean or median"
            end if
            where (ieee_is_nan(x(:, j))) x(:, j) = replacement
        end do
    end subroutine impute_columns

    subroutine scale_matrix(x, center, scale, means, scales, mask)
        real(dp), intent(inout) :: x(:, :) !! Observation-by-variable matrix scaled in place.
        logical, intent(in), optional :: center !! If true, subtract complete-case column means; defaults to true.
        logical, intent(in), optional :: scale !! If true, divide by sample standard deviations; defaults to true.
        real(dp), intent(out), optional :: means(:) !! Returned centering constants; length must equal the number of columns.
        real(dp), intent(out), optional :: scales(:) !! Returned scale constants; length must equal the number of columns.
        logical, intent(in), optional :: mask(:) !! Per-column selector; false columns are left unchanged.
        logical :: do_center
        logical :: do_scale
        logical :: selected
        integer :: i
        integer :: j
        integer :: n
        real(dp) :: mu
        real(dp) :: ss
        real(dp) :: sd

        do_center = .true.
        if (present(center)) do_center = center
        do_scale = .true.
        if (present(scale)) do_scale = scale
        if (present(means)) then
            if (size(means) /= size(x, 2)) error stop "scale_matrix: means has wrong length"
            means = 0.0_dp
        end if
        if (present(scales)) then
            if (size(scales) /= size(x, 2)) error stop "scale_matrix: scales has wrong length"
            scales = 1.0_dp
        end if
        if (present(mask)) then
            if (size(mask) /= size(x, 2)) error stop "scale_matrix: mask has wrong length"
        end if

        do j = 1, size(x, 2)
            selected = .true.
            if (present(mask)) selected = mask(j)
            if (.not. selected) cycle
            call complete_mean(x(:, j), mu, n)
            if (n == 0) cycle
            ss = 0.0_dp
            do i = 1, size(x, 1)
                if (ieee_is_nan(x(i, j))) cycle
                ss = ss + (x(i, j) - mu)**2
            end do
            if (n > 1) then
                sd = sqrt(ss / real(n - 1, dp))
            else
                sd = 0.0_dp
            end if
            if (present(means)) means(j) = mu
            if (present(scales)) scales(j) = sd
            if (do_center) then
                do i = 1, size(x, 1)
                    if (.not. ieee_is_nan(x(i, j))) x(i, j) = x(i, j) - mu
                end do
            end if
            if (do_scale) then
                if (sd <= 0.0_dp) error stop "scale_matrix: cannot scale a constant column"
                do i = 1, size(x, 1)
                    if (.not. ieee_is_nan(x(i, j))) x(i, j) = x(i, j) / sd
                end do
            end if
        end do
    end subroutine scale_matrix

    function interpolate_nd(points, data, dims, axes, method) result(values)
        real(dp), intent(in) :: points(:, :) !! Query points; rows are points and columns correspond to array dimensions.
        real(dp), intent(in) :: data(:) !! Flattened Fortran-column-major data array with product(dims) elements.
        integer, intent(in) :: dims(:) !! Length of each array dimension; product must equal size(data).
        type(interpolation_axis), intent(in) :: axes(:) !! Ordered coordinate values for each dimension; axis lengths equal dims.
        character(len=*), intent(in), optional :: method !! Interpolation rule: "linear" or "constant"; default linear.
        real(dp), allocatable :: values(:)
        integer, allocatable :: left(:)
        real(dp), allocatable :: dx(:)
        real(dp), allocatable :: width(:)
        integer :: npoint
        integer :: nd
        integer :: i
        integer :: j
        integer :: corner
        integer :: index
        integer :: bit
        real(dp) :: weight
        character(len=:), allocatable :: use_method

        nd = size(dims)
        if (size(points, 2) /= nd .or. size(axes) /= nd) error stop "interpolate_nd: dimension mismatch"
        if (product(dims) /= size(data)) error stop "interpolate_nd: data size does not match dims"
        do j = 1, nd
            if (size(axes(j)%values) /= dims(j)) error stop "interpolate_nd: axis length does not match dims"
            if (dims(j) < 2) error stop "interpolate_nd: every dimension must contain at least two grid points"
            if (any(axes(j)%values(2:) < axes(j)%values(:dims(j) - 1))) error stop "interpolate_nd: axes must be sorted"
        end do
        use_method = "linear"
        if (present(method)) use_method = trim(adjustl(method))
        if (use_method /= "linear" .and. use_method /= "constant") error stop "interpolate_nd: invalid method"

        npoint = size(points, 1)
        allocate(values(npoint), left(nd), dx(nd), width(nd))
        do i = 1, npoint
            do j = 1, nd
                if (points(i, j) < axes(j)%values(1) .or. points(i, j) > axes(j)%values(dims(j))) then
                    error stop "interpolate_nd: extrapolation is not allowed"
                end if
                left(j) = 1
                do while (left(j) < dims(j) .and. axes(j)%values(left(j) + 1) <= points(i, j))
                    left(j) = left(j) + 1
                end do
                if (left(j) == dims(j)) left(j) = left(j) - 1
                dx(j) = points(i, j) - axes(j)%values(left(j))
                width(j) = axes(j)%values(left(j) + 1) - axes(j)%values(left(j))
                if (width(j) <= 0.0_dp) error stop "interpolate_nd: duplicate axis coordinates are not allowed"
            end do
            if (use_method == "constant") then
                values(i) = data(flat_index(left, dims))
            else
                values(i) = 0.0_dp
                do corner = 0, 2**nd - 1
                    index = 1
                    weight = 1.0_dp
                    do j = 1, nd
                        bit = modulo(corner / 2**(j - 1), 2)
                        index = index + (left(j) - 1 + bit) * stride_for_dimension(dims, j)
                        if (bit == 0) then
                            weight = weight * (1.0_dp - dx(j) / width(j))
                        else
                            weight = weight * dx(j) / width(j)
                        end if
                    end do
                    values(i) = values(i) + weight * data(index)
                end do
            end if
        end do
    end function interpolate_nd

    subroutine rwiener(end_time, frequency, rng, path)
        real(dp), intent(in) :: end_time !! Simulated time horizon in the same time units as 1/frequency; must be positive.
        integer, intent(in) :: frequency !! Number of Wiener increments per unit time; must be positive.
        type(rng_state), intent(inout) :: rng !! Mutable pseudo-random-number generator state used for Gaussian increments.
        real(dp), allocatable, intent(out) :: path(:) !! Cumulative Wiener path of length nint(end_time*frequency).
        integer :: n
        integer :: i
        real(dp) :: cumulative

        if (end_time <= 0.0_dp .or. frequency <= 0) error stop "rwiener: invalid horizon or frequency"
        n = nint(end_time * real(frequency, dp))
        allocate(path(n))
        cumulative = 0.0_dp
        do i = 1, n
            cumulative = cumulative + rng_normal(rng) / sqrt(real(frequency, dp))
            path(i) = cumulative
        end do
    end subroutine rwiener

    subroutine rbridge(end_time, frequency, rng, path)
        real(dp), intent(in) :: end_time !! Simulated bridge horizon; must be positive.
        integer, intent(in) :: frequency !! Number of increments per unit time; must be positive.
        type(rng_state), intent(inout) :: rng !! Mutable pseudo-random-number generator state used for Gaussian increments.
        real(dp), allocatable, intent(out) :: path(:) !! Brownian bridge samples on the regular grid excluding time zero.
        real(dp), allocatable :: w(:)
        integer :: i
        integer :: n
        real(dp) :: endpoint
        real(dp) :: time_i

        call rwiener(end_time, frequency, rng, w)
        n = size(w)
        endpoint = w(n)
        allocate(path(n))
        do i = 1, n
            time_i = real(i, dp) / real(frequency, dp)
            path(i) = w(i) - (time_i / end_time) * endpoint
        end do
    end subroutine rbridge

    subroutine complete_mean(x, mean_value, n_complete)
        real(dp), intent(in) :: x(:) !! Vector whose finite non-NaN observations are averaged.
        real(dp), intent(out) :: mean_value !! Arithmetic mean of non-NaN observations, or zero if none are present.
        integer, intent(out) :: n_complete !! Number of non-NaN observations contributing to the mean.
        integer :: i

        mean_value = 0.0_dp
        n_complete = 0
        do i = 1, size(x)
            if (ieee_is_nan(x(i))) cycle
            mean_value = mean_value + x(i)
            n_complete = n_complete + 1
        end do
        if (n_complete > 0) mean_value = mean_value / real(n_complete, dp)
    end subroutine complete_mean

    function complete_column_mean(x) result(value)
        real(dp), intent(in) :: x(:) !! Column values whose non-NaN arithmetic mean is requested.
        real(dp) :: value
        integer :: n

        call complete_mean(x, value, n)
        if (n == 0) value = ieee_value(0.0_dp, ieee_quiet_nan)
    end function complete_column_mean

    function complete_column_median(x) result(value)
        real(dp), intent(in) :: x(:) !! Column values whose non-NaN median is requested.
        real(dp) :: value
        real(dp), allocatable :: work(:)
        integer :: i
        integer :: n
        integer :: j
        real(dp) :: key

        n = count(.not. ieee_is_nan(x))
        if (n == 0) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        allocate(work(n))
        j = 0
        do i = 1, size(x)
            if (ieee_is_nan(x(i))) cycle
            j = j + 1
            work(j) = x(i)
        end do
        do i = 2, n
            key = work(i)
            j = i - 1
            do while (j >= 1)
                if (work(j) <= key) exit
                work(j + 1) = work(j)
                j = j - 1
            end do
            work(j + 1) = key
        end do
        if (modulo(n, 2) == 1) then
            value = work((n + 1) / 2)
        else
            value = 0.5_dp * (work(n / 2) + work(n / 2 + 1))
        end if
    end function complete_column_median

    recursive subroutine enumerate_permutations(work, first, values, row)
        integer, intent(inout) :: work(:) !! Current permutation workspace modified temporarily during recursive enumeration.
        integer, intent(in) :: first !! One-based position whose value is being chosen in this recursion level.
        integer, intent(out) :: values(:, :) !! Output matrix receiving one permutation per row.
        integer, intent(inout) :: row !! Number of output permutations already written; incremented by this routine.
        integer :: i
        integer :: tmp

        if (first == size(work)) then
            row = row + 1
            values(row, :) = work
            return
        end if
        do i = first, size(work)
            tmp = work(first)
            work(first) = work(i)
            work(i) = tmp
            call enumerate_permutations(work, first + 1, values, row)
            tmp = work(first)
            work(first) = work(i)
            work(i) = tmp
        end do
    end subroutine enumerate_permutations

    pure function factorial_int(n) result(value)
        integer, intent(in) :: n !! Nonnegative integer whose factorial is requested.
        integer :: value
        integer :: i

        value = 1
        do i = 2, n
            value = value * i
        end do
    end function factorial_int

    pure function flat_index(indices, dims) result(index)
        integer, intent(in) :: indices(:) !! One-based multidimensional subscripts, one per dimension.
        integer, intent(in) :: dims(:) !! Array dimensions corresponding to indices.
        integer :: index
        integer :: j
        integer :: stride

        index = 1
        stride = 1
        do j = 1, size(dims)
            index = index + (indices(j) - 1) * stride
            stride = stride * dims(j)
        end do
    end function flat_index

    pure function stride_for_dimension(dims, dimension) result(stride)
        integer, intent(in) :: dims(:) !! Full array dimension vector.
        integer, intent(in) :: dimension !! One-based dimension whose Fortran column-major stride is requested.
        integer :: stride
        integer :: j

        stride = 1
        do j = 1, dimension - 1
            stride = stride * dims(j)
        end do
    end function stride_for_dimension

end module e1071_utils
