module proxy_gower
    use proxy_kinds, only: dp
    use proxy_ieee, only: proxy_nan, proxy_is_missing
    implicit none
    private

    integer, parameter, public :: proxy_gower_logical = 1
    integer, parameter, public :: proxy_gower_factor = 2
    integer, parameter, public :: proxy_gower_metric = 3
    integer, parameter, public :: proxy_gower_ordinal = 4

    public :: gower_auto_similarity, gower_cross_similarity, gower_pair_similarity

contains

    subroutine gower_auto_similarity(x, types, similarity, weights, ranges, minima)
        real(dp), intent(in) :: x(:, :) !! Observation-by-variable matrix; logical/factor/ordinal variables are numerically
        !! encoded according to `types`, with NaN for missing values.
        integer, intent(in) :: types(:) !! Variable type codes of length `size(x,2)` using `proxy_gower_logical`, `factor`,
        !! `metric`, or `ordinal`.
        real(dp), allocatable, intent(out) :: similarity(:, :) !! Allocated symmetric `size(x,1)` square matrix of Gower
        !! similarities with unit diagonal when defined.
        real(dp), intent(in), optional :: weights(:) !! Optional variable weights recycled to the number of columns; defaults to
        !! equal weights of one.
        real(dp), intent(in), optional :: ranges(:) !! Optional metric scaling ranges recycled across metric/ordinal columns;
        !! zero ranges are treated as one.
        real(dp), intent(in), optional :: minima(:) !! Optional metric scaling minima recycled across metric/ordinal columns;
        !! defaults to observed column minima.
        real(dp), allocatable :: z(:, :)
        real(dp), allocatable :: w(:)
        integer :: i
        integer :: j

        call prepare_gower_single(x, types, z, w, weights, ranges, minima)
        allocate(similarity(size(x, 1), size(x, 1)))
        do j = 1, size(x, 1)
            do i = 1, size(x, 1)
                similarity(i, j) = gower_pair_similarity(z(i, :), z(j, :), types, w)
            end do
        end do
    end subroutine gower_auto_similarity

    subroutine gower_cross_similarity(x, y, types, similarity, pairwise, weights, ranges_x, minima_x, ranges_y, minima_y)
        real(dp), intent(in) :: x(:, :) !! First observation-by-variable matrix, using numeric encodings specified by `types` and
        !! NaN for missing data.
        real(dp), intent(in) :: y(:, :) !! Second observation-by-variable matrix with the same number and interpretation of
        !! variables as `x`.
        integer, intent(in) :: types(:) !! Variable type codes of length `size(x,2)` shared by `x` and `y`.
        real(dp), allocatable, intent(out) :: similarity(:, :) !! Cross-similarity matrix, or an `n-by-1` matrix when `pairwise`
        !! is true.
        logical, intent(in), optional :: pairwise !! If true, compare corresponding rows only and require equal row counts;
        !! default is the full Cartesian cross-matrix.
        real(dp), intent(in), optional :: weights(:) !! Optional variable weights recycled to all variables; defaults to one for
        !! every variable.
        real(dp), intent(in), optional :: ranges_x(:) !! Optional scaling ranges for metric/ordinal columns of `x`; absent values
        !! use the combined `x`/`y` range.
        real(dp), intent(in), optional :: minima_x(:) !! Optional scaling minima for metric/ordinal columns of `x`; absent values
        !! use combined minima.
        real(dp), intent(in), optional :: ranges_y(:) !! Optional scaling ranges for metric/ordinal columns of `y`; absent values
        !! normally reuse the combined range used for `x`.
        real(dp), intent(in), optional :: minima_y(:) !! Optional scaling minima for metric/ordinal columns of `y`; absent values
        !! normally reuse the combined minimum used for `x`.
        real(dp), allocatable :: zx(:, :)
        real(dp), allocatable :: zy(:, :)
        real(dp), allocatable :: w(:)
        logical :: do_pairwise
        integer :: i
        integer :: j

        do_pairwise = .false.
        if (present(pairwise)) do_pairwise = pairwise
        call prepare_gower_pair(x, y, types, zx, zy, w, weights, ranges_x, minima_x, ranges_y, minima_y)
        if (do_pairwise) then
            if (size(x, 1) /= size(y, 1)) then
                allocate(similarity(0, 0))
                return
            end if
            allocate(similarity(size(x, 1), 1))
            do i = 1, size(x, 1)
                similarity(i, 1) = gower_pair_similarity(zx(i, :), zy(i, :), types, w)
            end do
        else
            allocate(similarity(size(x, 1), size(y, 1)))
            do j = 1, size(y, 1)
                do i = 1, size(x, 1)
                    similarity(i, j) = gower_pair_similarity(zx(i, :), zy(j, :), types, w)
                end do
            end do
        end if
    end subroutine gower_cross_similarity

    pure function gower_pair_similarity(x, y, types, weights) result(similarity)
        real(dp), intent(in) :: x(:) !! First pre-scaled Gower observation; metric/ordinal components should already be on the
        !! intended scale.
        real(dp), intent(in) :: y(:) !! Second pre-scaled observation aligned variable-by-variable with `x`.
        integer, intent(in) :: types(:) !! Variable type codes identifying logical, factor, metric, and ordinal components.
        real(dp), intent(in) :: weights(:) !! Nonnegative variable weights; only variables eligible for the denominator
        !! contribute their weight.
        real(dp) :: similarity
        real(dp) :: denominator
        real(dp) :: numerator
        real(dp) :: score
        integer :: k
        logical :: x_true
        logical :: y_true

        if (size(y) /= size(x) .or. size(types) /= size(x) .or. size(weights) /= size(x)) then
            similarity = proxy_nan()
            return
        end if
        numerator = 0.0_dp
        denominator = 0.0_dp
        do k = 1, size(x)
            if (proxy_is_missing(x(k)) .or. proxy_is_missing(y(k))) cycle
            select case (types(k))
            case (proxy_gower_logical)
                x_true = abs(x(k)) > tiny(1.0_dp)
                y_true = abs(y(k)) > tiny(1.0_dp)
                if (.not. (x_true .or. y_true)) cycle
                denominator = denominator + weights(k)
                if (x_true .and. y_true) numerator = numerator + weights(k)
            case (proxy_gower_factor)
                denominator = denominator + weights(k)
                if (x(k) <= y(k) .and. x(k) >= y(k)) numerator = numerator + weights(k)
            case (proxy_gower_metric, proxy_gower_ordinal)
                score = 1.0_dp - abs(x(k) - y(k))
                denominator = denominator + weights(k)
                numerator = numerator + weights(k) * score
            case default
                similarity = proxy_nan()
                return
            end select
        end do
        if (abs(denominator) <= tiny(1.0_dp)) then
            similarity = proxy_nan()
        else
            similarity = numerator / denominator
        end if
    end function gower_pair_similarity

    subroutine prepare_gower_single(x, types, z, weights_out, weights, ranges, minima)
        real(dp), intent(in) :: x(:, :) !! Raw observation matrix to prepare for an auto-Gower calculation.
        integer, intent(in) :: types(:) !! Type code for every variable/column in `x`.
        real(dp), allocatable, intent(out) :: z(:, :) !! Allocated transformed matrix with ordinal conversion and metric scaling
        !! applied.
        real(dp), allocatable, intent(out) :: weights_out(:) !! Allocated effective variable weights after recycling optional
        !! input weights.
        real(dp), intent(in), optional :: weights(:) !! Optional variable weights recycled across columns; defaults to one.
        real(dp), intent(in), optional :: ranges(:) !! Optional scaling ranges recycled across metric/ordinal columns.
        real(dp), intent(in), optional :: minima(:) !! Optional scaling minima recycled across metric/ordinal columns.
        integer :: k
        integer :: m_index
        real(dp) :: col_min
        real(dp) :: col_max
        real(dp) :: scale

        allocate(z(size(x, 1), size(x, 2)), weights_out(size(x, 2)))
        z = x
        call recycle_weights(weights, weights_out)
        call transform_ordinals(z, types)
        m_index = 0
        do k = 1, size(x, 2)
            if (types(k) /= proxy_gower_metric .and. types(k) /= proxy_gower_ordinal) cycle
            m_index = m_index + 1
            call finite_minmax(z(:, k), col_min, col_max)
            if (present(minima)) col_min = minima(1 + modulo(m_index - 1, size(minima)))
            scale = col_max - col_min
            if (present(ranges)) scale = ranges(1 + modulo(m_index - 1, size(ranges)))
            if (abs(scale) <= tiny(1.0_dp)) scale = 1.0_dp
            call scale_column(z(:, k), col_min, scale)
        end do
    end subroutine prepare_gower_single

    subroutine prepare_gower_pair(x, y, types, zx, zy, weights_out, weights, ranges_x, minima_x, ranges_y, minima_y)
        real(dp), intent(in) :: x(:, :) !! Raw first observation matrix for a cross-Gower calculation.
        real(dp), intent(in) :: y(:, :) !! Raw second observation matrix with the same columns and encodings as `x`.
        integer, intent(in) :: types(:) !! Shared variable type code for every column.
        real(dp), allocatable, intent(out) :: zx(:, :) !! Allocated transformed/scaled copy of `x`.
        real(dp), allocatable, intent(out) :: zy(:, :) !! Allocated transformed/scaled copy of `y`.
        real(dp), allocatable, intent(out) :: weights_out(:) !! Effective recycled variable weights.
        real(dp), intent(in), optional :: weights(:) !! Optional variable weights recycled across all columns.
        real(dp), intent(in), optional :: ranges_x(:) !! Optional metric/ordinal scaling ranges applied to `x`.
        real(dp), intent(in), optional :: minima_x(:) !! Optional metric/ordinal scaling minima applied to `x`.
        real(dp), intent(in), optional :: ranges_y(:) !! Optional metric/ordinal scaling ranges applied to `y`; defaults to `x`
        !! scaling.
        real(dp), intent(in), optional :: minima_y(:) !! Optional metric/ordinal scaling minima applied to `y`; defaults to `x`
        !! scaling.
        integer :: k
        integer :: m_index
        real(dp) :: max_x
        real(dp) :: max_y
        real(dp) :: min_x
        real(dp) :: min_y
        real(dp) :: rx
        real(dp) :: ry

        if (size(x, 2) /= size(y, 2) .or. size(types) /= size(x, 2)) then
            allocate(zx(0, 0), zy(0, 0), weights_out(0))
            return
        end if
        allocate(zx(size(x, 1), size(x, 2)), zy(size(y, 1), size(y, 2)), weights_out(size(x, 2)))
        zx = x
        zy = y
        call recycle_weights(weights, weights_out)
        call transform_ordinals(zx, types)
        call transform_ordinals(zy, types)
        m_index = 0
        do k = 1, size(x, 2)
            if (types(k) /= proxy_gower_metric .and. types(k) /= proxy_gower_ordinal) cycle
            m_index = m_index + 1
            call finite_minmax(zx(:, k), min_x, max_x)
            call finite_minmax(zy(:, k), min_y, max_y)
            min_x = min(min_x, min_y)
            max_x = max(max_x, max_y)
            min_y = min_x
            max_y = max_x
            if (present(minima_x)) min_x = minima_x(1 + modulo(m_index - 1, size(minima_x)))
            if (present(minima_y)) min_y = minima_y(1 + modulo(m_index - 1, size(minima_y)))
            rx = max_x - min_x
            ry = max_y - min_y
            if (present(ranges_x)) rx = ranges_x(1 + modulo(m_index - 1, size(ranges_x)))
            if (present(ranges_y)) ry = ranges_y(1 + modulo(m_index - 1, size(ranges_y)))
            if (.not. present(ranges_y)) ry = rx
            if (.not. present(minima_y)) min_y = min_x
            if (abs(rx) <= tiny(1.0_dp)) rx = 1.0_dp
            if (abs(ry) <= tiny(1.0_dp)) ry = 1.0_dp
            call scale_column(zx(:, k), min_x, rx)
            call scale_column(zy(:, k), min_y, ry)
        end do
    end subroutine prepare_gower_pair

    subroutine transform_ordinals(z, types)
        real(dp), intent(inout) :: z(:, :) !! Matrix whose ordinal-coded columns are transformed in place to
        !! `(code-1)/(max(code)-1)`.
        integer, intent(in) :: types(:) !! Type code for each column; only `proxy_gower_ordinal` columns are modified.
        integer :: i
        integer :: k
        real(dp) :: col_min
        real(dp) :: col_max
        real(dp) :: den

        do k = 1, min(size(z, 2), size(types))
            if (types(k) /= proxy_gower_ordinal) cycle
            call finite_minmax(z(:, k), col_min, col_max)
            den = col_max - 1.0_dp
            if (abs(den) <= tiny(1.0_dp)) den = 1.0_dp
            do i = 1, size(z, 1)
                if (.not. proxy_is_missing(z(i, k))) z(i, k) = (z(i, k) - 1.0_dp) / den
            end do
        end do
    end subroutine transform_ordinals

    pure subroutine finite_minmax(values, minimum, maximum)
        real(dp), intent(in) :: values(:) !! Numeric column containing optional NaN missing values.
        real(dp), intent(out) :: minimum !! Smallest nonmissing value, or zero when every value is missing.
        real(dp), intent(out) :: maximum !! Largest nonmissing value, or zero when every value is missing.
        integer :: i
        logical :: found

        minimum = 0.0_dp
        maximum = 0.0_dp
        found = .false.
        do i = 1, size(values)
            if (proxy_is_missing(values(i))) cycle
            if (.not. found) then
                minimum = values(i)
                maximum = values(i)
                found = .true.
            else
                minimum = min(minimum, values(i))
                maximum = max(maximum, values(i))
            end if
        end do
    end subroutine finite_minmax

    pure subroutine scale_column(values, minimum, range)
        real(dp), intent(inout) :: values(:) !! Metric/ordinal column to be scaled in place; NaN entries remain missing.
        real(dp), intent(in) :: minimum !! Minimum to subtract from every observed component.
        real(dp), intent(in) :: range !! Positive or negative divisor used for scaling; caller ensures it is nonzero.
        integer :: i

        do i = 1, size(values)
            if (.not. proxy_is_missing(values(i))) values(i) = (values(i) - minimum) / range
        end do
    end subroutine scale_column

    pure subroutine recycle_weights(weights, output)
        real(dp), intent(in), optional :: weights(:) !! Optional source weights; values are recycled modulo their length when
        !! fewer than `output` entries are provided.
        real(dp), intent(out) :: output(:) !! Effective weight vector, set to ones when `weights` is absent or empty.
        integer :: i

        output = 1.0_dp
        if (.not. present(weights)) return
        if (size(weights) == 0) return
        do i = 1, size(output)
            output(i) = weights(1 + modulo(i - 1, size(weights)))
        end do
    end subroutine recycle_weights

end module proxy_gower
