module proxy_numeric_measures
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_positive_inf, ieee_negative_inf
    use proxy_kinds, only: dp
    use proxy_ieee, only: proxy_nan, proxy_is_missing
    implicit none
    private

    public :: euclidean_distance, mahalanobis_distance, bhjattacharyya_distance
    public :: manhattan_distance, supremum_distance, minkowski_distance, canberra_distance
    public :: wave_hedges_distance, divergence_distance, kullback_leibler_distance
    public :: bray_curtis_distance, soergel_distance, podani_distance
    public :: chord_distance, geodesic_distance, whittaker_distance, hellinger_distance
    public :: fuzzy_jaccard_distance, cosine_similarity, angular_similarity
    public :: extended_jaccard_similarity, extended_dice_similarity, correlation_similarity
    public :: mutual_information_similarity

contains

    pure function euclidean_distance(x, y) result(distance)
        real(dp), intent(in) :: x(:) !! First numeric observation vector; NaN components are omitted pairwise and compensated as
        !! in proxy's C implementation.
        real(dp), intent(in) :: y(:) !! Second numeric observation vector; must have the same number of components as `x`.
        real(dp) :: distance
        real(dp) :: delta
        integer :: count
        integer :: i
        integer :: n

        n = size(x)
        if (size(y) /= n) then
            distance = proxy_nan()
            return
        end if
        distance = 0.0_dp
        count = 0
        do i = 1, n
            if (.not. proxy_is_missing(x(i)) .and. .not. proxy_is_missing(y(i))) then
                delta = x(i) - y(i)
                if (.not. proxy_is_missing(delta)) then
                    distance = distance + delta * delta
                    count = count + 1
                end if
            end if
        end do
        if (count == 0) then
            distance = proxy_nan()
            return
        end if
        if (count /= n) distance = distance * real(n, dp) / real(count, dp)
        distance = sqrt(distance)
    end function euclidean_distance

    function mahalanobis_distance(x, y, covariance) result(distance)
        real(dp), intent(in) :: x(:) !! First numeric observation vector; all components participate in the quadratic form.
        real(dp), intent(in) :: y(:) !! Second numeric observation vector; must have the same length as `x`.
        real(dp), intent(in) :: covariance(:, :) !! Symmetric covariance matrix with shape `(size(x),size(x))`; it must be
        !! nonsingular.
        real(dp) :: distance
        real(dp), allocatable :: a(:, :)
        real(dp), allocatable :: b(:)
        integer :: info

        if (size(y) /= size(x) .or. size(covariance, 1) /= size(x) .or. size(covariance, 2) /= size(x)) then
            distance = proxy_nan()
            return
        end if
        allocate(a(size(x), size(x)), b(size(x)))
        a = covariance
        b = x - y
        call solve_linear_system(a, b, info)
        if (info /= 0) then
            distance = proxy_nan()
        else
            distance = sqrt(max(0.0_dp, dot_product(x - y, b)))
        end if
    end function mahalanobis_distance

    pure function bhjattacharyya_distance(x, y) result(distance)
        real(dp), intent(in) :: x(:) !! First nonnegative abundance or probability vector; negative values yield NaN through the
        !! square root.
        real(dp), intent(in) :: y(:) !! Second nonnegative vector with the same length as `x`.
        real(dp) :: distance

        if (size(y) /= size(x)) then
            distance = proxy_nan()
            return
        end if
        if (any(x < 0.0_dp) .or. any(y < 0.0_dp)) then
            distance = proxy_nan()
            return
        end if
        distance = sqrt(sum((sqrt(x) - sqrt(y))**2))
    end function bhjattacharyya_distance

    pure function manhattan_distance(x, y) result(distance)
        real(dp), intent(in) :: x(:) !! First numeric observation vector; missing components are omitted and the sum is rescaled
        !! for their exclusion.
        real(dp), intent(in) :: y(:) !! Second numeric observation vector; must have the same number of components as `x`.
        real(dp) :: distance
        real(dp) :: delta
        integer :: count
        integer :: i
        integer :: n

        n = size(x)
        if (size(y) /= n) then
            distance = proxy_nan()
            return
        end if
        distance = 0.0_dp
        count = 0
        do i = 1, n
            if (.not. proxy_is_missing(x(i)) .and. .not. proxy_is_missing(y(i))) then
                delta = abs(x(i) - y(i))
                if (.not. proxy_is_missing(delta)) then
                    distance = distance + delta
                    count = count + 1
                end if
            end if
        end do
        if (count == 0) then
            distance = proxy_nan()
            return
        end if
        if (count /= n) distance = distance * real(n, dp) / real(count, dp)
    end function manhattan_distance

    pure function supremum_distance(x, y) result(distance)
        real(dp), intent(in) :: x(:) !! First numeric observation vector; NaN components are omitted pairwise without compensation.
        real(dp), intent(in) :: y(:) !! Second numeric observation vector; must have the same number of components as `x`.
        real(dp) :: distance
        real(dp) :: delta
        integer :: count
        integer :: i

        if (size(y) /= size(x)) then
            distance = proxy_nan()
            return
        end if
        distance = -huge(1.0_dp)
        count = 0
        do i = 1, size(x)
            if (.not. proxy_is_missing(x(i)) .and. .not. proxy_is_missing(y(i))) then
                delta = abs(x(i) - y(i))
                if (.not. proxy_is_missing(delta)) then
                    distance = max(distance, delta)
                    count = count + 1
                end if
            end if
        end do
        if (count == 0) distance = proxy_nan()
    end function supremum_distance

    pure function minkowski_distance(x, y, p) result(distance)
        real(dp), intent(in) :: x(:) !! First numeric observation vector; missing components are omitted and compensated exactly
        !! as in proxy's C fast path.
        real(dp), intent(in) :: y(:) !! Second numeric observation vector; must have the same number of components as `x`.
        real(dp), intent(in) :: p !! Positive Minkowski exponent; exact positive infinity retains the upstream C fast-path
        !! result of one when any component pair is valid.
        real(dp) :: distance
        real(dp) :: delta
        integer :: count
        integer :: i
        integer :: n

        n = size(x)
        if (size(y) /= n .or. p <= 0.0_dp) then
            distance = proxy_nan()
            return
        end if
        if (.not. ieee_is_finite(p)) then
            count = 0
            do i = 1, n
                if (.not. proxy_is_missing(x(i)) .and. .not. proxy_is_missing(y(i))) then
                    delta = abs(x(i) - y(i))
                    if (.not. proxy_is_missing(delta)) count = count + 1
                end if
            end do
            if (count == 0) then
                distance = proxy_nan()
            else
                distance = 1.0_dp
            end if
            return
        end if
        distance = 0.0_dp
        count = 0
        do i = 1, n
            if (.not. proxy_is_missing(x(i)) .and. .not. proxy_is_missing(y(i))) then
                delta = abs(x(i) - y(i))
                if (.not. proxy_is_missing(delta)) then
                    distance = distance + delta**p
                    count = count + 1
                end if
            end if
        end do
        if (count == 0) then
            distance = proxy_nan()
            return
        end if
        if (count /= n) distance = distance * real(n, dp) / real(count, dp)
        distance = distance**(1.0_dp / p)
    end function minkowski_distance

    pure function canberra_distance(x, y) result(distance)
        real(dp), intent(in) :: x(:) !! First numeric observation vector; undefined zero-over-zero components are excluded and
        !! remaining terms are compensated.
        real(dp), intent(in) :: y(:) !! Second numeric observation vector; must have the same number of components as `x`.
        real(dp) :: distance
        real(dp) :: den
        real(dp) :: diff
        real(dp) :: term
        integer :: count
        integer :: i
        integer :: n

        n = size(x)
        if (size(y) /= n) then
            distance = proxy_nan()
            return
        end if
        distance = 0.0_dp
        count = 0
        do i = 1, n
            if (.not. proxy_is_missing(x(i)) .and. .not. proxy_is_missing(y(i))) then
                den = abs(x(i) + y(i))
                diff = abs(x(i) - y(i))
                if (den > tiny(1.0_dp) .or. diff > tiny(1.0_dp)) then
                    if (.not. ieee_is_finite(diff) .and. .not. ieee_is_finite(den) .and. diff <= den .and. diff >= den) then
                        term = 1.0_dp
                    else if (den <= tiny(1.0_dp)) then
                        term = ieee_value(0.0_dp, ieee_positive_inf)
                    else
                        term = diff / den
                    end if
                    if (.not. proxy_is_missing(term)) then
                        distance = distance + term
                        count = count + 1
                    end if
                end if
            end if
        end do
        if (count == 0) then
            distance = proxy_nan()
            return
        end if
        if (count /= n) distance = distance * real(n, dp) / real(count, dp)
    end function canberra_distance

    pure function wave_hedges_distance(x, y) result(distance)
        real(dp), intent(in) :: x(:) !! First numeric vector supplied to proxy's literal Wave/Hedges R formula.
        real(dp), intent(in) :: y(:) !! Second numeric vector; R `min`/`max` operate globally across both vectors.
        real(dp) :: distance
        real(dp) :: maximum
        real(dp) :: minimum

        if (size(y) /= size(x) .or. any(proxy_is_missing(x)) .or. any(proxy_is_missing(y))) then
            distance = proxy_nan()
            return
        end if
        minimum = min(minval(x), minval(y))
        maximum = max(maxval(x), maxval(y))
        if (abs(maximum) <= tiny(1.0_dp)) then
            if (abs(minimum) <= tiny(1.0_dp)) then
                distance = proxy_nan()
            else if (minimum > 0.0_dp) then
                distance = ieee_value(0.0_dp, ieee_negative_inf)
            else
                distance = ieee_value(0.0_dp, ieee_positive_inf)
            end if
        else
            distance = 1.0_dp - minimum / maximum
        end if
    end function wave_hedges_distance

    pure function divergence_distance(x, y) result(distance)
        real(dp), intent(in) :: x(:) !! First numeric vector in the divergence formula; NaN terms are omitted like R `!is.nan`.
        real(dp), intent(in) :: y(:) !! Second numeric vector aligned with `x`; nonzero-over-zero terms contribute positive
        !! infinity.
        real(dp) :: distance
        real(dp) :: den
        real(dp) :: num
        real(dp) :: term
        integer :: i

        if (size(y) /= size(x)) then
            distance = proxy_nan()
            return
        end if
        distance = 0.0_dp
        do i = 1, size(x)
            if (proxy_is_missing(x(i)) .or. proxy_is_missing(y(i))) cycle
            num = (x(i) - y(i))**2
            den = (x(i) + y(i))**2
            if (den <= tiny(1.0_dp)) then
                if (num <= tiny(1.0_dp)) cycle
                term = ieee_value(0.0_dp, ieee_positive_inf)
            else
                term = num / den
            end if
            if (.not. proxy_is_missing(term)) distance = distance + term
        end do
    end function divergence_distance

    pure function kullback_leibler_distance(x, y) result(distance)
        real(dp), intent(in) :: x(:) !! First vector normalized by its sum before proxy's literal Kullback-Leibler expression.
        real(dp), intent(in) :: y(:) !! Second vector normalized by its sum; zero/invalid components retain R NaN/Inf behavior.
        real(dp) :: distance
        real(dp) :: px
        real(dp) :: qy
        real(dp) :: sx
        real(dp) :: sy
        integer :: i

        if (size(y) /= size(x) .or. any(proxy_is_missing(x)) .or. any(proxy_is_missing(y))) then
            distance = proxy_nan()
            return
        end if
        sx = sum(x)
        sy = sum(y)
        if (abs(sx) <= tiny(1.0_dp) .or. abs(sy) <= tiny(1.0_dp)) then
            distance = proxy_nan()
            return
        end if
        distance = 0.0_dp
        do i = 1, size(x)
            px = x(i) / sx
            qy = y(i) / sy
            if (abs(px) <= tiny(1.0_dp)) then
                distance = proxy_nan()
                return
            end if
            if (abs(qy) <= tiny(1.0_dp)) then
                if (px > 0.0_dp) then
                    distance = ieee_value(0.0_dp, ieee_positive_inf)
                else
                    distance = proxy_nan()
                end if
                return
            end if
            if (px / qy <= 0.0_dp) then
                distance = proxy_nan()
                return
            end if
            distance = distance + px * log(px / qy)
        end do
    end function kullback_leibler_distance

    pure function bray_curtis_distance(x, y) result(distance)
        real(dp), intent(in) :: x(:) !! First abundance vector in the Bray-Curtis dissimilarity.
        real(dp), intent(in) :: y(:) !! Second abundance vector with the same length as `x`.
        real(dp) :: distance
        real(dp) :: den

        if (size(y) /= size(x)) then
            distance = proxy_nan()
            return
        end if
        den = sum(x + y)
        if (abs(den) <= tiny(1.0_dp)) then
            if (sum(abs(x - y)) <= tiny(1.0_dp)) then
                distance = proxy_nan()
            else if (den < 0.0_dp) then
                distance = ieee_value(0.0_dp, ieee_negative_inf)
            else
                distance = ieee_value(0.0_dp, ieee_positive_inf)
            end if
        else
            distance = sum(abs(x - y)) / den
        end if
    end function bray_curtis_distance

    pure function soergel_distance(x, y) result(distance)
        real(dp), intent(in) :: x(:) !! First numeric vector supplied to proxy's literal Soergel R expression.
        real(dp), intent(in) :: y(:) !! Second vector; its global maximum is combined with `x` by R `max`, not componentwise `pmax`.
        real(dp) :: distance
        real(dp) :: den
        real(dp) :: num

        if (size(y) /= size(x) .or. any(proxy_is_missing(x)) .or. any(proxy_is_missing(y))) then
            distance = proxy_nan()
            return
        end if
        den = max(maxval(x), maxval(y))
        num = sum(abs(x - y))
        if (abs(den) <= tiny(1.0_dp)) then
            if (num <= tiny(1.0_dp)) then
                distance = proxy_nan()
            else if (den < 0.0_dp) then
                distance = ieee_value(0.0_dp, ieee_negative_inf)
            else
                distance = ieee_value(0.0_dp, ieee_positive_inf)
            end if
        else
            distance = num / den
        end if
    end function soergel_distance

    pure function podani_distance(x, y) result(distance)
        real(dp), intent(in) :: x(:) !! First rank vector; exact ties and zero scores have the special Podani discordance meaning
        !! used by proxy.
        real(dp), intent(in) :: y(:) !! Second rank vector with the same length and ranking scale as `x`.
        real(dp) :: distance
        integer :: a
        integer :: b
        integer :: c
        integer :: d
        integer :: i
        integer :: j
        integer :: n
        integer :: z
        logical :: xeq
        logical :: yeq

        n = size(x)
        if (size(y) /= n .or. n < 2) then
            distance = proxy_nan()
            return
        end if
        a = 0
        b = 0
        c = 0
        d = 0
        do i = 1, n - 1
            do j = i + 1, n
                if ((x(i) < x(j) .and. y(i) < y(j)) .or. (x(i) > x(j) .and. y(i) > y(j))) a = a + 1
                if ((x(i) < x(j) .and. y(i) > y(j)) .or. (x(i) > x(j) .and. y(i) < y(j))) b = b + 1
                xeq = x(i) <= x(j) .and. x(i) >= x(j)
                yeq = y(i) <= y(j) .and. y(i) >= y(j)
                if (xeq .and. yeq) then
                    if ((abs(x(i)) <= tiny(1.0_dp) .and. abs(y(i)) <= tiny(1.0_dp)) .or. &
                        (x(i) > 0.0_dp .and. y(i) > 0.0_dp)) c = c + 1
                end if
                z = 0
                if (abs(x(i)) <= tiny(1.0_dp)) z = z + 1
                if (abs(x(j)) <= tiny(1.0_dp)) z = z + 1
                if (abs(y(i)) <= tiny(1.0_dp)) z = z + 1
                if (abs(y(j)) <= tiny(1.0_dp)) z = z + 1
                if ((xeq .or. yeq) .and. z > 0 .and. z < 4) d = d + 1
            end do
        end do
        distance = 1.0_dp - 2.0_dp * real(a - b + c - d, dp) / real(n * (n - 1), dp)
    end function podani_distance

    pure function chord_distance(x, y) result(distance)
        real(dp), intent(in) :: x(:) !! First numeric vector; the distance is based on its cosine with `y`.
        real(dp), intent(in) :: y(:) !! Second numeric vector with the same length as `x`.
        real(dp) :: distance
        real(dp) :: c

        c = raw_cosine(x, y)
        if (proxy_is_missing(c) .or. c > 1.0_dp) then
            distance = proxy_nan()
        else
            distance = sqrt(2.0_dp * (1.0_dp - c))
        end if
    end function chord_distance

    pure function geodesic_distance(x, y) result(distance)
        real(dp), intent(in) :: x(:) !! First numeric vector; the result is the angle in radians between this vector and `y`.
        real(dp), intent(in) :: y(:) !! Second numeric vector with the same length as `x`.
        real(dp) :: distance
        real(dp) :: c

        c = raw_cosine(x, y)
        if (proxy_is_missing(c) .or. c < -1.0_dp .or. c > 1.0_dp) then
            distance = proxy_nan()
        else
            distance = acos(c)
        end if
    end function geodesic_distance

    pure function whittaker_distance(x, y) result(distance)
        real(dp), intent(in) :: x(:) !! First nonnegative abundance vector; it is normalized by its total abundance.
        real(dp), intent(in) :: y(:) !! Second nonnegative abundance vector with the same length as `x`.
        real(dp) :: distance
        real(dp) :: sx
        real(dp) :: sy

        if (size(y) /= size(x)) then
            distance = proxy_nan()
            return
        end if
        sx = sum(x)
        sy = sum(y)
        if (abs(sx) <= tiny(1.0_dp) .or. abs(sy) <= tiny(1.0_dp)) then
            distance = proxy_nan()
        else
            distance = sum(abs(x / sx - y / sy)) / 2.0_dp
        end if
    end function whittaker_distance

    pure function hellinger_distance(x, y) result(distance)
        real(dp), intent(in) :: x(:) !! First nonnegative abundance vector; values are normalized to a probability vector.
        real(dp), intent(in) :: y(:) !! Second nonnegative abundance vector with the same length as `x`.
        real(dp) :: distance
        real(dp) :: sx
        real(dp) :: sy

        if (size(y) /= size(x)) then
            distance = proxy_nan()
            return
        end if
        if (any(x < 0.0_dp) .or. any(y < 0.0_dp)) then
            distance = proxy_nan()
            return
        end if
        sx = sum(x)
        sy = sum(y)
        if (sx <= tiny(1.0_dp) .or. sy <= tiny(1.0_dp)) then
            distance = proxy_nan()
        else
            distance = sqrt(sum((sqrt(x / sx) - sqrt(y / sy))**2))
        end if
    end function hellinger_distance

    pure function fuzzy_jaccard_distance(x, y) result(distance)
        real(dp), intent(in) :: x(:) !! First fuzzy-membership vector, normally restricted to `[0,1]`; nonfinite pairs are omitted.
        real(dp), intent(in) :: y(:) !! Second fuzzy-membership vector with the same length as `x`.
        real(dp) :: distance
        real(dp) :: smax
        real(dp) :: smin
        integer :: count
        integer :: i

        if (size(y) /= size(x)) then
            distance = proxy_nan()
            return
        end if
        smax = 0.0_dp
        smin = 0.0_dp
        count = 0
        do i = 1, size(x)
            if (ieee_is_finite(x(i)) .and. ieee_is_finite(y(i))) then
                smax = smax + max(x(i), y(i))
                smin = smin + min(x(i), y(i))
                count = count + 1
            end if
        end do
        if (count == 0 .or. .not. ieee_is_finite(smin)) then
            distance = proxy_nan()
        else if (abs(smax) <= tiny(1.0_dp)) then
            distance = 0.0_dp
        else
            distance = 1.0_dp - smin / smax
        end if
    end function fuzzy_jaccard_distance

    pure function cosine_similarity(x, y) result(similarity)
        real(dp), intent(in) :: x(:) !! First numeric vector; NaN component pairs are omitted as in proxy's optimized cosine
        !! implementation.
        real(dp), intent(in) :: y(:) !! Second numeric vector with the same length as `x`.
        real(dp) :: similarity
        real(dp) :: xx
        real(dp) :: xy
        real(dp) :: yy
        integer :: count
        integer :: i

        if (size(y) /= size(x)) then
            similarity = proxy_nan()
            return
        end if
        xx = 0.0_dp
        xy = 0.0_dp
        yy = 0.0_dp
        count = 0
        do i = 1, size(x)
            if (.not. proxy_is_missing(x(i)) .and. .not. proxy_is_missing(y(i))) then
                if (.not. proxy_is_missing(x(i) * y(i))) then
                    xy = xy + x(i) * y(i)
                    xx = xx + x(i) * x(i)
                    yy = yy + y(i) * y(i)
                    count = count + 1
                end if
            end if
        end do
        if (count == 0 .or. .not. ieee_is_finite(xy)) then
            similarity = proxy_nan()
        else if (xx <= tiny(1.0_dp) .and. yy <= tiny(1.0_dp)) then
            similarity = 1.0_dp
        else if (xx <= tiny(1.0_dp) .or. yy <= tiny(1.0_dp)) then
            similarity = 0.0_dp
        else
            similarity = xy / sqrt(xx * yy)
        end if
    end function cosine_similarity

    pure function angular_similarity(x, y) result(similarity)
        real(dp), intent(in) :: x(:) !! First numeric vector used to compute the angular similarity.
        real(dp), intent(in) :: y(:) !! Second numeric vector with the same length as `x`.
        real(dp) :: similarity
        real(dp) :: c
        real(dp), parameter :: pi = acos(-1.0_dp)

        c = raw_cosine(x, y)
        if (proxy_is_missing(c) .or. c < -1.0_dp .or. c > 1.0_dp) then
            similarity = proxy_nan()
        else
            similarity = 1.0_dp - acos(c) / pi
        end if
    end function angular_similarity

    pure function extended_jaccard_similarity(x, y) result(similarity)
        real(dp), intent(in) :: x(:) !! First numeric vector in the extended Jaccard/Tanimoto coefficient; NaN pairs are omitted.
        real(dp), intent(in) :: y(:) !! Second numeric vector with the same length as `x`.
        real(dp) :: similarity

        similarity = extended_binary_similarity(x, y, 1.0_dp)
    end function extended_jaccard_similarity

    pure function extended_dice_similarity(x, y) result(similarity)
        real(dp), intent(in) :: x(:) !! First numeric vector in the extended Dice/Sorensen coefficient; NaN pairs are omitted.
        real(dp), intent(in) :: y(:) !! Second numeric vector with the same length as `x`.
        real(dp) :: similarity

        similarity = extended_binary_similarity(x, y, 2.0_dp)
    end function extended_dice_similarity

    pure function correlation_similarity(x, y) result(similarity)
        real(dp), intent(in) :: x(:) !! First numeric vector; proxy's R correlation measure uses ordinary centered values and
        !! propagates missing data.
        real(dp), intent(in) :: y(:) !! Second numeric vector with the same length as `x`.
        real(dp) :: similarity
        real(dp) :: mx
        real(dp) :: my
        real(dp) :: xx
        real(dp) :: xy
        real(dp) :: yy

        if (size(y) /= size(x) .or. size(x) == 0) then
            similarity = proxy_nan()
            return
        end if
        mx = sum(x) / real(size(x), dp)
        my = sum(y) / real(size(y), dp)
        xy = dot_product(x - mx, y - my)
        xx = dot_product(x - mx, x - mx)
        yy = dot_product(y - my, y - my)
        if (xx <= tiny(1.0_dp) .or. yy <= tiny(1.0_dp)) then
            similarity = proxy_nan()
        else
            similarity = xy / sqrt(xx * yy)
        end if
    end function correlation_similarity

    pure function mutual_information_similarity(x, y) result(similarity)
        real(dp), intent(in) :: x(:) !! First vector interpreted as Boolean by zero/nonzero status; NaN pairs are omitted.
        real(dp), intent(in) :: y(:) !! Second Boolean-coded vector with the same length as `x`.
        real(dp) :: similarity
        integer :: a
        integer :: b
        integer :: c
        integer :: d
        integer :: n

        call binary_counts_real(x, y, a, b, c, d, n)
        if (n == 0) then
            similarity = proxy_nan()
        else
            similarity = mutual_from_counts(a, b, c, d, n)
            if (n /= size(x)) similarity = similarity * real(size(x), dp) / real(n, dp)
        end if
    end function mutual_information_similarity

    pure function raw_cosine(x, y) result(value)
        real(dp), intent(in) :: x(:) !! First vector for an uncorrected dot-product cosine used by R-level formulas.
        real(dp), intent(in) :: y(:) !! Second vector with the same length as `x`.
        real(dp) :: value
        real(dp) :: den

        if (size(y) /= size(x)) then
            value = proxy_nan()
            return
        end if
        den = sqrt(dot_product(x, x) * dot_product(y, y))
        if (den <= tiny(1.0_dp)) then
            value = proxy_nan()
        else
            value = dot_product(x, y) / den
        end if
    end function raw_cosine

    pure function extended_binary_similarity(x, y, factor) result(similarity)
        real(dp), intent(in) :: x(:) !! First numeric vector in an extended set similarity; NaN component pairs are omitted.
        real(dp), intent(in) :: y(:) !! Second numeric vector with the same length as `x`.
        real(dp), intent(in) :: factor !! Denominator convention: `1` gives extended Jaccard and `2` gives extended Dice.
        real(dp) :: similarity
        real(dp) :: den
        real(dp) :: dist
        real(dp) :: prod
        real(dp) :: xy
        integer :: count
        integer :: i

        if (size(y) /= size(x)) then
            similarity = proxy_nan()
            return
        end if
        dist = 0.0_dp
        xy = 0.0_dp
        count = 0
        do i = 1, size(x)
            if (.not. proxy_is_missing(x(i)) .and. .not. proxy_is_missing(y(i))) then
                dist = dist + (x(i) - y(i))**2
                prod = x(i) * y(i)
                if (.not. proxy_is_missing(prod)) then
                    xy = xy + prod
                    count = count + 1
                end if
            end if
        end do
        if (count == 0 .or. .not. ieee_is_finite(xy)) then
            similarity = proxy_nan()
            return
        end if
        den = dist / factor + xy
        if (abs(den) <= tiny(1.0_dp)) then
            similarity = 1.0_dp
        else
            similarity = xy / den
        end if
    end function extended_binary_similarity

    pure subroutine binary_counts_real(x, y, a, b, c, d, n)
        real(dp), intent(in) :: x(:) !! First numeric vector interpreted as FALSE for zero and TRUE for nonzero; NaN values are
        !! missing.
        real(dp), intent(in) :: y(:) !! Second numeric Boolean-coded vector with the same length as `x`.
        integer, intent(out) :: a !! Number of valid component pairs where both values are nonzero/TRUE.
        integer, intent(out) :: b !! Number of valid pairs where `x` is nonzero/TRUE and `y` is zero/FALSE.
        integer, intent(out) :: c !! Number of valid pairs where `x` is zero/FALSE and `y` is nonzero/TRUE.
        integer, intent(out) :: d !! Number of valid component pairs where both values are zero/FALSE.
        integer, intent(out) :: n !! Number of valid nonmissing component pairs contributing to the contingency table.
        integer :: i
        logical :: tx
        logical :: ty

        a = 0
        b = 0
        c = 0
        d = 0
        n = 0
        do i = 1, min(size(x), size(y))
            if (proxy_is_missing(x(i)) .or. proxy_is_missing(y(i))) cycle
            tx = abs(x(i)) > tiny(1.0_dp)
            ty = abs(y(i)) > tiny(1.0_dp)
            if (tx .and. ty) then
                a = a + 1
            else if (tx) then
                b = b + 1
            else if (ty) then
                c = c + 1
            else
                d = d + 1
            end if
            n = n + 1
        end do
    end subroutine binary_counts_real

    pure function mutual_from_counts(a, b, c, d, n) result(value)
        integer, intent(in) :: a !! Count of TRUE/TRUE observations in the 2-by-2 table.
        integer, intent(in) :: b !! Count of TRUE/FALSE observations in the 2-by-2 table.
        integer, intent(in) :: c !! Count of FALSE/TRUE observations in the 2-by-2 table.
        integer, intent(in) :: d !! Count of FALSE/FALSE observations in the 2-by-2 table.
        integer, intent(in) :: n !! Total number of valid observations in the table.
        real(dp) :: value
        integer :: cx
        integer :: cy

        cx = a + b
        cy = a + c
        if (n <= 0 .or. cx == 0 .or. cy == 0 .or. cx == n .or. cy == n) then
            value = 0.0_dp
            return
        end if
        value = 0.0_dp
        if (a > 0) value = value + mi_term(a, cx, cy, n)
        if (b > 0) value = value + mi_term(b, cx, n - cy, n)
        if (c > 0) value = value + mi_term(c, n - cx, cy, n)
        if (d > 0) value = value + mi_term(d, n - cx, n - cy, n)
    end function mutual_from_counts

    pure function mi_term(cell, row_total, col_total, n) result(value)
        integer, intent(in) :: cell !! Positive cell count contributing one mutual-information term.
        integer, intent(in) :: row_total !! Marginal count for the cell's row.
        integer, intent(in) :: col_total !! Marginal count for the cell's column.
        integer, intent(in) :: n !! Total number of valid observations.
        real(dp) :: value

        value = real(cell, dp) / real(n, dp) * &
                log(real(cell * n, dp) / real(row_total * col_total, dp))
    end function mi_term

    subroutine solve_linear_system(a, b, info)
        real(dp), intent(inout) :: a(:, :) !! Square coefficient matrix; overwritten by Gaussian elimination factors.
        real(dp), intent(inout) :: b(:) !! Right-hand side on entry and solution vector on successful return.
        integer, intent(out) :: info !! Zero on success; nonzero when dimensions are invalid or a numerically singular pivot is
        !! found.
        integer :: i
        integer :: j
        integer :: k
        integer :: n
        integer :: pivot
        real(dp) :: factor
        real(dp) :: max_abs
        real(dp) :: temp

        n = size(b)
        if (size(a, 1) /= n .or. size(a, 2) /= n) then
            info = 1
            return
        end if
        info = 0
        do k = 1, n - 1
            pivot = k
            max_abs = abs(a(k, k))
            do i = k + 1, n
                if (abs(a(i, k)) > max_abs) then
                    pivot = i
                    max_abs = abs(a(i, k))
                end if
            end do
            if (max_abs <= tiny(1.0_dp)) then
                info = k
                return
            end if
            if (pivot /= k) then
                do j = k, n
                    temp = a(k, j)
                    a(k, j) = a(pivot, j)
                    a(pivot, j) = temp
                end do
                temp = b(k)
                b(k) = b(pivot)
                b(pivot) = temp
            end if
            do i = k + 1, n
                factor = a(i, k) / a(k, k)
                a(i, k) = 0.0_dp
                do j = k + 1, n
                    a(i, j) = a(i, j) - factor * a(k, j)
                end do
                b(i) = b(i) - factor * b(k)
            end do
        end do
        if (abs(a(n, n)) <= tiny(1.0_dp)) then
            info = n
            return
        end if
        do i = n, 1, -1
            if (i < n) b(i) = b(i) - dot_product(a(i, i + 1:n), b(i + 1:n))
            b(i) = b(i) / a(i, i)
        end do
    end subroutine solve_linear_system

end module proxy_numeric_measures
