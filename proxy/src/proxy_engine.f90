module proxy_engine
    use proxy_kinds, only: dp
    use proxy_ieee, only: proxy_nan
    use proxy_utils, only: proxy_convert_default, proxy_convert_one_minus, &
                           proxy_simil_to_dist, proxy_dist_to_simil, proxy_normalize_name
    use proxy_numeric_measures
    use proxy_binary_measures, only: binary_counts_numeric, binary_similarity_from_counts
    use proxy_nominal_measures, only: nominal_similarity
    use proxy_gower, only: gower_auto_similarity, gower_cross_similarity, proxy_gower_metric
    use proxy_registry, only: proxy_numeric_callback, proxy_binary_callback, &
                              proxy_lookup_numeric, proxy_lookup_binary
    implicit none
    private

    integer, parameter :: category_numeric = 1
    integer, parameter :: category_binary = 2
    integer, parameter :: category_nominal = 3
    integer, parameter :: category_gower = 4

    public :: proxy_dist_auto, proxy_dist_cross, proxy_dist_pairwise
    public :: proxy_simil_auto, proxy_simil_cross, proxy_simil_pairwise
    public :: proxy_measure_value, proxy_builtin_names

contains

    subroutine proxy_dist_auto(x, method, distance, p, covariance, status)
        real(dp), intent(in) :: x(:, :) !! Observation-by-variable numeric matrix; rows are objects and NaN values represent
        !! missing components.
        character(len=*), intent(in) :: method !! Registered or built-in proximity name; similarity methods are converted to
        !! distances using proxy rules.
        real(dp), allocatable, intent(out) :: distance(:, :) !! Allocated symmetric square distance matrix for the rows of `x`.
        real(dp), intent(in), optional :: p(:) !! Optional method parameter vector; Minkowski uses `p(1)` as its positive
        !! exponent and custom callbacks receive the full vector.
        real(dp), intent(in), optional :: covariance(:, :) !! Optional covariance matrix for Mahalanobis distance; otherwise it
        !! is estimated from `x`.
        integer, intent(out), optional :: status !! Zero on success; nonzero identifies an unknown method, invalid shape, or
        !! unavailable Mahalanobis covariance.
        real(dp), allocatable :: cov(:, :)
        integer :: i
        integer :: j
        integer :: info

        if (present(status)) status = 0
        if (is_gower(method)) then
            call gower_as_distance_auto(x, distance)
            return
        end if
        if (is_mahalanobis(method)) then
            call choose_auto_covariance(x, covariance, cov, info)
            if (info /= 0) then
                allocate(distance(0, 0))
                if (present(status)) status = info
                return
            end if
        else
            allocate(cov(0, 0))
        end if
        allocate(distance(size(x, 1), size(x, 1)))
        do j = 1, size(x, 1)
            do i = 1, size(x, 1)
                distance(i, j) = proxy_measure_value(x(i, :), x(j, :), method, .true., p, cov, info)
                if (info /= 0 .and. present(status)) status = info
            end do
        end do
    end subroutine proxy_dist_auto

    subroutine proxy_dist_cross(x, y, method, distance, p, covariance, status)
        real(dp), intent(in) :: x(:, :) !! First observation-by-variable matrix defining rows of the returned cross-distance matrix.
        real(dp), intent(in) :: y(:, :) !! Second observation-by-variable matrix defining columns; it must have the same number
        !! of variables as `x`.
        character(len=*), intent(in) :: method !! Built-in or registered method name; similarity methods are converted to
        !! distances automatically.
        real(dp), allocatable, intent(out) :: distance(:, :) !! Allocated `size(x,1)` by `size(y,1)` cross-distance matrix.
        real(dp), intent(in), optional :: p(:) !! Optional method parameter vector forwarded to built-in/custom methods;
        !! Minkowski uses its first value.
        real(dp), intent(in), optional :: covariance(:, :) !! Optional Mahalanobis covariance/cross-covariance matrix; otherwise
        !! proxy-style `cov(x,y)` is estimated when row counts agree.
        integer, intent(out), optional :: status !! Zero on success; nonzero for nonconforming matrices, unknown methods, or
        !! failed covariance estimation.
        real(dp), allocatable :: cov(:, :)
        integer :: i
        integer :: j
        integer :: info

        if (present(status)) status = 0
        if (size(x, 2) /= size(y, 2)) then
            allocate(distance(0, 0))
            if (present(status)) status = 3
            return
        end if
        if (is_gower(method)) then
            call gower_as_distance_cross(x, y, distance)
            return
        end if
        if (is_mahalanobis(method)) then
            call choose_cross_covariance(x, y, covariance, cov, info)
            if (info /= 0) then
                allocate(distance(0, 0))
                if (present(status)) status = info
                return
            end if
        else
            allocate(cov(0, 0))
        end if
        allocate(distance(size(x, 1), size(y, 1)))
        do j = 1, size(y, 1)
            do i = 1, size(x, 1)
                distance(i, j) = proxy_measure_value(x(i, :), y(j, :), method, .true., p, cov, info)
                if (info /= 0 .and. present(status)) status = info
            end do
        end do
    end subroutine proxy_dist_cross

    subroutine proxy_dist_pairwise(x, y, method, distance, p, covariance, status)
        real(dp), intent(in) :: x(:, :) !! First observation matrix for row-by-row comparisons; its row and column counts must
        !! match `y`.
        real(dp), intent(in) :: y(:, :) !! Second observation matrix paired rowwise with `x`.
        character(len=*), intent(in) :: method !! Built-in or registered method name, converted to a distance if it natively
        !! returns similarity.
        real(dp), allocatable, intent(out) :: distance(:) !! Allocated vector where element `i` compares rows `x(i,:)` and `y(i,:)`.
        real(dp), intent(in), optional :: p(:) !! Optional parameter vector for the selected method.
        real(dp), intent(in), optional :: covariance(:, :) !! Optional Mahalanobis covariance/cross-covariance matrix; otherwise
        !! it is estimated from paired `x` and `y`.
        integer, intent(out), optional :: status !! Zero on success; nonzero for shape errors, unknown methods, or covariance
        !! failure.
        real(dp), allocatable :: cov(:, :)
        real(dp), allocatable :: tmp(:, :)
        integer :: i
        integer :: info

        if (present(status)) status = 0
        if (size(x, 1) /= size(y, 1) .or. size(x, 2) /= size(y, 2)) then
            allocate(distance(0))
            if (present(status)) status = 3
            return
        end if
        if (is_gower(method)) then
            call gower_cross_similarity(x, y, spread(proxy_gower_metric, 1, size(x, 2)), tmp, pairwise=.true.)
            allocate(distance(size(x, 1)))
            distance = 1.0_dp - abs(tmp(:, 1))
            return
        end if
        if (is_mahalanobis(method)) then
            call choose_cross_covariance(x, y, covariance, cov, info)
            if (info /= 0) then
                allocate(distance(0))
                if (present(status)) status = info
                return
            end if
        else
            allocate(cov(0, 0))
        end if
        allocate(distance(size(x, 1)))
        do i = 1, size(x, 1)
            distance(i) = proxy_measure_value(x(i, :), y(i, :), method, .true., p, cov, info)
            if (info /= 0 .and. present(status)) status = info
        end do
    end subroutine proxy_dist_pairwise

    subroutine proxy_simil_auto(x, method, similarity, p, covariance, status)
        real(dp), intent(in) :: x(:, :) !! Observation-by-variable matrix whose rows are compared with every other row.
        character(len=*), intent(in) :: method !! Built-in or registered method name; distance methods are converted to
        !! similarities using proxy rules.
        real(dp), allocatable, intent(out) :: similarity(:, :) !! Allocated symmetric similarity matrix for all rows of `x`.
        real(dp), intent(in), optional :: p(:) !! Optional method parameter vector; custom methods receive it unchanged.
        real(dp), intent(in), optional :: covariance(:, :) !! Optional Mahalanobis covariance matrix; otherwise it is estimated
        !! from `x` when needed.
        integer, intent(out), optional :: status !! Zero on success; nonzero on invalid input, method lookup failure, or
        !! covariance failure.
        real(dp), allocatable :: cov(:, :)
        integer :: i
        integer :: j
        integer :: info

        if (present(status)) status = 0
        if (is_gower(method)) then
            call gower_auto_similarity(x, spread(proxy_gower_metric, 1, size(x, 2)), similarity)
            return
        end if
        if (is_mahalanobis(method)) then
            call choose_auto_covariance(x, covariance, cov, info)
            if (info /= 0) then
                allocate(similarity(0, 0))
                if (present(status)) status = info
                return
            end if
        else
            allocate(cov(0, 0))
        end if
        allocate(similarity(size(x, 1), size(x, 1)))
        do j = 1, size(x, 1)
            do i = 1, size(x, 1)
                similarity(i, j) = proxy_measure_value(x(i, :), x(j, :), method, .false., p, cov, info)
                if (info /= 0 .and. present(status)) status = info
            end do
        end do
    end subroutine proxy_simil_auto

    subroutine proxy_simil_cross(x, y, method, similarity, p, covariance, status)
        real(dp), intent(in) :: x(:, :) !! First observation matrix defining rows of the returned similarity matrix.
        real(dp), intent(in) :: y(:, :) !! Second observation matrix defining columns; it must have the same number of variables
        !! as `x`.
        character(len=*), intent(in) :: method !! Built-in or registered proximity name; distance methods are converted to
        !! similarities.
        real(dp), allocatable, intent(out) :: similarity(:, :) !! Allocated cross-similarity matrix with shape
        !! `(size(x,1),size(y,1))`.
        real(dp), intent(in), optional :: p(:) !! Optional method parameter vector forwarded to the selected method.
        real(dp), intent(in), optional :: covariance(:, :) !! Optional covariance matrix for Mahalanobis calculations; otherwise
        !! cross-covariance is estimated when possible.
        integer, intent(out), optional :: status !! Zero on success; nonzero for method, dimension, or covariance errors.
        real(dp), allocatable :: cov(:, :)
        integer :: i
        integer :: j
        integer :: info

        if (present(status)) status = 0
        if (size(x, 2) /= size(y, 2)) then
            allocate(similarity(0, 0))
            if (present(status)) status = 3
            return
        end if
        if (is_gower(method)) then
            call gower_cross_similarity(x, y, spread(proxy_gower_metric, 1, size(x, 2)), similarity)
            return
        end if
        if (is_mahalanobis(method)) then
            call choose_cross_covariance(x, y, covariance, cov, info)
            if (info /= 0) then
                allocate(similarity(0, 0))
                if (present(status)) status = info
                return
            end if
        else
            allocate(cov(0, 0))
        end if
        allocate(similarity(size(x, 1), size(y, 1)))
        do j = 1, size(y, 1)
            do i = 1, size(x, 1)
                similarity(i, j) = proxy_measure_value(x(i, :), y(j, :), method, .false., p, cov, info)
                if (info /= 0 .and. present(status)) status = info
            end do
        end do
    end subroutine proxy_simil_cross

    subroutine proxy_simil_pairwise(x, y, method, similarity, p, covariance, status)
        real(dp), intent(in) :: x(:, :) !! First matrix for corresponding-row similarity comparisons.
        real(dp), intent(in) :: y(:, :) !! Second matrix paired rowwise with `x`; row and column counts must match.
        character(len=*), intent(in) :: method !! Built-in or registered proximity method name.
        real(dp), allocatable, intent(out) :: similarity(:) !! Allocated vector of pairwise row similarities.
        real(dp), intent(in), optional :: p(:) !! Optional parameter vector forwarded to the selected method.
        real(dp), intent(in), optional :: covariance(:, :) !! Optional Mahalanobis covariance/cross-covariance matrix.
        integer, intent(out), optional :: status !! Zero on success; nonzero for shape, lookup, or covariance errors.
        real(dp), allocatable :: cov(:, :)
        real(dp), allocatable :: tmp(:, :)
        integer :: i
        integer :: info

        if (present(status)) status = 0
        if (size(x, 1) /= size(y, 1) .or. size(x, 2) /= size(y, 2)) then
            allocate(similarity(0))
            if (present(status)) status = 3
            return
        end if
        if (is_gower(method)) then
            call gower_cross_similarity(x, y, spread(proxy_gower_metric, 1, size(x, 2)), tmp, pairwise=.true.)
            allocate(similarity(size(x, 1)))
            similarity = tmp(:, 1)
            return
        end if
        if (is_mahalanobis(method)) then
            call choose_cross_covariance(x, y, covariance, cov, info)
            if (info /= 0) then
                allocate(similarity(0))
                if (present(status)) status = info
                return
            end if
        else
            allocate(cov(0, 0))
        end if
        allocate(similarity(size(x, 1)))
        do i = 1, size(x, 1)
            similarity(i) = proxy_measure_value(x(i, :), y(i, :), method, .false., p, cov, info)
            if (info /= 0 .and. present(status)) status = info
        end do
    end subroutine proxy_simil_pairwise

    function proxy_measure_value(x, y, method, want_distance, p, covariance, status) result(value)
        real(dp), intent(in) :: x(:) !! First observation vector passed to the selected built-in or registered method.
        real(dp), intent(in) :: y(:) !! Second observation vector; it must have the same number of components as `x` for built-in
        !! numeric methods.
        character(len=*), intent(in) :: method !! Built-in alias or custom registry name.
        logical, intent(in) :: want_distance !! True requests distance output; false requests similarity output, applying the
        !! method's registered conversion when necessary.
        real(dp), intent(in), optional :: p(:) !! Optional method parameter vector; Minkowski requires `p(1)` and custom numeric
        !! callbacks receive the entire vector.
        real(dp), intent(in), optional :: covariance(:, :) !! Optional covariance matrix used only by Mahalanobis distance.
        integer, intent(out) :: status !! Zero when the method succeeds; one for unknown method and other positive values for
        !! invalid method-specific input.
        real(dp) :: value
        real(dp) :: raw
        logical :: found
        logical :: is_distance
        integer :: conversion
        integer :: category
        integer :: a
        integer :: b
        integer :: c
        integer :: d
        integer :: n
        procedure(proxy_numeric_callback), pointer :: custom_numeric
        procedure(proxy_binary_callback), pointer :: custom_binary

        status = 0
        call builtin_info(method, found, is_distance, conversion, category)
        if (found) then
            select case (category)
            case (category_numeric)
                raw = numeric_builtin(x, y, method, p, covariance, status)
            case (category_binary)
                call binary_counts_numeric(x, y, a, b, c, d, n)
                raw = binary_similarity_from_counts(method, a, b, c, d, n)
            case (category_nominal)
                raw = nominal_similarity(x, y, method)
            case default
                raw = proxy_nan()
                status = 1
            end select
        else
            call proxy_lookup_numeric(method, custom_numeric, found, is_distance, conversion)
            if (found) then
                if (present(p)) then
                    raw = custom_numeric(x, y, p)
                else
                    raw = custom_numeric(x, y)
                end if
            else
                call proxy_lookup_binary(method, custom_binary, found, is_distance, conversion)
                if (.not. found) then
                    value = proxy_nan()
                    status = 1
                    return
                end if
                call binary_counts_numeric(x, y, a, b, c, d, n)
                raw = custom_binary(a, b, c, d, n)
            end if
        end if

        if (status /= 0) then
            value = raw
        else if (want_distance .and. .not. is_distance) then
            value = proxy_simil_to_dist(raw, conversion)
        else if (.not. want_distance .and. is_distance) then
            value = proxy_dist_to_simil(raw, conversion)
        else
            value = raw
        end if
    end function proxy_measure_value

    function numeric_builtin(x, y, method, p, covariance, status) result(value)
        real(dp), intent(in) :: x(:) !! First numeric observation vector for built-in metric/similarity evaluation.
        real(dp), intent(in) :: y(:) !! Second numeric observation vector aligned with `x`.
        character(len=*), intent(in) :: method !! Normal or aliased built-in numeric method name.
        real(dp), intent(in), optional :: p(:) !! Optional parameter vector; only Minkowski currently requires it.
        real(dp), intent(in), optional :: covariance(:, :) !! Optional covariance matrix required for Mahalanobis evaluation at
        !! this level.
        integer, intent(out) :: status !! Zero on success; positive on a missing required parameter or unknown numeric method.
        real(dp) :: value
        character(len=:), allocatable :: key

        status = 0
        key = proxy_normalize_name(method)
        select case (key)
        case ('euclidean', 'l2')
            value = euclidean_distance(x, y)
        case ('mahalanobis')
            if (.not. present(covariance) .or. size(covariance) == 0) then
                value = proxy_nan()
                status = 4
            else
                value = mahalanobis_distance(x, y, covariance)
            end if
        case ('bhjattacharyya')
            value = bhjattacharyya_distance(x, y)
        case ('manhattan', 'cityblock', 'l1', 'taxi')
            value = manhattan_distance(x, y)
        case ('supremum', 'max', 'maximum', 'tschebyscheff', 'chebyshev')
            value = supremum_distance(x, y)
        case ('minkowski', 'lp')
            if (.not. present(p) .or. size(p) == 0) then
                value = proxy_nan()
                status = 5
            else
                value = minkowski_distance(x, y, p(1))
            end if
        case ('canberra')
            value = canberra_distance(x, y)
        case ('wave', 'hedges', 'wavehedges')
            value = wave_hedges_distance(x, y)
        case ('divergence')
            value = divergence_distance(x, y)
        case ('kullback', 'leibler', 'kullbackleibler')
            value = kullback_leibler_distance(x, y)
        case ('bray', 'curtis', 'braycurtis')
            value = bray_curtis_distance(x, y)
        case ('soergel')
            value = soergel_distance(x, y)
        case ('podani', 'discordance')
            value = podani_distance(x, y)
        case ('chord')
            value = chord_distance(x, y)
        case ('geodesic')
            value = geodesic_distance(x, y)
        case ('whittaker')
            value = whittaker_distance(x, y)
        case ('hellinger')
            value = hellinger_distance(x, y)
        case ('fjaccard', 'fuzzyjaccard')
            value = fuzzy_jaccard_distance(x, y)
        case ('cosine')
            value = cosine_similarity(x, y)
        case ('angular')
            value = angular_similarity(x, y)
        case ('ejaccard', 'extendedjaccard')
            value = extended_jaccard_similarity(x, y)
        case ('edice', 'extendeddice', 'esorensen')
            value = extended_dice_similarity(x, y)
        case ('correlation')
            value = correlation_similarity(x, y)
        case ('mutual')
            value = mutual_information_similarity(x, y)
        case default
            value = proxy_nan()
            status = 1
        end select
    end function numeric_builtin

    subroutine builtin_info(method, found, is_distance, conversion, category)
        character(len=*), intent(in) :: method !! Built-in method name or alias to classify.
        logical, intent(out) :: found !! True when the normalized name belongs to a built-in method.
        logical, intent(out) :: is_distance !! True for native distances/dissimilarities and false for native similarities.
        integer, intent(out) :: conversion !! Conversion rule used when callers request the opposite proximity type.
        integer, intent(out) :: category !! Internal dispatch category for numeric, binary, nominal, or Gower methods.
        character(len=:), allocatable :: key

        key = proxy_normalize_name(method)
        found = .true.
        conversion = proxy_convert_default
        category = category_numeric
        select case (key)
        case ('euclidean', 'l2', 'mahalanobis', 'bhjattacharyya', 'manhattan', 'cityblock', 'l1', 'taxi', &
              'supremum', 'max', 'maximum', 'tschebyscheff', 'chebyshev', 'minkowski', 'lp', 'canberra', &
              'wave', 'hedges', 'wavehedges', 'divergence', 'kullback', 'leibler', 'kullbackleibler', &
              'bray', 'curtis', 'braycurtis', 'soergel', 'podani', 'discordance', 'chord', 'geodesic', &
              'whittaker', 'hellinger', 'fjaccard', 'fuzzyjaccard')
            is_distance = .true.
        case ('cosine', 'angular')
            is_distance = .false.
            conversion = proxy_convert_one_minus
        case ('ejaccard', 'extendedjaccard', 'edice', 'extendeddice', 'esorensen', 'correlation')
            is_distance = .false.
        case ('jaccard', 'binary', 'reyssac', 'roux', 'kulczynski1', 'kulczynski2', 'mountford', 'fager', &
              'mcgowan', 'fagermcgowan', 'russel', 'rao', 'russelrao', 'simplematching', 'sokalmichener', &
              'hamman', 'faith', 'tanimoto', 'rogers', 'rogerstanimoto', 'dice', 'czekanowski', 'sorensen', &
              'phi', 'stiles', 'michael', 'mozley', 'margalef', 'mozleymargalef', 'yule', 'yule2', &
              'ochiai', 'simpson', 'braunblanquet')
            is_distance = .false.
            category = category_binary
        case ('chisquared', 'phisquared', 'tschuprow', 'cramer', 'pearson', 'contingency')
            is_distance = .false.
            category = category_nominal
        case ('gower')
            is_distance = .false.
            category = category_gower
        case ('mutual')
            is_distance = .false.
        case default
            found = .false.
            is_distance = .false.
            category = 0
        end select
    end subroutine builtin_info

    subroutine proxy_builtin_names(names)
        character(len=24), allocatable, intent(out) :: names(:) !! Allocated canonical names of all translated built-in proxy
        !! methods, excluding aliases.

        allocate(names(50))
        names = [character(len=24) :: &
            'Euclidean', 'Mahalanobis', 'Bhjattacharyya', 'Manhattan', 'supremum', 'Minkowski', 'Canberra', &
            'WaveHedges', 'divergence', 'KullbackLeibler', 'BrayCurtis', 'Soergel', 'Podani', 'Chord', &
            'Geodesic', 'Whittaker', 'Hellinger', 'fJaccard', 'Jaccard', 'Kulczynski1', 'Kulczynski2', &
            'Mountford', 'FagerMcGowan', 'RusselRao', 'simple matching', 'Hamman', 'Faith', 'RogersTanimoto', &
            'Dice', 'Phi', 'Stiles', 'Michael', 'MozleyMargalef', 'Yule', 'Yule2', 'Ochiai', 'Simpson', &
            'Braun-Blanquet', 'cosine', 'angular', 'eJaccard', 'eDice', 'correlation', 'Chi-squared', &
            'Phi-squared', 'Tschuprow', 'Cramer', 'Pearson', 'Gower', 'mutual']
    end subroutine proxy_builtin_names

    subroutine choose_auto_covariance(x, supplied, covariance, status)
        real(dp), intent(in) :: x(:, :) !! Observation matrix whose sample covariance is estimated when `supplied` is absent.
        real(dp), intent(in), optional :: supplied(:, :) !! Optional caller-provided covariance matrix, copied unchanged when
        !! conformable.
        real(dp), allocatable, intent(out) :: covariance(:, :) !! Allocated covariance matrix used by Mahalanobis calculations.
        integer, intent(out) :: status !! Zero on success; nonzero when the supplied matrix is nonconforming or too few rows
        !! exist for estimation.

        if (present(supplied)) then
            if (size(supplied, 1) /= size(x, 2) .or. size(supplied, 2) /= size(x, 2)) then
                allocate(covariance(0, 0))
                status = 4
            else
                covariance = supplied
                status = 0
            end if
        else
            call covariance_matrix(x, covariance, status)
        end if
    end subroutine choose_auto_covariance

    subroutine choose_cross_covariance(x, y, supplied, covariance, status)
        real(dp), intent(in) :: x(:, :) !! First observation matrix used in proxy-style cross-covariance estimation.
        real(dp), intent(in) :: y(:, :) !! Second observation matrix; automatic cross-covariance requires the same row count as `x`.
        real(dp), intent(in), optional :: supplied(:, :) !! Optional caller-provided square covariance/cross-covariance matrix.
        real(dp), allocatable, intent(out) :: covariance(:, :) !! Allocated matrix passed to the Mahalanobis evaluator.
        integer, intent(out) :: status !! Zero on success; nonzero for dimension mismatch or insufficient rows.

        if (present(supplied)) then
            if (size(supplied, 1) /= size(x, 2) .or. size(supplied, 2) /= size(x, 2)) then
                allocate(covariance(0, 0))
                status = 4
            else
                covariance = supplied
                status = 0
            end if
        else
            call cross_covariance_matrix(x, y, covariance, status)
        end if
    end subroutine choose_cross_covariance

    subroutine covariance_matrix(x, covariance, status)
        real(dp), intent(in) :: x(:, :) !! Observation-by-variable matrix used to estimate the ordinary unbiased sample
        !! covariance matrix.
        real(dp), allocatable, intent(out) :: covariance(:, :) !! Allocated `p` by `p` sample covariance matrix.
        integer, intent(out) :: status !! Zero on success; nonzero when fewer than two observations are available.
        real(dp), allocatable :: mean(:)
        real(dp), allocatable :: centered(:, :)
        integer :: n

        n = size(x, 1)
        if (n < 2) then
            allocate(covariance(0, 0))
            status = 4
            return
        end if
        allocate(mean(size(x, 2)), centered(n, size(x, 2)))
        mean = sum(x, dim=1) / real(n, dp)
        centered = x - spread(mean, 1, n)
        covariance = matmul(transpose(centered), centered) / real(n - 1, dp)
        status = 0
    end subroutine covariance_matrix

    subroutine cross_covariance_matrix(x, y, covariance, status)
        real(dp), intent(in) :: x(:, :) !! First observation matrix in the cross-covariance estimate.
        real(dp), intent(in) :: y(:, :) !! Second matrix; rows are paired with `x` and columns must conform.
        real(dp), allocatable, intent(out) :: covariance(:, :) !! Allocated cross-covariance matrix analogous to R `cov(x,y)`.
        integer, intent(out) :: status !! Zero on success; nonzero when row/column counts do not conform or fewer than two rows
        !! are present.
        real(dp), allocatable :: mx(:)
        real(dp), allocatable :: my(:)
        real(dp), allocatable :: cx(:, :)
        real(dp), allocatable :: cy(:, :)
        integer :: n

        n = size(x, 1)
        if (n < 2 .or. size(y, 1) /= n .or. size(y, 2) /= size(x, 2)) then
            allocate(covariance(0, 0))
            status = 4
            return
        end if
        allocate(mx(size(x, 2)), my(size(y, 2)), cx(n, size(x, 2)), cy(n, size(y, 2)))
        mx = sum(x, dim=1) / real(n, dp)
        my = sum(y, dim=1) / real(n, dp)
        cx = x - spread(mx, 1, n)
        cy = y - spread(my, 1, n)
        covariance = matmul(transpose(cx), cy) / real(n - 1, dp)
        status = 0
    end subroutine cross_covariance_matrix

    pure function is_gower(method) result(answer)
        character(len=*), intent(in) :: method !! Method name to test for the Gower special preprocessing path.
        logical :: answer

        answer = proxy_normalize_name(method) == 'gower'
    end function is_gower

    pure function is_mahalanobis(method) result(answer)
        character(len=*), intent(in) :: method !! Method name to test for automatic covariance preprocessing.
        logical :: answer

        answer = proxy_normalize_name(method) == 'mahalanobis'
    end function is_mahalanobis

    subroutine gower_as_distance_auto(x, distance)
        real(dp), intent(in) :: x(:, :) !! Numeric matrix treated as all-metric variables for explicit Gower distance requests.
        real(dp), allocatable, intent(out) :: distance(:, :) !! Allocated Gower dissimilarity matrix obtained via proxy's
        !! `1-abs(similarity)` conversion.
        real(dp), allocatable :: similarity(:, :)

        call gower_auto_similarity(x, spread(proxy_gower_metric, 1, size(x, 2)), similarity)
        distance = 1.0_dp - abs(similarity)
    end subroutine gower_as_distance_auto

    subroutine gower_as_distance_cross(x, y, distance)
        real(dp), intent(in) :: x(:, :) !! First all-metric numeric matrix for Gower cross-distance evaluation.
        real(dp), intent(in) :: y(:, :) !! Second all-metric numeric matrix with columns conforming to `x`.
        real(dp), allocatable, intent(out) :: distance(:, :) !! Allocated cross-dissimilarity matrix using the standard proxy
        !! similarity conversion.
        real(dp), allocatable :: similarity(:, :)

        call gower_cross_similarity(x, y, spread(proxy_gower_metric, 1, size(x, 2)), similarity)
        distance = 1.0_dp - abs(similarity)
    end subroutine gower_as_distance_cross

end module proxy_engine
