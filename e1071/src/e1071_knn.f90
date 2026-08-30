module e1071_knn
    use e1071_kinds, only: dp
    use e1071_rng, only: rng_state, rng_integer
    use proxy, only: proxy_dist_cross
    implicit none
    private

    type, public :: gknn_model
        real(dp), allocatable :: x(:, :)
        real(dp), allocatable :: y_real(:)
        integer, allocatable :: y_class(:)
        integer, allocatable :: class_labels(:)
        real(dp), allocatable :: center(:)
        real(dp), allocatable :: scale(:)
        logical, allocatable :: scale_mask(:)
        integer :: k = 1
        logical :: use_all = .true.
        logical :: regression = .false.
        character(len=32) :: method = "euclidean"
    end type gknn_model

    public :: gknn_fit_regression, gknn_fit_classification
    public :: gknn_predict_regression, gknn_predict_classification

contains

    subroutine gknn_fit_regression(x, y, model, k, method, scale, scale_mask, use_all)
        real(dp), intent(in) :: x(:, :) !! Training observation-by-variable matrix.
        real(dp), intent(in) :: y(:) !! Numeric response values, one per training observation.
        type(gknn_model), intent(out) :: model !! Stored scaled training set and generalized-kNN regression settings.
        integer, intent(in), optional :: k !! Number of nearest neighbors before tie expansion; defaults to one.
        character(len=*), intent(in), optional :: method !! proxy distance name; defaults to Euclidean.
        logical, intent(in), optional :: scale !! If true, standardize selected numeric columns; defaults to true.
        logical, intent(in), optional :: scale_mask(:) !! Optional per-column standardization selector; default selects all columns.
        logical, intent(in), optional :: use_all !! If true, include every observation tied at the kth distance; defaults to true.

        if (size(y) /= size(x, 1)) error stop "gknn_fit_regression: x/y row mismatch"
        model%regression = .true.
        model%y_real = y
        call initialize_gknn_x(x, model, k, method, scale, scale_mask, use_all)
    end subroutine gknn_fit_regression

    subroutine gknn_fit_classification(x, y, model, k, method, scale, scale_mask, use_all)
        real(dp), intent(in) :: x(:, :) !! Training observation-by-variable matrix.
        integer, intent(in) :: y(:) !! Integer class labels, one per training observation.
        type(gknn_model), intent(out) :: model !! Stored scaled training set, class labels, and generalized-kNN settings.
        integer, intent(in), optional :: k !! Number of nearest neighbors before tie expansion; defaults to one.
        character(len=*), intent(in), optional :: method !! proxy distance name; defaults to Euclidean.
        logical, intent(in), optional :: scale !! If true, standardize selected columns; defaults to true.
        logical, intent(in), optional :: scale_mask(:) !! Optional per-column standardization selector; default selects all columns.
        logical, intent(in), optional :: use_all !! If true, include every observation tied at the kth distance; defaults to true.

        if (size(y) /= size(x, 1)) error stop "gknn_fit_classification: x/y row mismatch"
        model%regression = .false.
        model%y_class = y
        call unique_labels(y, model%class_labels)
        call initialize_gknn_x(x, model, k, method, scale, scale_mask, use_all)
    end subroutine gknn_fit_classification

    subroutine gknn_predict_regression(model, newdata, prediction, rng)
        type(gknn_model), intent(in) :: model !! Fitted regression gknn model containing training observations and responses.
        real(dp), intent(in) :: newdata(:, :) !! New observation-by-variable matrix with the fitted model's column count.
        real(dp), allocatable, intent(out) :: prediction(:) !! Mean response among selected nearest neighbors for each new row.
        type(rng_state), intent(inout), optional :: rng !! Optional RNG used only when use_all is false and the kth distance is
        !! tied.
        real(dp), allocatable :: query(:, :)
        real(dp), allocatable :: distance(:, :)
        integer, allocatable :: neighbors(:)
        integer :: i
        integer :: status

        if (.not. model%regression) error stop "gknn_predict_regression: classification model supplied"
        call scaled_query(model, newdata, query)
        call proxy_dist_cross(model%x, query, trim(model%method), distance, status=status)
        if (status /= 0) error stop "gknn_predict_regression: proxy distance evaluation failed"
        allocate(prediction(size(query, 1)))
        do i = 1, size(query, 1)
            call select_neighbors(distance(:, i), model%k, model%use_all, neighbors, rng)
            prediction(i) = sum(model%y_real(neighbors)) / real(size(neighbors), dp)
        end do
    end subroutine gknn_predict_regression

    subroutine gknn_predict_classification(model, newdata, class, probability, votes, rng)
        type(gknn_model), intent(in) :: model !! Fitted classification gknn model containing training observations and labels.
        real(dp), intent(in) :: newdata(:, :) !! New observation-by-variable matrix with the fitted model's column count.
        integer, allocatable, intent(out) :: class(:) !! Predicted original integer class label for each new observation.
        real(dp), allocatable, intent(out), optional :: probability(:, :) !! Optional class proportions by row and class_labels
        !! order.
        integer, allocatable, intent(out), optional :: votes(:, :) !! Optional raw neighbor counts by row and class_labels order.
        type(rng_state), intent(inout), optional :: rng !! Optional RNG for kth-distance and class-vote tie breaking.
        real(dp), allocatable :: query(:, :)
        real(dp), allocatable :: distance(:, :)
        integer, allocatable :: neighbors(:)
        integer, allocatable :: count_class(:)
        integer :: i
        integer :: j
        integer :: c
        integer :: maxvote
        integer :: ntied
        integer :: pick
        integer :: status

        if (model%regression) error stop "gknn_predict_classification: regression model supplied"
        call scaled_query(model, newdata, query)
        call proxy_dist_cross(model%x, query, trim(model%method), distance, status=status)
        if (status /= 0) error stop "gknn_predict_classification: proxy distance evaluation failed"
        allocate(class(size(query, 1)), count_class(size(model%class_labels)))
        if (present(probability)) allocate(probability(size(query, 1), size(model%class_labels)))
        if (present(votes)) allocate(votes(size(query, 1), size(model%class_labels)))
        do i = 1, size(query, 1)
            call select_neighbors(distance(:, i), model%k, model%use_all, neighbors, rng)
            count_class = 0
            do j = 1, size(neighbors)
                c = find_label(model%class_labels, model%y_class(neighbors(j)))
                count_class(c) = count_class(c) + 1
            end do
            maxvote = maxval(count_class)
            ntied = count(count_class == maxvote)
            if (ntied == 1 .or. .not. present(rng)) then
                c = maxloc(count_class, dim=1)
            else
                pick = rng_integer(rng, ntied)
                c = nth_true(count_class == maxvote, pick)
            end if
            class(i) = model%class_labels(c)
            if (present(probability)) probability(i, :) = real(count_class, dp) / real(size(neighbors), dp)
            if (present(votes)) votes(i, :) = count_class
        end do
    end subroutine gknn_predict_classification

    subroutine initialize_gknn_x(x, model, k, method, scale, scale_mask, use_all)
        real(dp), intent(in) :: x(:, :) !! Training predictor matrix copied and optionally standardized into model storage.
        type(gknn_model), intent(inout) :: model !! gknn model whose common predictor and control fields are initialized.
        integer, intent(in), optional :: k !! Number of nearest neighbors before tie handling; defaults to one.
        character(len=*), intent(in), optional :: method !! proxy distance name stored in the model; defaults to Euclidean.
        logical, intent(in), optional :: scale !! If true, standardize selected columns; defaults to true.
        logical, intent(in), optional :: scale_mask(:) !! Optional per-column standardization selector.
        logical, intent(in), optional :: use_all !! If true, retain all kth-distance ties; defaults to true.
        logical :: do_scale
        integer :: j
        real(dp) :: mean_value
        real(dp) :: ss

        model%k = 1
        if (present(k)) model%k = k
        if (model%k < 1 .or. model%k > size(x, 1)) error stop "gknn: invalid k"
        model%method = "euclidean"
        if (present(method)) model%method = trim(adjustl(method))
        model%use_all = .true.
        if (present(use_all)) model%use_all = use_all
        do_scale = .true.
        if (present(scale)) do_scale = scale
        allocate(model%scale_mask(size(x, 2)), model%center(size(x, 2)), model%scale(size(x, 2)))
        model%scale_mask = do_scale
        if (present(scale_mask)) then
            if (size(scale_mask) /= size(x, 2)) error stop "gknn: scale_mask has wrong length"
            if (do_scale) model%scale_mask = scale_mask
        end if
        model%x = x
        model%center = 0.0_dp
        model%scale = 1.0_dp
        do j = 1, size(x, 2)
            if (.not. model%scale_mask(j)) cycle
            mean_value = sum(model%x(:, j)) / real(size(model%x, 1), dp)
            ss = sum((model%x(:, j) - mean_value)**2)
            if (size(model%x, 1) <= 1 .or. ss <= 0.0_dp) error stop "gknn: cannot scale a constant column"
            model%center(j) = mean_value
            model%scale(j) = sqrt(ss / real(size(model%x, 1) - 1, dp))
            model%x(:, j) = (model%x(:, j) - model%center(j)) / model%scale(j)
        end do
    end subroutine initialize_gknn_x

    subroutine scaled_query(model, newdata, query)
        type(gknn_model), intent(in) :: model !! Fitted model supplying expected column count and scaling constants.
        real(dp), intent(in) :: newdata(:, :) !! New observation matrix to validate and transform.
        real(dp), allocatable, intent(out) :: query(:, :) !! Allocated copy of newdata transformed by model scaling.
        integer :: j

        if (size(newdata, 2) /= size(model%x, 2)) error stop "gknn prediction: variable count mismatch"
        query = newdata
        do j = 1, size(query, 2)
            if (model%scale_mask(j)) query(:, j) = (query(:, j) - model%center(j)) / model%scale(j)
        end do
    end subroutine scaled_query

    subroutine select_neighbors(distance, k, use_all, neighbors, rng)
        real(dp), intent(in) :: distance(:) !! Distances from every training observation to one query observation.
        integer, intent(in) :: k !! Nominal neighbor count; must be within the distance-vector length.
        logical, intent(in) :: use_all !! If true, include every training row tied with the kth nearest distance.
        integer, allocatable, intent(out) :: neighbors(:) !! Selected one-based training-row indices.
        type(rng_state), intent(inout), optional :: rng !! Optional RNG for choosing one kth-distance tie when use_all is false.
        real(dp), allocatable :: sorted(:)
        integer, allocatable :: order(:)
        integer, allocatable :: tied(:)
        real(dp) :: key
        real(dp) :: kth
        integer :: ikey
        integer :: i
        integer :: j
        integer :: ntie
        integer :: pick

        if (k < 1 .or. k > size(distance)) error stop "select_neighbors: invalid k"
        sorted = distance
        allocate(order(size(distance)))
        order = [(i, i = 1, size(distance))]
        do i = 2, size(sorted)
            key = sorted(i)
            ikey = order(i)
            j = i - 1
            do while (j >= 1)
                if (sorted(j) <= key) exit
                sorted(j + 1) = sorted(j)
                order(j + 1) = order(j)
                j = j - 1
            end do
            sorted(j + 1) = key
            order(j + 1) = ikey
        end do
        kth = sorted(k)
        if (use_all) then
            allocate(neighbors(count(distance <= kth)))
            j = 0
            do i = 1, size(distance)
                if (distance(i) <= kth) then
                    j = j + 1
                    neighbors(j) = i
                end if
            end do
        else if (k == 1 .or. sorted(k - 1) < kth) then
            ntie = count(abs(distance - kth) <= 0.0_dp)
            if (ntie <= 1 .or. .not. present(rng)) then
                allocate(neighbors(k))
                neighbors = order(:k)
            else
                allocate(tied(ntie))
                j = 0
                do i = 1, size(distance)
                    if (abs(distance(i) - kth) <= 0.0_dp) then
                        j = j + 1
                        tied(j) = i
                    end if
                end do
                pick = rng_integer(rng, ntie)
                allocate(neighbors(k))
                if (k > 1) neighbors(:k - 1) = order(:k - 1)
                neighbors(k) = tied(pick)
            end if
        else
            allocate(neighbors(k))
            neighbors = order(:k)
        end if
    end subroutine select_neighbors

    subroutine unique_labels(y, labels)
        integer, intent(in) :: y(:) !! Integer class labels from which sorted unique values are extracted.
        integer, allocatable, intent(out) :: labels(:) !! Sorted distinct labels present in y.
        integer, allocatable :: work(:)
        integer :: i
        integer :: j
        integer :: n
        integer :: key

        work = y
        do i = 2, size(work)
            key = work(i)
            j = i - 1
            do while (j >= 1)
                if (work(j) <= key) exit
                work(j + 1) = work(j)
                j = j - 1
            end do
            work(j + 1) = key
        end do
        n = 1
        do i = 2, size(work)
            if (work(i) /= work(i - 1)) n = n + 1
        end do
        allocate(labels(n))
        labels(1) = work(1)
        j = 1
        do i = 2, size(work)
            if (work(i) /= work(i - 1)) then
                j = j + 1
                labels(j) = work(i)
            end if
        end do
    end subroutine unique_labels

    pure function find_label(labels, value) result(index)
        integer, intent(in) :: labels(:) !! Distinct class labels to search.
        integer, intent(in) :: value !! Original class label whose one-based position is requested.
        integer :: index
        integer :: i

        index = 0
        do i = 1, size(labels)
            if (labels(i) == value) then
                index = i
                return
            end if
        end do
    end function find_label

    pure function nth_true(mask, n) result(index)
        logical, intent(in) :: mask(:) !! Logical selection vector containing at least n true entries.
        integer, intent(in) :: n !! One-based rank among true entries to return.
        integer :: index
        integer :: i
        integer :: seen

        seen = 0
        index = 0
        do i = 1, size(mask)
            if (.not. mask(i)) cycle
            seen = seen + 1
            if (seen == n) then
                index = i
                return
            end if
        end do
    end function nth_true

end module e1071_knn
