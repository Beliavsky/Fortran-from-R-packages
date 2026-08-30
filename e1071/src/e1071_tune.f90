module e1071_tune
    use e1071_kinds, only: dp
    use e1071_svm, only: svm_options, svm_model, svm_fit_classification, svm_fit_regression, &
                         svm_predict_classification, svm_predict_regression
    use e1071_knn, only: gknn_model, gknn_fit_classification, gknn_fit_regression, &
                         gknn_predict_classification, gknn_predict_regression
    implicit none
    private

    type, public :: tune_result
        real(dp), allocatable :: parameter1(:)
        real(dp), allocatable :: parameter2(:)
        real(dp), allocatable :: error(:)
        real(dp), allocatable :: dispersion(:)
        integer :: best_index = 0
        real(dp) :: best_performance = huge(1.0_dp)
    end type tune_result

    public :: tune_svm_classification, tune_svm_regression
    public :: tune_gknn_classification, tune_gknn_regression

contains

    subroutine tune_svm_classification(x, y, costs, gammas, folds, result, base_options, fold_id)
        real(dp), intent(in) :: x(:, :) !! Classification predictor matrix used by every parameter combination.
        integer, intent(in) :: y(:) !! Integer class labels corresponding rowwise to x.
        real(dp), intent(in) :: costs(:) !! Positive C values crossed with every value in `gammas`.
        real(dp), intent(in) :: gammas(:) !! Nonnegative kernel gamma values crossed with every value in `costs`.
        integer, intent(in) :: folds !! Number of validation folds; must be between two and the number of observations.
        type(tune_result), intent(out) :: result !! Grid parameters, mean error, fold dispersion, and best combination index.
        type(svm_options), intent(in), optional :: base_options !! SVM controls copied before cost and gamma are overwritten.
        integer, intent(in), optional :: fold_id(:) !! Optional one-based fold labels; default uses deterministic cyclic folds.
        type(svm_options) :: options
        integer, allocatable :: prediction(:)
        real(dp), allocatable :: fold_error(:)
        integer :: icost
        integer :: igamma
        integer :: grid

        call validate_tune_inputs(size(x, 1), folds, fold_id)
        if (size(y) /= size(x, 1)) error stop "tune_svm_classification: x/y mismatch"
        if (size(costs) < 1 .or. size(gammas) < 1) error stop "tune_svm_classification: empty parameter grid"
        if (any(costs <= 0.0_dp) .or. any(gammas < 0.0_dp)) error stop "tune_svm_classification: invalid grid value"
        allocate(result%parameter1(size(costs) * size(gammas)), result%parameter2(size(costs) * size(gammas)))
        allocate(result%error(size(result%parameter1)), result%dispersion(size(result%parameter1)))
        allocate(fold_error(folds))
        options = svm_options()
        if (present(base_options)) options = base_options
        grid = 0
        do igamma = 1, size(gammas)
            do icost = 1, size(costs)
                grid = grid + 1
                options%cost = costs(icost)
                options%gamma = gammas(igamma)
                call cv_svm_classification(x, y, folds, options, prediction, fold_error, fold_id)
                result%parameter1(grid) = costs(icost)
                result%parameter2(grid) = gammas(igamma)
                result%error(grid) = sum(fold_error) / real(folds, dp)
                result%dispersion(grid) = sample_sd(fold_error)
            end do
        end do
        call select_best(result)
    end subroutine tune_svm_classification

    subroutine tune_svm_regression(x, y, costs, gammas, folds, result, base_options, fold_id)
        real(dp), intent(in) :: x(:, :) !! Regression predictor matrix used by every parameter combination.
        real(dp), intent(in) :: y(:) !! Numeric response vector corresponding rowwise to x.
        real(dp), intent(in) :: costs(:) !! Positive C values crossed with every value in `gammas`.
        real(dp), intent(in) :: gammas(:) !! Nonnegative kernel gamma values crossed with every value in `costs`.
        integer, intent(in) :: folds !! Number of validation folds; must be between two and the number of observations.
        type(tune_result), intent(out) :: result !! Grid parameters, mean squared errors, dispersion, and best combination index.
        type(svm_options), intent(in), optional :: base_options !! SVM controls copied before cost and gamma are overwritten.
        integer, intent(in), optional :: fold_id(:) !! Optional one-based fold labels; default uses deterministic cyclic folds.
        type(svm_options) :: options
        real(dp), allocatable :: prediction(:)
        real(dp), allocatable :: fold_error(:)
        integer :: icost
        integer :: igamma
        integer :: grid

        call validate_tune_inputs(size(x, 1), folds, fold_id)
        if (size(y) /= size(x, 1)) error stop "tune_svm_regression: x/y mismatch"
        if (size(costs) < 1 .or. size(gammas) < 1) error stop "tune_svm_regression: empty parameter grid"
        if (any(costs <= 0.0_dp) .or. any(gammas < 0.0_dp)) error stop "tune_svm_regression: invalid grid value"
        allocate(result%parameter1(size(costs) * size(gammas)), result%parameter2(size(costs) * size(gammas)))
        allocate(result%error(size(result%parameter1)), result%dispersion(size(result%parameter1)))
        allocate(fold_error(folds))
        options = svm_options()
        if (present(base_options)) options = base_options
        grid = 0
        do igamma = 1, size(gammas)
            do icost = 1, size(costs)
                grid = grid + 1
                options%cost = costs(icost)
                options%gamma = gammas(igamma)
                call cv_svm_regression(x, y, folds, options, prediction, fold_error, fold_id)
                result%parameter1(grid) = costs(icost)
                result%parameter2(grid) = gammas(igamma)
                result%error(grid) = sum(fold_error) / real(folds, dp)
                result%dispersion(grid) = sample_sd(fold_error)
            end do
        end do
        call select_best(result)
    end subroutine tune_svm_regression

    subroutine tune_gknn_classification(x, y, k_values, folds, result, method, fold_id)
        real(dp), intent(in) :: x(:, :) !! Classification predictor matrix used to assess each candidate k.
        integer, intent(in) :: y(:) !! Integer class labels corresponding rowwise to x.
        integer, intent(in) :: k_values(:) !! Positive nearest-neighbor counts to evaluate.
        integer, intent(in) :: folds !! Number of cross-validation folds.
        type(tune_result), intent(out) :: result !! Candidate k values in parameter1, errors, dispersion, and best index.
        character(len=*), intent(in), optional :: method !! proxy distance name used by all generalized-kNN fits; default Euclidean.
        integer, intent(in), optional :: fold_id(:) !! Optional one-based fold labels; default uses deterministic cyclic folds.
        real(dp), allocatable :: fold_error(:)
        integer :: i

        call validate_tune_inputs(size(x, 1), folds, fold_id)
        if (size(y) /= size(x, 1)) error stop "tune_gknn_classification: x/y mismatch"
        if (size(k_values) < 1 .or. any(k_values < 1)) error stop "tune_gknn_classification: invalid k grid"
        allocate(result%parameter1(size(k_values)), result%parameter2(0), result%error(size(k_values)))
        allocate(result%dispersion(size(k_values)), fold_error(folds))
        do i = 1, size(k_values)
            call cv_gknn_classification(x, y, k_values(i), folds, fold_error, method, fold_id)
            result%parameter1(i) = real(k_values(i), dp)
            result%error(i) = sum(fold_error) / real(folds, dp)
            result%dispersion(i) = sample_sd(fold_error)
        end do
        call select_best(result)
    end subroutine tune_gknn_classification

    subroutine tune_gknn_regression(x, y, k_values, folds, result, method, fold_id)
        real(dp), intent(in) :: x(:, :) !! Regression predictor matrix used to assess each candidate k.
        real(dp), intent(in) :: y(:) !! Numeric response vector corresponding rowwise to x.
        integer, intent(in) :: k_values(:) !! Positive nearest-neighbor counts to evaluate.
        integer, intent(in) :: folds !! Number of cross-validation folds.
        type(tune_result), intent(out) :: result !! Candidate k values in parameter1, errors, dispersion, and best index.
        character(len=*), intent(in), optional :: method !! proxy distance name used by all generalized-kNN fits; default Euclidean.
        integer, intent(in), optional :: fold_id(:) !! Optional one-based fold labels; default uses deterministic cyclic folds.
        real(dp), allocatable :: fold_error(:)
        integer :: i

        call validate_tune_inputs(size(x, 1), folds, fold_id)
        if (size(y) /= size(x, 1)) error stop "tune_gknn_regression: x/y mismatch"
        if (size(k_values) < 1 .or. any(k_values < 1)) error stop "tune_gknn_regression: invalid k grid"
        allocate(result%parameter1(size(k_values)), result%parameter2(0), result%error(size(k_values)))
        allocate(result%dispersion(size(k_values)), fold_error(folds))
        do i = 1, size(k_values)
            call cv_gknn_regression(x, y, k_values(i), folds, fold_error, method, fold_id)
            result%parameter1(i) = real(k_values(i), dp)
            result%error(i) = sum(fold_error) / real(folds, dp)
            result%dispersion(i) = sample_sd(fold_error)
        end do
        call select_best(result)
    end subroutine tune_gknn_regression

    subroutine cv_svm_classification(x, y, folds, options, prediction, fold_error, fold_id)
        real(dp), intent(in) :: x(:, :) !! Full predictor matrix partitioned into training and validation rows.
        integer, intent(in) :: y(:) !! Full class-label vector.
        integer, intent(in) :: folds !! Number of folds represented by `fold_id` or cyclic assignment.
        type(svm_options), intent(in) :: options !! SVM controls used independently in every fold.
        integer, allocatable, intent(out) :: prediction(:) !! Out-of-fold predicted label for every observation.
        real(dp), intent(out) :: fold_error(:) !! Classification error fraction for each fold; shape must equal folds.
        integer, intent(in), optional :: fold_id(:) !! Optional explicit one-based fold labels.
        real(dp), allocatable :: train_x(:, :)
        real(dp), allocatable :: test_x(:, :)
        integer, allocatable :: train_y(:)
        integer, allocatable :: test_index(:)
        integer, allocatable :: pred(:)
        type(svm_model) :: fit
        integer :: fold

        allocate(prediction(size(y)))
        do fold = 1, folds
            call split_classification(x, y, fold, folds, train_x, train_y, test_x, test_index, fold_id)
            call svm_fit_classification(train_x, train_y, fit, options)
            call svm_predict_classification(fit, test_x, pred)
            prediction(test_index) = pred
            fold_error(fold) = real(count(prediction(test_index) /= y(test_index)), dp) / real(size(test_index), dp)
        end do
    end subroutine cv_svm_classification

    subroutine cv_svm_regression(x, y, folds, options, prediction, fold_error, fold_id)
        real(dp), intent(in) :: x(:, :) !! Full predictor matrix partitioned into training and validation rows.
        real(dp), intent(in) :: y(:) !! Full numeric response vector.
        integer, intent(in) :: folds !! Number of folds represented by `fold_id` or cyclic assignment.
        type(svm_options), intent(in) :: options !! SVM controls used independently in every fold.
        real(dp), allocatable, intent(out) :: prediction(:) !! Out-of-fold numeric prediction for every observation.
        real(dp), intent(out) :: fold_error(:) !! Mean squared error for each fold; shape must equal folds.
        integer, intent(in), optional :: fold_id(:) !! Optional explicit one-based fold labels.
        real(dp), allocatable :: train_x(:, :)
        real(dp), allocatable :: test_x(:, :)
        real(dp), allocatable :: train_y(:)
        real(dp), allocatable :: pred(:)
        integer, allocatable :: test_index(:)
        type(svm_model) :: fit
        integer :: fold

        allocate(prediction(size(y)))
        do fold = 1, folds
            call split_regression(x, y, fold, folds, train_x, train_y, test_x, test_index, fold_id)
            call svm_fit_regression(train_x, train_y, fit, options)
            call svm_predict_regression(fit, test_x, pred)
            prediction(test_index) = pred
            fold_error(fold) = sum((pred - y(test_index))**2) / real(size(test_index), dp)
        end do
    end subroutine cv_svm_regression

    subroutine cv_gknn_classification(x, y, k, folds, fold_error, method, fold_id)
        real(dp), intent(in) :: x(:, :) !! Full classification predictor matrix partitioned foldwise.
        integer, intent(in) :: y(:) !! Full integer class-label vector.
        integer, intent(in) :: k !! Candidate nearest-neighbor count.
        integer, intent(in) :: folds !! Number of cross-validation folds.
        real(dp), intent(out) :: fold_error(:) !! Classification error fraction for each fold.
        character(len=*), intent(in), optional :: method !! proxy distance name used by each fold model.
        integer, intent(in), optional :: fold_id(:) !! Optional explicit one-based fold labels.
        real(dp), allocatable :: train_x(:, :)
        real(dp), allocatable :: test_x(:, :)
        integer, allocatable :: train_y(:)
        integer, allocatable :: test_index(:)
        integer, allocatable :: pred(:)
        type(gknn_model) :: fit
        integer :: fold

        do fold = 1, folds
            call split_classification(x, y, fold, folds, train_x, train_y, test_x, test_index, fold_id)
            if (k > size(train_y)) error stop "tune_gknn_classification: k exceeds fold training size"
            call gknn_fit_classification(train_x, train_y, fit, k=k, method=optional_method(method))
            call gknn_predict_classification(fit, test_x, pred)
            fold_error(fold) = real(count(pred /= y(test_index)), dp) / real(size(test_index), dp)
        end do
    end subroutine cv_gknn_classification

    subroutine cv_gknn_regression(x, y, k, folds, fold_error, method, fold_id)
        real(dp), intent(in) :: x(:, :) !! Full regression predictor matrix partitioned foldwise.
        real(dp), intent(in) :: y(:) !! Full numeric response vector.
        integer, intent(in) :: k !! Candidate nearest-neighbor count.
        integer, intent(in) :: folds !! Number of cross-validation folds.
        real(dp), intent(out) :: fold_error(:) !! Mean squared prediction error for each fold.
        character(len=*), intent(in), optional :: method !! proxy distance name used by each fold model.
        integer, intent(in), optional :: fold_id(:) !! Optional explicit one-based fold labels.
        real(dp), allocatable :: train_x(:, :)
        real(dp), allocatable :: test_x(:, :)
        real(dp), allocatable :: train_y(:)
        real(dp), allocatable :: pred(:)
        integer, allocatable :: test_index(:)
        type(gknn_model) :: fit
        integer :: fold

        do fold = 1, folds
            call split_regression(x, y, fold, folds, train_x, train_y, test_x, test_index, fold_id)
            if (k > size(train_y)) error stop "tune_gknn_regression: k exceeds fold training size"
            call gknn_fit_regression(train_x, train_y, fit, k=k, method=optional_method(method))
            call gknn_predict_regression(fit, test_x, pred)
            fold_error(fold) = sum((pred - y(test_index))**2) / real(size(test_index), dp)
        end do
    end subroutine cv_gknn_regression

    subroutine split_classification(x, y, fold, folds, train_x, train_y, test_x, test_index, fold_id)
        real(dp), intent(in) :: x(:, :) !! Full predictor matrix to split.
        integer, intent(in) :: y(:) !! Full integer response vector to split in parallel with x.
        integer, intent(in) :: fold !! One-based validation fold selected for this split.
        integer, intent(in) :: folds !! Total number of folds used by cyclic fallback assignment.
        real(dp), allocatable, intent(out) :: train_x(:, :) !! Predictor rows outside the selected fold.
        integer, allocatable, intent(out) :: train_y(:) !! Class labels outside the selected fold.
        real(dp), allocatable, intent(out) :: test_x(:, :) !! Predictor rows inside the selected fold.
        integer, allocatable, intent(out) :: test_index(:) !! Original row indices in the selected validation fold.
        integer, intent(in), optional :: fold_id(:) !! Optional explicit one-based fold labels.
        integer :: i
        integer :: itrain
        integer :: itest

        call allocate_split(x, fold, folds, train_x, test_x, test_index, fold_id)
        allocate(train_y(size(train_x, 1)))
        itrain = 0
        itest = 0
        do i = 1, size(y)
            if (assigned_fold(i, folds, fold_id) == fold) then
                itest = itest + 1
            else
                itrain = itrain + 1
                train_y(itrain) = y(i)
            end if
        end do
    end subroutine split_classification

    subroutine split_regression(x, y, fold, folds, train_x, train_y, test_x, test_index, fold_id)
        real(dp), intent(in) :: x(:, :) !! Full predictor matrix to split.
        real(dp), intent(in) :: y(:) !! Full numeric response vector to split in parallel with x.
        integer, intent(in) :: fold !! One-based validation fold selected for this split.
        integer, intent(in) :: folds !! Total number of folds used by cyclic fallback assignment.
        real(dp), allocatable, intent(out) :: train_x(:, :) !! Predictor rows outside the selected fold.
        real(dp), allocatable, intent(out) :: train_y(:) !! Numeric responses outside the selected fold.
        real(dp), allocatable, intent(out) :: test_x(:, :) !! Predictor rows inside the selected fold.
        integer, allocatable, intent(out) :: test_index(:) !! Original row indices in the selected validation fold.
        integer, intent(in), optional :: fold_id(:) !! Optional explicit one-based fold labels.
        integer :: i
        integer :: itrain

        call allocate_split(x, fold, folds, train_x, test_x, test_index, fold_id)
        allocate(train_y(size(train_x, 1)))
        itrain = 0
        do i = 1, size(y)
            if (assigned_fold(i, folds, fold_id) /= fold) then
                itrain = itrain + 1
                train_y(itrain) = y(i)
            end if
        end do
    end subroutine split_regression

    subroutine allocate_split(x, fold, folds, train_x, test_x, test_index, fold_id)
        real(dp), intent(in) :: x(:, :) !! Full predictor matrix whose rows are partitioned.
        integer, intent(in) :: fold !! One-based validation fold selected for this split.
        integer, intent(in) :: folds !! Total number of folds used by cyclic fallback assignment.
        real(dp), allocatable, intent(out) :: train_x(:, :) !! Predictor rows not belonging to the selected fold.
        real(dp), allocatable, intent(out) :: test_x(:, :) !! Predictor rows belonging to the selected fold.
        integer, allocatable, intent(out) :: test_index(:) !! Original row indices of validation rows.
        integer, intent(in), optional :: fold_id(:) !! Optional explicit one-based fold labels.
        integer :: ntest
        integer :: i
        integer :: itrain
        integer :: itest

        ntest = 0
        do i = 1, size(x, 1)
            if (assigned_fold(i, folds, fold_id) == fold) ntest = ntest + 1
        end do
        if (ntest < 1 .or. ntest >= size(x, 1)) error stop "allocate_split: empty training or validation set"
        allocate(train_x(size(x, 1) - ntest, size(x, 2)), test_x(ntest, size(x, 2)), test_index(ntest))
        itrain = 0
        itest = 0
        do i = 1, size(x, 1)
            if (assigned_fold(i, folds, fold_id) == fold) then
                itest = itest + 1
                test_x(itest, :) = x(i, :)
                test_index(itest) = i
            else
                itrain = itrain + 1
                train_x(itrain, :) = x(i, :)
            end if
        end do
    end subroutine allocate_split

    pure function assigned_fold(row, folds, fold_id) result(fold)
        integer, intent(in) :: row !! One-based observation index whose fold label is requested.
        integer, intent(in) :: folds !! Positive fold count used for cyclic fallback assignment.
        integer, intent(in), optional :: fold_id(:) !! Optional explicit fold labels overriding cyclic assignment.
        integer :: fold

        if (present(fold_id)) then
            fold = fold_id(row)
        else
            fold = modulo(row - 1, folds) + 1
        end if
    end function assigned_fold

    subroutine validate_tune_inputs(n, folds, fold_id)
        integer, intent(in) :: n !! Number of observations available for resampling.
        integer, intent(in) :: folds !! Requested fold count, which must lie between two and n.
        integer, intent(in), optional :: fold_id(:) !! Optional explicit fold labels with length n and values 1..folds.

        if (folds < 2 .or. folds > n) error stop "tune: invalid fold count"
        if (present(fold_id)) then
            if (size(fold_id) /= n) error stop "tune: fold_id length mismatch"
            if (any(fold_id < 1) .or. any(fold_id > folds)) error stop "tune: fold_id value out of range"
        end if
    end subroutine validate_tune_inputs

    subroutine select_best(result)
        type(tune_result), intent(inout) :: result !! Tuning table updated with its minimum-error row and performance value.

        result%best_index = minloc(result%error, dim=1)
        result%best_performance = result%error(result%best_index)
    end subroutine select_best

    pure function sample_sd(x) result(value)
        real(dp), intent(in) :: x(:) !! Fold-level performance values whose sample standard deviation is returned.
        real(dp) :: value
        real(dp) :: mean_value

        if (size(x) <= 1) then
            value = 0.0_dp
            return
        end if
        mean_value = sum(x) / real(size(x), dp)
        value = sqrt(sum((x - mean_value)**2) / real(size(x) - 1, dp))
    end function sample_sd

    function optional_method(method) result(value)
        character(len=*), intent(in), optional :: method !! Optional proxy method name; absence maps to the e1071 Euclidean default.
        character(len=32) :: value

        value = "euclidean"
        if (present(method)) value = trim(adjustl(method))
    end function optional_method

end module e1071_tune
