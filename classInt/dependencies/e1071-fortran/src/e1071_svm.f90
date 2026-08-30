module e1071_svm
    use e1071_kinds, only: dp
    use e1071_rng, only: rng_state, rng_seed, rng_integer
    use e1071_svm_solver, only: solve_dual, initialize_nu_svc_alpha
    implicit none
    private

    integer, parameter, public :: svm_c_classification = 0
    integer, parameter, public :: svm_nu_classification = 1
    integer, parameter, public :: svm_one_classification = 2
    integer, parameter, public :: svm_eps_regression = 3
    integer, parameter, public :: svm_nu_regression = 4
    integer, parameter, public :: svm_linear = 0
    integer, parameter, public :: svm_polynomial = 1
    integer, parameter, public :: svm_radial = 2
    integer, parameter, public :: svm_sigmoid = 3

    type, public :: svm_options
        integer :: svm_type = svm_c_classification
        integer :: kernel = svm_radial
        integer :: degree = 3
        real(dp) :: gamma = -1.0_dp
        real(dp) :: coef0 = 0.0_dp
        real(dp) :: cost = 1.0_dp
        real(dp) :: nu = 0.5_dp
        real(dp) :: tolerance = 0.001_dp
        real(dp) :: epsilon = 0.1_dp
        logical :: scale = .true.
        logical :: probability = .false.
        integer :: max_iterations = 0
        integer :: probability_seed = 1
    end type svm_options

    type, public :: svm_pair_model
        real(dp), allocatable :: support_vectors(:, :)
        real(dp), allocatable :: coefficients(:)
        real(dp) :: rho = 0.0_dp
        real(dp) :: prob_a = 0.0_dp
        real(dp) :: prob_b = 0.0_dp
        integer :: positive_label = 1
        integer :: negative_label = -1
    end type svm_pair_model

    type, public :: svm_model
        type(svm_options) :: options
        type(svm_pair_model), allocatable :: pair(:)
        integer, allocatable :: class_labels(:)
        real(dp), allocatable :: center(:)
        real(dp), allocatable :: scale_value(:)
        logical, allocatable :: scale_mask(:)
        real(dp) :: y_center = 0.0_dp
        real(dp) :: y_scale = 1.0_dp
        integer :: n_features = 0
        logical :: probability_fitted = .false.
        real(dp) :: svr_probability_sigma = 0.0_dp
    end type svm_model

    public :: svm_fit_classification, svm_fit_regression, svm_fit_one_class
    public :: svm_predict_classification, svm_predict_regression, svm_predict_one_class
    public :: svm_decision_values, svm_linear_coefficients
    public :: svm_cross_validate_classification, svm_cross_validate_regression

contains

    subroutine svm_fit_classification(x, y, model, options, class_weight_labels, class_weights, scale_mask)
        real(dp), intent(in) :: x(:, :) !! Training observation-by-variable matrix for C- or nu-classification.
        integer, intent(in) :: y(:) !! Integer class labels; at least two distinct values are required.
        type(svm_model), intent(out) :: model !! Fitted one-vs-one SVM model with support vectors stored in scaled predictor space.
        type(svm_options), intent(in), optional :: options !! SVM controls; type must be C- or nu-classification and defaults to
        !! C-SVC.
        integer, intent(in), optional :: class_weight_labels(:) !! Original class labels receiving multiplicative C weights for
        !! C-SVC.
        real(dp), intent(in), optional :: class_weights(:) !! Positive C multipliers paired with class_weight_labels.
        logical, intent(in), optional :: scale_mask(:) !! Per-column scaling selector; by default every column follows
        !! options%scale.
        type(svm_options) :: opt
        real(dp), allocatable :: sx(:, :)
        integer, allocatable :: index(:)
        integer, allocatable :: binary_y(:)
        real(dp), allocatable :: pair_x(:, :)
        real(dp) :: cp
        real(dp) :: cn
        integer :: nclass
        integer :: npair
        integer :: a
        integer :: b
        integer :: p
        integer :: i
        integer :: n

        if (size(y) /= size(x, 1)) error stop "svm_fit_classification: x/y row mismatch"
        opt = svm_options()
        if (present(options)) opt = options
        if (opt%svm_type /= svm_c_classification .and. opt%svm_type /= svm_nu_classification) then
            error stop "svm_fit_classification: invalid svm_type"
        end if
        call validate_options(opt, size(x, 2))
        call unique_integer_labels(y, model%class_labels)
        nclass = size(model%class_labels)
        if (nclass < 2) error stop "svm_fit_classification: at least two classes are required"
        call prepare_predictors_fit(x, opt%scale, scale_mask, sx, model%center, model%scale_value, model%scale_mask)
        model%n_features = size(x, 2)
        model%options = opt
        if (model%options%gamma < 0.0_dp) model%options%gamma = 1.0_dp / real(size(x, 2), dp)
        npair = nclass * (nclass - 1) / 2
        allocate(model%pair(npair))
        p = 0
        do a = 1, nclass - 1
            do b = a + 1, nclass
                p = p + 1
                n = count(y == model%class_labels(a)) + count(y == model%class_labels(b))
                allocate(index(n), binary_y(n), pair_x(n, size(x, 2)))
                i = 0
                call collect_pair_rows(sx, y, model%class_labels(a), model%class_labels(b), pair_x, binary_y, index, i)
                cp = opt%cost * class_weight(model%class_labels(a), class_weight_labels, class_weights)
                cn = opt%cost * class_weight(model%class_labels(b), class_weight_labels, class_weights)
                call train_binary_classifier(pair_x, binary_y, model%pair(p), model%options, cp, cn)
                model%pair(p)%positive_label = model%class_labels(a)
                model%pair(p)%negative_label = model%class_labels(b)
                deallocate(index, binary_y, pair_x)
            end do
        end do
        model%probability_fitted = opt%probability
    end subroutine svm_fit_classification

    subroutine svm_fit_regression(x, y, model, options, scale_mask)
        real(dp), intent(in) :: x(:, :) !! Training observation-by-variable matrix for epsilon- or nu-SVR.
        real(dp), intent(in) :: y(:) !! Numeric response vector with one value per training row.
        type(svm_model), intent(out) :: model !! Fitted regression SVM with one support-vector expansion.
        type(svm_options), intent(in), optional :: options !! SVM controls; type must be eps-regression or nu-regression.
        logical, intent(in), optional :: scale_mask(:) !! Per-column scaling selector; default follows options%scale for all
        !! columns.
        type(svm_options) :: opt
        real(dp), allocatable :: sx(:, :)
        real(dp), allocatable :: sy(:)
        real(dp) :: ss

        if (size(y) /= size(x, 1)) error stop "svm_fit_regression: x/y row mismatch"
        opt = svm_options(svm_type=svm_eps_regression)
        if (present(options)) opt = options
        if (opt%svm_type /= svm_eps_regression .and. opt%svm_type /= svm_nu_regression) then
            error stop "svm_fit_regression: invalid svm_type"
        end if
        call validate_options(opt, size(x, 2))
        call prepare_predictors_fit(x, opt%scale, scale_mask, sx, model%center, model%scale_value, model%scale_mask)
        model%n_features = size(x, 2)
        model%options = opt
        if (model%options%gamma < 0.0_dp) model%options%gamma = 1.0_dp / real(size(x, 2), dp)
        sy = y
        model%y_center = 0.0_dp
        model%y_scale = 1.0_dp
        if (any(model%scale_mask)) then
            model%y_center = sum(y) / real(size(y), dp)
            ss = sum((y - model%y_center)**2)
            if (size(y) > 1 .and. ss > 0.0_dp) then
                model%y_scale = sqrt(ss / real(size(y) - 1, dp))
                sy = (y - model%y_center) / model%y_scale
            end if
        end if
        allocate(model%pair(1))
        call train_regressor(sx, sy, model%pair(1), model%options)
        if (model%options%probability) then
            call svr_probability(sx, sy, model%options, model%svr_probability_sigma)
            model%probability_fitted = .true.
        end if
    end subroutine svm_fit_regression

    subroutine svm_fit_one_class(x, model, options, scale_mask)
        real(dp), intent(in) :: x(:, :) !! Training observations defining the target one-class support.
        type(svm_model), intent(out) :: model !! Fitted one-class SVM model with support-vector expansion and threshold rho.
        type(svm_options), intent(in), optional :: options !! SVM controls; svm_type is forced/required to one-classification.
        logical, intent(in), optional :: scale_mask(:) !! Per-column scaling selector; default follows options%scale for all
        !! columns.
        type(svm_options) :: opt
        real(dp), allocatable :: sx(:, :)

        opt = svm_options(svm_type=svm_one_classification)
        if (present(options)) opt = options
        if (opt%svm_type /= svm_one_classification) error stop "svm_fit_one_class: invalid svm_type"
        call validate_options(opt, size(x, 2))
        call prepare_predictors_fit(x, opt%scale, scale_mask, sx, model%center, model%scale_value, model%scale_mask)
        model%n_features = size(x, 2)
        model%options = opt
        if (model%options%gamma < 0.0_dp) model%options%gamma = 1.0_dp / real(size(x, 2), dp)
        allocate(model%pair(1))
        call train_one_class_model(sx, model%pair(1), model%options)
    end subroutine svm_fit_one_class

    subroutine svm_predict_classification(model, x, class, decision, probability)
        type(svm_model), intent(in) :: model !! Fitted one-vs-one classification SVM.
        real(dp), intent(in) :: x(:, :) !! New observation-by-variable matrix in original unscaled units.
        integer, allocatable, intent(out) :: class(:) !! Predicted original integer class label for each row.
        real(dp), allocatable, intent(out), optional :: decision(:, :) !! Optional pairwise decision values by row and pair order.
        real(dp), allocatable, intent(out), optional :: probability(:, :) !! Optional coupled class probabilities when
        !! probability fitting was requested.
        real(dp), allocatable :: sx(:, :)
        real(dp), allocatable :: dec(:, :)
        integer, allocatable :: votes(:)
        real(dp), allocatable :: r(:, :)
        real(dp), allocatable :: pclass(:)
        real(dp) :: pij
        integer :: nclass
        integer :: npair
        integer :: row
        integer :: pair_index
        integer :: a
        integer :: b
        integer :: winner

        if (model%options%svm_type /= svm_c_classification .and. model%options%svm_type /= svm_nu_classification) then
            error stop "svm_predict_classification: non-classification model supplied"
        end if
        call prepare_predictors_apply(model, x, sx)
        nclass = size(model%class_labels)
        npair = size(model%pair)
        allocate(class(size(x, 1)), dec(size(x, 1), npair), votes(nclass))
        if (present(probability)) then
            if (.not. model%probability_fitted) error stop "svm_predict_classification: probability model was not fitted"
            allocate(probability(size(x, 1), nclass), r(nclass, nclass), pclass(nclass))
        end if
        do row = 1, size(x, 1)
            votes = 0
            if (present(probability)) then
                r = 0.0_dp
                do a = 1, nclass
                    r(a, a) = 0.0_dp
                end do
            end if
            pair_index = 0
            do a = 1, nclass - 1
                do b = a + 1, nclass
                    pair_index = pair_index + 1
                    dec(row, pair_index) = pair_decision(model%pair(pair_index), sx(row, :), model%options)
                    if (dec(row, pair_index) > 0.0_dp) then
                        votes(a) = votes(a) + 1
                    else
                        votes(b) = votes(b) + 1
                    end if
                    if (present(probability)) then
                        pij = sigmoid_predict(dec(row, pair_index), model%pair(pair_index)%prob_a, model%pair(pair_index)%prob_b)
                        pij = max(1.0e-7_dp, min(1.0_dp - 1.0e-7_dp, pij))
                        r(a, b) = pij
                        r(b, a) = 1.0_dp - pij
                    end if
                end do
            end do
            winner = maxloc(votes, dim=1)
            class(row) = model%class_labels(winner)
            if (present(probability)) then
                call multiclass_probability(r, pclass)
                probability(row, :) = pclass
            end if
        end do
        if (present(decision)) then
            allocate(decision(size(dec, 1), size(dec, 2)))
            decision = dec
        end if
    end subroutine svm_predict_classification

    subroutine svm_predict_regression(model, x, prediction, decision)
        type(svm_model), intent(in) :: model !! Fitted epsilon- or nu-regression SVM.
        real(dp), intent(in) :: x(:, :) !! New observation-by-variable matrix in original unscaled units.
        real(dp), allocatable, intent(out) :: prediction(:) !! Predicted numeric responses transformed back to original response
        !! units.
        real(dp), allocatable, intent(out), optional :: decision(:) !! Optional raw scaled-space decision values before response
        !! unscaling.
        real(dp), allocatable :: sx(:, :)
        real(dp) :: raw
        integer :: i

        if (model%options%svm_type /= svm_eps_regression .and. model%options%svm_type /= svm_nu_regression) then
            error stop "svm_predict_regression: non-regression model supplied"
        end if
        call prepare_predictors_apply(model, x, sx)
        allocate(prediction(size(x, 1)))
        if (present(decision)) allocate(decision(size(x, 1)))
        do i = 1, size(x, 1)
            raw = pair_decision(model%pair(1), sx(i, :), model%options)
            prediction(i) = raw * model%y_scale + model%y_center
            if (present(decision)) decision(i) = raw
        end do
    end subroutine svm_predict_regression

    subroutine svm_predict_one_class(model, x, accepted, decision)
        type(svm_model), intent(in) :: model !! Fitted one-class SVM model.
        real(dp), intent(in) :: x(:, :) !! New observation-by-variable matrix in original unscaled units.
        logical, allocatable, intent(out) :: accepted(:) !! True where the one-class decision function is nonnegative.
        real(dp), allocatable, intent(out), optional :: decision(:) !! Optional raw signed distance-like decision values.
        real(dp), allocatable :: sx(:, :)
        real(dp) :: raw
        integer :: i

        if (model%options%svm_type /= svm_one_classification) error stop "svm_predict_one_class: wrong model type"
        call prepare_predictors_apply(model, x, sx)
        allocate(accepted(size(x, 1)))
        if (present(decision)) allocate(decision(size(x, 1)))
        do i = 1, size(x, 1)
            raw = pair_decision(model%pair(1), sx(i, :), model%options)
            accepted(i) = raw >= 0.0_dp
            if (present(decision)) decision(i) = raw
        end do
    end subroutine svm_predict_one_class

    subroutine svm_decision_values(model, x, values)
        type(svm_model), intent(in) :: model !! Fitted SVM whose raw support-vector decision functions are requested.
        real(dp), intent(in) :: x(:, :) !! New observation-by-variable matrix in original unscaled units.
        real(dp), allocatable, intent(out) :: values(:, :) !! Raw decision matrix; one column per pair, or one column for
        !! regression/one-class.
        real(dp), allocatable :: sx(:, :)
        integer :: i
        integer :: p

        call prepare_predictors_apply(model, x, sx)
        allocate(values(size(x, 1), size(model%pair)))
        do p = 1, size(model%pair)
            do i = 1, size(x, 1)
                values(i, p) = pair_decision(model%pair(p), sx(i, :), model%options)
            end do
        end do
    end subroutine svm_decision_values

    subroutine svm_linear_coefficients(model, intercept, coefficients)
        type(svm_model), intent(in) :: model !! Linear-kernel binary/regression/one-class model whose primal coefficients are
        !! requested.
        real(dp), intent(out) :: intercept !! Decision-function intercept in the model's scaled predictor space, equal to -rho.
        real(dp), allocatable, intent(out) :: coefficients(:) !! Linear primal weight vector in scaled predictor units.
        integer :: i

        if (model%options%kernel /= svm_linear) error stop "svm_linear_coefficients: model kernel is not linear"
        if (size(model%pair) /= 1) error stop "svm_linear_coefficients: only binary or single-expansion models are supported"
        allocate(coefficients(model%n_features))
        coefficients = 0.0_dp
        do i = 1, size(model%pair(1)%coefficients)
            coefficients = coefficients + model%pair(1)%coefficients(i) * model%pair(1)%support_vectors(i, :)
        end do
        intercept = -model%pair(1)%rho
    end subroutine svm_linear_coefficients

    subroutine svm_cross_validate_classification(x, y, folds, prediction, options, fold_id)
        real(dp), intent(in) :: x(:, :) !! Full classification predictor matrix divided into training/validation folds.
        integer, intent(in) :: y(:) !! Full integer class-label vector corresponding rowwise to x.
        integer, intent(in) :: folds !! Number of cross-validation folds, between two and the number of observations.
        integer, allocatable, intent(out) :: prediction(:) !! Out-of-fold predicted class label for every row.
        type(svm_options), intent(in), optional :: options !! Classification SVM controls reused independently within every fold.
        integer, intent(in), optional :: fold_id(:) !! Optional one-based fold label for each row; otherwise cyclic folds are used.
        real(dp), allocatable :: train_x(:, :)
        real(dp), allocatable :: test_x(:, :)
        integer, allocatable :: train_y(:)
        integer, allocatable :: pred(:)
        integer, allocatable :: test_index(:)
        type(svm_model) :: fit
        type(svm_options) :: opt
        integer :: fold
        integer :: ntrain
        integer :: ntest
        integer :: i
        integer :: ti
        integer :: vi

        if (size(y) /= size(x, 1)) error stop "svm_cross_validate_classification: x/y mismatch"
        if (folds < 2 .or. folds > size(x, 1)) error stop "svm_cross_validate_classification: invalid folds"
        if (present(fold_id)) then
            if (size(fold_id) /= size(y)) error stop "svm_cross_validate_classification: fold_id length mismatch"
            if (any(fold_id < 1) .or. any(fold_id > folds)) error stop "svm_cross_validate_classification: invalid fold_id"
        end if
        opt = svm_options()
        if (present(options)) opt = options
        allocate(prediction(size(y)))
        do fold = 1, folds
            if (present(fold_id)) then
                ntest = count(fold_id == fold)
            else
                ntest = count([(modulo(i - 1, folds) + 1 == fold, i = 1, size(y))])
            end if
            ntrain = size(y) - ntest
            if (ntest == 0 .or. ntrain == 0) error stop "svm_cross_validate_classification: empty fold"
            allocate(train_x(ntrain, size(x, 2)), train_y(ntrain), test_x(ntest, size(x, 2)), test_index(ntest))
            ti = 0
            vi = 0
            do i = 1, size(y)
                if (row_fold(i, folds, fold_id) == fold) then
                    vi = vi + 1
                    test_x(vi, :) = x(i, :)
                    test_index(vi) = i
                else
                    ti = ti + 1
                    train_x(ti, :) = x(i, :)
                    train_y(ti) = y(i)
                end if
            end do
            call svm_fit_classification(train_x, train_y, fit, opt)
            call svm_predict_classification(fit, test_x, pred)
            do i = 1, ntest
                prediction(test_index(i)) = pred(i)
            end do
            deallocate(train_x, train_y, test_x, test_index, pred)
        end do
    end subroutine svm_cross_validate_classification

    subroutine svm_cross_validate_regression(x, y, folds, prediction, options, fold_id)
        real(dp), intent(in) :: x(:, :) !! Full regression predictor matrix divided into training/validation folds.
        real(dp), intent(in) :: y(:) !! Full numeric response vector corresponding rowwise to x.
        integer, intent(in) :: folds !! Number of cross-validation folds, between two and the number of observations.
        real(dp), allocatable, intent(out) :: prediction(:) !! Out-of-fold numeric prediction for every observation.
        type(svm_options), intent(in), optional :: options !! Regression SVM controls reused independently within every fold.
        integer, intent(in), optional :: fold_id(:) !! Optional one-based fold label for each row; otherwise cyclic folds are used.
        real(dp), allocatable :: train_x(:, :)
        real(dp), allocatable :: test_x(:, :)
        real(dp), allocatable :: train_y(:)
        real(dp), allocatable :: pred(:)
        integer, allocatable :: test_index(:)
        type(svm_model) :: fit
        type(svm_options) :: opt
        integer :: fold
        integer :: ntrain
        integer :: ntest
        integer :: i
        integer :: ti
        integer :: vi

        if (size(y) /= size(x, 1)) error stop "svm_cross_validate_regression: x/y mismatch"
        if (folds < 2 .or. folds > size(x, 1)) error stop "svm_cross_validate_regression: invalid folds"
        if (present(fold_id)) then
            if (size(fold_id) /= size(y)) error stop "svm_cross_validate_regression: fold_id length mismatch"
            if (any(fold_id < 1) .or. any(fold_id > folds)) error stop "svm_cross_validate_regression: invalid fold_id"
        end if
        opt = svm_options(svm_type=svm_eps_regression)
        if (present(options)) opt = options
        allocate(prediction(size(y)))
        do fold = 1, folds
            if (present(fold_id)) then
                ntest = count(fold_id == fold)
            else
                ntest = count([(modulo(i - 1, folds) + 1 == fold, i = 1, size(y))])
            end if
            ntrain = size(y) - ntest
            if (ntest == 0 .or. ntrain == 0) error stop "svm_cross_validate_regression: empty fold"
            allocate(train_x(ntrain, size(x, 2)), train_y(ntrain), test_x(ntest, size(x, 2)), test_index(ntest))
            ti = 0
            vi = 0
            do i = 1, size(y)
                if (row_fold(i, folds, fold_id) == fold) then
                    vi = vi + 1
                    test_x(vi, :) = x(i, :)
                    test_index(vi) = i
                else
                    ti = ti + 1
                    train_x(ti, :) = x(i, :)
                    train_y(ti) = y(i)
                end if
            end do
            call svm_fit_regression(train_x, train_y, fit, opt)
            call svm_predict_regression(fit, test_x, pred)
            do i = 1, ntest
                prediction(test_index(i)) = pred(i)
            end do
            deallocate(train_x, train_y, test_x, test_index, pred)
        end do
    end subroutine svm_cross_validate_regression

    recursive subroutine train_binary_classifier(x, y, pair, options, cp, cn)
        real(dp), intent(in) :: x(:, :) !! Two-class predictor matrix already transformed into fitted scaling space.
        integer, intent(in) :: y(:) !! Binary signs encoded +1 for the pair's positive class and -1 for the negative class.
        type(svm_pair_model), intent(out) :: pair !! Fitted support-vector expansion and optional sigmoid calibration.
        type(svm_options), intent(in) :: options !! Validated classification options with resolved gamma.
        real(dp), intent(in) :: cp !! Positive-class upper bound C after class weighting.
        real(dp), intent(in) :: cn !! Negative-class upper bound C after class weighting.
        real(dp), allocatable :: kernel(:, :)
        real(dp), allocatable :: q(:, :)
        real(dp), allocatable :: alpha(:)
        real(dp), allocatable :: p(:)
        integer, allocatable :: sign_y(:)
        real(dp) :: rho
        real(dp) :: r
        integer :: i
        integer :: nsv

        if (count(y > 0) == 0 .or. count(y < 0) == 0) error stop "train_binary_classifier: both signs are required"
        call kernel_matrix(x, x, options, kernel)
        allocate(q(size(y), size(y)), alpha(size(y)), p(size(y)), sign_y(size(y)))
        sign_y = merge(1, -1, y > 0)
        do i = 1, size(y)
            q(i, :) = real(sign_y(i), dp) * real(sign_y, dp) * kernel(i, :)
        end do
        if (options%svm_type == svm_c_classification) then
            alpha = 0.0_dp
            p = -1.0_dp
            call solve_dual(q, p, sign_y, alpha, cp, cn, options%tolerance, .false., options%max_iterations, rho, r)
            alpha = alpha * real(sign_y, dp)
        else
            call initialize_nu_svc_alpha(sign_y, options%nu, alpha)
            p = 0.0_dp
            call solve_dual(q, p, sign_y, alpha, 1.0_dp, 1.0_dp, options%tolerance, .true., options%max_iterations, rho, r)
            if (r <= 0.0_dp) error stop "train_binary_classifier: invalid nu-SVC scaling factor"
            alpha = alpha * real(sign_y, dp) / r
            rho = rho / r
        end if
        nsv = count(abs(alpha) > 1.0e-10_dp)
        allocate(pair%support_vectors(nsv, size(x, 2)), pair%coefficients(nsv))
        nsv = 0
        do i = 1, size(alpha)
            if (abs(alpha(i)) <= 1.0e-10_dp) cycle
            nsv = nsv + 1
            pair%support_vectors(nsv, :) = x(i, :)
            pair%coefficients(nsv) = alpha(i)
        end do
        pair%rho = rho
        if (options%probability) then
            call binary_svc_probability(x, sign_y, options, cp, cn, pair%prob_a, pair%prob_b)
        end if
    end subroutine train_binary_classifier

    recursive subroutine train_regressor(x, y, pair, options)
        real(dp), intent(in) :: x(:, :) !! Scaled training predictor matrix for support-vector regression.
        real(dp), intent(in) :: y(:) !! Training response in the model's scaled response units.
        type(svm_pair_model), intent(out) :: pair !! Fitted SVR support-vector expansion and rho.
        type(svm_options), intent(in) :: options !! Validated epsilon- or nu-SVR options with resolved gamma.
        real(dp), allocatable :: kernel(:, :)
        real(dp), allocatable :: q(:, :)
        real(dp), allocatable :: alpha2(:)
        real(dp), allocatable :: alpha(:)
        real(dp), allocatable :: p(:)
        integer, allocatable :: signs(:)
        real(dp) :: rho
        real(dp) :: r
        real(dp) :: sum_alpha
        integer :: n
        integer :: i
        integer :: j
        integer :: ii
        integer :: jj
        integer :: nsv

        n = size(y)
        call kernel_matrix(x, x, options, kernel)
        allocate(q(2 * n, 2 * n), alpha2(2 * n), alpha(n), p(2 * n), signs(2 * n))
        signs(:n) = 1
        signs(n + 1:) = -1
        do i = 1, 2 * n
            ii = modulo(i - 1, n) + 1
            do j = 1, 2 * n
                jj = modulo(j - 1, n) + 1
                q(i, j) = real(signs(i) * signs(j), dp) * kernel(ii, jj)
            end do
        end do
        if (options%svm_type == svm_eps_regression) then
            alpha2 = 0.0_dp
            p(:n) = options%epsilon - y
            p(n + 1:) = options%epsilon + y
            call solve_dual(q, p, signs, alpha2, options%cost, options%cost, options%tolerance, .false., &
                            options%max_iterations, rho, r)
        else
            sum_alpha = options%cost * options%nu * real(n, dp) / 2.0_dp
            do i = 1, n
                alpha2(i) = min(sum_alpha, options%cost)
                alpha2(i + n) = alpha2(i)
                sum_alpha = sum_alpha - alpha2(i)
            end do
            p(:n) = -y
            p(n + 1:) = y
            call solve_dual(q, p, signs, alpha2, options%cost, options%cost, options%tolerance, .true., &
                            options%max_iterations, rho, r)
        end if
        alpha = alpha2(:n) - alpha2(n + 1:)
        nsv = count(abs(alpha) > 1.0e-10_dp)
        allocate(pair%support_vectors(nsv, size(x, 2)), pair%coefficients(nsv))
        nsv = 0
        do i = 1, n
            if (abs(alpha(i)) <= 1.0e-10_dp) cycle
            nsv = nsv + 1
            pair%support_vectors(nsv, :) = x(i, :)
            pair%coefficients(nsv) = alpha(i)
        end do
        pair%rho = rho
    end subroutine train_regressor

    subroutine train_one_class_model(x, pair, options)
        real(dp), intent(in) :: x(:, :) !! Scaled observations defining the target one-class distribution.
        type(svm_pair_model), intent(out) :: pair !! Fitted one-class support-vector coefficients and threshold rho.
        type(svm_options), intent(in) :: options !! Validated one-class SVM controls with resolved gamma and nu.
        real(dp), allocatable :: kernel(:, :)
        real(dp), allocatable :: alpha(:)
        real(dp), allocatable :: p(:)
        integer, allocatable :: signs(:)
        real(dp) :: rho
        real(dp) :: r
        real(dp) :: remaining
        integer :: n
        integer :: i
        integer :: full
        integer :: nsv

        n = size(x, 1)
        call kernel_matrix(x, x, options, kernel)
        allocate(alpha(n), p(n), signs(n))
        alpha = 0.0_dp
        p = 0.0_dp
        signs = 1
        full = int(options%nu * real(n, dp))
        if (full > 0) alpha(:min(full, n)) = 1.0_dp
        remaining = options%nu * real(n, dp) - real(full, dp)
        if (full < n .and. remaining > 0.0_dp) alpha(full + 1) = remaining
        call solve_dual(kernel, p, signs, alpha, 1.0_dp, 1.0_dp, options%tolerance, .false., &
                        options%max_iterations, rho, r)
        nsv = count(alpha > 1.0e-10_dp)
        allocate(pair%support_vectors(nsv, size(x, 2)), pair%coefficients(nsv))
        nsv = 0
        do i = 1, n
            if (alpha(i) <= 1.0e-10_dp) cycle
            nsv = nsv + 1
            pair%support_vectors(nsv, :) = x(i, :)
            pair%coefficients(nsv) = alpha(i)
        end do
        pair%rho = rho
    end subroutine train_one_class_model

    subroutine kernel_matrix(x, y, options, kernel)
        real(dp), intent(in) :: x(:, :) !! First observation-by-variable matrix defining kernel rows.
        real(dp), intent(in) :: y(:, :) !! Second observation-by-variable matrix defining kernel columns.
        type(svm_options), intent(in) :: options !! Kernel code and parameters; gamma must already be resolved.
        real(dp), allocatable, intent(out) :: kernel(:, :) !! Allocated pairwise kernel matrix with shape nrow(x) by nrow(y).
        integer :: i
        integer :: j

        if (size(x, 2) /= size(y, 2)) error stop "kernel_matrix: variable count mismatch"
        allocate(kernel(size(x, 1), size(y, 1)))
        do i = 1, size(x, 1)
            do j = 1, size(y, 1)
                kernel(i, j) = kernel_value(x(i, :), y(j, :), options)
            end do
        end do
    end subroutine kernel_matrix

    pure function kernel_value(x, y, options) result(value)
        real(dp), intent(in) :: x(:) !! First feature vector in fitted scaling space.
        real(dp), intent(in) :: y(:) !! Second feature vector with the same length as x.
        type(svm_options), intent(in) :: options !! Kernel code and polynomial/RBF/sigmoid parameters.
        real(dp) :: value
        real(dp) :: dot

        dot = dot_product(x, y)
        select case (options%kernel)
        case (svm_linear)
            value = dot
        case (svm_polynomial)
            value = (options%gamma * dot + options%coef0)**options%degree
        case (svm_radial)
            value = exp(-options%gamma * sum((x - y)**2))
        case (svm_sigmoid)
            value = tanh(options%gamma * dot + options%coef0)
        case default
            value = 0.0_dp
        end select
    end function kernel_value

    function pair_decision(pair, x, options) result(value)
        type(svm_pair_model), intent(in) :: pair !! Support-vector expansion for one binary decision function.
        real(dp), intent(in) :: x(:) !! Feature vector in the same scaled space as pair support vectors.
        type(svm_options), intent(in) :: options !! Kernel parameters used by the fitted pair model.
        real(dp) :: value
        integer :: i

        value = -pair%rho
        do i = 1, size(pair%coefficients)
            value = value + pair%coefficients(i) * kernel_value(pair%support_vectors(i, :), x, options)
        end do
    end function pair_decision

    subroutine prepare_predictors_fit(x, do_scale, requested_mask, scaled, center, scale_value, scale_mask)
        real(dp), intent(in) :: x(:, :) !! Original training predictor matrix to copy and optionally standardize.
        logical, intent(in) :: do_scale !! Global scaling switch corresponding to e1071's scale argument.
        logical, intent(in), optional :: requested_mask(:) !! Optional per-column selector applied when do_scale is true.
        real(dp), allocatable, intent(out) :: scaled(:, :) !! Allocated training matrix in fitted scaling space.
        real(dp), allocatable, intent(out) :: center(:) !! Per-column means used for scaling; zero for unscaled columns.
        real(dp), allocatable, intent(out) :: scale_value(:) !! Per-column sample standard deviations; one for unscaled columns.
        logical, allocatable, intent(out) :: scale_mask(:) !! Final per-column scaling selector stored with the model.
        real(dp) :: mean_value
        real(dp) :: ss
        integer :: j

        allocate(center(size(x, 2)), scale_value(size(x, 2)), scale_mask(size(x, 2)))
        center = 0.0_dp
        scale_value = 1.0_dp
        scale_mask = do_scale
        if (present(requested_mask)) then
            if (size(requested_mask) /= size(x, 2)) error stop "svm: scale_mask has wrong length"
            if (do_scale) scale_mask = requested_mask
        end if
        scaled = x
        do j = 1, size(x, 2)
            if (.not. scale_mask(j)) cycle
            mean_value = sum(x(:, j)) / real(size(x, 1), dp)
            ss = sum((x(:, j) - mean_value)**2)
            if (size(x, 1) <= 1 .or. ss <= 0.0_dp) error stop "svm: cannot scale a constant predictor"
            center(j) = mean_value
            scale_value(j) = sqrt(ss / real(size(x, 1) - 1, dp))
            scaled(:, j) = (scaled(:, j) - center(j)) / scale_value(j)
        end do
    end subroutine prepare_predictors_fit

    subroutine prepare_predictors_apply(model, x, scaled)
        type(svm_model), intent(in) :: model !! Fitted model providing expected feature count and scaling constants.
        real(dp), intent(in) :: x(:, :) !! New predictor matrix in original units.
        real(dp), allocatable, intent(out) :: scaled(:, :) !! Allocated predictor matrix transformed into fitted scaling space.
        integer :: j

        if (size(x, 2) /= model%n_features) error stop "svm prediction: variable count mismatch"
        scaled = x
        do j = 1, size(x, 2)
            if (model%scale_mask(j)) scaled(:, j) = (scaled(:, j) - model%center(j)) / model%scale_value(j)
        end do
    end subroutine prepare_predictors_apply

    subroutine collect_pair_rows(x, y, label_a, label_b, pair_x, binary_y, index, nfound)
        real(dp), intent(in) :: x(:, :) !! Full scaled predictor matrix.
        integer, intent(in) :: y(:) !! Full original class-label vector corresponding rowwise to x.
        integer, intent(in) :: label_a !! Original label mapped to +1 in the binary subproblem.
        integer, intent(in) :: label_b !! Original label mapped to -1 in the binary subproblem.
        real(dp), intent(out) :: pair_x(:, :) !! Preallocated pair-specific predictor matrix receiving selected rows.
        integer, intent(out) :: binary_y(:) !! Preallocated +1/-1 label vector for the binary subproblem.
        integer, intent(out) :: index(:) !! Preallocated original-row indices corresponding to pair_x.
        integer, intent(inout) :: nfound !! Running selected-row count; normally supplied as zero and returned at pair size.
        integer :: i

        do i = 1, size(y)
            if (y(i) /= label_a .and. y(i) /= label_b) cycle
            nfound = nfound + 1
            pair_x(nfound, :) = x(i, :)
            index(nfound) = i
            binary_y(nfound) = merge(1, -1, y(i) == label_a)
        end do
    end subroutine collect_pair_rows

    function class_weight(label, labels, weights) result(value)
        integer, intent(in) :: label !! Original class label whose multiplicative C weight is requested.
        integer, intent(in), optional :: labels(:) !! Labels associated one-to-one with optional weights.
        real(dp), intent(in), optional :: weights(:) !! Positive class-weight multipliers; unspecified classes use one.
        real(dp) :: value
        integer :: i

        value = 1.0_dp
        if (.not. present(labels) .and. .not. present(weights)) return
        if (.not. present(labels) .or. .not. present(weights)) then
            error stop "svm: class weight labels and weights must both be supplied"
        end if
        if (size(labels) /= size(weights)) error stop "svm: class weight arrays have different lengths"
        if (any(weights <= 0.0_dp)) error stop "svm: class weights must be positive"
        do i = 1, size(labels)
            if (labels(i) == label) then
                value = weights(i)
                return
            end if
        end do
    end function class_weight

    subroutine validate_options(options, nfeatures)
        type(svm_options), intent(in) :: options !! SVM controls to validate before fitting.
        integer, intent(in) :: nfeatures !! Positive predictor count used to validate kernel defaults and model shape.

        if (nfeatures < 1) error stop "svm: at least one predictor is required"
        if (options%kernel < svm_linear .or. options%kernel > svm_sigmoid) error stop "svm: invalid kernel"
        if (options%degree < 0) error stop "svm: polynomial degree must be nonnegative"
        if (abs(options%gamma) <= 0.0_dp .and. options%kernel == svm_radial) error stop "svm: radial gamma must be positive"
        if (options%gamma < 0.0_dp .and. abs(options%gamma + 1.0_dp) > 0.0_dp) then
            error stop "svm: gamma must be nonnegative or -1 for default"
        end if
        if (options%cost <= 0.0_dp) error stop "svm: cost must be positive"
        if (options%nu <= 0.0_dp .or. options%nu > 1.0_dp) error stop "svm: nu must be in (0,1]"
        if (options%tolerance <= 0.0_dp) error stop "svm: tolerance must be positive"
        if (options%epsilon < 0.0_dp) error stop "svm: epsilon must be nonnegative"
        if (options%max_iterations < 0) error stop "svm: max_iterations must be nonnegative"
        if (options%probability_seed <= 0) error stop "svm: probability_seed must be positive"
    end subroutine validate_options

    subroutine unique_integer_labels(y, labels)
        integer, intent(in) :: y(:) !! Integer class-label vector from which sorted distinct labels are extracted.
        integer, allocatable, intent(out) :: labels(:) !! Sorted distinct labels occurring in y.
        integer, allocatable :: work(:)
        integer :: i
        integer :: j
        integer :: n
        integer :: key

        if (size(y) == 0) error stop "svm: empty class vector"
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
    end subroutine unique_integer_labels


    subroutine binary_svc_probability(x, y, options, cp, cn, prob_a, prob_b)
        real(dp), intent(in) :: x(:, :) !! Binary training predictors in the same scaled space used by the final SVM fit.
        integer, intent(in) :: y(:) !! Binary labels encoded +1 and -1 for LIBSVM-style sigmoid calibration.
        type(svm_options), intent(in) :: options !! SVM controls; probability is disabled in each calibration-fold subfit.
        real(dp), intent(in) :: cp !! Positive-class C bound including any e1071 class-weight multiplier.
        real(dp), intent(in) :: cn !! Negative-class C bound including any e1071 class-weight multiplier.
        real(dp), intent(out) :: prob_a !! First Platt-sigmoid coefficient fitted to cross-validated decision values.
        real(dp), intent(out) :: prob_b !! Second Platt-sigmoid coefficient fitted to cross-validated decision values.
        integer, allocatable :: permutation(:)
        real(dp), allocatable :: train_x(:, :)
        integer, allocatable :: train_y(:)
        real(dp), allocatable :: decision(:)
        type(svm_pair_model) :: fold_model
        type(svm_options) :: fold_options
        integer :: fold
        integer :: first
        integer :: last
        integer :: ntrain
        integer :: i
        integer :: j
        integer :: k
        integer :: npositive
        integer :: nnegative

        allocate(decision(size(y)))
        call shuffled_indices(size(y), options%probability_seed, permutation)
        fold_options = options
        fold_options%probability = .false.
        do fold = 1, 5
            first = (fold - 1) * size(y) / 5 + 1
            last = fold * size(y) / 5
            if (last < first) cycle
            ntrain = size(y) - (last - first + 1)
            if (ntrain == 0) then
                decision(permutation(first:last)) = 0.0_dp
                cycle
            end if
            allocate(train_x(ntrain, size(x, 2)), train_y(ntrain))
            k = 0
            do j = 1, size(y)
                if (j >= first .and. j <= last) cycle
                k = k + 1
                train_x(k, :) = x(permutation(j), :)
                train_y(k) = y(permutation(j))
            end do
            npositive = count(train_y > 0)
            nnegative = count(train_y < 0)
            if (npositive == 0 .and. nnegative == 0) then
                decision(permutation(first:last)) = 0.0_dp
            else if (nnegative == 0) then
                decision(permutation(first:last)) = 1.0_dp
            else if (npositive == 0) then
                decision(permutation(first:last)) = -1.0_dp
            else
                call train_binary_classifier(train_x, train_y, fold_model, fold_options, cp, cn)
                do i = first, last
                    decision(permutation(i)) = pair_decision(fold_model, x(permutation(i), :), fold_options)
                end do
            end if
            deallocate(train_x, train_y)
        end do
        call sigmoid_train(decision, y, prob_a, prob_b)
    end subroutine binary_svc_probability

    subroutine svr_probability(x, y, options, sigma)
        real(dp), intent(in) :: x(:, :) !! Scaled regression predictors used to estimate the SVR residual distribution.
        real(dp), intent(in) :: y(:) !! Scaled regression response corresponding rowwise to x.
        type(svm_options), intent(in) :: options !! SVR controls; probability is disabled in the five fold subfits.
        real(dp), intent(out) :: sigma !! LIBSVM Laplace residual scale after five-standard-deviation outlier removal.
        integer, allocatable :: permutation(:)
        real(dp), allocatable :: prediction(:)
        real(dp), allocatable :: residual(:)
        real(dp), allocatable :: train_x(:, :)
        real(dp), allocatable :: train_y(:)
        type(svm_pair_model) :: fold_model
        type(svm_options) :: fold_options
        integer :: folds
        integer :: fold
        integer :: first
        integer :: last
        integer :: ntrain
        integer :: i
        integer :: j
        integer :: k
        integer :: excluded
        real(dp) :: mae
        real(dp) :: sd_laplace

        if (size(y) /= size(x, 1)) error stop "svr_probability: x/y mismatch"
        if (size(y) < 2) then
            sigma = 0.0_dp
            return
        end if
        folds = min(5, size(y))
        allocate(prediction(size(y)), residual(size(y)))
        call shuffled_indices(size(y), options%probability_seed, permutation)
        fold_options = options
        fold_options%probability = .false.
        do fold = 1, folds
            first = (fold - 1) * size(y) / folds + 1
            last = fold * size(y) / folds
            ntrain = size(y) - (last - first + 1)
            allocate(train_x(ntrain, size(x, 2)), train_y(ntrain))
            k = 0
            do j = 1, size(y)
                if (j >= first .and. j <= last) cycle
                k = k + 1
                train_x(k, :) = x(permutation(j), :)
                train_y(k) = y(permutation(j))
            end do
            call train_regressor(train_x, train_y, fold_model, fold_options)
            do i = first, last
                prediction(permutation(i)) = pair_decision(fold_model, x(permutation(i), :), fold_options)
            end do
            deallocate(train_x, train_y)
        end do
        residual = y - prediction
        mae = sum(abs(residual)) / real(size(residual), dp)
        sd_laplace = sqrt(2.0_dp) * mae
        excluded = count(abs(residual) > 5.0_dp * sd_laplace)
        if (excluded >= size(residual)) then
            sigma = mae
        else
            sigma = sum(abs(residual), mask=abs(residual) <= 5.0_dp * sd_laplace) / &
                    real(size(residual) - excluded, dp)
        end if
    end subroutine svr_probability

    subroutine shuffled_indices(n, seed, permutation)
        integer, intent(in) :: n !! Number of one-based indices to place in a Fisher-Yates random permutation.
        integer, intent(in) :: seed !! Positive deterministic seed used for native probability-calibration resampling.
        integer, allocatable, intent(out) :: permutation(:) !! Random permutation of the integers one through n.
        type(rng_state) :: state
        integer :: i
        integer :: j
        integer :: tmp

        allocate(permutation(n))
        permutation = [(i, i = 1, n)]
        call rng_seed(state, seed)
        do i = 1, n
            j = i + rng_integer(state, n - i + 1) - 1
            tmp = permutation(i)
            permutation(i) = permutation(j)
            permutation(j) = tmp
        end do
    end subroutine shuffled_indices

    subroutine sigmoid_train(decision, labels, a, b)
        real(dp), intent(in) :: decision(:) !! Binary SVM decision values used for Platt sigmoid calibration.
        integer, intent(in) :: labels(:) !! Corresponding +1/-1 labels used to define calibrated positive probabilities.
        real(dp), intent(out) :: a !! Fitted sigmoid slope in 1/(1+exp(A*f+B)).
        real(dp), intent(out) :: b !! Fitted sigmoid intercept in 1/(1+exp(A*f+B)).
        real(dp), allocatable :: target(:)
        real(dp) :: prior1
        real(dp) :: prior0
        real(dp) :: hi_target
        real(dp) :: lo_target
        real(dp) :: fval
        real(dp) :: fapb
        real(dp) :: p
        real(dp) :: q
        real(dp) :: h11
        real(dp) :: h22
        real(dp) :: h21
        real(dp) :: g1
        real(dp) :: g2
        real(dp) :: d1
        real(dp) :: d2
        real(dp) :: determinant
        real(dp) :: da
        real(dp) :: db
        real(dp) :: gd
        real(dp) :: step
        real(dp) :: new_a
        real(dp) :: new_b
        real(dp) :: new_f
        integer :: i
        integer :: iter

        prior1 = real(count(labels > 0), dp)
        prior0 = real(count(labels <= 0), dp)
        hi_target = (prior1 + 1.0_dp) / (prior1 + 2.0_dp)
        lo_target = 1.0_dp / (prior0 + 2.0_dp)
        allocate(target(size(labels)))
        where (labels > 0)
            target = hi_target
        elsewhere
            target = lo_target
        end where
        a = 0.0_dp
        b = log((prior0 + 1.0_dp) / (prior1 + 1.0_dp))
        fval = sigmoid_objective(decision, target, a, b)
        do iter = 1, 100
            h11 = 1.0e-12_dp
            h22 = 1.0e-12_dp
            h21 = 0.0_dp
            g1 = 0.0_dp
            g2 = 0.0_dp
            do i = 1, size(decision)
                fapb = decision(i) * a + b
                if (fapb >= 0.0_dp) then
                    p = exp(-fapb) / (1.0_dp + exp(-fapb))
                    q = 1.0_dp / (1.0_dp + exp(-fapb))
                else
                    p = 1.0_dp / (1.0_dp + exp(fapb))
                    q = exp(fapb) / (1.0_dp + exp(fapb))
                end if
                d2 = p * q
                h11 = h11 + decision(i)**2 * d2
                h22 = h22 + d2
                h21 = h21 + decision(i) * d2
                d1 = target(i) - p
                g1 = g1 + decision(i) * d1
                g2 = g2 + d1
            end do
            if (abs(g1) < 1.0e-5_dp .and. abs(g2) < 1.0e-5_dp) exit
            determinant = h11 * h22 - h21 * h21
            da = -(h22 * g1 - h21 * g2) / determinant
            db = -(-h21 * g1 + h11 * g2) / determinant
            gd = g1 * da + g2 * db
            step = 1.0_dp
            do while (step >= 1.0e-10_dp)
                new_a = a + step * da
                new_b = b + step * db
                new_f = sigmoid_objective(decision, target, new_a, new_b)
                if (new_f < fval + 0.0001_dp * step * gd) then
                    a = new_a
                    b = new_b
                    fval = new_f
                    exit
                end if
                step = 0.5_dp * step
            end do
            if (step < 1.0e-10_dp) exit
        end do
    end subroutine sigmoid_train

    pure function sigmoid_objective(decision, target, a, b) result(value)
        real(dp), intent(in) :: decision(:) !! Binary decision values used by the calibration log-likelihood.
        real(dp), intent(in) :: target(:) !! Smoothed 0/1 calibration targets paired with decision values.
        real(dp), intent(in) :: a !! Candidate sigmoid slope.
        real(dp), intent(in) :: b !! Candidate sigmoid intercept.
        real(dp) :: value
        real(dp) :: fapb
        integer :: i

        value = 0.0_dp
        do i = 1, size(decision)
            fapb = decision(i) * a + b
            if (fapb >= 0.0_dp) then
                value = value + target(i) * fapb + log(1.0_dp + exp(-fapb))
            else
                value = value + (target(i) - 1.0_dp) * fapb + log(1.0_dp + exp(fapb))
            end if
        end do
    end function sigmoid_objective

    pure function sigmoid_predict(decision, a, b) result(value)
        real(dp), intent(in) :: decision !! Binary SVM decision value.
        real(dp), intent(in) :: a !! Fitted Platt sigmoid slope.
        real(dp), intent(in) :: b !! Fitted Platt sigmoid intercept.
        real(dp) :: value
        real(dp) :: z

        z = decision * a + b
        if (z >= 0.0_dp) then
            value = exp(-z) / (1.0_dp + exp(-z))
        else
            value = 1.0_dp / (1.0_dp + exp(z))
        end if
    end function sigmoid_predict

    subroutine multiclass_probability(r, probability)
        real(dp), intent(in) :: r(:, :) !! Pairwise class probabilities where r(i,j) estimates P(class i | i versus j).
        real(dp), intent(out) :: probability(:) !! Coupled class probabilities summing approximately to one.
        real(dp), allocatable :: q(:, :)
        real(dp), allocatable :: qp(:)
        real(dp) :: pqp
        real(dp) :: max_error
        real(dp) :: error
        real(dp) :: diff
        real(dp) :: eps
        integer :: k
        integer :: t
        integer :: j
        integer :: iter
        integer :: max_iter

        k = size(probability)
        if (size(r, 1) /= k .or. size(r, 2) /= k) error stop "multiclass_probability: shape mismatch"
        allocate(q(k, k), qp(k))
        q = 0.0_dp
        probability = 1.0_dp / real(k, dp)
        do t = 1, k
            do j = 1, t - 1
                q(t, t) = q(t, t) + r(j, t)**2
                q(t, j) = q(j, t)
            end do
            do j = t + 1, k
                q(t, t) = q(t, t) + r(j, t)**2
                q(t, j) = -r(j, t) * r(t, j)
            end do
        end do
        eps = 0.005_dp / real(k, dp)
        max_iter = max(100, k)
        do iter = 1, max_iter
            qp = matmul(q, probability)
            pqp = dot_product(probability, qp)
            max_error = 0.0_dp
            do t = 1, k
                error = abs(qp(t) - pqp)
                max_error = max(max_error, error)
            end do
            if (max_error < eps) exit
            do t = 1, k
                if (q(t, t) <= tiny(1.0_dp)) cycle
                diff = (-qp(t) + pqp) / q(t, t)
                probability(t) = probability(t) + diff
                pqp = (pqp + diff * (diff * q(t, t) + 2.0_dp * qp(t))) / (1.0_dp + diff)**2
                do j = 1, k
                    qp(j) = (qp(j) + diff * q(t, j)) / (1.0_dp + diff)
                    probability(j) = probability(j) / (1.0_dp + diff)
                end do
            end do
        end do
    end subroutine multiclass_probability

    pure function row_fold(row, folds, fold_id) result(fold)
        integer, intent(in) :: row !! One-based observation row whose cross-validation fold is requested.
        integer, intent(in) :: folds !! Total number of cyclic folds when fold_id is absent.
        integer, intent(in), optional :: fold_id(:) !! Optional explicit one-based fold assignment for all observations.
        integer :: fold

        if (present(fold_id)) then
            fold = fold_id(row)
        else
            fold = modulo(row - 1, folds) + 1
        end if
    end function row_fold

end module e1071_svm
