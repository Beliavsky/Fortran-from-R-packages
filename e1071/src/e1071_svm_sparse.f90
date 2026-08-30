module e1071_svm_sparse
    use e1071_kinds, only: dp
    use e1071_sparse, only: matrix_csr, csr_to_dense
    use e1071_svm, only: svm_model, svm_options, svm_pair_model, svm_fit_classification, svm_fit_regression, &
                         svm_predict_classification, svm_predict_regression, svm_predict_one_class, &
                         svm_c_classification, svm_nu_classification, svm_one_classification, svm_eps_regression, &
                         svm_nu_regression, svm_linear, svm_polynomial, svm_radial, svm_sigmoid
    use e1071_svm_solver, only: solve_dual, initialize_nu_svc_alpha
    implicit none
    private

    public :: svm_fit_classification_csr, svm_fit_regression_csr, svm_fit_one_class_csr
    public :: svm_predict_classification_csr, svm_predict_regression_csr, svm_predict_one_class_csr

contains

    subroutine svm_fit_classification_csr(x, y, model, options, class_weight_labels, class_weights, scale_mask)
        type(matrix_csr), intent(in) :: x !! Sparse CSR training predictors with observations stored by row.
        integer, intent(in) :: y(:) !! Integer class labels corresponding one-to-one with CSR rows.
        type(svm_model), intent(out) :: model !! Fitted classification SVM using sparse kernel construction.
        type(svm_options), intent(in), optional :: options !! Optional SVM controls; sparse inputs always disable scaling.
        integer, intent(in), optional :: class_weight_labels(:) !! Optional class labels receiving C multipliers.
        real(dp), intent(in), optional :: class_weights(:) !! Positive C multipliers paired with class_weight_labels.
        logical, intent(in), optional :: scale_mask(:) !! Compatibility argument validated for shape but ignored for sparse input.
        type(svm_options) :: sparse_options
        type(svm_model) :: calibration_model
        type(svm_options) :: calibration_options
        integer, allocatable :: row_index(:)
        integer, allocatable :: binary_y(:)
        real(dp), allocatable :: dense_pair(:, :)
        integer, allocatable :: calibration_labels(:)
        real(dp), allocatable :: calibration_weights(:)
        integer :: nclass
        integer :: npair
        integer :: pair_index
        integer :: a
        integer :: b
        integer :: i
        integer :: n
        real(dp) :: cp
        real(dp) :: cn

        if (size(y) /= x%nrow) error stop "svm_fit_classification_csr: x/y row mismatch"
        call validate_sparse_scale_argument(x, scale_mask)
        sparse_options = svm_options()
        if (present(options)) sparse_options = options
        sparse_options%scale = .false.
        call validate_sparse_options(sparse_options, x%ncol)
        if (sparse_options%svm_type /= svm_c_classification .and. &
            sparse_options%svm_type /= svm_nu_classification) then
            error stop "svm_fit_classification_csr: invalid svm_type"
        end if
        call unique_integer_labels_sparse(y, model%class_labels)
        nclass = size(model%class_labels)
        if (nclass < 2) error stop "svm_fit_classification_csr: at least two classes are required"
        model%n_features = x%ncol
        model%options = sparse_options
        if (model%options%gamma < 0.0_dp) model%options%gamma = 1.0_dp / real(x%ncol, dp)
        call set_identity_scaling(model)
        npair = nclass * (nclass - 1) / 2
        allocate(model%pair(npair))
        pair_index = 0
        do a = 1, nclass - 1
            do b = a + 1, nclass
                pair_index = pair_index + 1
                n = count(y == model%class_labels(a)) + count(y == model%class_labels(b))
                allocate(row_index(n), binary_y(n))
                i = 0
                call collect_sparse_pair_rows(y, model%class_labels(a), model%class_labels(b), row_index, binary_y, i)
                cp = sparse_options%cost * sparse_class_weight(model%class_labels(a), class_weight_labels, class_weights)
                cn = sparse_options%cost * sparse_class_weight(model%class_labels(b), class_weight_labels, class_weights)
                call train_binary_classifier_csr(x, row_index, binary_y, model%pair(pair_index), model%options, cp, cn)
                model%pair(pair_index)%positive_label = model%class_labels(a)
                model%pair(pair_index)%negative_label = model%class_labels(b)
                if (model%options%probability) then
                    call csr_selected_rows_to_dense(x, row_index, dense_pair)
                    calibration_options = model%options
                    calibration_options%scale = .false.
                    allocate(calibration_labels(2), calibration_weights(2))
                    calibration_labels = [1, -1]
                    calibration_weights = [cp / model%options%cost, cn / model%options%cost]
                    call svm_fit_classification(dense_pair, binary_y, calibration_model, calibration_options, &
                                                calibration_labels, calibration_weights)
                    model%pair(pair_index)%prob_a = calibration_model%pair(1)%prob_a
                    model%pair(pair_index)%prob_b = calibration_model%pair(1)%prob_b
                    deallocate(dense_pair, calibration_labels, calibration_weights)
                end if
                deallocate(row_index, binary_y)
            end do
        end do
        model%probability_fitted = model%options%probability
    end subroutine svm_fit_classification_csr

    subroutine svm_fit_regression_csr(x, y, model, options, scale_mask)
        type(matrix_csr), intent(in) :: x !! Sparse CSR training predictors with observations stored by row.
        real(dp), intent(in) :: y(:) !! Numeric response vector corresponding one-to-one with CSR rows.
        type(svm_model), intent(out) :: model !! Fitted sparse-input regression SVM with unscaled response semantics.
        type(svm_options), intent(in), optional :: options !! Optional epsilon- or nu-SVR controls; scaling is disabled.
        logical, intent(in), optional :: scale_mask(:) !! Compatibility argument validated for shape but ignored for sparse input.
        type(svm_options) :: sparse_options
        type(svm_model) :: calibration_model
        real(dp), allocatable :: dense(:, :)

        if (size(y) /= x%nrow) error stop "svm_fit_regression_csr: x/y row mismatch"
        call validate_sparse_scale_argument(x, scale_mask)
        sparse_options = svm_options(svm_type=svm_eps_regression)
        if (present(options)) sparse_options = options
        sparse_options%scale = .false.
        call validate_sparse_options(sparse_options, x%ncol)
        if (sparse_options%svm_type /= svm_eps_regression .and. sparse_options%svm_type /= svm_nu_regression) then
            error stop "svm_fit_regression_csr: invalid svm_type"
        end if
        model%n_features = x%ncol
        model%options = sparse_options
        if (model%options%gamma < 0.0_dp) model%options%gamma = 1.0_dp / real(x%ncol, dp)
        model%y_center = 0.0_dp
        model%y_scale = 1.0_dp
        call set_identity_scaling(model)
        allocate(model%pair(1))
        call train_regressor_csr(x, y, model%pair(1), model%options)
        if (model%options%probability) then
            call csr_to_dense(x, dense)
            call svm_fit_regression(dense, y, calibration_model, model%options)
            model%svr_probability_sigma = calibration_model%svr_probability_sigma
            model%probability_fitted = .true.
        end if
    end subroutine svm_fit_regression_csr

    subroutine svm_fit_one_class_csr(x, model, options, scale_mask)
        type(matrix_csr), intent(in) :: x !! Sparse CSR observations defining the target one-class support.
        type(svm_model), intent(out) :: model !! Fitted sparse-input one-class SVM with unscaled predictor semantics.
        type(svm_options), intent(in), optional :: options !! Optional one-class SVM controls; sparse scaling is disabled.
        logical, intent(in), optional :: scale_mask(:) !! Compatibility argument validated for shape but ignored for sparse input.
        type(svm_options) :: sparse_options

        call validate_sparse_scale_argument(x, scale_mask)
        sparse_options = svm_options(svm_type=svm_one_classification)
        if (present(options)) sparse_options = options
        sparse_options%scale = .false.
        call validate_sparse_options(sparse_options, x%ncol)
        if (sparse_options%svm_type /= svm_one_classification) then
            error stop "svm_fit_one_class_csr: invalid svm_type"
        end if
        model%n_features = x%ncol
        model%options = sparse_options
        if (model%options%gamma < 0.0_dp) model%options%gamma = 1.0_dp / real(x%ncol, dp)
        call set_identity_scaling(model)
        allocate(model%pair(1))
        call train_one_class_csr(x, model%pair(1), model%options)
    end subroutine svm_fit_one_class_csr

    subroutine svm_predict_classification_csr(model, x, class, decision, probability)
        type(svm_model), intent(in) :: model !! Fitted classification SVM used for sparse-row prediction.
        type(matrix_csr), intent(in) :: x !! Sparse CSR prediction matrix with the model's fitted feature count.
        integer, allocatable, intent(out) :: class(:) !! Predicted original integer class label for each sparse row.
        real(dp), allocatable, intent(out), optional :: decision(:, :) !! Optional pairwise decision values by row.
        real(dp), allocatable, intent(out), optional :: probability(:, :) !! Optional coupled class probabilities by row.
        real(dp), allocatable :: row_matrix(:, :)
        integer, allocatable :: row_class(:)
        real(dp), allocatable :: row_decision(:, :)
        real(dp), allocatable :: row_probability(:, :)
        integer :: i

        call validate_prediction_shape(model, x)
        allocate(class(x%nrow), row_matrix(1, x%ncol))
        if (present(decision)) allocate(decision(x%nrow, size(model%pair)))
        if (present(probability)) allocate(probability(x%nrow, size(model%class_labels)))
        do i = 1, x%nrow
            call csr_row_to_dense(x, i, row_matrix(1, :))
            if (present(probability)) then
                call svm_predict_classification(model, row_matrix, row_class, row_decision, row_probability)
                probability(i, :) = row_probability(1, :)
                if (present(decision)) decision(i, :) = row_decision(1, :)
            else if (present(decision)) then
                call svm_predict_classification(model, row_matrix, row_class, row_decision)
                decision(i, :) = row_decision(1, :)
            else
                call svm_predict_classification(model, row_matrix, row_class)
            end if
            class(i) = row_class(1)
        end do
    end subroutine svm_predict_classification_csr

    subroutine svm_predict_regression_csr(model, x, prediction, decision)
        type(svm_model), intent(in) :: model !! Fitted regression SVM used for sparse-row prediction.
        type(matrix_csr), intent(in) :: x !! Sparse CSR prediction matrix with the model's fitted feature count.
        real(dp), allocatable, intent(out) :: prediction(:) !! Predicted numeric response for each sparse row.
        real(dp), allocatable, intent(out), optional :: decision(:) !! Optional raw scaled-space SVR decision values.
        real(dp), allocatable :: row_matrix(:, :)
        real(dp), allocatable :: row_prediction(:)
        real(dp), allocatable :: row_decision(:)
        integer :: i

        call validate_prediction_shape(model, x)
        allocate(prediction(x%nrow), row_matrix(1, x%ncol))
        if (present(decision)) allocate(decision(x%nrow))
        do i = 1, x%nrow
            call csr_row_to_dense(x, i, row_matrix(1, :))
            if (present(decision)) then
                call svm_predict_regression(model, row_matrix, row_prediction, row_decision)
                decision(i) = row_decision(1)
            else
                call svm_predict_regression(model, row_matrix, row_prediction)
            end if
            prediction(i) = row_prediction(1)
        end do
    end subroutine svm_predict_regression_csr

    subroutine svm_predict_one_class_csr(model, x, accepted, decision)
        type(svm_model), intent(in) :: model !! Fitted one-class SVM used for sparse-row prediction.
        type(matrix_csr), intent(in) :: x !! Sparse CSR prediction matrix with the model's fitted feature count.
        logical, allocatable, intent(out) :: accepted(:) !! True for sparse rows on the accepted side of the learned boundary.
        real(dp), allocatable, intent(out), optional :: decision(:) !! Optional signed one-class decision value for each row.
        real(dp), allocatable :: row_matrix(:, :)
        logical, allocatable :: row_accepted(:)
        real(dp), allocatable :: row_decision(:)
        integer :: i

        call validate_prediction_shape(model, x)
        allocate(accepted(x%nrow), row_matrix(1, x%ncol))
        if (present(decision)) allocate(decision(x%nrow))
        do i = 1, x%nrow
            call csr_row_to_dense(x, i, row_matrix(1, :))
            if (present(decision)) then
                call svm_predict_one_class(model, row_matrix, row_accepted, row_decision)
                decision(i) = row_decision(1)
            else
                call svm_predict_one_class(model, row_matrix, row_accepted)
            end if
            accepted(i) = row_accepted(1)
        end do
    end subroutine svm_predict_one_class_csr

    subroutine train_binary_classifier_csr(x, rows, y, pair, options, cp, cn)
        type(matrix_csr), intent(in) :: x !! Full sparse predictor matrix containing the selected pair rows.
        integer, intent(in) :: rows(:) !! One-based original row indices forming the binary subproblem.
        integer, intent(in) :: y(:) !! Binary signs encoded +1 for the positive class and -1 for the negative class.
        type(svm_pair_model), intent(out) :: pair !! Fitted support-vector expansion materialized only for selected SV rows.
        type(svm_options), intent(in) :: options !! Validated sparse classification controls with resolved gamma.
        real(dp), intent(in) :: cp !! Positive-class dual upper bound after class weighting.
        real(dp), intent(in) :: cn !! Negative-class dual upper bound after class weighting.
        real(dp), allocatable :: kernel(:, :)
        real(dp), allocatable :: q(:, :)
        real(dp), allocatable :: alpha(:)
        real(dp), allocatable :: linear(:)
        integer, allocatable :: signs(:)
        real(dp) :: rho
        real(dp) :: r
        integer :: i
        integer :: nsv

        if (size(rows) /= size(y)) error stop "train_binary_classifier_csr: row/label mismatch"
        if (count(y > 0) == 0 .or. count(y < 0) == 0) error stop "train_binary_classifier_csr: both signs required"
        call kernel_matrix_csr_subset(x, rows, options, kernel)
        allocate(q(size(y), size(y)), alpha(size(y)), linear(size(y)), signs(size(y)))
        signs = merge(1, -1, y > 0)
        do i = 1, size(y)
            q(i, :) = real(signs(i), dp) * real(signs, dp) * kernel(i, :)
        end do
        if (options%svm_type == svm_c_classification) then
            alpha = 0.0_dp
            linear = -1.0_dp
            call solve_dual(q, linear, signs, alpha, cp, cn, options%tolerance, .false., &
                            options%max_iterations, rho, r)
            alpha = alpha * real(signs, dp)
        else
            call initialize_nu_svc_alpha(signs, options%nu, alpha)
            linear = 0.0_dp
            call solve_dual(q, linear, signs, alpha, 1.0_dp, 1.0_dp, options%tolerance, .true., &
                            options%max_iterations, rho, r)
            if (r <= 0.0_dp) error stop "train_binary_classifier_csr: invalid nu-SVC scaling factor"
            alpha = alpha * real(signs, dp) / r
            rho = rho / r
        end if
        nsv = count(abs(alpha) > 1.0e-10_dp)
        allocate(pair%support_vectors(nsv, x%ncol), pair%coefficients(nsv))
        nsv = 0
        do i = 1, size(alpha)
            if (abs(alpha(i)) <= 1.0e-10_dp) cycle
            nsv = nsv + 1
            call csr_row_to_dense(x, rows(i), pair%support_vectors(nsv, :))
            pair%coefficients(nsv) = alpha(i)
        end do
        pair%rho = rho
    end subroutine train_binary_classifier_csr

    subroutine train_regressor_csr(x, y, pair, options)
        type(matrix_csr), intent(in) :: x !! Sparse regression predictors whose kernel is constructed without dense expansion.
        real(dp), intent(in) :: y(:) !! Unscaled numeric response vector, one value per sparse row.
        type(svm_pair_model), intent(out) :: pair !! Fitted SVR support-vector expansion materialized only for selected SV rows.
        type(svm_options), intent(in) :: options !! Validated epsilon- or nu-SVR controls with resolved gamma.
        real(dp), allocatable :: kernel(:, :)
        real(dp), allocatable :: q(:, :)
        real(dp), allocatable :: alpha2(:)
        real(dp), allocatable :: alpha(:)
        real(dp), allocatable :: linear(:)
        integer, allocatable :: signs(:)
        integer, allocatable :: rows(:)
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
        allocate(rows(n))
        rows = [(i, i=1, n)]
        call kernel_matrix_csr_subset(x, rows, options, kernel)
        allocate(q(2 * n, 2 * n), alpha2(2 * n), alpha(n), linear(2 * n), signs(2 * n))
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
            linear(:n) = options%epsilon - y
            linear(n + 1:) = options%epsilon + y
            call solve_dual(q, linear, signs, alpha2, options%cost, options%cost, options%tolerance, .false., &
                            options%max_iterations, rho, r)
        else
            sum_alpha = options%cost * options%nu * real(n, dp) / 2.0_dp
            do i = 1, n
                alpha2(i) = min(sum_alpha, options%cost)
                alpha2(i + n) = alpha2(i)
                sum_alpha = sum_alpha - alpha2(i)
            end do
            linear(:n) = -y
            linear(n + 1:) = y
            call solve_dual(q, linear, signs, alpha2, options%cost, options%cost, options%tolerance, .true., &
                            options%max_iterations, rho, r)
        end if
        alpha = alpha2(:n) - alpha2(n + 1:)
        nsv = count(abs(alpha) > 1.0e-10_dp)
        allocate(pair%support_vectors(nsv, x%ncol), pair%coefficients(nsv))
        nsv = 0
        do i = 1, n
            if (abs(alpha(i)) <= 1.0e-10_dp) cycle
            nsv = nsv + 1
            call csr_row_to_dense(x, i, pair%support_vectors(nsv, :))
            pair%coefficients(nsv) = alpha(i)
        end do
        pair%rho = rho
    end subroutine train_regressor_csr

    subroutine train_one_class_csr(x, pair, options)
        type(matrix_csr), intent(in) :: x !! Sparse observations defining the target one-class distribution.
        type(svm_pair_model), intent(out) :: pair !! Fitted one-class support-vector coefficients and threshold.
        type(svm_options), intent(in) :: options !! Validated one-class controls with resolved gamma and nu.
        real(dp), allocatable :: kernel(:, :)
        real(dp), allocatable :: alpha(:)
        real(dp), allocatable :: linear(:)
        integer, allocatable :: signs(:)
        integer, allocatable :: rows(:)
        real(dp) :: rho
        real(dp) :: r
        real(dp) :: remaining
        integer :: n
        integer :: full
        integer :: nsv
        integer :: i

        n = x%nrow
        allocate(rows(n))
        rows = [(i, i=1, n)]
        call kernel_matrix_csr_subset(x, rows, options, kernel)
        allocate(alpha(n), linear(n), signs(n))
        alpha = 0.0_dp
        linear = 0.0_dp
        signs = 1
        full = int(options%nu * real(n, dp))
        if (full > 0) alpha(:min(full, n)) = 1.0_dp
        remaining = options%nu * real(n, dp) - real(full, dp)
        if (full < n .and. remaining > 0.0_dp) alpha(full + 1) = remaining
        call solve_dual(kernel, linear, signs, alpha, 1.0_dp, 1.0_dp, options%tolerance, .false., &
                        options%max_iterations, rho, r)
        nsv = count(alpha > 1.0e-10_dp)
        allocate(pair%support_vectors(nsv, x%ncol), pair%coefficients(nsv))
        nsv = 0
        do i = 1, n
            if (alpha(i) <= 1.0e-10_dp) cycle
            nsv = nsv + 1
            call csr_row_to_dense(x, i, pair%support_vectors(nsv, :))
            pair%coefficients(nsv) = alpha(i)
        end do
        pair%rho = rho
    end subroutine train_one_class_csr

    subroutine kernel_matrix_csr_subset(x, rows, options, kernel)
        type(matrix_csr), intent(in) :: x !! Sparse matrix supplying rows used by the kernel subproblem.
        integer, intent(in) :: rows(:) !! One-based original row indices defining the kernel matrix order.
        type(svm_options), intent(in) :: options !! Kernel family and resolved gamma/degree/coef0 parameters.
        real(dp), allocatable, intent(out) :: kernel(:, :) !! Dense symmetric kernel matrix for the selected sparse rows.
        real(dp), allocatable :: norm2(:)
        real(dp) :: dot
        real(dp) :: distance2
        integer :: i
        integer :: j

        allocate(kernel(size(rows), size(rows)), norm2(size(rows)))
        do i = 1, size(rows)
            norm2(i) = csr_row_norm2(x, rows(i))
        end do
        do i = 1, size(rows)
            do j = i, size(rows)
                dot = csr_row_dot(x, rows(i), rows(j))
                select case (options%kernel)
                case (svm_linear)
                    kernel(i, j) = dot
                case (svm_polynomial)
                    kernel(i, j) = (options%gamma * dot + options%coef0)**options%degree
                case (svm_radial)
                    distance2 = max(0.0_dp, norm2(i) + norm2(j) - 2.0_dp * dot)
                    kernel(i, j) = exp(-options%gamma * distance2)
                case (svm_sigmoid)
                    kernel(i, j) = tanh(options%gamma * dot + options%coef0)
                case default
                    error stop "kernel_matrix_csr_subset: invalid kernel"
                end select
                kernel(j, i) = kernel(i, j)
            end do
        end do
    end subroutine kernel_matrix_csr_subset

    function csr_row_dot(x, row_a, row_b) result(value)
        type(matrix_csr), intent(in) :: x !! Sparse CSR matrix containing both rows in the inner product.
        integer, intent(in) :: row_a !! One-based first row index in the sparse inner product.
        integer, intent(in) :: row_b !! One-based second row index in the sparse inner product.
        real(dp) :: value
        integer :: ka
        integer :: kb

        value = 0.0_dp
        do ka = x%row_pointer(row_a), x%row_pointer(row_a + 1) - 1
            do kb = x%row_pointer(row_b), x%row_pointer(row_b + 1) - 1
                if (x%col_index(ka) == x%col_index(kb)) value = value + x%values(ka) * x%values(kb)
            end do
        end do
    end function csr_row_dot

    function csr_row_norm2(x, row) result(value)
        type(matrix_csr), intent(in) :: x !! Sparse CSR matrix containing the row whose squared Euclidean norm is requested.
        integer, intent(in) :: row !! One-based row index whose nonzero values are squared and summed.
        real(dp) :: value
        integer :: k

        value = 0.0_dp
        do k = x%row_pointer(row), x%row_pointer(row + 1) - 1
            value = value + x%values(k) * x%values(k)
        end do
    end function csr_row_norm2

    subroutine csr_selected_rows_to_dense(x, rows, dense)
        type(matrix_csr), intent(in) :: x !! Sparse matrix containing the selected rows to materialize.
        integer, intent(in) :: rows(:) !! One-based row indices copied into output order.
        real(dp), allocatable, intent(out) :: dense(:, :) !! Dense selected-row matrix used only by probability calibration.
        integer :: i

        allocate(dense(size(rows), x%ncol))
        do i = 1, size(rows)
            call csr_row_to_dense(x, rows(i), dense(i, :))
        end do
    end subroutine csr_selected_rows_to_dense

    subroutine collect_sparse_pair_rows(y, label_a, label_b, rows, binary_y, nfound)
        integer, intent(in) :: y(:) !! Full original class-label vector corresponding rowwise to the sparse matrix.
        integer, intent(in) :: label_a !! Original label mapped to +1 in the pairwise binary subproblem.
        integer, intent(in) :: label_b !! Original label mapped to -1 in the pairwise binary subproblem.
        integer, intent(out) :: rows(:) !! Original sparse-row indices selected for the two requested classes.
        integer, intent(out) :: binary_y(:) !! Pairwise +1/-1 labels in the same order as `rows`.
        integer, intent(inout) :: nfound !! Running selected-row count, normally supplied as zero.
        integer :: i

        do i = 1, size(y)
            if (y(i) /= label_a .and. y(i) /= label_b) cycle
            nfound = nfound + 1
            rows(nfound) = i
            binary_y(nfound) = merge(1, -1, y(i) == label_a)
        end do
    end subroutine collect_sparse_pair_rows

    subroutine unique_integer_labels_sparse(y, labels)
        integer, intent(in) :: y(:) !! Integer class vector from which sorted distinct labels are extracted.
        integer, allocatable, intent(out) :: labels(:) !! Sorted distinct labels occurring in `y`.
        integer, allocatable :: work(:)
        integer :: i
        integer :: j
        integer :: n
        integer :: key

        if (size(y) == 0) error stop "svm sparse fit: empty class vector"
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
    end subroutine unique_integer_labels_sparse

    function sparse_class_weight(label, labels, weights) result(value)
        integer, intent(in) :: label !! Original class label whose multiplicative C weight is requested.
        integer, intent(in), optional :: labels(:) !! Labels associated one-to-one with optional class weights.
        real(dp), intent(in), optional :: weights(:) !! Positive class-weight multipliers; omitted classes use one.
        real(dp) :: value
        integer :: i

        value = 1.0_dp
        if (.not. present(labels) .and. .not. present(weights)) return
        if (.not. present(labels) .or. .not. present(weights)) then
            error stop "svm sparse fit: class weight labels and weights must both be supplied"
        end if
        if (size(labels) /= size(weights)) error stop "svm sparse fit: class weight array length mismatch"
        if (any(weights <= 0.0_dp)) error stop "svm sparse fit: class weights must be positive"
        do i = 1, size(labels)
            if (labels(i) == label) then
                value = weights(i)
                return
            end if
        end do
    end function sparse_class_weight

    subroutine set_identity_scaling(model)
        type(svm_model), intent(inout) :: model !! Sparse-input model receiving identity predictor scaling metadata.

        allocate(model%center(model%n_features), model%scale_value(model%n_features), model%scale_mask(model%n_features))
        model%center = 0.0_dp
        model%scale_value = 1.0_dp
        model%scale_mask = .false.
    end subroutine set_identity_scaling

    subroutine validate_sparse_options(options, nfeatures)
        type(svm_options), intent(in) :: options !! Sparse SVM controls checked against the dense engine's parameter rules.
        integer, intent(in) :: nfeatures !! Positive predictor count used to validate default gamma and model shape.

        if (nfeatures < 1) error stop "svm sparse fit: at least one predictor is required"
        if (options%kernel < svm_linear .or. options%kernel > svm_sigmoid) error stop "svm sparse fit: invalid kernel"
        if (options%degree < 0) error stop "svm sparse fit: polynomial degree must be nonnegative"
        if (abs(options%gamma) <= 0.0_dp .and. options%kernel == svm_radial) then
            error stop "svm sparse fit: radial gamma must be positive"
        end if
        if (options%gamma < 0.0_dp .and. abs(options%gamma + 1.0_dp) > 0.0_dp) then
            error stop "svm sparse fit: gamma must be nonnegative or -1 for default"
        end if
        if (options%cost <= 0.0_dp) error stop "svm sparse fit: cost must be positive"
        if (options%nu <= 0.0_dp .or. options%nu > 1.0_dp) error stop "svm sparse fit: nu must be in (0,1]"
        if (options%tolerance <= 0.0_dp) error stop "svm sparse fit: tolerance must be positive"
        if (options%epsilon < 0.0_dp) error stop "svm sparse fit: epsilon must be nonnegative"
        if (options%max_iterations < 0) error stop "svm sparse fit: max_iterations must be nonnegative"
        if (options%probability_seed <= 0) error stop "svm sparse fit: probability_seed must be positive"
    end subroutine validate_sparse_options

    subroutine validate_sparse_scale_argument(x, scale_mask)
        type(matrix_csr), intent(in) :: x !! Sparse matrix providing the feature count used for scale-mask validation.
        logical, intent(in), optional :: scale_mask(:) !! Optional compatibility mask whose values are intentionally ignored.

        if (present(scale_mask)) then
            if (size(scale_mask) /= x%ncol) error stop "svm sparse fit: scale_mask length mismatch"
        end if
    end subroutine validate_sparse_scale_argument

    subroutine validate_prediction_shape(model, x)
        type(svm_model), intent(in) :: model !! Fitted model providing the required feature count.
        type(matrix_csr), intent(in) :: x !! Sparse prediction matrix whose feature count is checked.

        if (x%ncol /= model%n_features) error stop "svm sparse prediction: variable count mismatch"
    end subroutine validate_prediction_shape

    subroutine csr_row_to_dense(x, row, values)
        type(matrix_csr), intent(in) :: x !! Sparse CSR matrix containing the row to expand.
        integer, intent(in) :: row !! One-based row index to expand; valid range is 1 through `x%nrow`.
        real(dp), intent(out) :: values(:) !! Dense feature vector of length `x%ncol` receiving the selected sparse row.
        integer :: k

        if (row < 1 .or. row > x%nrow) error stop "csr_row_to_dense: row index out of range"
        if (size(values) /= x%ncol) error stop "csr_row_to_dense: output feature count mismatch"
        values = 0.0_dp
        do k = x%row_pointer(row), x%row_pointer(row + 1) - 1
            values(x%col_index(k)) = x%values(k)
        end do
    end subroutine csr_row_to_dense

end module e1071_svm_sparse
