program test_svm
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use e1071
    implicit none

    real(dp) :: x(6, 2)
    integer :: y(6)
    real(dp) :: xr(10, 1)
    real(dp) :: yr(10)
    real(dp) :: one_x(5, 2)
    real(dp) :: one_query(2, 2)
    integer, allocatable :: class(:)
    integer, allocatable :: cv_class(:)
    logical, allocatable :: accepted(:)
    real(dp), allocatable :: decision(:, :)
    real(dp), allocatable :: probability(:, :)
    real(dp), allocatable :: prediction(:)
    real(dp), allocatable :: cv_prediction(:)
    real(dp), allocatable :: dense_prediction(:)
    real(dp), allocatable :: sparse_prediction(:)
    integer, allocatable :: dense_class(:)
    integer, allocatable :: sparse_class(:)
    real(dp), allocatable :: dense_decision(:, :)
    real(dp), allocatable :: sparse_decision(:, :)
    type(svm_options) :: options
    type(svm_model) :: model
    type(svm_model) :: dense_model
    type(svm_model) :: sparse_model
    type(tune_result) :: tune_fit
    type(matrix_csr) :: sparse_x
    integer :: i

    x = reshape([0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, &
                 3.0_dp, 3.0_dp, 3.0_dp, 4.0_dp, 4.0_dp, 3.0_dp], [6, 2], order=[2, 1])
    y = [-1, -1, -1, 1, 1, 1]
    options = svm_options(svm_type=svm_c_classification, kernel=svm_radial, gamma=0.5_dp, &
                          cost=1.0_dp, scale=.false.)
    call svm_fit_classification(x, y, model, options)
    call svm_predict_classification(model, x, class, decision=decision)
    if (any(class /= y)) error stop "C-SVC classification failed"
    if (maxval(abs(abs(decision(:, 1)) - 1.0_dp)) > 0.02_dp) error stop "C-SVC LIBSVM parity fixture failed"

    options = svm_options(svm_type=svm_c_classification, kernel=svm_polynomial, degree=2, gamma=0.5_dp, &
                          coef0=1.0_dp, cost=1.0_dp, scale=.false.)
    call svm_fit_classification(x, y, model, options)
    call svm_predict_classification(model, x, class, decision=decision)
    if (.not. all(ieee_is_finite(decision))) error stop "polynomial SVM kernel failed"

    options = svm_options(svm_type=svm_c_classification, kernel=svm_sigmoid, gamma=0.05_dp, &
                          coef0=0.0_dp, cost=1.0_dp, scale=.false.)
    call svm_fit_classification(x, y, model, options)
    call svm_predict_classification(model, x, class, decision=decision)
    if (.not. all(ieee_is_finite(decision))) error stop "sigmoid SVM kernel failed"

    options = svm_options(svm_type=svm_c_classification, kernel=svm_radial, gamma=0.5_dp, &
                          cost=1.0_dp, scale=.false.)
    call dense_to_csr(x, sparse_x)
    call svm_fit_classification(x, y, dense_model, options)
    call svm_predict_classification(dense_model, x, dense_class, decision=dense_decision)
    call svm_fit_classification_csr(sparse_x, y, sparse_model, options)
    call svm_predict_classification_csr(sparse_model, sparse_x, sparse_class, decision=sparse_decision)
    if (any(dense_class /= sparse_class)) error stop "sparse C-SVC class parity failed"
    if (maxval(abs(dense_decision - sparse_decision)) > 1.0e-10_dp) then
        error stop "sparse C-SVC decision parity failed"
    end if
    options%scale = .true.
    call svm_fit_classification_csr(sparse_x, y, model, options)
    if (any(model%scale_mask)) error stop "sparse SVM must disable predictor scaling like e1071"
    call svm_predict_classification_csr(model, sparse_x, class, decision=decision)
    if (any(class /= y)) error stop "sparse unscaled C-SVC parity failed"
    options%scale = .false.

    options%probability = .true.
    options%probability_seed = 123
    call svm_fit_classification(x, y, model, options)
    call svm_predict_classification(model, x, class, probability=probability)
    if (maxval(abs(sum(probability, dim=2) - 1.0_dp)) > 1.0e-8_dp) error stop "SVC probability coupling failed"
    if (.not. all(ieee_is_finite(probability))) error stop "SVC probability calibration failed"
    options%probability = .false.

    options%svm_type = svm_nu_classification
    options%nu = 0.4_dp
    call svm_fit_classification(x, y, model, options)
    call svm_predict_classification(model, x, class)
    if (any(class /= y)) error stop "nu-SVC classification failed"

    do i = 1, 10
        xr(i, 1) = real(i - 1, dp) / 3.0_dp
        yr(i) = 1.0_dp + 2.0_dp * xr(i, 1)
    end do
    options = svm_options(svm_type=svm_eps_regression, kernel=svm_linear, cost=10.0_dp, epsilon=0.01_dp, scale=.false.)
    call svm_fit_regression(xr, yr, dense_model, options)
    call svm_predict_regression(dense_model, xr, prediction)
    if (sqrt(sum((prediction - yr)**2) / real(size(yr), dp)) > 0.08_dp) error stop "epsilon-SVR failed"
    dense_prediction = prediction
    call dense_to_csr(xr, sparse_x)
    call svm_fit_regression_csr(sparse_x, yr, sparse_model, options)
    call svm_predict_regression_csr(sparse_model, sparse_x, sparse_prediction)
    if (maxval(abs(dense_prediction - sparse_prediction)) > 1.0e-10_dp) then
        error stop "sparse epsilon-SVR decision parity failed"
    end if
    options%scale = .true.
    call svm_fit_regression_csr(sparse_x, yr, model, options)
    if (any(model%scale_mask)) error stop "sparse SVR must disable predictor scaling like e1071"
    if (abs(model%y_center) > 0.0_dp .or. abs(model%y_scale - 1.0_dp) > 0.0_dp) then
        error stop "sparse SVR must disable response scaling like e1071"
    end if
    call svm_predict_regression_csr(model, sparse_x, prediction)
    if (sqrt(sum((prediction - yr)**2) / real(size(yr), dp)) > 0.08_dp) error stop "sparse epsilon-SVR failed"
    options%scale = .false.

    options%probability = .true.
    options%probability_seed = 321
    call svm_fit_regression(xr, yr, model, options)
    if (.not. model%probability_fitted) error stop "SVR probability flag failed"
    if (.not. ieee_is_finite(model%svr_probability_sigma)) error stop "SVR probability scale failed"
    if (model%svr_probability_sigma < 0.0_dp) error stop "SVR probability scale is negative"
    options%probability = .false.

    options%svm_type = svm_nu_regression
    options%nu = 0.5_dp
    call svm_fit_regression(xr, yr, model, options)
    call svm_predict_regression(model, xr, prediction)
    if (.not. all(ieee_is_finite(prediction))) error stop "nu-SVR produced nonfinite values"
    if (sqrt(sum((prediction - yr)**2) / real(size(yr), dp)) > 0.2_dp) error stop "nu-SVR failed"

    one_x = reshape([-0.1_dp, 0.0_dp, 0.1_dp, 0.0_dp, 0.0_dp, -0.1_dp, &
                     0.1_dp, 0.1_dp, -0.1_dp, -0.1_dp], [5, 2], order=[2, 1])
    one_query(1, :) = [0.0_dp, 0.0_dp]
    one_query(2, :) = [5.0_dp, 5.0_dp]
    options = svm_options(svm_type=svm_one_classification, kernel=svm_radial, gamma=1.0_dp, nu=0.2_dp, scale=.false.)
    call svm_fit_one_class(one_x, dense_model, options)
    call svm_predict_one_class(dense_model, one_query, accepted)
    if (accepted(2)) error stop "one-class SVM failed to reject distant point"
    call dense_to_csr(one_x, sparse_x)
    call svm_fit_one_class_csr(sparse_x, sparse_model, options)
    call dense_to_csr(one_query, sparse_x)
    call svm_predict_one_class_csr(sparse_model, sparse_x, accepted)
    if (accepted(2)) error stop "sparse one-class SVM failed to reject distant point"

    options = svm_options(svm_type=svm_c_classification, kernel=svm_radial, gamma=0.5_dp, scale=.false.)
    call svm_cross_validate_classification(x, y, 3, cv_class, options)
    if (count(cv_class == y) < 5) error stop "SVM classification CV failed"
    call tune_svm_classification(x, y, [0.5_dp, 1.0_dp], [0.25_dp, 0.5_dp], 3, tune_fit, options)
    if (tune_fit%best_index < 1 .or. .not. ieee_is_finite(tune_fit%best_performance)) error stop "tune.svm failed"

    options = svm_options(svm_type=svm_eps_regression, kernel=svm_linear, cost=10.0_dp, epsilon=0.01_dp, scale=.false.)
    call svm_cross_validate_regression(xr, yr, 5, cv_prediction, options)
    if (.not. all(ieee_is_finite(cv_prediction))) error stop "SVM regression CV failed"

    print '(a)', "test_svm: PASS"
end program test_svm
