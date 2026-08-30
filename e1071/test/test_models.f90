program test_models
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use e1071
    implicit none

    real(dp) :: x(8, 2)
    real(dp) :: xn(8, 1)
    integer :: xc(8, 1)
    integer :: y(8)
    integer :: lca_x(12, 3)
    real(dp) :: ica_x(20, 2)
    real(dp) :: iw(2, 2)
    integer :: table(2, 2)
    integer, allocatable :: cls(:)
    integer, allocatable :: mapping(:)
    integer, allocatable :: controls(:)
    real(dp), allocatable :: pred(:)
    real(dp) :: match_x(4, 2)
    integer :: match_types(2)
    logical :: is_case(4)
    logical :: is_control(4)
    type(naive_bayes_model) :: nb
    type(gknn_model) :: knn
    type(lca_model) :: lca_fit
    type(ica_model) :: ica_fit_result
    type(class_agreement_result) :: agreement
    type(agreement_matrix_result) :: agreement_matrix
    type(bclust_model) :: bag
    type(rng_state) :: rng
    integer :: i

    x = reshape([0.0_dp, 0.0_dp, 0.2_dp, 0.1_dp, -0.1_dp, 0.2_dp, 0.1_dp, -0.2_dp, &
                 5.0_dp, 5.0_dp, 5.2_dp, 5.1_dp, 4.9_dp, 5.2_dp, 5.1_dp, 4.8_dp], [8, 2], order=[2, 1])
    y = [1, 1, 1, 1, 2, 2, 2, 2]
    xn(:, 1) = x(:, 1)
    xc(:, 1) = [1, 1, 1, 1, 2, 2, 2, 2]
    call naive_bayes_fit(y, nb, x_numeric=xn, x_categorical=xc, categorical_levels=[2])
    call naive_bayes_predict(nb, cls, x_numeric=xn, x_categorical=xc)
    if (count(cls == y) < 8) error stop "naive Bayes training predictions failed"

    call gknn_fit_classification(x, y, knn, k=1, scale=.false.)
    call gknn_predict_classification(knn, x, cls)
    if (count(cls == y) /= 8) error stop "gknn classification failed"
    call gknn_fit_regression(x, real(y, dp), knn, k=2, scale=.false.)
    call gknn_predict_regression(knn, x, pred)
    if (.not. all(ieee_is_finite(pred))) error stop "gknn regression failed"

    lca_x = reshape([0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, &
                     1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1], [12, 3], order=[2, 1])
    call rng_seed(rng, 12345)
    call lca_fit_data(lca_x, 2, 100, rng, lca_fit)
    if (.not. ieee_is_finite(lca_fit%log_likelihood)) error stop "LCA likelihood failed"
    if (abs(sum(lca_fit%class_probability) - 1.0_dp) > 1.0e-10_dp) error stop "LCA class probabilities failed"

    do i = 1, 20
        ica_x(i, 1) = sin(0.3_dp * real(i, dp))
        ica_x(i, 2) = cos(0.7_dp * real(i, dp))
    end do
    iw = 0.0_dp
    iw(1, 1) = 1.0_dp
    iw(2, 2) = 1.0_dp
    call rng_seed(rng, 6789)
    call ica_fit(ica_x, 0.01_dp, rng, ica_fit_result, epochs=5, ncomp=2, initial_weights=iw)
    if (size(ica_fit_result%projection, 1) /= 20) error stop "ICA projection shape failed"
    if (.not. all(ieee_is_finite(ica_fit_result%weights))) error stop "ICA weights failed"

    table = reshape([30, 2, 3, 25], [2, 2])
    agreement = class_agreement(table)
    if (agreement%diagonal <= 0.8_dp) error stop "classAgreement failed"
    call match_classes(table, mapping, method=match_exact)
    if (any(mapping /= [1, 2])) error stop "matchClasses failed"

    call compare_matched_classes(reshape([1, 1, 2, 2, 2, 2, 1, 1], [4, 2]), agreement_matrix)
    if (agreement_matrix%diagonal(1, 2) < 0.99_dp) error stop "compareMatchedClasses failed"

    match_x = reshape([0.0_dp, 0.1_dp, 1.0_dp, 0.9_dp, 1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp], [4, 2])
    match_types = [proxy_gower_metric, proxy_gower_factor]
    is_case = [.true., .false., .true., .false.]
    is_control = [.false., .true., .false., .true.]
    call match_controls_gower(match_x, match_types, is_case, is_control, controls)
    if (any(controls /= [2, 4])) error stop "mixed Gower matchControls failed"

    call rng_seed(rng, 42)
    call bclust_fit(x, 2, rng, bag, iter_base=3, base_centers=2, maxcluster=2, resample=.false.)
    if (size(bag%centers, 1) /= 2) error stop "bclust center count failed"
    if (size(bag%cluster) /= 8) error stop "bclust labels failed"

    print '(a)', "test_models: PASS"
end program test_models
