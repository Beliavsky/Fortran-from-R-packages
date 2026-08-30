program svm_io_example
    use e1071, only: dp, svm_model, svm_options, svm_c_classification, svm_linear, &
                     svm_fit_classification, svm_predict_classification, svm_write_libsvm, svm_read_libsvm
    implicit none

    real(dp) :: x(6, 2)
    integer :: y(6)
    integer, allocatable :: prediction(:)
    type(svm_options) :: options
    type(svm_model) :: fitted
    type(svm_model) :: loaded
    integer :: unit
    integer :: ios

    x(1, :) = [-2.0_dp, -1.0_dp]
    x(2, :) = [-1.5_dp, -1.2_dp]
    x(3, :) = [-1.8_dp, -0.8_dp]
    x(4, :) = [1.8_dp, 0.9_dp]
    x(5, :) = [2.1_dp, 1.2_dp]
    x(6, :) = [1.6_dp, 0.8_dp]
    y = [-1, -1, -1, 1, 1, 1]

    options = svm_options(svm_type=svm_c_classification, kernel=svm_linear, scale=.false.)
    call svm_fit_classification(x, y, fitted, options)
    call svm_write_libsvm(fitted, 'svm_io_example.model')
    call svm_read_libsvm('svm_io_example.model', loaded, n_features=2)
    call svm_predict_classification(loaded, x, prediction)
    print '(a,*(1x,i0))', 'reloaded predictions:', prediction

    open(newunit=unit, file='svm_io_example.model', status='old', action='readwrite', iostat=ios)
    if (ios == 0) close(unit, status='delete')
end program svm_io_example
