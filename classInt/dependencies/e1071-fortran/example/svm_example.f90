program svm_example
    use e1071
    implicit none

    real(dp) :: x(6, 2)
    integer :: y(6)
    integer, allocatable :: prediction(:)
    type(svm_options) :: options
    type(svm_model) :: model

    x = reshape([0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, &
                 3.0_dp, 3.0_dp, 3.0_dp, 4.0_dp, 4.0_dp, 3.0_dp], [6, 2], order=[2, 1])
    y = [-1, -1, -1, 1, 1, 1]

    options = svm_options(kernel=svm_radial, gamma=0.5_dp, scale=.false.)
    call svm_fit_classification(x, y, model, options)
    call svm_predict_classification(model, x, prediction)

    print '(a,*(1x,i0))', "training predictions:", prediction
end program svm_example
