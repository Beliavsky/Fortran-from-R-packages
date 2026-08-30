program test_svm_io
    use e1071, only: dp, svm_model, svm_options, svm_c_classification, svm_eps_regression, &
                     svm_one_classification, svm_linear, svm_radial, svm_fit_classification, &
                     svm_fit_regression, svm_fit_one_class, svm_predict_classification, &
                     svm_predict_regression, svm_predict_one_class, svm_write_libsvm, svm_read_libsvm
    implicit none

    call test_external_libsvm_import()
    call test_multiclass_round_trip()
    call test_probability_round_trip()
    call test_regression_round_trip()
    call test_one_class_round_trip()
    print '(a)', 'test_svm_io: PASS'

contains

    subroutine test_external_libsvm_import()
        type(svm_model) :: model
        real(dp) :: x(9, 2)
        integer, allocatable :: prediction(:)
        integer :: expected(9)

        x = reshape([ &
            -2.0_dp, -2.0_dp, -2.2_dp, -1.8_dp, -1.8_dp, -2.2_dp, &
             2.0_dp, -2.0_dp,  2.2_dp, -1.8_dp,  1.8_dp, -2.2_dp, &
             0.0_dp,  2.0_dp,  0.2_dp,  2.2_dp, -0.2_dp,  1.8_dp], [9, 2], order=[2, 1])
        expected = [10, 10, 10, 20, 20, 20, 30, 30, 30]
        call svm_read_libsvm('test/data/libsvm_multiclass.model', model, n_features=2)
        call svm_predict_classification(model, x, prediction)
        if (any(prediction /= expected)) error stop 'external LIBSVM model import mismatch'
    end subroutine test_external_libsvm_import

    subroutine test_multiclass_round_trip()
        type(svm_model) :: fitted
        type(svm_model) :: loaded
        type(svm_options) :: options
        real(dp) :: x(9, 2)
        integer :: y(9)
        integer, allocatable :: class_before(:)
        integer, allocatable :: class_after(:)
        real(dp), allocatable :: decision_before(:, :)
        real(dp), allocatable :: decision_after(:, :)

        x = reshape([ &
            -2.0_dp, -2.0_dp, -2.2_dp, -1.8_dp, -1.8_dp, -2.2_dp, &
             2.0_dp, -2.0_dp,  2.2_dp, -1.8_dp,  1.8_dp, -2.2_dp, &
             0.0_dp,  2.0_dp,  0.2_dp,  2.2_dp, -0.2_dp,  1.8_dp], [9, 2], order=[2, 1])
        y = [10, 10, 10, 20, 20, 20, 30, 30, 30]
        options = svm_options(svm_type=svm_c_classification, kernel=svm_linear, scale=.false.)
        call svm_fit_classification(x, y, fitted, options)
        call svm_predict_classification(fitted, x, class_before, decision_before)
        call svm_write_libsvm(fitted, 'test_fortran_multiclass.model')
        call svm_read_libsvm('test_fortran_multiclass.model', loaded, n_features=2)
        call svm_predict_classification(loaded, x, class_after, decision_after)
        if (any(class_before /= class_after)) error stop 'multiclass serialization class mismatch'
        if (maxval(abs(decision_before - decision_after)) > 5.0e-7_dp) then
            error stop 'multiclass serialization decision mismatch'
        end if
        call delete_file('test_fortran_multiclass.model')
    end subroutine test_multiclass_round_trip

    subroutine test_probability_round_trip()
        type(svm_model) :: fitted
        type(svm_model) :: loaded
        type(svm_options) :: options
        real(dp) :: x(12, 2)
        integer :: y(12)
        integer, allocatable :: before(:)
        integer, allocatable :: after(:)
        real(dp), allocatable :: probability_before(:, :)
        real(dp), allocatable :: probability_after(:, :)
        integer :: i

        do i = 1, 6
            x(i, :) = [-2.0_dp - 0.1_dp * real(i, dp), -1.0_dp + 0.05_dp * real(i, dp)]
            x(i + 6, :) = [2.0_dp + 0.1_dp * real(i, dp), 1.0_dp - 0.05_dp * real(i, dp)]
        end do
        y(:6) = -1
        y(7:) = 1
        options = svm_options(svm_type=svm_c_classification, kernel=svm_radial, gamma=0.4_dp, &
                              scale=.false., probability=.true., probability_seed=91)
        call svm_fit_classification(x, y, fitted, options)
        call svm_predict_classification(fitted, x, before, probability=probability_before)
        call svm_write_libsvm(fitted, 'test_fortran_probability.model')
        call svm_read_libsvm('test_fortran_probability.model', loaded, n_features=2)
        call svm_predict_classification(loaded, x, after, probability=probability_after)
        if (any(before /= after)) error stop 'probability serialization class mismatch'
        if (maxval(abs(probability_before - probability_after)) > 2.0e-7_dp) then
            error stop 'probability serialization value mismatch'
        end if
        call delete_file('test_fortran_probability.model')
    end subroutine test_probability_round_trip

    subroutine test_regression_round_trip()
        type(svm_model) :: fitted
        type(svm_model) :: loaded
        type(svm_options) :: options
        real(dp) :: x(8, 3)
        real(dp) :: y(8)
        real(dp), allocatable :: before(:)
        real(dp), allocatable :: after(:)
        logical :: scale_mask(3)
        integer :: i

        scale_mask = [.true., .true., .false.]
        do i = 1, 8
            x(i, 1) = real(i - 4, dp)
            x(i, 2) = real(modulo(i, 3) - 1, dp)
            x(i, 3) = 0.0_dp
        end do
        y = 0.75_dp * x(:, 1) - 0.4_dp * x(:, 2) + 2.0_dp
        options = svm_options(svm_type=svm_eps_regression, kernel=svm_radial, gamma=0.3_dp, &
                              cost=5.0_dp, epsilon=0.01_dp, scale=.true., probability=.true.)
        call svm_fit_regression(x, y, fitted, options, scale_mask)
        call svm_predict_regression(fitted, x, before)
        call svm_write_libsvm(fitted, 'test_fortran_svr.model', 'test_fortran_svr.scale', 'test_fortran_svr.yscale')
        call svm_read_libsvm('test_fortran_svr.model', loaded, 'test_fortran_svr.scale', &
                             'test_fortran_svr.yscale')
        call svm_predict_regression(loaded, x, after)
        if (maxval(abs(before - after)) > 2.0e-6_dp) error stop 'SVR serialization prediction mismatch'
        if (.not. loaded%probability_fitted) error stop 'SVR probability scale was not restored'
        if (abs(loaded%svr_probability_sigma - fitted%svr_probability_sigma) > 1.0e-12_dp) then
            error stop 'SVR probability scale mismatch'
        end if
        call delete_file('test_fortran_svr.model')
        call delete_file('test_fortran_svr.scale')
        call delete_file('test_fortran_svr.yscale')
    end subroutine test_regression_round_trip

    subroutine test_one_class_round_trip()
        type(svm_model) :: fitted
        type(svm_model) :: loaded
        type(svm_options) :: options
        real(dp) :: x(6, 2)
        logical, allocatable :: before(:)
        logical, allocatable :: after(:)
        real(dp), allocatable :: decision_before(:)
        real(dp), allocatable :: decision_after(:)

        x = reshape([0.0_dp, 0.0_dp, 0.1_dp, 0.0_dp, -0.1_dp, 0.0_dp, &
                     0.0_dp, 0.1_dp, 0.0_dp, -0.1_dp, 0.05_dp, 0.05_dp], [6, 2], order=[2, 1])
        options = svm_options(svm_type=svm_one_classification, kernel=svm_radial, gamma=1.0_dp, scale=.false.)
        call svm_fit_one_class(x, fitted, options)
        call svm_predict_one_class(fitted, x, before, decision_before)
        call svm_write_libsvm(fitted, 'test_fortran_oneclass.model')
        call svm_read_libsvm('test_fortran_oneclass.model', loaded, n_features=2)
        call svm_predict_one_class(loaded, x, after, decision_after)
        if (any(before .neqv. after)) error stop 'one-class serialization classification mismatch'
        if (maxval(abs(decision_before - decision_after)) > 5.0e-7_dp) then
            error stop 'one-class serialization decision mismatch'
        end if
        call delete_file('test_fortran_oneclass.model')
    end subroutine test_one_class_round_trip

    subroutine delete_file(filename)
        character(len=*), intent(in) :: filename !! Temporary test file removed after serialization checks complete.
        integer :: unit
        integer :: ios

        open(newunit=unit, file=filename, status='old', action='readwrite', iostat=ios)
        if (ios == 0) close(unit, status='delete')
    end subroutine delete_file

end program test_svm_io
