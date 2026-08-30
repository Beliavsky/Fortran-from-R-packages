program test_proxy_all_methods
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use proxy, only: dp, proxy_dist_auto, proxy_dist_cross, proxy_dist_pairwise
    implicit none

    character(len=24), parameter :: numeric_methods(*) = [character(len=24) :: &
        'Euclidean', 'Mahalanobis', 'Bhjattacharyya', 'Manhattan', 'supremum', &
        'Minkowski', 'Canberra', 'WaveHedges', 'divergence', 'KullbackLeibler', &
        'BrayCurtis', 'Soergel', 'Podani', 'Chord', 'Geodesic', 'Whittaker', &
        'Hellinger', 'fJaccard', 'cosine', 'angular', 'eJaccard', 'eDice', &
        'correlation', 'Gower', 'mutual']
    character(len=24), parameter :: binary_methods(*) = [character(len=24) :: &
        'Jaccard', 'Kulczynski1', 'Kulczynski2', 'Mountford', 'FagerMcGowan', &
        'RusselRao', 'simple matching', 'Hamman', 'Faith', 'RogersTanimoto', &
        'Dice', 'Phi', 'Stiles', 'Michael', 'MozleyMargalef', 'Yule', 'Yule2', &
        'Ochiai', 'Simpson', 'Braun-Blanquet']
    character(len=24), parameter :: nominal_methods(*) = [character(len=24) :: &
        'Chi-squared', 'Phi-squared', 'Tschuprow', 'Cramer', 'Pearson']
    real(dp) :: x(6, 4)
    real(dp) :: b(6, 6)
    real(dp) :: n(6, 6)
    integer :: i

    x(1, :) = [0.10_dp, 0.20_dp, 0.30_dp, 0.40_dp]
    x(2, :) = [0.20_dp, 0.40_dp, 0.10_dp, 0.30_dp]
    x(3, :) = [0.80_dp, 0.10_dp, 0.60_dp, 0.20_dp]
    x(4, :) = [0.60_dp, 0.90_dp, 0.20_dp, 0.70_dp]
    x(5, :) = [0.35_dp, 0.55_dp, 0.75_dp, 0.25_dp]
    x(6, :) = [0.90_dp, 0.65_dp, 0.45_dp, 0.85_dp]

    b(1, :) = [1, 0, 1, 0, 1, 0]
    b(2, :) = [1, 1, 0, 0, 1, 0]
    b(3, :) = [0, 1, 1, 0, 0, 1]
    b(4, :) = [1, 0, 0, 1, 0, 1]
    b(5, :) = [0, 1, 0, 1, 1, 0]
    b(6, :) = [1, 1, 1, 0, 0, 0]

    n(1, :) = [1, 1, 2, 2, 3, 3]
    n(2, :) = [1, 2, 1, 2, 3, 2]
    n(3, :) = [2, 1, 2, 3, 1, 3]
    n(4, :) = [2, 2, 3, 1, 2, 1]
    n(5, :) = [3, 3, 1, 2, 1, 2]
    n(6, :) = [3, 2, 3, 1, 2, 3]

    do i = 1, size(numeric_methods)
        if (trim(numeric_methods(i)) == 'Minkowski') then
            call check_method(x, numeric_methods(i), [1.5_dp])
        else
            call check_method(x, numeric_methods(i))
        end if
    end do
    do i = 1, size(binary_methods)
        call check_method(b, binary_methods(i))
    end do
    do i = 1, size(nominal_methods)
        call check_method(n, nominal_methods(i))
    end do

contains

    subroutine check_method(data, method, p)
        real(dp), intent(in) :: data(:, :) !! Deterministic test matrix appropriate for the selected method family.
        character(len=*), intent(in) :: method !! Canonical proxy method name being checked for auto/cross/pairwise consistency.
        real(dp), intent(in), optional :: p(:) !! Optional method parameters; used here for the Minkowski exponent.
        real(dp), allocatable :: auto(:, :)
        real(dp), allocatable :: cross(:, :)
        real(dp), allocatable :: pair(:)
        integer :: i
        integer :: j
        integer :: status

        if (present(p)) then
            call proxy_dist_auto(data, method, auto, p=p, status=status)
            call assert_true(status == 0, trim(method)//' auto status')
            call proxy_dist_cross(data, data, method, cross, p=p, status=status)
            call assert_true(status == 0, trim(method)//' cross status')
            call proxy_dist_pairwise(data, data, method, pair, p=p, status=status)
        else
            call proxy_dist_auto(data, method, auto, status=status)
            call assert_true(status == 0, trim(method)//' auto status')
            call proxy_dist_cross(data, data, method, cross, status=status)
            call assert_true(status == 0, trim(method)//' cross status')
            call proxy_dist_pairwise(data, data, method, pair, status=status)
        end if
        call assert_true(status == 0, trim(method)//' pairwise status')
        call assert_true(size(auto, 1) == size(data, 1), trim(method)//' auto shape')
        do j = 1, size(data, 1)
            do i = 1, size(data, 1)
                call assert_true(same_value(auto(i, j), cross(i, j)), trim(method)//' auto-cross parity')
            end do
            call assert_true(same_value(auto(j, j), pair(j)), trim(method)//' diagonal-pairwise parity')
        end do
    end subroutine check_method

    pure function same_value(a, b) result(equal)
        real(dp), intent(in) :: a !! First floating-point result in a parity comparison; NaN is considered equal to NaN.
        real(dp), intent(in) :: b !! Second result; finite values use a scaled absolute tolerance.
        !! Infinities are compared by sign/value.
        logical :: equal

        if (ieee_is_nan(a) .or. ieee_is_nan(b)) then
            equal = ieee_is_nan(a) .and. ieee_is_nan(b)
        else if (abs(a) >= huge(1.0_dp) / 2.0_dp .or. abs(b) >= huge(1.0_dp) / 2.0_dp) then
            equal = (a > 0.0_dp .and. b > 0.0_dp) .or. (a < 0.0_dp .and. b < 0.0_dp)
        else
            equal = abs(a - b) <= 1.0e-11_dp * max(1.0_dp, abs(a), abs(b))
        end if
    end function same_value

    subroutine assert_true(condition, label)
        logical, intent(in) :: condition !! Boolean condition that must be true for this all-method consistency test.
        character(len=*), intent(in) :: label !! Method-specific diagnostic label printed if the assertion fails.

        if (.not. condition) then
            write (*, '(a)') 'FAIL '//trim(label)
            error stop 1
        end if
    end subroutine assert_true

end program test_proxy_all_methods
