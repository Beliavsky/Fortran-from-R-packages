! SPDX-License-Identifier: BSD-3-Clause
module metrics_information_retrieval
    use metrics_kinds, only : dp, metrics_success, metrics_invalid_size, metrics_invalid_argument
    use metrics_utils, only : quiet_nan
    implicit none
    private

    type, public :: integer_vector
        integer, allocatable :: values(:)
    end type integer_vector

    type, public :: real_vector
        real(dp), allocatable :: values(:)
    end type real_vector

    type, public :: string_vector
        character(len=:), allocatable :: values(:)
    end type string_vector

    public :: f1, apk, mapk

    interface f1
        module procedure f1_integer
        module procedure f1_real
        module procedure f1_character
    end interface f1

    interface apk
        module procedure apk_integer
        module procedure apk_real
        module procedure apk_character
    end interface apk

    interface mapk
        module procedure mapk_integer
        module procedure mapk_real
        module procedure mapk_character
    end interface mapk

contains

    real(dp) function f1_integer(actual, predicted) result(value)
        integer, intent(in) :: actual(:), predicted(:)
        integer :: tp, fp, fn
        call counts_integer(actual, predicted, tp, fp, fn)
        value = f1_from_counts(tp, fp, fn)
    end function f1_integer

    real(dp) function f1_real(actual, predicted) result(value)
        real(dp), intent(in) :: actual(:), predicted(:)
        integer :: tp, fp, fn
        call counts_real(actual, predicted, tp, fp, fn)
        value = f1_from_counts(tp, fp, fn)
    end function f1_real

    real(dp) function f1_character(actual, predicted) result(value)
        character(len=*), intent(in) :: actual(:), predicted(:)
        integer :: tp, fp, fn
        call counts_character(actual, predicted, tp, fp, fn)
        value = f1_from_counts(tp, fp, fn)
    end function f1_character

    real(dp) function apk_integer(k, actual, predicted, stat) result(value)
        integer, intent(in) :: k, actual(:), predicted(:)
        integer, intent(out), optional :: stat
        integer :: i, limit, count_hits
        real(dp) :: score

        if (k < 0) then
            value = quiet_nan()
            if (present(stat)) stat = metrics_invalid_argument
            return
        end if
        score = 0.0_dp
        count_hits = 0
        limit = min(k, size(predicted))
        do i = 1, limit
            if (any(actual == predicted(i)) .and. .not. any(predicted(:i - 1) == predicted(i))) then
                count_hits = count_hits + 1
                score = score + real(count_hits, dp) / real(i, dp)
            end if
        end do
        if (min(size(actual), k) == 0) then
            value = quiet_nan()
        else
            value = score / real(min(size(actual), k), dp)
        end if
        if (present(stat)) stat = metrics_success
    end function apk_integer

    real(dp) function apk_real(k, actual, predicted, stat) result(value)
        integer, intent(in) :: k
        real(dp), intent(in) :: actual(:), predicted(:)
        integer, intent(out), optional :: stat
        integer :: i, limit, count_hits
        real(dp) :: score

        if (k < 0) then
            value = quiet_nan()
            if (present(stat)) stat = metrics_invalid_argument
            return
        end if
        score = 0.0_dp
        count_hits = 0
        limit = min(k, size(predicted))
        do i = 1, limit
            if (any(actual == predicted(i)) .and. .not. any(predicted(:i - 1) == predicted(i))) then
                count_hits = count_hits + 1
                score = score + real(count_hits, dp) / real(i, dp)
            end if
        end do
        if (min(size(actual), k) == 0) then
            value = quiet_nan()
        else
            value = score / real(min(size(actual), k), dp)
        end if
        if (present(stat)) stat = metrics_success
    end function apk_real

    real(dp) function apk_character(k, actual, predicted, stat) result(value)
        integer, intent(in) :: k
        character(len=*), intent(in) :: actual(:), predicted(:)
        integer, intent(out), optional :: stat
        integer :: i, limit, count_hits
        real(dp) :: score

        if (k < 0) then
            value = quiet_nan()
            if (present(stat)) stat = metrics_invalid_argument
            return
        end if
        score = 0.0_dp
        count_hits = 0
        limit = min(k, size(predicted))
        do i = 1, limit
            if (any(actual == predicted(i)) .and. .not. any(predicted(:i - 1) == predicted(i))) then
                count_hits = count_hits + 1
                score = score + real(count_hits, dp) / real(i, dp)
            end if
        end do
        if (min(size(actual), k) == 0) then
            value = quiet_nan()
        else
            value = score / real(min(size(actual), k), dp)
        end if
        if (present(stat)) stat = metrics_success
    end function apk_character

    real(dp) function mapk_integer(k, actual, predicted, stat) result(value)
        integer, intent(in) :: k
        type(integer_vector), intent(in) :: actual(:), predicted(:)
        integer, intent(out), optional :: stat
        integer :: i

        if (size(actual) == 0 .or. size(predicted) == 0) then
            value = 0.0_dp
            if (present(stat)) stat = metrics_success
            return
        end if
        if (size(actual) /= size(predicted)) then
            value = quiet_nan()
            if (present(stat)) stat = metrics_invalid_size
            return
        end if
        value = 0.0_dp
        do i = 1, size(actual)
            value = value + apk_integer(k, actual(i)%values, predicted(i)%values)
        end do
        value = value / real(size(actual), dp)
        if (present(stat)) stat = metrics_success
    end function mapk_integer

    real(dp) function mapk_real(k, actual, predicted, stat) result(value)
        integer, intent(in) :: k
        type(real_vector), intent(in) :: actual(:), predicted(:)
        integer, intent(out), optional :: stat
        integer :: i

        if (size(actual) == 0 .or. size(predicted) == 0) then
            value = 0.0_dp
            if (present(stat)) stat = metrics_success
            return
        end if
        if (size(actual) /= size(predicted)) then
            value = quiet_nan()
            if (present(stat)) stat = metrics_invalid_size
            return
        end if
        value = 0.0_dp
        do i = 1, size(actual)
            value = value + apk_real(k, actual(i)%values, predicted(i)%values)
        end do
        value = value / real(size(actual), dp)
        if (present(stat)) stat = metrics_success
    end function mapk_real

    real(dp) function mapk_character(k, actual, predicted, stat) result(value)
        integer, intent(in) :: k
        type(string_vector), intent(in) :: actual(:), predicted(:)
        integer, intent(out), optional :: stat
        integer :: i

        if (size(actual) == 0 .or. size(predicted) == 0) then
            value = 0.0_dp
            if (present(stat)) stat = metrics_success
            return
        end if
        if (size(actual) /= size(predicted)) then
            value = quiet_nan()
            if (present(stat)) stat = metrics_invalid_size
            return
        end if
        value = 0.0_dp
        do i = 1, size(actual)
            value = value + apk_character(k, actual(i)%values, predicted(i)%values)
        end do
        value = value / real(size(actual), dp)
        if (present(stat)) stat = metrics_success
    end function mapk_character

    pure real(dp) function f1_from_counts(tp, fp, fn) result(value)
        integer, intent(in) :: tp, fp, fn
        real(dp) :: precision_value, recall_value

        if (tp == 0) then
            value = 0.0_dp
        else
            precision_value = real(tp, dp) / real(tp + fp, dp)
            recall_value = real(tp, dp) / real(tp + fn, dp)
            value = 2.0_dp * precision_value * recall_value / (precision_value + recall_value)
        end if
    end function f1_from_counts

    subroutine counts_integer(actual, predicted, tp, fp, fn)
        integer, intent(in) :: actual(:), predicted(:)
        integer, intent(out) :: tp, fp, fn
        integer :: i

        tp = 0
        fp = 0
        fn = 0
        do i = 1, size(predicted)
            if (first_integer(predicted, i)) then
                if (any(actual == predicted(i))) then
                    tp = tp + 1
                else
                    fp = fp + 1
                end if
            end if
        end do
        do i = 1, size(actual)
            if (first_integer(actual, i) .and. .not. any(predicted == actual(i))) fn = fn + 1
        end do
    end subroutine counts_integer

    subroutine counts_real(actual, predicted, tp, fp, fn)
        real(dp), intent(in) :: actual(:), predicted(:)
        integer, intent(out) :: tp, fp, fn
        integer :: i

        tp = 0
        fp = 0
        fn = 0
        do i = 1, size(predicted)
            if (first_real(predicted, i)) then
                if (any(actual == predicted(i))) then
                    tp = tp + 1
                else
                    fp = fp + 1
                end if
            end if
        end do
        do i = 1, size(actual)
            if (first_real(actual, i) .and. .not. any(predicted == actual(i))) fn = fn + 1
        end do
    end subroutine counts_real

    subroutine counts_character(actual, predicted, tp, fp, fn)
        character(len=*), intent(in) :: actual(:), predicted(:)
        integer, intent(out) :: tp, fp, fn
        integer :: i

        tp = 0
        fp = 0
        fn = 0
        do i = 1, size(predicted)
            if (first_character(predicted, i)) then
                if (any(actual == predicted(i))) then
                    tp = tp + 1
                else
                    fp = fp + 1
                end if
            end if
        end do
        do i = 1, size(actual)
            if (first_character(actual, i) .and. .not. any(predicted == actual(i))) fn = fn + 1
        end do
    end subroutine counts_character

    logical function first_integer(x, i) result(first)
        integer, intent(in) :: x(:)
        integer, intent(in) :: i
        if (i == 1) then
            first = .true.
        else
            first = .not. any(x(:i - 1) == x(i))
        end if
    end function first_integer

    logical function first_real(x, i) result(first)
        real(dp), intent(in) :: x(:)
        integer, intent(in) :: i
        if (i == 1) then
            first = .true.
        else
            first = .not. any(x(:i - 1) == x(i))
        end if
    end function first_real

    logical function first_character(x, i) result(first)
        character(len=*), intent(in) :: x(:)
        integer, intent(in) :: i
        if (i == 1) then
            first = .true.
        else
            first = .not. any(x(:i - 1) == x(i))
        end if
    end function first_character

end module metrics_information_retrieval
