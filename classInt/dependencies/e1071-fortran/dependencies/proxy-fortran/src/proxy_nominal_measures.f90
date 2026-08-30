module proxy_nominal_measures
    use proxy_kinds, only: dp
    use proxy_ieee, only: proxy_nan, proxy_is_missing
    implicit none
    private

    public :: nominal_similarity

contains

    pure function nominal_similarity(x, y, method) result(similarity)
        real(dp), intent(in) :: x(:) !! First nominal variable encoded by numeric category labels; NaN entries are omitted pairwise.
        real(dp), intent(in) :: y(:) !! Second nominal variable encoded by numeric category labels; must align componentwise with
        !! `x`.
        character(len=*), intent(in) :: method !! Association measure: `Chi-squared`, `Phi-squared`, `Tschuprow`, `Cramer`, or
        !! `Pearson`/`contingency`.
        real(dp) :: similarity
        integer, allocatable :: table(:, :)
        real(dp) :: chi
        integer :: nr
        integer :: nc
        integer :: n
        character(len=:), allocatable :: key

        call contingency_table(x, y, table, nr, nc, n)
        if (n == 0 .or. nr == 0 .or. nc == 0) then
            similarity = proxy_nan()
            return
        end if
        chi = chi_square(table)
        key = normalize(method)
        select case (key)
        case ('chisquared')
            similarity = chi
        case ('phisquared')
            similarity = chi / real(n, dp)
        case ('tschuprow')
            if (nr <= 1 .or. nc <= 1) then
                similarity = proxy_nan()
            else
                similarity = sqrt(chi / real(n, dp) / sqrt(real((nr - 1) * (nc - 1), dp)))
            end if
        case ('cramer')
            if (min(nr - 1, nc - 1) <= 0) then
                similarity = proxy_nan()
            else
                similarity = sqrt(chi / real(n, dp) / real(min(nr - 1, nc - 1), dp))
            end if
        case ('pearson', 'contingency')
            similarity = sqrt(chi / (real(n, dp) + chi))
        case default
            similarity = proxy_nan()
        end select
    end function nominal_similarity

    pure subroutine contingency_table(x, y, table, nr, nc, nvalid)
        real(dp), intent(in) :: x(:) !! First numeric category-label vector; NaN pairs are ignored.
        real(dp), intent(in) :: y(:) !! Second category-label vector aligned with `x`.
        integer, allocatable, intent(out) :: table(:, :) !! Allocated observed contingency table over category levels present
        !! among valid pairs.
        integer, intent(out) :: nr !! Number of distinct valid category labels found in `x`.
        integer, intent(out) :: nc !! Number of distinct valid category labels found in `y`.
        integer, intent(out) :: nvalid !! Number of nonmissing aligned observations entered into the table.
        real(dp), allocatable :: xr(:)
        real(dp), allocatable :: yc(:)
        integer :: i
        integer :: ir
        integer :: ic

        allocate(xr(size(x)), yc(size(y)))
        nr = 0
        nc = 0
        nvalid = 0
        do i = 1, min(size(x), size(y))
            if (proxy_is_missing(x(i)) .or. proxy_is_missing(y(i))) cycle
            call find_or_append(xr, nr, x(i), ir)
            call find_or_append(yc, nc, y(i), ic)
            nvalid = nvalid + 1
        end do
        allocate(table(nr, nc))
        table = 0
        do i = 1, min(size(x), size(y))
            if (proxy_is_missing(x(i)) .or. proxy_is_missing(y(i))) cycle
            ir = find_value(xr, nr, x(i))
            ic = find_value(yc, nc, y(i))
            table(ir, ic) = table(ir, ic) + 1
        end do
    end subroutine contingency_table

    pure subroutine find_or_append(values, n, value, index)
        real(dp), intent(inout) :: values(:) !! Workspace containing distinct category labels in positions `1:n`.
        integer, intent(inout) :: n !! Number of currently stored distinct labels; incremented when `value` is new.
        real(dp), intent(in) :: value !! Category label to locate or append using exact numeric equality semantics.
        integer, intent(out) :: index !! One-based position of `value` in the distinct-label workspace after the operation.

        index = find_value(values, n, value)
        if (index == 0) then
            n = n + 1
            values(n) = value
            index = n
        end if
    end subroutine find_or_append

    pure function find_value(values, n, value) result(index)
        real(dp), intent(in) :: values(:) !! Array whose first `n` elements hold distinct numeric category labels.
        integer, intent(in) :: n !! Number of valid labels at the start of `values`.
        real(dp), intent(in) :: value !! Numeric category label to locate using exact equality semantics.
        integer :: index
        integer :: i

        index = 0
        do i = 1, n
            if (values(i) <= value .and. values(i) >= value) then
                index = i
                return
            end if
        end do
    end function find_value

    pure function chi_square(table) result(value)
        integer, intent(in) :: table(:, :) !! Observed contingency counts; dimensions are the numbers of category levels in the
        !! two variables.
        real(dp) :: value
        integer, allocatable :: rows(:)
        integer, allocatable :: cols(:)
        real(dp) :: expected
        integer :: i
        integer :: j
        integer :: n

        allocate(rows(size(table, 1)), cols(size(table, 2)))
        rows = sum(table, dim=2)
        cols = sum(table, dim=1)
        n = sum(rows)
        value = 0.0_dp
        if (n <= 0) return
        do i = 1, size(table, 1)
            do j = 1, size(table, 2)
                expected = real(rows(i) * cols(j), dp) / real(n, dp)
                if (expected > tiny(1.0_dp)) then
                    value = value + (real(table(i, j), dp) - expected)**2 / expected
                end if
            end do
        end do
    end function chi_square

    pure function normalize(text) result(key)
        character(len=*), intent(in) :: text !! Measure name; ASCII case, spaces, underscores, and hyphens are normalized for
        !! lookup.
        character(len=:), allocatable :: key
        character(len=len_trim(text)) :: work
        integer :: c
        integer :: i
        integer :: j

        work = ''
        j = 0
        do i = 1, len_trim(text)
            c = iachar(text(i:i))
            if (c >= iachar('A') .and. c <= iachar('Z')) c = c + 32
            if ((c >= iachar('a') .and. c <= iachar('z')) .or. (c >= iachar('0') .and. c <= iachar('9'))) then
                j = j + 1
                work(j:j) = achar(c)
            end if
        end do
        key = work(:j)
    end function normalize

end module proxy_nominal_measures
