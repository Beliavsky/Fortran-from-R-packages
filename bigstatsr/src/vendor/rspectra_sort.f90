! SPDX-License-Identifier: MPL-2.0
! Modern Fortran translation of the RSpectra computational interface.

module rspectra_sort
    use rspectra_kinds, only: dp
    implicit none
    private
    public :: sort_real_pairs, sort_complex_pairs

contains

    subroutine sort_real_pairs(values, vectors, which, sigma, shift_mode)
        real(dp), intent(inout) :: values(:)
        real(dp), intent(inout), optional :: vectors(:,:)
        character(len=*), intent(in) :: which
        real(dp), intent(in), optional :: sigma
        logical, intent(in), optional :: shift_mode
        integer :: i, j, n, best
        real(dp) :: tmp, sbest, scur
        real(dp), allocatable :: vtmp(:)
        logical :: shifted
        n = size(values)
        shifted = .false.
        if (present(shift_mode)) shifted = shift_mode
        if (present(vectors)) allocate(vtmp(size(vectors, 1)))
        if (trim(which) == 'BE' .and. .not. shifted) then
            call reorder_both_ends(values, vectors)
            return
        end if
        do i = 1, n - 1
            best = i
            sbest = real_score(values(i), which, sigma, shifted)
            do j = i + 1, n
                scur = real_score(values(j), which, sigma, shifted)
                if (scur < sbest) then
                    best = j
                    sbest = scur
                end if
            end do
            if (best /= i) then
                tmp = values(i)
                values(i) = values(best)
                values(best) = tmp
                if (present(vectors)) then
                    vtmp = vectors(:, i)
                    vectors(:, i) = vectors(:, best)
                    vectors(:, best) = vtmp
                end if
            end if
        end do
    end subroutine sort_real_pairs

    subroutine reorder_both_ends(values, vectors)
        real(dp), intent(inout) :: values(:)
        real(dp), intent(inout), optional :: vectors(:,:)
        integer :: i, j, n, best, hi, lo, src
        real(dp) :: tmp
        real(dp), allocatable :: vals(:), vecs(:,:)
        n = size(values)
        ! First sort algebraically from largest to smallest.
        do i = 1, n - 1
            best = i
            do j = i + 1, n
                if (values(j) > values(best)) best = j
            end do
            if (best /= i) then
                tmp = values(i)
                values(i) = values(best)
                values(best) = tmp
                if (present(vectors)) then
                    allocate(vals(size(vectors, 1)))
                    vals = vectors(:, i)
                    vectors(:, i) = vectors(:, best)
                    vectors(:, best) = vals
                    deallocate(vals)
                end if
            end if
        end do
        allocate(vals(n))
        if (present(vectors)) then
            allocate(vecs(size(vectors, 1), n))
            hi = 1
            lo = n
            do i = 1, n
                if (mod(i, 2) == 1) then
                    src = hi
                    hi = hi + 1
                else
                    src = lo
                    lo = lo - 1
                end if
                vals(i) = values(src)
                vecs(:, i) = vectors(:, src)
            end do
            values = vals
            vectors = vecs
        else
            hi = 1
            lo = n
            do i = 1, n
                if (mod(i, 2) == 1) then
                    src = hi
                    hi = hi + 1
                else
                    src = lo
                    lo = lo - 1
                end if
                vals(i) = values(src)
            end do
            values = vals
        end if
    end subroutine reorder_both_ends

    real(dp) function real_score(x, which, sigma, shifted) result(s)
        real(dp), intent(in) :: x
        character(len=*), intent(in) :: which
        real(dp), intent(in), optional :: sigma
        logical, intent(in) :: shifted
        real(dp) :: y
        y = x
        if (shifted .and. present(sigma)) then
            if (abs(x - sigma) > tiny(1.0_dp)) then
                y = 1.0_dp / (x - sigma)
            else
                y = huge(1.0_dp)
            end if
        end if
        select case (trim(which))
        case ('LM')
            s = -abs(y)
        case ('SM')
            s = abs(y)
        case ('LA', 'LR')
            s = -y
        case ('SA', 'SR')
            s = y
        case default
            s = -abs(y)
        end select
    end function real_score

    subroutine sort_complex_pairs(values, vectors, which, sigma, shift_mode)
        complex(dp), intent(inout) :: values(:)
        complex(dp), intent(inout), optional :: vectors(:,:)
        character(len=*), intent(in) :: which
        complex(dp), intent(in), optional :: sigma
        logical, intent(in), optional :: shift_mode
        integer :: i, j, n, best
        real(dp) :: sbest, scur
        complex(dp) :: tmp
        complex(dp), allocatable :: vtmp(:)
        logical :: shifted
        n = size(values)
        shifted = .false.
        if (present(shift_mode)) shifted = shift_mode
        if (present(vectors)) allocate(vtmp(size(vectors, 1)))
        do i = 1, n - 1
            best = i
            sbest = complex_score(values(i), which, sigma, shifted)
            do j = i + 1, n
                scur = complex_score(values(j), which, sigma, shifted)
                if (scur < sbest) then
                    best = j
                    sbest = scur
                end if
            end do
            if (best /= i) then
                tmp = values(i)
                values(i) = values(best)
                values(best) = tmp
                if (present(vectors)) then
                    vtmp = vectors(:, i)
                    vectors(:, i) = vectors(:, best)
                    vectors(:, best) = vtmp
                end if
            end if
        end do
    end subroutine sort_complex_pairs

    real(dp) function complex_score(x, which, sigma, shifted) result(s)
        complex(dp), intent(in) :: x
        character(len=*), intent(in) :: which
        complex(dp), intent(in), optional :: sigma
        logical, intent(in) :: shifted
        complex(dp) :: y
        y = x
        if (shifted .and. present(sigma)) then
            if (abs(x - sigma) > tiny(1.0_dp)) then
                y = 1.0_dp / (x - sigma)
            else
                y = cmplx(huge(1.0_dp), 0.0_dp, dp)
            end if
        end if
        select case (trim(which))
        case ('LM')
            s = -abs(y)
        case ('SM')
            s = abs(y)
        case ('LR')
            s = -real(y, dp)
        case ('SR')
            s = real(y, dp)
        case ('LI')
            s = -abs(aimag(y))
        case ('SI')
            s = abs(aimag(y))
        case default
            s = -abs(y)
        end select
    end function complex_score

end module rspectra_sort
