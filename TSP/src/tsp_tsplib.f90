! Modern Fortran translation of computational code from TSP 1.2.7.
! Original Copyright (C) Michael Hahsler and Kurt Hornik.
! SPDX-License-Identifier: GPL-3.0-only
! See LICENSE, COPYING, and UPSTREAM.md for provenance and licensing.

module tsp_tsplib
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use tsp_kinds, only : dp
    use tsp_core, only : replace_infinite
    implicit none
    private

    type, public :: tsplib_instance
        character(len=16) :: problem_type = "TSP"
        character(len=24) :: edge_weight_type = "EXPLICIT"
        integer :: dimension = 0
        real(dp), allocatable :: cost(:,:)
        real(dp), allocatable :: coords(:,:)
    end type tsplib_instance

    public :: read_tsplib, write_tsplib_tsp, write_tsplib_atsp, write_tsplib_etsp
    public :: tsplib_att_distance, tsplib_geo_distance

contains

    subroutine read_tsplib(filename, instance, precision, ierr)
        character(len=*), intent(in) :: filename
        type(tsplib_instance), intent(out) :: instance
        integer, intent(in), optional :: precision
        integer, intent(out), optional :: ierr
        character(len=1024) :: line
        character(len=128) :: key, value
        character(len=24) :: weight_format
        integer :: unit, ios, p, nvals, expected, ncoord, ndim
        real(dp), allocatable :: vals(:), rowvals(:)
        logical :: in_weights, in_coords

        if (present(ierr)) ierr = 0
        p = 0
        if (present(precision)) p = precision
        weight_format = ""
        in_weights = .false.
        in_coords = .false.
        nvals = 0
        ncoord = 0
        open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            if (present(ierr)) ierr = 1
            return
        end if

        do
            read(unit,'(A)',iostat=ios) line
            if (ios /= 0) exit
            if (.not. in_weights .and. .not. in_coords) then
                if (index(adjustl(to_upper(line)), 'EDGE_WEIGHT_SECTION') == 1) then
                    if (instance%dimension <= 0) then
                        if (present(ierr)) ierr = 2
                        close(unit)
                        return
                    end if
                    expected = expected_weight_count(instance%dimension, trim(weight_format))
                    if (expected <= 0) then
                        if (present(ierr)) ierr = 3
                        close(unit)
                        return
                    end if
                    allocate(vals(expected))
                    in_weights = .true.
                    cycle
                else if (index(adjustl(to_upper(line)), 'NODE_COORD_SECTION') == 1) then
                    if (instance%dimension <= 0) then
                        if (present(ierr)) ierr = 2
                        close(unit)
                        return
                    end if
                    if (trim(instance%edge_weight_type) == 'EUC_3D') then
                        ndim = 3
                    else
                        ndim = 2
                    end if
                    allocate(instance%coords(instance%dimension,ndim))
                    in_coords = .true.
                    cycle
                end if
                call parse_metadata(line, key, value)
                select case(trim(key))
                case('TYPE')
                    if (index(value,'ATSP') == 1 .or. index(value,'ATS') == 1) then
                        instance%problem_type = 'ATSP'
                    else if (index(value,'TSP') == 1) then
                        instance%problem_type = 'TSP'
                    end if
                case('DIMENSION')
                    read(value,*,iostat=ios) instance%dimension
                case('EDGE_WEIGHT_TYPE')
                    instance%edge_weight_type = trim(value)
                case('EDGE_WEIGHT_FORMAT')
                    weight_format = trim(value)
                end select
            else if (in_weights) then
                if (index(to_upper(line),'EOF') > 0) exit
                call append_real_tokens(line, vals, nvals)
                if (nvals >= size(vals)) exit
            else if (in_coords) then
                if (ncoord >= instance%dimension) exit
                allocate(rowvals(4))
                rowvals = 0.0_dp
                call parse_coord_line(line, rowvals, ios)
                if (ios == 0) then
                    ncoord = ncoord + 1
                    instance%coords(ncoord,:) = rowvals(2:1+size(instance%coords,2))
                end if
                deallocate(rowvals)
                if (ncoord >= instance%dimension) exit
            end if
        end do
        close(unit)

        if (in_weights) then
            if (nvals < size(vals)) then
                if (present(ierr)) ierr = 4
                return
            end if
            if (p /= 0) vals = vals / (10.0_dp**p)
            call explicit_to_matrix(vals, instance%dimension, trim(weight_format), &
                trim(instance%problem_type) == 'ATSP', instance%cost, ios)
            if (ios /= 0 .and. present(ierr)) ierr = 5
        else if (in_coords) then
            if (ncoord /= instance%dimension) then
                if (present(ierr)) ierr = 6
                return
            end if
            select case(trim(instance%edge_weight_type))
            case('ATT')
                call tsplib_att_distance(instance%coords, instance%cost)
            case('GEO')
                call tsplib_geo_distance(instance%coords, instance%cost)
            case('EUC_2D','EUC_3D')
                instance%problem_type = 'ETSP'
            case default
                if (present(ierr)) ierr = 7
            end select
        else
            if (present(ierr)) ierr = 8
        end if
    end subroutine read_tsplib

    subroutine parse_metadata(line, key, value)
        character(len=*), intent(in) :: line
        character(len=*), intent(out) :: key, value
        integer :: k
        key = ''
        value = ''
        k = index(line, ':')
        if (k <= 0) return
        key = trim(to_upper(adjustl(line(:k-1))))
        value = trim(to_upper(adjustl(line(k+1:))))
    end subroutine parse_metadata

    pure function to_upper(s) result(out)
        character(len=*), intent(in) :: s
        character(len=len(s)) :: out
        integer :: i, c
        out = s
        do i = 1, len(s)
            c = iachar(out(i:i))
            if (c >= iachar('a') .and. c <= iachar('z')) out(i:i) = achar(c - 32)
        end do
    end function to_upper

    integer function expected_weight_count(n, format) result(k)
        integer, intent(in) :: n
        character(len=*), intent(in) :: format
        select case(trim(format))
        case('FULL_MATRIX')
            k = n*n
        case('UPPER_ROW','LOWER_COL','UPPER_COL','LOWER_ROW')
            k = n*(n-1)/2
        case('UPPER_DIAG_ROW','LOWER_DIAG_COL','UPPER_DIAG_COL','LOWER_DIAG_ROW')
            k = n*(n+1)/2
        case default
            k = -1
        end select
    end function expected_weight_count

    subroutine append_real_tokens(line, vals, nvals)
        character(len=*), intent(in) :: line
        real(dp), intent(inout) :: vals(:)
        integer, intent(inout) :: nvals
        character(len=128) :: token
        integer :: i, j, n, ios
        real(dp) :: x

        n = len_trim(line)
        i = 1
        do while (i <= n .and. nvals < size(vals))
            do while (i <= n .and. (line(i:i) == ' ' .or. line(i:i) == achar(9)))
                i = i + 1
            end do
            if (i > n) exit
            j = i
            do while (j <= n .and. line(j:j) /= ' ' .and. line(j:j) /= achar(9))
                j = j + 1
            end do
            token = ''
            token = line(i:j-1)
            if (trim(to_upper(token)) == 'EOF') exit
            read(token,*,iostat=ios) x
            if (ios == 0) then
                nvals = nvals + 1
                vals(nvals) = x
            end if
            i = j + 1
        end do
    end subroutine append_real_tokens

    subroutine parse_coord_line(line, vals, ierr)
        character(len=*), intent(in) :: line
        real(dp), intent(out) :: vals(:)
        integer, intent(out) :: ierr
        character(len=128) :: token
        integer :: i, j, n, count, ios

        vals = 0.0_dp
        ierr = 0
        count = 0
        n = len_trim(line)
        i = 1
        do while (i <= n .and. count < size(vals))
            do while (i <= n .and. (line(i:i) == ' ' .or. line(i:i) == achar(9)))
                i = i + 1
            end do
            if (i > n) exit
            j = i
            do while (j <= n .and. line(j:j) /= ' ' .and. line(j:j) /= achar(9))
                j = j + 1
            end do
            token = ''
            token = line(i:j-1)
            read(token,*,iostat=ios) vals(count+1)
            if (ios /= 0) then
                ierr = 1
                return
            end if
            count = count + 1
            i = j + 1
        end do
        if (count < 3) ierr = 1
    end subroutine parse_coord_line

    subroutine explicit_to_matrix(vals, n, format, asymmetric, cost, ierr)
        real(dp), intent(in) :: vals(:)
        integer, intent(in) :: n
        character(len=*), intent(in) :: format
        logical, intent(in) :: asymmetric
        real(dp), allocatable, intent(out) :: cost(:,:)
        integer, intent(out) :: ierr
        integer :: i, j, k

        ierr = 0
        allocate(cost(n,n), source=0.0_dp)
        k = 0
        if (asymmetric) then
            if (trim(format) /= 'FULL_MATRIX') then
                ierr = 1
                return
            end if
            do i = 1, n
                do j = 1, n
                    k = k + 1
                    cost(i,j) = vals(k)
                end do
            end do
            return
        end if

        select case(trim(format))
        case('FULL_MATRIX')
            do i = 1, n
                do j = 1, n
                    k = k + 1
                    cost(i,j) = vals(k)
                end do
            end do
        case('UPPER_ROW')
            do i = 1, n-1
                do j = i+1, n
                    k = k + 1
                    cost(i,j) = vals(k); cost(j,i) = vals(k)
                end do
            end do
        case('LOWER_ROW')
            do i = 2, n
                do j = 1, i-1
                    k = k + 1
                    cost(i,j) = vals(k); cost(j,i) = vals(k)
                end do
            end do
        case('UPPER_COL')
            do j = 2, n
                do i = 1, j-1
                    k = k + 1
                    cost(i,j) = vals(k); cost(j,i) = vals(k)
                end do
            end do
        case('LOWER_COL')
            do j = 1, n-1
                do i = j+1, n
                    k = k + 1
                    cost(i,j) = vals(k); cost(j,i) = vals(k)
                end do
            end do
        case('UPPER_DIAG_ROW')
            do i = 1, n
                do j = i, n
                    k = k + 1
                    cost(i,j) = vals(k); cost(j,i) = vals(k)
                end do
            end do
        case('LOWER_DIAG_ROW')
            do i = 1, n
                do j = 1, i
                    k = k + 1
                    cost(i,j) = vals(k); cost(j,i) = vals(k)
                end do
            end do
        case('UPPER_DIAG_COL')
            do j = 1, n
                do i = 1, j
                    k = k + 1
                    cost(i,j) = vals(k); cost(j,i) = vals(k)
                end do
            end do
        case('LOWER_DIAG_COL')
            do j = 1, n
                do i = j, n
                    k = k + 1
                    cost(i,j) = vals(k); cost(j,i) = vals(k)
                end do
            end do
        case default
            ierr = 1
        end select
    end subroutine explicit_to_matrix

    subroutine tsplib_att_distance(coords, cost)
        real(dp), intent(in) :: coords(:,:)
        real(dp), allocatable, intent(out) :: cost(:,:)
        integer :: i, j, n
        real(dp) :: rij, tij, d2

        n = size(coords,1)
        allocate(cost(n,n), source=0.0_dp)
        do i = 1, n-1
            do j = i+1, n
                d2 = sum((coords(i,1:2) - coords(j,1:2))**2)
                rij = sqrt(d2 / 10.0_dp)
                tij = anint(rij)
                if (tij < rij) tij = tij + 1.0_dp
                cost(i,j) = tij
                cost(j,i) = tij
            end do
        end do
    end subroutine tsplib_att_distance

    subroutine tsplib_geo_distance(coords, cost)
        real(dp), intent(in) :: coords(:,:)
        real(dp), allocatable, intent(out) :: cost(:,:)
        real(dp), allocatable :: lat(:), lon(:)
        real(dp) :: q1, q2, q3, z
        integer :: i, j, n
        real(dp), parameter :: pi = acos(-1.0_dp)

        n = size(coords,1)
        allocate(lat(n), lon(n), cost(n,n), source=0.0_dp)
        do i = 1, n
            lat(i) = pi * (anint(coords(i,1)) + 5.0_dp * (coords(i,1) - anint(coords(i,1))) / 3.0_dp) / 180.0_dp
            lon(i) = pi * (anint(coords(i,2)) + 5.0_dp * (coords(i,2) - anint(coords(i,2))) / 3.0_dp) / 180.0_dp
        end do
        do i = 1, n-1
            do j = i+1, n
                q1 = cos(lon(i) - lon(j))
                q2 = cos(lat(i) - lat(j))
                q3 = cos(lat(i) + lat(j))
                z = 0.5_dp * ((1.0_dp + q1) * q2 - (1.0_dp - q1) * q3)
                z = max(-1.0_dp, min(1.0_dp, z))
                cost(i,j) = real(int(6378.388_dp * acos(z) + 1.0_dp), dp)
                cost(j,i) = cost(i,j)
            end do
        end do
    end subroutine tsplib_geo_distance

    subroutine write_tsplib_tsp(filename, cost, precision, positive_inf, negative_inf, ierr)
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: cost(:,:)
        integer, intent(in), optional :: precision
        real(dp), intent(in), optional :: positive_inf, negative_inf
        integer, intent(out), optional :: ierr
        real(dp), allocatable :: x(:,:)
        integer :: p, unit, ios, i, j
        integer(kind=8) :: w

        if (present(ierr)) ierr = 0
        p = 6; if (present(precision)) p = precision
        call replace_for_write(cost, x, positive_inf, negative_inf, ios)
        if (ios /= 0) then
            if (present(ierr)) ierr = 1
            return
        end if
        open(newunit=unit,file=filename,status='replace',action='write',iostat=ios)
        if (ios /= 0) then
            if (present(ierr)) ierr = 2
            return
        end if
        write(unit,'(A)') 'NAME: TSP'
        write(unit,'(A)') 'COMMENT: Generated by tsp-fortran'
        write(unit,'(A)') 'TYPE: TSP'
        write(unit,'(A,I0)') 'DIMENSION: ', size(cost,1)
        write(unit,'(A)') 'EDGE_WEIGHT_TYPE: EXPLICIT'
        write(unit,'(A)') 'EDGE_WEIGHT_FORMAT: UPPER_ROW'
        write(unit,'(A)') 'EDGE_WEIGHT_SECTION'
        do i = 1, size(cost,1)-1
            do j = i+1, size(cost,1)
                w = int(x(i,j) * 10.0_dp**p, kind=8)
                write(unit,'(I0)') w
            end do
        end do
        write(unit,'(A)') 'EOF'
        close(unit)
    end subroutine write_tsplib_tsp

    subroutine write_tsplib_atsp(filename, cost, precision, positive_inf, negative_inf, ierr)
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: cost(:,:)
        integer, intent(in), optional :: precision
        real(dp), intent(in), optional :: positive_inf, negative_inf
        integer, intent(out), optional :: ierr
        real(dp), allocatable :: x(:,:)
        integer :: p, unit, ios, i, j
        integer(kind=8) :: w

        if (present(ierr)) ierr = 0
        p = 6; if (present(precision)) p = precision
        call replace_for_write(cost, x, positive_inf, negative_inf, ios)
        if (ios /= 0) then
            if (present(ierr)) ierr = 1
            return
        end if
        open(newunit=unit,file=filename,status='replace',action='write',iostat=ios)
        if (ios /= 0) then
            if (present(ierr)) ierr = 2
            return
        end if
        write(unit,'(A)') 'NAME: ATSP'
        write(unit,'(A)') 'COMMENT: Generated by tsp-fortran'
        write(unit,'(A)') 'TYPE: ATSP'
        write(unit,'(A,I0)') 'DIMENSION: ', size(cost,1)
        write(unit,'(A)') 'EDGE_WEIGHT_TYPE: EXPLICIT'
        write(unit,'(A)') 'EDGE_WEIGHT_FORMAT: FULL_MATRIX'
        write(unit,'(A)') 'EDGE_WEIGHT_SECTION'
        do i = 1, size(cost,1)
            do j = 1, size(cost,2)
                w = int(x(i,j) * 10.0_dp**p, kind=8)
                write(unit,'(I0)') w
            end do
        end do
        write(unit,'(A)') 'EOF'
        close(unit)
    end subroutine write_tsplib_atsp

    subroutine write_tsplib_etsp(filename, coords, precision, ierr)
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: coords(:,:)
        integer, intent(in), optional :: precision
        integer, intent(out), optional :: ierr
        integer :: p, unit, ios, i
        character(len=64) :: fmt

        if (present(ierr)) ierr = 0
        p = 6; if (present(precision)) p = precision
        if (size(coords,2) /= 2 .and. size(coords,2) /= 3) then
            if (present(ierr)) ierr = 1
            return
        end if
        if (any(.not. ieee_is_finite(coords))) then
            if (present(ierr)) ierr = 2
            return
        end if
        open(newunit=unit,file=filename,status='replace',action='write',iostat=ios)
        if (ios /= 0) then
            if (present(ierr)) ierr = 3
            return
        end if
        write(unit,'(A)') 'NAME: ETSP'
        write(unit,'(A)') 'COMMENT: Generated by tsp-fortran'
        write(unit,'(A)') 'TYPE: TSP'
        write(unit,'(A,I0)') 'DIMENSION: ', size(coords,1)
        if (size(coords,2) == 2) then
            write(unit,'(A)') 'EDGE_WEIGHT_TYPE: EUC_2D'
        else
            write(unit,'(A)') 'EDGE_WEIGHT_TYPE: EUC_3D'
        end if
        write(unit,'(A)') 'NODE_COORD_SECTION'
        write(fmt,'("(I0,1X,",I0,"(ES",I0,".",I0,",1X))")') size(coords,2), p+8, p
        do i = 1, size(coords,1)
            if (size(coords,2) == 2) then
                write(unit,*) i, coords(i,1), coords(i,2)
            else
                write(unit,*) i, coords(i,1), coords(i,2), coords(i,3)
            end if
        end do
        write(unit,'(A)') 'EOF'
        close(unit)
    end subroutine write_tsplib_etsp

    subroutine replace_for_write(cost, x, positive_inf, negative_inf, ierr)
        real(dp), intent(in) :: cost(:,:)
        real(dp), allocatable, intent(out) :: x(:,:)
        real(dp), intent(in), optional :: positive_inf, negative_inf
        integer, intent(out) :: ierr
        if (present(positive_inf) .and. present(negative_inf)) then
            call replace_infinite(cost, x, positive_inf, negative_inf, ierr)
        else if (present(positive_inf)) then
            call replace_infinite(cost, x, positive_value=positive_inf, ierr=ierr)
        else if (present(negative_inf)) then
            call replace_infinite(cost, x, negative_value=negative_inf, ierr=ierr)
        else
            call replace_infinite(cost, x, ierr=ierr)
        end if
    end subroutine replace_for_write

end module tsp_tsplib
