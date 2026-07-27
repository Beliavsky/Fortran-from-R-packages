! Part of the modern Fortran translation of longmemo 1.1-4.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original longmemo authors retain copyright; see ORIGINAL_PACKAGE.txt.
! SPDX-License-Identifier: GPL-2.0-or-later

module longmemo_io
    use longmemo_kinds, only : dp
    implicit none
    private

    public :: read_index_value_csv

contains

    subroutine read_index_value_csv(filename, values)
        character(len=*), intent(in) :: filename
        real(dp), allocatable, intent(out) :: values(:)
        character(len=1024) :: line
        integer :: unit, ios, n, index_value, i
        real(dp) :: value

        open(newunit=unit, file=filename, status="old", action="read", iostat=ios)
        if (ios /= 0) error stop "read_index_value_csv: cannot open input file"

        read(unit, '(a)', iostat=ios) line
        if (ios /= 0) error stop "read_index_value_csv: missing header"
        n = 0
        do
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) > 0) n = n + 1
        end do
        rewind(unit)
        read(unit, '(a)') line
        allocate(values(n))
        i = 0
        do
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) == 0) cycle
            read(line, *, iostat=ios) index_value, value
            if (ios /= 0) error stop "read_index_value_csv: malformed data row"
            i = i + 1
            values(i) = value
        end do
        close(unit)
    end subroutine read_index_value_csv

end module longmemo_io
