module sparsem_hb
    use sparsem_kinds, only: dp
    use sparsem_types, only: hb_matrix, csr_matrix, csc_matrix, ssc_matrix
    use sparsem_core, only: csr_from_csc, csr_from_ssc
    implicit none
    private

    public :: read_matrix_hb, hb_model_matrix, hb_model_response
    public :: hb_model_guess, hb_model_xexact

contains

    function read_matrix_hb(filename, info) result(data)
        character(len=*), intent(in) :: filename !! Path to a Harwell-Boeing file containing a real assembled sparse matrix.
        integer, intent(out), optional :: info !! Zero on success; negative values report I/O, format, or unsupported-storage errors.
        type(hb_matrix) :: data
        character(len=256) :: header(4)
        character(len=256) :: rhs_header
        character(len=32) :: ptrfmt
        character(len=32) :: indfmt
        character(len=32) :: valfmt
        character(len=32) :: rhsfmt
        character(len=3) :: matrix_type
        character(len=3) :: rhs_type
        real(dp), allocatable :: raw(:)
        integer :: unit
        integer :: ios
        integer :: stat
        integer :: ptr_lines
        integer :: ind_lines
        integer :: val_lines
        integer :: rhs_lines
        integer :: nz
        integer :: nrhs
        integer :: i

        if (present(info)) info = 0
        open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            if (present(info)) info = -1
            return
        end if

        do i = 1, 4
            read(unit, '(a)', iostat=ios) header(i)
            if (ios /= 0) then
                close(unit)
                if (present(info)) info = -2
                return
            end if
        end do

        data%title = trim(header(1)(1:72))
        data%key = trim(header(1)(73:80))
        call parse_integer_field(header(2), 15, 28, ptr_lines, stat)
        if (stat /= 0) goto 900
        call parse_integer_field(header(2), 29, 42, ind_lines, stat)
        if (stat /= 0) goto 900
        call parse_integer_field(header(2), 43, 56, val_lines, stat)
        if (stat /= 0) goto 900
        call parse_integer_field(header(2), 57, 70, rhs_lines, stat)
        if (stat /= 0) rhs_lines = 0

        matrix_type = upper_ascii(header(3)(1:3))
        if (matrix_type(1:1) /= 'R' .or. matrix_type(3:3) /= 'A') then
            close(unit)
            if (present(info)) info = -3
            return
        end if
        select case (matrix_type(2:2))
        case ('S')
            data%symmetric = .true.
        case ('R', 'U')
            data%symmetric = .false.
        case default
            close(unit)
            if (present(info)) info = -3
            return
        end select

        call parse_integer_field(header(3), 15, 28, data%nrow, stat)
        if (stat /= 0) goto 900
        call parse_integer_field(header(3), 29, 42, data%ncol, stat)
        if (stat /= 0) goto 900
        call parse_integer_field(header(3), 43, 56, nz, stat)
        if (stat /= 0) goto 900
        if (data%nrow < 0 .or. data%ncol < 0 .or. nz < 0) goto 900
        if (data%symmetric .and. data%nrow /= data%ncol) goto 900

        ptrfmt = trim(upper_ascii(header(4)(1:16)))
        indfmt = trim(upper_ascii(header(4)(17:32)))
        valfmt = trim(upper_ascii(header(4)(33:52)))
        rhsfmt = trim(upper_ascii(header(4)(53:72)))

        nrhs = 0
        if (rhs_lines > 0) then
            read(unit, '(a)', iostat=ios) rhs_header
            if (ios /= 0) then
                close(unit)
                if (present(info)) info = -2
                return
            end if
            rhs_type = upper_ascii(rhs_header(1:3))
            data%rhs_mode = rhs_type(1:1)
            if (data%rhs_mode /= 'F') then
                close(unit)
                if (present(info)) info = -6
                return
            end if
            if (rhs_type(2:2) /= ' ' .and. rhs_type(2:2) /= 'G') then
                close(unit)
                if (present(info)) info = -6
                return
            end if
            if (rhs_type(3:3) /= ' ' .and. rhs_type(3:3) /= 'X') then
                close(unit)
                if (present(info)) info = -6
                return
            end if
            call parse_integer_field(rhs_header, 15, 28, nrhs, stat)
            if (stat /= 0 .or. nrhs < 1) then
                close(unit)
                if (present(info)) info = -6
                return
            end if
        end if

        allocate(data%col_ptr(data%ncol + 1))
        allocate(data%row_ind(nz))
        allocate(data%values(nz))
        call read_hb_integers(unit, data%ncol + 1, ptrfmt, data%col_ptr, stat)
        if (stat /= 0) goto 910
        call read_hb_integers(unit, nz, indfmt, data%row_ind, stat)
        if (stat /= 0) goto 910
        call read_hb_reals(unit, nz, valfmt, data%values, stat)
        if (stat /= 0) goto 910

        if (.not. valid_hb_structure(data, nz)) then
            close(unit)
            if (present(info)) info = -7
            return
        end if

        if (rhs_lines > 0) then
            allocate(raw(data%nrow * nrhs))
            call read_hb_reals(unit, size(raw), rhsfmt, raw, stat)
            if (stat /= 0) goto 910
            allocate(data%rhs(data%nrow, nrhs))
            data%rhs = reshape(raw, [data%nrow, nrhs])
            data%has_rhs = .true.

            if (rhs_type(2:2) == 'G') then
                call read_hb_reals(unit, size(raw), rhsfmt, raw, stat)
                if (stat /= 0) goto 910
                allocate(data%guess(data%nrow, nrhs))
                data%guess = reshape(raw, [data%nrow, nrhs])
                data%has_guess = .true.
            end if
            if (rhs_type(3:3) == 'X') then
                call read_hb_reals(unit, size(raw), rhsfmt, raw, stat)
                if (stat /= 0) goto 910
                allocate(data%xexact(data%nrow, nrhs))
                data%xexact = reshape(raw, [data%nrow, nrhs])
                data%has_xexact = .true.
            end if
        end if

        close(unit)
        if (present(info)) info = 0
        return

900     continue
        close(unit)
        if (present(info)) info = -2
        return

910     continue
        close(unit)
        if (present(info)) info = -5
    end function read_matrix_hb

    function hb_model_matrix(data, info) result(x)
        type(hb_matrix), intent(in) :: data !! Harwell-Boeing object whose sparse coefficient matrix is requested.
        integer, intent(out), optional :: info !! Zero on success; -1 when the stored sparse structure is inconsistent.
        type(csr_matrix) :: x
        type(csc_matrix) :: csc
        type(ssc_matrix) :: ssc
        integer :: nz

        if (present(info)) info = 0
        nz = size(data%values)
        if (.not. valid_hb_structure(data, nz)) then
            if (present(info)) info = -1
            x%nrow = 0
            x%ncol = 0
            allocate(x%values(0))
            allocate(x%col_ind(0))
            allocate(x%row_ptr(1))
            x%row_ptr = 1
            return
        end if
        if (data%symmetric) then
            ssc%n = data%nrow
            ssc%values = data%values
            ssc%row_ind = data%row_ind
            ssc%col_ptr = data%col_ptr
            x = csr_from_ssc(ssc)
        else
            csc%nrow = data%nrow
            csc%ncol = data%ncol
            csc%values = data%values
            csc%row_ind = data%row_ind
            csc%col_ptr = data%col_ptr
            x = csr_from_csc(csc)
        end if
    end function hb_model_matrix

    function hb_model_response(data, info) result(y)
        type(hb_matrix), intent(in) :: data !! Harwell-Boeing object whose full-storage right-hand side is requested.
        integer, intent(out), optional :: info !! Zero on success; -1 when no supported right-hand side is present.
        real(dp), allocatable :: y(:,:)

        if (data%has_rhs .and. allocated(data%rhs) .and. data%rhs_mode == 'F') then
            y = data%rhs
            if (present(info)) info = 0
        else
            allocate(y(0, 0))
            if (present(info)) info = -1
        end if
    end function hb_model_response

    function hb_model_guess(data, info) result(x)
        type(hb_matrix), intent(in) :: data !! Harwell-Boeing object whose optional initial-solution guess is requested.
        integer, intent(out), optional :: info !! Zero on success; -1 when no initial guess is stored.
        real(dp), allocatable :: x(:,:)

        if (data%has_guess .and. allocated(data%guess)) then
            x = data%guess
            if (present(info)) info = 0
        else
            allocate(x(0, 0))
            if (present(info)) info = -1
        end if
    end function hb_model_guess

    function hb_model_xexact(data, info) result(x)
        type(hb_matrix), intent(in) :: data !! Harwell-Boeing object whose optional exact solution is requested.
        integer, intent(out), optional :: info !! Zero on success; -1 when no exact solution is stored.
        real(dp), allocatable :: x(:,:)

        if (data%has_xexact .and. allocated(data%xexact)) then
            x = data%xexact
            if (present(info)) info = 0
        else
            allocate(x(0, 0))
            if (present(info)) info = -1
        end if
    end function hb_model_xexact

    subroutine parse_integer_field(line, first, last, value, info)
        character(len=*), intent(in) :: line !! Fixed-width header record containing the integer field.
        integer, intent(in) :: first !! One-based first column of the field.
        integer, intent(in) :: last !! One-based last column of the field.
        integer, intent(out) :: value !! Parsed integer value.
        integer, intent(out) :: info !! Zero on success and one when the field is blank or invalid.
        character(len=64) :: field
        integer :: ios

        field = ' '
        if (first < 1 .or. last < first .or. first > len(line)) then
            value = 0
            info = 1
            return
        end if
        field = line(first:min(last, len(line)))
        read(field, *, iostat=ios) value
        if (ios /= 0) then
            value = 0
            info = 1
        else
            info = 0
        end if
    end subroutine parse_integer_field

    subroutine parse_hb_format(fmt, nper, width, code, info)
        character(len=*), intent(in) :: fmt !! Harwell-Boeing Fortran format such as `(16I5)` or `(10F7.1)`.
        integer, intent(out) :: nper !! Number of values represented on each record.
        integer, intent(out) :: width !! Fixed field width in characters.
        character(len=1), intent(out) :: code !! Conversion descriptor: I, D, E, F, or G.
        integer, intent(out) :: info !! Zero on success and one for an unsupported or malformed descriptor.
        character(len=:), allocatable :: f
        integer :: pos
        integer :: first_digit
        integer :: last_digit
        integer :: i
        integer :: ios

        f = trim(upper_ascii(fmt))
        pos = 0
        code = ' '
        do i = 1, len(f)
            if (index('IDEFG', f(i:i)) > 0) then
                pos = i
                code = f(i:i)
                exit
            end if
        end do
        if (pos == 0) then
            info = 1
            nper = 0
            width = 0
            return
        end if

        last_digit = pos - 1
        do while (last_digit >= 1 .and. (f(last_digit:last_digit) < '0' .or. f(last_digit:last_digit) > '9'))
            last_digit = last_digit - 1
        end do
        first_digit = last_digit
        do while (first_digit > 1 .and. f(first_digit - 1:first_digit - 1) >= '0' .and. &
            f(first_digit - 1:first_digit - 1) <= '9')
            first_digit = first_digit - 1
        end do
        if (last_digit < first_digit) then
            info = 1
            nper = 0
            width = 0
            return
        end if
        read(f(first_digit:last_digit), *, iostat=ios) nper
        if (ios /= 0 .or. nper < 1) then
            info = 1
            width = 0
            return
        end if

        first_digit = pos + 1
        do while (first_digit <= len(f) .and. (f(first_digit:first_digit) < '0' .or. f(first_digit:first_digit) > '9'))
            first_digit = first_digit + 1
        end do
        last_digit = first_digit
        do while (last_digit <= len(f) .and. f(last_digit:last_digit) >= '0' .and. f(last_digit:last_digit) <= '9')
            last_digit = last_digit + 1
        end do
        last_digit = last_digit - 1
        if (last_digit < first_digit) then
            info = 1
            width = 0
            return
        end if
        read(f(first_digit:last_digit), *, iostat=ios) width
        if (ios /= 0 .or. width < 1) then
            info = 1
            return
        end if
        info = 0
    end subroutine parse_hb_format

    subroutine read_hb_integers(unit, nvals, fmt, values, info)
        integer, intent(in) :: unit !! Open formatted unit positioned at the first integer data record.
        integer, intent(in) :: nvals !! Number of integer values to consume from consecutive fixed-width records.
        character(len=*), intent(in) :: fmt !! Harwell-Boeing integer field format.
        integer, intent(out) :: values(:) !! Destination array of length at least nvals.
        integer, intent(out) :: info !! Zero on success and one for malformed format, I/O failure, or invalid numeric text.
        character(len=512) :: line
        character(len=1) :: code
        character(len=128) :: field
        integer :: nper
        integer :: width
        integer :: stat
        integer :: ios
        integer :: count
        integer :: j
        integer :: first
        integer :: last

        info = 1
        if (nvals < 0 .or. size(values) < nvals) return
        if (nvals == 0) then
            info = 0
            return
        end if
        call parse_hb_format(fmt, nper, width, code, stat)
        if (stat /= 0 .or. code /= 'I') return
        count = 0
        do while (count < nvals)
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) return
            do j = 1, nper
                if (count >= nvals) exit
                first = (j - 1) * width + 1
                last = j * width
                if (first > len(line)) return
                field = ' '
                field = line(first:min(last, len(line)))
                read(field, *, iostat=ios) values(count + 1)
                if (ios /= 0) return
                count = count + 1
            end do
        end do
        info = 0
    end subroutine read_hb_integers

    subroutine read_hb_reals(unit, nvals, fmt, values, info)
        integer, intent(in) :: unit !! Open formatted unit positioned at the first real-valued data record.
        integer, intent(in) :: nvals !! Number of real values to consume from consecutive fixed-width records.
        character(len=*), intent(in) :: fmt !! Harwell-Boeing real field format.
        real(dp), intent(out) :: values(:) !! Destination array of length at least nvals.
        integer, intent(out) :: info !! Zero on success and one for malformed format, I/O failure, or invalid numeric text.
        character(len=512) :: line
        character(len=1) :: code
        character(len=128) :: field
        integer :: nper
        integer :: width
        integer :: stat
        integer :: ios
        integer :: count
        integer :: j
        integer :: first
        integer :: last
        integer :: k

        info = 1
        if (nvals < 0 .or. size(values) < nvals) return
        if (nvals == 0) then
            info = 0
            return
        end if
        call parse_hb_format(fmt, nper, width, code, stat)
        if (stat /= 0 .or. index('DEFG', code) == 0) return
        count = 0
        do while (count < nvals)
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) return
            do j = 1, nper
                if (count >= nvals) exit
                first = (j - 1) * width + 1
                last = j * width
                if (first > len(line)) return
                field = ' '
                field = line(first:min(last, len(line)))
                do k = 1, len(field)
                    if (field(k:k) == 'D' .or. field(k:k) == 'd') field(k:k) = 'E'
                end do
                read(field, *, iostat=ios) values(count + 1)
                if (ios /= 0) return
                count = count + 1
            end do
        end do
        info = 0
    end subroutine read_hb_reals

    pure logical function valid_hb_structure(data, nz) result(ok)
        type(hb_matrix), intent(in) :: data !! Parsed Harwell-Boeing sparse matrix metadata and compressed-column arrays.
        integer, intent(in) :: nz !! Expected number of stored sparse values and row indices.
        integer :: j

        ok = .false.
        if (.not. allocated(data%values) .or. .not. allocated(data%row_ind) .or. .not. allocated(data%col_ptr)) return
        if (size(data%values) /= nz .or. size(data%row_ind) /= nz) return
        if (size(data%col_ptr) /= data%ncol + 1) return
        if (data%nrow < 0 .or. data%ncol < 0) return
        if (data%col_ptr(1) /= 1 .or. data%col_ptr(data%ncol + 1) /= nz + 1) return
        do j = 1, data%ncol
            if (data%col_ptr(j + 1) < data%col_ptr(j)) return
        end do
        if (nz > 0) then
            if (any(data%row_ind < 1) .or. any(data%row_ind > data%nrow)) return
        end if
        if (data%symmetric .and. data%nrow /= data%ncol) return
        ok = .true.
    end function valid_hb_structure

    pure function upper_ascii(text) result(out)
        character(len=*), intent(in) :: text !! ASCII text to convert to uppercase without locale-dependent processing.
        character(len=len(text)) :: out
        integer :: i
        integer :: c

        out = text
        do i = 1, len(text)
            c = iachar(text(i:i))
            if (c >= iachar('a') .and. c <= iachar('z')) out(i:i) = achar(c - iachar('a') + iachar('A'))
        end do
    end function upper_ascii

end module sparsem_hb
