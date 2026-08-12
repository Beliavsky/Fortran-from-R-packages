! SPDX-License-Identifier: GPL-2.0-or-later
module linprog_mps
    use linprog_types
    use linprog_solver, only: solveLP
    implicit none
    private

    interface writeMps
        module procedure writeMps_arrays
        module procedure writeMps_model
    end interface writeMps

    public :: readMps, writeMps

contains

    subroutine readMps(filename, model, solve, maximum)
        character(len=*), intent(in) :: filename
        type(mps_model), intent(out) :: model
        logical, intent(in), optional :: solve, maximum
        character(len=512) :: line
        character(len=64) :: tok(8), objname, pname
        character(len=64), allocatable :: vars(:), rows(:)
        character(len=1), allocatable :: rtype(:)
        integer :: ntok, unit, ios, section, nbase, nup, i, j, iv, ir, ibound
        logical :: do_solve, do_max
        type(linprog_control) :: ctl

        do_solve = .false.
        do_max = .false.
        if (present(solve)) do_solve = solve
        if (present(maximum)) do_max = maximum
        objname = ''
        pname = ''
        allocate(vars(0), rows(0), rtype(0))
        nbase = 0
        nup = 0
        section = 0
        open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) error stop "readMps: cannot open file"
        do
            read(unit,'(a)',iostat=ios) line
            if (ios /= 0) exit
            call split_tokens(line, tok, ntok)
            if (ntok == 0) cycle
            if (tok(1)(1:1) == '*') cycle
            if (trim(tok(1)) == 'NAME') then
                section = 1
                if (ntok >= 2) pname = tok(2)
                cycle
            else if (ntok == 1) then
                select case (trim(tok(1)))
                case ('ROWS')
                    section = 2
                    cycle
                case ('COLUMNS')
                    section = 3
                    cycle
                case ('RHS')
                    section = 4
                    cycle
                case ('BOUNDS')
                    section = 5
                    cycle
                case ('ENDATA')
                    exit
                end select
            end if
            select case (section)
            case (2)
                if (ntok < 2) cycle
                select case (tok(1)(1:1))
                case ('N')
                    if (len_trim(objname) == 0) objname = tok(2)
                case ('L','G')
                    call append_row(rows, rtype, tok(2), tok(1)(1:1))
                    nbase = nbase + 1
                case ('E')
                    close(unit)
                    error stop "readMps: equality constraints are not implemented by original linprog parser"
                case default
                    close(unit)
                    error stop "readMps: invalid ROWS type"
                end select
            case (3)
                call append_unique(vars, tok(1))
            case (5)
                if (trim(tok(1)) == 'UP') then
                    nup = nup + 1
                else if (trim(tok(1)) == 'LO' .or. trim(tok(1)) == 'FX' .or. &
                    trim(tok(1)) == 'FR') then
                    close(unit)
                    error stop "readMps: LO, FX, and FR bounds are not implemented by original linprog parser"
                else
                    close(unit)
                    error stop "readMps: unsupported BOUNDS record"
                end if
            end select
        end do
        close(unit)

        if (len_trim(pname) == 0) pname = 'LP'
        allocate(character(len=max(1,len_trim(pname))) :: model%name)
        model%name = trim(pname)
        allocate(model%cvec(size(vars)), model%bvec(nbase+nup))
        allocate(model%amat(nbase+nup,size(vars)))
        allocate(character(len=64) :: model%var_names(size(vars)))
        allocate(character(len=64) :: model%con_names(nbase+nup))
        model%cvec = 0.0_dp
        model%bvec = 0.0_dp
        model%amat = 0.0_dp
        model%var_names = vars
        model%con_names = ''
        if (nbase > 0) model%con_names(1:nbase) = rows

        section = 0
        ir = 0
        ibound = nbase
        open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) error stop "readMps: cannot reopen file"
        do
            read(unit,'(a)',iostat=ios) line
            if (ios /= 0) exit
            call split_tokens(line, tok, ntok)
            if (ntok == 0) cycle
            if (tok(1)(1:1) == '*') cycle
            if (trim(tok(1)) == 'NAME') then
                section = 1
                cycle
            else if (ntok == 1) then
                select case (trim(tok(1)))
                case ('ROWS')
                    section = 2
                    cycle
                case ('COLUMNS')
                    section = 3
                    cycle
                case ('RHS')
                    section = 4
                    cycle
                case ('BOUNDS')
                    section = 5
                    cycle
                case ('ENDATA')
                    exit
                end select
            end if
            select case (section)
            case (3)
                iv = find_name(vars, tok(1))
                if (iv == 0) error stop "readMps: internal variable lookup failure"
                j = 2
                do while (j+1 <= ntok)
                    if (trim(tok(j)) == trim(objname)) then
                        read(tok(j+1),*,iostat=ios) model%cvec(iv)
                    else
                        ir = find_name(rows, tok(j))
                        if (ir == 0) error stop "readMps: undefined constraint in COLUMNS"
                        read(tok(j+1),*,iostat=ios) model%amat(ir,iv)
                    end if
                    if (ios /= 0) error stop "readMps: invalid numeric value in COLUMNS"
                    j = j + 2
                end do
            case (4)
                j = 2
                do while (j+1 <= ntok)
                    ir = find_name(rows, tok(j))
                    if (ir == 0) error stop "readMps: undefined constraint in RHS"
                    read(tok(j+1),*,iostat=ios) model%bvec(ir)
                    if (ios /= 0) error stop "readMps: invalid numeric value in RHS"
                    j = j + 2
                end do
            case (5)
                if (trim(tok(1)) == 'UP') then
                    if (ntok < 4) error stop "readMps: malformed UP bound"
                    iv = find_name(vars, tok(3))
                    if (iv == 0) error stop "readMps: undefined variable in BOUNDS"
                    ibound = ibound + 1
                    model%amat(ibound,iv) = 1.0_dp
                    read(tok(4),*,iostat=ios) model%bvec(ibound)
                    if (ios /= 0) error stop "readMps: invalid UP bound"
                    model%con_names(ibound) = 'UP'//trim(tok(3))
                end if
            end select
        end do
        close(unit)

        do i = 1, nbase
            if (rtype(i) == 'G') then
                model%bvec(i) = -model%bvec(i)
                model%amat(i,:) = -model%amat(i,:)
            end if
        end do

        if (do_solve) then
            ctl = linprog_control()
            ctl%maximum = do_max
            call solveLP(model%cvec, model%bvec, model%amat, model%result, ctl)
            model%has_result = .true.
        end if
    end subroutine readMps

    subroutine writeMps_model(filename, model)
        character(len=*), intent(in) :: filename
        type(mps_model), intent(in) :: model
        call writeMps_arrays(filename, model%cvec, model%bvec, model%amat, model%name, &
            model%var_names, model%con_names)
    end subroutine writeMps_model

    subroutine writeMps_arrays(filename, cvec, bvec, amat, name, var_names, con_names)
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: cvec(:), bvec(:), amat(:,:)
        character(len=*), intent(in), optional :: name
        character(len=*), intent(in), optional :: var_names(:), con_names(:)
        character(len=8), allocatable :: vname(:), rname(:)
        character(len=64) :: pname
        character(len=256) :: line
        integer :: nvar, ncon, unit, i, j

        nvar = size(cvec)
        ncon = size(bvec)
        if (size(amat,1) /= ncon .or. size(amat,2) /= nvar) then
            error stop "writeMps: inconsistent dimensions"
        end if
        pname = 'LP'
        if (present(name)) pname = trim(name)
        allocate(vname(nvar), rname(ncon))
        do i = 1, nvar
            if (present(var_names)) then
                vname(i) = clean_label(var_names(i), 'C', i, vname(1:max(0,i-1)))
            else
                vname(i) = clean_label('', 'C', i, vname(1:max(0,i-1)))
            end if
        end do
        do i = 1, ncon
            if (present(con_names)) then
                rname(i) = clean_label(con_names(i), 'R', i, rname(1:max(0,i-1)))
            else
                rname(i) = clean_label('', 'R', i, rname(1:max(0,i-1)))
            end if
        end do

        open(newunit=unit, file=filename, status='replace', action='write')
        write(unit,'(a)') 'NAME          '//trim(pname)
        write(unit,'(a)') 'ROWS'
        write(unit,'(a)') ' N  obj'
        do i = 1, ncon
            write(unit,'(a)') ' L  '//trim(rname(i))
        end do
        write(unit,'(a)') 'COLUMNS'
        do i = 1, nvar
            write(line,'(a,1x,a,1x,es24.16)') '    '//trim(vname(i)), 'obj', cvec(i)
            write(unit,'(a)') trim(line)
            do j = 1, ncon
                if (abs(amat(j,i)) > 0.0_dp) then
                    write(line,'(a,1x,a,1x,es24.16)') '    '//trim(vname(i)), trim(rname(j)), amat(j,i)
                    write(unit,'(a)') trim(line)
                end if
            end do
        end do
        write(unit,'(a)') 'RHS'
        do i = 1, ncon
            write(line,'(a,1x,a,1x,es24.16)') '    RHS', trim(rname(i)), bvec(i)
            write(unit,'(a)') trim(line)
        end do
        write(unit,'(a)') 'ENDATA'
        close(unit)
    end subroutine writeMps_arrays

    subroutine split_tokens(line, tok, ntok)
        character(len=*), intent(in) :: line
        character(len=64), intent(out) :: tok(:)
        integer, intent(out) :: ntok
        integer :: i, n, first, last
        ntok = 0
        tok = ''
        n = len_trim(line)
        i = 1
        do while (i <= n .and. ntok < size(tok))
            do while (i <= n .and. (line(i:i) == ' ' .or. line(i:i) == achar(9)))
                i = i + 1
            end do
            if (i > n) exit
            first = i
            do while (i <= n .and. line(i:i) /= ' ' .and. line(i:i) /= achar(9))
                i = i + 1
            end do
            last = i-1
            ntok = ntok+1
            tok(ntok) = line(first:min(last,first+63))
        end do
    end subroutine split_tokens

    subroutine append_unique(a, value)
        character(len=64), allocatable, intent(inout) :: a(:)
        character(len=*), intent(in) :: value
        character(len=64), allocatable :: tmp(:)
        integer :: n
        if (find_name(a, value) > 0) return
        n = size(a)
        allocate(tmp(n+1))
        if (n > 0) tmp(1:n) = a
        tmp(n+1) = trim(value)
        call move_alloc(tmp, a)
    end subroutine append_unique

    subroutine append_row(rows, rtype, name, typ)
        character(len=64), allocatable, intent(inout) :: rows(:)
        character(len=1), allocatable, intent(inout) :: rtype(:)
        character(len=*), intent(in) :: name, typ
        character(len=64), allocatable :: tr(:)
        character(len=1), allocatable :: tt(:)
        integer :: n
        n = size(rows)
        allocate(tr(n+1), tt(n+1))
        if (n > 0) then
            tr(1:n) = rows
            tt(1:n) = rtype
        end if
        tr(n+1) = trim(name)
        tt(n+1) = typ(1:1)
        call move_alloc(tr, rows)
        call move_alloc(tt, rtype)
    end subroutine append_row

    pure integer function find_name(a, value) result(idx)
        character(len=*), intent(in) :: a(:), value
        integer :: i
        idx = 0
        do i = 1, size(a)
            if (trim(a(i)) == trim(value)) then
                idx = i
                return
            end if
        end do
    end function find_name

    function clean_label(label, prefix, idx, previous) result(out)
        character(len=*), intent(in) :: label, prefix
        integer, intent(in) :: idx
        character(len=*), intent(in) :: previous(:)
        character(len=8) :: out, base, candidate
        character(len=16) :: num
        integer :: i, j, n
        logical :: duplicate

        base = ''
        if (len_trim(label) == 0) then
            write(num,'(i0)') idx
            base = trim(prefix)//'_'//trim(num)
        else
            j = 0
            do i = 1, len_trim(label)
                if (label(i:i) /= ' ') then
                    j = j+1
                    if (j <= 8) base(j:j) = label(i:i)
                end if
            end do
        end if
        candidate = base
        n = 2
        do
            duplicate = .false.
            do i = 1, size(previous)
                if (trim(previous(i)) == trim(candidate)) duplicate = .true.
            end do
            if (.not. duplicate) exit
            write(num,'(i0)') n
            candidate = ''
            j = max(1, 7-len_trim(num))
            candidate = base(1:min(j,len_trim(base)))//'_'//trim(num)
            n = n+1
        end do
        out = candidate
    end function clean_label

end module linprog_mps
