! SPDX-License-Identifier: Artistic-2.0
! Derived from the R matlab package; see COPYRIGHTS and upstream/.
module matlab_path
    use matlab_types, only : fileparts_result
    implicit none
    private

    public :: filesep
    public :: pathsep
    public :: fullfile
    public :: fileparts
    public :: strcmp

contains

    function filesep() result(sep)
        character(len=1) :: sep
        character(len=32) :: os
        integer :: status

        call get_environment_variable('OS', os, status=status)
        if (status == 0 .and. index(os, 'Windows') > 0) then
            sep = achar(92)
        else
            sep = '/'
        end if
    end function filesep

    function pathsep() result(sep)
        character(len=1) :: sep

        if (filesep() == achar(92)) then
            sep = ';'
        else
            sep = ':'
        end if
    end function pathsep

    function fullfile(a, b, c) result(path)
        character(len=*), intent(in) :: a
        character(len=*), intent(in), optional :: b
        character(len=*), intent(in), optional :: c
        character(len=:), allocatable :: path

        path = trim(a)
        if (present(b)) path = join_path(path, trim(b))
        if (present(c)) path = join_path(path, trim(c))
    end function fullfile

    function fileparts(pathname) result(parts)
        character(len=*), intent(in) :: pathname
        type(fileparts_result) :: parts
        character(len=:), allocatable :: path, fname
        integer :: last_sep, last_dot, n
        logical :: trailing

        path = trim(pathname)
        n = len(path)
        trailing = n > 0 .and. is_sep(path(n:n))
        last_sep = last_separator(path)

        if (trailing) then
            parts%pathstr = path(:n - 1)
            parts%name = ''
            parts%ext = ''
            parts%versn = ''
            return
        end if

        if (last_sep > 0) then
            parts%pathstr = path(:last_sep - 1)
            fname = path(last_sep + 1:)
        else
            parts%pathstr = ''
            fname = path
        end if

        if (fname == '.') then
            parts%name = ''
            parts%ext = '.'
        else if (fname == '..') then
            parts%name = '.'
            parts%ext = '.'
        else
            last_dot = scan(fname, '.', back=.true.)
            if (last_dot == 1) then
                parts%name = ''
                parts%ext = fname
            else if (last_dot > 1) then
                parts%name = fname(:last_dot - 1)
                parts%ext = fname(last_dot:)
            else
                parts%name = fname
                parts%ext = ''
            end if
        end if
        parts%versn = ''
    end function fileparts

    function strcmp(s, t) result(equal)
        character(len=*), intent(in) :: s
        character(len=*), intent(in) :: t
        logical :: equal

        equal = s == t
    end function strcmp

    pure function join_path(a, b) result(path)
        character(len=*), intent(in) :: a
        character(len=*), intent(in) :: b
        character(len=:), allocatable :: path
        character(len=1) :: sep

        sep = '/'
        if (len_trim(a) == 0) then
            path = b
        else if (len_trim(b) == 0) then
            path = a
        else if (is_sep(a(len_trim(a):len_trim(a)))) then
            path = trim(a) // trim(b)
        else
            path = trim(a) // sep // trim(b)
        end if
    end function join_path

    pure function last_separator(path) result(pos)
        character(len=*), intent(in) :: path
        integer :: pos, i

        pos = 0
        do i = len_trim(path), 1, -1
            if (is_sep(path(i:i))) then
                pos = i
                return
            end if
        end do
    end function last_separator

    pure function is_sep(ch) result(answer)
        character(len=1), intent(in) :: ch
        logical :: answer

        answer = ch == '/' .or. ch == achar(92)
    end function is_sep
end module matlab_path
