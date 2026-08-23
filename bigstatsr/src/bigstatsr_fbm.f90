! SPDX-License-Identifier: GPL-3.0-only
module bigstatsr_fbm
    use iso_fortran_env, only: int8, int32, int64
    use bigstatsr_kinds, only: dp
    implicit none
    private

    type, public :: fbm_real
        character(len=:), allocatable :: filename
        integer :: nrow = 0
        integer :: ncol = 0
    contains
        procedure :: get => fbm_get
        procedure :: set => fbm_set
        procedure :: read_col => fbm_read_col
        procedure :: write_col => fbm_write_col
        procedure :: read_cols => fbm_read_cols
        procedure :: write_cols => fbm_write_cols
        procedure :: to_array => fbm_to_array
    end type fbm_real

    type, public :: fbm_code256
        character(len=:), allocatable :: filename
        integer :: nrow = 0
        integer :: ncol = 0
    contains
        procedure :: read_col => code_read_col
        procedure :: write_col => code_write_col
    end type fbm_code256

    public :: create_fbm, attach_fbm, create_code256, attach_code256
    public :: fbm_from_array, fbm_copy, fbm_increment, fbm_transpose

contains

    function attach_fbm(filename, nrow, ncol) result(x)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: nrow, ncol
        type(fbm_real) :: x
        x%filename = trim(filename)
        x%nrow = nrow
        x%ncol = ncol
    end function attach_fbm

    function create_fbm(filename, nrow, ncol, init) result(x)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: nrow, ncol
        real(dp), intent(in), optional :: init
        type(fbm_real) :: x
        real(dp), allocatable :: buf(:)
        real(dp) :: val
        integer :: unit, ios, left, m
        val = 0.0_dp
        if (present(init)) val = init
        x = attach_fbm(filename, nrow, ncol)
        open(newunit=unit, file=x%filename, access='stream', form='unformatted', &
             status='replace', action='write', iostat=ios)
        if (ios /= 0) error stop 'create_fbm: cannot create backing file'
        allocate(buf(min(max(1,nrow*ncol), 1048576)))
        buf = val
        left = nrow*ncol
        do while (left > 0)
            m = min(left, size(buf))
            write(unit) buf(1:m)
            left = left - m
        end do
        close(unit)
    end function create_fbm

    subroutine fbm_from_array(x, a)
        type(fbm_real), intent(in) :: x
        real(dp), intent(in) :: a(:,:)
        integer :: unit, ios
        if (size(a,1) /= x%nrow .or. size(a,2) /= x%ncol) &
            error stop 'fbm_from_array: dimension mismatch'
        open(newunit=unit, file=x%filename, access='stream', form='unformatted', &
             status='old', action='write', iostat=ios)
        if (ios /= 0) error stop 'fbm_from_array: cannot open backing file'
        write(unit, pos=1) a
        close(unit)
    end subroutine fbm_from_array

    function fbm_get(self, i, j) result(value)
        class(fbm_real), intent(in) :: self
        integer, intent(in) :: i, j
        real(dp) :: value, dummy
        integer :: unit, ios, reclen
        integer(int64) :: pos
        if (i < 1 .or. i > self%nrow .or. j < 1 .or. j > self%ncol) &
            error stop 'fbm_get: index out of bounds'
        inquire(iolength=reclen) dummy
        pos = int((j-1)*self%nrow + i - 1, int64) * int(reclen,int64) + 1_int64
        open(newunit=unit, file=self%filename, access='stream', form='unformatted', &
             status='old', action='read', iostat=ios)
        if (ios /= 0) error stop 'fbm_get: cannot open backing file'
        read(unit, pos=pos, iostat=ios) value
        close(unit)
        if (ios /= 0) error stop 'fbm_get: read failed'
    end function fbm_get

    subroutine fbm_set(self, i, j, value)
        class(fbm_real), intent(in) :: self
        integer, intent(in) :: i, j
        real(dp), intent(in) :: value
        real(dp) :: dummy
        integer :: unit, ios, reclen
        integer(int64) :: pos
        if (i < 1 .or. i > self%nrow .or. j < 1 .or. j > self%ncol) &
            error stop 'fbm_set: index out of bounds'
        inquire(iolength=reclen) dummy
        pos = int((j-1)*self%nrow + i - 1, int64) * int(reclen,int64) + 1_int64
        open(newunit=unit, file=self%filename, access='stream', form='unformatted', &
             status='old', action='write', iostat=ios)
        if (ios /= 0) error stop 'fbm_set: cannot open backing file'
        write(unit, pos=pos, iostat=ios) value
        close(unit)
        if (ios /= 0) error stop 'fbm_set: write failed'
    end subroutine fbm_set

    subroutine fbm_read_col(self, j, col, rows)
        class(fbm_real), intent(in) :: self
        integer, intent(in) :: j
        real(dp), intent(out) :: col(:)
        integer, intent(in), optional :: rows(:)
        real(dp), allocatable :: full(:)
        real(dp) :: dummy
        integer :: unit, ios, reclen
        integer(int64) :: pos
        if (j < 1 .or. j > self%ncol) error stop 'fbm_read_col: bad column'
        inquire(iolength=reclen) dummy
        pos = int((j-1)*self%nrow, int64) * int(reclen,int64) + 1_int64
        open(newunit=unit, file=self%filename, access='stream', form='unformatted', &
             status='old', action='read', iostat=ios)
        if (ios /= 0) error stop 'fbm_read_col: cannot open backing file'
        if (present(rows)) then
            if (size(col) /= size(rows)) error stop 'fbm_read_col: output mismatch'
            allocate(full(self%nrow))
            read(unit, pos=pos, iostat=ios) full
            if (ios == 0) col = full(rows)
        else
            if (size(col) /= self%nrow) error stop 'fbm_read_col: output mismatch'
            read(unit, pos=pos, iostat=ios) col
        end if
        close(unit)
        if (ios /= 0) error stop 'fbm_read_col: read failed'
    end subroutine fbm_read_col

    subroutine fbm_write_col(self, j, col)
        class(fbm_real), intent(in) :: self
        integer, intent(in) :: j
        real(dp), intent(in) :: col(:)
        real(dp) :: dummy
        integer :: unit, ios, reclen
        integer(int64) :: pos
        if (size(col) /= self%nrow) error stop 'fbm_write_col: input mismatch'
        inquire(iolength=reclen) dummy
        pos = int((j-1)*self%nrow, int64) * int(reclen,int64) + 1_int64
        open(newunit=unit, file=self%filename, access='stream', form='unformatted', &
             status='old', action='write', iostat=ios)
        if (ios /= 0) error stop 'fbm_write_col: cannot open backing file'
        write(unit, pos=pos, iostat=ios) col
        close(unit)
        if (ios /= 0) error stop 'fbm_write_col: write failed'
    end subroutine fbm_write_col

    subroutine fbm_read_cols(self, j1, j2, a, rows)
        class(fbm_real), intent(in) :: self
        integer, intent(in) :: j1, j2
        real(dp), intent(out) :: a(:,:)
        integer, intent(in), optional :: rows(:)
        integer :: j
        if (j1 < 1 .or. j2 > self%ncol .or. j2 < j1) error stop 'fbm_read_cols: bad range'
        if (size(a,2) /= j2-j1+1) error stop 'fbm_read_cols: column mismatch'
        do j = j1, j2
            if (present(rows)) then
                call self%read_col(j, a(:,j-j1+1), rows)
            else
                call self%read_col(j, a(:,j-j1+1))
            end if
        end do
    end subroutine fbm_read_cols

    subroutine fbm_write_cols(self, j1, a)
        class(fbm_real), intent(in) :: self
        integer, intent(in) :: j1
        real(dp), intent(in) :: a(:,:)
        integer :: j
        if (size(a,1) /= self%nrow) error stop 'fbm_write_cols: row mismatch'
        if (j1 < 1 .or. j1 + size(a,2)-1 > self%ncol) error stop 'fbm_write_cols: range mismatch'
        do j = 1, size(a,2)
            call self%write_col(j1+j-1, a(:,j))
        end do
    end subroutine fbm_write_cols

    function fbm_to_array(self) result(a)
        class(fbm_real), intent(in) :: self
        real(dp), allocatable :: a(:,:)
        integer :: unit, ios
        allocate(a(self%nrow,self%ncol))
        open(newunit=unit, file=self%filename, access='stream', form='unformatted', &
             status='old', action='read', iostat=ios)
        if (ios /= 0) error stop 'fbm_to_array: cannot open backing file'
        read(unit, pos=1, iostat=ios) a
        close(unit)
        if (ios /= 0) error stop 'fbm_to_array: read failed'
    end function fbm_to_array

    function create_code256(filename, nrow, ncol, init) result(x)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: nrow, ncol
        integer, intent(in), optional :: init
        type(fbm_code256) :: x
        integer(int8), allocatable :: buf(:)
        integer(int8) :: val
        integer :: unit, ios, left, m
        val = 0_int8
        if (present(init)) val = int(iand(init,255), int8)
        x = attach_code256(filename,nrow,ncol)
        open(newunit=unit, file=x%filename, access='stream', form='unformatted', &
             status='replace', action='write', iostat=ios)
        if (ios /= 0) error stop 'create_code256: cannot create file'
        allocate(buf(min(max(1,nrow*ncol),1048576)))
        buf = val
        left = nrow*ncol
        do while (left > 0)
            m = min(left,size(buf))
            write(unit) buf(1:m)
            left = left-m
        end do
        close(unit)
    end function create_code256

    function attach_code256(filename,nrow,ncol) result(x)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: nrow,ncol
        type(fbm_code256) :: x
        x%filename=trim(filename)
        x%nrow=nrow
        x%ncol=ncol
    end function attach_code256

    subroutine code_read_col(self,j,col)
        class(fbm_code256), intent(in) :: self
        integer, intent(in) :: j
        integer, intent(out) :: col(:)
        integer(int8), allocatable :: raw(:)
        integer :: unit,ios,i,reclen
        integer(int8) :: dummy
        integer(int64) :: pos
        if (size(col) /= self%nrow) error stop 'code_read_col: size mismatch'
        allocate(raw(self%nrow))
        inquire(iolength=reclen) dummy
        pos=int((j-1)*self%nrow,int64)*int(reclen,int64)+1_int64
        open(newunit=unit,file=self%filename,access='stream',form='unformatted',status='old',action='read',iostat=ios)
        if (ios /= 0) error stop 'code_read_col: open failed'
        read(unit,pos=pos,iostat=ios) raw
        close(unit)
        if (ios /= 0) error stop 'code_read_col: read failed'
        do i=1,self%nrow
            col(i)=iand(int(raw(i),int32),255_int32)
        end do
    end subroutine code_read_col

    subroutine code_write_col(self,j,col)
        class(fbm_code256), intent(in) :: self
        integer, intent(in) :: j
        integer, intent(in) :: col(:)
        integer(int8), allocatable :: raw(:)
        integer :: unit,ios,i,reclen
        integer(int8) :: dummy
        integer(int64) :: pos
        if (size(col) /= self%nrow) error stop 'code_write_col: size mismatch'
        allocate(raw(self%nrow))
        do i=1,self%nrow
            if (col(i)<0 .or. col(i)>255) error stop 'code_write_col: value outside 0..255'
            raw(i)=int(iand(col(i),255),int8)
        end do
        inquire(iolength=reclen) dummy
        pos=int((j-1)*self%nrow,int64)*int(reclen,int64)+1_int64
        open(newunit=unit,file=self%filename,access='stream',form='unformatted',status='old',action='write',iostat=ios)
        if (ios /= 0) error stop 'code_write_col: open failed'
        write(unit,pos=pos,iostat=ios) raw
        close(unit)
        if (ios /= 0) error stop 'code_write_col: write failed'
    end subroutine code_write_col

    subroutine fbm_copy(src,dst)
        type(fbm_real), intent(in) :: src,dst
        real(dp), allocatable :: col(:)
        integer :: j
        if (src%nrow/=dst%nrow .or. src%ncol/=dst%ncol) error stop 'fbm_copy: dimension mismatch'
        allocate(col(src%nrow))
        do j=1,src%ncol
            call src%read_col(j,col)
            call dst%write_col(j,col)
        end do
    end subroutine fbm_copy

    subroutine fbm_increment(x,rows,cols,value)
        type(fbm_real), intent(in) :: x
        integer, intent(in) :: rows(:),cols(:)
        real(dp), intent(in) :: value
        real(dp), allocatable :: col(:)
        integer :: j
        allocate(col(x%nrow))
        do j=1,size(cols)
            call x%read_col(cols(j),col)
            col(rows)=col(rows)+value
            call x%write_col(cols(j),col)
        end do
    end subroutine fbm_increment

    subroutine fbm_transpose(src,dst,block_cols)
        type(fbm_real), intent(in) :: src,dst
        integer, intent(in), optional :: block_cols
        integer :: b,i1,i2,i,j
        real(dp), allocatable :: incol(:), out(:,:)
        if (dst%nrow/=src%ncol .or. dst%ncol/=src%nrow) error stop 'fbm_transpose: dimensions'
        b=64
        if (present(block_cols)) b=max(1,block_cols)
        allocate(incol(src%nrow))
        do i1=1,src%nrow,b
            i2=min(src%nrow,i1+b-1)
            allocate(out(src%ncol,i2-i1+1))
            do j=1,src%ncol
                call src%read_col(j,incol)
                out(j,:)=incol(i1:i2)
            end do
            do i=1,i2-i1+1
                call dst%write_col(i1+i-1,out(:,i))
            end do
            deallocate(out)
        end do
    end subroutine fbm_transpose

end module bigstatsr_fbm
