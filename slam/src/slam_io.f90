! SPDX-License-Identifier: GPL-2.0-only
!
! Translation of slam's CLUTO and MC sparse matrix readers/writers.
module slam_io
    use iso_fortran_env, only : int64
    use slam_kinds, only : dp
    use slam_utils, only : argsort_int64
    use slam_stm, only : simple_triplet_matrix, make_stm
    implicit none
    private

    public :: read_stm_cluto
    public :: write_stm_cluto
    public :: read_stm_mc
    public :: write_stm_mc

contains

    function read_stm_cluto(filename) result(x)
        character(len=*), intent(in) :: filename
        type(simple_triplet_matrix) :: x
        integer :: u, ios, nr, nc, nz, r, k, nt, p
        integer, allocatable :: ii(:), jj(:)
        real(dp), allocatable :: vv(:), tok(:)
        character(len=:), allocatable :: line

        allocate(character(len=1048576) :: line)
        open(newunit=u,file=filename,status='old',action='read',iostat=ios)
        if(ios/=0) error stop "read_stm_cluto: cannot open file"
        read(u,*,iostat=ios) nr,nc,nz
        if(ios/=0 .or. nr<0 .or. nc<0 .or. nz<0) error stop "read_stm_cluto: invalid header"
        allocate(ii(nz),jj(nz),vv(nz),tok(max(1,2*nc)))
        k=0
        do r=1,nr
            read(u,'(A)',iostat=ios) line
            if(ios/=0) error stop "read_stm_cluto: unexpected end of file"
            call parse_real_tokens(line,tok,nt)
            if(mod(nt,2)/=0) error stop "read_stm_cluto: odd number of row tokens"
            do p=1,nt,2
                k=k+1
                if(k>nz) error stop "read_stm_cluto: more nonzeros than header"
                ii(k)=r
                jj(k)=int(tok(p))
                vv(k)=tok(p+1)
            end do
        end do
        close(u)
        if(k/=nz) error stop "read_stm_cluto: nonzero count differs from header"
        x=make_stm(ii,jj,vv,nr,nc)
    end function read_stm_cluto

    subroutine write_stm_cluto(x,filename)
        type(simple_triplet_matrix), intent(in) :: x
        character(len=*), intent(in) :: filename
        integer :: u,ios,r,k

        open(newunit=u,file=filename,status='replace',action='write',iostat=ios)
        if(ios/=0) error stop "write_stm_cluto: cannot open file"
        write(u,'(i0,1x,i0,1x,i0)') x%nrow,x%ncol,x%nnz()
        do r=1,x%nrow
            do k=1,x%nnz()
                if(x%i(k)==r) write(u,'(i0,1x,es24.16,1x)',advance='no') x%j(k),x%v(k)
            end do
            write(u,'()')
        end do
        close(u)
    end subroutine write_stm_cluto

    function read_stm_mc(prefix,scalingtype) result(x)
        character(len=*), intent(in) :: prefix
        character(len=*), intent(in), optional :: scalingtype
        type(simple_triplet_matrix) :: x
        character(len=:), allocatable :: nzfile,scale
        integer :: u,ios,nr,nc,nz,k,c
        integer, allocatable :: rows(:),ptr(:),ii(:),jj(:)
        real(dp), allocatable :: vv(:)

        open(newunit=u,file=trim(prefix)//'_dim',status='old',action='read',iostat=ios)
        if(ios/=0) error stop "read_stm_mc: cannot open _dim"
        read(u,*,iostat=ios) nr,nc,nz
        close(u)
        if(ios/=0) error stop "read_stm_mc: invalid _dim"
        allocate(rows(nz),ptr(nc+1),ii(nz),jj(nz),vv(nz))
        open(newunit=u,file=trim(prefix)//'_row_ccs',status='old',action='read',iostat=ios)
        if(ios/=0) error stop "read_stm_mc: cannot open _row_ccs"
        do k=1,nz
            read(u,*,iostat=ios) rows(k)
            if(ios/=0) error stop "read_stm_mc: invalid _row_ccs"
        end do
        close(u)
        open(newunit=u,file=trim(prefix)//'_col_ccs',status='old',action='read',iostat=ios)
        if(ios/=0) error stop "read_stm_mc: cannot open _col_ccs"
        do k=1,nc+1
            read(u,*,iostat=ios) ptr(k)
            if(ios/=0) error stop "read_stm_mc: invalid _col_ccs"
        end do
        close(u)
        if(present(scalingtype)) then
            scale=trim(scalingtype)
            nzfile=trim(prefix)//'_'//scale//'_nz'
            if(.not.file_exists(nzfile)) error stop "read_stm_mc: requested scaling file not found"
        else
            call locate_mc_nz_file(prefix,nzfile,scale)
        end if
        open(newunit=u,file=nzfile,status='old',action='read',iostat=ios)
        if(ios/=0) error stop "read_stm_mc: cannot open nonzero file"
        do k=1,nz
            read(u,*,iostat=ios) vv(k)
            if(ios/=0) error stop "read_stm_mc: invalid nonzero file"
        end do
        close(u)
        ii=rows+1
        k=0
        do c=1,nc
            do while(k < ptr(c+1))
                k=k+1
                if(k>nz) error stop "read_stm_mc: invalid column pointers"
                jj(k)=c
            end do
        end do
        if(k/=nz) error stop "read_stm_mc: invalid column pointers"
        x=make_stm(ii,jj,vv,nr,nc)
    end function read_stm_mc

    subroutine write_stm_mc(x, prefix)
        type(simple_triplet_matrix), intent(in) :: x
        character(len=*), intent(in) :: prefix
        type(simple_triplet_matrix) :: t
        integer(int64), allocatable :: key(:)
        integer, allocatable :: ord(:), counts(:),ptr(:)
        integer :: u,ios,k

        t=x%transpose()
        allocate(key(t%nnz()),counts(t%ncol),ptr(0:t%ncol))
        key=int(t%i,int64)+int(t%j-1,int64)*int(max(1,t%nrow),int64)
        call argsort_int64(key,ord)
        open(newunit=u,file=trim(prefix)//'_dim',status='replace',action='write',iostat=ios)
        if(ios/=0) error stop "write_stm_mc: cannot open _dim"
        write(u,'(i0,1x,i0,1x,i0)') t%nrow,t%ncol,t%nnz()
        close(u)
        open(newunit=u,file=trim(prefix)//'_row_ccs',status='replace',action='write',iostat=ios)
        if(ios/=0) error stop "write_stm_mc: cannot open _row_ccs"
        do k=1,t%nnz()
            write(u,'(i0)') t%i(ord(k))-1
        end do
        close(u)
        counts=0
        do k=1,t%nnz()
            counts(t%j(ord(k)))=counts(t%j(ord(k)))+1
        end do
        ptr=0
        do k=1,t%ncol
            ptr(k)=ptr(k-1)+counts(k)
        end do
        open(newunit=u,file=trim(prefix)//'_col_ccs',status='replace',action='write',iostat=ios)
        if(ios/=0) error stop "write_stm_mc: cannot open _col_ccs"
        do k=0,t%ncol
            write(u,'(i0)') ptr(k)
        end do
        close(u)
        open(newunit=u,file=trim(prefix)//'_tfn_nz',status='replace',action='write',iostat=ios)
        if(ios/=0) error stop "write_stm_mc: cannot open _tfn_nz"
        do k=1,t%nnz()
            write(u,'(es24.16)') t%v(ord(k))
        end do
        close(u)
    end subroutine write_stm_mc

    subroutine parse_real_tokens(line,values,n)
        character(len=*), intent(in) :: line
        real(dp), intent(out) :: values(:)
        integer, intent(out) :: n
        integer :: p,q,l,ios
        character(len=128) :: token
        character :: ch

        n=0; l=len_trim(line); p=1
        do while(p<=l)
            do while(p<=l)
                ch=line(p:p)
                if(ch/=' ' .and. ch/=char(9)) exit
                p=p+1
            end do
            if(p>l) exit
            q=p
            do while(q<=l)
                ch=line(q:q)
                if(ch==' ' .or. ch==char(9)) exit
                q=q+1
            end do
            if(q-p>len(token)) error stop "parse_real_tokens: token too long"
            token=' '
            token(1:q-p)=line(p:q-1)
            n=n+1
            if(n>size(values)) error stop "parse_real_tokens: too many tokens"
            read(token,*,iostat=ios) values(n)
            if(ios/=0) error stop "parse_real_tokens: invalid numeric token"
            p=q+1
        end do
    end subroutine parse_real_tokens

    subroutine locate_mc_nz_file(prefix,filename,scale)
        character(len=*), intent(in) :: prefix
        character(len=:), allocatable, intent(out) :: filename,scale
        character(len=1), parameter :: a1(2)=[character(len=1)::'t','l']
        character(len=1), parameter :: a2(4)=[character(len=1)::'x','f','e','1']
        character(len=1), parameter :: a3(3)=[character(len=1)::'x','n','1']
        character(len=1), parameter :: a4(2)=[character(len=1)::' ','i']
        character(len=4) :: s
        character(len=:), allocatable :: f
        logical :: exists
        integer :: i,j,k,l
        do i=1,2; do j=1,4; do k=1,3; do l=1,2
            if(l==1) then
                s= a1(i)//a2(j)//a3(k)//' '
                scale=trim(s(:3))
            else
                s= a1(i)//a2(j)//a3(k)//a4(l)
                scale=trim(s)
            end if
            f=trim(prefix)//'_'//scale//'_nz'
            inquire(file=f,exist=exists)
            if(exists) then
                filename=f
                return
            end if
        end do; end do; end do; end do
        error stop "read_stm_mc: no recognized nonzero file found"
    end subroutine locate_mc_nz_file

    logical function file_exists(filename)
        character(len=*),intent(in)::filename
        inquire(file=filename,exist=file_exists)
    end function file_exists


end module slam_io
