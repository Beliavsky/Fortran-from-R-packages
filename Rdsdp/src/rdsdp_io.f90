! SDPA sparse-format I/O compatible with the Rdsdp reader/writer conventions.
! v0.2.0 preserves sparse SDP matrices without dense n*n*m materialization.
! DSDP copyright/license: see licenses/DSDP-LICENSE.
module rdsdp_io
   use rdsdp_kinds, only : dp
   use rdsdp_types, only : dsdp_problem, dsdp_sdp_block, dsdp_lp_block, dsdp_solution, &
      dsdp_data_dense, dsdp_data_sparse
   use rdsdp_data, only : get_data_dense
   implicit none
   private
   public :: read_sdpa, write_sdpa, flatten_primal

contains

   subroutine read_sdpa(filename,prob)
      character(len=*), intent(in) :: filename
      type(dsdp_problem), intent(out) :: prob
      integer :: u,ios,m,nb,k,nmat,nblk,row,col,add,p
      integer, allocatable :: bs(:),cnt0(:),cnta(:,:),pos0(:),posa(:,:)
      integer, allocatable :: bs2(:)
      real(dp), allocatable :: b(:),b2(:)
      real(dp) :: val
      character(len=16384) :: line

      open(newunit=u,file=filename,status='old',action='read',iostat=ios)
      if (ios/=0) error stop 'read_sdpa: cannot open file'
      call read_header(u,m,nb,bs,b)
      allocate(cnt0(nb),cnta(nb,m)); cnt0=0; cnta=0
      do
         read(u,'(A)',iostat=ios) line
         if (ios<0) exit
         if (ios>0) error stop 'read_sdpa: I/O failure'
         if (skip_line(line)) cycle
         call sanitize_delimiters(line)
         read(line,*,iostat=ios) nmat,nblk,row,col,val
         if (ios/=0) cycle
         if (nblk<1 .or. nblk>nb) error stop 'read_sdpa: block index out of range'
         if (nmat<0 .or. nmat>m) error stop 'read_sdpa: matrix index out of range'
         if (bs(nblk)>0 .and. abs(val)>tiny(1.0_dp)) then
            add=1+merge(1,0,row/=col)
            if (nmat==0) then; cnt0(nblk)=cnt0(nblk)+add; else; cnta(nblk,nmat)=cnta(nblk,nmat)+add; end if
         end if
      end do
      close(u)

      prob%m=m; allocate(prob%b(m),prob%block(nb)); prob%b=b
      do k=1,nb
         if (bs(k)>0) then
            prob%block(k)%category=dsdp_sdp_block; prob%block(k)%n=bs(k)
            prob%block(k)%c_storage=dsdp_data_sparse
            allocate(prob%block(k)%a_storage(m),prob%block(k)%a_sparse(m),prob%block(k)%a_lowrank(m))
            prob%block(k)%a_storage=dsdp_data_sparse
            prob%block(k)%c_sparse%n=bs(k); prob%block(k)%c_sparse%nnz=cnt0(k)
            allocate(prob%block(k)%c_sparse%row(cnt0(k)),prob%block(k)%c_sparse%col(cnt0(k)), &
                     prob%block(k)%c_sparse%val(cnt0(k)))
            do nmat=1,m
               prob%block(k)%a_sparse(nmat)%n=bs(k); prob%block(k)%a_sparse(nmat)%nnz=cnta(k,nmat)
               allocate(prob%block(k)%a_sparse(nmat)%row(cnta(k,nmat)), &
                        prob%block(k)%a_sparse(nmat)%col(cnta(k,nmat)), &
                        prob%block(k)%a_sparse(nmat)%val(cnta(k,nmat)))
            end do
         else if (bs(k)<0) then
            prob%block(k)%category=dsdp_lp_block; prob%block(k)%n=-bs(k)
            allocate(prob%block(k)%cdiag(-bs(k)),prob%block(k)%adiag(-bs(k),m))
            prob%block(k)%cdiag=0.0_dp; prob%block(k)%adiag=0.0_dp
         else
            error stop 'read_sdpa: zero block size not supported'
         end if
      end do

      allocate(pos0(nb),posa(nb,m)); pos0=0; posa=0
      open(newunit=u,file=filename,status='old',action='read',iostat=ios)
      if (ios/=0) error stop 'read_sdpa: cannot reopen file'
      call read_header(u,nmat,nblk,bs2,b2)
      if (nmat/=m .or. nblk/=nb) error stop 'read_sdpa: header changed while reading'
      do
         read(u,'(A)',iostat=ios) line
         if (ios<0) exit
         if (ios>0) error stop 'read_sdpa: I/O failure'
         if (skip_line(line)) cycle
         call sanitize_delimiters(line)
         read(line,*,iostat=ios) nmat,nblk,row,col,val
         if (ios/=0 .or. abs(val)<=tiny(1.0_dp)) cycle
         if (nmat==0) val=-val  ! Rreadsdpa.c convention
         if (prob%block(nblk)%category==dsdp_sdp_block) then
            if (row<1 .or. row>prob%block(nblk)%n .or. col<1 .or. col>prob%block(nblk)%n) &
               error stop 'read_sdpa: SDP index out of range'
            if (nmat==0) then
               p=pos0(nblk)+1; prob%block(nblk)%c_sparse%row(p)=row; prob%block(nblk)%c_sparse%col(p)=col; &
                  prob%block(nblk)%c_sparse%val(p)=val; pos0(nblk)=p
               if (row/=col) then
                  p=pos0(nblk)+1; prob%block(nblk)%c_sparse%row(p)=col; prob%block(nblk)%c_sparse%col(p)=row; &
                     prob%block(nblk)%c_sparse%val(p)=val; pos0(nblk)=p
               end if
            else
               p=posa(nblk,nmat)+1; prob%block(nblk)%a_sparse(nmat)%row(p)=row; &
                  prob%block(nblk)%a_sparse(nmat)%col(p)=col; prob%block(nblk)%a_sparse(nmat)%val(p)=val; posa(nblk,nmat)=p
               if (row/=col) then
                  p=posa(nblk,nmat)+1; prob%block(nblk)%a_sparse(nmat)%row(p)=col; &
                     prob%block(nblk)%a_sparse(nmat)%col(p)=row; prob%block(nblk)%a_sparse(nmat)%val(p)=val; posa(nblk,nmat)=p
               end if
            end if
         else
            if (row<1 .or. row>prob%block(nblk)%n) error stop 'read_sdpa: LP index out of range'
            if (nmat==0) then; prob%block(nblk)%cdiag(row)=prob%block(nblk)%cdiag(row)+val
            else; prob%block(nblk)%adiag(row,nmat)=prob%block(nblk)%adiag(row,nmat)+val; end if
         end if
      end do
      close(u)
   end subroutine read_sdpa

   subroutine write_sdpa(filename,prob)
      character(len=*), intent(in) :: filename
      type(dsdp_problem), intent(in) :: prob
      integer :: u,k,i,j,q,e
      real(dp), allocatable :: a(:,:)
      open(newunit=u,file=filename,status='replace',action='write')
      write(u,'(i0)') prob%m
      write(u,'(i0)') size(prob%block)
      do k=1,size(prob%block)
         if (prob%block(k)%category==dsdp_sdp_block) then
            write(u,'(i0,1x)',advance='no') prob%block(k)%n
         else
            write(u,'(i0,1x)',advance='no') -prob%block(k)%n
         end if
      end do
      write(u,*)
      write(u,'(*(es24.16,1x))') prob%b
      do k=1,size(prob%block)
         if (prob%block(k)%category==dsdp_sdp_block) then
            if (prob%block(k)%c_storage==dsdp_data_sparse) then
               do e=1,prob%block(k)%c_sparse%nnz
                  i=prob%block(k)%c_sparse%row(e); j=prob%block(k)%c_sparse%col(e)
                  if (i<=j .and. abs(prob%block(k)%c_sparse%val(e))>tiny(1.0_dp)) &
                     write(u,'(4(i0,1x),es24.16)') 0,k,i,j,-prob%block(k)%c_sparse%val(e)
               end do
            else
               call get_data_dense(prob%block(k),0,a)
               do i=1,prob%block(k)%n; do j=i,prob%block(k)%n
                  if (abs(a(i,j))>tiny(1.0_dp)) write(u,'(4(i0,1x),es24.16)') 0,k,i,j,-a(i,j)
               end do; end do
            end if
            do q=1,prob%m
               if (allocated(prob%block(k)%a_storage) .and. prob%block(k)%a_storage(q)==dsdp_data_sparse) then
                  do e=1,prob%block(k)%a_sparse(q)%nnz
                     i=prob%block(k)%a_sparse(q)%row(e); j=prob%block(k)%a_sparse(q)%col(e)
                     if (i<=j .and. abs(prob%block(k)%a_sparse(q)%val(e))>tiny(1.0_dp)) &
                        write(u,'(4(i0,1x),es24.16)') q,k,i,j,prob%block(k)%a_sparse(q)%val(e)
                  end do
               else
                  call get_data_dense(prob%block(k),q,a)
                  do i=1,prob%block(k)%n; do j=i,prob%block(k)%n
                     if (abs(a(i,j))>tiny(1.0_dp)) write(u,'(4(i0,1x),es24.16)') q,k,i,j,a(i,j)
                  end do; end do
               end if
            end do
         else
            do i=1,prob%block(k)%n
               if (abs(prob%block(k)%cdiag(i))>tiny(1.0_dp)) write(u,'(4(i0,1x),es24.16)') 0,k,i,i,-prob%block(k)%cdiag(i)
            end do
            do q=1,prob%m
               do i=1,prob%block(k)%n
                  if (abs(prob%block(k)%adiag(i,q))>tiny(1.0_dp)) write(u,'(4(i0,1x),es24.16)') &
                     q,k,i,i,prob%block(k)%adiag(i,q)
               end do
            end do
         end if
      end do
      close(u)
   end subroutine write_sdpa

   subroutine flatten_primal(sol,x)
      type(dsdp_solution), intent(in) :: sol
      real(dp), allocatable, intent(out) :: x(:)
      integer :: nall,k,i,j,p
      nall=0
      do k=1,size(sol%x)
         if (sol%x(k)%category==dsdp_sdp_block) then; nall=nall+sol%x(k)%n*sol%x(k)%n
         else; nall=nall+sol%x(k)%n; end if
      end do
      allocate(x(nall)); p=0
      do k=1,size(sol%x)
         if (sol%x(k)%category==dsdp_sdp_block) then
            do i=1,sol%x(k)%n; do j=1,sol%x(k)%n
               p=p+1; x(p)=sol%x(k)%x(i,j)
            end do; end do
         else
            x(p+1:p+sol%x(k)%n)=sol%x(k)%xdiag; p=p+sol%x(k)%n
         end if
      end do
   end subroutine flatten_primal

   subroutine read_header(u,m,nb,bs,b)
      integer, intent(in) :: u
      integer, intent(out) :: m,nb
      integer, allocatable, intent(out) :: bs(:)
      real(dp), allocatable, intent(out) :: b(:)
      integer :: ios
      character(len=16384) :: line
      logical :: got
      call next_data_line(u,line,got); if (.not.got) error stop 'read_sdpa: missing m'
      read(line,*,iostat=ios) m; if (ios/=0 .or. m<1) error stop 'read_sdpa: invalid m'
      call next_data_line(u,line,got); if (.not.got) error stop 'read_sdpa: missing block count'
      read(line,*,iostat=ios) nb; if (ios/=0 .or. nb<1) error stop 'read_sdpa: invalid block count'
      allocate(bs(nb))
      call next_data_line(u,line,got); if (.not.got) error stop 'read_sdpa: missing block sizes'
      call sanitize_delimiters(line); read(line,*,iostat=ios) bs
      if (ios/=0) error stop 'read_sdpa: invalid block sizes'
      allocate(b(m))
      call next_data_line(u,line,got); if (.not.got) error stop 'read_sdpa: missing objective vector b'
      call sanitize_delimiters(line); read(line,*,iostat=ios) b
      if (ios/=0) error stop 'read_sdpa: invalid objective vector b'
   end subroutine read_header

   subroutine next_data_line(u,line,got)
      integer, intent(in) :: u
      character(len=*), intent(out) :: line
      logical, intent(out) :: got
      integer :: ios
      got=.false.; line=''
      do
         read(u,'(A)',iostat=ios) line
         if (ios/=0) return
         if (.not.skip_line(line)) then; got=.true.; return; end if
      end do
   end subroutine next_data_line

   logical function skip_line(line) result(skip)
      character(len=*), intent(in) :: line
      character(len=:), allocatable :: t
      t=adjustl(line); skip=len_trim(t)==0
      if (.not.skip) skip=t(1:1)=='*' .or. t(1:1)=='"' .or. t(1:1)=='#' .or. t(1:1)=='%'
   end function skip_line

   subroutine sanitize_delimiters(line)
      character(len=*), intent(inout) :: line
      integer :: i
      do i=1,len_trim(line)
         select case(line(i:i)); case('{','}','(',')',',',';'); line(i:i)=' '; end select
      end do
   end subroutine sanitize_delimiters

end module rdsdp_io
