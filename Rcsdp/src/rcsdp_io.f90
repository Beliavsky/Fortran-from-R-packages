! SDPA sparse problem/solution I/O translated from CSDP readprob/writeprob/writesol.
! See LICENSE (CPL-1.0).
module rcsdp_io
   use rcsdp_kinds, only : dp
   use rcsdp_types, only : csdp_problem, csdp_solution, csdp_matrix, csdp_diag
   use rcsdp_problem_mod, only : init_problem, set_sparse_a_block
   use rcsdp_block_ops, only : zero_mat
   implicit none
   private
   public :: read_sdpa_sparse, write_sdpa_sparse, read_sdpa_solution, write_sdpa_solution

   type :: entry_buffer
      integer :: n=0
      integer, allocatable :: i(:),j(:)
      real(dp), allocatable :: v(:)
   contains
      procedure :: append => buffer_append
   end type entry_buffer

   type :: token_stream
      integer :: unit = -1
      character(len=:), allocatable :: line
      integer :: pos = 1
      integer :: nchar = 0
   end type token_stream

contains

   subroutine read_sdpa_sparse(filename, prob, info)
      character(len=*), intent(in) :: filename
      type(csdp_problem), intent(out) :: prob
      integer, intent(out) :: info
      type(token_stream) :: ts
      integer :: m,nb,k,matno,blkno,ii,jj,ios
      integer, allocatable :: signed_sizes(:),cats(:),sizes(:)
      real(dp) :: val
      type(entry_buffer), allocatable :: buffers(:,:)

      info=0
      allocate(character(len=1048576) :: ts%line)
      open(newunit=ts%unit,file=filename,status='old',action='read',iostat=ios)
      if (ios/=0) then
         info=1; return
      end if
      call next_int(ts,m,ios); if (ios/=0) goto 900
      call next_int(ts,nb,ios); if (ios/=0) goto 900
      allocate(signed_sizes(nb),cats(nb),sizes(nb))
      do k=1,nb
         call next_int(ts,signed_sizes(k),ios); if (ios/=0) goto 900
         if (signed_sizes(k)<0) then
            cats(k)=csdp_diag; sizes(k)=-signed_sizes(k)
         else
            cats(k)=csdp_matrix; sizes(k)=signed_sizes(k)
         end if
      end do
      call init_problem(prob,cats,sizes,m)
      do k=1,m
         call next_real(ts,prob%b(k),ios); if (ios/=0) goto 900
      end do
      allocate(buffers(0:m,nb))
      do
         call next_int(ts,matno,ios)
         if (ios<0) exit
         if (ios/=0) goto 900
         call next_int(ts,blkno,ios); if (ios/=0) goto 900
         call next_int(ts,ii,ios); if (ios/=0) goto 900
         call next_int(ts,jj,ios); if (ios/=0) goto 900
         call next_real(ts,val,ios); if (ios/=0) goto 900
         if (matno<0 .or. matno>m .or. blkno<1 .or. blkno>nb) then
            info=2; goto 900
         end if
         call buffers(matno,blkno)%append(ii,jj,val)
      end do

      do blkno=1,nb
         if (buffers(0,blkno)%n>0) call load_c_buffer(prob,blkno,buffers(0,blkno))
      end do
      do matno=1,m
         do blkno=1,nb
            if (buffers(matno,blkno)%n>0) then
               call set_sparse_a_block(prob,matno,blkno,buffers(matno,blkno)%i(1:buffers(matno,blkno)%n), &
                  buffers(matno,blkno)%j(1:buffers(matno,blkno)%n),buffers(matno,blkno)%v(1:buffers(matno,blkno)%n))
            end if
         end do
      end do
      close(ts%unit)
      return
900   continue
      info=3
      close(ts%unit)
   end subroutine read_sdpa_sparse

   subroutine write_sdpa_sparse(filename,prob,info)
      character(len=*), intent(in) :: filename
      type(csdp_problem), intent(in) :: prob
      integer, intent(out) :: info
      integer :: u,ios,k,i,j,c,b,e,ib
      info=0
      open(newunit=u,file=filename,status='replace',action='write',iostat=ios)
      if (ios/=0) then; info=1; return; end if
      write(u,'(i0)') size(prob%b)
      write(u,'(i0)') size(prob%c%block)
      do k=1,size(prob%c%block)
         if (prob%c%block(k)%category==csdp_diag) then
            write(u,'(i0,1x)',advance='no') -prob%c%block(k)%n
         else
            write(u,'(i0,1x)',advance='no') prob%c%block(k)%n
         end if
      end do
      write(u,*)
      do k=1,size(prob%b)
         write(u,'(es25.17,1x)',advance='no') prob%b(k)
      end do
      write(u,*)
      do k=1,size(prob%c%block)
         if (prob%c%block(k)%category==csdp_diag) then
            do i=1,prob%c%block(k)%n
               if (abs(prob%c%block(k)%diag(i))>tiny(1.0_dp)) write(u,'(4(i0,1x),es25.17)') 0,k,i,i,prob%c%block(k)%diag(i)
            end do
         else
            do j=1,prob%c%block(k)%n
               do i=1,j
                  if (abs(prob%c%block(k)%mat(i,j))>tiny(1.0_dp)) write(u,'(4(i0,1x),es25.17)') 0,k,i,j,prob%c%block(k)%mat(i,j)
               end do
            end do
         end if
      end do
      do c=1,size(prob%a)
         do b=1,size(prob%a(c)%block)
            ib=prob%a(c)%block(b)%blocknum
            do e=1,prob%a(c)%block(b)%nnz()
               write(u,'(4(i0,1x),es25.17)') c,ib,prob%a(c)%block(b)%i(e),prob%a(c)%block(b)%j(e),prob%a(c)%block(b)%v(e)
            end do
         end do
      end do
      close(u)
   end subroutine write_sdpa_sparse

   subroutine read_sdpa_solution(filename,prob,sol,info)
      character(len=*), intent(in) :: filename
      type(csdp_problem), intent(in) :: prob
      type(csdp_solution), intent(out) :: sol
      integer, intent(out) :: info
      type(token_stream) :: ts
      integer :: ios,k,which,blk,ii,jj
      real(dp) :: val
      sol%x=prob%c; sol%z=prob%c
      call zero_mat(sol%x); call zero_mat(sol%z)
      allocate(sol%y(size(prob%b))); sol%y=0.0_dp
      info=0
      allocate(character(len=1048576) :: ts%line)
      open(newunit=ts%unit,file=filename,status='old',action='read',iostat=ios)
      if (ios/=0) then; info=1; return; end if
      do k=1,size(sol%y)
         call next_real(ts,sol%y(k),ios)
         if (ios/=0) then; info=2; close(ts%unit); return; end if
      end do
      do
         call next_int(ts,which,ios)
         if (ios<0) exit
         if (ios/=0) then; info=3; exit; end if
         call next_int(ts,blk,ios); if (ios/=0) then; info=3; exit; end if
         call next_int(ts,ii,ios); if (ios/=0) then; info=3; exit; end if
         call next_int(ts,jj,ios); if (ios/=0) then; info=3; exit; end if
         call next_real(ts,val,ios); if (ios/=0) then; info=3; exit; end if
         if (which==1) then
            call put_solution_entry(sol%z,blk,ii,jj,val)
         else if (which==2) then
            call put_solution_entry(sol%x,blk,ii,jj,val)
         else
            info=4; exit
         end if
      end do
      close(ts%unit)
   end subroutine read_sdpa_solution

   subroutine put_solution_entry(a,blk,ii,jj,val)
      use rcsdp_types, only : csdp_block_matrix
      type(csdp_block_matrix), intent(inout) :: a
      integer, intent(in) :: blk,ii,jj
      real(dp), intent(in) :: val
      if (blk<1 .or. blk>size(a%block)) error stop 'put_solution_entry: invalid block'
      if (a%block(blk)%category==csdp_diag) then
         a%block(blk)%diag(ii)=a%block(blk)%diag(ii)+val
      else
         a%block(blk)%mat(ii,jj)=a%block(blk)%mat(ii,jj)+val
         if (ii/=jj) a%block(blk)%mat(jj,ii)=a%block(blk)%mat(jj,ii)+val
      end if
   end subroutine put_solution_entry

   subroutine write_sdpa_solution(filename,sol,info)
      character(len=*), intent(in) :: filename
      type(csdp_solution), intent(in) :: sol
      integer, intent(out) :: info
      integer :: u,ios,k,i,j
      info=0
      open(newunit=u,file=filename,status='replace',action='write',iostat=ios)
      if (ios/=0) then; info=1; return; end if
      do i=1,size(sol%y)
         write(u,'(es25.17,1x)',advance='no') sol%y(i)
      end do
      write(u,*)
      do k=1,size(sol%z%block)
         if (sol%z%block(k)%category==csdp_diag) then
            do i=1,sol%z%block(k)%n
               if (abs(sol%z%block(k)%diag(i))>tiny(1.0_dp)) write(u,'(4(i0,1x),es25.17)') 1,k,i,i,sol%z%block(k)%diag(i)
            end do
         else
            do j=1,sol%z%block(k)%n
               do i=1,j
                  if (abs(sol%z%block(k)%mat(i,j))>tiny(1.0_dp)) write(u,'(4(i0,1x),es25.17)') 1,k,i,j,sol%z%block(k)%mat(i,j)
               end do
            end do
         end if
      end do
      do k=1,size(sol%x%block)
         if (sol%x%block(k)%category==csdp_diag) then
            do i=1,sol%x%block(k)%n
               if (abs(sol%x%block(k)%diag(i))>tiny(1.0_dp)) write(u,'(4(i0,1x),es25.17)') 2,k,i,i,sol%x%block(k)%diag(i)
            end do
         else
            do j=1,sol%x%block(k)%n
               do i=1,j
                  if (abs(sol%x%block(k)%mat(i,j))>tiny(1.0_dp)) write(u,'(4(i0,1x),es25.17)') 2,k,i,j,sol%x%block(k)%mat(i,j)
               end do
            end do
         end if
      end do
      close(u)
   end subroutine write_sdpa_solution

   subroutine buffer_append(this,ii,jj,vv)
      class(entry_buffer), intent(inout) :: this
      integer, intent(in) :: ii,jj
      real(dp), intent(in) :: vv
      integer, allocatable :: it(:),jt(:)
      real(dp), allocatable :: vt(:)
      integer :: cap,newcap
      if (.not.allocated(this%i)) then
         allocate(this%i(8),this%j(8),this%v(8))
      else if (this%n==size(this%i)) then
         cap=size(this%i); newcap=max(8,2*cap)
         allocate(it(newcap),jt(newcap),vt(newcap))
         it(1:cap)=this%i; jt(1:cap)=this%j; vt(1:cap)=this%v
         call move_alloc(it,this%i); call move_alloc(jt,this%j); call move_alloc(vt,this%v)
      end if
      this%n=this%n+1
      this%i(this%n)=ii; this%j(this%n)=jj; this%v(this%n)=vv
   end subroutine buffer_append

   subroutine load_c_buffer(prob,blk,b)
      type(csdp_problem), intent(inout) :: prob
      integer, intent(in) :: blk
      type(entry_buffer), intent(in) :: b
      integer :: e,i,j
      if (prob%c%block(blk)%category==csdp_diag) then
         do e=1,b%n
            i=b%i(e)
            prob%c%block(blk)%diag(i)=prob%c%block(blk)%diag(i)+b%v(e)
         end do
      else
         do e=1,b%n
            i=b%i(e); j=b%j(e)
            prob%c%block(blk)%mat(i,j)=prob%c%block(blk)%mat(i,j)+b%v(e)
            if (i/=j) prob%c%block(blk)%mat(j,i)=prob%c%block(blk)%mat(j,i)+b%v(e)
         end do
      end if
   end subroutine load_c_buffer

   subroutine next_int(ts,v,ios)
      type(token_stream), intent(inout) :: ts
      integer, intent(out) :: v
      integer, intent(out) :: ios
      character(len=128) :: tok
      call next_token(ts,tok,ios)
      if (ios==0) read(tok,*,iostat=ios) v
   end subroutine next_int

   subroutine next_real(ts,v,ios)
      type(token_stream), intent(inout) :: ts
      real(dp), intent(out) :: v
      integer, intent(out) :: ios
      character(len=128) :: tok
      call next_token(ts,tok,ios)
      if (ios==0) read(tok,*,iostat=ios) v
   end subroutine next_real

   subroutine next_token(ts,tok,ios)
      type(token_stream), intent(inout) :: ts
      character(len=*), intent(out) :: tok
      integer, intent(out) :: ios
      integer :: s,e,p
      character :: ch
      tok=''; ios=0
      do
         if (ts%pos>ts%nchar) then
            read(ts%unit,'(a)',iostat=ios) ts%line
            if (ios/=0) return
            ts%nchar=len_trim(ts%line); ts%pos=1
            p=1
            do while(p<=ts%nchar .and. ts%line(p:p)==' '); p=p+1; end do
            if (p>ts%nchar) cycle
            if (ts%line(p:p)=='*' .or. ts%line(p:p)=='"') cycle
         end if
         do while(ts%pos<=ts%nchar)
            ch=ts%line(ts%pos:ts%pos)
            if (index(' '//achar(9)//',{}()[]',ch)==0) exit
            ts%pos=ts%pos+1
         end do
         if (ts%pos>ts%nchar) cycle
         s=ts%pos; e=s
         do while(e<=ts%nchar)
            ch=ts%line(e:e)
            if (index(' '//achar(9)//',{}()[]',ch)>0) exit
            e=e+1
         end do
         if (e-s>len(tok)) then
            ios=2; return
         end if
         tok(1:e-s)=ts%line(s:e-1)
         ts%pos=e
         return
      end do
   end subroutine next_token

end module rcsdp_io
