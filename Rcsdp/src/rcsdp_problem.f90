! Problem construction and sparse constraint operators for Rcsdp/CSDP.
! See LICENSE (CPL-1.0).
module rcsdp_problem_mod
   use rcsdp_kinds, only : dp
   use rcsdp_types, only : csdp_problem, csdp_block_matrix, csdp_constraint, csdp_sparse_block, &
      csdp_matrix, csdp_diag
   use rcsdp_block_ops, only : allocate_block_matrix, zero_mat, mat_mult, trace_prod
   implicit none
   private
   public :: init_problem, set_c_matrix_block, set_c_diag_block
   public :: set_a_matrix_block, set_a_diag_block, set_sparse_a_block
   public :: dense_constraint, build_dense_constraints, op_a, op_at, op_o
   public :: validate_problem

contains

   subroutine init_problem(prob, categories, sizes, m)
      type(csdp_problem), intent(out) :: prob
      integer, intent(in) :: categories(:), sizes(:)
      integer, intent(in) :: m
      integer :: i
      call allocate_block_matrix(categories,sizes,prob%c)
      allocate(prob%a(m),prob%b(m))
      prob%b = 0.0_dp
      do i = 1, m
         allocate(prob%a(i)%block(0))
      end do
   end subroutine init_problem

   subroutine set_c_matrix_block(prob, iblock, x)
      type(csdp_problem), intent(inout) :: prob
      integer, intent(in) :: iblock
      real(dp), intent(in) :: x(:,:)
      if (prob%c%block(iblock)%category /= csdp_matrix) error stop 'set_c_matrix_block: block is not matrix'
      if (size(x,1) /= prob%c%block(iblock)%n .or. size(x,2) /= prob%c%block(iblock)%n) &
         error stop 'set_c_matrix_block: wrong shape'
      prob%c%block(iblock)%mat = 0.5_dp*(x+transpose(x))
   end subroutine set_c_matrix_block

   subroutine set_c_diag_block(prob, iblock, x)
      type(csdp_problem), intent(inout) :: prob
      integer, intent(in) :: iblock
      real(dp), intent(in) :: x(:)
      if (prob%c%block(iblock)%category /= csdp_diag) error stop 'set_c_diag_block: block is not diagonal'
      if (size(x) /= prob%c%block(iblock)%n) error stop 'set_c_diag_block: wrong shape'
      prob%c%block(iblock)%diag = x
   end subroutine set_c_diag_block

   subroutine set_a_matrix_block(prob, icon, iblock, x, tol)
      type(csdp_problem), intent(inout) :: prob
      integer, intent(in) :: icon, iblock
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(in), optional :: tol
      integer :: i, j, n, nnz, p
      real(dp) :: eps
      integer, allocatable :: ii(:), jj(:)
      real(dp), allocatable :: vv(:)
      if (prob%c%block(iblock)%category /= csdp_matrix) error stop 'set_a_matrix_block: block is not matrix'
      n = prob%c%block(iblock)%n
      if (size(x,1) /= n .or. size(x,2) /= n) error stop 'set_a_matrix_block: wrong shape'
      eps = 0.0_dp
      if (present(tol)) eps = tol
      nnz = 0
      do j = 1, n
         do i = 1, j
            if (abs(x(i,j)) > eps) nnz = nnz + 1
         end do
      end do
      allocate(ii(nnz),jj(nnz),vv(nnz))
      p = 0
      do j = 1, n
         do i = 1, j
            if (abs(x(i,j)) > eps) then
               p = p + 1
               ii(p) = i
               jj(p) = j
               vv(p) = x(i,j)
            end if
         end do
      end do
      call set_sparse_a_block(prob,icon,iblock,ii,jj,vv)
   end subroutine set_a_matrix_block

   subroutine set_a_diag_block(prob, icon, iblock, x, tol)
      type(csdp_problem), intent(inout) :: prob
      integer, intent(in) :: icon, iblock
      real(dp), intent(in) :: x(:)
      real(dp), intent(in), optional :: tol
      integer :: i, nnz, p
      real(dp) :: eps
      integer, allocatable :: ii(:), jj(:)
      real(dp), allocatable :: vv(:)
      if (prob%c%block(iblock)%category /= csdp_diag) error stop 'set_a_diag_block: block is not diagonal'
      if (size(x) /= prob%c%block(iblock)%n) error stop 'set_a_diag_block: wrong shape'
      eps = 0.0_dp
      if (present(tol)) eps = tol
      nnz = count(abs(x) > eps)
      allocate(ii(nnz),jj(nnz),vv(nnz))
      p = 0
      do i = 1, size(x)
         if (abs(x(i)) > eps) then
            p = p + 1
            ii(p)=i
            jj(p)=i
            vv(p)=x(i)
         end if
      end do
      call set_sparse_a_block(prob,icon,iblock,ii,jj,vv)
   end subroutine set_a_diag_block

   subroutine set_sparse_a_block(prob, icon, iblock, ii, jj, vv)
      type(csdp_problem), intent(inout) :: prob
      integer, intent(in) :: icon, iblock
      integer, intent(in) :: ii(:), jj(:)
      real(dp), intent(in) :: vv(:)
      type(csdp_sparse_block) :: sb
      type(csdp_sparse_block), allocatable :: tmp(:)
      integer :: p, loc
      if (size(ii) /= size(jj) .or. size(ii) /= size(vv)) error stop 'set_sparse_a_block: inconsistent arrays'
      if (icon < 1 .or. icon > size(prob%a)) error stop 'set_sparse_a_block: invalid constraint'
      if (iblock < 1 .or. iblock > size(prob%c%block)) error stop 'set_sparse_a_block: invalid block'
      do p=1,size(ii)
         if (ii(p)<1 .or. jj(p)<1 .or. ii(p)>prob%c%block(iblock)%n .or. jj(p)>prob%c%block(iblock)%n) &
            error stop 'set_sparse_a_block: index out of range'
      end do
      sb%blocknum=iblock
      sb%n=prob%c%block(iblock)%n
      allocate(sb%i(size(ii)),sb%j(size(jj)),sb%v(size(vv)))
      sb%i=ii; sb%j=jj; sb%v=vv
      loc=0
      do p=1,size(prob%a(icon)%block)
         if (prob%a(icon)%block(p)%blocknum==iblock) then
            loc=p
            exit
         end if
      end do
      if (loc>0) then
         prob%a(icon)%block(loc)=sb
      else
         allocate(tmp(size(prob%a(icon)%block)+1))
         if (size(prob%a(icon)%block)>0) tmp(1:size(prob%a(icon)%block))=prob%a(icon)%block
         tmp(size(tmp))=sb
         call move_alloc(tmp,prob%a(icon)%block)
      end if
   end subroutine set_sparse_a_block

   subroutine dense_constraint(prob, icon, a)
      type(csdp_problem), intent(in) :: prob
      integer, intent(in) :: icon
      type(csdp_block_matrix), intent(out) :: a
      integer :: b, e, ib, i, j
      a = prob%c
      call zero_mat(a)
      do b = 1, size(prob%a(icon)%block)
         ib=prob%a(icon)%block(b)%blocknum
         do e=1,prob%a(icon)%block(b)%nnz()
            i=prob%a(icon)%block(b)%i(e)
            j=prob%a(icon)%block(b)%j(e)
            if (a%block(ib)%category==csdp_diag) then
               a%block(ib)%diag(i)=a%block(ib)%diag(i)+prob%a(icon)%block(b)%v(e)
            else
               a%block(ib)%mat(i,j)=a%block(ib)%mat(i,j)+prob%a(icon)%block(b)%v(e)
               if (i/=j) a%block(ib)%mat(j,i)=a%block(ib)%mat(j,i)+prob%a(icon)%block(b)%v(e)
            end if
         end do
      end do
   end subroutine dense_constraint

   subroutine build_dense_constraints(prob, adense)
      type(csdp_problem), intent(in) :: prob
      type(csdp_block_matrix), allocatable, intent(out) :: adense(:)
      integer :: i
      allocate(adense(size(prob%a)))
      do i=1,size(prob%a)
         call dense_constraint(prob,i,adense(i))
      end do
   end subroutine build_dense_constraints

   subroutine op_a(prob, x, result, adense)
      type(csdp_problem), intent(in) :: prob
      type(csdp_block_matrix), intent(in) :: x
      real(dp), intent(out) :: result(:)
      type(csdp_block_matrix), intent(in), optional :: adense(:)
      type(csdp_block_matrix) :: ai
      integer :: i
      if (size(result)/=size(prob%a)) error stop 'op_a: result has wrong size'
      if (present(adense)) then
         do i=1,size(prob%a)
            result(i)=trace_prod(adense(i),x)
         end do
      else
         do i=1,size(prob%a)
            call dense_constraint(prob,i,ai)
            result(i)=trace_prod(ai,x)
         end do
      end if
   end subroutine op_a

   subroutine op_at(prob, y, result, adense)
      type(csdp_problem), intent(in) :: prob
      real(dp), intent(in) :: y(:)
      type(csdp_block_matrix), intent(out) :: result
      type(csdp_block_matrix), intent(in), optional :: adense(:)
      type(csdp_block_matrix) :: ai
      integer :: i, k
      if (size(y)/=size(prob%a)) error stop 'op_at: y has wrong size'
      result=prob%c
      call zero_mat(result)
      if (present(adense)) then
         do i=1,size(y)
            do k=1,size(result%block)
               if (result%block(k)%category==csdp_diag) then
                  result%block(k)%diag=result%block(k)%diag+y(i)*adense(i)%block(k)%diag
               else
                  result%block(k)%mat=result%block(k)%mat+y(i)*adense(i)%block(k)%mat
               end if
            end do
         end do
      else
         do i=1,size(y)
            call dense_constraint(prob,i,ai)
            do k=1,size(result%block)
               if (result%block(k)%category==csdp_diag) then
                  result%block(k)%diag=result%block(k)%diag+y(i)*ai%block(k)%diag
               else
                  result%block(k)%mat=result%block(k)%mat+y(i)*ai%block(k)%mat
               end if
            end do
         end do
      end if
   end subroutine op_at

   subroutine op_o(prob, zi, x, o, adense)
      type(csdp_problem), intent(in) :: prob
      type(csdp_block_matrix), intent(in) :: zi, x
      real(dp), intent(out) :: o(:,:)
      type(csdp_block_matrix), intent(in) :: adense(:)
      type(csdp_block_matrix) :: w1,w2
      real(dp), allocatable :: col(:)
      integer :: j
      if (size(o,1)/=size(prob%a) .or. size(o,2)/=size(prob%a)) error stop 'op_o: wrong shape'
      allocate(col(size(prob%a)))
      w1=prob%c; w2=prob%c
      do j=1,size(prob%a)
         call mat_mult(1.0_dp,0.0_dp,zi,adense(j),w1)
         call mat_mult(1.0_dp,0.0_dp,w1,x,w2)
         call op_a(prob,w2,col,adense)
         o(:,j)=col
      end do
      o=0.5_dp*(o+transpose(o))
   end subroutine op_o

   subroutine validate_problem(prob, ok, message)
      type(csdp_problem), intent(in) :: prob
      logical, intent(out) :: ok
      character(len=:), allocatable, intent(out) :: message
      integer :: i,b,e,ib
      ok=.false.
      if (.not. allocated(prob%c%block)) then
         message='C has no blocks'; return
      end if
      if (.not. allocated(prob%a) .or. .not. allocated(prob%b)) then
         message='constraints or b are missing'; return
      end if
      if (size(prob%a)/=size(prob%b)) then
         message='number of constraints differs from size of b'; return
      end if
      do i=1,size(prob%a)
         do b=1,size(prob%a(i)%block)
            ib=prob%a(i)%block(b)%blocknum
            if (ib<1 .or. ib>size(prob%c%block)) then
               message='constraint block number out of range'; return
            end if
            do e=1,prob%a(i)%block(b)%nnz()
               if (prob%a(i)%block(b)%i(e)<1 .or. prob%a(i)%block(b)%i(e)>prob%c%block(ib)%n .or. &
                   prob%a(i)%block(b)%j(e)<1 .or. prob%a(i)%block(b)%j(e)>prob%c%block(ib)%n) then
                  message='constraint index out of range'; return
               end if
            end do
         end do
      end do
      ok=.true.
      message='ok'
   end subroutine validate_problem

end module rcsdp_problem_mod
