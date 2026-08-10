! Simple symmetric triplet matrix support corresponding to Rcsdp's sparse.R.
! See LICENSE (CPL-1.0).
module rcsdp_triplet
   use rcsdp_kinds, only : dp
   implicit none
   private
   public :: simple_triplet_sym_matrix, triplet_from_dense, triplet_to_dense
   public :: triplet_zero, triplet_diag

   type :: simple_triplet_sym_matrix
      integer :: n = 0
      integer, allocatable :: i(:), j(:)
      real(dp), allocatable :: v(:)
   contains
      procedure :: nnz => triplet_nnz
   end type simple_triplet_sym_matrix

contains

   pure integer function triplet_nnz(this) result(n)
      class(simple_triplet_sym_matrix), intent(in) :: this
      if (allocated(this%v)) then
         n=size(this%v)
      else
         n=0
      end if
   end function triplet_nnz

   function triplet_zero(n) result(t)
      integer, intent(in) :: n
      type(simple_triplet_sym_matrix) :: t
      t%n=n
      allocate(t%i(0),t%j(0),t%v(0))
   end function triplet_zero

   function triplet_diag(x,n) result(t)
      real(dp), intent(in) :: x(:)
      integer, intent(in), optional :: n
      type(simple_triplet_sym_matrix) :: t
      integer :: m,k
      if (present(n)) then
         m=n
      else
         m=size(x)
      end if
      t%n=m
      allocate(t%i(m),t%j(m),t%v(m))
      do k=1,m
         t%i(k)=k; t%j(k)=k; t%v(k)=x(1+mod(k-1,size(x)))
      end do
   end function triplet_diag

   function triplet_from_dense(a,tol,check_symmetric) result(t)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(in), optional :: tol
      logical, intent(in), optional :: check_symmetric
      type(simple_triplet_sym_matrix) :: t
      real(dp) :: eps
      integer :: i,j,p,n,nnz
      logical :: check
      n=size(a,1)
      if (size(a,2)/=n) error stop 'triplet_from_dense: matrix must be square'
      eps=0.0_dp; if (present(tol)) eps=tol
      check=.false.; if (present(check_symmetric)) check=check_symmetric
      if (check .and. maxval(abs(a-transpose(a)))>eps) error stop 'triplet_from_dense: matrix is not symmetric'
      nnz=0
      do j=1,n
         do i=j,n
            if (abs(a(i,j))>eps) nnz=nnz+1
         end do
      end do
      t%n=n
      allocate(t%i(nnz),t%j(nnz),t%v(nnz))
      p=0
      do j=1,n
         do i=j,n
            if (abs(a(i,j))>eps) then
               p=p+1; t%i(p)=i; t%j(p)=j; t%v(p)=a(i,j)
            end if
         end do
      end do
   end function triplet_from_dense

   function triplet_to_dense(t) result(a)
      type(simple_triplet_sym_matrix), intent(in) :: t
      real(dp), allocatable :: a(:,:)
      integer :: k,i,j
      allocate(a(t%n,t%n)); a=0.0_dp
      do k=1,t%nnz()
         i=t%i(k); j=t%j(k)
         if (i<1 .or. i>t%n .or. j<1 .or. j>t%n) error stop 'triplet_to_dense: index out of range'
         a(i,j)=a(i,j)+t%v(k)
         if (i/=j) a(j,i)=a(j,i)+t%v(k)
      end do
   end function triplet_to_dense

end module rcsdp_triplet
