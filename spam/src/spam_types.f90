module spam_types
use spam_kinds, only: dp
implicit none
private
public :: csr_matrix, spam_chol, eigen_result, mle_result

type :: csr_matrix
   integer :: nrow = 0
   integer :: ncol = 0
   real(dp), allocatable :: entries(:)
   integer, allocatable :: colindices(:)
   integer, allocatable :: rowpointers(:)
contains
   procedure :: nnz => csr_nnz
   procedure :: valid => csr_valid
end type csr_matrix

type :: spam_chol
   integer :: n = 0
   integer :: nsuper = 0
   integer :: nnz_a = 0
   real(dp), allocatable :: entries(:)
   integer, allocatable :: colindices(:)
   integer, allocatable :: colpointers(:)
   integer, allocatable :: rowpointers(:)
   integer, allocatable :: pivot(:)
   integer, allocatable :: invpivot(:)
   integer, allocatable :: supernodes(:)
   integer, allocatable :: snmember(:)
   integer :: cache_kb = 512
   integer :: info = 0
end type spam_chol

type :: eigen_result
   real(dp), allocatable :: values(:)
   real(dp), allocatable :: vectors(:,:)
   real(dp), allocatable :: imag_values(:)
   integer :: nconv = 0
   integer :: niter = 0
   integer :: info = 0
end type eigen_result

type :: mle_result
   real(dp), allocatable :: par(:)
   real(dp) :: objective = huge(1.0_dp)
   integer :: convergence = 1
   integer :: iterations = 0
end type mle_result

contains
integer function csr_nnz(self) result(n)
class(csr_matrix), intent(in) :: self
if (allocated(self%entries)) then
   n = size(self%entries)
else
   n = 0
end if
end function csr_nnz

logical function csr_valid(self) result(ok)
class(csr_matrix), intent(in) :: self
integer :: i, k
ok = self%nrow >= 0 .and. self%ncol >= 0
if (.not. allocated(self%rowpointers)) then
   ok = ok .and. self%nrow == 0
   return
end if
if (size(self%rowpointers) /= self%nrow + 1) then
   ok = .false.
   return
end if
if (self%rowpointers(1) /= 1) then
   ok = .false.
   return
end if
if (.not. allocated(self%entries) .or. .not. allocated(self%colindices)) then
   ok = self%rowpointers(self%nrow+1) == 1
   return
end if
if (size(self%entries) /= size(self%colindices)) then
   ok = .false.
   return
end if
if (self%rowpointers(self%nrow+1) /= size(self%entries)+1) then
   ok = .false.
   return
end if
if (any(self%rowpointers(2:) < self%rowpointers(:self%nrow))) then
   ok = .false.
   return
end if
if (size(self%colindices)>0) then
   if (minval(self%colindices)<1 .or. maxval(self%colindices)>self%ncol) then
      ok=.false.
      return
   end if
end if
! spam rows are maintained in strictly increasing column order.
do i=1,self%nrow
   do k=self%rowpointers(i)+1,self%rowpointers(i+1)-1
      if (self%colindices(k) <= self%colindices(k-1)) then
         ok=.false.
         return
      end if
   end do
end do
end function csr_valid
end module spam_types
