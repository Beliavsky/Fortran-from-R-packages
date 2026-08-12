module qap_io
   use qap_kinds, only : dp
   use qap_types, only : qap_problem_t
   implicit none
   private
   public :: read_qaplib

contains

   subroutine read_qaplib(filename, problem, status, message)
      character(len=*), intent(in) :: filename
      type(qap_problem_t), intent(out) :: problem
      integer, intent(out), optional :: status
      character(len=:), allocatable, intent(out), optional :: message

      integer :: u, ios, n, i, j
      real(dp), allocatable :: flat(:)
      character(len=:), allocatable :: sln
      logical :: exists
      integer :: nsol

      if (present(status)) status = 0
      if (present(message)) message = ""

      open(newunit=u, file=filename, status="old", action="read", iostat=ios)
      if (ios /= 0) then
         call fail("read_qaplib: cannot open problem file")
         return
      end if
      read(u, *, iostat=ios) n
      if (ios /= 0 .or. n <= 0) then
         close(u)
         call fail("read_qaplib: invalid problem dimension")
         return
      end if
      allocate(flat(2*n*n))
      read(u, *, iostat=ios) flat
      close(u)
      if (ios /= 0) then
         call fail("read_qaplib: incomplete problem data")
         return
      end if

      allocate(problem%A(n,n), problem%B(n,n), problem%solution(0))
      do i = 1, n
         do j = 1, n
            problem%A(i,j) = flat((i-1)*n+j)
            problem%B(i,j) = flat(n*n+(i-1)*n+j)
         end do
      end do

      sln = solution_filename(filename)
      inquire(file=sln, exist=exists)
      if (.not. exists) return

      open(newunit=u, file=sln, status="old", action="read", iostat=ios)
      if (ios /= 0) return
      read(u, *, iostat=ios) nsol, problem%opt
      if (ios /= 0 .or. nsol /= n) then
         close(u)
         problem%opt = huge(1.0_dp)
         return
      end if
      deallocate(problem%solution)
      allocate(problem%solution(n))
      read(u, *, iostat=ios) problem%solution
      close(u)
      if (ios /= 0) then
         deallocate(problem%solution)
         allocate(problem%solution(0))
         problem%opt = huge(1.0_dp)
         return
      end if
      problem%has_solution = .true.

   contains

      subroutine fail(text)
         character(len=*), intent(in) :: text
         if (present(status)) then
            status = 1
            if (present(message)) message = text
         else
            error stop text
         end if
      end subroutine fail

   end subroutine read_qaplib

   function solution_filename(filename) result(sln)
      character(len=*), intent(in) :: filename
      character(len=:), allocatable :: sln
      integer :: p

      p = index(filename, ".dat", back=.true.)
      if (p > 0 .and. p + 3 == len_trim(filename)) then
         sln = filename(:p-1) // ".sln"
      else
         sln = trim(filename) // ".sln"
      end if
   end function solution_filename

end module qap_io
