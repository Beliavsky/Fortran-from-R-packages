module gsl_qrng
  use iso_c_binding, only: c_ptr, c_null_ptr, c_associated, c_int, c_size_t, c_double, c_char, c_f_pointer
  implicit none
  private

  integer, parameter, public :: qrng_niederreiter_2 = 0
  integer, parameter, public :: qrng_sobol = 1

  type, public :: gsl_qrng_type
    type(c_ptr) :: handle = c_null_ptr
    integer :: dimension = 0
  end type gsl_qrng_type

  public :: qrng_alloc, qrng_free, qrng_clone, qrng_init, qrng_get
  public :: qrng_size, qrng_is_allocated, qrng_get_n, qrng_name

  interface
    function c_qrng_alloc(type_id, dim) bind(C, name='fgsl_qrng_alloc') result(p)
      import :: c_ptr, c_int
      integer(c_int), value :: type_id
      integer(c_int), value :: dim
      type(c_ptr) :: p
    end function
    subroutine c_qrng_free(p) bind(C, name='fgsl_qrng_free')
      import :: c_ptr
      type(c_ptr), value :: p
    end subroutine
    function c_qrng_clone(p) bind(C, name='fgsl_qrng_clone') result(q)
      import :: c_ptr
      type(c_ptr), value :: p
      type(c_ptr) :: q
    end function
    subroutine c_qrng_init(p) bind(C, name='fgsl_qrng_init')
      import :: c_ptr
      type(c_ptr), value :: p
    end subroutine
    function c_qrng_name(p) bind(C, name='fgsl_qrng_name') result(q)
      import :: c_ptr
      type(c_ptr), value :: p
      type(c_ptr) :: q
    end function
    function c_qrng_size(p) bind(C, name='fgsl_qrng_size') result(n)
      import :: c_ptr, c_size_t
      type(c_ptr), value :: p
      integer(c_size_t) :: n
    end function
    function c_qrng_get(p, x) bind(C, name='fgsl_qrng_get') result(status)
      import :: c_ptr, c_double, c_int
      type(c_ptr), value :: p
      real(c_double) :: x(*)
      integer(c_int) :: status
    end function
  end interface

contains

  function qrng_alloc(type_id, dimension) result(q)
    integer, intent(in) :: type_id, dimension
    type(gsl_qrng_type) :: q
    q%handle = c_qrng_alloc(int(type_id,c_int), int(dimension,c_int))
    q%dimension = dimension
  end function qrng_alloc

  subroutine qrng_free(q)
    type(gsl_qrng_type), intent(inout) :: q
    if (c_associated(q%handle)) call c_qrng_free(q%handle)
    q%handle = c_null_ptr
    q%dimension = 0
  end subroutine qrng_free


  logical function qrng_is_allocated(q) result(ok)
    type(gsl_qrng_type), intent(in) :: q
    ok = c_associated(q%handle)
  end function qrng_is_allocated

  function qrng_clone(q) result(r)
    type(gsl_qrng_type), intent(in) :: q
    type(gsl_qrng_type) :: r
    r%handle = c_qrng_clone(q%handle)
    r%dimension = q%dimension
  end function qrng_clone

  subroutine qrng_init(q)
    type(gsl_qrng_type), intent(in) :: q
    call c_qrng_init(q%handle)
  end subroutine qrng_init

  integer(c_size_t) function qrng_size(q) result(n)
    type(gsl_qrng_type), intent(in) :: q
    n = c_qrng_size(q%handle)
  end function qrng_size

  integer function qrng_get(q, x) result(status)
    type(gsl_qrng_type), intent(in) :: q
    real(c_double), intent(out) :: x(:)
    if (size(x) /= q%dimension) then
      status = -1
      return
    end if
    status = int(c_qrng_get(q%handle, x), kind(status))
  end function qrng_get

  integer function qrng_get_n(q, x) result(status)
    type(gsl_qrng_type), intent(in) :: q
    real(c_double), intent(out) :: x(:,:)
    integer :: i, st
    if (size(x,1) /= q%dimension) then
      status = -1
      return
    end if
    do i = 1, size(x,2)
      st = int(c_qrng_get(q%handle, x(:,i)), kind(st))
      if (st /= 0) then
        status = st
        return
      end if
    end do
    status = 0
  end function qrng_get_n

  function qrng_name(q) result(name)
    type(gsl_qrng_type), intent(in) :: q
    character(len=:), allocatable :: name
    type(c_ptr) :: p
    character(kind=c_char), pointer :: chars(:)
    integer :: i, n
    p = c_qrng_name(q%handle)
    if (.not. c_associated(p)) then
      name = ''
      return
    end if
    call c_f_pointer(p, chars, [256])
    n = 0
    do i = 1, 256
      if (iachar(chars(i)) == 0) exit
      n = n + 1
    end do
    allocate(character(len=n) :: name)
    do i = 1, n
      name(i:i) = chars(i)
    end do
  end function qrng_name

end module gsl_qrng
