module gsl_rng
  use iso_c_binding, only: c_ptr, c_null_ptr, c_associated, c_int, c_long, c_double, c_char, c_f_pointer
  implicit none
  private

  integer, parameter, public :: rng_mt19937 = 0
  integer, parameter, public :: rng_ranlxs0 = 1
  integer, parameter, public :: rng_ranlxs1 = 2
  integer, parameter, public :: rng_ranlxs2 = 3
  integer, parameter, public :: rng_ranlxd1 = 4
  integer, parameter, public :: rng_ranlxd2 = 5
  integer, parameter, public :: rng_ranlux = 6
  integer, parameter, public :: rng_ranlux389 = 7
  integer, parameter, public :: rng_cmrg = 8
  integer, parameter, public :: rng_mrg = 9
  integer, parameter, public :: rng_taus = 10
  integer, parameter, public :: rng_taus2 = 11
  integer, parameter, public :: rng_gfsr4 = 12
  integer, parameter, public :: rng_minstd = 13

  type, public :: gsl_rng_type
    type(c_ptr) :: handle = c_null_ptr
  end type gsl_rng_type

  public :: rng_alloc, rng_free, rng_set, rng_clone, rng_get
  public :: rng_uniform, rng_uniform_pos, rng_uniform_int
  public :: rng_min, rng_max, rng_is_allocated, rng_name
  public :: rng_get_array, rng_uniform_array, rng_uniform_pos_array, rng_uniform_int_array

  interface
    function c_rng_alloc(type_id) bind(C, name='fgsl_rng_alloc') result(p)
      import :: c_ptr, c_int
      integer(c_int), value :: type_id
      type(c_ptr) :: p
    end function
    subroutine c_rng_free(p) bind(C, name='fgsl_rng_free')
      import :: c_ptr
      type(c_ptr), value :: p
    end subroutine
    subroutine c_rng_set(p, seed) bind(C, name='fgsl_rng_set')
      import :: c_ptr, c_long
      type(c_ptr), value :: p
      integer(c_long), value :: seed
    end subroutine
    function c_rng_clone(p) bind(C, name='fgsl_rng_clone') result(q)
      import :: c_ptr
      type(c_ptr), value :: p
      type(c_ptr) :: q
    end function
    function c_rng_name(p) bind(C, name='fgsl_rng_name') result(q)
      import :: c_ptr
      type(c_ptr), value :: p
      type(c_ptr) :: q
    end function
    function c_rng_get(p) bind(C, name='fgsl_rng_get') result(x)
      import :: c_ptr, c_long
      type(c_ptr), value :: p
      integer(c_long) :: x
    end function
    function c_rng_uniform(p) bind(C, name='fgsl_rng_uniform') result(x)
      import :: c_ptr, c_double
      type(c_ptr), value :: p
      real(c_double) :: x
    end function
    function c_rng_uniform_pos(p) bind(C, name='fgsl_rng_uniform_pos') result(x)
      import :: c_ptr, c_double
      type(c_ptr), value :: p
      real(c_double) :: x
    end function
    function c_rng_uniform_int(p, n) bind(C, name='fgsl_rng_uniform_int') result(x)
      import :: c_ptr, c_long
      type(c_ptr), value :: p
      integer(c_long), value :: n
      integer(c_long) :: x
    end function
    function c_rng_min(p) bind(C, name='fgsl_rng_min') result(x)
      import :: c_ptr, c_long
      type(c_ptr), value :: p
      integer(c_long) :: x
    end function
    function c_rng_max(p) bind(C, name='fgsl_rng_max') result(x)
      import :: c_ptr, c_long
      type(c_ptr), value :: p
      integer(c_long) :: x
    end function
  end interface

contains

  function rng_alloc(type_id) result(r)
    integer, intent(in), optional :: type_id
    type(gsl_rng_type) :: r
    integer(c_int) :: k
    k = rng_mt19937
    if (present(type_id)) k = int(type_id, c_int)
    r%handle = c_rng_alloc(k)
  end function rng_alloc

  subroutine rng_free(r)
    type(gsl_rng_type), intent(inout) :: r
    if (c_associated(r%handle)) call c_rng_free(r%handle)
    r%handle = c_null_ptr
  end subroutine rng_free


  logical function rng_is_allocated(r) result(ok)
    type(gsl_rng_type), intent(in) :: r
    ok = c_associated(r%handle)
  end function rng_is_allocated

  subroutine rng_set(r, seed)
    type(gsl_rng_type), intent(in) :: r
    integer(c_long), intent(in) :: seed
    call c_rng_set(r%handle, seed)
  end subroutine rng_set

  function rng_clone(r) result(q)
    type(gsl_rng_type), intent(in) :: r
    type(gsl_rng_type) :: q
    q%handle = c_rng_clone(r%handle)
  end function rng_clone

  integer(c_long) function rng_get(r) result(x)
    type(gsl_rng_type), intent(in) :: r
    x = c_rng_get(r%handle)
  end function rng_get

  real(c_double) function rng_uniform(r) result(x)
    type(gsl_rng_type), intent(in) :: r
    x = c_rng_uniform(r%handle)
  end function rng_uniform

  real(c_double) function rng_uniform_pos(r) result(x)
    type(gsl_rng_type), intent(in) :: r
    x = c_rng_uniform_pos(r%handle)
  end function rng_uniform_pos

  integer(c_long) function rng_uniform_int(r, n) result(x)
    type(gsl_rng_type), intent(in) :: r
    integer(c_long), intent(in) :: n
    x = c_rng_uniform_int(r%handle, n)
  end function rng_uniform_int

  integer(c_long) function rng_min(r) result(x)
    type(gsl_rng_type), intent(in) :: r
    x = c_rng_min(r%handle)
  end function rng_min

  integer(c_long) function rng_max(r) result(x)
    type(gsl_rng_type), intent(in) :: r
    x = c_rng_max(r%handle)
  end function rng_max

  function rng_name(r) result(name)
    type(gsl_rng_type), intent(in) :: r
    character(len=:), allocatable :: name
    type(c_ptr) :: p
    character(kind=c_char), pointer :: chars(:)
    integer :: i, n
    p = c_rng_name(r%handle)
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
  end function rng_name

  subroutine rng_get_array(r, x)
    type(gsl_rng_type), intent(in) :: r
    integer(c_long), intent(out) :: x(:)
    integer :: i
    do i = 1, size(x)
      x(i) = c_rng_get(r%handle)
    end do
  end subroutine rng_get_array

  subroutine rng_uniform_array(r, x)
    type(gsl_rng_type), intent(in) :: r
    real(c_double), intent(out) :: x(:)
    integer :: i
    do i = 1, size(x)
      x(i) = c_rng_uniform(r%handle)
    end do
  end subroutine rng_uniform_array

  subroutine rng_uniform_pos_array(r, x)
    type(gsl_rng_type), intent(in) :: r
    real(c_double), intent(out) :: x(:)
    integer :: i
    do i = 1, size(x)
      x(i) = c_rng_uniform_pos(r%handle)
    end do
  end subroutine rng_uniform_pos_array

  subroutine rng_uniform_int_array(r, n, x)
    type(gsl_rng_type), intent(in) :: r
    integer(c_long), intent(in) :: n
    integer(c_long), intent(out) :: x(:)
    integer :: i
    do i = 1, size(x)
      x(i) = c_rng_uniform_int(r%handle, n)
    end do
  end subroutine rng_uniform_int_array

end module gsl_rng
