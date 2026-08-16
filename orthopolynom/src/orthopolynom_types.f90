module orthopolynom_types
  use polynom, only : dp, polynomial_t
  implicit none
  private

  type, public :: recurrence_t
    real(dp), allocatable :: c(:)
    real(dp), allocatable :: d(:)
    real(dp), allocatable :: e(:)
    real(dp), allocatable :: f(:)
  contains
    procedure :: size => recurrence_size
  end type recurrence_t

  type, public :: monic_recurrence_t
    real(dp), allocatable :: a(:)
    real(dp), allocatable :: b(:)
  contains
    procedure :: size => monic_recurrence_size
  end type monic_recurrence_t

  type, public :: real_vector_t
    real(dp), allocatable :: value(:)
  end type real_vector_t

  type, public :: real_vector_list_t
    type(real_vector_t), allocatable :: item(:)
  contains
    procedure :: size => real_vector_list_size
  end type real_vector_list_t

  type, public :: real_matrix_t
    real(dp), allocatable :: value(:,:)
  end type real_matrix_t

  type, public :: real_matrix_list_t
    type(real_matrix_t), allocatable :: item(:)
  contains
    procedure :: size => real_matrix_list_size
  end type real_matrix_list_t

  type, public :: polynomial_function_t
    type(polynomial_t) :: polynomial
  contains
    procedure :: evaluate_scalar => polynomial_function_evaluate_scalar
    procedure :: evaluate_vector => polynomial_function_evaluate_vector
    generic :: evaluate => evaluate_scalar, evaluate_vector
  end type polynomial_function_t

  type, public :: polynomial_function_list_t
    type(polynomial_function_t), allocatable :: item(:)
  contains
    procedure :: size => polynomial_function_list_size
  end type polynomial_function_list_t

contains

  integer function recurrence_size(self)
    class(recurrence_t), intent(in) :: self
    if (allocated(self%c)) then
      recurrence_size = size(self%c)
    else
      recurrence_size = 0
    end if
  end function recurrence_size

  integer function monic_recurrence_size(self)
    class(monic_recurrence_t), intent(in) :: self
    if (allocated(self%a)) then
      monic_recurrence_size = size(self%a)
    else
      monic_recurrence_size = 0
    end if
  end function monic_recurrence_size

  integer function real_vector_list_size(self)
    class(real_vector_list_t), intent(in) :: self
    if (allocated(self%item)) then
      real_vector_list_size = size(self%item)
    else
      real_vector_list_size = 0
    end if
  end function real_vector_list_size

  integer function real_matrix_list_size(self)
    class(real_matrix_list_t), intent(in) :: self
    if (allocated(self%item)) then
      real_matrix_list_size = size(self%item)
    else
      real_matrix_list_size = 0
    end if
  end function real_matrix_list_size

  integer function polynomial_function_list_size(self)
    class(polynomial_function_list_t), intent(in) :: self
    if (allocated(self%item)) then
      polynomial_function_list_size = size(self%item)
    else
      polynomial_function_list_size = 0
    end if
  end function polynomial_function_list_size

  pure real(dp) function polynomial_function_evaluate_scalar(self, x)
    class(polynomial_function_t), intent(in) :: self
    real(dp), intent(in) :: x
    polynomial_function_evaluate_scalar = self%polynomial%evaluate(x)
  end function polynomial_function_evaluate_scalar

  pure function polynomial_function_evaluate_vector(self, x) result(y)
    class(polynomial_function_t), intent(in) :: self
    real(dp), intent(in) :: x(:)
    real(dp) :: y(size(x))
    y = self%polynomial%evaluate(x)
  end function polynomial_function_evaluate_vector

end module orthopolynom_types
