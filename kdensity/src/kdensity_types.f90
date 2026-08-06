module kdensity_types
  use kdensity_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: kd_ok=0, kd_invalid_input=1, kd_unsupported=2
  integer, parameter, public :: kd_integration_failed=3, kd_optimization_failed=4

  abstract interface
    pure function start_density_interface(x, parameters) result(value)
      import dp
      real(dp), intent(in) :: x
      real(dp), intent(in) :: parameters(:)
      real(dp) :: value
    end function start_density_interface
    subroutine start_estimator_interface(data, parameters, status)
      import dp
      real(dp), intent(in) :: data(:)
      real(dp), allocatable, intent(out) :: parameters(:)
      integer, intent(out) :: status
    end subroutine start_estimator_interface
    pure function kernel_interface(y, x, h) result(value)
      import dp
      real(dp), intent(in) :: y, x, h
      real(dp) :: value
    end function kernel_interface
  end interface

  type, public :: kd_start
    character(len=32) :: name='uniform'
    real(dp) :: support(2)=[-huge(1.0_dp),huge(1.0_dp)]
    procedure(start_density_interface), pointer, nopass :: density => null()
    procedure(start_estimator_interface), pointer, nopass :: estimator => null()
  end type kd_start

  type, public :: kd_kernel
    character(len=32) :: name='gaussian'
    real(dp) :: support(2)=[-huge(1.0_dp),huge(1.0_dp)]
    logical :: has_sd=.true.
    real(dp) :: sd=1.0_dp
    procedure(kernel_interface), pointer, nopass :: evaluate => null()
  end type kd_kernel

  type, public :: kdensity_options
    character(len=32) :: kernel='auto'
    character(len=32) :: start='uniform'
    character(len=32) :: bandwidth='auto'
    real(dp) :: bw=-1.0_dp
    real(dp) :: adjust=1.0_dp
    real(dp) :: support(2)=[-huge(1.0_dp),huge(1.0_dp)]
    logical :: support_supplied=.false.
    logical :: normalized=.true.
    real(dp) :: integration_tolerance=1.0e-7_dp
  end type kdensity_options

  type, public :: kdensity_fit
    real(dp), allocatable :: data(:)
    real(dp), allocatable :: parameters(:)
    type(kd_start) :: start
    type(kd_kernel) :: kernel
    real(dp) :: support(2)=[-huge(1.0_dp),huge(1.0_dp)]
    real(dp) :: bw=0.0_dp
    real(dp) :: adjust=1.0_dp
    real(dp) :: h=0.0_dp
    real(dp) :: normalization=1.0_dp
    real(dp) :: parametric_loglik=0.0_dp
    logical :: normalized=.true.
    integer :: status=kd_ok
    character(len=256) :: message='ok'
  contains
    procedure :: pdf => kdensity_pdf_scalar
    procedure :: pdf_vector => kdensity_pdf_vector
    procedure :: start_pdf => kdensity_start_pdf
  end type kdensity_fit

contains

  pure function kdensity_start_pdf(self,y) result(value)
    class(kdensity_fit), intent(in) :: self
    real(dp), intent(in) :: y
    real(dp) :: value
    if (associated(self%start%density)) then
      value=self%start%density(y,self%parameters)
    else
      value=0.0_dp
    end if
  end function kdensity_start_pdf

  pure function kdensity_pdf_scalar(self,y) result(value)
    class(kdensity_fit), intent(in) :: self
    real(dp), intent(in) :: y
    real(dp) :: value, sy, denom
    integer :: i
    if (self%status /= kd_ok .or. .not. associated(self%start%density)) then
      value=0.0_dp; return
    end if
    if (y < self%support(1) .or. y > self%support(2)) then
      value=0.0_dp; return
    end if
    if (self%h >= 0.5_dp*huge(1.0_dp)) then
      value=self%start%density(y,self%parameters); return
    end if
    sy=self%start%density(y,self%parameters)
    if (sy <= 0.0_dp .or. self%h <= 0.0_dp) then
      value=0.0_dp; return
    end if
    value=0.0_dp
    do i=1,size(self%data)
      denom=self%start%density(self%data(i),self%parameters)
      if (denom > tiny(1.0_dp)) then
        value=value+self%kernel%evaluate(y,self%data(i),self%h)*sy/denom
      end if
    end do
    value=value/(real(size(self%data),dp)*self%h*self%normalization)
  end function kdensity_pdf_scalar

  pure function kdensity_pdf_vector(self,y) result(value)
    class(kdensity_fit), intent(in) :: self
    real(dp), intent(in) :: y(:)
    real(dp) :: value(size(y))
    integer :: i
    do i=1,size(y)
      value(i)=self%pdf(y(i))
    end do
  end function kdensity_pdf_vector

end module kdensity_types
