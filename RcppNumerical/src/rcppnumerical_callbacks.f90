module rcppnumerical_callbacks
  use rcppnumerical_kinds, only : dp
  implicit none
  private

  abstract interface
    function scalar_function_interface(x, user_data) result(f)
      import dp
      real(dp), intent(in) :: x
      class(*), intent(inout), optional :: user_data
      real(dp) :: f
    end function scalar_function_interface

    function multivariate_function_interface(x, user_data) result(f)
      import dp
      real(dp), intent(in) :: x(:)
      class(*), intent(inout), optional :: user_data
      real(dp) :: f
    end function multivariate_function_interface

    subroutine objective_gradient_interface(x, f, g, user_data)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: f
      real(dp), intent(out) :: g(:)
      class(*), intent(inout), optional :: user_data
    end subroutine objective_gradient_interface
  end interface

  public :: scalar_function_interface
  public :: multivariate_function_interface
  public :: objective_gradient_interface
end module rcppnumerical_callbacks
