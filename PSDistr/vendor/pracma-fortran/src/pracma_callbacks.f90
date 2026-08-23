! SPDX-License-Identifier: GPL-3.0-or-later
module pracma_callbacks
   use pracma_kinds, only : dp
   implicit none
   private
   public :: scalar_function, scalar_derivative, objective_function, vector_function, vector_field
   public :: complex_scalar_function, complex_objective_function, complex_vector_function
   public :: bivariate_function, trivariate_function, vector_curve, regression_model

   abstract interface
      function scalar_function(x) result(y)
         import dp
         real(dp), intent(in) :: x
         real(dp) :: y
      end function scalar_function


      function scalar_derivative(x) result(y)
         import dp
         real(dp), intent(in) :: x
         real(dp) :: y
      end function scalar_derivative


      function objective_function(x) result(y)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp) :: y
      end function objective_function

      subroutine vector_function(x, y)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: y(:)
      end subroutine vector_function

      subroutine vector_field(t, y, dydt)
         import dp
         real(dp), intent(in) :: t
         real(dp), intent(in) :: y(:)
         real(dp), intent(out) :: dydt(:)
      end subroutine vector_field


      function bivariate_function(x, y) result(z)
         import dp
         real(dp), intent(in) :: x, y
         real(dp) :: z
      end function bivariate_function

      function trivariate_function(x, y, z) result(v)
         import dp
         real(dp), intent(in) :: x, y, z
         real(dp) :: v
      end function trivariate_function

      subroutine vector_curve(t, x)
         import dp
         real(dp), intent(in) :: t
         real(dp), intent(out) :: x(:)
      end subroutine vector_curve


      subroutine regression_model(parameters, x, y)
         import dp
         real(dp), intent(in) :: parameters(:)
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: y(:)
      end subroutine regression_model

      function complex_scalar_function(x) result(y)
         import dp
         complex(dp), intent(in) :: x
         complex(dp) :: y
      end function complex_scalar_function


      function complex_objective_function(x) result(y)
         import dp
         complex(dp), intent(in) :: x(:)
         complex(dp) :: y
      end function complex_objective_function

      subroutine complex_vector_function(x, y)
         import dp
         complex(dp), intent(in) :: x(:)
         complex(dp), intent(out) :: y(:)
      end subroutine complex_vector_function
   end interface
end module pracma_callbacks
