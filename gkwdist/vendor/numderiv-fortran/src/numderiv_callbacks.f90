! SPDX-License-Identifier: GPL-2.0-or-later
module numderiv_callbacks
   use numderiv_kinds, only : dp
   implicit none
   private

   public :: scalar_real_function, vector_real_function
   public :: scalar_complex_function, vector_complex_function

   abstract interface
      function scalar_real_function(x) result(value)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp) :: value
      end function scalar_real_function

      function vector_real_function(x) result(value)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp), allocatable :: value(:)
      end function vector_real_function

      function scalar_complex_function(x) result(value)
         import dp
         complex(dp), intent(in) :: x(:)
         complex(dp) :: value
      end function scalar_complex_function

      function vector_complex_function(x) result(value)
         import dp
         complex(dp), intent(in) :: x(:)
         complex(dp), allocatable :: value(:)
      end function vector_complex_function
   end interface

end module numderiv_callbacks
