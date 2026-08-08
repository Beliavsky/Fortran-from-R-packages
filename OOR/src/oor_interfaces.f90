! Upstream OOR license declaration: LGPL (version unspecified).
module oor_interfaces
   use oor_kinds, only : dp
   implicit none
   private
   public :: scalar_objective, vector_objective

   abstract interface
      function scalar_objective(x) result(fx)
         import :: dp
         real(dp), intent(in) :: x
         real(dp) :: fx
      end function scalar_objective

      function vector_objective(x) result(fx)
         import :: dp
         real(dp), intent(in) :: x(:)
         real(dp) :: fx
      end function vector_objective
   end interface
end module oor_interfaces
