module ld_interfaces
use ld_kinds, only: dp
implicit none
private
public :: log_target_iface, vector_func_iface, scalar_stat_iface
abstract interface
   function log_target_iface(theta) result(lp)
      import dp
      real(dp), intent(in) :: theta(:)
      real(dp) :: lp
   end function log_target_iface
   subroutine vector_func_iface(theta, y)
      import dp
      real(dp), intent(in) :: theta(:)
      real(dp), intent(out) :: y(:)
   end subroutine vector_func_iface
   function scalar_stat_iface(x, w) result(v)
      import dp
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(in) :: w(:)
      real(dp) :: v
   end function scalar_stat_iface
end interface
end module ld_interfaces
