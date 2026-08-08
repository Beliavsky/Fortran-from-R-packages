! SPDX-License-Identifier: Apache-2.0
program generic_quadratic
  use psqn
  implicit none

  type(psqn_element_spec), allocatable :: specs(:)
  type(psqn_info) :: info
  real(dp) :: x(3)

  allocate(specs(2))
  specs(1)%idx = [1, 2]
  specs(2)%idx = [2, 3]
  x = 0.0_dp

  call psqn_optimize_generic(x, specs, element, info)
  print '(a,3f12.6)', 'x = ', x
  print '(a,es14.6)', 'f = ', info%value
  print '(a,i0)', 'status = ', info%info

contains

  subroutine element(i, z, f, g, comp_grad)
    integer, intent(in) :: i
    real(dp), intent(in) :: z(:)
    real(dp), intent(out) :: f
    real(dp), intent(out) :: g(:)
    logical, intent(in) :: comp_grad

    select case (i)
    case (1)
      f = 0.5_dp * ((z(1) - 1.0_dp)**2 + (z(2) + 2.0_dp)**2)
      if (comp_grad) g = [z(1) - 1.0_dp, z(2) + 2.0_dp]
    case (2)
      f = 0.5_dp * ((z(1) + 2.0_dp)**2 + (z(2) - 3.0_dp)**2)
      if (comp_grad) g = [z(1) + 2.0_dp, z(2) - 3.0_dp]
    case default
      f = huge(1.0_dp)
      if (comp_grad) g = 0.0_dp
    end select
  end subroutine element

end program generic_quadratic
