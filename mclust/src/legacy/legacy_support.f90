! Compatibility replacements for R's historical d1mach/i1mach/intpr helpers.
double precision function d1mach(i)
  implicit none
  integer, intent(in) :: i
  select case(i)
  case(1); d1mach = tiny(1.0d0)
  case(2); d1mach = huge(1.0d0)
  case(3); d1mach = 0.5d0*epsilon(1.0d0)
  case(4); d1mach = epsilon(1.0d0)
  case(5); d1mach = log10(2.0d0)
  case default; d1mach = 0.0d0
  end select
end function d1mach

integer function i1mach(i)
  implicit none
  integer, intent(in) :: i
  select case(i)
  case(9); i1mach = huge(1)
  case default; i1mach = 0
  end select
end function i1mach

subroutine intpr(label, nchar, x, n)
  implicit none
  character(len=*), intent(in) :: label
  integer, intent(in) :: nchar, n
  integer, intent(in) :: x(*)
  return
end subroutine intpr
