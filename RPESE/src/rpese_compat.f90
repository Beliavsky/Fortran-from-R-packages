! SPDX-License-Identifier: GPL-3.0-or-later
module rpese_compat
  use rpese_types, only : se_if_iid, se_if_cor, se_if_cor_adapt, se_if_cor_pw, &
    se_boot_iid, se_boot_cor, fit_exponential, fit_gamma, frequency_all, &
    frequency_decimate, frequency_truncate
  implicit none
  private
  public :: se_method_from_name, fitting_method_from_name, frequency_mode_from_name
contains
  integer function se_method_from_name(name) result(method)
    character(len=*), intent(in) :: name
    character(len=:), allocatable :: key
    key = lower_case(trim(adjustl(name)))
    select case (key)
    case ('ifiid')
      method = se_if_iid
    case ('ifcor')
      method = se_if_cor
    case ('ifcoradapt')
      method = se_if_cor_adapt
    case ('ifcorpw')
      method = se_if_cor_pw
    case ('bootiid')
      method = se_boot_iid
    case ('bootcor')
      method = se_boot_cor
    case default
      method = 0
    end select
  end function se_method_from_name

  integer function fitting_method_from_name(name) result(method)
    character(len=*), intent(in) :: name
    character(len=:), allocatable :: key
    key = lower_case(trim(adjustl(name)))
    select case (key)
    case ('exponential')
      method = fit_exponential
    case ('gamma')
      method = fit_gamma
    case default
      method = 0
    end select
  end function fitting_method_from_name

  integer function frequency_mode_from_name(name) result(mode)
    character(len=*), intent(in) :: name
    character(len=:), allocatable :: key
    key = lower_case(trim(adjustl(name)))
    select case (key)
    case ('all')
      mode = frequency_all
    case ('decimate')
      mode = frequency_decimate
    case ('truncate')
      mode = frequency_truncate
    case default
      mode = 0
    end select
  end function frequency_mode_from_name

  pure function lower_case(text) result(lower)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: lower
    integer :: i, code
    lower = text
    do i = 1, len(text)
      code = iachar(text(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code + 32)
    end do
  end function lower_case
end module rpese_compat
