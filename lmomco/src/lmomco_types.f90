module lmomco_types
   use lmomco_kinds, only : dp
   implicit none
   private

   type, public :: lmomco_params
      character(len=8) :: family = ''
      real(dp) :: p(5) = 0.0_dp
      integer :: npar = 0
   end type lmomco_params

   public :: make_params, family_npar

contains

   pure function lower_string(s) result(t)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: t
      integer :: i, k
      t = s
      do i = 1, len(s)
         k = iachar(t(i:i))
         if (k >= iachar('A') .and. k <= iachar('Z')) t(i:i) = achar(k + 32)
      end do
   end function lower_string

   pure integer function family_npar(family) result(n)
      character(len=*), intent(in) :: family
      character(len=:), allocatable :: f
      f = trim(lower_string(family))
      select case (f)
      case ('aep4','gld','kap','smd','gdd')
         if (f == 'gdd') then
            n = 4
         else
            n = 4
         end if
      case ('gev','glo','gno','gov','gpa','ln3','pdq3','pdq4','pe3','st3','tri','wei','gep')
         n = 3
      case ('wak')
         n = 5
      case ('cau','emu','exp','gam','gum','kmu','kur','lap','lmrq','nor','ray','revgum','rice','sla')
         n = 2
      case ('texp')
         n = 3
      case default
         n = 0
      end select
   end function family_npar

   pure function make_params(family, values) result(par)
      character(len=*), intent(in) :: family
      real(dp), intent(in) :: values(:)
      type(lmomco_params) :: par
      integer :: n
      n = min(size(values), 5)
      par%family = lower_string(adjustl(family))
      par%npar = n
      if (n > 0) par%p(1:n) = values(1:n)
   end function make_params

end module lmomco_types
