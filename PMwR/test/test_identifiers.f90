program test_identifiers
   use pmwr, only : dp, valid_isin, valid_sedol, quote32_from_string, quote32_components, format_quote32
   implicit none
   real(dp) :: q
   integer :: h, t, f, status
   character(len=32) :: text

   if (.not. valid_isin('US0378331005')) error stop "valid ISIN rejected"
   if (valid_isin('US0378331006')) error stop "invalid ISIN accepted"
   if (.not. valid_sedol('0263494')) error stop "valid SEDOL rejected"
   if (valid_sedol('0263495')) error stop "invalid SEDOL accepted"

   call quote32_from_string('99-162', q, status)
   if (status /= 0) error stop "quote parse status"
   call assert_close(q, 99.5078125_dp, 1.0e-12_dp, "quote value")
   call quote32_components(q, h, t, f)
   if (h /= 99 .or. t /= 16 .or. f /= 1) error stop "quote components"
   text = format_quote32(q)
   if (trim(text) /= '99-162') error stop "quote format"

   print *, "test_identifiers: PASS"
contains
   subroutine assert_close(a,b,tol,label)
      real(dp), intent(in) :: a,b,tol
      character(len=*), intent(in) :: label
      if (abs(a-b) > tol) then
         print *, trim(label), a, b
         error stop 1
      end if
   end subroutine assert_close
end program test_identifiers
