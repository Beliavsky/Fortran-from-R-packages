program test_backend
   use highs
   implicit none
   character(len=:), allocatable :: version
   if (.not. highs_backend_available()) then
      print *, "test_backend: SKIP (run scripts/build_backend first)"
      stop
   end if
   version = highs_backend_version()
   if (len_trim(version) == 0 .or. version == "unavailable") error stop "version unavailable"
   if (highs_infinity() < 1.0e20_dp) error stop "invalid infinity"
   print *, "test_backend: PASS; HiGHS " // trim(version)
end program test_backend
