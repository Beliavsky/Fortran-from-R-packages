! SPDX-License-Identifier: GPL-2.0-only
program mi_example
   use mitools, only : dp, mi_combine, mi_result, mi_summary, mitools_success
   implicit none
   real(dp) :: estimates_in(5)
   real(dp) :: variances(5)
   real(dp), allocatable :: estimates(:)
   real(dp), allocatable :: lower(:)
   real(dp), allocatable :: missinfo(:)
   real(dp), allocatable :: se(:)
   real(dp), allocatable :: upper(:)
   type(mi_result) :: result
   integer :: status

   estimates_in = [0.92_dp, 1.05_dp, 0.98_dp, 1.10_dp, 0.95_dp]
   variances = [0.040_dp, 0.042_dp, 0.041_dp, 0.039_dp, 0.040_dp]

   call mi_combine(estimates_in, variances, result, status)
   if (status /= mitools_success) error stop "MI combination failed"
   call mi_summary(result, 0.05_dp, estimates, se, lower, upper, missinfo, status)
   if (status /= mitools_success) error stop "MI summary failed"

   print '(a,f8.4)', "combined estimate: ", estimates(1)
   print '(a,f8.4)', "standard error:    ", se(1)
   print '(a,2f9.4)', "95% interval:      ", lower(1), upper(1)
   print '(a,f8.4)', "missing information:", missinfo(1)
end program mi_example
