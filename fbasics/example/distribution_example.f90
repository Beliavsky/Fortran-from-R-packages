! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from the GPL-2-or-later fBasics package.
program distribution_example
  use fbasics
  implicit none
  real(dp) :: p,x
  integer :: i
  write(*,'(a)')' p        Normal q       Student(6) q      NIG density'
  do i=1,9
    p=real(i,dp)/10.0_dp
    x=normal_quantile(p)
    write(*,'(f5.2,3(2x,es14.6))')p,x,student_quantile(p,6.0_dp),dnig(x,2.0_dp,0.3_dp,1.0_dp,0.0_dp)
  end do
end program distribution_example
