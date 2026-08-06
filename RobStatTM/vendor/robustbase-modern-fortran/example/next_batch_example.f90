! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of robustbase.
! It may be redistributed and/or modified under GPL version 2 or later.
program next_batch_example
   use robustbase
   implicit none
   integer,parameter::n=30
   real(dp)::x(n,2),y(n)
   type(detmcd_result)::mcd
   type(fast_lts_result)::lts
   integer::i
   do i=1,n
      x(i,1)=1.0_dp
      x(i,2)=-1.0_dp+2.0_dp*real(i-1,dp)/real(n-1,dp)
      y(i)=1.0_dp+2.0_dp*x(i,2)+0.02_dp*sin(real(i,dp))
   end do
   y(n-4:n)=y(n-4:n)+8.0_dp
   call cov_detmcd(x(:,2:2),mcd,alpha=0.75_dp)
   call fast_lts_regression(x,y,lts,alpha=0.5_dp,sampling='best',n_starts=100)
   write(*,'(a,1x,i0)')'detMCD best start:',mcd%best_start
   write(*,'(a,2(1x,f10.5))')'FAST-LTS coefficients:',lts%coefficients
end program next_batch_example
