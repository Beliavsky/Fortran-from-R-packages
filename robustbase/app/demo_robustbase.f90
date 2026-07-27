! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of robustbase.
! It may be redistributed and/or modified under GPL version 2 or later.
program demo_robustbase
   use robustbase
   implicit none
   integer,parameter::n=50
   real(dp)::x1(n),x(n,2),y(n),mu,s,se,lf,uf,mc
   integer::i,it
   type(robust_cov_result)::mcd
   type(robust_regression_result)::fit
   do i=1,n
      x1(i)=-2.0_dp+4.0_dp*real(i-1,dp)/real(n-1,dp)
      x(i,1)=x1(i)
      x(i,2)=0.7_dp*x1(i)+sin(real(i,dp))
      y(i)=1.0_dp+2.0_dp*x1(i)+0.05_dp*cos(real(i,dp))
   end do
   x(n-3:n,1)=x(n-3:n,1)+15.0_dp
   x(n-3:n,2)=x(n-3:n,2)-12.0_dp
   y(n-5:n)=y(n-5:n)+20.0_dp
   call huber_location(y,mu,s,it,standard_error=se)
   call adjusted_boxplot_stats(y,lf,uf,mc)
   call cov_mcd(x,mcd,n_starts=100)
   call mm_regression(reshape([spread(1.0_dp,1,n),x1],[n,2]),y,fit,n_starts=100)
   write(*,'(a,f12.6)')'Qn scale: ',qn_scale(y)
   write(*,'(a,f12.6)')'Sn scale: ',sn_scale(y)
   write(*,'(a,f12.6,a,i0)')'Huber location: ',mu,' iterations: ',it
   write(*,'(a,f12.6)')'Medcouple: ',mc
   write(*,'(a,2f12.6)')'Adjusted fences: ',lf,uf
   write(*,'(a,2f12.6)')'MCD center: ',mcd%center
   write(*,'(a,2f12.6)')'MM coefficients: ',fit%coefficients
end program demo_robustbase
