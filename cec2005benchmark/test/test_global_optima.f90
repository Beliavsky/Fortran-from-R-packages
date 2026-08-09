program test_global_optima
   use cec2005benchmark
   implicit none
   integer, parameter :: dims(4) = [2,10,30,50]
   integer :: fid, u, i, n, ios, k
   real(dp) :: opt(25,100), bias(25), x(50), got, err, worst

   open(newunit=u, file='data/global_optima.txt', status='old', action='read')
   do i = 1, 25
      read(u,*) opt(i,:)
   end do
   close(u)
   open(newunit=u, file='data/fbias_data.txt', status='old', action='read')
   read(u,*) bias
   close(u)

   worst = 0.0_dp
   do fid = 1, 25
      do k = 1, 4
         n = dims(k)
         x(1:n) = opt(fid,1:n)
         if (fid == 5) then
            x(1:(n+3)/4) = -100.0_dp
            x(max((3*n)/4,1):n) = 100.0_dp
         end if
         if (fid == 8) then
            do i = 1, n/2
               x(2*i-1) = -32.0_dp
            end do
         end if
         if (fid == 20) then
            do i = 1, n/2
               x(2*i) = 5.0_dp
            end do
         end if
         got = cec2005_eval(fid,x(1:n),'data',.false.,ios)
         if (ios /= 0) error stop 'failed to initialize optimum case'
         err = abs(got-bias(fid))/max(1.0_dp,abs(bias(fid)))
         worst = max(worst,err)
         if (err > 1.0e-8_dp) error stop 'global-optimum mismatch'
      end do
   end do
   print '(a,es12.4)', 'PASS test_global_optima; worst relative error = ', worst
end program test_global_optima
