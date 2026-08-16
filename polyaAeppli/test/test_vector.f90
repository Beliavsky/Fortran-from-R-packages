program test_vector
   use polya_aeppli
   implicit none
   real(dp) :: x(4),lambda(2),prob(1),d(4),pv(4),pp(4)
   integer :: qv(4),i,rv(1000),fails

   fails = 0
   x = [0.0_dp,1.0_dp,2.0_dp,3.0_dp]
   lambda = [2.0_dp,4.0_dp]
   prob = [0.2_dp]

   call d_polya_aeppli_vec(x,lambda,prob,d)
   do i = 1, 4
      if (abs(d(i)-d_polya_aeppli(x(i),lambda(mod(i-1,2)+1),prob(1))) > &
          1.0e-15_dp) fails=fails+1
   end do

   pp = [0.1_dp,0.25_dp,0.5_dp,0.9_dp]
   call q_polya_aeppli_vec(pp,lambda,prob,qv)
   do i = 1, 4
      if (qv(i) /= q_polya_aeppli(pp(i),lambda(mod(i-1,2)+1),prob(1))) fails=fails+1
   end do

   call p_polya_aeppli_vec(x,lambda,prob,pv)
   do i = 1, 4
      if (abs(pv(i)-p_polya_aeppli(x(i),lambda(mod(i-1,2)+1),prob(1))) > &
          1.0e-15_dp) fails=fails+1
   end do

   call set_polya_aeppli_seed(9)
   call r_polya_aeppli_vec(rv,lambda,prob)
   if (any(rv < 0)) fails=fails+1

   if (fails /= 0) then
      print *, "test_vector: FAIL", fails
      error stop 1
   end if
   print *, "test_vector: PASS"
end program test_vector
