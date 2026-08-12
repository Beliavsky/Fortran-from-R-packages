program basic_qap
   use qap
   implicit none
   real(dp) :: A(4,4), B(4,4)
   integer :: p(4)

   A = reshape([ &
      0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, &
      1.0_dp, 0.0_dp, 4.0_dp, 2.0_dp, &
      2.0_dp, 4.0_dp, 0.0_dp, 1.0_dp, &
      3.0_dp, 2.0_dp, 1.0_dp, 0.0_dp], [4,4])
   B = reshape([ &
      0.0_dp, 3.0_dp, 1.0_dp, 2.0_dp, &
      3.0_dp, 0.0_dp, 2.0_dp, 4.0_dp, &
      1.0_dp, 2.0_dp, 0.0_dp, 3.0_dp, &
      2.0_dp, 4.0_dp, 3.0_dp, 0.0_dp], [4,4])
   p = [2,4,1,3]
   write(*,'(a,f10.2)') 'QAP objective = ', qap_obj(A,B,p)
end program basic_qap
