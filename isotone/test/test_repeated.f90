program test_repeated
   use isotone
   implicit none
   real(dp) :: y(4,3), z(4), w(4,3)
   type(gpava_repeated_result) :: r
   y=reshape([3.0_dp,2.0_dp,5.0_dp,4.0_dp, &
              1.0_dp,4.0_dp,1.0_dp,6.0_dp, &
              2.0_dp,3.0_dp,2.0_dp,5.0_dp],[4,3])
   z=[1.0_dp,1.0_dp,2.0_dp,3.0_dp];w=1.0_dp
   call gpava_fit_repeated(y,r,z=z,weights=w,ties=GPAVA_SECONDARY)
   if(r%status/=0) error stop 'repeated gpava failed'
   if(abs(r%x(1)-r%x(2))>1.0e-13_dp) error stop 'secondary tie mismatch'
   if(r%x(2)>r%x(3)+1.0e-12_dp .or. r%x(3)>r%x(4)+1.0e-12_dp) then
      error stop 'repeated fit not monotone'
   end if
   print *, 'test_repeated: PASS'
end program
