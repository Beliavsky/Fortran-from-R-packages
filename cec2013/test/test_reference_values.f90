program test_reference_values
   use cec2013
   implicit none
   type(cec2013_context) :: ctx
   real(dp) :: x(10), got, scale
   real(dp), parameter :: expected(28) = [ &
      2.04381181288061198e4_dp, 1.18806586518990350e9_dp, &
      8.71893064996279936e17_dp, 1.59475167307608098e8_dp, &
      7.51043450590324210e4_dp, 2.16025043548001213e3_dp, &
      2.53498324137073010e6_dp, -6.78341479026762613e2_dp, &
      -5.80448134144029154e2_dp, 2.31850026728556850e3_dp, &
      -1.40712914648453875e1_dp, 5.09780854500377245e-1_dp, &
      7.50413064575145086e1_dp, 3.17177371624759553e3_dp, &
      3.95537577674006798e3_dp, 2.24663962281437023e2_dp, &
      6.05566011385872230e2_dp, 7.08496912082925974e2_dp, &
      2.74451324858295033e5_dp, 6.05000000000000114e2_dp, &
      2.06418123360306436e3_dp, 4.88015960333976545e3_dp, &
      4.59297591646461660e3_dp, 1.65900417848366328e3_dp, &
      1.41802095617497753e3_dp, 1.44894802159504288e4_dp, &
      2.70257214167002758e3_dp, 3.38727853830248250e3_dp ]
   integer :: i, status

   x = [-3.70_dp, 5.01_dp, -6.32_dp, 7.63_dp, -8.94_dp, &
         10.25_dp, -11.56_dp, 12.87_dp, -14.18_dp, 15.49_dp]
   call ctx%init(10, 'data', status)
   if (status /= CEC2013_OK) error stop 'context initialization failed'

   do i = 1, 28
      got = cec2013_evaluate(ctx, i, x, status)
      if (status /= CEC2013_OK) error stop 'evaluation failed'
      scale = max(1.0_dp, abs(expected(i)))
      if (abs(got-expected(i)) > 5.0e-12_dp*scale) then
         write(*,'(a,i0,2(1x,es24.16))') 'reference mismatch for problem', i, got, expected(i)
         error stop 1
      end if
   end do
   print *, 'PASS test_reference_values'
end program test_reference_values
