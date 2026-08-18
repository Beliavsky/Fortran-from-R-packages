program test_family2
   use frbinom
   implicit none
   real(dp), allocatable :: pmf(:)
   integer :: fails,status

   fails = 0
   call frbinom2_pmf_table(10,0.8_dp,0.2_dp,0.1_dp,.false.,pmf,status)
   if (status /= 0) fails=fails+1
   if (abs(pmf(0)-0.7652273598531224_dp) > 3.0e-14_dp) fails=fails+1
   if (abs(pmf(1)-0.1280874279081523_dp) > 3.0e-14_dp) fails=fails+1
   if (abs(pmf(2)-0.06635901093834663_dp) > 3.0e-14_dp) fails=fails+1
   if (abs(pmf(4)-0.009378033386844543_dp) > 3.0e-14_dp) fails=fails+1
   if (abs(sum(pmf)-1.0_dp) > 3.0e-14_dp) fails=fails+1

   call frbinom2_pmf_table(10,0.8_dp,0.2_dp,0.1_dp,.true.,pmf,status)
   if (abs(pmf(0)-0.3758898608917548_dp) > 3.0e-14_dp) fails=fails+1
   if (abs(pmf(1)-0.30343153851279236_dp) > 3.0e-14_dp) fails=fails+1
   if (abs(pmf(2)-0.186384831290287_dp) > 3.0e-14_dp) fails=fails+1
   if (abs(pmf(4)-0.033181353401811164_dp) > 3.0e-14_dp) fails=fails+1
   if (abs(sum(pmf)-1.0_dp) > 3.0e-14_dp) fails=fails+1

   ! Family-II size=1 special cases.
   call frbinom2_pmf_table(1,0.8_dp,0.2_dp,0.1_dp,.false.,pmf,status)
   if (abs(pmf(1)-0.1_dp) > 1.0e-15_dp) fails=fails+1
   call frbinom2_pmf_table(1,0.8_dp,0.2_dp,0.1_dp,.true.,pmf,status)
   if (abs(pmf(1)-0.2_dp) > 1.0e-15_dp) fails=fails+1

   if (fails /= 0) then
      print *, "test_family2: FAIL", fails
      error stop 1
   end if
   print *, "test_family2: PASS"
end program test_family2
