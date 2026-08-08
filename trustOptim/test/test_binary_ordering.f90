program test_binary_ordering
   use trustoptim
   use trustoptim_binary
   implicit none
   type(binary_data) :: dat
   type(binary_priors) :: pri
   real(dp) :: pcol(8), prow(8), gcol(8), grow(8)
   real(dp) :: fcol, frow
   integer :: n, k, i, j

   n = 3
   k = 2
   allocate(dat%y(n), dat%x(k,n), pri%inv_sigma(k,k), pri%inv_omega(k,k))
   dat%y = [2.0_dp, 7.0_dp, 4.0_dp]
   dat%x = reshape([1.0_dp,0.5_dp, -0.2_dp,1.3_dp, 0.7_dp,-0.9_dp],[k,n])
   dat%trials = 9
   pri%inv_sigma = 0.0_dp
   pri%inv_omega = 0.0_dp
   pri%inv_sigma(1,1) = 1.0_dp
   pri%inv_sigma(2,2) = 2.0_dp
   pri%inv_omega(1,1) = 0.3_dp
   pri%inv_omega(2,2) = 0.4_dp

   pcol = [0.1_dp,0.2_dp, 0.3_dp,0.4_dp, 0.5_dp,0.6_dp, -0.2_dp,0.7_dp]
   do j = 1, n
      do i = 1, k
         prow((i-1)*n+j) = pcol((j-1)*k+i)
      end do
   end do
   prow(n*k+1:) = pcol(n*k+1:)

   fcol = binary_value(pcol, dat, pri, .false.)
   frow = binary_value(prow, dat, pri, .true.)
   if (abs(fcol-frow) > 1.0e-12_dp) error stop 'order.row objective mismatch'
   call binary_gradient(pcol, dat, pri, gcol, .false.)
   call binary_gradient(prow, dat, pri, grow, .true.)
   do j = 1, n
      do i = 1, k
         if (abs(grow((i-1)*n+j)-gcol((j-1)*k+i)) > 1.0e-12_dp) then
            error stop 'order.row gradient mismatch'
         end if
      end do
   end do
   if (maxval(abs(grow(n*k+1:)-gcol(n*k+1:))) > 1.0e-12_dp) then
      error stop 'order.row mu-gradient mismatch'
   end if
   write(*,*) 'PASS binary order.row parameter layout'
end program test_binary_ordering
