program test_pde
   use ragtop
   implicit none
   type(grid_spec) :: grid
   type(market_spec) :: market
   type(instrument_spec) :: option
   type(option_value) :: exact
   real(dp) :: drift(5), sub(4), diag(5), super(4)
   real(dp) :: x_true(5), rhs(5), x(5), pde_price
   integer :: status, i

   call construct_implicit_grid_structure(5.0_dp,100,33.0_dp,33.0_dp, &
      0.01_dp,0.5_dp,0.25_dp,2.0_dp,grid,status)
   call assert_int(status,ragtop_ok,'grid status')
   call assert_close(grid%dt,0.05_dp,1.0e-14_dp,'grid dt')
   call assert_close(grid%dz,0.447213595499958_dp,1.0e-14_dp,'grid dz')
   call assert_close(grid%z0,-0.575_dp,1.0e-14_dp,'grid center')
   call assert_int(grid%n_space,11,'grid size')

   drift = 0.0011_dp*[2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp]
   call construct_tridiagonals(0.5_dp,0.05_dp,drift,sub,diag,super,status)
   call assert_vector(diag,[1.0022_dp,1.0125_dp,1.0125_dp,1.0125_dp, &
                            1.0066_dp],1.0e-12_dp,'diagonal')
   call assert_vector(super,[-0.0022_dp,-0.0079_dp,-0.00845_dp,-0.009_dp], &
                      1.0e-12_dp,'superdiagonal')
   call assert_vector(sub,[-0.0046_dp,-0.00405_dp,-0.0035_dp,-0.0066_dp], &
                      1.0e-12_dp,'subdiagonal')

   x_true = [0.2_dp,-0.5_dp,1.0_dp,0.7_dp,-0.1_dp]
   rhs = diag*x_true
   do i = 1, 4
      rhs(i) = rhs(i)+super(i)*x_true(i+1)
      rhs(i+1) = rhs(i+1)+sub(i)*x_true(i)
   end do
   call pde_matrix_solve(sub,diag,super,rhs,x,status)
   call assert_vector(x,x_true,1.0e-12_dp,'tridiagonal solve')

   market%short_rate = 0.0_dp
   market%default_intensity = 0.07_dp
   market%volatility = 0.5_dp
   option = EuropeanOption(1.0_dp,90.0_dp,put_option)
   pde_price = find_present_value(100.0_dp,100,option,market, &
                                  std_devs_width=4.0_dp)
   exact = black_scholes_on_term_structures(put_option,100.0_dp,90.0_dp, &
                                             1.0_dp,market)
   call assert_close(pde_price,exact%price,0.03_dp,'European PDE')

   print '(a)', 'test_pde: PASS'
contains
   subroutine assert_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: label
      if (abs(actual-expected) > tolerance) then
         write(*,'(a,2es24.14)') trim(label)//' failed: ',actual,expected
         error stop 1
      end if
   end subroutine assert_close
   subroutine assert_vector(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual(:), expected(:), tolerance
      character(len=*), intent(in) :: label
      if (maxval(abs(actual-expected)) > tolerance) then
         write(*,'(a,es24.14)') trim(label)//' max error: ', &
                                maxval(abs(actual-expected))
         error stop 1
      end if
   end subroutine assert_vector
   subroutine assert_int(actual, expected, label)
      integer, intent(in) :: actual, expected
      character(len=*), intent(in) :: label
      if (actual /= expected) error stop trim(label)
   end subroutine assert_int
end program test_pde
