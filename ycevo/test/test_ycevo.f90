program test_ycevo
   use ycevo, only : dp, ycevo_success, bond_panel_t, yield_curve_t
   use ycevo, only : epanechnikov, epanechnikov_quantile, calc_epaker_weights
   use ycevo, only : get_yield_at, generate_yield, calc_dbar, estimate_yield
   use ycevo, only : loess_predict, linear_interpolate, bilinear_interpolate
   implicit none

   integer :: failures

   failures = 0
   call test_kernel(failures)
   call test_yield_model(failures)
   call test_dbar_reference(failures)
   call test_zero_coupon_estimator(failures)
   call test_coupon_estimator(failures)
   call test_prediction(failures)

   if (failures > 0) then
      write(*, '(a,i0)') 'FAILED tests: ', failures
      error stop 1
   end if
   write(*, '(a)') 'All ycevo tests passed.'

contains

   subroutine test_kernel(failures)
      integer, intent(inout) :: failures
      real(dp), allocatable :: w(:, :)
      real(dp) :: gamma(11), grid(2), bandwidth(2), expected(11,2)
      integer :: i

      gamma = [(real(i-1,dp)/10.0_dp, i=1,11)]
      grid = [0.2_dp, 0.4_dp]
      bandwidth = [0.2_dp, 0.4_dp]
      call calc_epaker_weights(gamma, grid, bandwidth, w)
      expected(:,1) = [0.0_dp, 0.5625_dp, 0.75_dp, 0.5625_dp, 0.0_dp, &
                       0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
      expected(:,2) = [0.0_dp, 0.328125_dp, 0.5625_dp, 0.703125_dp, 0.75_dp, &
                       0.703125_dp, 0.5625_dp, 0.328125_dp, 0.0_dp, 0.0_dp, 0.0_dp]
      call assert_matrix_close('Epanechnikov weights', w, expected, 1.0e-14_dp, failures)
      call assert_close('kernel centre', epanechnikov(0.0_dp), 0.75_dp, 1.0e-15_dp, failures)
      call assert_close('kernel quantile median', &
         epanechnikov_quantile(0.5_dp, 4.0_dp, 2.0_dp), 4.0_dp, 1.0e-15_dp, failures)
   end subroutine test_kernel

   subroutine test_yield_model(failures)
      integer, intent(inout) :: failures
      real(dp), parameter :: expected(10) = [ &
         0.0337921613338867_dp, 0.0322856630629355_dp, 0.0309356840928623_dp, &
         0.0296248349480087_dp, 0.0282357261527160_dp, 0.0266509682313258_dp, &
         0.0247531717081795_dp, 0.0224249471076185_dp, 0.0195489049539844_dp, &
         0.0160076557716185_dp]
      real(dp), allocatable :: generated(:, :)
      real(dp) :: actual(10)
      integer :: i

      do i = 1, 10
         actual(i) = get_yield_at(real(i,dp)/10.0_dp, 1.0_dp)
      end do
      call assert_vector_close('get_yield_at R reference', actual, expected, 2.0e-14_dp, failures)
      call generate_yield(generated)
      call assert_close('generate_yield shape value', generated(18,6), &
         get_yield_at(0.5_dp, 5.0_dp), 2.0e-14_dp, failures)
   end subroutine test_yield_model

   subroutine test_dbar_reference(failures)
      integer, intent(inout) :: failures
      type(bond_panel_t) :: panel
      real(dp), allocatable :: num(:), den(:)
      real(dp) :: tau(4), ht(4), discount(10), price_by_id(4)
      real(dp), parameter :: expected_num(4) = [ &
         5678.94832046269_dp, 5705.41441494833_dp, &
         5625.45900788638_dp, 5538.31638436822_dp]
      real(dp), parameter :: expected_den(4) = [ &
         5626.6875_dp, 5739.1875_dp, 5738.625_dp, 5738.0625_dp]
      integer :: status, i

      panel%nday = 1
      panel%day = [1,1,1,1,1,1,1,1,1,1]
      panel%id = [1,2,2,3,3,3,4,4,4,4]
      panel%tupq = [180,360,180,540,360,180,720,540,360,180]
      panel%cashflow = [100.0_dp,101.0_dp,1.0_dp,101.0_dp,1.0_dp,1.0_dp, &
                        101.0_dp,1.0_dp,1.0_dp,1.0_dp]
      do i = 1, 10
         discount(i) = exp(-get_yield_at(0.0_dp, real(panel%tupq(i),dp)/365.0_dp)* &
                            real(panel%tupq(i),dp)/365.0_dp)
      end do
      price_by_id(1) = panel%cashflow(1)*discount(1)
      price_by_id(2) = sum(panel%cashflow(2:3)*discount(2:3))
      price_by_id(3) = sum(panel%cashflow(4:6)*discount(4:6))
      price_by_id(4) = sum(panel%cashflow(7:10)*discount(7:10))
      allocate(panel%price(10))
      panel%price = [price_by_id(1), price_by_id(2), price_by_id(2), &
                     price_by_id(3), price_by_id(3), price_by_id(3), &
                     price_by_id(4), price_by_id(4), price_by_id(4), price_by_id(4)]
      tau = [180.0_dp,360.0_dp,540.0_dp,720.0_dp]/365.0_dp
      ht = 180.0_dp/365.0_dp
      call calc_dbar(panel, 1.0_dp, 1.0_dp, tau, ht, num, den, status=status)
      call assert_true('calc_dbar status', status == ycevo_success, failures)
      call assert_vector_close('calc_dbar numerator R reference', num, expected_num, 2.0e-10_dp, failures)
      call assert_vector_close('calc_dbar denominator R reference', den, expected_den, 2.0e-12_dp, failures)
   end subroutine test_dbar_reference

   subroutine test_zero_coupon_estimator(failures)
      integer, intent(inout) :: failures
      type(bond_panel_t) :: panel
      type(yield_curve_t) :: curve
      real(dp) :: tau(3), ht(3), expected(3), interest(1)
      integer :: status
      character(len=256) :: message

      panel%nday = 1
      panel%day = [1,1,1]
      panel%id = [1,2,3]
      panel%tupq = [365,730,1095]
      panel%cashflow = [100.0_dp,100.0_dp,100.0_dp]
      panel%price = [95.0_dp,90.0_dp,85.0_dp]
      tau = [1.0_dp,2.0_dp,3.0_dp]
      ht = [0.49_dp,0.49_dp,0.49_dp]
      expected = [0.95_dp,0.90_dp,0.85_dp]
      call estimate_yield(panel, 1.0_dp, 0.5_dp, tau, ht, curve, status, message)
      call assert_true('zero coupon status', status == ycevo_success, failures, message)
      if (status == ycevo_success) then
         call assert_vector_close('zero coupon discount', curve%discount, expected, 2.0e-14_dp, failures)
      end if

      interest = [4.0_dp]
      call estimate_yield(panel, 1.0_dp, 0.5_dp, tau, ht, curve, status, message, &
                          interest=interest, rgrid=4.0_dp, hr=0.5_dp)
      call assert_true('covariate estimator status', status == ycevo_success, failures, message)
   end subroutine test_zero_coupon_estimator


   subroutine test_coupon_estimator(failures)
      integer, intent(inout) :: failures
      type(bond_panel_t) :: panel
      type(yield_curve_t) :: curve
      real(dp) :: tau(3), ht(3), expected(3), p1, p2, p3
      integer :: status
      character(len=256) :: message

      expected = [0.96_dp, 0.90_dp, 0.84_dp]
      p1 = 100.0_dp*expected(1)
      p2 = 5.0_dp*expected(1) + 105.0_dp*expected(2)
      p3 = 5.0_dp*expected(1) + 5.0_dp*expected(2) + 105.0_dp*expected(3)
      panel%nday = 1
      panel%day = [1,1,1,1,1,1]
      panel%id = [1,2,2,3,3,3]
      panel%tupq = [365,365,730,365,730,1095]
      panel%cashflow = [100.0_dp,5.0_dp,105.0_dp,5.0_dp,5.0_dp,105.0_dp]
      panel%price = [p1,p2,p2,p3,p3,p3]
      tau = [1.0_dp,2.0_dp,3.0_dp]
      ht = 0.49_dp
      call estimate_yield(panel, 1.0_dp, 0.5_dp, tau, ht, curve, status, message)
      call assert_true('coupon estimator status', status == ycevo_success, failures, message)
      if (status == ycevo_success) then
         expected = [0.981638356239861_dp, 0.920379756358323_dp, &
                     0.860642210264351_dp]
         call assert_vector_close('coupon estimator R-kernel result', curve%discount, expected, &
                                  3.0e-14_dp, failures)
      end if
   end subroutine test_coupon_estimator

   subroutine test_prediction(failures)
      integer, intent(inout) :: failures
      real(dp) :: x(11), y(11), xout(3), expected(3), grid(2,2)
      real(dp), allocatable :: actual(:)
      integer :: i, status

      do i = 1, 11
         x(i) = real(i-1,dp)/10.0_dp
         y(i) = 1.0_dp + 2.0_dp*x(i) + 3.0_dp*x(i)**2
      end do
      xout = [0.15_dp,0.45_dp,0.85_dp]
      expected = 1.0_dp + 2.0_dp*xout + 3.0_dp*xout**2
      call loess_predict(x, y, xout, actual, status)
      call assert_true('loess status', status == ycevo_success, failures)
      call assert_vector_close('local quadratic reproduction', actual, expected, 2.0e-11_dp, failures)
      call assert_close('linear interpolation', linear_interpolate([0.0_dp,1.0_dp], &
         [2.0_dp,4.0_dp], 0.25_dp), 2.5_dp, 1.0e-15_dp, failures)
      grid = reshape([0.0_dp,1.0_dp,2.0_dp,3.0_dp], [2,2])
      call assert_close('bilinear interpolation', bilinear_interpolate([0.0_dp,1.0_dp], &
         [0.0_dp,1.0_dp], grid, 0.5_dp, 0.5_dp), 1.5_dp, 1.0e-15_dp, failures)
   end subroutine test_prediction

   subroutine assert_close(name, actual, expected, tolerance, failures)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: actual, expected, tolerance
      integer, intent(inout) :: failures

      if (abs(actual - expected) > tolerance*max(1.0_dp,abs(expected))) then
         write(*, '(a,2(1x,es24.16))') 'FAIL '//trim(name)//':', actual, expected
         failures = failures + 1
      end if
   end subroutine assert_close

   subroutine assert_vector_close(name, actual, expected, tolerance, failures)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: actual(:), expected(:), tolerance
      integer, intent(inout) :: failures
      integer :: i

      if (size(actual) /= size(expected)) then
         write(*, '(a)') 'FAIL '//trim(name)//': size mismatch'
         failures = failures + 1
         return
      end if
      do i = 1, size(actual)
         if (abs(actual(i)-expected(i)) > tolerance*max(1.0_dp,abs(expected(i)))) then
            write(*, '(a,i0,2(1x,es24.16))') 'FAIL '//trim(name)//' at ', i, actual(i), expected(i)
            failures = failures + 1
            return
         end if
      end do
   end subroutine assert_vector_close

   subroutine assert_matrix_close(name, actual, expected, tolerance, failures)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: actual(:,:), expected(:,:), tolerance
      integer, intent(inout) :: failures
      integer :: i, j

      if (any(shape(actual) /= shape(expected))) then
         write(*, '(a)') 'FAIL '//trim(name)//': shape mismatch'
         failures = failures + 1
         return
      end if
      do j = 1, size(actual,2)
         do i = 1, size(actual,1)
            if (abs(actual(i,j)-expected(i,j)) > tolerance*max(1.0_dp,abs(expected(i,j)))) then
               write(*, '(a,2(i0,1x),2(es24.16,1x))') 'FAIL '//trim(name)//' at ', &
                  i, j, actual(i,j), expected(i,j)
               failures = failures + 1
               return
            end if
         end do
      end do
   end subroutine assert_matrix_close

   subroutine assert_true(name, condition, failures, detail)
      character(len=*), intent(in) :: name
      logical, intent(in) :: condition
      integer, intent(inout) :: failures
      character(len=*), intent(in), optional :: detail

      if (.not. condition) then
         if (present(detail)) then
            write(*, '(a,a)') 'FAIL '//trim(name)//': ', trim(detail)
         else
            write(*, '(a)') 'FAIL '//trim(name)
         end if
         failures = failures + 1
      end if
   end subroutine assert_true

end program test_ycevo
