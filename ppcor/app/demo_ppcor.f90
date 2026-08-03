program demo_ppcor
   use ppcor, only : dp, ppcor_result, ppcor_test_result, pcor, spcor, &
                     pcor_test, ppcor_pearson, ppcor_spearman, ppcor_kendall, &
                     method_name
   implicit none
   real(dp) :: data(40,4)
   type(ppcor_result) :: result
   type(ppcor_test_result) :: pair
   integer :: method, k
   integer, parameter :: methods(2) = [ppcor_pearson, ppcor_spearman]

   call make_data(data)

   print '(a)', 'ppcor-fortran demonstration'
   print '(a)', 'Variables are x, y, z1, z2; rows are observations.'
   print '(a)', ''
   do k = 1, size(methods)
      method = methods(k)
      call pcor(data, result, method)
      print '(a,a)', trim(method_name(method)), ' partial-correlation matrix:'
      call print_matrix(result%estimate)
      print '(a,f10.6,a,es11.4)', 'pcor(x,y | z1,z2) = ', result%estimate(1,2), &
            ', p = ', result%p_value(1,2)
      print '(a)', ''
   end do

   call pcor(data, result, ppcor_kendall)
   print '(a,f10.6,a,es11.4)', 'Kendall partial x,y = ', result%estimate(1,2), &
         ', p = ', result%p_value(1,2)

   call spcor(data, result, ppcor_pearson)
   print '(a)', ''
   print '(a)', 'Pearson semi-partial matrix (directional):'
   call print_matrix(result%estimate)

   call pcor_test(data(:,1), data(:,2), data(:,3:4), pair, ppcor_pearson)
   print '(a)', ''
   print '(a,f10.6,a,f10.6,a,es11.4)', 'Pairwise pcor.test: r = ', pair%estimate, &
         ', t = ', pair%statistic, ', p = ', pair%p_value

contains

   subroutine make_data(a)
      real(dp), intent(out) :: a(:,:)
      integer :: i
      real(dp) :: t, z1, z2, x, y
      do i = 1, size(a,1)
         t = real(i,dp)
         z1 = sin(0.21_dp*t)
         z2 = cos(0.13_dp*t) + 0.02_dp*t
         x = 0.7_dp*z1 - 0.35_dp*z2 + 0.55_dp*sin(0.91_dp*t)
         y = 0.25_dp*z1 + 0.45_dp*z2 + 0.5_dp*x + 0.65_dp*cos(0.57_dp*t)
         a(i,:) = [x, y, z1, z2]
      end do
   end subroutine make_data

   subroutine print_matrix(a)
      real(dp), intent(in) :: a(:,:)
      integer :: i
      do i = 1, size(a,1)
         print '(*(f11.6,1x))', a(i,:)
      end do
   end subroutine print_matrix

end program demo_ppcor
