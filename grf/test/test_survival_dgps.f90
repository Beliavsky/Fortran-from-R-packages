program test_survival_dgps
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
   use grf, only : dp, generate_causal_survival_data
   implicit none

   integer, parameter :: n = 36
   integer, parameter :: p = 5
   character(len=7), parameter :: designs(6) = [character(len=7) :: &
      'simple1', 'type1', 'type2', 'type3', 'type4', 'type5']
   real(dp), parameter :: default_horizon(6) = [1.0_dp, 1.5_dp, 2.0_dp, 15.0_dp, 3.0_dp, 2.0_dp]
   real(dp), allocatable :: x(:,:)
   real(dp), allocatable :: x_copy(:,:)
   real(dp), allocatable :: y(:)
   real(dp), allocatable :: y_copy(:)
   real(dp), allocatable :: cate(:)
   real(dp), allocatable :: cate_copy(:)
   real(dp), allocatable :: cate_probability(:)
   real(dp), allocatable :: cate_probability_copy(:)
   real(dp), allocatable :: cate_sign(:)
   real(dp), allocatable :: supplied_x(:,:)
   integer, allocatable :: treatment(:)
   integer, allocatable :: treatment_copy(:)
   integer, allocatable :: event(:)
   integer, allocatable :: event_copy(:)
   integer :: failures
   integer :: info
   integer :: i
   integer :: j
   integer :: k

   failures = 0
   do k = 1, size(designs)
      call generate_causal_survival_data(n, p, x, y, treatment, event, cate, &
                                         cate_probability, info, seed=2100 + k, &
                                         dgp=trim(designs(k)), n_mc=250, cate_sign=cate_sign)
      call check(info == 0, trim(designs(k)) // ' accepted', failures)
      if (info /= 0) cycle
      call check(all(shape(x) == [n,p]) .and. size(y) == n, &
                 trim(designs(k)) // ' output shapes', failures)
      call check(all(ieee_is_finite(x)) .and. all(ieee_is_finite(y)), &
                 trim(designs(k)) // ' finite observations', failures)
      call check(all(y >= 0.0_dp) .and. all(y <= default_horizon(k)), &
                 trim(designs(k)) // ' observed-time range', failures)
      call check(all(treatment == 0 .or. treatment == 1) .and. &
                 all(event == 0 .or. event == 1), &
                 trim(designs(k)) // ' binary indicators', failures)
      call check(all(ieee_is_finite(cate)) .and. all(ieee_is_finite(cate_probability)), &
                 trim(designs(k)) // ' finite target effects', failures)
      call check(all(cate_probability >= -1.0_dp) .and. all(cate_probability <= 1.0_dp), &
                 trim(designs(k)) // ' probability-effect range', failures)
      call check(all(ieee_is_finite(cate_sign) .or. ieee_is_nan(cate_sign)), &
                 trim(designs(k)) // ' effect signs', failures)
   end do

   call generate_causal_survival_data(n, p, x, y, treatment, event, cate, &
                                      cate_probability, info, seed=765, dgp='type2', &
                                      n_mc=180, cate_sign=cate_sign)
   x_copy = x
   y_copy = y
   treatment_copy = treatment
   event_copy = event
   cate_copy = cate
   cate_probability_copy = cate_probability
   call generate_causal_survival_data(n, p, x, y, treatment, event, cate, &
                                      cate_probability, info, seed=765, dgp='type2', &
                                      n_mc=180, cate_sign=cate_sign)
   call check(maxval(abs(x - x_copy)) <= 0.0_dp .and. &
              maxval(abs(y - y_copy)) <= 0.0_dp, &
              'repeatable survival observations', failures)
   call check(all(treatment == treatment_copy) .and. all(event == event_copy), &
              'repeatable survival indicators', failures)
   call check(maxval(abs(cate - cate_copy)) <= 0.0_dp .and. &
              maxval(abs(cate_probability - cate_probability_copy)) <= 0.0_dp, &
              'repeatable survival targets', failures)

   allocate(supplied_x(4,p))
   do j = 1, p
      do i = 1, size(supplied_x,1)
         supplied_x(i,j) = real(i + j,dp) / 12.0_dp
      end do
   end do
   call generate_causal_survival_data(1, 1, x, y, treatment, event, cate, &
                                      cate_probability, info, seed=55, dgp='type1', &
                                      n_mc=100, x_input=supplied_x)
   call check(info == 0 .and. all(shape(x) == shape(supplied_x)), &
              'caller-supplied predictors', failures)
   call check(maxval(abs(x - supplied_x)) <= 0.0_dp, &
              'supplied predictors retained', failures)

   call generate_causal_survival_data(n, p, x, y, treatment, event, cate, &
                                      cate_probability, info, seed=91, rho=0.65_dp, &
                                      n_mc=20)
   call check(info == 0 .and. all(x >= 0.0_dp) .and. all(x <= 1.0_dp), &
              'correlated uniform-marginal predictors', failures)

   call generate_causal_survival_data(10, 5, x, y, treatment, event, cate, &
                                      cate_probability, info, dgp='unknown')
   call check(info == -3, 'unknown survival selector validation', failures)
   call generate_causal_survival_data(10, 4, x, y, treatment, event, cate, &
                                      cate_probability, info, dgp='type1')
   call check(info == -2, 'survival predictor-count validation', failures)
   call generate_causal_survival_data(10, 5, x, y, treatment, event, cate, &
                                      cate_probability, info, y_max=0.0_dp)
   call check(info == -4, 'survival horizon validation', failures)
   call generate_causal_survival_data(10, 5, x, y, treatment, event, cate, &
                                      cate_probability, info, rho=1.1_dp)
   call check(info == -5, 'survival correlation validation', failures)
   call generate_causal_survival_data(10, 5, x, y, treatment, event, cate, &
                                      cate_probability, info, n_mc=0)
   call check(info == -6, 'survival Monte Carlo validation', failures)

   if (failures /= 0) error stop 'GRF causal-survival DGP tests failed'
   print '(a)', 'All GRF causal-survival simulation DGP tests passed.'

contains

   subroutine check(condition, label, failures)
      logical, intent(in) :: condition !! Assertion result.
      character(len=*), intent(in) :: label !! Human-readable assertion label.
      integer, intent(inout) :: failures !! Running failure count.

      if (.not. condition) then
         failures = failures + 1
         print '(a,a)', 'FAIL: ', trim(label)
      end if
   end subroutine check

end program test_survival_dgps
