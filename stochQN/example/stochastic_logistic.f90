program stochastic_logistic
   use stochqn_kinds, only : dp
   use stochqn_logistic, only : stochastic_logistic_t, logistic_olbfgs
   implicit none

   type(stochastic_logistic_t) :: model
   real(dp) :: x(12, 2), y(12)
   real(dp), allocatable :: coefficients(:), probability(:)
   integer, allocatable :: classification(:)
   integer :: i, stat

   do i = 1, 12
      x(i, 1) = -2.0_dp + 4.0_dp * real(i - 1, dp) / 11.0_dp
      x(i, 2) = sin(real(i, dp))
      y(i) = merge(1.0_dp, 0.0_dp, 1.5_dp * x(i, 1) - 0.3_dp * x(i, 2) > 0.0_dp)
   end do

   call model%initialize(2, optimizer_kind=logistic_olbfgs, fit_intercept=.true., &
                         lambda=1.0e-4_dp, initial_step=0.2_dp, stat=stat)
   if (stat /= 0) error stop 'Could not initialize model.'

   do i = 1, 150
      call model%partial_fit(x, y, stat=stat)
      if (stat /= 0) error stop 'Partial fit failed.'
   end do

   coefficients = model%get_coefficients()
   probability = model%predict_probability(x)
   classification = model%predict_class(x)

   print '(a,*(f12.6,1x))', 'coefficients: ', coefficients
   print '(a,i0,a,i0)', 'correct classifications: ', count(classification == int(y)), &
                        ' of ', size(y)
   print '(a,*(f8.4,1x))', 'probabilities: ', probability
end program stochastic_logistic
