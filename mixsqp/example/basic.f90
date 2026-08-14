program basic
  use mixsqp
  implicit none
  real(dp),allocatable::xdat(:),s(:),L(:,:)
  type(mixsqp_result)::fit
  type(mixsqp_control)::ctl
  call set_seed(1)
  call simulate_mix_data(1000,10,xdat,s,L)
  ctl=mixsqp_default_control();ctl%verbose=.false.
  call fit_mixsqp(L,fit,control=ctl)
  print '(a,a)', 'status: ', trim(fit%status_message)
  print '(a,es14.6)', 'objective: ', fit%value
  print '(a,10(f8.5,1x))', 'weights: ', fit%x
end program
