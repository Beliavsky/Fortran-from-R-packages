program test_models
  use ghyp_kinds, only : dp
  use tsgarch
  use test_support
  implicit none
  character(len=12),parameter::models(8)=[character(len=12)::'garch','gjrgarch','aparch','egarch','fgarch','cgarch','igarch','ewma']
  character(len=8),parameter::dists(10)=[character(len=8)::'norm','std','snorm','sstd','ged','sged','nig','gh','jsu','ghst']
  real(dp)::y(180)
  type(garch_spec)::spec
  type(garch_parameters)::par
  type(garch_filter_result)::f
  integer::i
  y=[(0.02_dp*sin(0.17_dp*real(i,dp))+0.01_dp*cos(0.071_dp*real(i*i,dp)),i=1,size(y))]
  do i=1,size(models)
    spec=standard_spec(trim(models(i)))
    par=standard_parameters(y,spec)
    f=filter_garch(y,spec,par)
    call assert_true(f%status==tsg_success,'filter failed for '//trim(models(i))//' '//trim(f%message))
    call assert_true(all(f%variance>0.0_dp),'nonpositive variance for '//trim(models(i)))
  end do
  spec=standard_spec('garch')
  do i=1,size(dists)
    spec%distribution=dists(i)
    par=standard_parameters(y,spec)
    par=initialize_parameters(y,spec)
    par%omega=0.02_dp
    par%alpha=0.06_dp
    par%beta=0.90_dp
    f=filter_garch(y,spec,par)
    call assert_true(f%status==tsg_success,'distribution filter failed for '//trim(dists(i)))
  end do
  call assert_true(index(model_equation(spec),'sigma')>0,'equation text missing')
  write(*,'(a)')'test_models: PASS'
end program test_models
