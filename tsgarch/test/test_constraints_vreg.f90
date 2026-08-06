program test_constraints_vreg
  use ghyp_kinds, only : dp
  use tsgarch
  use test_support
  implicit none
  real(dp)::y(100),vreg(100,1)
  real(dp),allocatable::packed(:)
  type(garch_spec)::spec
  type(garch_parameters)::par,copy
  type(garch_filter_result)::f
  integer::i,status
  y=[(0.02_dp*sin(0.1_dp*real(i,dp)),i=1,100)]
  vreg(:,1)=[(sin(0.03_dp*real(i,dp)),i=1,100)]
  spec=standard_spec('garch','norm')
  spec%variance_targeting=.true.
  par=initialize_parameters(y,spec,1)
  par%alpha=0.05_dp
  par%beta=0.90_dp
  par%xi(1)=1.0e-6_dp
  f=filter_garch(y,spec,par,vreg)
  if(f%status/=tsg_success)write(*,'(a)')trim(f%message)
  call assert_true(f%status==tsg_success,'variance-regressor filter failed')
  call pack_parameters(spec,par,packed)
  call unpack_parameters(spec,packed,par,copy,status)
  call assert_true(status==tsg_success.and.abs(copy%xi(1)-par%xi(1))<1e-12_dp,'pack/unpack with vreg')
  spec=standard_spec('igarch','norm')
  par=standard_parameters(y,spec)
  call pack_parameters(spec,par,packed)
  call unpack_parameters(spec,packed,par,copy,status)
  call assert_true(abs(sum(copy%alpha)+sum(copy%beta)-1.0_dp)<1e-12_dp,'IGARCH equality packing')
  write(*,'(a)')'test_constraints_vreg: PASS'
end program test_constraints_vreg
