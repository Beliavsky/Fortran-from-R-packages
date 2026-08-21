program test_imputation
  use compositions
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
  implicit none
  real(dp) :: comp(4,3),pred(4,3),cov(3,3),dl(4,3),x(4,2),ratio
  integer :: mt(4,3),i
  integer, allocatable :: mtc(:,:)
  real(dp), allocatable :: dlc(:,:)
  real(dp) :: cv(1,6)
  type(acomp_imputation_result) :: imp,fit,em
  cv(1,:)=[1.0_dp,0.0_dp,-0.2_dp,ieee_value(0.0_dp,ieee_quiet_nan), &
    ieee_value(0.0_dp,ieee_positive_inf),ieee_value(0.0_dp,ieee_negative_inf)]
  call classify_missingness(cv,mtc,dlc,0.1_dp)
  if(any(mtc(1,:)/=[mt_observed,mt_bdl,mt_bdl,mt_mar,mt_error,mt_sz])) error stop 'missing classification'
  if(abs(dlc(1,2)-0.1_dp)>1.0e-15_dp.or.abs(dlc(1,3)-0.2_dp)>1.0e-15_dp) error stop 'detection limits'
  comp=reshape([0.5_dp,0.3_dp,0.2_dp, 0.4_dp,0.4_dp,0.2_dp, &
                0.3_dp,0.5_dp,0.2_dp, 0.2_dp,0.5_dp,0.3_dp],[4,3],order=[2,1])
  x(:,1)=1.0_dp; x(:,2)=[-1.0_dp,-0.3_dp,0.4_dp,1.0_dp]
  pred=0.0_dp
  do i=1,4
    pred(i,:)=[0.2_dp*x(i,2),-0.1_dp*x(i,2),-0.1_dp*x(i,2)]
  end do
  cov=reshape([0.20_dp,-0.10_dp,-0.10_dp,-0.10_dp,0.18_dp,-0.08_dp, &
               -0.10_dp,-0.08_dp,0.18_dp],[3,3])
  mt=mt_observed; dl=0.0_dp
  mt(2,1)=mt_mar; comp(2,1)=0.0_dp
  mt(3,1)=mt_bdl; dl(3,1)=0.20_dp; comp(3,1)=0.0_dp
  imp=impute_acomp_conditional(comp,pred,cov,mt,dl,nsim=400,seed=17)
  if(any(abs(sum(imp%composition,dim=2)-1.0_dp)>1.0e-12_dp)) error stop 'imputation closure'
  ratio=imp%composition(2,2)/imp%composition(2,3)
  if(abs(ratio-(comp(2,2)/comp(2,3)))>1.0e-10_dp) error stop 'observed ratio not preserved'
  if(imp%composition(3,1)/imp%composition(3,3)>dl(3,1)/0.2_dp+1.0e-10_dp) error stop 'BDL constraint'
  fit=fit_acomp_projection(reshape([0.5_dp,0.3_dp,0.2_dp,0.4_dp,0.4_dp,0.2_dp, &
       0.3_dp,0.5_dp,0.2_dp,0.2_dp,0.5_dp,0.3_dp],[4,3],order=[2,1]),x,0*mt,source_compatible=.false.)
  if(.not.fit%ok.or.any(abs(sum(fit%composition,dim=2)-1.0_dp)>1.0e-12_dp)) error stop 'projection fit'
  em=fit_acomp_em(comp,x,mt,dl,steps=2,seed=19)
  if(.not.em%ok.or.any(abs(sum(em%composition,dim=2)-1.0_dp)>1.0e-12_dp)) error stop 'EM imputation'
  print *, 'test_imputation: PASS'
end program
