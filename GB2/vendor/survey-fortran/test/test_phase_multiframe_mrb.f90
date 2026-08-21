program test_phase_multiframe_mrb
  use survey_kinds, only : dp
  use survey_types, only : survey_design_t, rep_design_t, multiframe_design_t, phase_variance_t, svystat_t
  use survey_design, only : make_design
  use survey_phase, only : dcheck_strat, dcheck_multi, dcheck_multi_subset, combine_dcheck, &
    twophase2_variance, multiphase_variance, multiphase_total, multiphase_mean, project_phase_calibration
  use survey_multiframe, only : make_multiframe_constant, make_multiframe_expected, multiframe_total, multiframe_mean
  use survey_mrb, only : make_mrb
  implicit none
  real(dp) :: x(3,1),dfull(3,3),d2(3,3),dcube(3,3,2),pw(3,2),fw(3),expected, tol
  real(dp) :: x1(2,1),x2(2,1),dc1(2,2),dc2(2,2),w1(2),w2(2),ov1(2,2),ov2(2,2)
  real(dp) :: dcs(4,4),dcm(4,4),dca(4,4),dcb(4,4),dcc(4,4),probs(4,2),dcsub(2,2)
  real(dp) :: qrx(3,1),cscale(3),calx(3,1),calres(3,1)
  integer :: ids2(4,2),str2(4,2)
  logical :: subset4(4)
  type(phase_variance_t) :: pv,pvm
  type(svystat_t) :: st,sm,mft,mfm
  type(multiframe_design_t) :: mf,mfe
  logical :: of1(2),of2(2)
  type(survey_design_t) :: d
  type(rep_design_t) :: rep
  integer :: cl(4,1),str(4,1),i,r
  real(dp) :: bw(4),ss(4,1),ps(4,1)
  integer,allocatable :: seed(:)

  tol=1e-11_dp
  call dcheck_strat([1,1,2,2],[0.5_dp,0.5_dp,0.25_dp,0.25_dp],dcs)
  call assert_close(dcs(1,1),0.5_dp,tol,'Dcheck strat diag')
  call assert_close(dcs(1,2),-0.5_dp,tol,'Dcheck strat offdiag')
  call assert_close(dcs(3,4),-0.75_dp,tol,'Dcheck strat second stratum')
  call assert_close(dcs(1,3),0.0_dp,tol,'Dcheck strat cross-stratum')
  ids2(:,1)=[1,2,3,4];ids2(:,2)=[11,12,13,14];str2=1;probs(:,1)=0.5_dp;probs(:,2)=0.25_dp
  call dcheck_multi(ids2,str2,probs,dcm)
  call dcheck_strat(str2(:,1),probs(:,1),dca);call dcheck_strat(str2(:,2),probs(:,2),dcb);call combine_dcheck(dca,dcb,dcc)
  if(maxval(abs(dcm-dcc))>tol) error stop 'Dcheck multistage combination'
  subset4=[.true.,.true.,.false.,.false.]
  call dcheck_multi_subset(ids2(:,1:1),str2(:,1:1),subset4,probs(:,1:1),dcsub)
  call assert_close(dcsub(1,2),-(1.0_dp-0.5_dp)/3.0_dp,tol,'Dcheck subset uses full stratum sample size')
  qrx(:,1)=[1.0_dp,2.0_dp,3.0_dp];cscale=[1.0_dp,1.5_dp,2.0_dp];calx(:,1)=2.5_dp*qrx(:,1)*cscale
  call project_phase_calibration(calx,qrx,cscale,calres)
  if(maxval(abs(calres))>1e-10_dp) error stop 'phase calibration projection'
  x(:,1)=[1.0_dp,2.0_dp,4.0_dp]
  dfull=0.0_dp;d2=0.0_dp
  do i=1,3;dfull(i,i)=1.0_dp;d2(i,i)=0.25_dp;end do
  call twophase2_variance(x,dfull,d2,pv)
  expected=sum(x(:,1)**2)
  call assert_close(pv%variance(1,1),expected,tol,'twophase full variance')
  call assert_close(pv%phase(1,1,1)+pv%phase(1,1,2),pv%variance(1,1),tol,'twophase decomposition')
  call assert_close(pv%phase(1,1,2),0.25_dp*expected,tol,'twophase phase2')

  pw(:,1)=[1.0_dp,2.0_dp,1.0_dp];pw(:,2)=[0.5_dp,1.0_dp,1.5_dp]
  dcube=0.0_dp
  do i=1,3;dcube(i,i,1)=1.0_dp;dcube(i,i,2)=0.5_dp;end do
  call multiphase_variance(x,pw,dcube,pvm)
  expected=sum((x(:,1)*pw(:,1))**2)+0.5_dp*sum((x(:,1)*pw(:,2))**2)
  call assert_close(pvm%variance(1,1),expected,tol,'multiphase summed variance')
  call assert_close(sum(pvm%phase(1,1,:)),expected,tol,'multiphase phase sum')
  fw=[2.0_dp,1.0_dp,0.5_dp]
  st=multiphase_total(x,fw,pw,dcube)
  call assert_close(st%estimate(1),6.0_dp,tol,'multiphase total')
  dcube(:,:,2)=dcube(:,:,1)
  sm=multiphase_mean(x,fw,pw,dcube)
  call assert_close(sm%estimate(1),6.0_dp/3.5_dp,tol,'multiphase mean')

  x1(:,1)=[10.0_dp,20.0_dp];x2(:,1)=[20.0_dp,30.0_dp];w1=1.0_dp;w2=1.0_dp
  of1=[.false.,.true.];of2=[.true.,.false.];dc1=0.0_dp;dc2=0.0_dp
  do i=1,2;dc1(i,i)=1.0_dp;dc2(i,i)=1.0_dp;end do
  call make_multiframe_constant(w1,w2,of1,of2,mf)
  mft=multiframe_total(x1,x2,mf,dc1,dc2);mfm=multiframe_mean(x1,x2,mf,dc1,dc2)
  call assert_close(mft%estimate(1),60.0_dp,tol,'multiframe constant total')
  call assert_close(mfm%estimate(1),20.0_dp,tol,'multiframe constant mean')
  call assert_close(mft%variance(1,1),1200.0_dp,tol,'multiframe total variance')
  call assert_close(mfm%variance(1,1),200.0_dp/9.0_dp,tol,'multiframe mean variance')

  ov1=reshape([2.0_dp,2.0_dp,0.0_dp,4.0_dp],[2,2]);ov2=reshape([2.0_dp,0.0_dp,4.0_dp,4.0_dp],[2,2])
  w1=2.0_dp;w2=4.0_dp
  call make_multiframe_expected(w1,w2,ov1,ov2,mfe,overlaps_are_weights=.true.)
  call assert_close(mfe%frame_weight1(1)*w1(1),2.0_dp,tol,'expected exclusive frame1')
  call assert_close(mfe%frame_weight1(2)*w1(2),4.0_dp/3.0_dp,tol,'expected overlap frame1')
  call assert_close(mfe%frame_weight2(1)*w2(1),4.0_dp/3.0_dp,tol,'expected overlap frame2')
  call assert_close(mfe%frame_weight2(2)*w2(2),4.0_dp,tol,'expected exclusive frame2')

  bw=1.0_dp;cl(:,1)=[1,2,3,4];str=1;ss=4.0_dp;ps=huge(1.0_dp)
  call make_design(bw,cl,d,strata=str,samp_size=ss,pop_size=ps)
  call random_seed(size=i);allocate(seed(i));seed=[(37*r+11,r=1,i)];call random_seed(put=seed)
  call make_mrb(d,20,rep)
  call assert_close(rep%scale,1.0_dp,tol,'MRB scale')
  do r=1,rep%r
    call assert_close(sum(rep%repweights(:,r)),4.0_dp,tol,'MRB replicate total weight')
    if(count(abs(rep%repweights(:,r)-2.0_dp)<tol)/=2.or.count(abs(rep%repweights(:,r))<tol)/=2) error stop 'MRB half-sample pattern'
  end do
  call assert_close(rep%rscales(1),1.0_dp/19.0_dp,tol,'MRB rscale')

  print '(a)','test_phase_multiframe_mrb: PASS'
contains
  subroutine assert_close(a,b,eps,msg)
    real(dp),intent(in)::a,b,eps;character(len=*),intent(in)::msg
    if(abs(a-b)>eps*max(1.0_dp,abs(a),abs(b))) then
      print '(a,2es24.14)',trim(msg)//' mismatch: ',a,b
      error stop 1
    end if
  end subroutine assert_close
end program test_phase_multiframe_mrb
