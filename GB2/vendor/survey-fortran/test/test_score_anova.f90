program test_score_anova
  use survey_kinds, only : dp
  use survey_types, only : survey_design_t, glm_result_t, model_test_t, FAMILY_GAUSSIAN, LINK_IDENTITY
  use survey_design, only : make_design
  use survey_glm, only : svy_glm
  use survey_score, only : glm_pseudoscore_test, glm_working_score_test, wald_term_test, misspecified_lrt_test
  use survey_special, only : f_survival
  implicit none
  integer,parameter :: n=40
  real(dp) :: x(n,2),y(n),w(n),ss(n,1),pops(n,1),beta0(2),m0(2,2),vr(2,2),tol
  integer :: cl(n,1),str(n,1),i
  type(survey_design_t) :: d
  type(glm_result_t) :: fit
  type(model_test_t) :: pst,wst,wt,lr
  tol=2e-8_dp
  w=1.0_dp;str=1;ss=real(n,dp);pops=huge(1.0_dp)
  do i=1,n
    cl(i,1)=i;x(i,1)=1.0_dp;x(i,2)=(real(i,dp)-20.5_dp)/10.0_dp
    y(i)=1.5_dp+0.8_dp*x(i,2)+0.25_dp*sin(1.7_dp*real(i,dp))
  end do
  call make_design(w,cl,d,strata=str,samp_size=ss,pop_size=pops)
  call svy_glm(x,y,d,fit,FAMILY_GAUSSIAN,LINK_IDENTITY)
  call wald_term_test(fit%coef,fit%vcov,[2],fit%df_residual,wt)
  if(wt%statistic<=0.0_dp.or.wt%p_value<0.0_dp.or.wt%p_value>1.0_dp) error stop 'Wald term test range'

  beta0=[sum(y)/real(n,dp),0.0_dp]
  call glm_pseudoscore_test(x,y,d,beta0,[2],pst,FAMILY_GAUSSIAN,LINK_IDENTITY)
  call glm_working_score_test(x,y,d,beta0,[2],wst,FAMILY_GAUSSIAN,LINK_IDENTITY,dispersion=sum((y-matmul(x,fit%coef))**2)/real(n-2,dp))
  if(pst%statistic<=0.0_dp.or.pst%p_value<0.0_dp.or.pst%p_value>1.0_dp.or.pst%df/=1) error stop 'pseudoscore test range'
  if(wst%statistic<=0.0_dp.or.wst%p_value<0.0_dp.or.wst%p_value>1.0_dp.or.wst%df/=1) error stop 'working score test range'
  if(pst%p_value>0.05_dp.or.wst%p_value>0.05_dp.or.wt%p_value>0.05_dp) error stop 'score/Wald should detect slope'

  m0=0.0_dp;vr=0.0_dp;m0(1,1)=1.0_dp;m0(2,2)=1.0_dp;vr(1,1)=2.0_dp;vr(2,2)=3.0_dp
  call misspecified_lrt_test(5.0_dp,m0,vr,30,lr,saddlepoint=.false.)
  if(size(lr%lambda)/=2) error stop 'LRT eigenvalue count'
  call sort2(lr%lambda)
  call assert_close(lr%lambda(1),2.0_dp,tol,'LRT lambda 1')
  call assert_close(lr%lambda(2),3.0_dp,tol,'LRT lambda 2')
  if(lr%p_value<0.0_dp.or.lr%p_value>1.0_dp) error stop 'LRT p range'

  call wald_term_test([2.0_dp],reshape([4.0_dp],[1,1]),[1],20,wt)
  call assert_close(wt%statistic,1.0_dp,tol,'scalar Wald Q')
  call assert_close(wt%p_value,f_survival(1.0_dp,1.0_dp,20.0_dp),tol,'scalar Wald F')
  print '(a)','test_score_anova: PASS'
contains
  subroutine assert_close(a,b,eps,msg)
    real(dp),intent(in)::a,b,eps;character(len=*),intent(in)::msg
    if(abs(a-b)>eps*max(1.0_dp,abs(a),abs(b))) then
      print '(a,2es24.14)',trim(msg)//' mismatch: ',a,b;error stop 1
    end if
  end subroutine assert_close
  subroutine sort2(a)
    real(dp),intent(inout)::a(:);real(dp)::t
    if(size(a)==2.and.a(1)>a(2))then;t=a(1);a(1)=a(2);a(2)=t;end if
  end subroutine sort2
end program test_score_anova
