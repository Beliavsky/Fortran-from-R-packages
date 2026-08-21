program test_ordinal_loglinear
  use survey, only : dp, survey_design_t, olr_result_t, loglin_result_t, loglin_test_t, &
    make_design, svy_olr, olr_predict_proba, svy_loglin_prob, svy_loglin_compare, &
    OLR_LOGISTIC, OLR_CAUCHIT, weighted_chisq_survival, chisq_survival
  implicit none
  integer, parameter :: n=1600
  real(dp) :: x(n,2), w(n), beta(2), zeta(2), eta, u, c1, c2
  integer :: y(n), cl(n,1), i, m
  type(survey_design_t) :: d
  type(olr_result_t) :: fit
  real(dp) :: pp(5,3), xx(5,2)
  real(dp) :: prob(4), vc(4,4), xind(4,3), xsat(4,4), lam(3), p0, p1
  type(loglin_result_t) :: lind, lsat
  type(loglin_test_t) :: lt

  beta=[0.7_dp,-0.4_dp]; zeta=[-0.5_dp,0.8_dp]
  do i=1,n
    x(i,1)=sin(0.37_dp*real(i,dp));x(i,2)=cos(0.13_dp*real(i,dp));w(i)=1.0_dp;cl(i,1)=i
    eta=dot_product(x(i,:),beta);c1=logistic(zeta(1)-eta);c2=logistic(zeta(2)-eta)
    u=modulo(real(i,dp)*0.6180339887498949_dp,1.0_dp)
    if(u<c1)then;y(i)=1;else if(u<c2)then;y(i)=2;else;y(i)=3;end if
  end do
  call make_design(w,cl,d)
  call svy_olr(x,y,d,fit,OLR_LOGISTIC,maxfun=12000)
  call assert_true(fit%converged,'ordinal convergence')
  call assert_close(maxval(abs(fit%coef-beta)),0.0_dp,0.12_dp,'ordinal beta recovery')
  call assert_close(maxval(abs(fit%zeta-zeta)),0.0_dp,0.12_dp,'ordinal cutpoint recovery')
  call assert_true(all([(fit%zeta(i)>fit%zeta(i-1),i=2,size(fit%zeta))]),'ordered cutpoints')
  xx=reshape([0.0_dp,0.0_dp, 0.4_dp,-0.2_dp, -0.3_dp,0.7_dp, 0.8_dp,0.3_dp, -0.9_dp,-0.1_dp],[5,2],order=[2,1])
  do m=OLR_LOGISTIC,OLR_CAUCHIT
    call olr_predict_proba(xx,beta,zeta,m,pp)
    do i=1,5
      call assert_close(sum(pp(i,:)),1.0_dp,2e-13_dp,'ordinal probability sum')
      call assert_true(all(pp(i,:)>=0.0_dp.and.pp(i,:)<=1.0_dp),'ordinal probability range')
    end do
  end do

  prob=[0.32_dp,0.18_dp,0.28_dp,0.22_dp]
  vc=-spread(prob,2,4)*spread(prob,1,4)/1000.0_dp
  do i=1,4;vc(i,i)=vc(i,i)+prob(i)/1000.0_dp;end do
  ! Cells are (A,B)=(0,0),(1,0),(0,1),(1,1).
  xind=reshape([1.0_dp,0.0_dp,0.0_dp, 1.0_dp,1.0_dp,0.0_dp, 1.0_dp,0.0_dp,1.0_dp, 1.0_dp,1.0_dp,1.0_dp],[4,3],order=[2,1])
  xsat(:,1:3)=xind;xsat(:,4)=[0.0_dp,0.0_dp,0.0_dp,1.0_dp]
  call svy_loglin_prob(prob,vc,xind,1000.0_dp,999,lind)
  call svy_loglin_prob(prob,vc,xsat,1000.0_dp,999,lsat)
  call assert_true(lind%converged.and.lsat%converged,'loglinear convergence')
  call assert_close(maxval(abs(lsat%fitted_prob-prob)),0.0_dp,2e-7_dp,'saturated fit')
  call assert_true(lsat%deviance<1e-8_dp,'saturated deviance')
  call svy_loglin_compare(prob,vc,xind,xsat,1000.0_dp,999,lt)
  call assert_true(lt%deviance>=0.0_dp.and.lt%score>=0.0_dp,'loglinear test statistics')
  call assert_true(lt%p_deviance_saddle>=0.0_dp.and.lt%p_deviance_saddle<=1.0_dp,'loglinear saddle p')
  call assert_true(size(lt%lambda)==1,'loglinear one-df comparison')

  lam=1.0_dp;p0=weighted_chisq_survival(5.0_dp,lam,.false.);p1=chisq_survival(5.0_dp,3.0_dp)
  call assert_close(p0,p1,2e-13_dp,'weighted chi-square equal-lambda reduction')
  print '(a)','test_ordinal_loglinear: PASS'
contains
  pure real(dp) function logistic(z) result(p)
    real(dp),intent(in)::z
    p=1.0_dp/(1.0_dp+exp(-z))
  end function logistic
  subroutine assert_true(ok,msg)
    logical,intent(in)::ok;character(*),intent(in)::msg
    if(.not.ok)then;write(*,'(a)')'FAIL: '//trim(msg);error stop 1;end if
  end subroutine assert_true
  subroutine assert_close(a,b,tol,msg)
    real(dp),intent(in)::a,b,tol;character(*),intent(in)::msg
    if(abs(a-b)>tol)then;write(*,'(a,2es24.14)')'FAIL: '//trim(msg)//' ',a,b;error stop 1;end if
  end subroutine assert_close
end program test_ordinal_loglinear
