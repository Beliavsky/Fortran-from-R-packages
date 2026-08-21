! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module survey_mle
  use survey_kinds, only : dp
  use survey_types, only : survey_design_t, mle_result_t
  use survey_taylor, only : svyrecvar
  use survey_linalg, only : sym_pinv
  use minqa_module, only : minqa_result_t, minqa_control_t, uobyqa, bobyqa
  use numderiv, only : hessian, nd_success
  implicit none
  private
  public :: svy_mle, observation_loglike, observation_score

  abstract interface
    function observation_loglike(theta,i) result(value)
      import dp
      real(dp), intent(in) :: theta(:)
      integer, intent(in) :: i
      real(dp) :: value
    end function observation_loglike
    subroutine observation_score(theta,i,score)
      import dp
      real(dp), intent(in) :: theta(:)
      integer, intent(in) :: i
      real(dp), intent(out) :: score(:)
    end subroutine observation_score
  end interface

  procedure(observation_loglike), pointer :: active_loglike => null()
  real(dp), allocatable :: active_weight(:)
contains
  subroutine svy_mle(loglike,design,start,result,score,lower,upper,maxfun)
    procedure(observation_loglike) :: loglike
    type(survey_design_t),intent(in)::design
    real(dp),intent(in)::start(:)
    type(mle_result_t),intent(out)::result
    procedure(observation_score),optional::score
    real(dp),intent(in),optional::lower(:),upper(:)
    integer,intent(in),optional::maxfun
    type(minqa_result_t)::opt
    type(minqa_control_t)::ctrl
    real(dp),allocatable::theta(:),hh(:,:),infoinv(:,:),srow(:,:),infl(:,:),lo(:),hi(:)
    integer::p,i,rank,istat
    character(:),allocatable::msg
    p=size(start);if(p<1)error stop 'svy_mle: empty parameter vector'
    if(associated(active_loglike)) error stop 'svy_mle: nested/concurrent calls are not supported'
    active_loglike=>loglike;allocate(active_weight,source=design%weight)
    theta=start;ctrl%maxfun=10000;if(present(maxfun))ctrl%maxfun=maxfun
    if(p==1)then
      call optimize_one(theta(1),opt,lower,upper,ctrl%maxfun)
    else if(present(lower).or.present(upper))then
      allocate(lo(p),hi(p));lo=-huge(1.0_dp)/4.0_dp;hi=huge(1.0_dp)/4.0_dp
      if(present(lower))then;if(size(lower)/=p)error stop 'svy_mle: lower size mismatch';lo=lower;end if
      if(present(upper))then;if(size(upper)/=p)error stop 'svy_mle: upper size mismatch';hi=upper;end if
      call bobyqa(negloglik_active,theta,opt,lo,hi,ctrl)
    else
      call uobyqa(negloglik_active,theta,opt,ctrl)
    end if
    theta=opt%x
    call hessian(total_loglik_active,theta,hh,status=istat,message=msg)
    if(istat/=nd_success) then
      call clear_active(); error stop 'svy_mle: numerical Hessian failed'
    end if
    allocate(infoinv(p,p));call sym_pinv(-hh,infoinv,rank,info=istat)
    allocate(result%par(p),result%model_vcov(p,p),result%vcov(p,p));result%par=theta;result%model_vcov=infoinv
    result%loglik=total_loglik_active(theta);result%evaluations=opt%evaluations;result%status=opt%status;result%converged=(opt%status==0)
    if(present(score))then
      allocate(srow(design%n,p),infl(design%n,p))
      do i=1,design%n;call score(theta,i,srow(i,:));srow(i,:)=srow(i,:)*design%weight(i);end do
      infl=matmul(srow,infoinv);result%vcov=svyrecvar(infl,design);allocate(result%influence,source=infl)
    else
      result%vcov=infoinv
    end if
    call clear_active()
  end subroutine svy_mle

  real(dp) function total_loglik_active(x) result(v)
    real(dp),intent(in)::x(:);integer::i
    if(.not.associated(active_loglike).or..not.allocated(active_weight)) error stop 'survey_mle: no active likelihood'
    v=0.0_dp;do i=1,size(active_weight);v=v+active_weight(i)*active_loglike(x,i);end do
  end function total_loglik_active

  real(dp) function negloglik_active(x) result(v)
    real(dp),intent(in)::x(:);v=-total_loglik_active(x)
  end function negloglik_active

  subroutine optimize_one(x,opt,lower,upper,maxfun)
    real(dp),intent(inout)::x
    type(minqa_result_t),intent(out)::opt
    real(dp),intent(in),optional::lower(:),upper(:)
    integer,intent(in)::maxfun
    real(dp)::h,g,curv,step,oldf,newf,trial
    integer::it
    oldf=negloglik_active([x]);opt%status=1;opt%evaluations=1
    do it=1,max(1,maxfun/5)
      h=1.0e-4_dp*max(1.0_dp,abs(x))
      g=(negloglik_active([x+h])-negloglik_active([x-h]))/(2*h)
      curv=(negloglik_active([x+h])-2*oldf+negloglik_active([x-h]))/(h*h);opt%evaluations=opt%evaluations+4
      if(curv<=tiny(1.0_dp))then;step=-sign(min(1.0_dp,abs(g)),g);else;step=-g/curv;end if
      step=max(-2.0_dp*max(1.0_dp,abs(x)),min(2.0_dp*max(1.0_dp,abs(x)),step));trial=x+step
      if(present(lower))then;if(size(lower)/=1)error stop 'svy_mle: lower size mismatch';trial=max(trial,lower(1));end if
      if(present(upper))then;if(size(upper)/=1)error stop 'svy_mle: upper size mismatch';trial=min(trial,upper(1));end if
      newf=negloglik_active([trial]);opt%evaluations=opt%evaluations+1
      do while(newf>oldf.and.abs(step)>1e-10_dp)
        step=step/2;trial=x+step
        if(present(lower))trial=max(trial,lower(1));if(present(upper))trial=min(trial,upper(1))
        newf=negloglik_active([trial]);opt%evaluations=opt%evaluations+1
      end do
      x=trial;if(abs(newf-oldf)<=1e-10_dp*(1+abs(oldf)))then;opt%status=0;exit;end if;oldf=newf
    end do
    allocate(opt%x(1));opt%x(1)=x;opt%fval=oldf;opt%raw_status=opt%status;opt%message='one-dimensional Newton/line-search'
  end subroutine optimize_one

  subroutine clear_active()
    nullify(active_loglike);if(allocated(active_weight))deallocate(active_weight)
  end subroutine clear_active
end module survey_mle
