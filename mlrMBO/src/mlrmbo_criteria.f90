module mlrmbo_criteria
  use mlrmbo_kinds, only : dp
  use mlrmbo_types, only : mbo_control, mbo_path, crit_mean, crit_se, crit_ei, crit_cb, crit_aei, crit_eqi, crit_adacb
  use mlrmbo_math, only : normal_pdf, normal_cdf, normal_quantile, variance_value
  use mlrmbo_gp, only : mbo_surrogates, predict_surrogate
  implicit none
  private
  public :: eval_single_criterion, expected_improvement
contains
  elemental real(dp) function expected_improvement(mu,se,best,se_threshold) result(ei)
    real(dp), intent(in) :: mu,se,best,se_threshold
    real(dp) :: d,z
    if(se<se_threshold) then
      ei=0.0_dp
    else
      d=best-mu; z=d/se
      ei=d*normal_cdf(z)+se*normal_pdf(z)
    end if
  end function expected_improvement

  subroutine eval_single_criterion(x,sur,path,control,iter,values,criterion,lambda_override)
    real(dp), intent(in) :: x(:,:)
    type(mbo_surrogates), intent(in) :: sur
    type(mbo_path), intent(in) :: path
    type(mbo_control), intent(in) :: control
    integer, intent(in) :: iter
    real(dp), allocatable, intent(out) :: values(:)
    integer, intent(in), optional :: criterion
    real(dp), intent(in), optional :: lambda_override
    real(dp), allocatable :: mu(:),se(:),train_mu(:),train_se(:)
    real(dp) :: best,d,z,tau,qmin,mq,sq,lambda,progress,sgn
    integer :: i,id
    id=control%infill_criterion; if(present(criterion)) id=criterion
    if(control%minimize(1)) then; sgn=1.0_dp; else; sgn=-1.0_dp; end if
    call predict_surrogate(sur%model(1),x,mu,se)
    mu=sgn*mu
    allocate(values(size(mu)))
    select case(id)
    case(crit_mean)
      values=mu
    case(crit_se)
      values=-se
    case(crit_ei)
      best=minval(sgn*path%y(:,1))
      do i=1,size(mu); values(i)=-expected_improvement(mu(i),se(i),best,control%se_threshold); end do
    case(crit_cb)
      lambda=control%cb_lambda; if(present(lambda_override)) lambda=lambda_override
      values=mu-lambda*se
    case(crit_adacb)
      if(control%max_iter>0) then
        progress=min(1.0_dp,max(0.0_dp,real(iter,dp)/real(control%max_iter,dp)))
      else
        progress=0.0_dp
      end if
      lambda=(1.0_dp-progress)*control%cb_lambda_start+progress*control%cb_lambda_end
      values=mu-lambda*se
    case(crit_aei)
      call predict_surrogate(sur%model(1),path%x,train_mu,train_se)
      train_mu=sgn*train_mu
      i=minloc(train_mu+train_se,dim=1); best=train_mu(i)
      tau=sqrt(max(variance_value(path%y(:,1)-sgn*train_mu),1.0e-12_dp))
      if(control%aei_use_nugget .and. sur%model(1)%covariance%nugget_flag) &
        tau=sqrt(max(sur%model(1)%covariance%nugget,1.0e-12_dp))
      do i=1,size(mu)
        if(se(i)<control%se_threshold) then
          values(i)=0.0_dp
        else
          d=best-mu(i); z=d/se(i)
          values(i)=-(d*normal_cdf(z)+se(i)*normal_pdf(z))* &
            (1.0_dp-tau/sqrt(tau*tau+se(i)*se(i)))
        end if
      end do
    case(crit_eqi)
      call predict_surrogate(sur%model(1),path%x,train_mu,train_se)
      train_mu=sgn*train_mu
      tau=sqrt(max(variance_value(path%y(:,1)-sgn*train_mu),1.0e-12_dp))
      qmin=minval(train_mu+normal_quantile(control%eqi_beta)*train_se)
      do i=1,size(mu)
        if(se(i)<control%se_threshold) then
          values(i)=0.0_dp
        else
          mq=mu(i)+normal_quantile(control%eqi_beta)*sqrt((tau*se(i)*se(i))/(tau+se(i)*se(i)))
          sq=se(i)*se(i)/sqrt(tau*tau+se(i)*se(i))
          if(sq<=1.0e-15_dp) then
            values(i)=0.0_dp
          else
            d=qmin-mq; z=d/sq
            values(i)=-sq*(z*normal_cdf(z)+normal_pdf(z))
          end if
        end if
      end do
    case default
      error stop 'eval_single_criterion: unsupported criterion'
    end select
  end subroutine eval_single_criterion
end module mlrmbo_criteria
