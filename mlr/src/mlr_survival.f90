module mlr_survival
  use mlr_kinds, only : dp
  use survival_types, only : coxph_result, concordance_result, survfit_result
  use survival_cox, only : coxph_fit
  use survival_stats, only : concordance_right
  use survival_nonparametric, only : kaplan_meier
  implicit none
  private
  public :: mlr_cox_model, fit_cox_learner, predict_cox_risk, measure_cindex, fit_kaplan_meier
  type :: mlr_cox_model
    type(coxph_result) :: fit
  end type mlr_cox_model
contains
  subroutine fit_cox_learner(time,status,x,model,method,weights)
    real(dp),intent(in)::time(:),x(:,:);integer,intent(in)::status(:);type(mlr_cox_model),intent(out)::model
    character(len=*),intent(in),optional::method;real(dp),intent(in),optional::weights(:)
    if(present(method))then
      if(present(weights))then
        call coxph_fit(time,status,x,model%fit,method=method,weights=weights)
      else
        call coxph_fit(time,status,x,model%fit,method=method)
      end if
    else
      if(present(weights))then
        call coxph_fit(time,status,x,model%fit,weights=weights)
      else
        call coxph_fit(time,status,x,model%fit)
      end if
    end if
  end subroutine

  subroutine predict_cox_risk(model,x,risk)
    type(mlr_cox_model),intent(in)::model;real(dp),intent(in)::x(:,:);real(dp),allocatable,intent(out)::risk(:)
    allocate(risk(size(x,1)));risk=matmul(x,model%fit%coef)
  end subroutine

  real(dp) function measure_cindex(time,status,risk) result(v)
    real(dp),intent(in)::time(:),risk(:);integer,intent(in)::status(:);type(concordance_result)::r
    call concordance_right(time,status,risk,r);v=r%cindex
  end function

  subroutine fit_kaplan_meier(time,status,fit,weights)
    real(dp),intent(in)::time(:);integer,intent(in)::status(:);type(survfit_result),intent(out)::fit
    real(dp),intent(in),optional::weights(:)
    if(present(weights))then;call kaplan_meier(time,status,fit,weights);else;call kaplan_meier(time,status,fit);end if
  end subroutine
end module mlr_survival
