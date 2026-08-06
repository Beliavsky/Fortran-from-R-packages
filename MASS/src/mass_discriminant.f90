! SPDX-License-Identifier: GPL-3.0-only
module mass_discriminant
  use rrcov_kinds, only : dp
  use rrcov_types, only : lda_model, qda_model, rrcov_success
  use rrcov_da, only : lda_classic_fit, lda_cov_fit, lda_predict_rr => lda_predict, &
    qda_classic_fit, qda_cov_fit, qda_predict_rr => qda_predict, confusion_matrix
  use mass_types, only : mass_success, mass_invalid_argument
  implicit none
  private
  public :: lda_model, qda_model
  public :: lda_fit, qda_fit, lda_predict, qda_predict, confusion_matrix
contains

  subroutine lda_fit(x,grouping,model,priors,method,nu,nsamp,seed)
    real(dp),intent(in)::x(:,:)
    integer,intent(in)::grouping(:)
    type(lda_model),intent(out)::model
    real(dp),intent(in),optional::priors(:),nu
    character(len=*),intent(in),optional::method
    integer,intent(in),optional::nsamp,seed
    character(len=16)::m
    m="moment";if(present(method))m=method
    if(present(nu))then
      if(nu<=0.0_dp)then;model%status=1;return;end if
    end if
    select case(trim(m))
    case("mve")
      call lda_cov_fit(x,grouping,model,"mve",priors=priors,nsamp=nsamp,seed=seed)
      model%method="MASS LDA (MVE covariance)"
    case("t")
      call lda_cov_fit(x,grouping,model,"mest",priors=priors,nsamp=nsamp,seed=seed)
      model%method="MASS LDA (robust t covariance)"
    case default
      call lda_classic_fit(x,grouping,model,priors)
      if(trim(m)=="mle") model%method="MASS LDA (maximum likelihood scaling)"
    end select
  end subroutine lda_fit

  subroutine qda_fit(x,grouping,model,priors,method,nsamp,seed)
    real(dp),intent(in)::x(:,:)
    integer,intent(in)::grouping(:)
    type(qda_model),intent(out)::model
    real(dp),intent(in),optional::priors(:)
    character(len=*),intent(in),optional::method
    integer,intent(in),optional::nsamp,seed
    character(len=16)::m
    m="moment";if(present(method))m=method
    if(trim(m)=="mve")then
      call qda_cov_fit(x,grouping,model,"mve",priors=priors,nsamp=nsamp,seed=seed)
      model%method="MASS QDA (MVE covariance)"
    else
      call qda_classic_fit(x,grouping,model,priors)
      if(trim(m)=="mle")model%method="MASS QDA (maximum likelihood scaling)"
    end if
  end subroutine qda_fit

  subroutine lda_predict(model,x,predicted,posterior,scores,status)
    type(lda_model),intent(in)::model
    real(dp),intent(in)::x(:,:)
    integer,allocatable,intent(out)::predicted(:)
    real(dp),allocatable,intent(out),optional::posterior(:,:),scores(:,:)
    integer,intent(out),optional::status
    integer::st
    call lda_predict_rr(model,x,predicted,posterior,scores,st)
    if(present(status))status=merge(mass_success,mass_invalid_argument,st==rrcov_success)
  end subroutine lda_predict

  subroutine qda_predict(model,x,predicted,posterior,scores,status)
    type(qda_model),intent(in)::model
    real(dp),intent(in)::x(:,:)
    integer,allocatable,intent(out)::predicted(:)
    real(dp),allocatable,intent(out),optional::posterior(:,:),scores(:,:)
    integer,intent(out),optional::status
    integer::st
    call qda_predict_rr(model,x,predicted,posterior,scores,st)
    if(present(status))status=merge(mass_success,mass_invalid_argument,st==rrcov_success)
  end subroutine qda_predict

end module mass_discriminant
