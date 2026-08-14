module mlr_evaluate
  use mlr_kinds, only : dp
  use mlr_types, only : resample_plan, metric_summary
  use mlr_utils, only : mean_dp
  implicit none
  private
  public :: regression_predictor, classification_predictor, regression_measure, classification_measure
  public :: resample_regression, resample_classification
  abstract interface
    subroutine regression_predictor(xtrain,ytrain,xtest,pred)
      import dp
      real(dp),intent(in)::xtrain(:,:),ytrain(:),xtest(:,:)
      real(dp),allocatable,intent(out)::pred(:)
    end subroutine
    subroutine classification_predictor(xtrain,ytrain,xtest,pred)
      import dp
      real(dp),intent(in)::xtrain(:,:),xtest(:,:)
      integer,intent(in)::ytrain(:)
      integer,allocatable,intent(out)::pred(:)
    end subroutine
    real(dp) function regression_measure(truth,pred)
      import dp
      real(dp),intent(in)::truth(:),pred(:)
    end function
    real(dp) function classification_measure(truth,pred)
      import dp
      integer,intent(in)::truth(:),pred(:)
    end function
  end interface
contains
  subroutine resample_regression(x,y,plan,predictor,measure,result)
    real(dp),intent(in)::x(:,:),y(:);type(resample_plan),intent(in)::plan
    procedure(regression_predictor)::predictor;procedure(regression_measure)::measure
    type(metric_summary),intent(out)::result
    real(dp),allocatable::pred(:),xt(:,:),xv(:,:),yt(:),yv(:);integer::i
    allocate(result%values(size(plan%test)))
    do i=1,size(plan%test)
      xt=x(plan%train(i)%idx,:);yt=y(plan%train(i)%idx);xv=x(plan%test(i)%idx,:);yv=y(plan%test(i)%idx)
      call predictor(xt,yt,xv,pred);result%values(i)=measure(yv,pred)
    end do
    result%mean=mean_dp(result%values)
    if(size(result%values)>1)then
      result%sd=sqrt(sum((result%values-result%mean)**2)/real(size(result%values)-1,dp))
    else
      result%sd=0.0_dp
    end if
  end subroutine

  subroutine resample_classification(x,y,plan,predictor,measure,result)
    real(dp),intent(in)::x(:,:);integer,intent(in)::y(:);type(resample_plan),intent(in)::plan
    procedure(classification_predictor)::predictor;procedure(classification_measure)::measure
    type(metric_summary),intent(out)::result
    real(dp),allocatable::xt(:,:),xv(:,:);integer,allocatable::yt(:),yv(:),pred(:);integer::i
    allocate(result%values(size(plan%test)))
    do i=1,size(plan%test)
      xt=x(plan%train(i)%idx,:);yt=y(plan%train(i)%idx);xv=x(plan%test(i)%idx,:);yv=y(plan%test(i)%idx)
      call predictor(xt,yt,xv,pred);result%values(i)=measure(yv,pred)
    end do
    result%mean=mean_dp(result%values)
    if(size(result%values)>1)then
      result%sd=sqrt(sum((result%values-result%mean)**2)/real(size(result%values)-1,dp))
    else
      result%sd=0.0_dp
    end if
  end subroutine
end module mlr_evaluate
