module mlr_metrics
  use mlr_kinds, only : dp
  use mlr_utils, only : mean_dp, median_dp, ranks_average
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  implicit none
  private
  public :: measure_sse, measure_mse, measure_rmse, measure_medse, measure_sae, measure_mae, measure_medae
  public :: measure_rsq, measure_expvar, measure_rrse, measure_rae, measure_mape, measure_msle, measure_rmsle
  public :: measure_spearman, measure_kendall
  public :: confusion_matrix, measure_mmce, measure_acc, measure_ber, measure_bac, measure_kappa, measure_wkappa
  public :: measure_tp, measure_tn, measure_fp, measure_fn, measure_tpr, measure_tnr, measure_fpr, measure_fnr
  public :: measure_ppv, measure_npv, measure_fdr, measure_mcc, measure_f1, measure_gmean, measure_auc
  public :: measure_brier, measure_brier_scaled, measure_logloss, measure_multiclass_brier
  public :: aggregate_mean, aggregate_sd, aggregate_rmse, aggregate_b632
contains
  pure real(dp) function measure_sse(truth, response) result(v)
    real(dp), intent(in) :: truth(:), response(:)
    v = sum((response-truth)**2)
  end function
  pure real(dp) function measure_mse(truth, response) result(v)
    real(dp), intent(in) :: truth(:), response(:)
    v = measure_sse(truth,response)/real(size(truth),dp)
  end function
  pure real(dp) function measure_rmse(truth, response) result(v)
    real(dp), intent(in) :: truth(:), response(:)
    v = sqrt(measure_mse(truth,response))
  end function
  real(dp) function measure_medse(truth, response) result(v)
    real(dp), intent(in) :: truth(:), response(:)
    v = median_dp((response-truth)**2)
  end function
  pure real(dp) function measure_sae(truth, response) result(v)
    real(dp), intent(in) :: truth(:), response(:)
    v = sum(abs(response-truth))
  end function
  pure real(dp) function measure_mae(truth, response) result(v)
    real(dp), intent(in) :: truth(:), response(:)
    v = measure_sae(truth,response)/real(size(truth),dp)
  end function
  real(dp) function measure_medae(truth, response) result(v)
    real(dp), intent(in) :: truth(:), response(:)
    v = median_dp(abs(response-truth))
  end function
  pure real(dp) function measure_rsq(truth, response) result(v)
    real(dp), intent(in) :: truth(:), response(:)
    real(dp) :: tss, m
    m = mean_dp(truth); tss = sum((truth-m)**2)
    if (tss <= tiny(1.0_dp)) then
      v = ieee_value(0.0_dp, ieee_quiet_nan)
    else
      v = 1.0_dp - measure_sse(truth,response)/tss
    end if
  end function
  pure real(dp) function measure_expvar(truth, response) result(v)
    real(dp), intent(in) :: truth(:), response(:)
    real(dp) :: tss, m
    m = mean_dp(truth); tss = sum((truth-m)**2)
    if (tss <= tiny(1.0_dp)) then
      v = ieee_value(0.0_dp, ieee_quiet_nan)
    else
      v = sum((response-m)**2)/tss
    end if
  end function
  pure real(dp) function measure_rrse(truth, response) result(v)
    real(dp), intent(in) :: truth(:), response(:)
    real(dp) :: tss
    tss = sum((truth-mean_dp(truth))**2)
    if (tss <= tiny(1.0_dp)) then
      v = ieee_value(0.0_dp, ieee_quiet_nan)
    else
      v = sqrt(measure_sse(truth,response)/tss)
    end if
  end function
  pure real(dp) function measure_rae(truth, response) result(v)
    real(dp), intent(in) :: truth(:), response(:)
    real(dp) :: d
    d = sum(abs(truth-mean_dp(truth)))
    if (d <= tiny(1.0_dp)) then
      v = ieee_value(0.0_dp, ieee_quiet_nan)
    else
      v = measure_sae(truth,response)/d
    end if
  end function
  pure real(dp) function measure_mape(truth, response) result(v)
    real(dp), intent(in) :: truth(:), response(:)
    integer :: i
    v = 0.0_dp
    do i=1,size(truth)
      if (abs(truth(i)) <= tiny(1.0_dp)) then
        v = ieee_value(0.0_dp, ieee_quiet_nan); return
      end if
      v = v + abs((truth(i)-response(i))/truth(i))
    end do
    v = v/real(size(truth),dp)
  end function
  pure real(dp) function measure_msle(truth, response) result(v)
    real(dp), intent(in) :: truth(:), response(:)
    if (any(truth < -1.0_dp) .or. any(response < -1.0_dp)) then
      v = ieee_value(0.0_dp, ieee_quiet_nan)
    else
      v = sum((log(response+1.0_dp)-log(truth+1.0_dp))**2)/real(size(truth),dp)
    end if
  end function
  pure real(dp) function measure_rmsle(truth, response) result(v)
    real(dp), intent(in) :: truth(:), response(:)
    v = sqrt(measure_msle(truth,response))
  end function

  subroutine confusion_matrix(truth, response, nclass, cm)
    integer, intent(in) :: truth(:), response(:), nclass
    integer, allocatable, intent(out) :: cm(:,:)
    integer :: i
    allocate(cm(nclass,nclass)); cm=0
    do i=1,size(truth)
      if (truth(i)>=1 .and. truth(i)<=nclass .and. response(i)>=1 .and. response(i)<=nclass) &
        cm(truth(i),response(i)) = cm(truth(i),response(i)) + 1
    end do
  end subroutine

  real(dp) function measure_mmce(truth,response) result(v)
    integer,intent(in)::truth(:),response(:)
    v=real(count(truth/=response),dp)/real(size(truth),dp)
  end function
  real(dp) function measure_acc(truth,response) result(v)
    integer,intent(in)::truth(:),response(:)
    v=1.0_dp-measure_mmce(truth,response)
  end function
  real(dp) function measure_bac(truth,response,nclass) result(v)
    integer,intent(in)::truth(:),response(:),nclass
    integer,allocatable::cm(:,:); integer::k,den
    call confusion_matrix(truth,response,nclass,cm); v=0.0_dp
    do k=1,nclass
      den=sum(cm(k,:)); if(den>0)v=v+real(cm(k,k),dp)/real(den,dp)
    end do
    v=v/real(nclass,dp)
  end function
  real(dp) function measure_ber(truth,response,nclass) result(v)
    integer,intent(in)::truth(:),response(:),nclass
    v=1.0_dp-measure_bac(truth,response,nclass)
  end function
  real(dp) function measure_kappa(truth,response,nclass) result(v)
    integer,intent(in)::truth(:),response(:),nclass
    integer,allocatable::cm(:,:); integer::n,k
    real(dp)::p0,pe
    call confusion_matrix(truth,response,nclass,cm); n=sum(cm)
    p0=0.0_dp; pe=0.0_dp
    do k=1,nclass
      p0=p0+real(cm(k,k),dp)
      pe=pe+real(sum(cm(k,:))*sum(cm(:,k)),dp)
    end do
    p0=p0/real(n,dp); pe=pe/real(n*n,dp)
    if(abs(1.0_dp-pe)<=tiny(1.0_dp))then
      v=ieee_value(0.0_dp,ieee_quiet_nan)
    else
      v=(p0-pe)/(1.0_dp-pe)
    end if
  end function
  real(dp) function measure_wkappa(truth,response,nclass) result(v)
    integer,intent(in)::truth(:),response(:),nclass
    integer,allocatable::cm(:,:); integer::i,j,n
    real(dp)::num,den,w,e
    call confusion_matrix(truth,response,nclass,cm); n=sum(cm); num=0.0_dp; den=0.0_dp
    do i=1,nclass; do j=1,nclass
      w=real((i-j)*(i-j),dp)
      num=num+w*real(cm(i,j),dp)/real(n,dp)
      e=real(sum(cm(i,:))*sum(cm(:,j)),dp)/real(n*n,dp)
      den=den+w*e
    end do; end do
    if(den<=tiny(1.0_dp))then; v=1.0_dp; else; v=1.0_dp-num/den; end if
  end function

  integer function measure_tp(truth,response,positive) result(v)
    integer,intent(in)::truth(:),response(:),positive
    v=count(truth==positive .and. response==positive)
  end function
  integer function measure_tn(truth,response,negative) result(v)
    integer,intent(in)::truth(:),response(:),negative
    v=count(truth==negative .and. response==negative)
  end function
  integer function measure_fp(truth,response,positive) result(v)
    integer,intent(in)::truth(:),response(:),positive
    v=count(truth/=response .and. response==positive)
  end function
  integer function measure_fn(truth,response,negative) result(v)
    integer,intent(in)::truth(:),response(:),negative
    v=count(truth/=response .and. response==negative)
  end function
  real(dp) function measure_tpr(truth,response,positive) result(v)
    integer,intent(in)::truth(:),response(:),positive
    v=real(measure_tp(truth,response,positive),dp)/real(count(truth==positive),dp)
  end function
  real(dp) function measure_tnr(truth,response,negative) result(v)
    integer,intent(in)::truth(:),response(:),negative
    v=real(measure_tn(truth,response,negative),dp)/real(count(truth==negative),dp)
  end function
  real(dp) function measure_fpr(truth,response,negative,positive) result(v)
    integer,intent(in)::truth(:),response(:),negative,positive
    v=real(measure_fp(truth,response,positive),dp)/real(count(truth==negative),dp)
  end function
  real(dp) function measure_fnr(truth,response,negative,positive) result(v)
    integer,intent(in)::truth(:),response(:),negative,positive
    v=real(measure_fn(truth,response,negative),dp)/real(count(truth==positive),dp)
  end function
  real(dp) function measure_ppv(truth,response,positive) result(v)
    integer,intent(in)::truth(:),response(:),positive
    integer::d; d=count(response==positive)
    if(d==0)then; v=0.0_dp; else; v=real(measure_tp(truth,response,positive),dp)/real(d,dp); end if
  end function
  real(dp) function measure_npv(truth,response,negative) result(v)
    integer,intent(in)::truth(:),response(:),negative
    integer::d; d=count(response==negative)
    if(d==0)then; v=0.0_dp; else; v=real(measure_tn(truth,response,negative),dp)/real(d,dp); end if
  end function
  real(dp) function measure_fdr(truth,response,positive) result(v)
    integer,intent(in)::truth(:),response(:),positive
    v=1.0_dp-measure_ppv(truth,response,positive)
  end function
  real(dp) function measure_mcc(truth,response,negative,positive) result(v)
    integer,intent(in)::truth(:),response(:),negative,positive
    real(dp)::tp,tn,fp,fn,den
    tp=real(measure_tp(truth,response,positive),dp); tn=real(measure_tn(truth,response,negative),dp)
    fp=real(measure_fp(truth,response,positive),dp); fn=real(measure_fn(truth,response,negative),dp)
    den=sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    if(den<=tiny(1.0_dp))then; v=0.0_dp; else; v=(tp*tn-fp*fn)/den; end if
  end function
  real(dp) function measure_f1(truth,response,positive) result(v)
    integer,intent(in)::truth(:),response(:),positive
    real(dp)::p,r
    p=measure_ppv(truth,response,positive); r=measure_tpr(truth,response,positive)
    if(p+r<=tiny(1.0_dp))then; v=0.0_dp; else; v=2.0_dp*p*r/(p+r); end if
  end function
  real(dp) function measure_gmean(truth,response,negative,positive) result(v)
    integer,intent(in)::truth(:),response(:),negative,positive
    v=sqrt(measure_tpr(truth,response,positive)*measure_tnr(truth,response,negative))
  end function

  real(dp) function measure_auc(probabilities, truth, positive) result(v)
    real(dp),intent(in)::probabilities(:)
    integer,intent(in)::truth(:),positive
    real(dp),allocatable::r(:); integer::npos,nneg,i
    allocate(r(size(probabilities))); call ranks_average(probabilities,r)
    npos=count(truth==positive); nneg=size(truth)-npos
    if(npos==0 .or. nneg==0)then
      v=ieee_value(0.0_dp,ieee_quiet_nan); return
    end if
    v=0.0_dp
    do i=1,size(truth); if(truth(i)==positive)v=v+r(i); end do
    v=(v-real(npos*(npos+1),dp)/2.0_dp)/real(npos*nneg,dp)
  end function
  real(dp) function measure_brier(probabilities,truth,positive) result(v)
    real(dp),intent(in)::probabilities(:); integer,intent(in)::truth(:),positive
    integer::i; v=0.0_dp
    do i=1,size(truth); v=v+(merge(1.0_dp,0.0_dp,truth(i)==positive)-probabilities(i))**2; end do
    v=v/real(size(truth),dp)
  end function
  real(dp) function measure_brier_scaled(probabilities,truth,positive) result(v)
    real(dp),intent(in)::probabilities(:); integer,intent(in)::truth(:),positive
    real(dp)::inc,bmax
    inc=mean_dp(probabilities); bmax=inc*(1.0_dp-inc)**2+(1.0_dp-inc)*inc**2
    if(bmax<=tiny(1.0_dp))then; v=ieee_value(0.0_dp,ieee_quiet_nan)
    else; v=1.0_dp-measure_brier(probabilities,truth,positive)/bmax; end if
  end function
  real(dp) function measure_logloss(prob,truth) result(v)
    real(dp),intent(in)::prob(:,:); integer,intent(in)::truth(:)
    integer::i; real(dp)::p
    v=0.0_dp
    do i=1,size(truth)
      p=max(min(prob(i,truth(i)),1.0_dp-epsilon(1.0_dp)),epsilon(1.0_dp)); v=v-log(p)
    end do
    v=v/real(size(truth),dp)
  end function
  real(dp) function measure_multiclass_brier(prob,truth) result(v)
    real(dp),intent(in)::prob(:,:); integer,intent(in)::truth(:)
    integer::i,k; real(dp)::y
    v=0.0_dp
    do i=1,size(truth); do k=1,size(prob,2)
      y=merge(1.0_dp,0.0_dp,truth(i)==k); v=v+(prob(i,k)-y)**2
    end do; end do
    v=v/real(size(truth),dp)
  end function

  real(dp) function measure_spearman(truth,response) result(v)
    real(dp),intent(in)::truth(:),response(:)
    real(dp),allocatable::a(:),b(:); real(dp)::ma,mb,den
    allocate(a(size(truth)),b(size(truth))); call ranks_average(truth,a); call ranks_average(response,b)
    ma=mean_dp(a); mb=mean_dp(b); den=sqrt(sum((a-ma)**2)*sum((b-mb)**2))
    if(den<=tiny(1.0_dp))then; v=ieee_value(0.0_dp,ieee_quiet_nan); else; v=sum((a-ma)*(b-mb))/den; end if
  end function
  real(dp) function measure_kendall(truth,response) result(v)
    real(dp),intent(in)::truth(:),response(:)
    integer::i,j,c,d,tx,ty; real(dp)::den
    c=0;d=0;tx=0;ty=0
    do i=1,size(truth)-1; do j=i+1,size(truth)
      if(abs(truth(i)-truth(j))<=32.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(truth(i)),abs(truth(j))))tx=tx+1
      if(abs(response(i)-response(j))<=32.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(response(i)),abs(response(j))))ty=ty+1
      if((truth(i)-truth(j))*(response(i)-response(j))>0.0_dp)c=c+1
      if((truth(i)-truth(j))*(response(i)-response(j))<0.0_dp)d=d+1
    end do; end do
    den=sqrt(real(c+d+tx,dp)*real(c+d+ty,dp))
    if(den<=tiny(1.0_dp))then; v=ieee_value(0.0_dp,ieee_quiet_nan); else; v=real(c-d,dp)/den; end if
  end function

  pure real(dp) function aggregate_mean(x) result(v)
    real(dp),intent(in)::x(:); v=mean_dp(x)
  end function
  pure real(dp) function aggregate_sd(x) result(v)
    real(dp),intent(in)::x(:)
    if(size(x)<2)then;v=0.0_dp;else;v=sqrt(sum((x-mean_dp(x))**2)/real(size(x)-1,dp));end if
  end function
  pure real(dp) function aggregate_rmse(x) result(v)
    real(dp),intent(in)::x(:); v=sqrt(sum(x*x)/real(size(x),dp))
  end function
  pure real(dp) function aggregate_b632(test_perf,train_perf) result(v)
    real(dp),intent(in)::test_perf(:),train_perf(:)
    v=mean_dp(0.632_dp*test_perf+0.368_dp*train_perf)
  end function
end module mlr_metrics
