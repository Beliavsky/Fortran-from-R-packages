! Modern Fortran translation of the computational core of DiceKriging 1.6.1.
! Upstream DiceKriging is distributed under GPL-2 | GPL-3.
! This translation is distributed under the same license choice; see
! LICENSE-GPL-2 and LICENSE-GPL-3 in the project root.
module dk_diagnostics
  use dk_kinds, only : dp
  use dk_model, only : km_model, km_prediction, km_predict, km_recompute
  implicit none
  private
  public :: cv_predict
contains
  subroutine cv_predict(model,fold,kind,mean,var,trend_reestimate)
    type(km_model),intent(in)::model
    integer,intent(in)::fold(:)
    character(len=*),intent(in)::kind
    real(dp),allocatable,intent(out)::mean(:),var(:)
    logical,intent(in),optional::trend_reestimate
    logical::tr
    integer::g,ng,nkeep,ntest,i,j,ii,it
    integer,allocatable::keep(:),test(:)
    type(km_model)::m2 = km_model()
    type(km_prediction)::pr
    real(dp),allocatable::xt(:,:),ft(:,:)
    tr=.true.;if(present(trend_reestimate))tr=trend_reestimate
    if(size(fold)/=model%n)error stop 'cv_predict: fold length mismatch'
    if(model%noise_flag)error stop 'cv_predict: noisy observations are not supported'
    ng=maxval(fold);allocate(mean(model%n),var(model%n));mean=0.0_dp;var=0.0_dp
    do g=1,ng
      nkeep=count(fold/=g);ntest=count(fold==g);allocate(keep(nkeep),test(ntest));ii=0;it=0
      do i=1,model%n
        if(fold(i)==g)then;it=it+1;test(it)=i;else;ii=ii+1;keep(ii)=i;end if
      end do
      m2=model;m2%n=nkeep;m2%x=model%x(keep,:);m2%y=model%y(keep);m2%f=model%f(keep,:)
      call km_recompute(m2,reestimate_trend=tr)
      xt=model%x(test,:);ft=model%f(test,:);call km_predict(m2,xt,ft,kind,pr,se_compute=.true.)
      do j=1,ntest;mean(test(j))=pr%mean(j);var(test(j))=pr%sd(j)**2;end do
      deallocate(keep,test,xt,ft)
    end do
  end subroutine cv_predict
end module dk_diagnostics
