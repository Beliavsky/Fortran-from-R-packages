! SPDX-License-Identifier: MIT
! PortfolioTesteR modern Fortran translation
module ptr_ml
  use ptr_kinds, only : dp
  use ptr_types, only : linear_model
  use ptr_utils, only : nan_dp, is_finite, finite_mean, finite_sd, solve_linear
  use ptr_utils, only : rank_vector, normalize_nonnegative, percentile
  implicit none
  private
  public :: panel_lag, make_labels, transform_scores, select_top_k_scores
  public :: weight_from_scores, combine_scores, fit_linear_model, predict_linear_model
  public :: rolling_fit_predict, ic_series, evaluate_scores, bucket_returns
  public :: coverage_by_date, panel_op, add_interaction

contains

  subroutine panel_lag(panel,k,out)
    real(dp),intent(in)::panel(:,:)
    integer,intent(in)::k
    real(dp),allocatable,intent(out)::out(:,:)
    allocate(out(size(panel,1),size(panel,2)));out=nan_dp()
    if(k>=0.and.k<size(panel,1))out(k+1:,:)=panel(:size(panel,1)-k,:)
  end subroutine panel_lag

  subroutine make_labels(prices,horizon,label_type,out)
    real(dp),intent(in)::prices(:,:)
    integer,intent(in)::horizon,label_type
    real(dp),allocatable,intent(out)::out(:,:)
    real(dp)::r
    integer::t,j
    allocate(out(size(prices,1),size(prices,2)));out=nan_dp()
    do j=1,size(prices,2)
      do t=1,size(prices,1)-horizon
        if(prices(t,j)<=0.0_dp.or.prices(t+horizon,j)<=0.0_dp)cycle
        select case(label_type)
        case(1);r=log(prices(t+horizon,j)/prices(t,j))
        case(2);r=prices(t+horizon,j)/prices(t,j)-1.0_dp
        case(3);r=merge(1.0_dp,-1.0_dp,prices(t+horizon,j)>prices(t,j))
        case default;r=log(prices(t+horizon,j)/prices(t,j))
        end select
        out(t,j)=r
      end do
    end do
  end subroutine make_labels

  subroutine panel_op(a,b,out,operation,fill)
    real(dp),intent(in)::a(:,:),b(:,:)
    real(dp),allocatable,intent(out)::out(:,:)
    character(len=*),intent(in)::operation
    real(dp),intent(in),optional::fill
    real(dp)::f,x,y
    integer::i,j
    f=nan_dp();if(present(fill))f=fill
    allocate(out(min(size(a,1),size(b,1)),min(size(a,2),size(b,2))));out=nan_dp()
    do j=1,size(out,2);do i=1,size(out,1)
      x=a(i,j);y=b(i,j);if(.not.is_finite(x))x=f;if(.not.is_finite(y))y=f
      if(.not.is_finite(x).or..not.is_finite(y))cycle
      select case(trim(operation))
      case('+');out(i,j)=x+y
      case('-');out(i,j)=x-y
      case('/');if(abs(y)>tiny(1.0_dp))out(i,j)=x/y
      case default;out(i,j)=x*y
      end select
    end do;end do
  end subroutine panel_op

  subroutine add_interaction(a,b,out)
    real(dp),intent(in)::a(:,:),b(:,:)
    real(dp),allocatable,intent(out)::out(:,:)
    call panel_op(a,b,out,'*')
  end subroutine add_interaction

  subroutine transform_scores(scores,out,method,robust)
    real(dp),intent(in)::scores(:,:)
    real(dp),allocatable,intent(out)::out(:,:)
    character(len=*),intent(in)::method
    logical,intent(in),optional::robust
    logical::use_robust
    real(dp),allocatable::r(:),absdev(:)
    real(dp)::mu,sd,med,mad
    integer::t
    use_robust=.false.;if(present(robust))use_robust=robust
    allocate(out(size(scores,1),size(scores,2)));out=nan_dp()
    do t=1,size(scores,1)
      select case(trim(method))
      case('rank')
        call rank_vector(scores(t,:),r,descending=.false.,normalize=.true.);out(t,:)=r
      case('zscore')
        if(use_robust)then
          med=percentile(scores(t,:),0.5_dp);allocate(absdev(size(scores,2)))
          absdev=abs(scores(t,:)-med);mad=percentile(absdev,0.5_dp);deallocate(absdev)
          if(is_finite(mad).and.mad>0.0_dp)then
            where(is_finite(scores(t,:)))out(t,:)=(scores(t,:)-med)/(1.4826_dp*mad)
          end if
        else
          mu=finite_mean(scores(t,:));sd=finite_sd(scores(t,:))
          if(is_finite(sd).and.sd>0.0_dp)then
            where(is_finite(scores(t,:)))out(t,:)=(scores(t,:)-mu)/sd
          end if
        end if
      case default
        out(t,:)=scores(t,:)
      end select
    end do
  end subroutine transform_scores

  subroutine select_top_k_scores(scores,k,mask,groups,max_per_group)
    real(dp),intent(in)::scores(:,:)
    integer,intent(in)::k
    real(dp),allocatable,intent(out)::mask(:,:)
    integer,intent(in),optional::groups(:),max_per_group
    logical,allocatable::used(:)
    integer,allocatable::gcount(:)
    integer::t,j,m,best,ng,limit
    real(dp)::bestv
    allocate(mask(size(scores,1),size(scores,2)),used(size(scores,2)));mask=0.0_dp
    ng=0;if(present(groups))ng=maxval(groups)
    allocate(gcount(max(1,ng)));limit=huge(1);if(present(max_per_group))limit=max_per_group
    do t=1,size(scores,1)
      used=.false.;gcount=0
      do m=1,min(k,count(is_finite(scores(t,:))))
        best=0;bestv=-huge(1.0_dp)
        do j=1,size(scores,2)
          if(used(j).or..not.is_finite(scores(t,j)))cycle
          if(present(groups))then
            if(groups(j)<1.or.groups(j)>ng)cycle
            if(gcount(groups(j))>=limit)cycle
          end if
          if(best==0.or.scores(t,j)>bestv)then;best=j;bestv=scores(t,j);end if
        end do
        if(best==0)exit
        used(best)=.true.;mask(t,best)=1.0_dp
        if(present(groups))gcount(groups(best))=gcount(groups(best))+1
      end do
    end do
  end subroutine select_top_k_scores

  subroutine weight_from_scores(scores,weights,method,temperature,floor_weight)
    real(dp),intent(in)::scores(:,:)
    real(dp),allocatable,intent(out)::weights(:,:)
    character(len=*),intent(in)::method
    real(dp),intent(in),optional::temperature,floor_weight
    real(dp),allocatable::row(:),r(:)
    real(dp)::temp,floorw,mx
    integer::t,j,n
    temp=1.0_dp;if(present(temperature))temp=max(1.0e-8_dp,temperature)
    floorw=0.0_dp;if(present(floor_weight))floorw=max(0.0_dp,floor_weight)
    allocate(weights(size(scores,1),size(scores,2)),row(size(scores,2)));weights=0.0_dp
    do t=1,size(scores,1)
      row=0.0_dp;n=count(is_finite(scores(t,:)));if(n==0)cycle
      select case(trim(method))
      case('equal')
        where(is_finite(scores(t,:)))row=1.0_dp
      case('rank')
        call rank_vector(scores(t,:),r,descending=.true.,normalize=.false.)
        do j=1,size(row);if(is_finite(r(j)))row(j)=real(n,dp)-r(j)+1.0_dp;end do
      case('linear')
        mx=minval(scores(t,:),mask=is_finite(scores(t,:)))
        where(is_finite(scores(t,:)))row=max(scores(t,:)-mx+floorw,0.0_dp)
      case default
        mx=maxval(scores(t,:),mask=is_finite(scores(t,:)))
        where(is_finite(scores(t,:)))row=exp((scores(t,:)-mx)*temp)+floorw
      end select
      call normalize_nonnegative(row);weights(t,:)=row
    end do
  end subroutine weight_from_scores

  subroutine combine_scores(score_cube,out,method,blend,trim_fraction)
    real(dp),intent(in)::score_cube(:,:,:)
    real(dp),allocatable,intent(out)::out(:,:)
    character(len=*),intent(in)::method
    real(dp),intent(in),optional::blend(:),trim_fraction
    real(dp),allocatable::b(:),v(:)
    real(dp)::s,trim_frac
    integer::t,j,k,n,lo,hi
    allocate(out(size(score_cube,1),size(score_cube,2)));out=nan_dp()
    allocate(b(size(score_cube,3)),v(size(score_cube,3)));b=1.0_dp
    if(present(blend))then;if(size(blend)==size(b))b=blend;end if
    if(sum(abs(b))>0.0_dp)b=b/sum(b)
    trim_frac=0.1_dp
    if(present(trim_fraction))trim_frac=max(0.0_dp,min(0.49_dp,trim_fraction))
    do t=1,size(out,1);do j=1,size(out,2)
      n=0
      do k=1,size(score_cube,3);if(is_finite(score_cube(t,j,k)))then;n=n+1;v(n)=score_cube(t,j,k);end if;end do
      if(n==0)cycle
      select case(trim(method))
      case('weighted')
        s=0.0_dp;out(t,j)=0.0_dp
        do k=1,size(score_cube,3)
          if(is_finite(score_cube(t,j,k)))then;out(t,j)=out(t,j)+b(k)*score_cube(t,j,k);s=s+b(k);end if
        end do
        if(s>0.0_dp)out(t,j)=out(t,j)/s
      case('trimmed_mean')
        call sort_local(v(:n));lo=1+int(trim_frac*real(n,dp));hi=n-int(trim_frac*real(n,dp))
        if(hi>=lo)out(t,j)=sum(v(lo:hi))/real(hi-lo+1,dp)
      case default
        out(t,j)=sum(v(:n))/real(n,dp)
      end select
    end do;end do
  contains
    subroutine sort_local(x)
      real(dp),intent(inout)::x(:)
      real(dp)::key
      integer::ii,jj
      do ii = 2, size(x)
        key = x(ii)
        jj = ii - 1
        do while (jj >= 1)
          if (x(jj) <= key) exit
          x(jj+1) = x(jj)
          jj = jj - 1
        end do
        x(jj+1) = key
      end do
    end subroutine sort_local
  end subroutine combine_scores

  subroutine fit_linear_model(x,y,model,lambda,standardize,status)
    real(dp),intent(in)::x(:,:),y(:)
    type(linear_model),intent(out)::model
    real(dp),intent(in),optional::lambda
    logical,intent(in),optional::standardize
    integer,intent(out),optional::status
    logical::std
    real(dp),allocatable::z(:,:),a(:,:),rhs(:),coef(:)
    real(dp)::lam
    integer::n,p,j,st
    n=size(x,1);p=size(x,2);lam=0.0_dp;if(present(lambda))lam=max(0.0_dp,lambda)
    std=.true.;if(present(standardize))std=standardize
    allocate(model%center(p),model%scale(p),z(n,p+1));model%center=0.0_dp;model%scale=1.0_dp
    z(:,1)=1.0_dp
    do j=1,p
      if(std)then;model%center(j)=finite_mean(x(:,j));model%scale(j)=finite_sd(x(:,j));end if
      if(.not.is_finite(model%scale(j)).or.model%scale(j)<=0.0_dp)model%scale(j)=1.0_dp
      z(:,j+1)=(x(:,j)-model%center(j))/model%scale(j)
    end do
    allocate(a(p+1,p+1),rhs(p+1));a=matmul(transpose(z),z);rhs=matmul(transpose(z),y)
    do j=2,p+1;a(j,j)=a(j,j)+lam;end do
    call solve_linear(a,rhs,coef,st)
    allocate(model%coef(p+1));model%coef=0.0_dp
    if(st==0)then;model%coef=coef;model%fitted=.true.;else;model%fitted=.false.;end if
    model%lambda=lam;if(present(status))status=st
  end subroutine fit_linear_model

  subroutine predict_linear_model(model,x,prediction)
    type(linear_model),intent(in)::model
    real(dp),intent(in)::x(:,:)
    real(dp),allocatable,intent(out)::prediction(:)
    integer::i,j
    allocate(prediction(size(x,1)));prediction=nan_dp()
    if(.not.model%fitted)return
    do i=1,size(x,1)
      if(count(is_finite(x(i,:)))<size(x,2))cycle
      prediction(i)=model%coef(1)
      do j=1,size(x,2);prediction(i)=prediction(i)+model%coef(j+1)*(x(i,j)-model%center(j))/model%scale(j);end do
    end do
  end subroutine predict_linear_model

  subroutine rolling_fit_predict(features,labels,is_periods,oos_periods,step,scores,lambda)
    real(dp),intent(in)::features(:,:,:),labels(:,:)
    integer,intent(in)::is_periods,oos_periods,step
    real(dp),allocatable,intent(out)::scores(:,:)
    real(dp),intent(in),optional::lambda
    real(dp),allocatable::xtrain(:,:),ytrain(:),xpred(:,:),pred(:)
    type(linear_model)::model
    integer::start,isend,oosend,t,j,ntrain,nrow,k,st
    allocate(scores(size(labels,1),size(labels,2)));scores=nan_dp();start=1
    do while(start+is_periods+oos_periods-1<=size(labels,1))
      isend=start+is_periods-1;oosend=isend+oos_periods;ntrain=0
      do t=start,isend;do j=1,size(labels,2)
        if(is_finite(labels(t,j)).and.count(is_finite(features(t,j,:)))==size(features,3))ntrain=ntrain+1
      end do;end do
      if(ntrain>size(features,3))then
        allocate(xtrain(ntrain,size(features,3)),ytrain(ntrain));k=0
        do t=start,isend;do j=1,size(labels,2)
          if(is_finite(labels(t,j)).and.count(is_finite(features(t,j,:)))==size(features,3))then
            k=k+1;xtrain(k,:)=features(t,j,:);ytrain(k)=labels(t,j)
          end if
        end do;end do
        call fit_linear_model(xtrain,ytrain,model,lambda,status=st)
        if(st==0)then
          nrow=(oosend-isend)*size(labels,2);allocate(xpred(nrow,size(features,3)));xpred=nan_dp();k=0
          do t=isend+1,oosend;do j=1,size(labels,2);k=k+1;xpred(k,:)=features(t,j,:);end do;end do
          call predict_linear_model(model,xpred,pred);k=0
          do t=isend+1,oosend;do j=1,size(labels,2);k=k+1;scores(t,j)=pred(k);end do;end do
          deallocate(xpred,pred)
        end if
        deallocate(xtrain,ytrain)
      end if
      start=start+step
    end do
  end subroutine rolling_fit_predict

  subroutine ic_series(scores,labels,ic,spearman)
    real(dp),intent(in)::scores(:,:),labels(:,:)
    real(dp),allocatable,intent(out)::ic(:)
    logical,intent(in),optional::spearman
    logical::spr
    real(dp),allocatable::x(:),y(:),rx(:),ry(:)
    integer::t,j,n,k
    spr=.true.;if(present(spearman))spr=spearman
    allocate(ic(size(scores,1)));ic=nan_dp()
    do t=1,size(scores,1)
      n=0;do j=1,size(scores,2);if(is_finite(scores(t,j)).and.is_finite(labels(t,j)))n=n+1;end do
      if(n<2)cycle;allocate(x(n),y(n));k=0
      do j = 1, size(scores,2)
        if (is_finite(scores(t,j)) .and. is_finite(labels(t,j))) then
          k = k + 1
          x(k) = scores(t,j)
          y(k) = labels(t,j)
        end if
      end do
      if(spr)then;call rank_vector(x,rx);call rank_vector(y,ry);ic(t)=corr_local(rx,ry);else;ic(t)=corr_local(x,y);end if
      deallocate(x,y);if(allocated(rx))deallocate(rx);if(allocated(ry))deallocate(ry)
    end do
  contains
    real(dp) function corr_local(a,b)
      real(dp),intent(in)::a(:),b(:)
      real(dp)::ma,mb,sa,sb,sab
      ma=sum(a)/real(size(a),dp);mb=sum(b)/real(size(b),dp)
      sa=sum((a-ma)**2);sb=sum((b-mb)**2);sab=sum((a-ma)*(b-mb))
      if(sa>0.0_dp.and.sb>0.0_dp)then;corr_local=sab/sqrt(sa*sb);else;corr_local=nan_dp();end if
    end function corr_local
  end subroutine ic_series

  subroutine evaluate_scores(scores,labels,mean_ic,top_return,bottom_return,top_fraction)
    real(dp),intent(in)::scores(:,:),labels(:,:)
    real(dp),intent(out)::mean_ic,top_return,bottom_return
    real(dp),intent(in),optional::top_fraction
    real(dp),allocatable::ic(:),r(:)
    real(dp)::frac
    integer::t,j,n,k,ntop
    frac=0.2_dp;if(present(top_fraction))frac=top_fraction
    call ic_series(scores,labels,ic,.true.);mean_ic=finite_mean(ic)
    top_return=0.0_dp;bottom_return=0.0_dp;n=0
    do t=1,size(scores,1)
      if(count(is_finite(scores(t,:)).and.is_finite(labels(t,:)))<2)cycle
      call rank_vector(scores(t,:), r, descending=.false., normalize=.false.)
      ntop = max(1, int(frac*real(count(is_finite(r)),dp)))
      do j=1,size(scores,2)
        if(.not.is_finite(r(j)).or..not.is_finite(labels(t,j)))cycle
        k=count(is_finite(r))
        if(r(j)>real(k-ntop,dp))top_return=top_return+labels(t,j)/real(ntop,dp)
        if(r(j)<=real(ntop,dp))bottom_return=bottom_return+labels(t,j)/real(ntop,dp)
      end do
      n=n+1
    end do
    if (n > 0) then
      top_return = top_return / real(n,dp)
      bottom_return = bottom_return / real(n,dp)
    else
      top_return = nan_dp()
      bottom_return = nan_dp()
    end if
  end subroutine evaluate_scores

  subroutine bucket_returns(scores,labels,n_buckets,means,counts)
    real(dp),intent(in)::scores(:,:),labels(:,:)
    integer,intent(in)::n_buckets
    real(dp),allocatable,intent(out)::means(:)
    integer,allocatable,intent(out)::counts(:)
    real(dp),allocatable::r(:)
    integer::t,j,n,b
    allocate(means(n_buckets),counts(n_buckets));means=0.0_dp;counts=0
    do t=1,size(scores,1)
      call rank_vector(scores(t,:),r,descending=.false.,normalize=.true.)
      do j=1,size(scores,2)
        if(.not.is_finite(r(j)).or..not.is_finite(labels(t,j)))cycle
        b=min(n_buckets,1+int(r(j)*real(n_buckets,dp)))
        means(b)=means(b)+labels(t,j);counts(b)=counts(b)+1
      end do
    end do
    do n=1,n_buckets;if(counts(n)>0)means(n)=means(n)/real(counts(n),dp);end do
  end subroutine bucket_returns

  subroutine coverage_by_date(panel,coverage)
    real(dp),intent(in)::panel(:,:)
    real(dp),allocatable,intent(out)::coverage(:)
    integer::t
    allocate(coverage(size(panel,1)))
    do t=1,size(panel,1);coverage(t)=real(count(is_finite(panel(t,:))),dp)/real(size(panel,2),dp);end do
  end subroutine coverage_by_date

end module ptr_ml
