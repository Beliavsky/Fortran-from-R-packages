! SPDX-License-Identifier: MIT
! PortfolioTesteR modern Fortran translation
module ptr_weighting
  use ptr_kinds, only : dp
  use ptr_utils, only : is_finite, nan_dp, normalize_nonnegative, rank_vector
  use ptr_utils, only : covariance_matrix, correlation_matrix, solve_linear
  implicit none
  private
  public :: weight_equally, weight_by_signal, weight_by_rank, weight_by_volatility
  public :: combine_weights, switch_weights, cap_exposure, cap_turnover
  public :: calculate_hrp_weights, calculate_erc_weights, calculate_max_div_weights
  public :: rolling_hrp_weights, rolling_risk_parity_weights

contains

  subroutine weight_equally(selection, weights)
    real(dp), intent(in) :: selection(:,:)
    real(dp), allocatable, intent(out) :: weights(:,:)
    integer :: t,n
    allocate(weights(size(selection,1),size(selection,2))); weights=0.0_dp
    do t=1,size(selection,1)
      n=count(is_finite(selection(t,:)) .and. selection(t,:)>0.0_dp)
      if(n>0) where(is_finite(selection(t,:)) .and. selection(t,:)>0.0_dp) weights(t,:)=1.0_dp/real(n,dp)
    end do
  end subroutine weight_equally

  subroutine weight_by_signal(selection, signal, weights, positive_only)
    real(dp), intent(in) :: selection(:,:), signal(:,:)
    real(dp), allocatable, intent(out) :: weights(:,:)
    logical, intent(in), optional :: positive_only
    logical :: pos
    real(dp), allocatable :: row(:)
    integer :: t,j
    pos=.true.; if(present(positive_only)) pos=positive_only
    allocate(weights(size(selection,1),size(selection,2)));weights=0.0_dp
    allocate(row(size(selection,2)))
    do t=1,size(selection,1)
      row=0.0_dp
      do j=1,size(selection,2)
        if(selection(t,j)>0.0_dp .and. is_finite(signal(t,j))) then
          if(pos) then
            row(j)=max(signal(t,j),0.0_dp)
          else
            row(j)=abs(signal(t,j))
          end if
        end if
      end do
      call normalize_nonnegative(row);weights(t,:)=row
    end do
  end subroutine weight_by_signal

  subroutine weight_by_rank(selection, signal, weights, exponential, ascending)
    real(dp), intent(in) :: selection(:,:), signal(:,:)
    real(dp), allocatable, intent(out) :: weights(:,:)
    logical, intent(in), optional :: exponential,ascending
    logical :: expo,asc
    real(dp),allocatable::v(:),r(:),row(:)
    integer::t,j,n
    expo=.false.;if(present(exponential))expo=exponential
    asc=.false.;if(present(ascending))asc=ascending
    allocate(weights(size(selection,1),size(selection,2)));weights=0.0_dp
    allocate(v(size(selection,2)),row(size(selection,2)))
    do t=1,size(selection,1)
      v=nan_dp();row=0.0_dp
      do j=1,size(selection,2)
        if(selection(t,j)>0.0_dp .and. is_finite(signal(t,j))) v(j)=signal(t,j)
      end do
      call rank_vector(v,r,descending=.not.asc)
      n=count(is_finite(r))
      do j=1,size(r)
        if(.not.is_finite(r(j)))cycle
        if(expo)then
          row(j)=0.5_dp**(r(j)-1.0_dp)
        else
          row(j)=real(n,dp)-r(j)+1.0_dp
        end if
      end do
      call normalize_nonnegative(row);weights(t,:)=row
    end do
  end subroutine weight_by_rank

  subroutine weight_by_volatility(selection, volatility, weights, floor_vol)
    real(dp),intent(in)::selection(:,:),volatility(:,:)
    real(dp),allocatable,intent(out)::weights(:,:)
    real(dp),intent(in),optional::floor_vol
    real(dp)::vf
    real(dp),allocatable::row(:)
    integer::t,j
    vf=1.0e-8_dp;if(present(floor_vol))vf=floor_vol
    allocate(weights(size(selection,1),size(selection,2)));weights=0.0_dp
    allocate(row(size(selection,2)))
    do t=1,size(selection,1)
      row=0.0_dp
      do j=1,size(selection,2)
        if (selection(t,j) > 0.0_dp .and. is_finite(volatility(t,j)) .and. &
            volatility(t,j) > 0.0_dp) then
          row(j) = 1.0_dp / max(vf, volatility(t,j))
        end if
      end do
      call normalize_nonnegative(row);weights(t,:)=row
    end do
  end subroutine weight_by_volatility

  subroutine combine_weights(weight_cube,weights,blend)
    real(dp),intent(in)::weight_cube(:,:,:)
    real(dp),allocatable,intent(out)::weights(:,:)
    real(dp),intent(in),optional::blend(:)
    real(dp),allocatable::b(:),row(:)
    integer::k,t
    allocate(weights(size(weight_cube,1),size(weight_cube,2)));weights=0.0_dp
    allocate(b(size(weight_cube,3)),row(size(weight_cube,2)))
    if(present(blend))then
      if(size(blend)==size(b))then;b=blend;else;b=1.0_dp;end if
    else;b=1.0_dp;end if
    call normalize_nonnegative(b)
    do k=1,size(weight_cube,3);weights=weights+b(k)*weight_cube(:,:,k);end do
    do t=1,size(weights,1);row=weights(t,:);call normalize_nonnegative(row);weights(t,:)=row;end do
  end subroutine combine_weights

  subroutine switch_weights(a,b,condition,out,partial_blend)
    real(dp),intent(in)::a(:,:),b(:,:),condition(:)
    real(dp),allocatable,intent(out)::out(:,:)
    real(dp),intent(in),optional::partial_blend
    real(dp)::q
    integer::t
    q=1.0_dp;if(present(partial_blend))q=max(0.0_dp,min(1.0_dp,partial_blend))
    allocate(out(size(a,1),size(a,2)));out=0.0_dp
    do t=1,size(a,1)
      if(condition(t)>0.0_dp)then;out(t,:)=(1.0_dp-q)*a(t,:)+q*b(t,:);else;out(t,:)=a(t,:);end if
      call normalize_nonnegative(out(t,:))
    end do
  end subroutine switch_weights

  subroutine cap_exposure(weights,max_per_symbol,out,groups,max_per_group,renormalize)
    real(dp),intent(in)::weights(:,:),max_per_symbol
    real(dp),allocatable,intent(out)::out(:,:)
    integer,intent(in),optional::groups(:)
    real(dp),intent(in),optional::max_per_group
    logical,intent(in),optional::renormalize
    logical::renorm
    integer::t,j,g,iter,ng
    real(dp)::s,excess,room,scale,group_cap,gsum
    renorm=.true.;if(present(renormalize))renorm=renormalize
    group_cap=1.0_dp;if(present(max_per_group))group_cap=max_per_group
    allocate(out(size(weights,1),size(weights,2)));out=weights
    where(.not.is_finite(out).or.out<0.0_dp)out=0.0_dp
    do t=1,size(out,1)
      call normalize_nonnegative(out(t,:))
      do iter=1,50
        excess=sum(max(out(t,:)-max_per_symbol,0.0_dp))
        where(out(t,:)>max_per_symbol)out(t,:)=max_per_symbol
        if(excess<=1.0e-13_dp)exit
        room=sum(max(max_per_symbol-out(t,:),0.0_dp))
        if(room<=tiny(1.0_dp))exit
        do j=1,size(out,2)
          if(out(t,j)<max_per_symbol)out(t,j)=out(t,j)+excess*(max_per_symbol-out(t,j))/room
        end do
      end do
      if (present(groups)) then
        if (size(groups) == size(out,2)) then
          ng=maxval(groups)
          do g=1,ng
            gsum=0.0_dp
            do j=1,size(out,2)
              if(groups(j)==g)gsum=gsum+out(t,j)
            end do
            if(gsum>group_cap .and. gsum>0.0_dp)then
              scale=group_cap/gsum
              do j=1,size(out,2)
                if(groups(j)==g)out(t,j)=out(t,j)*scale
              end do
            end if
          end do
        end if
      end if
      if(renorm)then;s=sum(out(t,:));if(s>0.0_dp)out(t,:)=out(t,:)/s;end if
    end do
  end subroutine cap_exposure

  subroutine cap_turnover(desired,max_turnover,executed)
    real(dp),intent(in)::desired(:,:),max_turnover
    real(dp),allocatable,intent(out)::executed(:,:)
    real(dp),allocatable::prev(:),row(:)
    real(dp)::turn,alpha
    integer::t
    allocate(executed(size(desired,1),size(desired,2)),prev(size(desired,2)),row(size(desired,2)))
    executed=0.0_dp;prev=0.0_dp
    do t=1,size(desired,1)
      row=desired(t,:);call normalize_nonnegative(row)
      turn=0.5_dp*sum(abs(row-prev))
      if(turn>max_turnover.and.turn>0.0_dp)then;alpha=max_turnover/turn;row=prev+alpha*(row-prev);end if
      call normalize_nonnegative(row);executed(t,:)=row;prev=row
    end do
  end subroutine cap_turnover

  subroutine calculate_hrp_weights(returns,weights,status)
    real(dp),intent(in)::returns(:,:)
    real(dp),allocatable,intent(out)::weights(:)
    integer,intent(out),optional::status
    real(dp),allocatable::cov(:,:),corr(:,:),dist(:,:),ivp(:),subcov(:,:)
    logical,allocatable::members(:,:),active(:)
    integer,allocatable::left(:),right(:),order(:)
    integer::n,node,a,b,besta,bestb,i,pos
    real(dp)::d,bestd,var1,var2,alpha
    n=size(returns,2);allocate(weights(n));weights=0.0_dp
    if(n==0)then;if(present(status))status=1;return;end if
    if(n==1)then;weights=1.0_dp;if(present(status))status=0;return;end if
    call covariance_matrix(returns,cov);call correlation_matrix(returns,corr)
    allocate(dist(n,n));dist=sqrt(max(0.0_dp,0.5_dp*(1.0_dp-corr)))
    allocate(members(n,2*n-1),active(2*n-1),left(2*n-1),right(2*n-1));members=.false.;active=.false.;left=0;right=0
    do i=1,n;members(i,i)=.true.;active(i)=.true.;end do
    do node=n+1,2*n-1
      bestd=huge(1.0_dp);besta=0;bestb=0
      do a=1,node-1
        if(.not.active(a))cycle
        do b=a+1,node-1
          if(.not.active(b))cycle
          d=cluster_distance(a,b)
          if(d<bestd)then;bestd=d;besta=a;bestb=b;end if
        end do
      end do
      left(node)=besta;right(node)=bestb;members(:,node)=members(:,besta).or.members(:,bestb)
      active(besta)=.false.;active(bestb)=.false.;active(node)=.true.
    end do
    allocate(order(n));pos=0;call traverse(2*n-1)
    weights=1.0_dp;call bisect(1,n)
    call normalize_nonnegative(weights)
    if(present(status))status=0
  contains
    real(dp) function cluster_distance(c1,c2)
      integer,intent(in)::c1,c2
      integer::ii,jj,cnt
      cluster_distance=0.0_dp;cnt=0
      do ii=1,n;if(.not.members(ii,c1))cycle
        do jj=1,n;if(.not.members(jj,c2))cycle
          cluster_distance=cluster_distance+dist(ii,jj);cnt=cnt+1
        end do
      end do
      if(cnt>0)cluster_distance=cluster_distance/real(cnt,dp)
    end function cluster_distance
    recursive subroutine traverse(nd)
      integer,intent(in)::nd
      if(nd<=n)then;pos=pos+1;order(pos)=nd;return;end if
      call traverse(left(nd));call traverse(right(nd))
    end subroutine traverse
    real(dp) function cluster_var(lo,hi)
      integer,intent(in)::lo,hi
      integer::m,ii,jj
      m=hi-lo+1;allocate(subcov(m,m),ivp(m))
      do ii=1,m
        do jj=1,m;subcov(ii,jj)=cov(order(lo+ii-1),order(lo+jj-1));end do
        ivp(ii)=1.0_dp/max(subcov(ii,ii),1.0e-12_dp)
      end do
      call normalize_nonnegative(ivp);cluster_var=dot_product(ivp,matmul(subcov,ivp))
      deallocate(subcov,ivp)
    end function cluster_var
    recursive subroutine bisect(lo,hi)
      integer,intent(in)::lo,hi
      integer::mid,ii
      if(lo>=hi)return
      mid=(lo+hi)/2;var1=cluster_var(lo,mid);var2=cluster_var(mid+1,hi)
      if(var1+var2>0.0_dp)then;alpha=var2/(var1+var2);else;alpha=0.5_dp;end if
      do ii=lo,mid;weights(order(ii))=weights(order(ii))*alpha;end do
      do ii=mid+1,hi;weights(order(ii))=weights(order(ii))*(1.0_dp-alpha);end do
      call bisect(lo,mid);call bisect(mid+1,hi)
    end subroutine bisect
  end subroutine calculate_hrp_weights

  subroutine calculate_erc_weights(cov,weights,status,max_iter,tol)
    real(dp),intent(in)::cov(:,:)
    real(dp),allocatable,intent(out)::weights(:)
    integer,intent(out),optional::status
    integer,intent(in),optional::max_iter
    real(dp),intent(in),optional::tol
    real(dp),allocatable::marg(:),rc(:)
    real(dp)::portvar,target,err,tolerance
    integer::n,iter,itmax
    n=size(cov,1);allocate(weights(n),marg(n),rc(n));weights=1.0_dp/real(max(1,n),dp)
    itmax=5000;if(present(max_iter))itmax=max_iter;tolerance=1.0e-9_dp;if(present(tol))tolerance=tol
    do iter=1,itmax
      marg=matmul(cov,weights);portvar=dot_product(weights,marg)
      if(portvar<=0.0_dp)exit
      rc=weights*marg;target=portvar/real(n,dp);err=maxval(abs(rc-target))/max(target,1.0e-16_dp)
      if(err<tolerance)exit
      where(rc>1.0e-18_dp)weights=weights*sqrt(target/rc)
      call normalize_nonnegative(weights)
    end do
    if(present(status))then;if(iter<=itmax)then;status=0;else;status=1;end if;end if
  end subroutine calculate_erc_weights

  subroutine calculate_max_div_weights(cov,weights,status)
    real(dp),intent(in)::cov(:,:)
    real(dp),allocatable,intent(out)::weights(:)
    integer,intent(out),optional::status
    real(dp),allocatable::vol(:),x(:),a(:,:)
    integer::i,n,st
    n=size(cov,1);allocate(vol(n),a(n,n));a=cov
    do i=1,n;vol(i)=sqrt(max(cov(i,i),0.0_dp));a(i,i)=a(i,i)+1.0e-10_dp;end do
    call solve_linear(a,vol,x,st)
    allocate(weights(n))
    if (st /= 0) then
      weights = 1.0_dp / real(n, dp)
    else
      weights = max(x, 0.0_dp)
      call normalize_nonnegative(weights)
    end if
    if(present(status))status=st
  end subroutine calculate_max_div_weights

  subroutine rolling_hrp_weights(prices,selection,lookback,weights)
    real(dp),intent(in)::prices(:,:),selection(:,:)
    integer,intent(in)::lookback
    real(dp),allocatable,intent(out)::weights(:,:)
    real(dp),allocatable::r(:,:),sub(:,:),w(:)
    integer,allocatable::idx(:)
    integer::t,j,k,nsel,st
    allocate(r(size(prices,1),size(prices,2)));r=nan_dp()
    do j = 1, size(prices,2)
      do t = 2, size(prices,1)
        if (prices(t,j) > 0.0_dp .and. prices(t-1,j) > 0.0_dp) then
          r(t,j) = log(prices(t,j) / prices(t-1,j))
        end if
      end do
    end do
    allocate(weights(size(prices,1),size(prices,2)));weights=0.0_dp
    do t=lookback+1,size(prices,1)
      nsel=count(selection(t,:)>0.0_dp);if(nsel==0)cycle;allocate(idx(nsel));k=0
      do j=1,size(prices,2);if(selection(t,j)>0.0_dp)then;k=k+1;idx(k)=j;end if;end do
      allocate(sub(lookback,nsel));do j=1,nsel;sub(:,j)=r(t-lookback+1:t,idx(j));end do
      call calculate_hrp_weights(sub,w,st);if(st==0)then;do j=1,nsel;weights(t,idx(j))=w(j);end do;end if
      deallocate(idx,sub,w)
    end do
  end subroutine rolling_hrp_weights

  subroutine rolling_risk_parity_weights(prices,selection,lookback,weights)
    real(dp),intent(in)::prices(:,:),selection(:,:)
    integer,intent(in)::lookback
    real(dp),allocatable,intent(out)::weights(:,:)
    real(dp),allocatable::r(:,:),sub(:,:),cov(:,:),w(:)
    integer,allocatable::idx(:)
    integer::t,j,k,nsel,st
    allocate(r(size(prices,1),size(prices,2)));r=nan_dp()
    do j = 1, size(prices,2)
      do t = 2, size(prices,1)
        if (prices(t,j) > 0.0_dp .and. prices(t-1,j) > 0.0_dp) then
          r(t,j) = log(prices(t,j) / prices(t-1,j))
        end if
      end do
    end do
    allocate(weights(size(prices,1),size(prices,2)));weights=0.0_dp
    do t=lookback+1,size(prices,1)
      nsel=count(selection(t,:)>0.0_dp);if(nsel==0)cycle;allocate(idx(nsel));k=0
      do j=1,size(prices,2);if(selection(t,j)>0.0_dp)then;k=k+1;idx(k)=j;end if;end do
      allocate(sub(lookback,nsel));do j=1,nsel;sub(:,j)=r(t-lookback+1:t,idx(j));end do
      call covariance_matrix(sub, cov)
      call calculate_erc_weights(cov, w, st)
      if (st == 0) then
        do j = 1, nsel
          weights(t,idx(j)) = w(j)
        end do
      end if
      deallocate(idx,sub,cov,w)
    end do
  end subroutine rolling_risk_parity_weights

end module ptr_weighting
