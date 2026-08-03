module ffp_views
  use ffp_kinds, only : dp
  implicit none
  private
  public :: view_on_mean, view_on_covariance, view_on_volatility, view_on_correlation, view_on_rank
contains
  subroutine view_on_mean(x,mean,aeq,beq)
    real(dp),intent(in)::x(:,:),mean(:); real(dp),intent(out)::aeq(:,:),beq(:)
    aeq=transpose(x); beq=mean
  end subroutine
  subroutine view_on_covariance(x,mean,sigma,aeq,beq)
    real(dp),intent(in)::x(:,:),mean(:),sigma(:,:); real(dp),intent(out)::aeq(:,:),beq(:)
    integer::i,j,row,k; real(dp),allocatable::sec(:,:)
    k=size(x,2); allocate(sec(k,k)); sec=sigma+spread(mean,2,k)*spread(mean,1,k); row=0
    do i=1,k; do j=i,k; row=row+1; aeq(row,:)=x(:,i)*x(:,j); beq(row)=sec(i,j); end do; end do
  end subroutine
  subroutine view_on_volatility(x,vol,aeq,beq)
    real(dp),intent(in)::x(:,:),vol(:); real(dp),intent(out)::aeq(:,:),beq(:)
    integer::j
    do j=1,size(x,2); aeq(j,:)=x(:,j)**2; beq(j)=vol(j)**2; end do
  end subroutine
  subroutine view_on_correlation(x,cor,aeq,beq)
    real(dp),intent(in)::x(:,:),cor(:,:); real(dp),intent(out)::aeq(:,:),beq(:)
    real(dp),allocatable::sd(:); integer::i,j,row,k
    k=size(x,2); allocate(sd(k)); do i=1,k; sd(i)=sqrt(sum(x(:,i)**2)/real(size(x,1),dp)); end do
    row=0; do i=1,k; do j=i+1,k; row=row+1; aeq(row,:)=x(:,i)*x(:,j); beq(row)=cor(i,j)*sd(i)*sd(j); end do; end do
  end subroutine
  subroutine view_on_rank(x,rank_target,aeq,beq)
    real(dp),intent(in)::x(:,:),rank_target(:); real(dp),intent(out)::aeq(:,:),beq(:)
    integer::i,j,n,k; real(dp)::r
    n=size(x,1); k=size(x,2)
    do j=1,k
      do i=1,n
        r=1.0_dp+real(count(x(:,j)<x(i,j)),dp)+0.5_dp*real(count(abs(x(:,j)-x(i,j)) <= 10.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(x(i,j))))-1,dp)
        aeq(j,i)=r
      end do
      beq(j)=rank_target(j)
    end do
  end subroutine
end module ffp_views
