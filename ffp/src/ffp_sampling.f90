module ffp_sampling
  use ffp_kinds, only : dp
  implicit none
  private
  public :: bootstrap_scenarios
contains
  subroutine bootstrap_scenarios(x,p,sample)
    real(dp),intent(in)::x(:,:),p(:); real(dp),intent(out)::sample(:,:)
    real(dp),allocatable::cdf(:),u(:); integer::i,j,n
    n=size(p); allocate(cdf(n),u(size(sample,1))); cdf(1)=p(1)
    do i=2,n; cdf(i)=cdf(i-1)+p(i); end do
    call random_number(u)
    do i=1,size(sample,1)
      j=1; do while(j<n .and. u(i)>cdf(j)); j=j+1; end do
      sample(i,:)=x(j,:)
    end do
  end subroutine
end module ffp_sampling
