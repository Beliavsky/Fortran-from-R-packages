program test_randomized_pareto
  use rmoo
  implicit none
  integer,parameter :: maxn=25,maxm=5
  real(dp) :: f(maxn,maxm)
  integer :: r1(maxn),r2(maxn)
  integer :: case_id,n,m

  call ga_seed(987654)
  do case_id=1,250
    n=5+mod(case_id,21)
    m=2+mod(case_id,4)
    call random_number(f(1:n,1:m))
    if(mod(case_id,7)==0)f(2,1:m)=f(1,1:m)
    call non_dominated_sort(f(1:n,1:m),r1(1:n))
    call brute_layers(f(1:n,1:m),r2(1:n))
    if(any(r1(1:n)/=r2(1:n)))then
      write(*,*)"FAIL randomized Pareto case",case_id
      write(*,*)"translated",r1(1:n)
      write(*,*)"reference ",r2(1:n)
      error stop 1
    end if
  end do
  print *,"test_randomized_pareto: PASS (250 cases)"
contains
  subroutine brute_layers(x,rank)
    real(dp),intent(in)::x(:,:)
    integer,intent(out)::rank(size(x,1))
    logical :: active(size(x,1)),isdom
    integer::layer,i,j,nleft
    rank=0;active=.true.;layer=1;nleft=size(x,1)
    do while(nleft>0)
      do i=1,size(x,1)
        if(.not.active(i))cycle
        isdom=.false.
        do j=1,size(x,1)
          if(i==j.or..not.active(j))cycle
          if(all(x(j,:)<=x(i,:)).and.any(x(j,:)<x(i,:)))then
            isdom=.true.;exit
          end if
        end do
        if(.not.isdom)rank(i)=-layer
      end do
      do i=1,size(x,1)
        if(rank(i)==-layer)then
          rank(i)=layer;active(i)=.false.;nleft=nleft-1
        end if
      end do
      layer=layer+1
    end do
  end subroutine brute_layers
end program test_randomized_pareto
