module mixsqp_em
  use mixsqp_kinds, only : dp
  use mixsqp_utils, only : mixobjective
  implicit none
  private
  public :: mixem_update, run_mixem
contains
  subroutine mixem_update(L,w,x)
    real(dp), intent(in) :: L(:,:), w(:)
    real(dp), intent(inout) :: x(:)
    real(dp), allocatable :: P(:,:)
    real(dp) :: rmax,rsum
    integer :: i,j
    real(dp), parameter :: tiny_p = 1.0e-15_dp
    allocate(P(size(L,1),size(L,2)))
    P = L
    do j=1,size(P,2)
      P(:,j)=P(:,j)*(x(j)+tiny_p)
    end do
    do i=1,size(P,1)
      rmax=maxval(P(i,:))
      if (rmax>0.0_dp) P(i,:)=P(i,:)/rmax
      P(i,:)=P(i,:)+tiny_p
      rsum=sum(P(i,:))
      P(i,:)=P(i,:)/rsum
    end do
    x=matmul(transpose(P),w)
  end subroutine mixem_update

  subroutine run_mixem(L,w,z,x,e,numiter,zero_threshold,obj,nnz,dmax)
    real(dp), intent(in) :: L(:,:), w(:), z(:), e(:), zero_threshold
    real(dp), intent(inout) :: x(:)
    integer, intent(in) :: numiter
    real(dp), intent(out) :: obj(numiter), dmax(numiter)
    integer, intent(out) :: nnz(numiter)
    real(dp), allocatable :: xold(:)
    integer :: i
    allocate(xold(size(x)))
    do i=1,numiter
      xold=x
      call mixem_update(L,w,x)
      obj(i)=mixobjective(L,x,w,z,e)
      nnz(i)=count(x>zero_threshold)
      dmax(i)=maxval(abs(x-xold))
    end do
  end subroutine run_mixem
end module mixsqp_em
