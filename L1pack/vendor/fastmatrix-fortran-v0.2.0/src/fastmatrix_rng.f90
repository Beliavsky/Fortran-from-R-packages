module fastmatrix_rng
  use fastmatrix_base, only: dp, normal_rand, chol_lower
  implicit none
  private
  public :: rmnorm, rsphere, rball
contains
  subroutine rmnorm(n,mean,sigma,y,info)
    integer,intent(in)::n
    real(dp),intent(in)::mean(:),sigma(:,:)
    real(dp),intent(out)::y(n,size(mean))
    integer,intent(out),optional::info
    real(dp)::l(size(mean),size(mean)),z(size(mean))
    integer::i,j,ier
    call chol_lower(sigma,l,ier)
    if(present(info))info=ier
    if(ier/=0)return
    do i=1,n
    do j=1,size(mean)
    z(j)=normal_rand()
    end do
    y(i,:)=mean+matmul(l,z)
    end do
  end subroutine
  subroutine rsphere(n,p,y)
    integer,intent(in)::n,p
    real(dp),intent(out)::y(n,p)
    real(dp)::s
    integer::i,j
    do i=1,n
    do j=1,p
    y(i,j)=normal_rand()
    end do
    s=sqrt(sum(y(i,:)**2))
    y(i,:)=y(i,:)/s
    end do
  end subroutine
  subroutine rball(n,p,y)
    integer,intent(in)::n,p
    real(dp),intent(out)::y(n,p)
    real(dp)::u
    integer::i
    call rsphere(n,p,y)
    do i=1,n
    call random_number(u)
    y(i,:)=y(i,:)*u**(1.0_dp/real(p,dp))
    end do
  end subroutine
end module
