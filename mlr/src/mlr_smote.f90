module mlr_smote
  use mlr_kinds, only : dp
  use mlr_rng, only : rng_state, rng_integer, rng_uniform
  implicit none
  private
  public :: smote_generate, nearest_neighbors
contains
  subroutine nearest_neighbors(x,k,nn)
    real(dp),intent(in)::x(:,:);integer,intent(in)::k
    integer,allocatable,intent(out)::nn(:,:)
    real(dp),allocatable::d(:);integer::i,j,t,best
    if(k<1.or.k>=size(x,1))error stop "nearest_neighbors: invalid k"
    allocate(nn(size(x,1),k),d(size(x,1)))
    do i=1,size(x,1)
      do j=1,size(x,1);d(j)=sum((x(i,:)-x(j,:))**2);end do;d(i)=huge(1.0_dp)
      do t=1,k
        best=minloc(d,dim=1);nn(i,t)=best;d(best)=huge(1.0_dp)
      end do
    end do
  end subroutine

  subroutine smote_generate(x,is_numeric,nn,nnew,rng,res)
    real(dp),intent(in)::x(:,:);logical,intent(in)::is_numeric(:);integer,intent(in)::nn(:,:),nnew
    type(rng_state),intent(inout)::rng;real(dp),allocatable,intent(out)::res(:,:)
    integer::i,j,j_sel,j_nn;real(dp)::lambda
    if(size(is_numeric)/=size(x,2).or.size(nn,1)/=size(x,1))error stop "smote_generate: dimensions"
    allocate(res(nnew,size(x,2)))
    do i=1,nnew
      j_sel=rng_integer(rng,1,size(x,1));j_nn=nn(j_sel,rng_integer(rng,1,size(nn,2)))
      lambda=rng_uniform(rng)
      do j=1,size(x,2)
        if(is_numeric(j))then
          res(i,j)=lambda*x(j_sel,j)+(1.0_dp-lambda)*x(j_nn,j)
        else
          if(lambda<0.5_dp)then;res(i,j)=x(j_sel,j);else;res(i,j)=x(j_nn,j);end if
        end if
      end do
    end do
  end subroutine
end module mlr_smote
