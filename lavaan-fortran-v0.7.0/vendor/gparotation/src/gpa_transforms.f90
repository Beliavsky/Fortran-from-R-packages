module gpa_transforms
  use gpa_kinds, only: dp
  use gpa_linalg, only: eye, inverse_matrix, cholesky_lower
  use gpa_rotation, only: rotation_result
  implicit none
  private
  public :: eiv_rotate, echelon_rotate
contains
  subroutine eiv_rotate(a, identity_rows, r)
    real(dp),intent(in)::a(:,:)
    integer,intent(in)::identity_rows(:)
    type(rotation_result),intent(out)::r
    integer::p,k,i,info
    real(dp),allocatable::a1(:,:),ainv(:,:)
    p=size(a,1)
    k=size(a,2)
    call init_r(r,p,k,.false.)
    if(size(identity_rows)/=k) then
    r%info=-1
    return
    end if
    allocate(a1(k,k),ainv(k,k))
    do i=1,k
      if(identity_rows(i)<1 .or. identity_rows(i)>p) then
      r%info=-2
      return
      end if
      a1(i,:)=a(identity_rows(i),:)
    end do
    call inverse_matrix(a1,ainv,info)
    if(info/=0) then
    r%info=info
    return
    end if
    r%loadings=matmul(a,ainv)
    do i=1,k
    r%loadings(identity_rows(i),:)=0.0_dp
    r%loadings(identity_rows(i),i)=1.0_dp
    end do
    r%phi=matmul(a1,transpose(a1))
    r%th=transpose(a1)
    r%method='eiv'
    r%converged=.true.
    r%info=0
  end subroutine eiv_rotate

  subroutine echelon_rotate(a, reference_rows, r)
    real(dp),intent(in)::a(:,:)
    integer,intent(in)::reference_rows(:)
    type(rotation_result),intent(out)::r
    integer::p,k,i,info
    real(dp),allocatable::a1(:,:),b1(:,:),ainv(:,:)
    p=size(a,1)
    k=size(a,2)
    call init_r(r,p,k,.true.)
    if(size(reference_rows)/=k) then
    r%info=-1
    return
    end if
    allocate(a1(k,k),b1(k,k),ainv(k,k))
    do i=1,k
      if(reference_rows(i)<1 .or. reference_rows(i)>p) then
      r%info=-2
      return
      end if
      a1(i,:)=a(reference_rows(i),:)
    end do
    call cholesky_lower(matmul(a1,transpose(a1)),b1,info)
    if(info/=0) then
    r%info=info
    return
    end if
    call inverse_matrix(a1,ainv,info)
    if(info/=0) then
    r%info=info
    return
    end if
    r%th=matmul(ainv,b1)
    r%loadings=matmul(a,r%th)
    do i=1,k
    r%loadings(reference_rows(i),:)=b1(i,:)
    end do
    r%phi=eye(k)
    r%method='echelon'
    r%converged=.true.
    r%info=0
  end subroutine echelon_rotate

  subroutine init_r(r,p,k,orth)
    type(rotation_result),intent(out)::r
    integer,intent(in)::p,k
    logical,intent(in)::orth
    allocate(r%loadings(p,k),r%th(k,k),r%phi(k,k),r%gq(p,k),r%table(1,4))
    r%loadings=0.0_dp
    r%th=eye(k)
    r%phi=eye(k)
    r%gq=0.0_dp
    r%table=0.0_dp
    r%orthogonal=orth
    r%converged=.false.
    r%info=0
    r%objective=0.0_dp
  end subroutine init_r
end module gpa_transforms
