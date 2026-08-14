module mixsqp_utils
  use mixsqp_kinds, only : dp
  use mixsqp_lapack, only : dgesdd
  implicit none
  private
  public :: mixobjective, normalize_likelihoods, normalize_loglikelihoods
  public :: normalize_rows_with_logscale, compute_grad_hessian
  public :: truncated_svd, logspace, set_seed
contains
  function mixobjective(L,x,w,z,e) result(f)
    real(dp), intent(in) :: L(:,:), x(:), w(:)
    real(dp), intent(in), optional :: z(:), e(:)
    real(dp) :: f
    real(dp), allocatable :: u(:)
    integer :: n
    n = size(L,1)
    allocate(u(n))
    u = matmul(L,x)
    if (present(e)) u = u + e
    if (any(u <= 0.0_dp)) then
      f = huge(1.0_dp)
    else if (present(z)) then
      f = -sum(w*(z + log(u)))
    else
      f = -sum(w*log(u))
    end if
  end function mixobjective

  subroutine normalize_likelihoods(A)
    real(dp), intent(inout) :: A(:,:)
    integer :: i
    real(dp) :: z
    do i = 1, size(A,1)
      z = maxval(A(i,:))
      if (z > 0.0_dp) A(i,:) = A(i,:)/z
    end do
  end subroutine normalize_likelihoods

  subroutine normalize_loglikelihoods(logA,A)
    real(dp), intent(in) :: logA(:,:)
    real(dp), intent(out) :: A(size(logA,1),size(logA,2))
    integer :: i
    real(dp) :: z
    do i = 1, size(logA,1)
      z = maxval(logA(i,:))
      A(i,:) = exp(logA(i,:) - z)
    end do
  end subroutine normalize_loglikelihoods

  subroutine normalize_rows_with_logscale(A,z)
    real(dp), intent(inout) :: A(:,:)
    real(dp), intent(out) :: z(size(A,1))
    integer :: i
    real(dp) :: s
    do i = 1, size(A,1)
      s = maxval(A(i,:))
      if (s > 0.0_dp) then
        A(i,:) = A(i,:)/s
        z(i) = log(s)
      else
        z(i) = 0.0_dp
      end if
    end do
  end subroutine normalize_rows_with_logscale

  subroutine compute_grad_hessian(L,w,x,e,g,H,U,V,usesvd)
    real(dp), intent(in) :: L(:,:), w(:), x(:), e(:)
    real(dp), intent(out) :: g(size(L,2)), H(size(L,2),size(L,2))
    real(dp), intent(in), optional :: U(:,:), V(:,:)
    logical, intent(in), optional :: usesvd
    logical :: us
    real(dp), allocatable :: uvec(:), Z(:,:)
    integer :: j
    us = .false.
    if (present(usesvd)) us = usesvd
    if (us) then
      if (.not. present(U) .or. .not. present(V)) error stop "SVD factors missing"
      allocate(uvec(size(L,1)), Z(size(U,1),size(U,2)))
      uvec = matmul(U,matmul(transpose(V),x)) + e
      g = -matmul(V,matmul(transpose(U),w/uvec))
      Z = U
      do j = 1, size(Z,2)
        Z(:,j) = Z(:,j)*(sqrt(w)/uvec)
      end do
      H = matmul(V,matmul(matmul(transpose(Z),Z),transpose(V)))
    else
      allocate(uvec(size(L,1)), Z(size(L,1),size(L,2)))
      uvec = matmul(L,x) + e
      g = -matmul(transpose(L),w/uvec)
      Z = L
      do j = 1, size(Z,2)
        Z(:,j) = Z(:,j)*(sqrt(w)/uvec)
      end do
      H = matmul(transpose(Z),Z)
    end if
  end subroutine compute_grad_hessian

  subroutine truncated_svd(X,tol,Uf,Vf,rank,ok)
    real(dp), intent(in) :: X(:,:), tol
    real(dp), allocatable, intent(out) :: Uf(:,:), Vf(:,:)
    integer, intent(out) :: rank
    logical, intent(out) :: ok
    integer :: m,n,r,lda,ldu,ldvt,lwork,info,k,i
    real(dp), allocatable :: A(:,:), s(:), U(:,:), VT(:,:), work(:)
    integer, allocatable :: iwork(:)
    real(dp) :: wk(1)
    m = size(X,1); n = size(X,2); r = min(m,n)
    ok = .false.; rank = 0
    if (r < 2) return
    lda = max(1,m); ldu = max(1,m); ldvt = max(1,r)
    allocate(A(m,n),s(r),U(m,r),VT(r,n),iwork(8*r))
    A = X
    lwork = -1
    allocate(work(1))
    call dgesdd('S',m,n,A,lda,s,U,ldu,VT,ldvt,work,lwork,iwork,info)
    if (info /= 0) return
    wk(1) = work(1)
    deallocate(work)
    lwork = max(1,int(wk(1)))
    allocate(work(lwork))
    A = X
    call dgesdd('S',m,n,A,lda,s,U,ldu,VT,ldvt,work,lwork,iwork,info)
    if (info /= 0) return
    k = count(s > tol)
    if (k < 2) k = min(2,r)
    rank = k
    allocate(Uf(m,k),Vf(n,k))
    do i = 1, k
      Uf(:,i) = U(:,i)*sqrt(max(s(i),0.0_dp))
      Vf(:,i) = VT(i,:)*sqrt(max(s(i),0.0_dp))
    end do
    ok = .true.
  end subroutine truncated_svd

  function logspace(x,y,n) result(v)
    real(dp), intent(in) :: x,y
    integer, intent(in) :: n
    real(dp), allocatable :: v(:)
    integer :: i
    allocate(v(n))
    if (n == 1) then
      v(1) = x
    else
      do i = 1, n
        v(i) = 2.0_dp**(log(x)/log(2.0_dp) + real(i-1,dp)*(log(y/x)/log(2.0_dp))/real(n-1,dp))
      end do
    end if
  end function logspace

  subroutine set_seed(seed)
    integer, intent(in) :: seed
    integer :: n,i
    integer, allocatable :: put(:)
    call random_seed(size=n)
    allocate(put(n))
    do i = 1,n
      put(i) = modulo(seed + 104729*i, huge(1)-1)
      if (put(i) == 0) put(i) = i
    end do
    call random_seed(put=put)
  end subroutine set_seed
end module mixsqp_utils
