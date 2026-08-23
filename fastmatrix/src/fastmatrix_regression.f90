module fastmatrix_regression
  use fastmatrix_base, only: dp, inverse_matrix, solve_linear
  use fastmatrix_linalg, only: cg_solve
  implicit none
  private
  type, public :: ols_result
    real(dp), allocatable :: coefficients(:), fitted(:), residuals(:), cov_unscaled(:,:)
    real(dp) :: rss=0.0_dp
    integer :: rank=0, iterations=0
  end type
  public :: ols_fit, ols_fit_cg, ridge_fit
contains
  subroutine ols_fit(x,y,res,info)
    real(dp),intent(in)::x(:,:),y(:)
    type(ols_result),intent(out)::res
    integer,intent(out),optional::info
    real(dp),allocatable::xtx(:,:),xty(:)
    integer::p,ier
    p=size(x,2)
    allocate(xtx(p,p),xty(p),res%coefficients(p),res%fitted(size(y)),res%residuals(size(y)),res%cov_unscaled(p,p))
    xtx=matmul(transpose(x),x)
    xty=matmul(transpose(x),y)
    call solve_linear(xtx,xty,res%coefficients,ier)
    call inverse_matrix(xtx,res%cov_unscaled,ier)
    res%fitted=matmul(x,res%coefficients)
    res%residuals=y-res%fitted
    res%rss=sum(res%residuals**2)
    res%rank=p
    if(present(info))info=ier
  end subroutine
  subroutine ols_fit_cg(x,y,res,tol,maxiter)
    real(dp),intent(in)::x(:,:),y(:)
    type(ols_result),intent(out)::res
    real(dp),intent(in),optional::tol
    integer,intent(in),optional::maxiter
    real(dp),allocatable::xtx(:,:),xty(:)
    integer::p,it
    p=size(x,2)
    allocate(xtx(p,p),xty(p),res%coefficients(p),res%fitted(size(y)),res%residuals(size(y)))
    xtx=matmul(transpose(x),x)
    xty=matmul(transpose(x),y)
    call cg_solve(xtx,xty,res%coefficients,tol,maxiter,it)
    res%iterations=it
    res%fitted=matmul(x,res%coefficients)
    res%residuals=y-res%fitted
    res%rss=sum(res%residuals**2)
    res%rank=p
  end subroutine
  subroutine ridge_fit(x,y,lambda,coef)
    real(dp),intent(in)::x(:,:),y(:),lambda
    real(dp),intent(out)::coef(:)
    real(dp)::a(size(coef),size(coef)),b(size(coef))
    integer::i,info
    a=matmul(transpose(x),x)
    do i=1,size(coef)
    a(i,i)=a(i,i)+lambda
    end do
    b=matmul(transpose(x),y)
    call solve_linear(a,b,coef,info)
  end subroutine
end module
