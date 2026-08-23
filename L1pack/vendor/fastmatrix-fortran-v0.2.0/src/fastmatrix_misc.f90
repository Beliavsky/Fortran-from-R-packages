module fastmatrix_misc
  use fastmatrix_base, only: dp, eye, solve_linear
  use fastmatrix_linalg, only: jacobi_eigen
  implicit none
  private
  public :: minkowski, floyd_warshall, bezier_curve, de_casteljau, krylov, array_mult
  public :: matrix_exp_sym, matrix_log_sym, matrix_power_sym, modified_cholesky
contains
  pure function minkowski(x,p) result(v)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: p
    real(dp) :: v, pp
    pp = 2.0_dp
    if (present(p)) pp = p
    if (pp > 1.0e100_dp) then
      v = maxval(abs(x))
    else
      v = sum(abs(x)**pp)**(1.0_dp/pp)
    end if
  end function

  subroutine floyd_warshall(a, dist)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: dist(:,:)
    integer :: i, j, k, n
    n = size(a,1)
    dist = a
    do k = 1, n
      do i = 1, n
        do j = 1, n
          dist(i,j) = min(dist(i,j), dist(i,k) + dist(k,j))
        end do
      end do
    end do
  end subroutine

  pure function de_casteljau(ctrl, t) result(pt)
    real(dp), intent(in) :: ctrl(:,:)
    real(dp), intent(in) :: t
    real(dp) :: pt(size(ctrl,2)), work(size(ctrl,1),size(ctrl,2))
    integer :: r, i, n
    n = size(ctrl,1)
    work = ctrl
    do r = 1, n-1
      do i = 1, n-r
        work(i,:) = (1.0_dp-t)*work(i,:) + t*work(i+1,:)
      end do
    end do
    pt = work(1,:)
  end function

  subroutine bezier_curve(ctrl, t, curve)
    real(dp), intent(in) :: ctrl(:,:), t(:)
    real(dp), intent(out) :: curve(size(t),size(ctrl,2))
    integer :: i
    do i=1,size(t)
      curve(i,:) = de_casteljau(ctrl,t(i))
    end do
  end subroutine

  subroutine krylov(a, x, k, kmat)
    real(dp), intent(in) :: a(:,:), x(:)
    integer, intent(in) :: k
    real(dp), intent(out) :: kmat(size(x),k)
    integer :: j
    kmat(:,1) = x
    do j=2,k
      kmat(:,j) = matmul(a,kmat(:,j-1))
    end do
  end subroutine

  pure function array_mult(a,b) result(c)
    real(dp), intent(in) :: a(:,:,:), b(:,:,:)
    real(dp) :: c(size(a,1),size(b,2),size(a,3))
    integer :: k
    do k=1,size(a,3)
      c(:,:,k)=matmul(a(:,:,k),b(:,:,k))
    end do
  end function

  subroutine matrix_exp_sym(a, f, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: f(:,:)
    integer, intent(out), optional :: info
    call spectral_function(a,1,f,info)
  end subroutine

  subroutine matrix_log_sym(a, f, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: f(:,:)
    integer, intent(out), optional :: info
    call spectral_function(a,2,f,info)
  end subroutine

  subroutine matrix_power_sym(a, power, f, info)
    real(dp), intent(in) :: a(:,:), power
    real(dp), intent(out) :: f(:,:)
    integer, intent(out), optional :: info
    real(dp) :: eval(size(a,1)), evec(size(a,1),size(a,1))
    integer :: i, ier
    call jacobi_eigen(a,eval,evec,info=ier)
    if (ier /= 0) then
      if (present(info)) info=ier
      f=0.0_dp
      return
    end if
    if (any(eval < 0.0_dp) .and. abs(power-nint(power)) > 1.0e-12_dp) then
      if (present(info)) info=-1
      f=0.0_dp
      return
    end if
    f=0.0_dp
    do i=1,size(eval)
      f=f+eval(i)**power*spread(evec(:,i),2,size(eval))*spread(evec(:,i),1,size(eval))
    end do
    if(present(info)) info=0
  end subroutine

  subroutine spectral_function(a, which, f, info)
    real(dp),intent(in)::a(:,:)
    integer,intent(in)::which
    real(dp),intent(out)::f(:,:)
    integer,intent(out),optional::info
    real(dp)::eval(size(a,1)),evec(size(a,1),size(a,1)),val
    integer::i,ier
    call jacobi_eigen(a,eval,evec,info=ier)
    if(ier/=0)then
      if(present(info))info=ier
      f=0.0_dp
      return
    end if
    if(which==2 .and. any(eval<=0.0_dp))then
      if(present(info))info=-1
      f=0.0_dp
      return
    end if
    f=0.0_dp
    do i=1,size(eval)
      if(which==1)then
        val=exp(eval(i))
      else
        val=log(eval(i))
      end if
      f=f+val*spread(evec(:,i),2,size(eval))*spread(evec(:,i),1,size(eval))
    end do
    if(present(info))info=0
  end subroutine

  subroutine modified_cholesky(a,l,d,delta,info)
    real(dp),intent(in)::a(:,:)
    real(dp),intent(out)::l(:,:),d(:)
    real(dp),intent(in),optional::delta
    integer,intent(out),optional::info
    real(dp)::s,del
    integer::n,i,j,k
    n=size(a,1)
    del=sqrt(epsilon(1.0_dp))*max(1.0_dp,maxval(abs(a)))
    if(present(delta))del=delta
    l=0.0_dp
    d=0.0_dp
    do i=1,n
      do j=1,i-1
        s=a(i,j)
        do k=1,j-1
          s=s-l(i,k)*d(k)*l(j,k)
        end do
        l(i,j)=s/d(j)
      end do
      s=a(i,i)
      do k=1,i-1
        s=s-l(i,k)**2*d(k)
      end do
      d(i)=max(abs(s),del)
      l(i,i)=1.0_dp
    end do
    if(present(info))info=0
  end subroutine
end module
