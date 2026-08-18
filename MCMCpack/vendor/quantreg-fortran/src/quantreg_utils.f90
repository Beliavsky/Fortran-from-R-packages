! SPDX-License-Identifier: GPL-2.0-or-later
module quantreg_utils
  use quantreg_kinds, only : dp
  use quantreg_types, only : rq_result
  use quantreg_dense, only : rq_fit_fnb
  use quantreg_select, only : kuantiles
  use quantreg_linalg, only : symmetric_inverse
  use quantreg_bootstrap, only : seed_rng
  implicit none
  private
  public :: rq_fit_pfn, recursive_least_squares, combinations, random_exponential
contains

  subroutine rq_fit_pfn(x, y, tau, result, seed, mm_factor, max_bad_fixups, eps)
    real(dp), intent(in) :: x(:,:), y(:), tau
    type(rq_result), intent(out) :: result
    integer, intent(in), optional :: seed, max_bad_fixups
    real(dp), intent(in), optional :: mm_factor, eps
    integer :: n, p, m, mbad, ifix, ibad, badfix, ns, i, j, k
    integer, allocatable :: idx(:), sample_idx(:)
    logical, allocatable :: sl(:), su(:), submask(:)
    real(dp) :: mmf, tol, loq, hiq, mcap, u
    real(dp), allocatable :: xs(:,:), ys(:), gram(:,:), invgram(:,:), band(:), r(:), stdr(:), q(:)
    real(dp), allocatable :: xx(:,:), yy(:), globx(:), ghibx(:)
    type(rq_result) :: fit

    n = size(x,1)
    p = size(x,2)
    result%tau = tau
    if (size(y) /= n .or. tau <= 0.0_dp .or. tau >= 1.0_dp) then
      result%info = -1
      return
    end if
    mmf = 0.8_dp
    if (present(mm_factor)) mmf = mm_factor
    mbad = 3
    if (present(max_bad_fixups)) mbad = max_bad_fixups
    tol = 1.0e-6_dp
    if (present(eps)) tol = eps
    if (present(seed)) call seed_rng(seed)

    m = nint(sqrt(real(p,dp)) * real(n,dp)**(2.0_dp/3.0_dp))
    m = max(p+2, min(n,m))
    allocate(idx(n), sl(n), su(n), submask(n), r(n), band(n), stdr(n), q(2))
    allocate(gram(p,p), invgram(p,p), globx(p), ghibx(p))
    ifix = 0
    ibad = 0

    do
      ibad = ibad + 1
      if (m >= n) then
        call rq_fit_fnb(x,y,tau,fit,eps=tol)
        exit
      end if
      do i=1,n
        idx(i)=i
      end do
      do i=n,2,-1
        call random_number(u)
        j=min(i,int(u*real(i,dp))+1)
        k=idx(i)
        idx(i)=idx(j)
        idx(j)=k
      end do
      allocate(sample_idx(m), xs(m,p), ys(m))
      sample_idx=idx(1:m)
      do i=1,m
        xs(i,:)=x(sample_idx(i),:)
        ys(i)=y(sample_idx(i))
      end do
      call rq_fit_fnb(xs,ys,tau,fit,eps=tol)
      if (fit%info /= 0) then
        deallocate(sample_idx,xs,ys)
        m=min(n,2*m)
        cycle
      end if
      gram=matmul(transpose(xs),xs)
      call symmetric_inverse(gram,invgram,k)
      if (k /= 0) then
        deallocate(sample_idx,xs,ys)
        m=min(n,2*m)
        cycle
      end if
      do i=1,n
        band(i)=sqrt(max(tol,dot_product(x(i,:),matmul(invgram,x(i,:)))))
      end do
      r=y-matmul(x,fit%coefficients)
      mcap=mmf*real(m,dp)
      loq=max(1.0_dp/real(n,dp),tau-mcap/(2.0_dp*real(n,dp)))
      hiq=min(tau+mcap/(2.0_dp*real(n,dp)),real(n-1,dp)/real(n,dp))
      stdr=r/max(tol,band)
      call kuantiles(stdr,[loq,hiq],q,7)
      sl=r < band*q(1)
      su=r > band*q(2)
      badfix=0

      do
        ifix=ifix+1
        submask=.not.(su .or. sl)
        ns=count(submask)
        allocate(xx(ns+merge(1,0,any(sl))+merge(1,0,any(su)),p))
        allocate(yy(size(xx,1)))
        k=0
        do i=1,n
          if (submask(i)) then
            k=k+1
            xx(k,:)=x(i,:)
            yy(k)=y(i)
          end if
        end do
        if (any(sl)) then
          globx=0.0_dp
          do i=1,n
            if (sl(i)) globx=globx+x(i,:)
          end do
          k=k+1
          xx(k,:)=globx
          yy(k)=sum(y,mask=sl)
        end if
        if (any(su)) then
          ghibx=0.0_dp
          do i=1,n
            if (su(i)) ghibx=ghibx+x(i,:)
          end do
          k=k+1
          xx(k,:)=ghibx
          yy(k)=sum(y,mask=su)
        end if
        call rq_fit_fnb(xx,yy,tau,fit,eps=tol)
        deallocate(xx,yy)
        if (fit%info /= 0) then
          badfix=badfix+1
          m=min(n,2*m)
          exit
        end if
        r=y-matmul(x,fit%coefficients)
        if (.not. any((r < 0.0_dp) .and. su) .and. .not. any((r > 0.0_dp) .and. sl)) exit
        if (count(((r < 0.0_dp) .and. su) .or. ((r > 0.0_dp) .and. sl)) > int(0.1_dp*mcap)) then
          badfix=badfix+1
          m=min(n,2*m)
          exit
        end if
        su=su .and. .not.((r < 0.0_dp) .and. su)
        sl=sl .and. .not.((r > 0.0_dp) .and. sl)
        if (badfix >= mbad) exit
      end do
      deallocate(sample_idx,xs,ys)
      if (fit%info == 0 .and. .not. any((r < 0.0_dp) .and. su) &
          .and. .not. any((r > 0.0_dp) .and. sl)) exit
      if (m >= n) then
        call rq_fit_fnb(x,y,tau,fit,eps=tol)
        exit
      end if
    end do
    result=fit
    result%iterations=fit%iterations+ifix+ibad
  end subroutine rq_fit_pfn

  subroutine recursive_least_squares(x, y, beta_path, covariance, info)
    real(dp), intent(in) :: x(:,:), y(:)
    real(dp), intent(out) :: beta_path(:,:), covariance(:,:)
    integer, intent(out) :: info
    real(dp), allocatable :: gram(:,:), ginv(:,:), rhs(:), ax(:)
    real(dp) :: f, r
    integer :: n,p,i,j

    n=size(x,1)
    p=size(x,2)
    if (size(y)/=n .or. size(beta_path,1)/=p .or. size(beta_path,2)/=n) then
      info=-1
      return
    end if
    allocate(gram(p,p),ginv(p,p),rhs(p),ax(p))
    gram=matmul(transpose(x(1:p,:)),x(1:p,:))
    call symmetric_inverse(gram,ginv,info)
    if (info/=0) return
    rhs=matmul(transpose(x(1:p,:)),y(1:p))
    beta_path(:,p)=matmul(ginv,rhs)
    do j=1,p-1
      beta_path(:,j)=beta_path(:,p)
    end do
    do i=p+1,n
      ax=matmul(ginv,x(i,:))
      f=1.0_dp+dot_product(x(i,:),ax)
      r=(y(i)-dot_product(x(i,:),beta_path(:,i-1)))/f
      beta_path(:,i)=beta_path(:,i-1)+r*ax
      ginv=ginv-outer(ax,ax)/f
    end do
    covariance=ginv
  end subroutine recursive_least_squares

  subroutine combinations(n, k, a, count_out)
    integer, intent(in) :: n,k
    integer, allocatable, intent(out) :: a(:,:)
    integer, intent(out) :: count_out
    integer :: total, i
    integer, allocatable :: c(:)

    if (k<0 .or. k>n) then
      allocate(a(0,0))
      count_out=0
      return
    end if
    total=binomial(n,k)
    allocate(a(k,total),c(k))
    if (k==0) then
      count_out=1
      return
    end if
    do i=1,k
      c(i)=i
    end do
    count_out=0
    do
      count_out=count_out+1
      a(:,count_out)=c
      i=k
      do while (i>=1)
        if (c(i)/=n-k+i) exit
        i=i-1
      end do
      if (i==0) exit
      c(i)=c(i)+1
      do while (i<k)
        i=i+1
        c(i)=c(i-1)+1
      end do
    end do
  end subroutine combinations

  subroutine random_exponential(x, rate)
    real(dp), intent(out) :: x(:)
    real(dp), intent(in) :: rate
    real(dp) :: u
    integer :: i
    do i=1,size(x)
      call random_number(u)
      u=max(u,tiny(1.0_dp))
      x(i)=-log(u)/rate
    end do
  end subroutine random_exponential

  pure function outer(x,y) result(a)
    real(dp), intent(in) :: x(:),y(:)
    real(dp) :: a(size(x),size(y))
    integer :: i
    do i=1,size(x)
      a(i,:)=x(i)*y
    end do
  end function outer

  pure integer function binomial(n,k) result(v)
    integer, intent(in) :: n,k
    integer :: i,kk
    kk=min(k,n-k)
    v=1
    do i=1,kk
      v=v*(n-kk+i)/i
    end do
  end function binomial
end module quantreg_utils
